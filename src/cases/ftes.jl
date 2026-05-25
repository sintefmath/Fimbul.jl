import Jutul.CutCellMeshes: PlaneCut, PolygonalSurface, cut_mesh

function expand_ftes_fracture_property(values, fracture_disc::Dict{Symbol, Any})
    if ismissing(values)
        return missing
    elseif !(values isa AbstractVector)
        return values
    end

    n_regular = fracture_disc[:num_regular_cells]
    n_total = fracture_disc[:num_cells]
    if length(values) == n_total
        return values
    elseif length(values) == 1
        return fill(only(values), n_total)
    elseif length(values) != fracture_disc[:num_input_fractures]
        error("FTES fracture property vector must have length 1, $(fracture_disc[:num_input_fractures]), or $n_total. Got $(length(values)).")
    end

    cut_numbers = fracture_disc[:cut_numbers]
    expanded = fill(sum(values)/length(values), n_total)
    expanded[1:n_regular] .= values[cut_numbers]
    return expanded
end

function expand_ftes_fracture_properties(properties, fracture_disc::Dict{Symbol, Any})
    return Dict{Symbol, Any}(k => expand_ftes_fracture_property(v, fracture_disc) for (k, v) in pairs(properties))
end

function to_symbol_dict(values)
    return Dict{Symbol, Any}(pairs(values))
end

function to_nested_symbol_dict(values)
    if values isa AbstractDict
        return Dict{Symbol, Any}(
            key => to_nested_symbol_dict(value) for (key, value) in pairs(to_symbol_dict(values))
        )
    else
        return values
    end
end

function normalize_ftes_layer_properties(properties)
    return Dict{Symbol, Any}(
        key => (value isa AbstractVector ? value : [value]) for (key, value) in pairs(properties)
    )
end

function ftes(discretization::AbstractDict, parameters::AbstractDict, controls::AbstractDict;
    info_level = 0,
    )

    discretization = to_nested_symbol_dict(discretization)
    parameters = to_nested_symbol_dict(parameters)
    controls = to_nested_symbol_dict(controls)

    reservoir_disc = discretization[:reservoir]
    well_disc = discretization[:well]
    matrix_disc = reservoir_disc[:matrix]
    fracture_disc = reservoir_disc[:fractures]

    matrix_mesh = matrix_disc[:mesh]
    layers = matrix_disc[:layers]
    geo = matrix_disc[:geometry]
    fracture_mesh = fracture_disc[:mesh]

    matrix_parameters = parameters[:reservoir][:matrix]
    fracture_parameters = parameters[:reservoir][:fractures]
    injector_parameters = parameters[:well][:injector]
    producer_parameters = parameters[:well][:producer]
    model_parameters = parameters[:model]
    state0_parameters = parameters[:state0]
    bc_parameters = parameters[:boundary_conditions]

    info_level > 0 && @info "Setting up matrix, fracture, and well domains"
    matrix_domain = layered_reservoir_domain(matrix_mesh, layers, (; pairs(matrix_parameters)...))
    # volumes = max.(matrix_domain[:volumes], 1e-6)
    # matrix_domain[:volumes, Cells()] = volumes
    # areas = max.(matrix_domain[:areas], 1e-6)
    # matrix_domain[:areas, Faces()] = areas

    fracture_properties = expand_ftes_fracture_properties(fracture_parameters, fracture_disc)
    fracture_domain = JutulDarcy.fracture_domain(
        fracture_mesh,
        matrix_domain;
        pairs(fracture_properties)...,
    )
    # volumes = max.(fracture_domain[:volumes], 1e-6)
    # fracture_domain[:volumes, Cells()] = volumes
    # areas = max.(fracture_domain[:areas], 1e-6)
    # fracture_domain[:areas, Faces()] = areas

    injector_disc = well_disc[:injector]
    well_inj = setup_well(
        matrix_domain,
        injector_disc[:cells];
        pairs(injector_parameters)...,
    )

    producer_disc = well_disc[:producer]
    well_prod = setup_well(
        matrix_domain,
        producer_disc[:cells];
        pairs(producer_parameters)...,
        neighborship = producer_disc[:neighborship],
        perforation_cells_well = producer_disc[:perforation_cells_well],
        well_cell_centers = producer_disc[:well_cell_centers],
    )
    wells = [well_inj, well_prod]

    info_level > 0 && @info "Setting up DFM model"
    model_system = model_parameters[:system]
    model_kwargs = (k => v for (k, v) in pairs(model_parameters) if k != :system)
    model = JutulDarcy.setup_fractured_reservoir_model(
        matrix_domain,
        fracture_domain,
        model_system;
        wells = wells,
        model_kwargs...,
    )

    info_level > 0 && @info "Setting up initial state, boundary conditions, and well controls"
    pressure_fn = state0_parameters[:pressure]
    temperature_fn = state0_parameters[:temperature]
    z = geo.cell_centroids[3, :]
    state0 = setup_reservoir_state(model; Pressure = pressure_fn(z), Temperature = temperature_fn(z))

    bc_type = bc_parameters[:type]
    bc_type == :flow || error("Unsupported FTES boundary condition type: $bc_type")
    bc_cells = geo.boundary_neighbors
    z_bc = geo.cell_centroids[3, bc_cells]
    bc_pressure = bc_parameters[:pressure](z_bc)
    bc_temperature = bc_parameters[:temperature](z_bc)
    bc = flow_boundary_condition(bc_cells, matrix_domain, bc_pressure, bc_temperature)

    charge_period = controls[:schedule][:charge_period]
    discharge_period = controls[:schedule][:discharge_period]
    default_rate = controls[:default_rate]
    rate_charge = controls[:charge][:Injector][:target]
    if ismissing(rate_charge)
        rate_charge = default_rate[:multiplier]*Fimbul.scaled_rate(
            fracture_domain,
            wells,
            charge_period;
            mean_well_coordinate = default_rate[:mean_well_coordinate],
        )
    end
    rate_discharge = controls[:discharge][:Injector][:target]
    if ismissing(rate_discharge)
        rate_discharge = rate_charge
    end

    ρ = reservoir_model(model).system.rho_ref[1]
    p_ref = pressure_fn(0.0)
    function setup_ftes_control(control_spec, target)
        control_type = control_spec[:type]
        if control_type == :rate
            return InjectorControl(
                TotalRateTarget(target),
                control_spec[:composition];
                density = ρ,
                temperature = control_spec[:temperature],
            )
        elseif control_type == :bhp_fraction
            return ProducerControl(BottomHolePressureTarget(p_ref*control_spec[:target]))
        else
            error("Unsupported FTES control type: $control_type")
        end
    end

    control_charge = Dict{Symbol, Any}(
        :Injector => setup_ftes_control(controls[:charge][:Injector], rate_charge),
        :Producer => setup_ftes_control(controls[:charge][:Producer], controls[:charge][:Producer][:target]),
    )
    control_discharge = Dict{Symbol, Any}(
        :Injector => setup_ftes_control(controls[:discharge][:Injector], rate_discharge),
        :Producer => setup_ftes_control(controls[:discharge][:Producer], controls[:discharge][:Producer][:target]),
    )

    forces_charge = setup_reservoir_forces(model, bc = bc, control = control_charge)
    forces_discharge = setup_reservoir_forces(model, bc = bc, control = control_discharge)
    forces_rest = setup_reservoir_forces(model, bc = bc)

    dt, forces, timestamps = make_utes_schedule(
        forces_charge,
        forces_discharge,
        forces_rest;
        charge_period = charge_period,
        discharge_period = discharge_period,
        controls[:schedule][:utes_schedule_args]...,
    )

    info = Dict{Symbol, Any}(
        :discretization => discretization,
        :parameters => parameters,
        :controls => controls,
        :timestamps => timestamps,
    )
    haskey(discretization, :well_coordinates) && (info[:well_coordinates] = discretization[:well_coordinates])
    haskey(discretization, :fractures) && (info[:fractures] = discretization[:fractures])

    return JutulCase(model, dt, forces; state0 = state0, input_data = info)

end

function ftes(well_coordinates::Vector{Matrix{Float64}}, fractures::AbstractDict;
    depths = nothing,
    matrix_properties = Dict{Symbol, Any}(),
    fracture_properties = Dict{Symbol, Any}(),
    rate_charge = missing,
    rate_discharge = missing,
    temperature_charge = convert_to_si(95.0, :Celsius),
    temperature_discharge = convert_to_si(20.0, :Celsius),
    producer_bhp_fraction = 0.1,
    charge_period = ["April", "November"],
    discharge_period = ["December", "March"],
    utes_schedule_args = NamedTuple(),
    mesh_args = NamedTuple(),
    info_level = 0,
    )

    if isempty(fracture_properties)
        fractures = to_symbol_dict(fractures)
        fracture_properties = Dict{Symbol, Any}(
            :aperture => fractures[:aperture],
            :permeability => get(fractures, :permeability, missing),
            :porosity => fractures[:porosity],
        )
    else
        fractures = to_symbol_dict(fractures)
        fracture_properties = to_symbol_dict(fracture_properties)
    end

    discretization = ftes_discretization(
        well_coordinates,
        fractures;
        depths = depths,
        info_level = info_level,
        mesh_args...,
    )
    parameters = ftes_parameters(
        matrix_properties = matrix_properties,
        fracture_properties = fracture_properties,
    )
    controls = ftes_controls(
        rate_charge = rate_charge,
        rate_discharge = rate_discharge,
        temperature_charge = temperature_charge,
        temperature_discharge = temperature_discharge,
        producer_bhp_fraction = producer_bhp_fraction,
        charge_period = charge_period,
        discharge_period = discharge_period,
        utes_schedule_args = utes_schedule_args,
    )

    return ftes(discretization, parameters, controls; info_level = info_level)

end

function ftes(wells, fractures=Union{NamedTuple, Int}; kwargs...)
    well_coordinates, fractures = normalize_ftes_inputs(wells, fractures)
    return ftes(well_coordinates, fractures; kwargs...)

end

function ftes_discretization(well_coordinates::Vector{Matrix{Float64}}, fractures::AbstractDict;
    depths = nothing, info_level=0, hxy_min=missing, hz=missing, mesh_args...)

    fractures = to_symbol_dict(fractures)

    # Make constraints from well coordinates
    info_level > 0 && @info "Setting up wells and making mesh"
    collars = hcat([x[1:2, 1] for x in well_coordinates]...)
    Δx_min, Δx_max = Fimbul.min_max_distance(collars)
    hxy_min = ifelse(ismissing(hxy_min), Δx_min/4, hxy_min)
    r_given = filter(!isinf, fractures[:radius])
    r_max = ifelse(isempty(r_given), 0.0, maximum(r_given))
    well_offset = max(Δx_max/2, r_max-Δx_max/2)
    well_outline = Fimbul.offset_boundary(collars, well_offset; n=24)
    well_outline = hcat(well_outline, well_outline[:, 1]) # Close the loop
    
    collars = [permutedims([x[1] x[2]]).+hxy_min/2 for x in eachcol(collars)]
    cell_constraints = [x for x in collars]
    push!(cell_constraints, well_outline)
    # Determine layers including the well depths
    well_depth = maximum(maximum(x[3, :]) for x in well_coordinates)
    if isnothing(depths)
        depths = [0.0, well_depth+1e-2, well_depth*1.25]
    else
        extra_depths = [0.0, well_depth+1e-2, well_depth*1.25]
        for d in extra_depths
            if all(!isapprox.(depths, d, atol=0.5))
                push!(depths, d)
            end
        end
        depths = sort(depths)
    end
    # Generate mesh
    num_fractures = length(fractures[:normal])
    hz = ifelse(ismissing(hz), diff(depths)./[num_fractures*3, 2], hz)
    matrix_mesh, layers, _ = extruded_mesh(cell_constraints, depths;
        hxy_min=hxy_min, hz=hz, offset=well_offset*4, offset_rel=missing, mesh_args...)
    # Add fractures as polygonal disk cuts (radius from fracture dict)
    n_poly = 16
    θ_poly = range(0.0, 2π, length = n_poly + 1)[1:n_poly]
    function fracture_plane_basis(n_hat)
        ref = abs(n_hat[1]) < 0.9 ? [1.0, 0.0, 0.0] : [0.0, 1.0, 0.0]
        u = ref .- dot(ref, n_hat) .* n_hat; u ./= norm(u)
        v = cross(n_hat, u); v ./= norm(v)
        return u, v
    end
    cuts = []
    # Compute domain radius
    info_level > 0 && @info "Adding fractures to the mesh"
    geo = tpfv_geometry(matrix_mesh)
    for (normal, center, r) in zip(fractures[:normal], fractures[:centers], fractures[:radius])
        if isinf(r)
            push!(cuts, PlaneCut(center, normal))
        else
            n_hat = normal ./ norm(normal)
            u, v = fracture_plane_basis(n_hat)
            polygon = [Jutul.SVector{3, Float64}(
                center[1] + r * (cos(θ) * u[1] + sin(θ) * v[1]),
                center[2] + r * (cos(θ) * u[2] + sin(θ) * v[2]),
                center[3] + r * (cos(θ) * u[3] + sin(θ) * v[3])) for θ in θ_poly]
            push!(cuts, PolygonalSurface([polygon]))
        end
    end
    matrix_mesh, info = cut_mesh(matrix_mesh, cuts; extra_out = true, min_cut_fraction = 0.0)
    fracture_faces = findall(info[:face_index] .== 0)
    layers = layers[info[:cell_index]]
    # Generate embedded mesh for fractures
    fracture_mesh = Jutul.EmbeddedMeshes.EmbeddedMesh(matrix_mesh, fracture_faces)

    geo_matrix = tpfv_geometry(matrix_mesh)
    geo_fractures = tpfv_geometry(fracture_mesh)
    n_regular_fracture_cells = length(fracture_mesh.parent_faces)
    cut_numbers = info[:cut_no][fracture_mesh.parent_faces]

    reservoir_disc = Dict{Symbol, Any}()
    reservoir_disc[:matrix] = Dict(
        :mesh => matrix_mesh, :layers => layers, :geometry => geo_matrix)
    reservoir_disc[:fractures] = Dict(
        :mesh => fracture_mesh,
        :geometry => geo_fractures,
        :cut_numbers => cut_numbers,
        :num_cells => number_of_cells(fracture_mesh),
        :num_regular_cells => n_regular_fracture_cells,
        :num_input_fractures => length(fractures[:normal]),
    )

    # Generate matrix and fracture domains
    cells_inj = Jutul.find_enclosing_cells(matrix_mesh, permutedims(well_coordinates[1]), n=1_000_000)

    x_prod = [permutedims(x) for x in well_coordinates[2:end]]
    connectivity = zeros(Int, length(x_prod)+1, 2)
    connectivity[2:end, 1] .= 1
    cells_prod, wcells, neighborship = Fimbul.get_well_neighborship(
        matrix_mesh, x_prod, connectivity, geo_matrix; top_node=true, n=1_000_000)
    well_cell_centers = hcat([0; 0; 0], geo_matrix.cell_centroids[:, cells_prod])

    well_disc = Dict{Symbol, Any}()
    well_disc[:injector] = Dict(:cells => cells_inj)
    well_disc[:producer] = Dict(
        :cells => cells_prod,
        :neighborship => neighborship,
        :well_cell_centers => well_cell_centers,
        :perforation_cells_well => wcells[2:end],
        )

    discretization = Dict{Symbol, Any}()
    discretization[:reservoir] = reservoir_disc
    discretization[:well] = well_disc
    discretization[:well_coordinates] = well_coordinates
    discretization[:fractures] = fractures

    return discretization

end

function ftes_discretization(wells, fractures=Union{NamedTuple, Int}; kwargs...)
    well_coordinates, fractures = normalize_ftes_inputs(wells, fractures)
    return ftes_discretization(well_coordinates, fractures; kwargs...)
end

function ftes_parameters(;
    matrix_properties = Dict{Symbol, Any}(),
    fracture_properties = Dict{Symbol, Any}(),
    )

    default_matrix_properties = Dict{Symbol, Any}(
        :permeability => 1e-4si_unit(:darcy),
        :porosity => 0.01,
        :rock_thermal_conductivity => 2.5*watt/(meter*Kelvin),
    )
    matrix_properties = to_symbol_dict(matrix_properties)
    matrix_properties = merge(default_matrix_properties, matrix_properties)

    default_fracture_properties = Dict{Symbol, Any}(
        :aperture => 1.0e-4,
        :permeability => missing,
        :porosity => 0.5,
    )
    fracture_properties = to_symbol_dict(fracture_properties)
    fracture_properties = merge(default_fracture_properties, fracture_properties)

    matrix_properties = normalize_ftes_layer_properties(matrix_properties)

    p0(z) = 20si_unit(:atm)
    T0(z) = convert_to_si(10.0, :Celsius)

    parameters = Dict{Symbol, Any}()
    parameters[:reservoir] = Dict{Symbol, Any}(
        :matrix => matrix_properties,
        :fractures => fracture_properties,
    )
    parameters[:well] = Dict{Symbol, Any}(
        :injector => Dict{Symbol, Any}(
            :name => :Injector,
            :radius => 75e-3,
            :simple_well => false,
        ),
        :producer => Dict{Symbol, Any}(
            :name => :Producer,
            :radius => 75e-3,
            :use_top_node => true,
            :simple_well => false,
        ),
    )
    parameters[:model] = Dict{Symbol, Any}(
        :system => :geothermal,
        :block_backend => true,
    )
    parameters[:state0] = Dict{Symbol, Any}(
        :pressure => p0,
        :temperature => T0,
    )
    parameters[:boundary_conditions] = Dict{Symbol, Any}(
        :type => :flow,
        :pressure => p0,
        :temperature => T0,
    )
    return parameters

end

function ftes_controls(;
    rate_charge = missing,
    rate_discharge = missing,
    temperature_charge = convert_to_si(95.0, :Celsius),
    temperature_discharge = convert_to_si(20.0, :Celsius),
    producer_bhp_fraction = 0.1,
    charge_period = ["April", "November"],
    discharge_period = ["December", "March"],
    utes_schedule_args = NamedTuple(),
    )

    if !ismissing(rate_charge) && ismissing(rate_discharge)
        rate_discharge = rate_charge
    end

    controls = Dict{Symbol, Any}()
    controls[:default_rate] = Dict{Symbol, Any}(
        :type => :scaled_rate,
        :multiplier => 20_000,
        :mean_well_coordinate => true,
    )
    controls[:charge] = Dict{Symbol, Any}(
        :Injector => Dict{Symbol, Any}(
            :type => :rate,
            :target => rate_charge,
            :composition => [1.0],
            :temperature => temperature_charge,
        ),
        :Producer => Dict{Symbol, Any}(
            :type => :bhp_fraction,
            :target => producer_bhp_fraction,
        ),
    )
    controls[:discharge] = Dict{Symbol, Any}(
        :Injector => Dict{Symbol, Any}(
            :type => :rate,
            :target => rate_discharge,
            :composition => [1.0],
            :temperature => temperature_discharge,
        ),
        :Producer => Dict{Symbol, Any}(
            :type => :bhp_fraction,
            :target => producer_bhp_fraction,
        ),
    )
    controls[:schedule] = Dict{Symbol, Any}(
        :charge_period => charge_period,
        :discharge_period => discharge_period,
        :utes_schedule_args => utes_schedule_args,
    )

    return controls

end

function normalize_ftes_inputs(wells, fractures)
    well_coordinates = wells
    if well_coordinates isa NamedTuple
        if !haskey(well_coordinates, :num_producers) || !haskey(well_coordinates, :radius) || !haskey(well_coordinates, :depth)
            error("Named tuple defining wells must contain the keys: :num_producers, :radius, and :depth")
        end
        well_coordinates = setup_ftes_well_coordinates(wells.num_producers, wells.radius, wells.depth)
    end
    well_coordinates isa Vector{Matrix{Float64}} || error("well_coordinates must be a Vector of Matrix{Float64}")

    if fractures isa Int
        z_min = minimum(minimum(x[3, :] for x in well_coordinates))
        z_max = maximum(maximum(x[3, :] for x in well_coordinates))
        Δz = z_max - z_min
        z_min += Δz/8
        z_max -= Δz/8
        fractures = (num=fractures, z_min=z_min, z_max=z_max)
    end
    if fractures isa NamedTuple
        if !haskey(fractures, :num) || !haskey(fractures, :z_min) || !haskey(fractures, :z_max)
            error("Named tuple defining fractures must contain the keys: :num, :z_min, and :z_max")
        end
        num_fractures = fractures.num
        z_min = fractures.z_min
        z_max = fractures.z_max
        fractures = (; (k => v for (k, v) in pairs(fractures) if k ∉ (:num, :z_min, :z_max))...)
        fractures = setup_ftes_fractures(num_fractures, z_min, z_max; fractures...)
    end
    fractures isa AbstractDict || error("fractures must be an AbstractDict with Symbol keys")
    fractures = to_symbol_dict(fractures)
    required_keys = [:normal, :centers, :radius, :aperture, :porosity]
    for key in required_keys
        haskey(fractures, key) || error("fractures must contain the key: $key")
    end

    return well_coordinates, fractures
end

function setup_ftes_well_coordinates(num_producers::Int, radius::Float64, depth::Float64)
    # Place producers in a circle around the injector
    Δθ = 2π/(num_producers)
    r = radius
    producer_coordinates = [[r*cos(i*Δθ), r*sin(i*Δθ), 0.0] for i in 0:num_producers-1]
    wc = pushfirst!(producer_coordinates, [0.0, 0.0, 0.0]) # Add injector at the center

    well_coordinates = Vector{Matrix{Float64}}(undef, length(wc))
    for (wno, x) in enumerate(wc)
        x = repeat(x, 1, 2)
        x[3, 2] = depth
        well_coordinates[wno] = x
    end
    well_coordinates = [wc for wc in well_coordinates]

    return well_coordinates
end

function setup_ftes_fractures(num::Int, z_min::Float64, z_max::Float64;
    uniform_spacing = false,
    strike::Union{Float64, Tuple{Float64, Float64}}=(0.0, 5.0),
    dip::Union{Float64, Tuple{Float64, Float64}}=(0.0, 5.0),
    radius::Union{Float64, Tuple{Float64, Float64}}=Inf,
    aperture::Union{Float64, Tuple{Float64, Float64}}=1.0*1e-4,
    porosity::Union{Float64, Tuple{Float64, Float64}}=0.5,
    boundary_or_center=[0.0, 0.0],
    )

    if !(strike isa Tuple{Float64, Float64})
        strike = (strike, 0.0)
    end
    if !(dip isa Tuple{Float64, Float64})
        dip = (dip, 0.0)
    end
    if !(radius isa Tuple{Float64, Float64})
        radius = (radius, 0.0)
    end
    if !(aperture isa Tuple{Float64, Float64})
        aperture = (aperture, 0.0)
    end
    if !(porosity isa Tuple{Float64, Float64})
        porosity = (porosity, 0.0)
    end
    # Set strike and dip angles
    strike = strike[1] .+ randn(num).*strike[2]
    dip = dip[1] .+ randn(num).*dip[2]
    normal = [strike_dip_to_normal(s, d) for (s, d) in zip(strike, dip)]
    # Set radius, aperture, and porosity
    radius = radius[1] .+ randn(num).*radius[2]
    aperture = aperture[1] .+ randn(num).*aperture[2]
    porosity = porosity[1] .+ randn(num).*porosity[2]

    if uniform_spacing
        z = range(z_min, z_max, length=num)
    else
        z = z_min .+ rand(num) .* (z_max - z_min)
    end

    if boundary_or_center isa Vector
        centers = [vcat(boundary_or_center, z[i]) for i in 1:num]
    else
        x_min = minimum(boundary_or_center[1,:])
        x_max = maximum(boundary_or_center[1,:])
        y_min = minimum(boundary_or_center[2,:])
        y_max = maximum(boundary_or_center[2,:])
        centers = Vector{Vector{Float64}}()
        while length(centers) < num
            x = x_min + rand() * (x_max - x_min)
            y = y_min + rand() * (y_max - y_min)
            xy = [x, y]
            if Fimbul.point_in_polygon(xy, boundary_or_center)
                push!(centers, [x, y, z[length(centers)+1]])
            end
        end
    end

    fractures = Dict{Symbol, Any}()
    fractures[:normal] = normal
    fractures[:centers] = centers
    fractures[:radius] = radius
    fractures[:aperture] = aperture
    fractures[:porosity] = porosity

    return fractures

end
