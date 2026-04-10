import Jutul.CutCellMeshes: PlaneCut, PolygonalSurface, cut_mesh

function ftes(well_coordinates::Vector{Matrix{Float64}}, fractures::Dict{Symbol, Any};
    depths = nothing,
    matrix_properties = NamedTuple(),
    fracture_properties = NamedTuple(),
    rate_charge = missing,
    rate_discharge = missing,
    temperature_charge = convert_to_si(95.0, :Celsius),
    temperature_discharge = convert_to_si(10.0, :Celsius),
    charge_period = ["April", "November"],
    discharge_period = ["December", "March"],
    utes_schedule_args = NamedTuple(),
    mesh_args = NamedTuple(),
    )

    @warn "This function is under development and currently requires the `dfm` \
    branch of JutulDarcy."

    # Make constraints from well coordinates
    collars = hcat([x[1:2, 1] for x in well_coordinates]...)
    Δx_min, Δx_max = Fimbul.min_max_distance(collars)
    hxy_min = Δx_min/3
    well_outline = Fimbul.offset_boundary(collars, Δx_max; n=24)
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
    hz = diff(depths)./[num_fractures*3, 2]
    matrix_mesh, layers, _ = extruded_mesh(cell_constraints, depths;
        hxy_min=hxy_min, hz=hz, offset_rel=2.5, mesh_args...)
    # Add fractures
    fracture_faces = Int[]
    for (i, (normal, center)) in enumerate(zip(fractures[:normal], fractures[:centers]))
        is_frac = falses(number_of_faces(matrix_mesh))
        is_frac[fracture_faces] .= true
        plane = PlaneCut(center, normal)
        matrix_mesh, info = cut_mesh(matrix_mesh, plane; min_cut_fraction = 0.0, extra_out=true)
        # Update fracture face vector
        face_index = filter(x->x>0, info[:face_index])
        fracture_faces = findall(is_frac[face_index])
        new_faces = findall(info[:face_index] .== 0)
        append!(fracture_faces, new_faces)
        # Update layer vector
        layers = layers[info[:cell_index]]
    end
    # Generate embedded mesh for fractures
    fracture_mesh = Jutul.EmbeddedMeshes.EmbeddedMesh(matrix_mesh, fracture_faces)

    geo = tpfv_geometry(matrix_mesh)
    # Generate matrix and fracture domains
    if isempty(matrix_properties)
        matrix_properties = (permeability=1e-4si_unit(:darcy), porosity=0.01)
    end
    matrix_domain = layered_reservoir_domain(matrix_mesh, layers, matrix_properties)
    # matrix_domain = reservoir_domain(matrix_mesh;
    #     permeability=marix_permeability, porosity=matrix_porosity)
    fracture_domain = JutulDarcy.fracture_domain(fracture_mesh;
        aperture=fractures[:aperture][1],
        porosity=fractures[:porosity][1],
        matrix_faces=fracture_faces)
    
    frac_args = (WI = missing,)
    cells = Jutul.find_enclosing_cells(matrix_mesh, permutedims(well_coordinates[1]), n=1_000_000)
    well_inj = setup_well(matrix_domain, cells, fracture_domain;
        name=:Injector, radius=75e-3, frac_args=frac_args)

    x_prod = [permutedims(x) for x in well_coordinates[2:end]]
    connectivity = zeros(Int, length(x_prod)+1, 2)
    connectivity[2:end, 1] .= 1
    println(connectivity)
    cells, wcells, neighborship = Fimbul.get_well_neighborship(
        matrix_mesh, x_prod, connectivity, geo; top_node=true, n=1_000_000)
    well_cell_centers = hcat([0; 0; 0], geo.cell_centroids[:, cells])

    well_prod = setup_well(matrix_domain, cells, fracture_domain;
        name=:Producer,
        radius=75e-3,
        frac_args=frac_args,
        neighborship=neighborship,
        perforation_cells_well=wcells[2:end],
        well_cell_centers=well_cell_centers,
        use_top_node=true,
        )
    
    wells = [well_inj, well_prod]
    
    model = JutulDarcy.setup_reservoir_model(matrix_domain, fracture_domain, :geothermal;
        wells=wells, block_backend=true)
    ρ = reservoir_model(model).system.rho_ref[1]
    dpdz = gravity_constant * ρ
    p0(z) = 20si_unit(:atm)# .+ dpdz.*z
    dTdz = 0.03si_unit(:Kelvin)/si_unit(:meter)
    T0(z) = convert_to_si(10.0, :Celsius)# .+ dTdz.*z
    z = geo.cell_centroids[3, :]
    # Set up initial state
    state0 = setup_reservoir_state(model; Pressure=p0(z), Temperature=T0(z))
    # Set up boundary conditions
    bc_cells = geo.boundary_neighbors
    bc_temperature = state0[:Reservoir][:Temperature][bc_cells]
    bc_pressure = state0[:Reservoir][:Pressure][bc_cells]
    bc = flow_boundary_condition(bc_cells, matrix_domain, bc_pressure, bc_temperature)

    if ismissing(rate_charge)
        rate_charge = 20_000*Fimbul.scaled_rate(
            fracture_domain, wells, charge_period; mean_well_coordinate=true)
    end
    if ismissing(rate_discharge)
        rate_discharge = rate_charge
    end
    ctrl_prod = ProducerControl(BottomHolePressureTarget(p0(0.0)*0.1))
    ctrl_inj = InjectorControl(TotalRateTarget(rate_charge), [1.0];
        density=ρ, temperature=temperature_charge)
    control_charge = Dict()
    for w in wells
        name = w.representation.name
        if name == :Injector
            control_charge[name] = ctrl_inj
        else
            control_charge[name] = ctrl_prod
        end
    end
    forces_charge = setup_reservoir_forces(model, bc=bc, control=control_charge)
    
    ctrl_inj = InjectorControl(TotalRateTarget(rate_discharge), [1.0];
        density=ρ, temperature=temperature_discharge)
    control_discharge = Dict()
    for w in wells
        name = w.representation.name
        if name == :Injector
            control_discharge[name] = ctrl_inj
        else
            control_discharge[name] = ctrl_prod
        end
    end
    forces_discharge = setup_reservoir_forces(model, bc=bc, control=control_discharge)

    forces_rest = setup_reservoir_forces(model, bc=bc)

    dt, forces = make_utes_schedule(forces_charge, forces_discharge, forces_rest;
        charge_period=charge_period, discharge_period=discharge_period,
        utes_schedule_args...)

    info = Dict{Symbol, Any}()
    info[:well_coordinates] = well_coordinates
    info[:fractures] = fractures
    case = JutulCase(model, dt, forces; state0=state0, input_data=info)

    return case

end

function ftes(wells, fractures=Union{NamedTuple, Int}; kwargs...)

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
        # Remove these keys from the named tuple before passing to the setup function
        fractures = (; (k => v for (k, v) in pairs(fractures) if k ∉ (:num, :z_min, :z_max))...)
        println(fractures)
        fractures = setup_ftes_fractures(num_fractures, z_min, z_max; fractures...)
    end
    fractures isa Dict{Symbol, Any} || error("fractures must be a Dict{Symbol, Any}")
    required_keys = [:normal, :centers, :radius, :aperture, :porosity]
    for key in required_keys
        haskey(fractures, key) || error("fractures must contain the key: $key")
    end

    return ftes(well_coordinates, fractures; kwargs...)

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

    z = z_min .+ rand(num) .* (z_max - z_min)

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
