"""
    btes(; <keyword arguments>)

Setup function for borehole thermal energy storage (BTES) system.

# Keyword arguments
- `field = missing`: Explicit wells for the system. If given, this overrides
  `num_wells`, `num_sectors` and `well_spacing`. The wells are given using
  the following structure:
  - A full field is represented as `[sector_1, sector_2, ..., sector_n]`.
  - `sector_k` contains all wells in sector `k`, represented as
    `[well_1, well_2, ..., well_nk]`.
  - Each `well_l` is a `3 x m` matrix containing the coordinates of the well
    trajectory.
  Wells within a sector are coupled in series, in the order they appear in
  `sector_k`: `well_1` is charged first and discharged last, `well_nk` is
  charged last and discharged first. A single well (a `3 x m` matrix) can
  also be passed directly as `field`, representing the special case of a
  field with a single sector containing a single well.
- `pattern = :sunflower`: Well placement pattern used when `field` is not
  given. One of `:sunflower`, `:rectangular`, `:circular` or `:polygonal`.
- `num_sides = 6`: Number of sides of the polygon when `pattern = :polygonal`.
- `num_wells = 48`: Number of wells in the BTES system.
- `num_sections = 6`: Number of sections in the BTES system. The system is
  divided into equal circle sectors, and all wells in each sector are coupled in series.
- `well_spacing = 5.0`: Horizontal spacing between wells [m].
- `depths = [0.0, 0.5, 50, 65]`: Depths delineating geological layers [m].
- `well_layers = [1, 2]`: Layers in which the wells are placed
- `density = [30, 2580, 2580]: Rock density in the layers [kg/m³].
- `thermal_conductivity = [0.034, 3.7, 3.7]: Thermal conductivity in the layers [W/(m⋅K)].
- `heat_capacity = [1500, 900, 900]`: Heat capacity in the layers [J/(kg⋅K)].
- `geothermal_gradient = 0.03 K/m`: Geothermal gradient [K/m].
- `temperature_charge = 90 °C/363.15 K`: Injection temperature during charging [K].
- `temperature_discharge = 10 °C/283.15 K`: Injection temperature during discharging [K].
- `rate_charge = 0.5 l/s`: Injection rate during charging [m³/s].
- `rate_discharge = rate_charge`: Injection rate during discharging [m³/s].
- `temperature_surface = 10 °C/283.15 K`: Temperature at the surface [K].
- `num_years = 5`: Number of years to run the simulation.
- `charge_period = ["June", "September"]`: Period during which the system is charged.
- `discharge_period = ["December", "March"]`: Period during which the system is discharged.
- `report_interval = 14 day`: Reporting interval for the simulation.
- `utes_schedule_args = NamedTuple()`: Additional arguments for the UTES schedule.
- `n_z = [3, 8, 3]`: Number of layers in the vertical direction for each layer.
- `n_xy = 3`: Number of layers in the horizontal direction for each layer.
- `mesh_args = NamedTuple()`: Additional arguments for the mesh generation.
"""
function btes(;
    field = missing,
    pattern = :sunflower,
    num_sides = 6,
    num_wells = 48,
    num_sectors = 6,
    well_spacing = 5.0,
    depths = [0.0, 0.5, 50, 65],
    well_layers = [1, 2],
    density = [30, 2580, 2580]*kilogram/meter^3,
    thermal_conductivity = [0.034, 3.7, 3.7]*watt/meter/Kelvin,
    heat_capacity = [1500, 900, 900]*joule/kilogram/Kelvin,
    geothermal_gradient = 0.03Kelvin/meter,
    temperature_charge = convert_to_si(90.0, :Celsius),
    temperature_discharge = convert_to_si(10.0, :Celsius),
    rate_charge = 0.5litre/second,
    rate_discharge = rate_charge,
    temperature_surface = convert_to_si(10.0, :Celsius),
    num_years = 4,
    charge_period = ["June", "September"],
    discharge_period = ["December", "March"],
    report_interval = 14day,
    utes_schedule_args = NamedTuple(),
    n_z = [3, 8, 3],
    n_xy = 3,
    mesh_args = NamedTuple(),
    )

    if ismissing(field)
        if pattern == :sunflower
            field = sunflower_pattern(num_wells, well_spacing; num_sectors = num_sectors, depths = depths)
        elseif pattern == :rectangular
            field = rectangular_pattern(num_wells, well_spacing; num_sectors = num_sectors, depths = depths)
        elseif pattern == :circular
            field = circular_pattern(num_wells, well_spacing; num_sectors = num_sectors, depths = depths)
        elseif pattern == :polygonal
            field = polygonal_pattern(num_wells, well_spacing, num_sides; num_sectors = num_sectors, depths = depths)
        else
            error("Unknown pattern: $pattern. Supported patterns are :sunflower, :rectangular, :circular and :polygonal.")
        end
    elseif field isa AbstractMatrix
        # Special case: a single well, given as a single 3 x m matrix
        field = [[field]]
    end

    well_coords_3d = vcat(field...)
    num_wells = length(well_coords_3d)

    # ## Create mesh
    # Use the (x,y) projection of each well trajectory as a mesh constraint
    well_coordinates = [wc[1:2, :] for wc in well_coords_3d]
    println(well_coordinates)
    hz = diff(depths)./n_z
    hxy = well_spacing/n_xy

    # ## Set up model
    # Set up reservoir domain with rock properties similar to that of granite,
    # with a styrofoam layer on top                
    domain, layers, metrics = layered_reservoir_domain(well_coordinates, depths,
        (
            rock_density = density,
            rock_thermal_conductivity = thermal_conductivity,
            rock_heat_capacity = heat_capacity
        );
        mesh_args = (; hxy_min = hxy, hz = hz, mesh_args...),
        permeability = 1e-6darcy,
        porosity = 0.01,
        component_heat_capacity = 4.278e3joule/kilogram/Kelvin,
    )
    mesh = physical_representation(domain)
    # Set up BTES wells
    hxy_min = metrics.hxy_min
    well_models = []
    well_names = Symbol[]
    nl = length(layers)
    geo = tpfv_geometry(mesh)

    for (wno, wc) in enumerate(well_coords_3d)
        name = Symbol("B$wno")
        println("Adding well $name ($wno/$num_wells)")
        cells = Jutul.find_enclosing_cells(mesh, permutedims(wc), n = 100)
        filter!(c -> layers[c] ∈ well_layers, cells)
        w_sup, w_ret = setup_btes_well(domain, cells, name=name, closed_loop_type=:u1)
        push!(well_models, w_sup, w_ret)
        push!(well_names, name)
    end

    # Make the model
    # Closed-loop wells must not be shut: The supply and return are hydraulically
    # connected, and shutting one of them gives a dead end.
    model = setup_reservoir_model(
        domain, :geothermal,
        wells = well_models,
        can_shut_wells = false,
    );

    # ## Set up initial state and boundary conditions
    geo = tpfv_geometry(mesh)
    z_bc = geo.boundary_centroids[3, :]
    bottom = map(v -> isapprox(v, maximum(z_bc)), z_bc)
    # Define pressure and temperature profiles
    rho = reservoir_model(model).system.rho_ref[1]
    dpdz = rho*gravity_constant
    dTdz = geothermal_gradient
    T = z -> temperature_surface .+ dTdz*z
    p = z -> 5atm .+ dpdz.*z
    # Set initial conditions
    z_cells = geo.cell_centroids[3, :]
    z_hat = z_cells .- minimum(z_cells)
    state0 = setup_reservoir_state(model,
        Pressure = p(z_hat),
        Temperature = T(z_hat)
    );
    # Set boundary conditions
    z_bc = z_bc[.!bottom]
    z_hat = z_bc .- minimum(z_bc)
    bc_cells = geo.boundary_neighbors[.!bottom]
    bc = flow_boundary_condition(bc_cells, domain, p(z_hat), T(z_hat));

    # Group supply well names by sector, preserving the order given in `field`
    wells_per_sector = Vector{Vector{Symbol}}()
    wtot = 0
    for sector in field
        sw = [Symbol(well_names[wtot + l], "_supply") for l in 1:length(sector)]
        wtot += length(sector)
        push!(wells_per_sector, sw)
    end
    control_charge, control_discharge, sectors = setup_controls(model, wells_per_sector,
        rate_charge, rate_discharge, temperature_charge, temperature_discharge);
    
    forces_charge = setup_reservoir_forces(model, control=control_charge, bc=bc)
    forces_discharge = setup_reservoir_forces(model, control=control_discharge, bc=bc);
    forces_rest = setup_reservoir_forces(model, bc=bc)
    # Make schedule
    dt, forces, timestamps = make_utes_schedule(
        forces_charge, forces_discharge, forces_rest;
        charge_period = charge_period,
        discharge_period = discharge_period,
        num_years = num_years,
        report_interval = report_interval,
        utes_schedule_args...,
    )

    # ## Useful case info
    info = Dict()
    info[:description] = "Borehole thermal energy storage (BTES) case set up using Fimbul.btes()"
    info[:sectors] = sectors
    info[:timestamps] = timestamps

    # ## Assemble and return model
    case = JutulCase(model, dt, forces, state0 = state0, input_data = info)
    return case

end

# ## Patterns
#
# Each pattern function takes the number of wells and the approximate
# spacing between neighboring wells, and returns a full field

function sunflower_pattern(num_wells, spacing; num_sectors = 6, depths = [0.0, 0.5, 50, 65])
    xy = fibonacci_pattern_2d(num_wells; spacing = spacing)
    return field_from_points(xy, num_sectors, depths)
end

function rectangular_pattern(num_wells, spacing; num_sectors = 6, depths = [0.0, 0.5, 50, 65], jitter = 0.003)
    nx = max(1, round(Int, sqrt(num_wells)))
    ny = ceil(Int, num_wells/nx)
    xs = ((0:nx-1) .- (nx-1)/2).*spacing
    ys = ((0:ny-1) .- (ny-1)/2).*spacing
    xy = Matrix{Float64}(undef, 2, nx*ny)
    k = 0
    for y in ys, x in xs
        k += 1
        xy[:, k] = [x, y]
    end
    xy = xy[:, 1:num_wells]
    # Perturb points slightly to avoid exact co-circularity, which causes an
    # ambiguous Delaunay triangulation (and hence irregular mesh cells) for a
    # perfectly regular grid.
    return field_from_points(xy, num_sectors, depths)
end

function circular_pattern(num_wells, spacing; num_sectors = 6, depths = [0.0, 0.5, 50, 65])
    points = Vector{Vector{Float64}}()
    push!(points, [0.0, 0.0])
    ring = 1
    while length(points) < num_wells
        r = ring*spacing
        n_ring = min(max(1, round(Int, 2π*r/spacing)), num_wells - length(points))
        for k in 0:n_ring-1
            θ = 2π*k/n_ring
            push!(points, [r*cos(θ), r*sin(θ)])
        end
        ring += 1
    end
    xy = hcat(points...)
    return field_from_points(xy, num_sectors, depths)
end

function field_from_points(xy::AbstractMatrix, num_sectors::Int, depths::AbstractVector)

    num_wells = size(xy, 2)
    well_coords = Vector{Matrix{Float64}}(undef, num_wells)
    for i in 1:num_wells
        x_top, y_top, z_top = xy[1, i], xy[2, i], 0.0 + 1e-3
        x_bottom, y_bottom, z_bottom = xy[1, i], xy[2, i], depths[end] - 1e-3
        well_coords[i] = permutedims([x_top y_top z_top; x_bottom y_bottom z_bottom])
    end
    sector_indices = group_into_sectors(xy, num_sectors)

    return [well_coords[idx] for idx in sector_indices]

end

function group_into_sectors(xy::AbstractMatrix, num_sectors::Int)

    n = size(xy, 2)
    r = sqrt.(xy[1,:].^2 .+ xy[2,:].^2)
    θ = atan.(xy[2,:], xy[1,:]) .+ π
    order_θ = sortperm(θ)
    wells_per_sector = div(n, num_sectors)
    rem = n - wells_per_sector*num_sectors
    wells_per_sector = fill(wells_per_sector, num_sectors)
    wells_per_sector[1:rem] .+= 1

    sector_indices = Vector{Vector{Int}}()
    wtot = 0
    for sno in 1:num_sectors
        idx = order_θ[(1:wells_per_sector[sno]) .+ wtot]
        wtot += wells_per_sector[sno]
        order_r = sortperm(r[idx])
        push!(sector_indices, idx[order_r])
    end

    return sector_indices

end

#=
function setup_controls(model, number_of_sectors,
    rate_charge, rate_discharge, temperature_charge, temperature_discharge)

    rho = reservoir_model(model).system.rho_ref[1]
    rate_target = TotalRateTarget(rate_charge)
    ctrl_charge = InjectorControl(rate_target, [1.0],
        density=rho, temperature=temperature_charge)
    rate_target = TotalRateTarget(rate_discharge)
    ctrl_discharge = InjectorControl(rate_target, [1.0],
        density=rho, temperature=temperature_discharge);
    # BHP control for return side
    bhp_target = BottomHolePressureTarget(1.0si_unit(:atm))
    ctrl_ret = ProducerControl(bhp_target);
    # Set up forces
    control_charge = Dict()
    control_discharge = Dict()

    msh = physical_representation(reservoir_model(model).data_domain)
    geo = tpfv_geometry(msh)
    xy = geo.cell_centroids[1:2,:]
    # map from (x,y) to polar coordiates
    r = sqrt.(xy[1,:].^2 .+ xy[2,:].^2)
    θ = atan.(xy[2,:], xy[1,:]) .+ π
    assigned = []
    # Closed-loop wells, each with a supply and a return part
    wells = [k for (k, m) in pairs(model.models) if is_closed_loop_well(m.data_domain)]
    get_supply = (well) -> first(closed_loop_well_names(model.models[well].data_domain))
    get_return = (well) -> last(closed_loop_well_names(model.models[well].data_domain))
    sectors = Dict()
    wells_per_sector = div(length(wells), number_of_sectors)
    rem = length(wells) - wells_per_sector*number_of_sectors
    wells_per_sector = fill(wells_per_sector, number_of_sectors)
    wells_per_sector[1:rem] .+= 1
    wmodels = [model.models[well] for well in wells]
    wc = [wm.domain.representation.perforations.reservoir[1] for wm in wmodels]
    θ = θ[wc]
    order_θ = sortperm(θ)
    wtot = 0
    for sno in 1:number_of_sectors
        wno = (1:wells_per_sector[sno]) .+ wtot
        wtot += wells_per_sector[sno]
        sw = wells[order_θ[wno]]
        well_radii = r[wc[order_θ[wno]]]
        order = sortperm(well_radii)
        sec_wells = Symbol[]
        for (k, wno) in enumerate(order)
            well_sup = get_supply(sw[wno])
            well_ret = get_return(sw[wno])
            @assert well_sup ∉ assigned
            @assert well_ret ∉ assigned
            if k == 1
                # Water is injected into innermost well during charging
                control_charge[well_sup] = ctrl_charge
                # Discharging runs from outer to inner
                well_prev = get_return(sw[order[k+1]])
                target = JutulDarcy.ReinjectionTarget([well_prev])
                ctrl = InjectorControl(target, [1.0],
                    density=rho, temperature=NaN; check=false)
                control_discharge[well_sup] = ctrl
            elseif k == length(sw)
                # Water is injected into outermost well during discharging
                control_discharge[well_sup] = ctrl_discharge
                # Charging runs from inner to outer
                well_prev = get_return(sw[order[k-1]])
                target = JutulDarcy.ReinjectionTarget([well_prev])
                ctrl = InjectorControl(target, [1.0],
                    density=rho, temperature=NaN; check=false)
                control_charge[well_sup] = ctrl
            else
                # Charging runs from inner to outer
                well_prev = get_return(sw[order[k-1]])
                target = JutulDarcy.ReinjectionTarget([well_prev])
                ctrl = InjectorControl(target, [1.0],
                    density=rho, temperature=NaN; check=false)
                control_charge[well_sup] = ctrl
                # Discharging runs from outer to inner
                well_prev = get_return(sw[order[k+1]])
                target = JutulDarcy.ReinjectionTarget([well_prev])
                ctrl = InjectorControl(target, [1.0],
                    density=rho, temperature=NaN; check=false)
                control_discharge[well_sup] = ctrl
            end
            control_charge[well_ret] = ctrl_ret
            control_discharge[well_ret] = ctrl_ret
            push!(assigned, well_sup, well_ret)
            push!(sec_wells, sw[wno])
        end
        sectors[Symbol("S$sno")] = sec_wells
    end

    @assert sort(assigned) == sort(well_symbols(model))

    return control_charge, control_discharge, sectors

end
=#

function setup_controls(model, wells_per_sector::AbstractVector{<:AbstractVector{Symbol}},
    rate_charge, rate_discharge, temperature_charge, temperature_discharge)

    rho = reservoir_model(model).system.rho_ref[1]
    rate_target = TotalRateTarget(rate_charge)
    ctrl_charge = InjectorControl(rate_target, [1.0],
        density=rho, temperature=temperature_charge)
    rate_target = TotalRateTarget(rate_discharge)
    ctrl_discharge = InjectorControl(rate_target, [1.0],
        density=rho, temperature=temperature_discharge);
    # BHP control for return side
    bhp_target = BottomHolePressureTarget(1.0si_unit(:atm))
    ctrl_ret = ProducerControl(bhp_target);
    # Set up forces
    control_charge = Dict()
    control_discharge = Dict()
    assigned = []
    get_return = (well) -> Symbol(replace(String(well), "_supply" => "_return"))
    sectors = Dict()

    for (sno, sw) in enumerate(wells_per_sector)
        sec_wells = Symbol[]
        for (k, well_sup) in enumerate(sw)
            well_ret = get_return(well_sup)
            @assert well_sup ∉ assigned
            @assert well_ret ∉ assigned
            if length(sw) == 1
                # Single well in sector: charge and discharge it directly
                control_charge[well_sup] = ctrl_charge
                control_discharge[well_sup] = ctrl_discharge
            elseif k == 1
                # Water is injected into innermost well during charging
                control_charge[well_sup] = ctrl_charge
                # Discharging runs from outer to inner
                well_prev = get_return(sw[k+1])
                target = JutulDarcy.ReinjectionTarget([well_prev])
                ctrl = InjectorControl(target, [1.0],
                    density=rho, temperature=NaN; check=false)
                control_discharge[well_sup] = ctrl
            elseif k == length(sw)
                # Water is injected into outermost well during discharging
                control_discharge[well_sup] = ctrl_discharge
                # Charging runs from inner to outer
                well_prev = get_return(sw[k-1])
                target = JutulDarcy.ReinjectionTarget([well_prev])
                ctrl = InjectorControl(target, [1.0],
                    density=rho, temperature=NaN; check=false)
                control_charge[well_sup] = ctrl
            else
                # Charging runs from inner to outer
                well_prev = get_return(sw[k-1])
                target = JutulDarcy.ReinjectionTarget([well_prev])
                ctrl = InjectorControl(target, [1.0],
                    density=rho, temperature=NaN; check=false)
                control_charge[well_sup] = ctrl
                # Discharging runs from outer to inner
                well_prev = get_return(sw[k+1])
                target = JutulDarcy.ReinjectionTarget([well_prev])
                ctrl = InjectorControl(target, [1.0],
                    density=rho, temperature=NaN; check=false)
                control_discharge[well_sup] = ctrl
            end
            control_charge[well_ret] = ctrl_ret
            control_discharge[well_ret] = ctrl_ret
            push!(assigned, well_sup, well_ret)
            push!(sec_wells, well_sup, well_ret)
        end
        sectors[Symbol("S$sno")] = sec_wells
    end

    @assert sort(assigned) == sort(well_symbols(model))

    return control_charge, control_discharge, sectors

end
