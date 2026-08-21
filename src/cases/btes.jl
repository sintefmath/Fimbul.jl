"""
    btes(pattern::Symbol = :sunflower; <keyword arguments>)
    btes(field; <keyword arguments>)

Setup function for borehole thermal energy storage (BTES) system.

The first form places the wells using one of the built-in patterns. The second
form takes an explicit well field, and does not accept `num_sides`, `num_wells`
or `num_sectors`, since the placement and grouping are given by `field` itself.

# Arguments
- `pattern = :sunflower`: Well placement pattern. One of `:sunflower`,
  `:rectangular`, `:circular` or `:polygonal`.
- `field`: Explicit wells for the system, given using the following structure:
  - A full field is represented as `[sector_1, sector_2, ..., sector_n]`.
  - `sector_k` contains all wells in sector `k`, represented as
    `[well_1, well_2, ..., well_nk]`.
  - Each `well_l` is a `3 x m` matrix containing the coordinates of the well
    trajectory.
  Wells within a sector are coupled in series, in the order they appear in
  `sector_k`: `well_1` is always charged first and `well_nk` charged last.
  Discharge order depends on `reversed_discharge`. A single well (a `3 x m`
  matrix) can also be passed directly as `field`, representing the special
  case of a field with a single sector containing a single well.

# Keyword arguments
- `num_sides = 6`: Number of sides of the polygon when `pattern = :polygonal`.
  Only used with the `pattern` form.
- `num_wells = 48`: Number of wells in the BTES system. Only used with the
  `pattern` form.
- `num_sectors = 6`: Number of sectors in the BTES system. The system is
  divided into equal circle sectors, and all wells in each sector are coupled in
  series. Only used with the `pattern` form.
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
- `topology = :sectors_parallel`: How the boreholes are connected.
  - `:sectors_parallel`: the sectors are operated in parallel, and the
    boreholes within each sector are coupled in series. Each sector receives
    `rate_charge`, so the field circulates `num_sectors*rate_charge` in total.
    Each borehole is its own pair of wells, `B<n>_supply` and `B<n>_return`.
  - `:sectors_series`: the sectors are coupled in series, and the boreholes
    within each sector are coupled in parallel off a shared wellhead. The field
    circulates `rate_charge` in total, or `rate_charge/wells_per_sector` per
    borehole on average. Each sector is a single pair of wells, `S<n>_supply`
    and `S<n>_return`, holding all its boreholes between a shared inlet and a
    shared outlet node. The split between the boreholes is resolved from
    wellbore friction, while the manifold itself is treated as resistance-free.
    Per-borehole results are reached through the `borehole` field of the well
    domain, since a sector rather than a borehole is the unit of well output.
- `reversed_discharge = false`: During charging, flow runs from the first to the
  last well in each sector, or from the first to the last sector when
  `topology = :sectors_series`. If `reversed_discharge = false`, discharge uses
  the same direction as charge. If `true`, discharge flow is reversed.
- `temperature_surface = 10 °C/283.15 K`: Temperature at the surface [K].
- `num_years = 4`: Number of years to run the simulation.
- `charge_period = ["June", "September"]`: Period during which the system is charged.
- `discharge_period = ["December", "March"]`: Period during which the system is discharged.
- `report_interval = 14 day`: Reporting interval for the simulation.
- `utes_schedule_args = NamedTuple()`: Additional arguments for the UTES schedule.
- `n_z = [3, 8, 3]`: Number of layers in the vertical direction for each layer.
- `n_xy = 3`: Number of layers in the horizontal direction for each layer.
- `mesh_args = NamedTuple()`: Additional arguments for the mesh generation.
"""
function btes(
    pattern::Symbol = :sunflower;
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
    reversed_discharge = false,
    topology = :sectors_parallel,
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

    return btes(field; well_spacing = well_spacing,depths = depths,well_layers = well_layers,topology = topology,
    density = density, thermal_conductivity = thermal_conductivity,heat_capacity = heat_capacity,
    geothermal_gradient = geothermal_gradient,temperature_charge = temperature_charge,temperature_discharge = temperature_discharge,
    rate_charge = rate_charge,rate_discharge = rate_discharge,reversed_discharge = reversed_discharge,temperature_surface = temperature_surface,num_years = num_years,charge_period = charge_period,
    discharge_period = discharge_period,report_interval = report_interval,utes_schedule_args = utes_schedule_args,n_z = n_z,n_xy = n_xy,mesh_args = mesh_args)

end

function btes(
    field::Vector{Vector{Matrix{Float64}}};
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
    reversed_discharge = false,
    topology = :sectors_parallel,
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

    topology in (:sectors_parallel, :sectors_series) ||
        error("Unknown topology: $topology. Supported topologies are \
            :sectors_parallel and :sectors_series.")

    well_coords_3d = vcat(field...)
    num_wells = length(well_coords_3d)

    # ## Create mesh
    # Use the (x,y) projection of each well trajectory as a mesh constraint
    well_coordinates = [wc[1:2, :] for wc in well_coords_3d]
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
    manifold_segments = Vector{Vector{Int}}()
    nl = length(layers)
    geo = tpfv_geometry(mesh)

    # Reservoir cells traversed by a well trajectory, top to bottom
    function trajectory_cells(wc)
        cells = Jutul.find_enclosing_cells(mesh, permutedims(wc), n = 100)
        filter!(c -> layers[c] ∈ well_layers, cells)
        return cells
    end

    if topology == :sectors_parallel
        # One pair of wells per borehole
        for (wno, wc) in enumerate(well_coords_3d)
            name = Symbol("B$wno")
            println("Adding well $name ($wno/$num_wells)")
            cells = trajectory_cells(wc)
            w_sup, w_ret = setup_btes_well(domain, cells, name=name, closed_loop_type=:u1)
            push!(well_models, w_sup, w_ret)
            push!(well_names, name)
        end
    else
        # One pair of wells per sector, with the boreholes of the sector hung in
        # parallel between a shared inlet and a shared outlet node
        for (sno, sector) in enumerate(field)
            name = Symbol("S$sno")
            println("Adding sector $name ($sno/$(length(field))) with \
                $(length(sector)) wells")
            cells = [trajectory_cells(wc) for wc in sector]
            disc = discretize_closed_loop_u1_manifold(cells, geo.cell_centroids)
            w_sup, w_ret = setup_closed_loop_well_u1(domain, disc.reservoir_cells;
                name                  = name,
                neighborship          = disc.neighborship,
                pipe_cells            = disc.pipe_cells,
                grout_cells           = disc.grout_cells,
                well_cell_centers     = disc.well_cell_centers,
                end_nodes             = disc.end_nodes,
                return_reservoir_cell = disc.return_reservoir_cell,
                tag                   = disc.tag,
                borehole              = disc.borehole)
            push!(well_models, w_sup, w_ret)
            push!(well_names, name)
            push!(manifold_segments, disc.manifold_segments)
        end
    end

    # Make the model
    model = setup_reservoir_model(
        domain, :geothermal,
        wells = well_models,
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

    if topology == :sectors_parallel
        # Group supply well names by sector, preserving the order given in `field`
        wells_per_sector = Vector{Vector{Symbol}}()
        wtot = 0
        for sector in field
            sw = [Symbol(well_names[wtot + l], "_supply") for l in 1:length(sector)]
            wtot += length(sector)
            push!(wells_per_sector, sw)
        end
    else
        # A single chain running through the sectors, in the order given in `field`
        wells_per_sector = [[Symbol(name, "_supply") for name in well_names]]
    end
    control_charge, control_discharge, sectors = setup_controls(model, wells_per_sector,
        rate_charge, rate_discharge, temperature_charge, temperature_discharge;
        reversed_discharge = reversed_discharge);
    if topology == :sectors_series
        # setup_controls saw one chain, and so returned a single group. Restore
        # one entry per sector.
        sectors = Dict(Symbol("S$sno") =>
            [Symbol(name, "_supply"), Symbol(name, "_return")]
            for (sno, name) in enumerate(well_names))
    end
    
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
    info[:topology] = topology

    # ## Set parameters
    parameters = nothing
    if topology == :sectors_series
        # The manifold is treated as resistance-free: its segment lengths, which
        # would otherwise be the distance from the wellhead to each borehole,
        # are set to zero so that node placement does not bias the flow split.
        parameters = setup_parameters(model)
        for (sno, segments) in enumerate(manifold_segments)
            parameters[Symbol("S$sno", "_supply")][:SegmentLength][segments] .= 0.0
        end
    end

    # ## Assemble and return model
    case = JutulCase(model, dt, forces,
        state0 = state0, parameters = parameters, input_data = info)
    return case

end

function btes(field::AbstractMatrix; kwargs...)
    # Special case: a single well, given as a single 3 x m matrix, representing
    # a field with a single sector containing a single well
    return btes([[Matrix{Float64}(field)]]; kwargs...)
end

# ## Patterns
#
# Each pattern function takes the number of wells and the approximate
# spacing between neighboring wells, and returns a full field

function sunflower_pattern(num_wells, spacing; num_sectors = 6, depths = [0.0, 0.5, 50, 65])
    xy = fibonacci_pattern_2d(num_wells; spacing = spacing)
    return field_from_points(:angular,xy, num_sectors, depths)
end

function rectangular_pattern(num_wells, spacing; num_sectors = 6, depths = [0.0, 0.5, 50, 65])
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
    return field_from_points(:cartesian,xy, num_sectors, depths)
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
    return field_from_points(:angular,xy, num_sectors, depths)
end

function polygonal_pattern(num_wells, spacing, num_sides; num_sectors = 6, depths = [0.0, 0.5, 50, 65])
    # Regular polygon with an area roughly matching num_wells points at the
    # given spacing
    R = sqrt(2*num_wells*spacing^2/(num_sides*sin(2π/num_sides)))
    θ = range(0.0, 2π, length = num_sides + 1)[1:num_sides]
    polygon = vcat((R*cos.(θ))', (R*sin.(θ))')

    # Build a rectangular grid covering the polygon with some margin, and
    # keep only the points that fall inside it
    margin = 1.2
    nxy = 2*ceil(Int, margin*R/spacing) + 1
    xs = ((0:nxy-1) .- (nxy-1)/2).*spacing
    xy = Matrix{Float64}(undef, 2, nxy^2)
    k = 0
    for y in xs, x in xs
        k += 1
        xy[:, k] = [x, y]
    end
    xy = xy[:, points_in_polygon(xy, polygon)]

    # Keep the num_wells points closest to the center
    r = sqrt.(xy[1,:].^2 .+ xy[2,:].^2)
    order = sortperm(r)[1:min(num_wells, size(xy, 2))]
    xy = xy[:, order]

    return field_from_points(:angular,xy, num_sectors, depths)
end

function setup_controls(model, wells_per_sector::AbstractVector{<:AbstractVector{Symbol}},
    rate_charge, rate_discharge, temperature_charge, temperature_discharge;
    reversed_discharge::Bool = false)

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
            else
                # Charging always runs from the first to the last well in the sector
                if k == 1
                    control_charge[well_sup] = ctrl_charge
                else
                    well_prev = get_return(sw[k-1])
                    target = JutulDarcy.ReinjectionTarget([well_prev])
                    control_charge[well_sup] = InjectorControl(target, [1.0],
                        density=rho, temperature=NaN; check=false)
                end
                # Discharging runs from first to last (reversed_discharge = false)
                # or from last to first (reversed_discharge = true)
                discharge_first = reversed_discharge ? length(sw) : 1
                if k == discharge_first
                    control_discharge[well_sup] = ctrl_discharge
                else
                    well_prev = get_return(sw[reversed_discharge ? k+1 : k-1])
                    target = JutulDarcy.ReinjectionTarget([well_prev])
                    control_discharge[well_sup] = InjectorControl(target, [1.0],
                        density=rho, temperature=NaN; check=false)
                end
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
