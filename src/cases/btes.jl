"""
    btes(pattern::Symbol = :sunflower; <keyword arguments>)
    btes(field; <keyword arguments>)

Setup function for borehole thermal energy storage (BTES) system.

The first form places the wells using one of the built-in patterns. The second
form takes an explicit well field, and does not accept the pattern arguments
(`num_wells`, `radius`, `extent`, `num_sectors`, `sector_grouping` and
`sector_dims`), since the placement and grouping are given by `field` itself.

# Arguments
- `pattern = :sunflower`: Well placement pattern.
  - `:sunflower`: Fibonacci spiral (see [`fibonacci_pattern_2d`](@ref)), with
    the outermost well at distance `radius` from the center. The well spacing
    is the median distance between neighboring wells.
  - `:circular`: a center well surrounded by concentric rings, with the
    outermost ring at distance `radius` from the center. The wells are
    numbered from the center outwards, counter-clockwise on each ring starting
    from the positive x axis.
  - `:rectangular`: rows of wells, with the outermost wells spanning at most
    `extent = (Lx, Ly)`. The last row is partly filled if needed. Without
    `extent`, the field is a whole number of sector blocks wide (see
    `sector_dims`), and close to square.
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
- `num_wells`, `well_spacing`: Number of wells, and horizontal spacing between
  neighboring wells [m]. Give at most one of them: the other follows from the
  field size (`radius` or `extent`). Without a field size, the missing one
  defaults to `num_wells = 48` or `well_spacing = 5.0`, and the field size
  follows. With the `field` form, `well_spacing = 5.0` only sets the horizontal
  mesh resolution.
- `radius`: Field size for `:sunflower` and `:circular` [m]. Only used with the
  `pattern` form.
- `extent`: Field size `(Lx, Ly)` for `:rectangular` [m]. Only used with the
  `pattern` form.
- `num_sectors = 6`: Number of sectors in the BTES system. Only used with the
  `pattern` form.
- `sector_grouping`: How the wells are grouped into sectors. Only used with the
  `pattern` form.
  - `:spiral` (default for `:sunflower`): bundles of adjacent spiral spines of
    the sunflower, with the wells of each sector ordered from the center
    outwards: the first well of each spine, counter-clockwise, then the second
    well of each spine, and so on (see `Fimbul.group_into_sectors_spiral`).
    The spines come in families of Fibonacci size `F` (2, 3, 5, 8, 13, 21,
    ...), and only families with at most `num_wells/2` spines are used. Whole
    spines divide evenly into sectors when `num_wells` is a multiple of
    `num_sectors` and `num_sectors` divides such a family: 2, 3, 5 or 8
    sectors (one spine per sector) for any multiple, 4 sectors (two spines of
    the 8-family) from 16 wells, and 7 sectors (three spines of the 21-family)
    from 42 wells. Otherwise, as for 6 sectors, the outermost wells on the
    boundary spines are moved to neighboring sectors, so that the sector sizes
    differ by at most one.
  - `:angular` (default for `:circular`): circle sectors, with the wells of
    each sector ordered from the center outwards.
  - `:sequential` (`:sunflower` and `:circular`): consecutive groups of wells
    in the order they are placed, from the center outwards, so that the first
    sector holds the wells closest to the center.
  - `:cartesian` (default and only option for `:rectangular`): blocks of
    `sector_dims` wells tiling the full rows of the field from the bottom up,
    and bands across the field for the remaining wells at the top (see
    `Fimbul.group_into_sectors_cartesian`). Within each sector, the wells are
    ordered column by column. For example, 45 wells in 5 sectors gives a
    6 x 8 field with four blocks of 3 x 3 wells, and a last sector of 9 wells
    in the two top rows.
- `sector_dims`: Dimensions `(nx, ny)` of the blocks of the `:rectangular`
  pattern, in number of wells along x and y. Defaults to the factor pair of
  `num_wells/num_sectors` closest to square. Only used with the `pattern`
  form.
- `depths = [0.0, 0.5, 50, 65]`: Depths delineating geological layers [m].
- `well_layers = [2]`: Layers in which the wells are perforated. Layer 1 is the
  insulation cover; the pipes pass through it but exchange nothing with it, so
  by default the boreholes start below it.
- `density = [30, 2580, 2580]: Rock density in the layers [kg/m³].
- `thermal_conductivity = [0.034, 3.7, 3.7]: Thermal conductivity in the layers [W/(m⋅K)].
- `heat_capacity = [1500, 900, 900]`: Heat capacity in the layers [J/(kg⋅K)].
- `permeability = 1e-6 darcy`: Rock permeability [m²].
- `geothermal_gradient = 0.03 K/m`: Geothermal gradient [K/m].
- `temperature_charge = 90 °C/363.15 K`: Injection temperature during charging [K].
- `temperature_discharge = 10 °C/283.15 K`: Injection temperature during discharging [K].
- `rate_charge = 0.5 l/s`: Injection rate during charging [m³/s].
- `rate_discharge = rate_charge`: Injection rate during discharging [m³/s].
- `topology = :sectors_parallel`: How the boreholes are connected.
  - `:sectors_parallel`: the sectors are operated in parallel, and the
    boreholes within each sector are coupled in series. Each sector receives
    `rate_charge`, so the field circulates `num_sectors*rate_charge` in total.
    How the boreholes are represented is set by `well_grouping`.
  - `:sectors_series`: the sectors are coupled in series, and the boreholes
    within each sector are coupled in parallel off a shared wellhead. The field
    circulates `rate_charge` in total, or `rate_charge/wells_per_sector` per
    borehole on average. Each sector is a single closed-loop well `S<n>`,
    controlled as `S<n>_supply` and `S<n>_return`, holding all its boreholes
    between a shared inlet and a shared outlet node. The split between the
    boreholes is resolved from wellbore friction, while the manifold itself is
    treated as resistance-free.
- `well_grouping = :sector`: Wells representing the boreholes when
  `topology = :sectors_parallel`.
  - `:sector`: each sector is a single closed-loop well `S<n>`, controlled as
    `S<n>_supply` and `S<n>_return`, with the boreholes connected in series
    within the well. The surface pipes between the boreholes are treated as
    resistance-free.
  - `:borehole`: each borehole is its own closed-loop well `B<n>`, controlled
    as `B<n>_supply` and `B<n>_return` (see [`closed_loop_well_names`](@ref)).
    Each borehole is fed by the return of the previous one with a
    `ReinjectionTarget`. Use this to control individual boreholes, e.g. to
    operate only some of them.
  With a single well per sector, per-borehole results are reached through the
  `borehole` field of the well domain.
- `reversed_discharge = false`: During charging, flow runs from the first to the
  last well in each sector, or from the first to the last sector when
  `topology = :sectors_series`. If `reversed_discharge = false`, discharge uses
  the same direction as charge. If `true`, discharge flow is reversed. A sector
  modelled as a single well (`well_grouping = :sector`) is then fed at its
  return part, so that the fluid passes the boreholes, and each U-tube, in
  reverse.
- `temperature_surface = 10 °C/283.15 K`: Temperature at the surface [K].
- `num_years = 4`: Number of years to run the simulation.
- `charge_period = ["June", "September"]`: Period during which the system is charged.
- `discharge_period = ["December", "March"]`: Period during which the system is discharged.
- `report_interval = 14 day`: Reporting interval for the simulation.
  During rest periods the supply parts of the closed-loop wells are put under
  [`ClosedLoopRestControl`](@ref) at the wellhead pressure, and the return parts
  are disabled.
- `utes_schedule_args = NamedTuple()`: Additional arguments for the UTES schedule.
- `n_z = [3, 8, 3]`: Number of layers in the vertical direction for each layer.
- `n_xy = 3`: Number of layers in the horizontal direction for each layer.
- `mesh_args = NamedTuple()`: Additional arguments for the mesh generation.
"""
function btes(
        pattern::Symbol = :sunflower;
        num_wells = missing,
        well_spacing = missing,
        radius = missing,
        extent = missing,
        sector_grouping = missing,
        sector_dims = missing,
        num_sectors = 6,
        depths = [0.0, 0.5, 50, 65],
        well_layers = [2],
        density = [30, 2580, 2580] * kilogram / meter^3,
        thermal_conductivity = [0.034, 3.7, 3.7] * watt / meter / Kelvin,
        heat_capacity = [1500, 900, 900] * joule / kilogram / Kelvin,
        permeability = 1.0e-6darcy,
        geothermal_gradient = 0.03Kelvin / meter,
        temperature_charge = convert_to_si(90.0, :Celsius),
        temperature_discharge = convert_to_si(10.0, :Celsius),
        rate_charge = 0.5litre / second,
        rate_discharge = rate_charge,
        reversed_discharge = false,
        topology = :sectors_parallel,
        well_grouping = :sector,
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

    field, well_spacing = btes_well_field(
        pattern;
        num_wells = num_wells, well_spacing = well_spacing,
        radius = radius, extent = extent, sector_grouping = sector_grouping,
        sector_dims = sector_dims, num_sectors = num_sectors, depths = depths
    )

    case = btes(
        field; well_spacing = well_spacing, depths = depths, well_layers = well_layers, topology = topology, well_grouping = well_grouping,
        density = density, thermal_conductivity = thermal_conductivity, heat_capacity = heat_capacity, permeability = permeability,
        geothermal_gradient = geothermal_gradient, temperature_charge = temperature_charge, temperature_discharge = temperature_discharge,
        rate_charge = rate_charge, rate_discharge = rate_discharge, reversed_discharge = reversed_discharge, temperature_surface = temperature_surface, num_years = num_years, charge_period = charge_period,
        discharge_period = discharge_period, report_interval = report_interval, utes_schedule_args = utes_schedule_args, n_z = n_z, n_xy = n_xy, mesh_args = mesh_args
    )
    case.input_data[:pattern] = pattern
    return case

end

function btes(
        field::Vector{Vector{Matrix{Float64}}};
        well_spacing = 5.0,
        depths = [0.0, 0.5, 50, 65],
        well_layers = [2],
        density = [30, 2580, 2580] * kilogram / meter^3,
        thermal_conductivity = [0.034, 3.7, 3.7] * watt / meter / Kelvin,
        heat_capacity = [1500, 900, 900] * joule / kilogram / Kelvin,
        permeability = 1.0e-6darcy,
        geothermal_gradient = 0.03Kelvin / meter,
        temperature_charge = convert_to_si(90.0, :Celsius),
        temperature_discharge = convert_to_si(10.0, :Celsius),
        rate_charge = 0.5litre / second,
        rate_discharge = rate_charge,
        reversed_discharge = false,
        topology = :sectors_parallel,
        well_grouping = :sector,
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
    well_grouping in (:sector, :borehole) ||
        error("Unknown well_grouping: $well_grouping. Supported values are \
            :sector and :borehole.")
    if topology == :sectors_series && well_grouping == :borehole
        error("well_grouping = :borehole is only supported with topology = \
            :sectors_parallel. With :sectors_series, the boreholes of a sector \
            share a single well.")
    end

    well_coords_3d = vcat(field...)
    num_wells = length(well_coords_3d)

    # ## Create mesh
    # Use the (x,y) projection of each well trajectory as a mesh constraint.
    # Repeated points are removed, so that a vertical well is a point
    # constraint: a curve of zero length would give a zero distance between
    # constraints, and with it no refinement around the wells.
    well_coordinates = [unique_consecutive_columns(wc[1:2, :]) for wc in well_coords_3d]
    hz = diff(depths) ./ n_z
    hxy = well_spacing / n_xy

    # ## Set up model
    # Set up reservoir domain with rock properties similar to that of granite,
    # with a styrofoam layer on top
    domain, layers, metrics = layered_reservoir_domain(
        well_coordinates, depths,
        (
            rock_density = density,
            rock_thermal_conductivity = thermal_conductivity,
            rock_heat_capacity = heat_capacity,
        );
        mesh_args = (; hxy_min = hxy, hz = hz, mesh_args...),
        permeability = permeability,
        porosity = 0.01,
        component_heat_capacity = 4.278e3joule / kilogram / Kelvin,
    )
    mesh = physical_representation(domain)
    # Set up BTES wells
    hxy_min = metrics.hxy_min
    well_models = []
    well_names = Symbol[]
    surface_segments_per_well = Vector{Vector{Int}}()
    nl = length(layers)
    geo = tpfv_geometry(mesh)

    # Reservoir cells traversed by a well trajectory, top to bottom. The well
    # positions are mesh nodes, so the trajectory is shifted half a cell away
    # from the field center to run through a cell rather than along its edges.
    function trajectory_cells(wc)
        xy = wc[1:2, 1]
        d = norm(xy)
        v = d > 0 ? xy ./ d : [1.0, 0.0]
        wc = copy(wc)
        wc[1:2, :] .+= (hxy_min / 2) .* v
        cells = Jutul.find_enclosing_cells(mesh, permutedims(wc), geometry = geo)
        filter!(c -> layers[c] ∈ well_layers, cells)
        return cells
    end

    if topology == :sectors_parallel && well_grouping == :borehole
        # One closed-loop well per borehole
        for (wno, wc) in enumerate(well_coords_3d)
            name = Symbol("B$wno")
            println("Adding well $name ($wno/$num_wells)")
            cells = trajectory_cells(wc)
            well = setup_btes_well(domain, cells, name = name, closed_loop_type = :u1)
            push!(well_models, well)
            push!(well_names, name)
        end
    else
        # One closed-loop well per sector, with the boreholes of the sector in
        # series, or hung in parallel between a shared inlet and outlet node
        for (sno, sector) in enumerate(field)
            name = Symbol("S$sno")
            println("Adding sector $name ($sno/$(length(field))) with \
                $(length(sector)) wells")
            cells = [trajectory_cells(wc) for wc in sector]
            if topology == :sectors_parallel
                disc = discretize_closed_loop_u1_series(cells, geo.cell_centroids)
                surface_segments = disc.connector_segments
            else
                disc = discretize_closed_loop_u1_manifold(cells, geo.cell_centroids)
                surface_segments = disc.manifold_segments
            end
            well = setup_btes_well(
                domain, disc.reservoir_cells;
                name = name,
                closed_loop_type = :u1,
                neighborship = disc.neighborship,
                pipe_cells = disc.pipe_cells,
                grout_cells = disc.grout_cells,
                well_cell_centers = disc.well_cell_centers,
                end_nodes = disc.end_nodes,
                return_reservoir_cell = disc.return_reservoir_cell,
                tag = disc.tag,
                borehole = disc.borehole
            )
            push!(well_models, well)
            push!(well_names, name)
            push!(surface_segments_per_well, surface_segments)
        end
    end

    # Make the model. The facility must not shut wells on its own: the supply
    # and return parts of a closed loop are hydraulically connected, so shutting
    # one of them gives a dead end, and the chained reinjection controls start
    # every period from zero rate, which is exactly when the facility would shut
    # them.
    model = setup_reservoir_model(
        domain, :geothermal,
        wells = well_models,
        can_shut_wells = false,
    )

    # ## Set up initial state and boundary conditions
    geo = tpfv_geometry(mesh)
    z_bc = geo.boundary_centroids[3, :]
    bottom = map(v -> isapprox(v, maximum(z_bc)), z_bc)
    # Define pressure and temperature profiles
    rho = reservoir_model(model).system.rho_ref[1]
    dpdz = rho * gravity_constant
    dTdz = geothermal_gradient
    T = z -> temperature_surface .+ dTdz * z
    p = z -> 5atm .+ dpdz .* z
    # Set initial conditions
    z_cells = geo.cell_centroids[3, :]
    z_hat = z_cells .- minimum(z_cells)
    state0 = setup_reservoir_state(
        model,
        Pressure = p(z_hat),
        Temperature = T(z_hat)
    )
    # Set boundary conditions
    z_bc = z_bc[.!bottom]
    z_hat = z_bc .- minimum(z_bc)
    bc_cells = geo.boundary_neighbors[.!bottom]
    bc = flow_boundary_condition(bc_cells, domain, p(z_hat), T(z_hat))

    if topology == :sectors_parallel && well_grouping == :borehole
        # Group closed-loop wells by sector, preserving the order given in `field`
        wells_per_sector = Vector{Vector{Symbol}}()
        wtot = 0
        for sector in field
            push!(wells_per_sector, well_names[wtot .+ (1:length(sector))])
            wtot += length(sector)
        end
    elseif topology == :sectors_parallel
        # One well per sector, each fed on its own
        wells_per_sector = [[name] for name in well_names]
    else
        # A single chain running through the sectors, in the order given in `field`
        wells_per_sector = [copy(well_names)]
    end
    control_charge, control_discharge, control_rest, sectors = setup_controls(
        model, wells_per_sector,
        rate_charge, rate_discharge, temperature_charge, temperature_discharge;
        reversed_discharge = reversed_discharge,
        reverse_single_wells = topology == :sectors_parallel && well_grouping == :sector
    )
    if topology == :sectors_series || well_grouping == :sector
        # One entry per sector, each holding the closed-loop well of the sector
        sectors = Dict(
            Symbol("S$sno") => [name]
                for (sno, name) in enumerate(well_names)
        )
    end

    forces_charge = setup_reservoir_forces(model, control = control_charge, bc = bc)
    forces_discharge = setup_reservoir_forces(model, control = control_discharge, bc = bc)
    forces_rest = setup_reservoir_forces(model, control = control_rest, bc = bc)
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
    info[:field] = field
    info[:well_spacing] = well_spacing

    # ## Set parameters
    parameters = nothing
    if !isempty(surface_segments_per_well)
        # The surface pipes (the manifold, or the connections between boreholes
        # in series) are treated as resistance-free: their segment lengths,
        # which would otherwise follow from the node positions, are set to zero
        # so that node placement does not affect the flow.
        parameters = setup_parameters(model)
        for (name, segments) in zip(well_names, surface_segments_per_well)
            parameters[name][:SegmentLength][segments] .= 0.0
        end
    end

    # ## Assemble and return model
    case = JutulCase(
        model, dt, forces,
        state0 = state0, parameters = parameters, input_data = info
    )
    return case

end

function btes(field::AbstractMatrix; kwargs...)
    # Special case: a single well, given as a single 3 x m matrix, representing
    # a field with a single sector containing a single well
    return btes([[Matrix{Float64}(field)]]; kwargs...)
end

# Columns of x without consecutive repetitions
function unique_consecutive_columns(x::AbstractMatrix; tol = 1.0e-10)
    keep = [1]
    for i in 2:size(x, 2)
        norm(x[:, i] - x[:, keep[end]]) > tol && push!(keep, i)
    end
    return x[:, keep]
end

# ## Patterns
#
# Each pattern is set by the number of wells `N`, the spacing `s` between
# neighboring wells and the field size, of which exactly two are given. The
# third is computed, and the full field is returned along with the spacing.

"""
    btes_well_field(pattern = :sunflower; num_wells = missing, well_spacing = missing,
        radius = missing, extent = missing, sector_grouping = missing,
        sector_dims = missing, num_sectors = 6, depths = [0.0, 0.5, 50, 65])

Well field for [`btes`](@ref) with a built-in pattern, returned as
`(field, well_spacing)`. See [`btes`](@ref) for the patterns and how the number
of wells, the spacing and the field size are set.
"""
function btes_well_field(
        pattern::Symbol = :sunflower;
        num_wells = missing,
        well_spacing = missing,
        radius = missing,
        extent = missing,
        sector_grouping = missing,
        sector_dims = missing,
        num_sectors = 6,
        depths = [0.0, 0.5, 50, 65],
    )
    ismissing(num_wells) || ismissing(well_spacing) ||
        error("Give either num_wells or well_spacing, not both. The other \
            follows from the field size (radius or extent).")
    if pattern in (:sunflower, :circular)
        ismissing(extent) || error("extent is only used with the :rectangular \
            pattern. Use radius to set the size of a :$pattern field.")
        field_size = radius
    elseif pattern == :rectangular
        ismissing(radius) || error("radius is not used with the :rectangular \
            pattern. Use extent = (Lx, Ly) to set the size of the field.")
        field_size = extent
    else
        error("Unknown pattern: $pattern. Supported patterns are :sunflower, \
            :rectangular and :circular.")
    end
    if ismissing(field_size)
        ismissing(num_wells) && (num_wells = 48)
        ismissing(well_spacing) && (well_spacing = 5.0)
    elseif ismissing(num_wells) && ismissing(well_spacing)
        well_spacing = 5.0
    end
    ismissing(num_wells) || num_wells >= 1 ||
        error("num_wells must be at least 1, got $num_wells.")
    ismissing(well_spacing) || well_spacing > 0 ||
        error("well_spacing must be positive, got $well_spacing.")

    allowed = Dict(
        :sunflower => (:spiral, :sequential),
        :circular => (:angular, :sequential), :rectangular => (:cartesian,)
    )[pattern]
    if ismissing(sector_grouping)
        sector_grouping = first(allowed)
    end
    sector_grouping in allowed || error("sector_grouping = :$sector_grouping \
        is not supported for the :$pattern pattern. Supported values are \
        $(join(":" .* string.(allowed), ", ")).")

    # The wells of the sunflower and circular patterns are placed from the
    # center outwards (counter-clockwise on each ring for :circular), which is
    # the order :sequential groups them in
    if pattern == :sunflower
        xy, well_spacing = sunflower_pattern(num_wells, well_spacing, field_size)
    elseif pattern == :circular
        xy, well_spacing = circular_pattern(num_wells, well_spacing, field_size)
    else
        xy, well_spacing = rectangular_pattern(
            num_wells, well_spacing, field_size;
            num_sectors = num_sectors, sector_dims = sector_dims
        )
    end
    division = sector_grouping
    ismissing(sector_dims) || pattern == :rectangular ||
        error("sector_dims is only used with the :rectangular pattern.")
    size(xy, 2) >= 1 || error("No wells fit in a :$pattern field of size \
        $field_size with spacing $well_spacing.")
    field = field_from_points(
        division, xy, min(num_sectors, size(xy, 2)), depths;
        sector_dims = sector_dims
    )
    return (field, Float64(well_spacing))
end

# Outermost radius of a sunflower pattern with N wells at spacing s. The
# distance between the innermost wells, which sets the spacing, scales as
# 1/sqrt(N) for a fixed outer radius.
function sunflower_radius(N, s)
    N == 1 && return 0.0
    xy = fibonacci_pattern_2d(N; spacing = s)
    return maximum(sqrt.(xy[1, :] .^ 2 .+ xy[2, :] .^ 2))
end

function sunflower_pattern(N, s, R)
    if ismissing(N)
        # Largest number of wells that fits within R. The radius grows about
        # as sqrt(N), which gives a first guess to correct from.
        fits(n) = sunflower_radius(n, s) <= R * (1 + 1.0e-10)
        N = max(1, floor(Int, (R / sunflower_radius(2, s))^2 * 2))
        while N > 1 && !fits(N)
            N -= 1
        end
        while fits(N + 1)
            N += 1
        end
    elseif ismissing(s)
        N == 1 && return (zeros(2, 1), Float64(R))
        s = R / sunflower_radius(N, 1.0)
    end
    N == 1 && return (zeros(2, 1), Float64(s))
    return (fibonacci_pattern_2d(N; spacing = s), Float64(s))
end

# Number of wells on ring k of the circular pattern, where the ring radius is k
# spacings. The count does not depend on the spacing itself.
circular_ring_size(k) = max(1, round(Int, 2π * k))

function circular_pattern(N, s, R)
    if ismissing(N)
        K = floor(Int, R / s * (1 + 1.0e-10))
        N = 1 + sum(circular_ring_size, 1:K; init = 0)
    elseif ismissing(s)
        # The number of rings is fixed by N, so the outer ring sets the spacing
        K, n = 0, 1
        while n < N
            K += 1
            n += circular_ring_size(K)
        end
        s = K == 0 ? Float64(R) : R / K
    end
    points = [[0.0, 0.0]]
    ring = 1
    while length(points) < N
        r = ring * s
        n_ring = min(circular_ring_size(ring), N - length(points))
        for k in 0:(n_ring - 1)
            θ = 2π * k / n_ring
            push!(points, [r * cos(θ), r * sin(θ)])
        end
        ring += 1
    end
    return (hcat(points...), Float64(s))
end

function rectangular_pattern(N, s, extent; num_sectors = 1, sector_dims = missing)
    if !ismissing(extent)
        length(extent) == 2 && all(>=(0), extent) ||
            error("extent must be (Lx, Ly) with Lx, Ly ≥ 0, got $extent.")
        Lx, Ly = extent
        counts(s) = (floor(Int, Lx / s * (1 + 1.0e-10)) + 1, floor(Int, Ly / s * (1 + 1.0e-10)) + 1)
        if ismissing(N)
            nx, ny = counts(s)
            N = nx * ny
        else
            # Largest spacing that fits N wells. The count only changes where
            # a whole number of spacings fits Lx or Ly, so those are the
            # candidates.
            candidates = filter(>(0), vcat(Lx ./ (1:N), Ly ./ (1:N)))
            candidates = filter(c -> prod(counts(c)) >= N, candidates)
            if isempty(candidates)
                N == 1 || error("extent $extent has no room for $N wells.")
                s = 1.0
            else
                s = maximum(candidates)
            end
            nx, ny = counts(s)
        end
    else
        # A whole number of sector blocks wide, so that the full rows can be
        # tiled with sectors
        num_sectors = clamp(num_sectors, 1, N)
        dims = rectangular_sector_dims(N ÷ num_sectors, sector_dims)
        nx = rectangular_grid_width(N, num_sectors, dims)
        ny = ceil(Int, N / nx)
    end
    xs = ((0:(nx - 1)) .- (nx - 1) / 2) .* s
    ys = ((0:(ny - 1)) .- (ny - 1) / 2) .* s
    xy = Matrix{Float64}(undef, 2, nx * ny)
    k = 0
    for y in ys, x in xs
        k += 1
        xy[:, k] = [x, y]
    end
    return (xy[:, 1:N], Float64(s))
end

# Controls for the closed-loop wells in `wells_per_sector`, where the wells of
# each sector are coupled in series in the given order
function setup_controls(
        model, wells_per_sector::AbstractVector{<:AbstractVector{Symbol}},
        rate_charge, rate_discharge, temperature_charge, temperature_discharge;
        reversed_discharge::Bool = false,
        reverse_single_wells::Bool = false
    )

    rho = reservoir_model(model).system.rho_ref[1]
    rate_target = TotalRateTarget(rate_charge)
    ctrl_charge = InjectorControl(
        rate_target, [1.0],
        density = rho, temperature = temperature_charge
    )
    rate_target = TotalRateTarget(rate_discharge)
    ctrl_discharge = InjectorControl(
        rate_target, [1.0],
        density = rho, temperature = temperature_discharge
    )
    # BHP control for return side
    bhp_target = BottomHolePressureTarget(1.0si_unit(:atm))
    ctrl_ret = ProducerControl(bhp_target)
    # Rest: the loops are not circulating. Holding the wellhead pressure lets the
    # fluid in each loop expand and contract freely with temperature, as an
    # expansion vessel would. Disabling every well instead seals the loops, and
    # a cooling loop then cannot converge (see ClosedLoopRestControl).
    ctrl_rest = ClosedLoopRestControl(bhp_target)
    # Set up forces
    control_charge = Dict()
    control_discharge = Dict()
    control_rest = Dict()
    assigned = []
    # Supply and return parts of a closed-loop well, controlled separately
    get_supply = (well) -> first(closed_loop_well_names(model.models[well].data_domain))
    get_return = (well) -> last(closed_loop_well_names(model.models[well].data_domain))
    sectors = Dict()

    for (sno, sw) in enumerate(wells_per_sector)
        sec_wells = Symbol[]
        for (k, well) in enumerate(sw)
            well_sup = get_supply(well)
            well_ret = get_return(well)
            @assert well_sup ∉ assigned
            @assert well_ret ∉ assigned
            control_charge[well_ret] = ctrl_ret
            control_discharge[well_ret] = ctrl_ret
            if length(sw) == 1
                # Single well in sector: charge and discharge it directly
                control_charge[well_sup] = ctrl_charge
                if reversed_discharge && reverse_single_wells
                    # Reverse the flow through the well: inject into the return
                    # part and produce from the supply part
                    control_discharge[well_ret] = ctrl_discharge
                    control_discharge[well_sup] = ctrl_ret
                else
                    control_discharge[well_sup] = ctrl_discharge
                end
            else
                # Charging always runs from the first to the last well in the sector
                if k == 1
                    control_charge[well_sup] = ctrl_charge
                else
                    well_prev = get_return(sw[k - 1])
                    target = JutulDarcy.ReinjectionTarget([well_prev])
                    control_charge[well_sup] = InjectorControl(
                        target, [1.0],
                        density = rho, temperature = NaN; check = false
                    )
                end
                # Discharging runs from first to last (reversed_discharge = false)
                # or from last to first (reversed_discharge = true)
                discharge_first = reversed_discharge ? length(sw) : 1
                if k == discharge_first
                    control_discharge[well_sup] = ctrl_discharge
                else
                    well_prev = get_return(sw[reversed_discharge ? k + 1 : k - 1])
                    target = JutulDarcy.ReinjectionTarget([well_prev])
                    control_discharge[well_sup] = InjectorControl(
                        target, [1.0],
                        density = rho, temperature = NaN; check = false
                    )
                end
            end
            control_rest[well_sup] = ctrl_rest
            control_rest[well_ret] = DisabledControl()
            push!(assigned, well_sup, well_ret)
            push!(sec_wells, well)
        end
        sectors[Symbol("S$sno")] = sec_wells
    end

    @assert sort(assigned) == sort(well_symbols(model))

    return control_charge, control_discharge, control_rest, sectors

end
