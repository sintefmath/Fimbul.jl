"""
    coaxial_well_branches(; kwargs...)

Create a multi-branch coaxial geothermal well system simulation case.

The system consists of a single shared vertical wellbore (the "trunk") from
which multiple branches diverge in a conical/fan-out geometry at depth. Each
branch is a coaxial closed-loop well that shares a common top node (wellhead).
Branch collars are arranged using a Fibonacci spiral pattern for optimal
spacing.

# Keyword Arguments

## Branch geometry
- `n_branches = 4`: Number of lateral branches.
- `branch_surface_spacing = 200.0`: Minimum inter-branch spacing at surface
  level [m]. Used to set the Fibonacci pattern radius.
- `trunk_depth = 2000.0`: Depth at which the trunk ends and branches begin [m].
- `branch_length = 1500.0`: Length of each branch (measured along trajectory) [m].
- `branch_dz = 400.0`: Additional vertical depth gained by each branch [m].
- `deviation_depth = 1500.0`: Depth at which branches start to deviate from the
  trunk [m].
- `trajectory_points = 20`: Number of points used to discretize each branch
  trajectory.

## Geological parameters
- `depths = [0.0, 500.0, 1500.0, 2000.0, 2500.0, 3000.0]`: Depth boundaries
  for geological layers [m].
- `permeability = [1e-3, 1e-3, 1e-3, 1e-2, 1e-3]*darcy`: Permeability per
  layer [darcy].
- `porosity = [0.01, 0.01, 0.01, 0.05, 0.01]`: Porosity per layer [-].
- `rock_thermal_conductivity`: Thermal conductivity per layer [W/(m·K)].
- `rock_heat_capacity`: Rock heat capacity per layer [J/(kg·K)].
- `rock_density`: Rock density per layer [kg/m³].

## Thermal parameters
- `temperature_surface = 10°C`: Surface temperature [K].
- `thermal_gradient = 0.03`: Geothermal gradient [K/m].

## Operational parameters
- `rate = 50 m³/h`: Total flow rate through the system [m³/s].
- `temperature_inj = 25°C`: Injection temperature [K].
- `num_years = 30`: Total simulation time [years].
- `report_interval = year/4`: Reporting interval [s].
- `schedule_args`: Additional keyword arguments passed to `make_schedule`.

## Mesh parameters
- `hz`: Vertical cell sizes per layer [m].
- `hxy_min = 25.0`: Minimum horizontal cell size [m].
- `hxy_max = 500.0`: Maximum horizontal cell size [m].
- `mesh_args`: Additional keyword arguments passed to `extruded_mesh`.

# Returns
A `JutulCase` for multi-branch coaxial geothermal simulation.
"""
function coaxial_well_branches(;
    # Branch geometry
    n_branches = 4,
    branch_surface_spacing = 200.0,
    trunk_depth = 2000.0,
    branch_length = 1500.0,
    branch_dz = 400.0,
    deviation_depth = 1500.0,
    trajectory_points = 20,
    # Geological parameters
    depths = [0.0, 500.0, 1500.0, 2000.0, 2500.0, 3000.0],
    permeability = [1e-3, 1e-3, 1e-3, 1e-2, 1e-3]*darcy,
    porosity = [0.01, 0.01, 0.01, 0.05, 0.01],
    rock_thermal_conductivity = [2.5, 2.5, 2.8, 3.5, 2.5]*watt/(meter*Kelvin),
    rock_heat_capacity = [900, 900, 900, 900, 900]*joule/(kilogram*Kelvin),
    rock_density = [2600, 2600, 2600, 2600, 2600]*kilogram/meter^3,
    # Thermal parameters
    temperature_surface = convert_to_si(10.0, :Celsius),
    thermal_gradient = 0.03Kelvin/meter,
    # Operational parameters
    rate = 50meter^3/hour,
    temperature_inj = convert_to_si(25.0, :Celsius),
    num_years = 30,
    report_interval = year/4,
    schedule_args = NamedTuple(),
    # Mesh parameters
    hz = missing,
    hxy_min = 25.0,
    hxy_max = 500.0,
    mesh_args = NamedTuple()
)
    @assert n_branches >= 2 "At least 2 branches are required"
    @assert deviation_depth < trunk_depth "deviation_depth must be less than trunk_depth"
    @assert trunk_depth > 0 "trunk_depth must be positive"
    @assert branch_length > 0 "branch_length must be positive"

    # ## Generate well trajectories
    well_coords, well_connectivity = generate_branch_trajectories(
        n_branches,
        branch_surface_spacing,
        trunk_depth,
        branch_length,
        branch_dz,
        deviation_depth,
        trajectory_points
    )

    # ## Set up vertical mesh sizing
    if ismissing(hz)
        num_layers = length(depths) - 1
        dz = diff(depths)
        hz = fill(250.0, num_layers)
        # Use finer resolution in the branch target layer
        for (i, d) in enumerate(dz)
            hz[i] = min(hz[i], d / 5)
        end
    end

    # ## Create mesh constraints from well trajectories
    constraints = get_branch_constraints(well_coords; hxy_min = hxy_min)

    # ## Create domain
    domain, layers, metrics = layered_reservoir_domain(constraints, depths;
        mesh_args = (;
            hz = hz,
            hxy_min = hxy_min,
            hxy_max = hxy_max,
            offset_rel = 1.0,
            dist_min_factor = 50.0,
            mesh_args...
        ),
        layer_properties = (
            permeability = permeability,
            porosity = porosity,
            rock_thermal_conductivity = rock_thermal_conductivity,
            rock_heat_capacity = rock_heat_capacity,
            rock_density = rock_density,
        )
    )

    # ## Set up wells using the AGS well setup pattern
    wells, section_info = setup_branch_wells(domain, well_coords, well_connectivity)

    model = setup_reservoir_model(
        domain, :geothermal; wells = wells)
    bc, state0 = set_dirichlet_bcs(model;
        pressure_surface = 10atm,
        temperature_surface = temperature_surface,
        geothermal_gradient = thermal_gradient,
    )

    # ## Set up controls
    rho = reservoir_model(model).system.rho_ref[1]
    # Injector: inject cooled water at specified rate
    rate_target = TotalRateTarget(rate)
    ctrl_supply = InjectorControl(
        rate_target, [1.0], density = rho, temperature = temperature_inj)
    # Producer: BHP-controlled
    bhp_target = BottomHolePressureTarget(10atm)
    ctrl_return = ProducerControl(bhp_target)
    control = Dict(
        :CoaxialBranch_supply => ctrl_supply,
        :CoaxialBranch_return => ctrl_return,
    )

    forces = setup_reservoir_forces(model, control = control, bc = bc)

    dt, forces = make_schedule(
        [forces], [(1,1), (1,1)];
        num_years = num_years,
        report_interval = report_interval,
        schedule_args...)

    # ## Additional case info
    info = Dict(
        :description => "Multi-branch coaxial geothermal well system",
        :well_coords => well_coords,
        :sections => section_info,
        :n_branches => n_branches,
    )

    # ## Create and return the complete simulation case
    case = JutulCase(model, dt, forces; state0 = state0, input_data = info)

    return case

end

"""
    generate_branch_trajectories(n_branches, surface_spacing, trunk_depth,
        branch_length, branch_dz, deviation_depth, trajectory_points)

Generate well trajectories for a multi-branch coaxial system.

All branches share a common vertical trunk from the surface. Below
`deviation_depth`, each branch smoothly curves outward in a conical fan-out
geometry. Branch collar positions at the bottom of the trunk are placed using
`fibonacci_pattern_2d` for a space-filling, low-interference arrangement.

# Returns
- `well_coords`: Vector of N×3 matrices (one per well section: trunk +
  branches + producer).
- `well_connectivity`: Matrix specifying section-to-section connectivity.
"""
function generate_branch_trajectories(
    n_branches,
    surface_spacing,
    trunk_depth,
    branch_length,
    branch_dz,
    deviation_depth,
    trajectory_points
)
    # Place branch endpoints using Fibonacci pattern
    collar_xy = fibonacci_pattern_2d(n_branches;
        spacing = surface_spacing)

    # Generate trunk trajectory (vertical from surface to deviation depth)
    trunk_step = 25.0
    n_trunk = max(Int(ceil(trunk_depth / trunk_step)), 2)
    trunk_z = range(0.0, trunk_depth, length = n_trunk + 1)
    trunk = zeros(length(trunk_z), 3)
    trunk[:, 3] .= collect(trunk_z)

    # Generate branch trajectories
    branches = Vector{Matrix{Float64}}()
    for b in 1:n_branches
        # Target endpoint of this branch
        x_end = collar_xy[1, b]
        y_end = collar_xy[2, b]
        z_end = trunk_depth + branch_dz

        # The branch starts at deviation_depth with (x,y) = (0,0) and smoothly
        # curves to (x_end, y_end) at z_end
        t = range(0.0, 1.0, length = trajectory_points)
        branch = zeros(trajectory_points, 3)
        for (i, ti) in enumerate(t)
            # Use a smooth cubic interpolation for x,y to avoid sharp kinks
            s = 3 * ti^2 - 2 * ti^3  # Hermite basis (smooth 0→1)
            branch[i, 1] = s * x_end
            branch[i, 2] = s * y_end
            branch[i, 3] = deviation_depth + ti * (z_end - deviation_depth)
        end
        push!(branches, branch)
    end

    # The producer well returns from the last point of each branch back to
    # the surface. We pick the first branch's endpoint as the return well
    # location for simplicity. The return well goes straight up.
    return_x = collar_xy[1, 1]
    return_y = collar_xy[2, 1]
    return_depth = trunk_depth + branch_dz
    n_return = max(Int(ceil(return_depth / trunk_step)), 2)
    return_z = range(return_depth, 0.0, length = n_return + 1)
    producer = zeros(length(return_z), 3)
    producer[:, 1] .= return_x
    producer[:, 2] .= return_y
    producer[:, 3] .= collect(return_z)

    # Assemble well_coords: [trunk, branch_1, ..., branch_n, producer]
    well_coords = Vector{Matrix{Float64}}()
    push!(well_coords, trunk)
    for branch in branches
        push!(well_coords, branch)
    end
    push!(well_coords, producer)

    # Connectivity: trunk connects to all branches, all branches connect to
    # the producer. Section indices: trunk = 1, branches = 2:(n+1), producer = n+2
    n_sections = n_branches + 2
    well_connectivity = zeros(Int, n_sections, 2)
    for b in 1:n_branches
        well_connectivity[b + 1, 1] = 1          # branch receives from trunk
        well_connectivity[b + 1, 2] = n_sections  # branch feeds into producer
    end

    return well_coords, well_connectivity

end

"""
    get_branch_constraints(well_coords; hxy_min)

Generate 2D mesh constraints from multi-branch well trajectories. Each unique
(x,y) footprint is added as a constraint for mesh refinement.
"""
function get_branch_constraints(well_coords; hxy_min)

    Δ = hxy_min / 2
    well_coords_2x = []
    for wc in well_coords
        wc_left = copy(wc)
        wc_left[:, 2] .-= Δ / 2
        push!(well_coords_2x, wc_left)
        wc_right = copy(wc)
        wc_right[:, 2] .+= Δ / 2
        push!(well_coords_2x, wc_right)
    end

    cell_constraints = Vector{Matrix{Float64}}()
    for wc in well_coords_2x
        cc_new = unique(wc[:, 1:2], dims = 1)'  # 2×N matrix
        if !isempty(cell_constraints)
            for cc in cell_constraints
                remove = falses(size(cc_new, 2))
                for (k, x) in enumerate(eachcol(cc_new))
                    if any(norm(x - y, 2) < Δ for y in eachcol(cc))
                        remove[k] = true
                    end
                end
                cc_new = cc_new[:, .!remove]
            end
        end
        if size(cc_new, 2) > 0
            push!(cell_constraints, cc_new)
        end
    end

    return cell_constraints

end

"""
    setup_branch_wells(domain, well_coords, well_connectivity)

Set up multi-segment wells for a branching well system. Follows the same
pattern as `setup_ags_wells` to create a supply well that traverses all
sections and a return well at the endpoint.
"""
function setup_branch_wells(domain, well_coords, well_connectivity)

    msh = physical_representation(domain)
    geo = tpfv_geometry(msh)
    wells = []
    directions = []
    wcells = []
    rcells = []
    neighbors = []
    cell_to_section = []
    wc0 = 0
    for (k, wc) in enumerate(well_coords)
        println("Processing well section $k/$(length(well_coords))")

        rc, extra = Jutul.find_enclosing_cells(
            msh, wc, n = 1000; geometry = geo, extra_out = true)
        dir = Vector.(extra[:direction] .* extra[:lengths])

        wc_idx = collect(1:length(rc)) .+ wc0
        wc0 += length(wc_idx)
        push!(rcells, rc)
        push!(wcells, wc_idx)
        push!(directions, dir)
        push!(neighbors, vcat(wc_idx[1:end-1]', wc_idx[2:end]'))
        push!(cell_to_section, fill(k, length(rc)))
        println("  Number of cells in section $k: ", length(rc))
    end

    connection_segments = []
    seg_no = 0
    for k in 1:length(well_coords)
        from = well_connectivity[k, 1]
        if from != 0
            ix = findfirst(rcells[from] .== rcells[k][1])
            if isnothing(ix)
                error("Could not find connection from well section $from to well $k")
            end
            in_seg = [wcells[from][ix], wcells[k][1]]
            neighbors[k] = hcat(in_seg, neighbors[k])
        end
        to = well_connectivity[k, 2]
        if to != 0
            ix = findfirst(rcells[to] .== rcells[k][end])
            if isnothing(ix)
                error("Could not find connection from well section $k to well $to")
            end
            out_seg = [wcells[k][end], wcells[to][ix]]
            neighbors[k] = hcat(neighbors[k], out_seg)
        end
        push!(connection_segments, [1, size(neighbors[k], 2)] .+ seg_no)
        seg_no += size(neighbors[k], 2)
    end

    rcells = vcat(rcells...)
    wcells = vcat(wcells...)
    directions = vcat(directions...)
    neighbors = hcat(neighbors...)
    cell_to_section = vcat(cell_to_section...)

    common_args = (
        simple_well = false,
        type = :closed_loop,
        radius = 150e-3,
        WI = 0.0
    )

    w_supply = setup_well(domain, rcells;
        neighborship = neighbors,
        perforation_cells_well = wcells,
        dir = directions,
        well_cell_centers = geo.cell_centroids[:, rcells],
        end_nodes = [length(wcells)],
        name = :CoaxialBranch_supply,
        common_args...
    )

    w_return = setup_well(domain, rcells[end];
        name = :CoaxialBranch_return,
        WIth = 0.0,
        common_args...
    )

    wells = [w_supply, w_return]

    section_info = (
        cell_to_section = cell_to_section,
        connection_segments = connection_segments
    )

    return wells, section_info

end
