"""
    coaxial_well_branches(; kwargs...)

Create a deep coaxial geothermal well simulation case.

The system consists of a single coaxial closed-loop well defined by a general
well trajectory given as an m×3 matrix. The well is set up using
`setup_closed_loop_well` with `closed_loop_type = :coaxial`.

# Keyword Arguments

## Well trajectory
- `well_trajectory`: An m×3 matrix defining the well path as (x, y, z) 
  coordinates [m]. Defaults to a vertical well from the surface to 2500 m
  depth.

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
- `rate = 50 m³/h`: Flow rate through the well [m³/s].
- `temperature_inj = 25°C`: Injection temperature [K].
- `num_years = 30`: Total simulation time [years].
- `report_interval = year/4`: Reporting interval [s].
- `schedule_args`: Additional keyword arguments passed to `make_schedule`.

## Mesh parameters
- `hz`: Vertical cell sizes per layer [m].
- `hxy_min = 25.0`: Minimum horizontal cell size [m].
- `hxy_max = 500.0`: Maximum horizontal cell size [m].
- `mesh_args`: Additional keyword arguments passed to `extruded_mesh`.

## Coaxial well parameters
- `well_name = :CoaxialWell`: Name prefix for the well.
- `well_args`: Additional keyword arguments passed to `setup_closed_loop_well`.

# Returns
A `JutulCase` for coaxial geothermal well simulation.
"""
function coaxial_well_branches(;
    # Well trajectory (m×3 matrix)
    well_trajectory = default_coaxial_trajectory(),
    # Geological parameters
    depths = [0.0, 500.0, 1500.0, 2000.0, 2500.0, 3000.0],
    permeability = [1e-3, 1e-3, 1e-3, 1e-2, 1e-3]*darcy,
    porosity = [0.01, 0.01, 0.01, 0.05, 0.01],
    rock_thermal_conductivity = [2.5, 2.5, 2.8, 3.5, 2.5]*watt/(meter*Kelvin),
    rock_heat_capacity = [900, 900, 900, 900, 900]*joule/(kilogram*Kelvin),
    rock_density = [2600, 2600, 2600, 2600]*kilogram/meter^3,
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
    hxy_min = 5.0,
    hxy_max = 10*hxy_min,
    mesh_args = NamedTuple(),
    # Well parameters
    well_name = :CoaxialWell,
    well_args = NamedTuple()
)
    @assert size(well_trajectory, 2) == 3 "well_trajectory must be an m×3 matrix"
    @assert size(well_trajectory, 1) >= 2 "well_trajectory must have at least 2 points"

    # ## Set up vertical mesh sizing
    if ismissing(hz)
        num_layers = length(depths) - 1
        dz = diff(depths)
        hz = fill(250.0, num_layers)
        for (i, d) in enumerate(dz)
            hz[i] = min(hz[i], d / 5)
        end
    end

    # ## Create mesh constraints from well trajectory
    constraints = get_well_constraints(well_trajectory; hxy_min = hxy_min)

    # return constraints

    # ## Create domain
    domain, layers, metrics = layered_reservoir_domain(constraints, depths,
        (
            permeability = permeability,
            porosity = porosity,
            rock_thermal_conductivity = rock_thermal_conductivity,
            rock_heat_capacity = rock_heat_capacity,
            rock_density = rock_density,
        );
        mesh_args = (;
            hz = hz,
            hxy_min = hxy_min,
            hxy_max = hxy_max,
            offset_rel = 1.0,
            mesh_args...
        )
    )

    # ## Set up the coaxial well
    mesh = physical_representation(domain)
    cells = Jutul.find_enclosing_cells(mesh, well_trajectory, n = 1_000)
    wells = setup_closed_loop_well(domain, cells;
        name = well_name,
        closed_loop_type = :coaxial,
        well_args...
    )

    model = setup_reservoir_model(
        domain, :geothermal; wells = wells)
    bc, state0, p, T = set_dirichlet_bcs(model;
        pressure_surface = 10atm,
        temperature_surface = temperature_surface,
        geothermal_gradient = thermal_gradient,
    )
    z = wells[1][:cell_centroids][3, :]
    state0[:CoaxialWell_supply][:Pressure] .= p(z)
    state0[:CoaxialWell_supply][:Temperature] .= T(z)

    # ## Set up controls
    rho = reservoir_model(model).system.rho_ref[1]
    supply_name = Symbol(well_name, "_supply")
    return_name = Symbol(well_name, "_return")
    
    rate_target = TotalRateTarget(rate)
    ctrl_inj = InjectorControl(
        rate_target, [1.0], density = rho, temperature = temperature_inj)
    bhp_target = BottomHolePressureTarget(5atm)
    ctrl_prod = ProducerControl(bhp_target)

    control = Dict(
        supply_name => ctrl_prod,
        return_name => ctrl_inj,
    )

    forces = setup_reservoir_forces(model, control = control, bc = bc)

    dt, forces = make_schedule(
        [forces], [(1,1), (1,1)];
        num_years = num_years,
        report_interval = report_interval,
        schedule_args...)

    # ## Additional case info
    info = Dict(
        :description => "Deep coaxial geothermal well system",
        :well_trajectory => well_trajectory,
    )

    # ## Create and return the complete simulation case
    case = JutulCase(model, dt, forces; state0 = state0, input_data = info)

    return case

end

"""
    default_coaxial_trajectory(; depth = 2500.0, step = 25.0)

Generate a default vertical well trajectory from the surface to the given
depth.

# Returns
An m×3 matrix of (x, y, z) coordinates [m], where z represents depth with
z = 0 at the surface.
"""
function default_coaxial_trajectory(; depth = 2500.0, step = 25.0)
    n = max(Int(ceil(depth / step)), 2)
    z = range(0.0, depth, length = n + 1)
    trajectory = zeros(length(z), 3)
    trajectory[:, 3] .= collect(z)
    return trajectory
end

"""
    get_well_constraints(well_trajectory; hxy_min)

Generate 2D mesh constraints from a well trajectory. The unique (x,y)
footprint of the trajectory is used as constraints for mesh refinement.
"""
function get_well_constraints(well_trajectory; hxy_min)

    Δ = hxy_min / 2
    well_coords_2x = []
    wc_left = copy(well_trajectory)
    wc_left[:, 1] .-= Δ / 2
    push!(well_coords_2x, wc_left)
    wc_right = copy(well_trajectory)
    wc_right[:, 1] .+= Δ / 2
    push!(well_coords_2x, wc_right)

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
