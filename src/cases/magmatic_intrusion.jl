"""
    magmatic_intrusion(; <keyword arguments>)

High-enthalpy geothermal system heated by a magmatic intrusion, based on the
production and reinjection scenarios of Yapparova et al. (2023).

The domain is a box with a cylindrical intrusion centered horizontally. The
background temperature follows a linear geothermal gradient, and the intrusion
is initialized at `temperature_intrusion` (or the background temperature, if
higher). The fluid is pure H₂O modelled with the pressure-enthalpy formulation
(`H2OSystem`), so the intrusion and its surroundings may be liquid, vapor,
two-phase or supercritical.

The simulation consists of two periods:
1. A natural-state period of `natural_state_time` with all wells shut, during
   which the system develops a hydrothermal convection pattern around the
   intrusion.
2. A production period of `num_years` where a producer is operated at constant
   bottom-hole pressure, optionally with a cold-water injector that is active
   for the first `injection_years`.

Initial pressure and boundary pressures are hydrostatic, computed by
integrating the density of water along the background temperature profile.

# Keyword arguments

## Geometry and mesh
- `dims = (48, 48, 80)`: Number of cells in x, y and z.
- `domain_size = (10000.0, 10000.0, 5000.0).*meter`: Domain extent in x, y and z.

## Geology and initial conditions
- `pressure_surface = 10.0bar`: Pressure at the top of the domain.
- `temperature_surface = convert_to_si(10.0, :Celsius)`: Temperature at the
  top of the domain.
- `geothermal_gradient = 0.0375Kelvin/meter`: Background geothermal gradient.
- `intrusion_radius = 1500.0meter`: Radius of the cylindrical intrusion.
- `intrusion_depths = [2000.0, 3000.0].*meter`: Top and bottom depth of the
  intrusion.
- `temperature_intrusion = convert_to_si(500.0, :Celsius)`: Initial
  temperature of the intrusion.
- `permeability = 1e-15`: Host-rock permeability [m²].
- `permeability_cap = 1e-16`: Cap-rock permeability [m²].
- `cap_thickness = 500.0meter`: Thickness of the low-permeability cap rock at
  the top of the domain.
- `porosity = 0.1`: Rock porosity.
- `rock_density = 2700.0kilogram/meter^3`: Rock density.
- `rock_heat_capacity = 880.0joule/kilogram/Kelvin`: Rock heat capacity.
- `rock_thermal_conductivity = 2.0watt/meter/Kelvin`: Rock thermal
  conductivity.
- `boundaries = :all`: Boundaries with fixed hydrostatic pressure and
  background temperature. Either `:all` or a subset of `[:top, :sides,
  :bottom]`. Boundaries that are not included are closed.

## Wells
- `producer_position = (-250.0, -250.0).*meter`: Horizontal producer position
  relative to the intrusion center.
- `producer_depths = [750.0, 1750.0].*meter`: Top and bottom depth of the
  producer's open interval.
- `injector_position = nothing`: Horizontal injector position relative to the
  intrusion center, e.g. `(-100.0, 250.0).*meter`. If `nothing`, no injector is
  added.
- `injector_depths = [2050.0, 2300.0].*meter`: Top and bottom depth of the
  injector's open interval.
- `bhp_producer = 50.0bar`: Producer bottom-hole pressure.
- `bhp_injector = 300.0bar`: Injector bottom-hole pressure.
- `temperature_injection = convert_to_si(80.0, :Celsius)`: Injection
  temperature.

## Schedule
- `natural_state_time = 2400.0year`: Duration of the natural-state period.
- `dt_natural_state = 100.0year`: Target timestep during the natural-state
  period.
- `num_years = 100`: Duration of the production period in years.
- `injection_years = num_years`: Number of years the injector is active.
- `report_interval = year/4`: Target timestep during the production period.

# Returns

A `JutulCase` ready for `simulate_reservoir`.

# Reference

Yapparova, A., Lamy-Chappuis, B., Scott, S. W., Gunnarsson, G., & Driesner, T.
(2023). Cold water injection near the magmatic heat source can enhance
production from high-enthalpy geothermal fields. *Geothermics*, 112, 102744.
https://doi.org/10.1016/j.geothermics.2023.102744
"""
function magmatic_intrusion(;
    dims = (48, 48, 80),
    domain_size = (10000.0, 10000.0, 5000.0).*meter,
    pressure_surface = 10.0*si_unit(:bar),
    temperature_surface = convert_to_si(10.0, :Celsius),
    geothermal_gradient = 0.0375Kelvin/meter,
    intrusion_radius = 1500.0meter,
    intrusion_depths = [2000.0, 3000.0].*meter,
    temperature_intrusion = convert_to_si(500.0, :Celsius),
    permeability = 1e-15,
    permeability_cap = 1e-16,
    cap_thickness = 500.0meter,
    porosity = 0.1,
    rock_density = 2700.0kilogram/meter^3,
    rock_heat_capacity = 880.0joule/kilogram/Kelvin,
    rock_thermal_conductivity = 2.0watt/meter/Kelvin,
    boundaries = :all,
    producer_position = (-250.0, -250.0).*meter,
    producer_depths = [750.0, 1750.0].*meter,
    injector_position = nothing,
    injector_depths = [2050.0, 2300.0].*meter,
    bhp_producer = 50.0*si_unit(:bar),
    bhp_injector = 300.0*si_unit(:bar),
    temperature_injection = convert_to_si(80.0, :Celsius),
    natural_state_time = 2400.0year,
    dt_natural_state = 100.0year,
    num_years = 100,
    injection_years = num_years,
    report_interval = year/4,
)
    if boundaries == :all
        boundaries = [:top, :sides, :bottom]
    end
    all(b -> b in (:top, :sides, :bottom), boundaries) ||
        throw(ArgumentError("boundaries must be :all or a subset of [:top, :sides, :bottom]"))

    sys = H2OSystem()
    tables = sys.pvt_tables

    # ## Mesh and rock properties
    g = CartesianMesh(dims, domain_size)
    geo = tpfv_geometry(g)
    x, y, z = (geo.cell_centroids[i, :] for i in 1:3)
    center = domain_size[1:2]./2

    T_background = depth -> temperature_surface + geothermal_gradient*depth
    in_intrusion = @. (intrusion_depths[1] <= z <= intrusion_depths[2]) &&
        (x - center[1])^2 + (y - center[2])^2 <= intrusion_radius^2

    perm = fill(permeability, length(z))
    perm[z .<= cap_thickness] .= permeability_cap

    domain = reservoir_domain(g;
        permeability = perm,
        porosity = porosity,
        rock_density = rock_density,
        rock_heat_capacity = rock_heat_capacity,
        rock_thermal_conductivity = rock_thermal_conductivity,
    )

    # ## Wells
    well_cells = (pos, depths) -> magmatic_intrusion_well_cells(
        g, center .+ pos, depths)
    wells = [setup_well(domain, well_cells(producer_position, producer_depths);
        name = :Producer, simple_well = true, use_top_node = true)]
    well_names = [:Producer]
    has_injector = !isnothing(injector_position)
    if has_injector
        push!(wells, setup_well(domain, well_cells(injector_position, injector_depths);
            name = :Injector, simple_well = true, use_top_node = true))
        push!(well_names, :Injector)
    end

    # ## Model
    model, parameters = setup_reservoir_model(domain, sys;
        wells = wells,
        thermal = true,
        block_backend = true,
        extra_out = true,
    )
    rmodel = reservoir_model(model)
    push!(rmodel.output_variables, :PhaseMassDensities, :PhaseViscosities)

    # ## Initial state
    p_hydrostatic = hydrostatic_pressure_h2o(tables, pressure_surface,
        T_background, domain_size[3])
    p0 = p_hydrostatic.(z)
    T0 = T_background.(z)
    T0[in_intrusion] .= max.(T0[in_intrusion], temperature_intrusion)
    state0 = setup_reservoir_state(model,
        Pressure = p0,
        Temperature = T0,
        Enthalpy = tables[:enthalpy].(p0, T0),
    )

    # ## Boundary conditions
    # Boundary fluxes do not include gravity, so boundary pressures are taken
    # at the boundary-cell centroids to be in equilibrium with the initial
    # state. Temperatures are taken at the boundary faces.
    bc_cells, bc_dir, bc_depth = Int[], Symbol[], Float64[]
    z_top, z_bottom = extrema(geo.boundary_centroids[3, :])
    for (f, c) in enumerate(geo.boundary_neighbors)
        dir = (:x, :y, :z)[argmax(abs.(geo.boundary_normals[:, f]))]
        zf = geo.boundary_centroids[3, f]
        if dir == :z
            kind = isapprox(zf, z_top) ? :top : :bottom
        else
            kind = :sides
        end
        kind in boundaries || continue
        push!(bc_cells, c)
        push!(bc_dir, dir)
        push!(bc_depth, zf)
    end
    bc_pressure = p_hydrostatic.(z[bc_cells])
    bc_temperature = T_background.(bc_depth)
    bc_enthalpy = tables[:enthalpy].(bc_pressure, bc_temperature)
    bc_density = tables[:density_mix].(bc_pressure, bc_enthalpy)
    bc = []
    for d in unique(bc_dir)
        ix = bc_dir .== d
        append!(bc, flow_boundary_condition(bc_cells[ix], domain,
            bc_pressure[ix], bc_temperature[ix];
            density = bc_density[ix], enthalpy = bc_enthalpy[ix], dir = d))
    end
    bc = [b for b in bc]
    isempty(bc) && (bc = nothing)

    # ## Forces
    ctrl_prod = ProducerControl(BottomHolePressureTarget(bhp_producer))
    ctrl_inj = InjectorControl(BottomHolePressureTarget(bhp_injector), [1.0];
        temperature = temperature_injection,
        enthalpy = tables[:enthalpy](bhp_injector, temperature_injection),
        density = tables[:density_mix](bhp_injector,
            tables[:enthalpy](bhp_injector, temperature_injection)),
        check = false,
    )
    shut = Dict{Symbol, Any}(w => DisabledControl() for w in well_names)
    forces_natural = setup_reservoir_forces(model; bc = bc, control = shut)
    control = Dict{Symbol, Any}(:Producer => ctrl_prod)
    has_injector && (control[:Injector] = ctrl_inj)
    forces_inj = setup_reservoir_forces(model; bc = bc, control = control)
    control_no_inj = Dict{Symbol, Any}(:Producer => ctrl_prod)
    has_injector && (control_no_inj[:Injector] = DisabledControl())
    forces_prod = setup_reservoir_forces(model; bc = bc, control = control_no_inj)

    # ## Schedule
    dt, forces = Float64[], []
    if natural_state_time > 0
        dt_ns = _benchmark_rampup_timesteps(Float64(natural_state_time),
            Float64(dt_natural_state))
        append!(dt, dt_ns)
        append!(forces, fill(forces_natural, length(dt_ns)))
    end
    if num_years > 0
        dt_prod = _benchmark_rampup_timesteps(Float64(num_years*year),
            Float64(report_interval))
        t_end = cumsum(dt_prod)
        append!(dt, dt_prod)
        append!(forces, [t <= injection_years*year + 1.0 ? forces_inj : forces_prod
            for t in t_end])
    end

    # ## Additional case info
    info = Dict{Symbol, Any}(
        :description => "Magmatic intrusion case set up using Fimbul.magmatic_intrusion()",
        :intrusion_cells => findall(in_intrusion),
        :natural_state_steps => count(f -> f === forces_natural, forces),
    )

    return JutulCase(model, dt, forces;
        state0 = state0,
        parameters = parameters,
        input_data = info,
    )
end

# Cells of a vertical well in a Cartesian mesh at horizontal position `xy`,
# open between `depths[1]` and `depths[2]`.
function magmatic_intrusion_well_cells(g, xy, depths)
    nx, ny, nz = g.dims
    Δ = g.deltas
    ijk = (v, d, n) -> clamp(Int(floor(v/d)) + 1, 1, n)
    i = ijk(xy[1], Δ[1], nx)
    j = ijk(xy[2], Δ[2], ny)
    k1 = ijk(depths[1], Δ[3], nz)
    k2 = ijk(depths[2], Δ[3], nz)
    return [(i, j, k) for k in k1:k2]
end

"""
    hydrostatic_pressure_h2o(tables, p_top, T, depth; dz = 1.0)

Compute a hydrostatic pressure profile for pure H₂O with temperature profile
`T(depth)` by integrating `dp/dz = ρ(p, T)g` from `p_top` at depth zero down to
`depth`. Returns a function of depth.
"""
function hydrostatic_pressure_h2o(tables, p_top, T, depth; dz = 1.0)
    rho = (p, z) -> tables[:density_mix](p, tables[:enthalpy](p, T(z)))
    n = max(ceil(Int, depth/dz), 1) + 1
    zv = collect(range(0.0, depth, length = n))
    pv = similar(zv)
    pv[1] = p_top
    for i in 2:n
        h = zv[i] - zv[i-1]
        rho_0 = rho(pv[i-1], zv[i-1])
        rho_1 = rho(pv[i-1] + rho_0*gravity_constant*h, zv[i])
        pv[i] = pv[i-1] + 0.5*(rho_0 + rho_1)*gravity_constant*h
    end
    return get_1d_interpolator(zv, pv)
end
