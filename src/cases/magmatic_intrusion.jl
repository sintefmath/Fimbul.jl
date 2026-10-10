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

Only wells that are active at some point in the schedule are included in the
model, since fluid can flow along the wellbore of a shut well. To compute the
natural state without wells, set `num_years = 0`, and use the resulting
reservoir state as `initial_state` for a production case with
`natural_state_time = 0`.

Initial pressure and boundary pressures are hydrostatic, computed by
integrating the density of water along the background temperature profile.

# Keyword arguments

## Geometry and mesh
- `dims = (48, 48, 80)`: Number of cells in x, y and z.
- `domain_size = (10000.0, 10000.0, 5000.0).*meter`: Domain extent in x, y and z.
- `refinement_margin = 250.0meter`: Distance outside the intrusion within which
  the mesh has uniform, fine cells.
- `cell_growth = 1.25`: Approximate size ratio between neighboring cells outside
  the refined region. The mesh is a tensor grid with the number of cells given
  by `dims`, with uniform cells in the refined region and cells growing
  geometrically towards the domain boundaries. Use `cell_growth = 1.0` for a
  uniform mesh.
- `cap_cells = 3`: Number of uniform cell layers in the cap rock for the tensor
  grid. Cells between the cap rock and the refined region are no larger than
  the cap-rock cells. If zero, the cap rock is not treated separately.

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
  intrusion center, e.g. `(-100.0, 250.0).*meter`. If `nothing` (or if
  `injection_years = 0`), no injector is added.
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
- `num_years = 100`: Duration of the production period in years. If zero, the
  model has no wells.
- `injection_years = num_years`: Number of years the injector is active.
- `report_interval = year/4`: Target timestep during the production period.

## Restart
- `initial_state = nothing`: Reservoir state to start from, e.g. the final
  state of a natural-state simulation. Must contain `:Pressure` and
  `:Enthalpy` (and optionally `:Temperature`) for all cells. If `nothing`, the
  hydrostatic initial state described above is used.

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
    refinement_margin = 250.0meter,
    cell_growth = 1.25,
    cap_cells = 3,
    pressure_surface = 10.0*si_unit(:bar),
    temperature_surface = convert_to_si(10.0, :Celsius),
    geothermal_gradient = 0.0375Kelvin/meter,
    intrusion_radius = 1500.0meter,
    intrusion_depths = [2000.0, 3000.0].*meter,
    temperature_intrusion = convert_to_si(900.0, :Celsius),
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
    initial_state = nothing,
)
    if boundaries == :all
        boundaries = [:top, :sides, :bottom]
    end
    all(b -> b in (:top, :sides, :bottom), boundaries) ||
        throw(ArgumentError("boundaries must be :all or a subset of [:top, :sides, :bottom]"))

    sys = H2OSystem()
    tables = sys.pvt_tables

    # ## Mesh and rock properties
    center = domain_size[1:2]./2
    if cell_growth > 1.0
        r = intrusion_radius + refinement_margin
        sizes = (
            magmatic_intrusion_cell_sizes(domain_size[1], dims[1],
                center[1] - r, center[1] + r; growth = cell_growth),
            magmatic_intrusion_cell_sizes(domain_size[2], dims[2],
                center[2] - r, center[2] + r; growth = cell_growth),
            magmatic_intrusion_cell_sizes(domain_size[3], dims[3],
                intrusion_depths[1] - refinement_margin,
                intrusion_depths[2] + refinement_margin;
                growth = cell_growth,
                top = cap_cells > 0 ? (cap_thickness, cap_cells) : nothing),
        )
        g = CartesianMesh(dims, sizes)
    else
        g = CartesianMesh(dims, domain_size)
    end
    geo = tpfv_geometry(g)
    x, y, z = (geo.cell_centroids[i, :] for i in 1:3)

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
    # Only wells that are active at some point in the schedule are added to the
    # model. Fluid can flow through the wellbore of a shut well, which would
    # disturb the natural state.
    well_cells = (pos, depths) -> magmatic_intrusion_well_cells(
        g, center .+ pos, depths)
    has_producer = num_years > 0
    has_injector = has_producer && !isnothing(injector_position) &&
        injection_years > 0
    wells, well_names = [], Symbol[]
    if has_producer
        push!(wells, setup_well(domain, well_cells(producer_position, producer_depths);
            name = :Producer, simple_well = true, use_top_node = true))
        push!(well_names, :Producer)
    end
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
    if isnothing(initial_state)
        p0 = p_hydrostatic.(z)
        T0 = T_background.(z)
        T0[in_intrusion] .= max.(T0[in_intrusion], temperature_intrusion)
        h0 = tables[:enthalpy].(p0, T0)
    else
        p0 = initial_state[:Pressure]
        h0 = initial_state[:Enthalpy]
        T0 = haskey(initial_state, :Temperature) ?
            initial_state[:Temperature] : tables[:temperature].(p0, h0)
    end
    state0 = setup_reservoir_state(model,
        Pressure = p0,
        Temperature = T0,
        Enthalpy = h0,
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
    if has_producer
        shut = Dict{Symbol, Any}(w => DisabledControl() for w in well_names)
        forces_natural = setup_reservoir_forces(model; bc = bc, control = shut)
        control = Dict{Symbol, Any}(:Producer => ctrl_prod)
        has_injector && (control[:Injector] = ctrl_inj)
        forces_inj = setup_reservoir_forces(model; bc = bc, control = control)
        control_no_inj = Dict{Symbol, Any}(:Producer => ctrl_prod)
        has_injector && (control_no_inj[:Injector] = DisabledControl())
        forces_prod = setup_reservoir_forces(model; bc = bc, control = control_no_inj)
    else
        forces_natural = setup_reservoir_forces(model; bc = bc)
    end

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
    function ijk(v, d)
        n = g.dims[d]
        Δ = g.deltas[d]
        faces = cumsum(Δ isa Number ? fill(Δ, n) : Δ)
        return clamp(searchsortedfirst(faces, v), 1, n)
    end
    i = ijk(xy[1], 1)
    j = ijk(xy[2], 2)
    k1 = ijk(depths[1], 3)
    k2 = ijk(depths[2], 3)
    return [(i, j, k) for k in k1:k2]
end

# Cell sizes for `n` cells on `[0, L]`: uniform cells in `[a, b]`, and cells
# growing geometrically by approximately `growth` towards both ends. If `top =
# (thickness, n_top)` is given, the top layer `[0, thickness]` gets `n_top`
# uniform cells, and cells between the top layer and `a` are no larger than
# the top-layer cells.
function magmatic_intrusion_cell_sizes(L, n, a, b; growth = 1.25, top = nothing)
    a, b = clamp(a, 0.0, L), clamp(b, 0.0, L)
    W = b - a
    W > 0 || throw(ArgumentError("Refined region must have positive length"))
    if isnothing(top)
        L_top, n_top, h_top = 0.0, 0, Inf
    else
        L_top, n_top = top
        L_top < a || throw(ArgumentError("Top layer must be above the refined region"))
        h_top = L_top/n_top
    end
    # Number of cells (real-valued) needed to grade outwards over length Lo,
    # starting from cells of size h, with cell sizes limited by h_max
    function count(Lo, h, h_max)
        Lo <= 0 && return 0.0
        h_max = max(h_max, h)
        k = floor(log(h_max/h)/log(growth))
        S = h*growth*(growth^k - 1)/(growth - 1)
        S >= Lo && return log(1 + Lo*(growth - 1)/(h*growth))/log(growth)
        return k + (Lo - S)/h_max
    end
    total = h -> W/h + count(a - L_top, h, h_top) + n_top + count(L - b, h, Inf)
    lo, hi = 1e-6*L, Float64(L)
    for _ in 1:200
        h = sqrt(lo*hi)
        total(h) > n ? (lo = h) : (hi = h)
    end
    h = sqrt(lo*hi)
    m_a = ceil(Int, count(a - L_top, h, h_top) - 0.1)
    m_b = round(Int, count(L - b, h, Inf))
    n_fine = n - m_a - m_b - n_top
    n_fine > 0 || throw(ArgumentError("Too few cells for the refined region"))
    h_fine = W/n_fine
    left = magmatic_intrusion_graded_sizes(a - L_top, m_a, h_fine, h_top)
    right = magmatic_intrusion_graded_sizes(L - b, m_b, h_fine, Inf)
    return vcat(fill(h_top, n_top), reverse(left), fill(h_fine, n_fine), right)
end

# Sizes of `m` cells over length `Lo`, growing geometrically from a neighboring
# cell of size `h0`, with cell sizes limited by `h_max` where possible.
function magmatic_intrusion_graded_sizes(Lo, m, h0, h_max = Inf)
    m == 0 && return Float64[]
    m*h_max <= Lo && return fill(Lo/m, m)
    sizes = r -> [min(h0*r^k, h_max) for k in 1:m]
    lo, hi = 1e-3, 10.0
    for _ in 1:200
        r = (lo + hi)/2
        sum(sizes(r)) > Lo ? (hi = r) : (lo = r)
    end
    return sizes((lo + hi)/2)
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
