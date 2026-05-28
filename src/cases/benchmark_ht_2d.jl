"""
    benchmark_ht_2d(; <keyword arguments>)

2D high-temperature geothermal benchmark derived from the upper-crust plume
examples in Weis et al. (2014).

The domain is a 9 km wide by 3 km deep vertical section with unit thickness.
Pressure is initialized hydrostatically from atmospheric pressure at the
surface, and temperature is initialized uniformly at 10 °C. The top boundary is
held at atmospheric pressure and 10 °C. The implemented benchmark cases are the
magmatic-fluid-source variants, where hot H₂O is injected into the bottom
center of the domain as a single-cell point source. A moderate-enthalpy source
case reproduces the single-phase plume benchmark, while a high-enthalpy source
case reproduces the two-phase plume benchmark. Reported horizontal coordinates
are centered so that the injection region is located at x = 0. The alternative
bottom heat-flux variant is not included here because heat-flux boundary
conditions are not yet supported in Fimbul.

## Keyword arguments

- `benchmark_case = :fluid_source`: Case selector. `:fluid_source` is kept as
    an alias for the moderate-enthalpy single-phase source case.
    `:single_phase_source`, `:two_phase_source`, and `:test` are supported.
- `nx = 90`: Number of cells in the horizontal direction.
- `nz = 60`: Number of cells in the vertical direction.
- `domain_width = 9000.0*meter`: Horizontal domain extent [m].
- `domain_depth = 3000.0*meter`: Vertical domain extent [m].
- `domain_thickness = 1.0*meter`: Out-of-plane thickness [m].
- `top_pressure = 1.0*atm`: Surface pressure [Pa].
- `top_temperature = convert_to_si(10.0, :Celsius)`: Surface temperature [K].
- `source_mass_rate = 0.04*kilogram/second`: Magmatic H₂O source rate [kg/s].
- `source_enthalpy = nothing`: Injected specific enthalpy [J/kg]. Defaults to
    `0.5e6*joule/kilogram` for the single-phase source case and
    `1.5e6*joule/kilogram` for the two-phase source case.
- `source_i = cld(nx, 2)`: Horizontal source-cell index.
- `source_k = nz`: Vertical source-cell index.
- `permeability = 1e-15`: Rock permeability [m²].
- `porosity = 0.1`: Rock porosity [-].
- `rock_density = 2700.0`: Rock density [kg/m³].
- `rock_heat_capacity = 880.0`: Rock heat capacity [J/(kg·K)].
- `rock_thermal_conductivity = 2.0`: Rock thermal conductivity [W/(m·K)].
- `total_time = 1000.0*year`: Total simulated time [s].
- `dt_target = 5.0*year`: Target report step after ramp-up [s].
- `enthalpy_tables = build_steam_tables_2ph()`: Pre-built (P,H) steam tables.

## Returns

A `JutulCase` ready for `simulate_reservoir`.

## Reference

Weis, P., Driesner, T., & Heinrich, C. A. (2014). Hydrothermal, multiphase
convection of H₂O-NaCl fluids from ambient to magmatic temperatures: A new
numerical scheme and benchmarks for code comparison. *Geofluids*, 14(3),
347–371. https://doi.org/10.1111/gfl.12080
"""
function benchmark_ht_2d(;
    benchmark_case :: Symbol = :fluid_source,
    nx :: Int = 90,
    nz :: Int = 60,
    domain_width :: Float64 = 9000.0*meter,
    domain_depth :: Float64 = 3000.0*meter,
    domain_thickness :: Float64 = 1.0*meter,
    top_pressure :: Float64 = 1.0*atm,
    top_temperature :: Float64 = convert_to_si(10.0, :Celsius),
    source_mass_rate :: Float64 = 0.04*kilogram/second,
    source_enthalpy :: Union{Nothing, Float64} = nothing,
    source_i :: Int = cld(nx, 2),
    source_k :: Int = nz,
    permeability :: Float64 = 1e-15,
    porosity :: Float64 = 0.1,
    rock_density :: Float64 = 2700.0,
    rock_heat_capacity :: Float64 = 880.0,
    rock_thermal_conductivity :: Float64 = 2.0,
    total_time :: Float64 = 5000.0*year,
    dt_target :: Float64 = 10.0*year,
    enthalpy_tables = build_steam_tables_2ph(),
)
    benchmark_case == :fluid_source && (benchmark_case = :single_phase_source)
    benchmark_case in (:single_phase_source, :two_phase_source, :test) ||
        throw(ArgumentError("benchmark_case must be :fluid_source, :single_phase_source, :two_phase_source, or :test"))
    1 <= source_i <= nx || throw(ArgumentError("source_i must satisfy 1 <= source_i <= nx"))
    1 <= source_k <= nz || throw(ArgumentError("source_k must satisfy 1 <= source_k <= nz"))

    if source_enthalpy === nothing
        source_enthalpy = ifelse(benchmark_case == :two_phase_source, 1.5e6*joule/kilogram, 0.5e6*joule/kilogram)
    end

    source_regime = benchmark_case == :two_phase_source ? :two_phase : :single_phase

    if benchmark_case == :test
        nx = min(nx, 20)
        nz = min(nz, 12)
        domain_width = 2000.0*meter
        domain_depth = 1000.0*meter
        source_i = cld(nx, 2)
        source_k = nz
        total_time = 50.0*year
        dt_target = min(dt_target, 1.0*year)
    end

    g = CartesianMesh((nx, 1, nz), (domain_width, domain_thickness, domain_depth))
    domain = reservoir_domain(g;
        permeability = permeability,
        porosity = porosity,
        rock_density = rock_density,
        rock_heat_capacity = rock_heat_capacity,
        rock_thermal_conductivity = rock_thermal_conductivity,
    )

    source_well = setup_well(
        domain,
        [(source_i, 1, source_k)];
        name = :MagmaticSource,
        simple_well = true,
        use_top_node = true,
    )

    sys = H2OSystem(enthalpy_tables)
    model, parameters = setup_reservoir_model(
        domain,
        sys;
        wells = [source_well],
        thermal = true,
        block_backend = true,
        extra_out = true,
    )

    rmodel = reservoir_model(model)
    push!(rmodel.output_variables, :PhaseMassDensities, :PhaseViscosities)

    geo = tpfv_geometry(physical_representation(domain))
    z_cells = geo.cell_centroids[3, :]
    z0 = minimum(z_cells)
    depth_cells = z_cells .- z0

    rho_ref = rmodel.system.rho_ref[1]
    dp_dz = rho_ref*gravity_constant
    pressure0 = top_pressure .+ dp_dz.*depth_cells
    enthalpy0 = enthalpy_tables[:enthalpy].(pressure0, top_temperature)

    state0 = setup_reservoir_state(model,
        Pressure = pressure0,
        Enthalpy = enthalpy0,
        Temperature = top_temperature,
    )

    z_bdr = geo.boundary_centroids[3, :]
    top = isapprox.(z_bdr, minimum(z_bdr))
    bc_cells = geo.boundary_neighbors[top]
    bc_pressure = fill(top_pressure, length(bc_cells))
    bc_temperature = fill(top_temperature, length(bc_cells))
    bc_enthalpy = enthalpy_tables[:enthalpy].(bc_pressure, bc_temperature)
    bc_density = enthalpy_tables[:density_mix].(bc_pressure, bc_enthalpy)
    bc = flow_boundary_condition(
        bc_cells,
        domain,
        bc_pressure,
        bc_temperature;
        density = bc_density,
        enthalpy = bc_enthalpy,
    )

    source_cell = cell_index(g, (source_i, 1, source_k))
    source_x = geo.cell_centroids[1, source_cell]
    source_pressure = pressure0[source_cell]
    source_temperature = enthalpy_tables[:temperature](source_pressure, source_enthalpy)
    source_saturation_vapor = enthalpy_tables[:saturation_vapor_ph](source_pressure, source_enthalpy)
    source_density_liquid = enthalpy_tables[:density_liquid_ph](source_pressure, source_enthalpy)
    source_density_vapor = enthalpy_tables[:density_vapor_ph](source_pressure, source_enthalpy)
    source_density = (1.0 - source_saturation_vapor)*source_density_liquid +
        source_saturation_vapor*source_density_vapor
    source_rate = source_mass_rate/source_density

    control = Dict(
        :MagmaticSource => InjectorControl(
            TotalRateTarget(source_rate),
            [1.0];
            density = source_density,
            temperature = source_temperature,
            enthalpy = source_enthalpy,
            check = false,
        ),
    )
    forces = setup_reservoir_forces(model; bc = bc, control = control)

    dt = _benchmark_rampup_timesteps(total_time, dt_target)
    info = Dict{Symbol, Any}(
        :description => "2D high-temperature geothermal benchmark case $(benchmark_case) (Weis et al. 2014)",
        :benchmark_case => benchmark_case,
        :source_regime => source_regime,
        :domain_width => domain_width,
        :domain_depth => domain_depth,
        :top_pressure => top_pressure,
        :top_temperature => top_temperature,
        :source_cell => source_cell,
        :x_coordinate_origin => source_x,
        :source_mass_rate => source_mass_rate,
        :source_enthalpy => source_enthalpy,
        :source_energy_rate => source_mass_rate*source_enthalpy,
        :total_time => total_time,
    )

    return JutulCase(model, dt, forces;
        state0 = state0,
        parameters = parameters,
        input_data = info,
    )
end
