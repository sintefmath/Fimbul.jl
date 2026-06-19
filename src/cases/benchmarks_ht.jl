const HYDROTHERM_BENCHMARKS_ROOT = joinpath(artifact"HYDROTHERMBenchmarks", "hydrotherm")
const HYDROTHERM_1D_ROOT = joinpath(HYDROTHERM_BENCHMARKS_ROOT, "benchmark_1d")
const HYDROTHERM_2D_ROOT = joinpath(HYDROTHERM_BENCHMARKS_ROOT, "benchmark_2d")

available_vertical_modes_ht_1d(case_symbol) = case_symbol == :e ? (false,) : (false, true)

hydrotherm_1d_profile_path(case_symbol, vertical; root = HYDROTHERM_1D_ROOT) =
    joinpath(root, "case_$(case_symbol)_$(vertical ? "vertical" : "horizontal").txt")

hydrotherm_1d_timestep_path(case_symbol, vertical; root = HYDROTHERM_1D_ROOT) =
    joinpath(root, "case_$(case_symbol)_$(vertical ? "vertical" : "horizontal")_timesteps.txt")

function load_hydrotherm_1d_table(case_symbol, vertical; root = HYDROTHERM_1D_ROOT)
    path = hydrotherm_1d_profile_path(case_symbol, vertical; root)
    isfile(path) || return nothing
    header = split(strip(replace(readline(path), "#" => "")))
    data = readdlm(path, comments = true)
    return Dict(name => vec(data[:, i]) for (i, name) in enumerate(header))
end

function load_hydrotherm_1d_timesteps(case_symbol, vertical; root = HYDROTHERM_1D_ROOT)
    path = hydrotherm_1d_timestep_path(case_symbol, vertical; root)
    isfile(path) || return nothing
    data = readdlm(path, comments = true)
    return data[:, 2] .* si_unit(:year)
end

function load_hydrotherm_1d_property(case_symbol, vertical, column; root = HYDROTHERM_1D_ROOT)
    table = load_hydrotherm_1d_table(case_symbol, vertical; root)
    table === nothing && return nothing
    haskey(table, column) || return nothing
    coordinate_column = vertical ? "depth_m" : "distance_m"
    return (coordinate_m = table[coordinate_column], values = table[column])
end

function load_hydrotherm_1d_phase_path(case_symbol; root = HYDROTHERM_1D_ROOT)
    table = load_hydrotherm_1d_table(case_symbol, false; root)
    table === nothing && return nothing
    return (pressure_mpa = table["pressure_mpa"], enthalpy_kj_per_kg = table["enthalpy_kj_per_kg"])
end

function replace_case_timesteps(case, timesteps; check_sum = false, rtol = 1e-3)
    timesteps === nothing && return case
    if check_sum
        isapprox(sum(timesteps), sum(case.dt); rtol = rtol) ||
            error("HYDROTHERM timestep sum $(sum(timesteps)) does not match case duration $(sum(case.dt))")
    end

    (; model, forces, state0, parameters, input_data) = case
    return JutulCase(model, timesteps, forces; state0 = state0, parameters = parameters, input_data = input_data)
end

hydrotherm_2d_case_name(case) = String(case.input_data[:benchmark_case])
hydrotherm_2d_path(case_name, name; root = HYDROTHERM_2D_ROOT) = joinpath(root, "$(case_name)_$(name).txt")

load_hydrotherm_2d_vector(path) = vec(readdlm(path, comments = true))
load_hydrotherm_2d_matrix(path) = permutedims(readdlm(path, comments = true))

function load_hydrotherm_2d_timesteps(case; root = HYDROTHERM_2D_ROOT)
    path = hydrotherm_2d_path(hydrotherm_2d_case_name(case), "timesteps"; root)
    data = readdlm(path, comments = true)
    return data[:, 2] .* si_unit(:year)
end

function load_hydrotherm_2d_reference(case; root = HYDROTHERM_2D_ROOT)
    case_name = hydrotherm_2d_case_name(case)
    x_m = load_hydrotherm_2d_vector(hydrotherm_2d_path(case_name, "x_m"; root))
    z_m = load_hydrotherm_2d_vector(hydrotherm_2d_path(case_name, "z_m"; root))
    pressure_mpa = load_hydrotherm_2d_matrix(hydrotherm_2d_path(case_name, "pressure_mpa"; root))
    temperature_c = load_hydrotherm_2d_matrix(hydrotherm_2d_path(case_name, "temperature_c"; root))
    liquid_saturation = load_hydrotherm_2d_matrix(hydrotherm_2d_path(case_name, "liquid_saturation"; root))
    vapor_saturation = map(liquid_saturation) do value
        isnan(value) ? NaN : 1.0 - value
    end

    x0 = x_m[cld(length(x_m), 2)]
    return (
        x_km = (x_m .- x0) ./ 1e3,
        depth_km = z_m ./ 1e3,
        pressure = pressure_mpa,
        temperature = temperature_c,
        vapor_saturation = vapor_saturation,
    )
end

"""
    benchmark_ht_1d(; <keyword arguments>)

1D high-temperature geothermal benchmark cases from Weis et al. (2014).

Five cases (`:a`–`:e`) cover different initial pressure/temperature conditions
and pressure gradients, spanning single-phase liquid, single-phase vapor, and
two-phase flow regimes. Each case is a 1D column with Dirichlet-like
pressure/enthalpy conditions at both ends, imposed through flow boundary
conditions. The initial state uses a linear pressure profile between the two
boundary pressures and a uniform producer-side temperature.

| Case | p_inj [MPa] | p_prod [MPa] | T_inj [°C] | T_prod [°C] | vertical? |
|:----:|:-----------:|:------------:|:----------:|:-----------:|:--------:|
| :a   | 50          | 25           | 350        | 150         | optional |
| :b   | 40          | 20           | 450        | 300         | optional |
| :c   | 15          | 1            | 500        | 350         | optional |
| :d   | 20          | 1            | 400        | 150         | optional |
| :e   | 4           | 1            | 300        | 150         | no       |

Case `:test` is also accepted as a lightweight internal sanity-check case.

## Grid orientation

- `vertical = false`: flow in x, grid `(nx, 1, 1)`. Injector at cell 1 (left),
  producer at cell nx (right).
- `vertical = true`:  flow in z, grid `(1, 1, nx)`. Injector at cell nx
  (bottom = deep = high pressure), producer at cell 1 (top = shallow).

## Keyword arguments

- `benchmark_case = :a`: Case selector. Benchmark cases are `:a`, `:b`, `:c`,
    `:d`, `:e`; `:test` is a lightweight internal case.
- `vertical = false`: Enable vertical flow (not available for case `:e`).
- `nx = 200`: Number of cells along the active flow direction.
- `domain_length = 2000.0*meter`: Length of the 1D domain [m].
- `cell_size = 10.0*meter`: Cell width in the two inactive directions [m].
- `dt_target = 5.0*year`: Target timestep size [s] after ramp-up.

The total simulation time is selected internally from the benchmark case and
orientation:

- `:a`: 250 years horizontal, 750 years vertical.
- `:b`: 120 years horizontal, 350 years vertical.
- `:c`: 1500 years.
- `:d`: 200 years horizontal, 1000 years vertical.
- `:e`: 2000 years.

## Returns

A `JutulCase` ready for `simulate_reservoir`.

## Reference

Weis, P., Driesner, T., & Heinrich, C. A. (2014). Hydrothermal, multiphase
convection of H₂O-NaCl fluids from ambient to magmatic temperatures: A new
numerical scheme and benchmarks for code comparison. *Geofluids*, 14(3),
347–371. https://doi.org/10.1111/gfl.12080
"""
function benchmark_ht_1d(;
    benchmark_case  :: Symbol  = :a,
    vertical        :: Bool    = false,
    nx              :: Int     = 200,
    domain_length   :: Float64 = 2000.0*meter,
    cell_size       :: Float64 = 10.0*meter,
    dt_target       :: Float64 = 5.0*year,
)

    benchmark_case in (:a, :b, :c, :d, :e, :test) ||
        throw(ArgumentError("benchmark_case must be one of :a, :b, :c, :d, :e"))
    (benchmark_case == :e && vertical) &&
        throw(ArgumentError("Case :e is only defined without vertical flow (vertical = false)"))

    # ── Case parameters ────────────────────────────────────────────────────
    K0 = 273.15  # 0 °C in Kelvin

    if benchmark_case == :a
        p_inj  = 50e6;  p_prod = 25e6
        T_inj  = K0 + 350;  T_prod = K0 + 150
    elseif benchmark_case == :b
        p_inj  = 40e6;  p_prod = 20e6
        T_inj  = K0 + 450;  T_prod = K0 + 300
    elseif benchmark_case == :c
        p_inj  = 15e6;  p_prod = 1e6
        T_inj  = K0 + 500;  T_prod = K0 + 350
    elseif benchmark_case == :d
        p_inj  = 20e6;  p_prod = 1e6
        T_inj  = K0 + 400;  T_prod = K0 + 150
    elseif benchmark_case == :e
        p_inj  = 4e6;   p_prod = 1e6
        T_inj  = K0 + 300;  T_prod = K0 + 150
     else
        error("Invalid benchmark_case: $benchmark_case")
    end

    sys    = H2OSystem()
    tables = sys.pvt_tables

    H_inj  = tables[:enthalpy](p_inj, T_inj)
    H_prod = tables[:enthalpy](p_prod, T_prod)

    # ── Grid ──────────────────────────────────────────────────────────────
    #
    # Horizontal 1D flow in x  →  (nx, 1, 1)
    #   cell 1  = left  (inlet  / high p, high T)
    #   cell nx = right (outlet / low  p, low  T)
    #
    # Vertical 1D flow in z  →  (1, 1, nx)
    #   cell 1  = top    (shallow / low  p, low  T)  ← producer
    #   cell nx = bottom (deep    / high p, high T)  ← injector
    #
    if vertical
        g         = CartesianMesh((1, 1, nx), (cell_size, cell_size, domain_length))
        cell_inj  = nx   # bottom cell = deep = high pressure
        cell_prod = 1    # top    cell = shallow = low pressure
        p0 = collect(range(p_prod, p_inj; length = nx))
    else
        g         = CartesianMesh((nx, 1, 1), (domain_length, cell_size, cell_size))
        cell_inj  = 1    # left  cell = inlet
        cell_prod = nx   # right cell = outlet
        p0 = collect(range(p_inj, p_prod; length = nx))
    end
    h0 = tables[:enthalpy].(p0, T_prod)

    domain = reservoir_domain(g;
        permeability              = 1e-15,   # m²  (≈ 10 μD)
        porosity                  = 0.1,
        rock_density              = 2700.0,  # kg/m³
        rock_heat_capacity        = 880.0,   # J/(kg·K)
        rock_thermal_conductivity = 2.0,     # W/(m·K)
    )

    # ── Model ─────────────────────────────────────────────────────────────
    model, parameters = setup_reservoir_model(domain, sys; block_backend = true, extra_out = true, thermal = true)
    kr = BrooksCoreyRelativePermeabilities(sys, 1.0, [0.3, 0.0], 1.0)
    model = replace_variables!(model, RelativePermeabilities = kr)

    # ── Initial state ──────────────────────────────────────────────────────
    # Initialize with a linear pressure profile between the boundary pressures
    # and a uniform producer-side temperature.
    state0 = setup_reservoir_state(model,
        Pressure = p0,
        Enthalpy = h0,
        Temperature = T_prod,
    )

    bc = flow_boundary_condition(
        [cell_inj, cell_prod], domain, [p_inj, p_prod], [T_inj, T_prod];
        enthalpy = [H_inj, H_prod],
    )
    forces = setup_reservoir_forces(model; bc = bc)

    # ── Timesteps ──────────────────────────────────────────────────────────
    if benchmark_case == :a
        total_time = ifelse(!vertical, 250.0*year, 750.0*year)
    elseif benchmark_case == :b
        total_time = ifelse(!vertical, 120.0*year, 350.0*year)
    elseif benchmark_case == :c
        total_time = 1500.0*year
    elseif benchmark_case == :d
        total_time = ifelse(!vertical, 200.0*year, 1000.0*year)
    else  # :e
        total_time = 2000.0*year
    end
    dt_vec = _benchmark_rampup_timesteps(total_time, dt_target)

    # ── Assemble case ──────────────────────────────────────────────────────
    info = Dict{Symbol, Any}(
        :description    => "1D high-temperature geothermal benchmark case $(benchmark_case) (Weis et al. 2014)",
        :benchmark_case => benchmark_case,
        :total_time     => total_time,
        :vertical       => vertical,
        :p_inj          => p_inj,
        :p_prod         => p_prod,
        :T_inj          => T_inj,
        :T_prod         => T_prod,
    )
    return JutulCase(model, dt_vec, forces; state0 = state0, parameters=parameters, input_data = info)
end

# ──────────────────────────────────────────────────────────────────────────────

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
- `total_time = 5000.0*year`: Total simulated time [s].
- `dt_target = 10.0*year`: Target timestep size after ramp-up [s].

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
)
    benchmark_case == :fluid_source && (benchmark_case = :single_phase_source)
    benchmark_case in (:single_phase_source, :two_phase_source) ||
        throw(ArgumentError("benchmark_case must be :fluid_source, :single_phase_source, or :two_phase_source"))
    1 <= source_i <= nx || throw(ArgumentError("source_i must satisfy 1 <= source_i <= nx"))
    1 <= source_k <= nz || throw(ArgumentError("source_k must satisfy 1 <= source_k <= nz"))

    if source_enthalpy === nothing
        source_enthalpy = ifelse(benchmark_case == :two_phase_source, 1.5e6*joule/kilogram, 0.5e6*joule/kilogram)
    end

    source_regime = benchmark_case == :two_phase_source ? :two_phase : :single_phase

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
        WIth = 0.0
    )

    sys    = H2OSystem()
    tables = sys.pvt_tables

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
    enthalpy0 = tables[:enthalpy].(pressure0, top_temperature)

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
    bc_enthalpy = tables[:enthalpy].(bc_pressure, bc_temperature)
    bc_density = tables[:density_mix].(bc_pressure, bc_enthalpy)
    bc = flow_boundary_condition(
        bc_cells,
        domain,
        bc_pressure,
        bc_temperature;
        density = bc_density,
        enthalpy = bc_enthalpy,
        dir = :z,
    )

    source_cell = cell_index(g, (source_i, 1, source_k))
    source_x = geo.cell_centroids[1, source_cell]
    source_pressure = pressure0[source_cell]
    source_temperature = tables[:temperature](source_pressure, source_enthalpy)
    source_saturation_vapor = tables[:saturation_vapor_ph](source_pressure, source_enthalpy)
    source_density_liquid = tables[:density_liquid_ph](source_pressure, source_enthalpy)
    source_density_vapor = tables[:density_vapor_ph](source_pressure, source_enthalpy)
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

# ── Private helper ─────────────────────────────────────────────────────────────
#
# Mimics MRST's rampupTimesteps: geometric ramp from a small initial step up
# to dt_target, followed by uniform steps to fill the remaining time.

function _benchmark_rampup_timesteps(
        total_time :: Float64,
        dt_target  :: Float64;
        factor     :: Float64 = 2.0,
        n_ramp     :: Int     = 8,
    )
    dt0  = dt_target / factor^n_ramp
    ramp = [dt0 * factor^i for i in 0:(n_ramp - 1)]
    ramp = min.(ramp, dt_target)

    t_ramp = sum(ramp)
    if t_ramp >= total_time
        # Ramp alone exceeds total time — trim to fit exactly
        cum = cumsum(ramp)
        keep = searchsortedfirst(cum, total_time)
        ramp = ramp[1:keep]
        ramp[end] -= (sum(ramp) - total_time)
        return ramp
    end

    remaining = total_time - t_ramp
    n_steady  = max(1, round(Int, remaining / dt_target))
    dt_steady = remaining / n_steady
    return [ramp; fill(dt_steady, n_steady)]
end
