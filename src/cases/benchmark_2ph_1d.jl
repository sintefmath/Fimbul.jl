"""
    benchmark_1d_geothermal(; <keyword arguments>)

1D two-phase geothermal benchmark cases from Weis et al. (2014), translated
from the MRST `benchmark_1d_geothermal` test suite and implemented using
`GeothermalTwoPhaseSystem` with IAPWS-IF97 steam tables.

Five cases (`:a`–`:e`) cover different initial pressure/temperature conditions
and pressure gradients. Each case is a 1D column with Dirichlet-like
pressure/enthalpy conditions at both ends, imposed via BHP-controlled simple
injector and producer wells.

| Case | p_inj [MPa] | p_prod [MPa] | T_inj [°C] | T_prod [°C] | gravity? |
|:----:|:-----------:|:------------:|:----------:|:-----------:|:--------:|
| :a   | 50          | 25           | 350        | 150         | optional |
| :b   | 40          | 20           | 450        | 300         | optional |
| :c   | 15          | 1            | 500        | 350         | optional |
| :d   | 20          | 1            | 400        | 150         | optional |
| :e   | 4           | 1            | 300        | 150         | no       |

## Grid orientation

- `gravity = false`: flow in x, grid `(nx, 1, 1)`. Injector at cell 1 (left),
  producer at cell nx (right).
- `gravity = true`:  flow in z, grid `(1, 1, nx)`. Injector at cell nx
  (bottom = deep = high pressure), producer at cell 1 (top = shallow).

## Keyword arguments

- `benchmark_case = :a`: Case selector, one of `:a`, `:b`, `:c`, `:d`, `:e`.
- `gravity = false`: Enable gravity-driven flow (not available for case `:e`).
- `nx = 200`: Number of cells along the flow direction (default matches MRST
  with `length = 2000 m` and `dx = 10 m`).
- `domain_length = 2000.0*meter`: Length of the 1D domain [m].
- `cell_size = 10.0*meter`: Cell width in the two inactive directions [m].
- `dt_target = 5*year`: Target timestep size [s] after ramp-up.
- `num_years = 1500`: Total simulation duration [years] (matches MRST default;
  the "interesting" observation times per case are a/b/d: shorter, c/e: 1500 yr).
- `enthalpy_tables = build_steam_tables_2ph()`: Pre-built (P,H) steam tables;
  pass your own to avoid rebuilding when calling the function multiple times.

## Returns

A `JutulCase` ready for `simulate_reservoir`.

## Reference

Weis, P., Driesner, T., & Heinrich, C. A. (2014). Hydrothermal, multiphase
convection of H₂O-NaCl fluids from ambient to magmatic temperatures: A new
numerical scheme and benchmarks for code comparison. *Geofluids*, 14(3),
347–371. https://doi.org/10.1111/gfl.12080
"""
function benchmark_2ph_1d(;
    benchmark_case  :: Symbol  = :a,
    gravity         :: Bool    = false,
    nx              :: Int     = 200,
    domain_length   :: Float64 = 2000.0*meter,
    cell_size       :: Float64 = 10.0*meter,
    dt_target       :: Float64 = 5.0*year,
    num_years       :: Int     = 1500,
    enthalpy_tables                = build_steam_tables_2ph(),
)

    benchmark_case in (:a, :b, :c, :d, :e) ||
        throw(ArgumentError("benchmark_case must be one of :a, :b, :c, :d, :e"))
    (benchmark_case == :e && gravity) &&
        throw(ArgumentError("Case :e is only defined without gravity (gravity = false)"))

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
    else  # :e
        p_inj  = 4e6;   p_prod = 1e6
        T_inj  = K0 + 300;  T_prod = K0 + 150
    end

    H_inj  = enthalpy_tables[:enthalpy](p_inj, T_inj)
    H_prod = enthalpy_tables[:enthalpy](p_prod, T_prod)
    # H_prod = enthalpy_tables[:H_pT](p_prod, T_prod)

    # ── Grid ──────────────────────────────────────────────────────────────
    #
    # No gravity: horizontal 1D flow in x  →  (nx, 1, 1)
    #   cell 1  = left  (inlet  / high p, high T)
    #   cell nx = right (outlet / low  p, low  T)
    #
    # With gravity: vertical 1D flow in z  →  (1, 1, nx)
    #   cell 1  = top    (shallow / low  p, low  T)  ← producer
    #   cell nx = bottom (deep    / high p, high T)  ← injector
    #
    if gravity
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
    h0 = maximum(enthalpy_tables[:enthalpy].(p0, T_prod))
    # h0 = enthalpy_tables[:enthalpy].(p0, T_prod)

    domain = reservoir_domain(g;
        permeability              = 1e-15,   # m²  (≈ 10 μD)
        porosity                  = 0.1,
        rock_density              = 2700.0,  # kg/m³
        rock_heat_capacity        = 880.0,   # J/(kg·K)
        rock_thermal_conductivity = 2.0,     # W/(m·K)
    )

    # ── Wells ──────────────────────────────────────────────────────────────
    WI = 1e-8
    WIth = 0.0
    well_inj  = setup_well(domain, [cell_inj];  name = :Injector, WI = WI, WIth = WIth, simple_well = false, use_top_node = true)
    well_prod = setup_well(domain, [cell_prod]; name = :Producer, WI = WI, WIth = WIth, simple_well = false, use_top_node = true)

    # ── Model ──────────────────────────────────────────────────────────────
    model = setup_reservoir_model_geothermal_2ph(
        domain;
        enthalpy_tables = enthalpy_tables,
        wells           = [well_inj, well_prod],
        block_backend=false,
    )

    # Add function handle for omputing enthalpy from (P, T)
    for k in keys(model.models)
        model.models[k].extra[:enthalpy] = enthalpy_tables[:enthalpy]
        push!(model.models[k].output_variables, :PhaseViscosities)
    end
    # H_prod = 1e6
    # H_inj  = 2e6

    # ── Initial state ──────────────────────────────────────────────────────
    # Uniform initial conditions at the producer-side state (low p, low T).
    # The MRST benchmark also initialises uniformly (not hydrostatic).
    state0 = setup_reservoir_state(model;
        Pressure = p0,
        Enthalpy = h0,
        Temperature = T_prod,
    )

    # ── Well controls ──────────────────────────────────────────────────────
    # Both ends are driven by fixed bottom-hole pressure (BHP), mimicking
    # the Dirichlet pressure / temperature boundary conditions in MRST.
    # rhoL_ref = first(JutulDarcy.reference_densities(reservoir_model(model).system))
    rho_inj = enthalpy_tables[:density_mix](p_inj, H_inj)
    rate = 1si"meter^3/day"
    ctrl_inj = InjectorControl(
        BottomHolePressureTarget(p_inj),
        # TotalRateTarget(rate),
        [1.0],
        density     = rho_inj,
        temperature = T_inj,
        enthalpy    = H_inj,
    )
    ctrl_prod = ProducerControl(BottomHolePressureTarget(p_prod))

    forces = setup_reservoir_forces(model;
        control = Dict(:Injector => ctrl_inj, :Producer => ctrl_prod),
    )
    forces = with_property_evaluators(model, forces)  # add (P,H)-dependent properties

    # ── Timesteps ──────────────────────────────────────────────────────────
    total_time = Float64(num_years) * year
    dt_vec     = _benchmark_rampup_timesteps(total_time, dt_target)

    if benchmark_case == :a
        plot_time = ifelse(!gravity, 250.0*year, 750.0*year)
    elseif benchmark_case == :b
        plot_time = ifelse(!gravity, 120.0*year, 350.0*year)
    elseif benchmark_case == :c
        plot_time = 1500.0*year
    elseif benchmark_case == :d
        plot_time = ifelse(!gravity, 200.0*year, 1000.0*year)
    else  # :e
        plot_time = 2000.0*year
    end

    # ── Assemble case ──────────────────────────────────────────────────────
    info = Dict{Symbol, Any}(
        :description    => "1D geothermal benchmark case $(benchmark_case) (Weis et al. 2014)",
        :benchmark_case => benchmark_case,
        :plot_time      => plot_time,
        :gravity        => gravity,
        :p_inj          => p_inj,
        :p_prod         => p_prod,
        :T_inj          => T_inj,
        :T_prod         => T_prod,
    )
    return JutulCase(model, dt_vec, forces; state0 = state0, input_data = info)
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
