"""
    Fimbul.build_steam_tables_2ph(; n_pressure = 100, n_enthalpy = 100,
        p_min = 1e4, p_max = 100e6, h_min = 8e4, h_max = 5e6) -> Dict{Symbol, Any}

Build interpolated pure-water steam tables from CoolProp for the two-phase
pressure-enthalpy formulation used by Fimbul.

All returned tables use SI units. Pressure is in Pa, enthalpy in J/kg,
temperature in K, density in kg/m^3, and viscosity in Pa*s.

Returned keys:
- `:temperature`: `(P, H) -> T`
- `:enthalpy`: `(P, T) -> H`
- `:density_mix`: `(P, H) -> value`
- `:enthalpy_liquid_sat`, `:enthalpy_vapor_sat`: `P -> value` saturation-line
    enthalpy tables
- `:density_liquid_ph`, `:density_vapor_ph`: `(P, H) -> value` per-phase tables
- `:viscosity_liquid_ph`, `:viscosity_vapor_ph`: `(P, H) -> value` per-phase tables
- `:enthalpy_liquid_ph`, `:enthalpy_vapor_ph`: `(P, H) -> value` per-phase tables
- `:saturation_vapor_ph`: `(P, H) -> S_v`

The per-phase `(P, H)` tables are defined over the full grid to keep secondary
variable updates branch-free in the OBL-style formulation.
"""
function Fimbul.build_steam_tables_2ph(;
        n_pressure::Int = 100,
        n_enthalpy::Int = 100,
        p_min::Float64  = 1e4,
        p_max::Float64  = 100e6,
        h_min::Float64  = 8e4,
        h_max::Float64  = 5e6,
    )
    p_min > 0     || throw(ArgumentError("p_min must be positive"))
    p_max > p_min || throw(ArgumentError("p_max must be greater than p_min"))

    p = collect(range(p_min, p_max; length = n_pressure))
    p_sat = collect(range(p_min, min(p_max, Fimbul.WATER_CRITICAL_PRESSURE * 0.999); length = n_pressure))

    h_l = [PropsSI("H", "P", p_i, "Q", 0.0, "Water") for p_i in p_sat]
    h_v = [PropsSI("H", "P", p_i, "Q", 1.0, "Water") for p_i in p_sat]
    ρ_l = [PropsSI("D", "P", p_i, "Q", 0.0, "Water") for p_i in p_sat]
    ρ_v = [PropsSI("D", "P", p_i, "Q", 1.0, "Water") for p_i in p_sat]
    μ_l = [PropsSI("V", "P", p_i, "Q", 0.0, "Water") for p_i in p_sat]
    μ_v = [PropsSI("V", "P", p_i, "Q", 1.0, "Water") for p_i in p_sat]

    h_l_subcritical = filter(!isnan, h_l)
    h_v_subcritical = filter(!isnan, h_v)

    if isnan(h_min)
        if isempty(h_l_subcritical)
            h_min = PropsSI("H", "P", Fimbul.WATER_CRITICAL_PRESSURE, "Q", 0.0, "Water")
        else
            h_min = minimum(h_l_subcritical)
        end
    end
    if isnan(h_max)
        if isempty(h_v_subcritical)
            h_max = PropsSI("H", "P", Fimbul.WATER_CRITICAL_PRESSURE, "Q", 1.0, "Water") + 200e3
        else
            h_max = maximum(h_v_subcritical) + 200e3
        end
    end

    h_min < h_max || throw(ArgumentError("h_min must be less than h_max"))

    h = collect(range(h_min, h_max; length = n_enthalpy))

    ρ = [PropsSI("D", "P", p_i, "H", h_j, "Water") for p_i in p, h_j in h]
    μ = [PropsSI("V", "P", p_i, "H", h_j, "Water") for p_i in p, h_j in h]
    T = [PropsSI("T", "P", p_i, "H", h_j, "Water") for p_i in p, h_j in h]

    T_range = collect(range(minimum(T), maximum(T), length = n_enthalpy))
    H = [PropsSI("H", "P", p_i, "T", T_j, "Water") for p_i in p, T_j in T_range]

    make_table_1d(data) = Jutul.get_1d_interpolator(p_sat, data)
    make_table_2d(data) = Jutul.get_2d_interpolator(p, h, data)

    h_l_itp = make_table_1d(h_l)
    h_v_itp = make_table_1d(h_v)
    ρ_l_itp = make_table_1d(ρ_l)
    ρ_v_itp = make_table_1d(ρ_v)
    μ_l_itp = make_table_1d(μ_l)
    μ_v_itp = make_table_1d(μ_v)

    ρ_liq_ph = similar(ρ)
    ρ_vap_ph = similar(ρ)
    μ_liq_ph = similar(μ)
    μ_vap_ph = similar(μ)
    H_liq_ph = similar(ρ)
    H_vap_ph = similar(ρ)
    S_vap_ph = similar(ρ)

    for (i, p_i) in enumerate(p)
        if p_i >= Fimbul.WATER_CRITICAL_PRESSURE
            for (j, h_j) in enumerate(h)
                ρ_liq_ph[i, j] = ρ[i, j]
                ρ_vap_ph[i, j] = ρ[i, j]
                μ_liq_ph[i, j] = μ[i, j]
                μ_vap_ph[i, j] = μ[i, j]
                H_liq_ph[i, j] = h_j
                H_vap_ph[i, j] = h_j
                S_vap_ph[i, j] = 0.0
            end
            continue
        end

        hl = h_l_itp(p_i)
        hv = h_v_itp(p_i)
        rl = ρ_l_itp(p_i)
        rv = ρ_v_itp(p_i)
        ml = μ_l_itp(p_i)
        mv = μ_v_itp(p_i)

        for (j, h_j) in enumerate(h)
            if h_j < hl
                ρ_liq_ph[i, j] = ρ[i, j]
                ρ_vap_ph[i, j] = rv
                μ_liq_ph[i, j] = μ[i, j]
                μ_vap_ph[i, j] = mv
                H_liq_ph[i, j] = h_j
                H_vap_ph[i, j] = hv
                S_vap_ph[i, j] = 0.0
            elseif h_j <= hv
                vapor_quality = (h_j - hl) / (hv - hl)
                @assert 0.0 <= vapor_quality <= 1.0
                vapor_saturation = clamp(
                    (vapor_quality / rv) / (vapor_quality / rv + (1 - vapor_quality) / rl),
                    0.0,
                    1.0,
                )
                ρ_liq_ph[i, j] = rl
                ρ_vap_ph[i, j] = rv
                μ_liq_ph[i, j] = ml
                μ_vap_ph[i, j] = mv
                H_liq_ph[i, j] = hl
                H_vap_ph[i, j] = hv
                S_vap_ph[i, j] = vapor_saturation
            else
                ρ_liq_ph[i, j] = rl
                ρ_vap_ph[i, j] = ρ[i, j]
                μ_liq_ph[i, j] = ml
                μ_vap_ph[i, j] = μ[i, j]
                H_liq_ph[i, j] = hl
                H_vap_ph[i, j] = h_j
                S_vap_ph[i, j] = 1.0
            end
        end
    end

    density_mix = make_table_2d(ρ)
    temperature = make_table_2d(T)
    enthalpy = Jutul.get_2d_interpolator(p, T_range, H)

    return Dict{Symbol, Any}(
        :density_mix         => density_mix,
        :temperature         => temperature,
        :enthalpy            => enthalpy,
        :enthalpy_liquid_sat => h_l_itp,
        :enthalpy_vapor_sat  => h_v_itp,
        :density_liquid_ph   => make_table_2d(ρ_liq_ph),
        :density_vapor_ph    => make_table_2d(ρ_vap_ph),
        :viscosity_liquid_ph => make_table_2d(μ_liq_ph),
        :viscosity_vapor_ph  => make_table_2d(μ_vap_ph),
        :enthalpy_liquid_ph  => make_table_2d(H_liq_ph),
        :enthalpy_vapor_ph   => make_table_2d(H_vap_ph),
        :saturation_vapor_ph => make_table_2d(S_vap_ph),
    )
end
