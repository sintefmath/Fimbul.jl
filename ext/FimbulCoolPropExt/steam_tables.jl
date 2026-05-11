const _PCRIT_WATER_CP = 22.064e6   # critical pressure of water [Pa] (IAPWS-IF97 / CoolProp)

function Fimbul.build_steam_tables_2ph(;
        n_pressure::Int = 100,
        n_enthalpy::Int = 100,
        p_min::Float64  = 1e5,     # Pa  (1 bar)
        p_max::Float64  = 100e6,   # Pa  (1000 bar)
        h_min::Float64  = NaN,     # J/kg, auto-computed if not provided
        h_max::Float64  = NaN,     # J/kg, auto-computed if not provided
    )
    p_min > 0     || throw(ArgumentError("p_min must be positive"))
    p_max > p_min || throw(ArgumentError("p_max must be greater than p_min"))

    # Pressure axis [Pa]
    p = collect(range(p_min, p_max; length = n_pressure))

    p_sub = collect(range(p_min, min(p_max, _PCRIT_WATER_CP * 0.999); length = n_pressure))

    # Saturation properties over subcritical pressures
    h_l = [PropsSI("H", "P", p_i, "Q", 0.0, "Water") for p_i in p_sub]
    h_v = [PropsSI("H", "P", p_i, "Q", 1.0, "Water") for p_i in p_sub]
    ρ_l = [PropsSI("D", "P", p_i, "Q", 0.0, "Water") for p_i in p_sub]
    ρ_v = [PropsSI("D", "P", p_i, "Q", 1.0, "Water") for p_i in p_sub]
    μ_l = [PropsSI("V", "P", p_i, "Q", 0.0, "Water") for p_i in p_sub]
    μ_v = [PropsSI("V", "P", p_i, "Q", 1.0, "Water") for p_i in p_sub]

    # Collect subcritical saturation enthalpies for auto-range computation
    hL_subcrit = filter(!isnan, h_l)
    hV_subcrit = filter(!isnan, h_v)

    # Enthalpy axis defaults
    if isnan(h_min)
        if isempty(hL_subcrit)
            h_min = PropsSI("H", "P", _PCRIT_WATER_CP, "Q", 0.0, "Water")
        else
            h_min = minimum(hL_subcrit)
        end
    end
    if isnan(h_max)
        if isempty(hV_subcrit)
            h_max = PropsSI("H", "P", _PCRIT_WATER_CP, "Q", 1.0, "Water") + 200e3
        else
            h_max = maximum(hV_subcrit) + 200e3
        end
    end

    h_min < h_max || throw(ArgumentError("h_min must be less than h_max"))

    h = collect(range(h_min, h_max; length = n_enthalpy))  # J/kg axis

    ρ = [PropsSI("D", "P", p_i, "H", h_j, "Water") for p_i in p, h_j in h]
    μ = [PropsSI("V", "P", p_i, "H", h_j, "Water") for p_i in p, h_j in h]
    T = [PropsSI("T", "P", p_i, "H", h_j, "Water") for p_i in p, h_j in h]

    T_range = collect(range(minimum(T), maximum(T), length = n_enthalpy))
    H = [PropsSI("H", "P", p_i, "T", T_j, "Water") for p_i in p, T_j in T_range]

    make_table_1d(data) = Jutul.get_1d_interpolator(p_sub, data; cap_endpoints = true)
    make_table_2d(data) = Jutul.get_2d_interpolator(p, h, data)

    return Dict{Symbol, Any}(
        :density_liquid   => make_table_1d(ρ_l),
        :density_vapor    => make_table_1d(ρ_v),
        :density_mix      => make_table_2d(ρ),
        :viscosity_liquid => make_table_1d(μ_l),
        :viscosity_vapor  => make_table_1d(μ_v),
        :viscosity_mix    => make_table_2d(μ),
        :enthalpy_liquid  => make_table_1d(h_l),
        :enthalpy_vapor   => make_table_1d(h_v),
        :temperature      => make_table_2d(T),
        :enthalpy         => Jutul.get_2d_interpolator(p, T_range, H),
    )
end
