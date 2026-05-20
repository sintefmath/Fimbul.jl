const _PCRIT_WATER_CP = 22.064e6   # critical pressure of water [Pa] (IAPWS-IF97 / CoolProp)

function Fimbul.build_steam_tables_2ph(;
        n_pressure::Int = 100,
        n_enthalpy::Int = 100,
        p_min::Float64  = 1e4,     # Pa  (0.1 bar)
        p_max::Float64  = 100e6,   # Pa  (1000 bar)
        h_min::Float64  = 8e4,     # J/kg, auto-computed if not provided
        h_max::Float64  = 5e6,     # J/kg, auto-computed if not provided
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

    make_table_1d(data) = Jutul.get_1d_interpolator(p_sub, data)
    make_table_2d(data) = Jutul.get_2d_interpolator(p, h, data)

    # ── Per-phase 2D lookup tables (OBL approach) ───────────────────────────
    # Each table is defined continuously over the full (P, H) grid so that all
    # secondary-variable updates are branch-free and AD-compatible.
    #
    # Fill rules per region:
    #   Supercritical (P ≥ Pcrit): liquid = vapor = single-phase value
    #   Subcooled (H < h_l(P)):    liquid = actual (P,H); vapor = sat-vapor ref
    #   Two-phase (h_l ≤ H ≤ h_v): liquid = sat-liquid; vapor = sat-vapor
    #   Superheated (H > h_v(P)):  liquid = sat-liquid ref; vapor = actual (P,H)

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
        if p_i >= _PCRIT_WATER_CP
            # Supercritical: single-phase fluid; both phase slots identical
            for (j, h_j) in enumerate(h)
                ρ_liq_ph[i, j] = ρ[i, j]
                ρ_vap_ph[i, j] = ρ[i, j]
                μ_liq_ph[i, j] = μ[i, j]
                μ_vap_ph[i, j] = μ[i, j]
                H_liq_ph[i, j] = h_j
                H_vap_ph[i, j] = h_j
                S_vap_ph[i, j] = 0.0
            end
        else
            hl = h_l_itp(p_i)
            hv = h_v_itp(p_i)
            rl = ρ_l_itp(p_i)
            rv = ρ_v_itp(p_i)
            ml = μ_l_itp(p_i)
            mv = μ_v_itp(p_i)
            for (j, h_j) in enumerate(h)
                if h_j < hl
                    # Subcooled liquid
                    ρ_liq_ph[i, j] = ρ[i, j];  ρ_vap_ph[i, j] = rv
                    μ_liq_ph[i, j] = μ[i, j];  μ_vap_ph[i, j] = mv
                    H_liq_ph[i, j] = h_j;       H_vap_ph[i, j] = hv
                    S_vap_ph[i, j] = 0.0
                elseif h_j <= hv
                    # Two-phase: saturation values; volume saturation from lever rule
                    x = (h_j - hl) / (hv - hl)
                    @assert 0.0 <= x <= 1.0
                    sv = clamp((x / rv) / (x / rv + (1 - x) / rl), 0.0, 1.0)
                    ρ_liq_ph[i, j] = rl;  ρ_vap_ph[i, j] = rv
                    μ_liq_ph[i, j] = ml;  μ_vap_ph[i, j] = mv
                    H_liq_ph[i, j] = hl;  H_vap_ph[i, j] = hv
                    S_vap_ph[i, j] = sv
                else
                    # Superheated vapor
                    ρ_liq_ph[i, j] = rl;       ρ_vap_ph[i, j] = ρ[i, j]
                    μ_liq_ph[i, j] = ml;       μ_vap_ph[i, j] = μ[i, j]
                    H_liq_ph[i, j] = hl;       H_vap_ph[i, j] = h_j
                    S_vap_ph[i, j] = 1.0
                end
            end
        end
    end

    return Dict{Symbol, Any}(
        :density_liquid      => make_table_1d(ρ_l),
        :density_vapor       => make_table_1d(ρ_v),
        :density_mix         => make_table_2d(ρ),
        :viscosity_liquid    => make_table_1d(μ_l),
        :viscosity_vapor     => make_table_1d(μ_v),
        :viscosity_mix       => make_table_2d(μ),
        :enthalpy_liquid     => make_table_1d(h_l),
        :enthalpy_vapor      => make_table_1d(h_v),
        :temperature         => make_table_2d(T),
        :enthalpy            => Jutul.get_2d_interpolator(p, T_range, H),
        :enthalpy_liquid_sat => make_table_1d(h_l),
        :enthalpy_vapor_sat  => make_table_1d(h_v),
        # Per-phase 2D tables (OBL)
        :density_liquid_ph   => make_table_2d(ρ_liq_ph),
        :density_vapor_ph    => make_table_2d(ρ_vap_ph),
        :viscosity_liquid_ph => make_table_2d(μ_liq_ph),
        :viscosity_vapor_ph  => make_table_2d(μ_vap_ph),
        :enthalpy_liquid_ph  => make_table_2d(H_liq_ph),
        :enthalpy_vapor_ph   => make_table_2d(H_vap_ph),
        :saturation_vapor_ph => make_table_2d(S_vap_ph),
    )
end
