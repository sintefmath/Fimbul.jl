"""
    Fimbul.build_steam_tables_h2o(; n_pressure = 100, n_enthalpy = 100,
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
function Fimbul.build_steam_tables_h2o(;
        n_pressure::Int = 100,
        n_enthalpy::Int = 100,
        p_min::Float64  = 1e4,
        p_max::Float64  = 100e6,
        h_min::Float64  = 1e3,
        h_max::Float64  = 5e6,
        info_level::Int = 0
    )
    # Check input validity
    p_min > 0     || throw(ArgumentError("p_min must be positive"))
    p_max > p_min || throw(ArgumentError("p_max must be greater than p_min"))
    h_min < h_max || throw(ArgumentError("h_min must be less than h_max"))

    # Saturation pressures relies on temperature, so we use this to sample the saturation envelope
    T_min = PropsSI("T", "P", p_min, "H", h_min, "Water")
    T_sat_max = PropsSI("T", "P", Fimbul.WATER_CRITICAL_PRESSURE, "H", Fimbul.WATER_CRITICAL_ENTHALPY, "Water")
    # Sample saturation and temperature pressures naively from min to max
    p_sat = collect(range(p_min, min(p_max, Fimbul.WATER_CRITICAL_PRESSURE * 0.999); length = n_pressure))
    T_sat = collect(range(T_min, T_sat_max; length = n_enthalpy))
    # Identify missing saturation pressures and temperatures along the line
    p_l_missing = sample_property("P", "T", T_sat, "Q", [0.0], info_level)
    p_v_missing = sample_property("P", "T", T_sat, "Q", [1.0], info_level)
    T_l_missing = sample_property("T", "P", p_sat, "Q", [0.0], info_level)
    T_v_missing = sample_property("T", "P", p_sat, "Q", [1.0], info_level)
    # Add missing saturation pressures and temperatures to the sampling points
    p_sat = sort(unique(vcat(p_sat, p_l_missing, p_v_missing)))
    T_sat = sort(unique(vcat(T_sat, T_l_missing, T_v_missing)))
    # Sample saturation enthalpies at the saturation pressures
    h_sat = sample_property("H", "P", p_sat, "Q", 0.0, info_level)
    # Sample saturation properties at the saturation pressures
    info_level > 0 && @info "Building saturation tables with \
    $n_pressure pressure points ($minimum(p_sat), $maximum(p_sat)) Pa"
    h_l = [PropsSI("H", "P", p_i, "Q", 0.0, "Water") for p_i in p_sat]
    h_v = [PropsSI("H", "P", p_i, "Q", 1.0, "Water") for p_i in p_sat]
    ρ_l = [PropsSI("D", "P", p_i, "Q", 0.0, "Water") for p_i in p_sat]
    ρ_v = [PropsSI("D", "P", p_i, "Q", 1.0, "Water") for p_i in p_sat]
    μ_l = [PropsSI("V", "P", p_i, "Q", 0.0, "Water") for p_i in p_sat]
    μ_v = [PropsSI("V", "P", p_i, "Q", 1.0, "Water") for p_i in p_sat]
    # Make 1D tables for saturation properties
    make_table_1d(data) = Jutul.get_1d_interpolator(p_sat, data; cap_endpoints = true)
    h_l_itp = make_table_1d(h_l)
    h_v_itp = make_table_1d(h_v)
    ρ_l_itp = make_table_1d(ρ_l)
    ρ_v_itp = make_table_1d(ρ_v)
    μ_l_itp = make_table_1d(μ_l)
    μ_v_itp = make_table_1d(μ_v)

    # Complete the grid with single-phase points if the max pressure and enthalpy exceed the critical point
    if p_max > Fimbul.WATER_CRITICAL_PRESSURE
        p_sc = collect(range(Fimbul.WATER_CRITICAL_PRESSURE, p_max; length = n_pressure))
        p = sort(unique(vcat(p_sat, p_sc)))
    end
    if h_max > Fimbul.WATER_CRITICAL_ENTHALPY
        h_sc = collect(range(Fimbul.WATER_CRITICAL_ENTHALPY, h_max; length = n_enthalpy))
        h = sort(unique(vcat(h_sat, h_sc)))
    end
    T_max = PropsSI("T", "P", p_max, "H", h_max, "Water")
    if T_max > Fimbul.WATER_CRITICAL_TEMPERATURE
        # T_sc = collect(range(Fimbul.WATER_CRITICAL_TEMPERATURE, T_max; length = n_enthalpy))
        n_temperature = length(h) - length(T_sat)
        T_sc = collect(range(Fimbul.WATER_CRITICAL_TEMPERATURE, T_max; length = n_temperature))
        T = sort(unique(vcat(T_sat, T_sc)))
    end
    info_level > 0 && @info "Building two-phase tables with \
    $(length(p)) pressure points ($(minimum(p)), $(maximum(p))) Pa and \
    $(length(h)) enthalpy points ($(minimum(h)), $(maximum(h))) J/kg."
    # Sample properties over the full grid
    ρ_tab = sample_property("D", "P", p, "H", h, info_level)
    μ_tab = sample_property("V", "P", p, "H", h, info_level)
    T_tab = sample_property("T", "P", p, "H", h, info_level)
    h_tab = sample_property("H", "P", p, "T", T, info_level)

    # Construct per-phase tables by identifying the phase at each point and
    # assigning the corresponding properties.
    ρ_liq_ph = similar(ρ_tab)
    ρ_vap_ph = similar(ρ_tab)
    μ_liq_ph = similar(μ_tab)
    μ_vap_ph = similar(μ_tab)
    h_liq_ph = similar(h_tab)
    h_vap_ph = similar(h_tab)
    S_vap_ph = similar(h_tab)

    for (i, p_i) in enumerate(p)
        if p_i >= Fimbul.WATER_CRITICAL_PRESSURE
            for (j, h_j) in enumerate(h)
                ρ_liq_ph[i, j] = ρ_tab[i, j]
                ρ_vap_ph[i, j] = ρ_tab[i, j]
                μ_liq_ph[i, j] = μ_tab[i, j]
                μ_vap_ph[i, j] = μ_tab[i, j]
                h_liq_ph[i, j] = h_j
                h_vap_ph[i, j] = h_j
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
                ρ_liq_ph[i, j] = ρ_tab[i, j]
                ρ_vap_ph[i, j] = rv
                μ_liq_ph[i, j] = μ_tab[i, j]
                μ_vap_ph[i, j] = mv
                h_liq_ph[i, j] = h_j
                h_vap_ph[i, j] = hv
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
                h_liq_ph[i, j] = hl
                h_vap_ph[i, j] = hv
                S_vap_ph[i, j] = vapor_saturation
            else
                ρ_liq_ph[i, j] = rl
                ρ_vap_ph[i, j] = ρ_tab[i, j]
                μ_liq_ph[i, j] = ml
                μ_vap_ph[i, j] = μ_tab[i, j]
                h_liq_ph[i, j] = hl
                h_vap_ph[i, j] = h_j
                S_vap_ph[i, j] = 1.0
            end
        end
    end

    make_table_2d(p, h, T_tab)
    make_table_2d(p, T, h_tab)
    h_l_itp
    h_v_itp
    make_table_2d(p, h, ρ_liq_ph)
    make_table_2d(p, h, ρ_vap_ph)
    make_table_2d(p, h, ρ_tab)
    make_table_2d(p, h, μ_liq_ph)
    make_table_2d(p, h, μ_vap_ph)
    make_table_2d(p, h, μ_tab)
    make_table_2d(p, h, h_liq_ph)
    make_table_2d(p, h, h_vap_ph)
    make_table_2d(p, h, S_vap_ph)

    return Dict{Symbol, Any}(
        :temperature         => make_table_2d(p, h, T_tab),
        :enthalpy            => make_table_2d(p, T, h_tab),
        :enthalpy_liquid_sat => h_l_itp,
        :enthalpy_vapor_sat  => h_v_itp,
        :density_liquid_ph   => make_table_2d(p, h, ρ_liq_ph),
        :density_vapor_ph    => make_table_2d(p, h, ρ_vap_ph),
        :density_mix         => make_table_2d(p, h, ρ_tab),
        :viscosity_liquid_ph => make_table_2d(p, h, μ_liq_ph),
        :viscosity_vapor_ph  => make_table_2d(p, h, μ_vap_ph),
        :viscosity_mix       => make_table_2d(p, h, μ_tab),
        :enthalpy_liquid_ph  => make_table_2d(p, h, h_liq_ph),
        :enthalpy_vapor_ph   => make_table_2d(p, h, h_vap_ph),
        :saturation_vapor_ph => make_table_2d(p, h, S_vap_ph),
    )
end

function sample_property(property, x_name, x, y_name, y, info_level::Int = 0)

    v = zeros(Float64, length(x), length(y))
    for (i, xi) in enumerate(x)
        for (j, yj) in enumerate(y)
            v_ij = NaN
            try
                v_ij = PropsSI(property, x_name, xi, y_name, yj, "Water")
            catch e
                info_level > 1 && @warn "PropsSI error for $property \
                $x_name = $xi, $y_name = $yj: $(e.msg). \
                Setting property to NaN." exception = e
            end
            v[i, j] = v_ij
        end
    end
    # Replace NaNs with nearest valid value in the same column
    for j in eachindex(y)
        column = view(v, :, j)
        @assert !all(isnan, column) "$property: All values are NaN at \
        $y_name = $(y[j]). This cannot generate table."
        for i in eachindex(x)
            if isnan(column[i])
                # Find nearest non-NaN value in the same column
                valid_indices = findall(!isnan, column)
                nearest_index = valid_indices[argmin(abs.(valid_indices .- i))]
                column[i] = column[nearest_index]
            end
        end
    end
    return v

end

function make_table_2d(x, y, v)
    # x, y, v = pad_table(x, y, v)
    return Jutul.get_2d_interpolator(x, y, v; cap_endpoints = true)
end

function pad_table(x, y, v, ϵ_xy = 1.0, ϵ_v=1e-3)

    x_min, x_max = minimum(x), maximum(x)
    y_min, y_max = minimum(y), maximum(y)

    Δx = x_max - x_min
    Δy = y_max - y_min
    x_padded = vcat(x_min - Δx * ϵ_xy, x, x_max + Δx * ϵ_xy)
    y_padded = vcat(y_min - Δy * ϵ_xy, y, y_max + Δy * ϵ_xy)

    v_padded = zeros(size(v) .+ 2)
    v_padded[2:end-1, 2:end-1] .= v
    v_padded[1, :] .= v_padded[2, :].* (1 - ϵ_v)
    v_padded[end, :] .= v_padded[end-1, :].* (1 + ϵ_v)
    v_padded[:, 1] .= v_padded[:, 2].* (1 - ϵ_v)
    v_padded[:, end] .= v_padded[:, end-1].* (1 + ϵ_v)

    return x_padded, y_padded, v_padded

end