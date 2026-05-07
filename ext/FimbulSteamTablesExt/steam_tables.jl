# (P, H) property tables for a two-phase (liquid + vapour) water system.
#
# SteamTables (IAPWS-IF97) uses P in MPa and H in kJ/kg.
# Fimbul/Jutul uses SI: P in Pa, H in J/kg, T in K, rho in kg/m³, mu in Pa·s,
# c_p in J/(kg·K).
#
# Region-3 of IAPWS-IF97 (P ≥ 16.529 MPa) is avoided because SteamTables.jl
# has a scoping bug in Region3_TPh / Region3_VPh. Above this threshold the
# code falls back to a safe (P, T) sweep via `SpecificH`, `SpecificV`,
# `SpecificCP`.

const _Pc_MPa = 22.064   # water critical pressure [MPa]

# ── Viscosity helpers ─────────────────────────────────────────────────────────

"""
    _water_viscosity(T) -> μ [Pa·s]

Vogel-Fulcher-Tammann (VFT / Andrade) viscosity approximation for liquid water.
`T` in K. Accurate to ~2 % for 25–300 °C, pressure-independent.
"""
_water_viscosity(T::Real) = 2.414e-5 * 10^(247.8 / (T - 140.0))

"""
    _steam_viscosity(T) -> μ [Pa·s]

Linear approximation for steam dynamic viscosity. `T` in K.
Fit to IAPWS data over 100–400 °C: μ ≈ −3.4×10⁻⁶ + 4.3×10⁻⁸ T.
Clipped to 8×10⁻⁶ Pa·s as a lower bound.
"""
_steam_viscosity(T::Real) = max(8e-6, -3.4e-6 + 4.3e-8 * T)

# ── Two-phase (liquid + vapour) (P, H) tables ─────────────────────────────────

"""
    build_ph_tables_2ph(; P_min, P_max, H_min, H_max, T_min, T_max,
                           nP, nH, nT) -> Dict{Symbol, Any}

Build bilinear (P, H) lookup tables for a two-phase (liquid + vapour) pure-
water system, suitable for use with `setup_reservoir_model_geothermal_2ph`.

Requires `SteamTables.jl` (loaded via the `FimbulSteamTablesExt` extension).

# Returned `Dict` callables (all SI: P [Pa], H [J/kg])
- `:T`        — `(P, H) → T [K]`
- `:rho`      — `(P, H) → (ρ_l, ρ_v) [kg/m³]`
- `:mu`       — `(P, H) → (μ_l, μ_v) [Pa·s]`
- `:S`        — `(P, H) → (S_l, S_v)` volume saturations (sum = 1)
- `:H_phases` — `(P, H) → (H_l, H_v) [J/kg]` per-phase specific enthalpy
- `:c_p`      — `(P, H) → (cp_l, cp_v) [J/(kg·K)]`
- `:H_pT`     — `(P, T) → H [J/kg]` liquid enthalpy (for T → H initialisation)

# Regions
- Sub-cooled liquid  (H < H_l_sat):  S_l = 1, S_v = 0
- Two-phase (H_l_sat ≤ H ≤ H_v_sat): S_v from vapour quality x = (H-Hf)/(Hg-Hf),
  then volume fraction S_v = (x/ρ_v) / (x/ρ_v + (1-x)/ρ_l)
- Superheated steam  (H > H_v_sat):  S_l = 0, S_v = 1
- Supercritical / high-P (P ≥ 16.529 MPa): single-phase via safe (P,T) sweep

# Keyword arguments
- `P_min` / `P_max`: pressure range [Pa], default 1 bar – 1000 bar
- `H_min` / `H_max`: specific enthalpy range [J/kg], default 50–3000 kJ/kg
- `T_min` / `T_max`: temperature range for H(P,T) table [K]
- `nP`, `nH`, `nT`: grid sizes
"""
function Fimbul.build_steam_tables_2ph(;
        P_min = 1e5,      # Pa  (1 bar)
        P_max = 100e6,    # Pa  (1000 bar)
        H_min = 50e3,     # J/kg
        H_max = 3000e3,   # J/kg  (superheated steam up to ~540 °C at 1 bar)
        T_min = 274.0,    # K
        T_max = 623.15,   # K   (350 °C)
        nP    = 80,
        nH    = 100,
        nT    = 60,
    )
    P_grid = collect(range(P_min, P_max; length = nP))   # [Pa]
    H_grid = collect(range(H_min, H_max; length = nH))   # [J/kg]

    T_arr    = Matrix{Float64}(undef, nP, nH)
    rhoL_arr = Matrix{Float64}(undef, nP, nH)
    rhoV_arr = Matrix{Float64}(undef, nP, nH)
    muL_arr  = Matrix{Float64}(undef, nP, nH)
    muV_arr  = Matrix{Float64}(undef, nP, nH)
    SL_arr   = Matrix{Float64}(undef, nP, nH)
    SV_arr   = Matrix{Float64}(undef, nP, nH)
    HL_arr   = Matrix{Float64}(undef, nP, nH)
    HV_arr   = Matrix{Float64}(undef, nP, nH)
    cpL_arr  = Matrix{Float64}(undef, nP, nH)
    cpV_arr  = Matrix{Float64}(undef, nP, nH)

    for (i, P) in enumerate(P_grid)
        P_mpa = P / 1e6

        # Avoid Region 3 of IAPWS-IF97 (P ≥ 16.529 MPa, T ≥ 623.15 K) where
        # SteamTables.jl has a scoping bug in Region3_TPh / Region3_VPh.
        subcrit = P_mpa < 16.529

        if subcrit
            T_sat    = SteamTables.Tsat(P_mpa)           # K
            Hf       = SteamTables.SatHL(T_sat) * 1e3    # J/kg, sat-liquid enthalpy
            Hg       = SteamTables.SatHV(T_sat) * 1e3    # J/kg, sat-vapour enthalpy
            rhoL_sat = SteamTables.SatDensL(T_sat)        # kg/m³
            rhoV_sat = SteamTables.SatDensV(T_sat)        # kg/m³
            muL_sat  = _water_viscosity(T_sat)
            muV_sat  = _steam_viscosity(T_sat)
            # cp diverges exactly at saturation near P_c; evaluate 0.5 K inside
            # each phase and cap at 50 kJ/(kg·K).
            T_liq = T_sat - 0.5
            T_vap = min(T_sat + 0.5, 646.0)   # stay below T_c = 647.096 K
            cpL_sat = min(SteamTables.SpecificCP(P_mpa, T_liq) * 1e3, 50e3)
            cpV_sat = min(SteamTables.SpecificCP(P_mpa, T_vap) * 1e3, 50e3)

            # Not used in the subcritical path
            T_sweep_K    = nothing
            H_sweep_kjkg = nothing
            rho_sweep    = nothing
            cp_sweep     = nothing
        else
            # High-pressure / supercritical: pre-build property sweeps over T
            # using the safe (P, T) path to avoid Region-3 failures.
            T_sat = NaN; Hf = Inf; Hg = -Inf
            rhoL_sat = rhoV_sat = NaN
            muL_sat = muV_sat = NaN
            cpL_sat = cpV_sat = NaN

            T_sweep_K    = collect(range(274.0, 750.0; length = 300))
            H_sweep_kjkg = [SteamTables.SpecificH(P_mpa, T) for T in T_sweep_K]
            rho_sweep    = [1.0 / SteamTables.SpecificV(P_mpa, T) for T in T_sweep_K]
            cp_sweep     = [min(SteamTables.SpecificCP(P_mpa, T) * 1e3, 50e3)
                            for T in T_sweep_K]
        end

        for (j, H) in enumerate(H_grid)
            H_kjkg = H / 1e3

            if subcrit && H < Hf
                # ── Sub-cooled liquid ──────────────────────────────────────
                T     = SteamTables.Temperature_Ph(P_mpa, H_kjkg)
                rho_l = 1.0 / SteamTables.SpecificV_Ph(P_mpa, H_kjkg)
                cp_l  = SteamTables.SpecificCP_Ph(P_mpa, H_kjkg) * 1e3
                T_arr[i, j]    = T
                rhoL_arr[i, j] = rho_l
                rhoV_arr[i, j] = rhoV_sat   # reference (S_v = 0)
                muL_arr[i, j]  = _water_viscosity(T)
                muV_arr[i, j]  = muV_sat
                SL_arr[i, j]   = 1.0
                SV_arr[i, j]   = 0.0
                HL_arr[i, j]   = H
                HV_arr[i, j]   = Hg
                cpL_arr[i, j]  = cp_l
                cpV_arr[i, j]  = cpV_sat

            elseif subcrit && H <= Hg
                # ── Two-phase (liquid + vapour) ────────────────────────────
                x  = clamp((H - Hf) / (Hg - Hf), 0.0, 1.0)   # vapour quality
                # Volume fraction of vapour: S_v = (x/ρ_v) / (x/ρ_v + (1-x)/ρ_l)
                sv = (x / rhoV_sat) / (x / rhoV_sat + (1.0 - x) / rhoL_sat)
                T_arr[i, j]    = T_sat
                rhoL_arr[i, j] = rhoL_sat
                rhoV_arr[i, j] = rhoV_sat
                muL_arr[i, j]  = muL_sat
                muV_arr[i, j]  = muV_sat
                SL_arr[i, j]   = 1.0 - sv
                SV_arr[i, j]   = sv
                HL_arr[i, j]   = Hf
                HV_arr[i, j]   = Hg
                cpL_arr[i, j]  = cpL_sat
                cpV_arr[i, j]  = cpV_sat

            elseif subcrit && H > Hg
                # ── Superheated steam ──────────────────────────────────────
                T     = SteamTables.Temperature_Ph(P_mpa, H_kjkg)
                rho_v = 1.0 / SteamTables.SpecificV_Ph(P_mpa, H_kjkg)
                cp_v  = SteamTables.SpecificCP_Ph(P_mpa, H_kjkg) * 1e3
                T_arr[i, j]    = T
                rhoL_arr[i, j] = rhoL_sat   # reference (S_l = 0)
                rhoV_arr[i, j] = rho_v
                muL_arr[i, j]  = muL_sat
                muV_arr[i, j]  = _steam_viscosity(T)
                SL_arr[i, j]   = 0.0
                SV_arr[i, j]   = 1.0
                HL_arr[i, j]   = Hf
                HV_arr[i, j]   = H
                cpL_arr[i, j]  = cpL_sat
                cpV_arr[i, j]  = cp_v

            else
                # ── Supercritical / high-P: T-sweep inversion via (P,T) path ──
                idx  = searchsortedlast(H_sweep_kjkg, H_kjkg)
                idx  = clamp(idx, 1, length(T_sweep_K) - 1)
                dH   = H_sweep_kjkg[idx + 1] - H_sweep_kjkg[idx]
                frac = dH > 0 ? clamp((H_kjkg - H_sweep_kjkg[idx]) / dH, 0.0, 1.0) : 0.0
                T    = T_sweep_K[idx]  + frac * (T_sweep_K[idx + 1]  - T_sweep_K[idx])
                rho  = rho_sweep[idx]  + frac * (rho_sweep[idx + 1]  - rho_sweep[idx])
                cp   = cp_sweep[idx]   + frac * (cp_sweep[idx + 1]   - cp_sweep[idx])
                mu   = T < 647.0 ? _water_viscosity(T) : _steam_viscosity(T)
                T_arr[i, j]    = T
                rhoL_arr[i, j] = rho
                rhoV_arr[i, j] = rho
                muL_arr[i, j]  = mu
                muV_arr[i, j]  = mu
                SL_arr[i, j]   = 1.0
                SV_arr[i, j]   = 0.0
                HL_arr[i, j]   = H
                HV_arr[i, j]   = H
                cpL_arr[i, j]  = cp
                cpV_arr[i, j]  = cp
            end
        end
    end

    T_itp    = Jutul.get_2d_interpolator(P_grid, H_grid, T_arr)
    rhoL_itp = Jutul.get_2d_interpolator(P_grid, H_grid, rhoL_arr)
    rhoV_itp = Jutul.get_2d_interpolator(P_grid, H_grid, rhoV_arr)
    muL_itp  = Jutul.get_2d_interpolator(P_grid, H_grid, muL_arr)
    muV_itp  = Jutul.get_2d_interpolator(P_grid, H_grid, muV_arr)
    SL_itp   = Jutul.get_2d_interpolator(P_grid, H_grid, SL_arr)
    SV_itp   = Jutul.get_2d_interpolator(P_grid, H_grid, SV_arr)
    HL_itp   = Jutul.get_2d_interpolator(P_grid, H_grid, HL_arr)
    HV_itp   = Jutul.get_2d_interpolator(P_grid, H_grid, HV_arr)
    cpL_itp  = Jutul.get_2d_interpolator(P_grid, H_grid, cpL_arr)
    cpV_itp  = Jutul.get_2d_interpolator(P_grid, H_grid, cpV_arr)

    # ── H(P, T) table for liquid-region initialisation ────────────────────────
    T_grid   = collect(range(T_min, T_max; length = nT))
    H_pT_arr = Matrix{Float64}(undef, nP, nT)
    for (i, P) in enumerate(P_grid)
        P_mpa   = P / 1e6
        T_sat_K = P_mpa < _Pc_MPa ? SteamTables.Tsat(P_mpa) : Inf
        for (j, T) in enumerate(T_grid)
            # Clamp to subcooled liquid region
            T_eval = isfinite(T_sat_K) ? min(T, T_sat_K - 0.5) : T
            H_pT_arr[i, j] = SteamTables.SpecificH(P_mpa, T_eval) * 1e3
        end
    end
    H_pT_itp = Jutul.get_2d_interpolator(P_grid, T_grid, H_pT_arr)

    tabs = Dict{Symbol, Any}()
    tabs[:T]        = (P, H) -> T_itp(P, H)
    tabs[:rho]      = (P, H) -> (rhoL_itp(P, H), rhoV_itp(P, H))
    tabs[:mu]       = (P, H) -> (muL_itp(P, H),  muV_itp(P, H))
    tabs[:S]        = (P, H) -> (SL_itp(P, H),   SV_itp(P, H))
    tabs[:H_phases] = (P, H) -> (HL_itp(P, H),   HV_itp(P, H))
    tabs[:c_p]      = (P, H) -> (cpL_itp(P, H),  cpV_itp(P, H))
    tabs[:H_pT]     = (P, T) -> H_pT_itp(P, T)
    return tabs
end
