## Variable types for the pressure-enthalpy (P,H) thermal formulation

"""
    Enthalpy{T} <: ScalarVariable

Specific enthalpy [J/kg] as a primary variable, for use in the pressure-enthalpy
(P,H) thermal formulation. The default bounds cover liquid water up through
superheated steam.
"""
Base.@kwdef struct Enthalpy{T} <: ScalarVariable
    min::T = 0.0
    max::T = 1e8
    max_rel::Union{T, Nothing} = nothing
    max_abs::Union{T, Nothing} = nothing
end

Jutul.default_value(model, ::Enthalpy) = 125e3  # ~30°C liquid water [J/kg]
Jutul.minimum_value(h::Enthalpy) = h.min
Jutul.maximum_value(h::Enthalpy) = h.max
Jutul.absolute_increment_limit(h::Enthalpy) = h.max_abs
Jutul.relative_increment_limit(h::Enthalpy) = h.max_rel

function JutulDarcy.model_is_thermal(model::SimulationModel)
    pvars = Jutul.get_primary_variables(model)
    return haskey(pvars, :Temperature) || haskey(pvars, :Enthalpy)
end


## Temperature from (P, H) ─────────────────────────────────────────────────────

"""
    TemperatureFromEnthalpy{T} <: ScalarVariable

Secondary variable that computes temperature [K] from pressure and specific
enthalpy via a `(P, H) → T` lookup table callable.
"""
struct TemperatureFromEnthalpy{T} <: ScalarVariable
    tab::T
end

Jutul.subvariable(v::TemperatureFromEnthalpy, map) = v

@jutul_secondary function update_temperature_from_enthalpy!(T_out, var::TemperatureFromEnthalpy, model, Pressure, Enthalpy, ix)
    for c in ix
        T_out[c] = var.tab(Pressure[c], Enthalpy[c])
    end
end

struct FluidInternalEnergyFromEnthalpy <: JutulDarcy.PhaseVariables end
@jutul_secondary function update_internal_energy_from_enthalpy!(U, var::FluidInternalEnergyFromEnthalpy, model, FluidEnthalpy, PhaseMassDensities, Pressure, ix)
    nph = size(PhaseMassDensities, 1)
    for ph in 1:nph
        for c in ix
            H = FluidEnthalpy[ph, c]
            p = Pressure[c]
            rho = PhaseMassDensities[ph, c]
            U[ph, c] = H - p./rho
        end
    end
end

## (P, H)-dependent vector variable ───────────────────────────────────────────

"""
    PressureEnthalpyDependentVariable{T, N} <: VectorVariables

Vector secondary variable backed by a `(P, H) → SVector/Vector` lookup table.
Used to represent pressure-enthalpy-dependent fluid properties (density,
viscosity, heat capacity) in the H-formulation. `N` is the number of returned
values per cell and is inferred at construction time by probing the table.
"""
struct PressureEnthalpyDependentVariable{T, N} <: VectorVariables
    tab::T
    function PressureEnthalpyDependentVariable(tab)
        N = length(tab(1e8, 125e3))
        new{typeof(tab), N}(tab)
    end
end

Jutul.subvariable(v::PressureEnthalpyDependentVariable, map) = v
Jutul.values_per_entity(model, ::PressureEnthalpyDependentVariable{T, N}) where {T, N} = N

@jutul_secondary function update_ph_dependent!(result, var::PressureEnthalpyDependentVariable{T, N}, model, Pressure, Enthalpy, ix) where {T, N}
    for c in ix
        vals = var.tab(Pressure[c], Enthalpy[c])
        for i in 1:N
            result[i, c] = vals[i]
        end
    end
    return result
end





## Saturations from (P, H) ────────────────────────────────────────────────────

"""
    SaturationsFromEnthalpy{T, N} <: VectorVariables

Secondary variable that computes per-phase volume saturations from pressure and
specific enthalpy via a `(P, H) → NTuple{N}` lookup table.  `N` is inferred at
construction time by probing the table at a representative state.

Used in the two-phase (liquid + vapour) enthalpy formulation to replace
`Saturations` as a primary variable; the state is fully determined by (P, H).
"""
struct SaturationsFromEnthalpy{T, N} <: VectorVariables
    tab::T
    function SaturationsFromEnthalpy(tab)
        N = length(tab(1e5, 500e3))
        new{typeof(tab), N}(tab)
    end
end

Jutul.subvariable(v::SaturationsFromEnthalpy, map) = v
Jutul.values_per_entity(model, ::SaturationsFromEnthalpy{T, N}) where {T, N} = N

@jutul_secondary function update_saturations_from_enthalpy!(S, var::SaturationsFromEnthalpy{T, N}, model, Pressure, Enthalpy, ix) where {T, N}
    for c in ix
        vals = var.tab(Pressure[c], Enthalpy[c])
        for i in 1:N
            S[i, c] = vals[i]
        end
    end
    return S
end

## Smooth per-phase 2D-table variables (OBL approach) ────────────────────────
# All four types read from pre-built (P, H) → scalar tables, one per phase.
# No branching on phase state → fully AD-compatible Jacobians.

"""
    LVPhaseDensity{TL, TV} <: PhaseVariables

Per-phase mass density [kg/m³] from smooth 2D (P,H) lookup tables.
"""
struct LVPhaseDensity{TL, TV} <: JutulDarcy.PhaseVariables
    tab_liq::TL
    tab_vap::TV
end

@jutul_secondary function update_lv_phase_density!(result, var::LVPhaseDensity, model, Pressure, Enthalpy, ix)
    for c in ix
        p = Pressure[c]
        h = Enthalpy[c]
        result[1, c] = var.tab_liq(p, h)
        result[2, c] = var.tab_vap(p, h)
    end
    return result
end

"""
    LVPhaseViscosity{TL, TV} <: PhaseVariables

Per-phase dynamic viscosity [Pa·s] from smooth 2D (P,H) lookup tables.
"""
struct LVPhaseViscosity{TL, TV} <: JutulDarcy.PhaseVariables
    tab_liq::TL
    tab_vap::TV
end

@jutul_secondary function update_lv_phase_viscosity!(result, var::LVPhaseViscosity, model, Pressure, Enthalpy, ix)
    for c in ix
        p = Pressure[c]
        h = Enthalpy[c]
        result[1, c] = var.tab_liq(p, h)
        result[2, c] = var.tab_vap(p, h)
    end
    return result
end

"""
    LVPhaseEnthalpy{TL, TV} <: PhaseVariables

Per-phase specific enthalpy [J/kg] from smooth 2D (P,H) lookup tables.
"""
struct LVPhaseEnthalpy{TL, TV} <: JutulDarcy.PhaseVariables
    tab_liq::TL
    tab_vap::TV
end

@jutul_secondary function update_lv_phase_enthalpy!(H_out, var::LVPhaseEnthalpy, model, Pressure, Enthalpy, Saturations, ix)
    for c in ix
        p = Pressure[c]
        h = Enthalpy[c]
        S_l = Saturations[1, c]
        ϵ = 1e-6
        if ϵ < S_l < 1 - ϵ
            H_out[1, c] = var.tab_liq(p, h)
            H_out[2, c] = var.tab_vap(p, h)
        else
            H_out[1, c] = h
            H_out[2, c] = h
        end
    end
    return H_out
end

"""
    LVPhaseSaturation{TV} <: PhaseVariables

Per-phase volume saturation [-] from a smooth 2D (P,H) → S_vapor lookup table.
Liquid saturation is `1 - S_vapor`.
"""
struct LVPhaseSaturation{TV} <: JutulDarcy.PhaseVariables
    tab_vap::TV
end

@jutul_secondary function update_lv_phase_saturation!(S, var::LVPhaseSaturation, model, Pressure, Enthalpy, ix)
    for c in ix
        p = Pressure[c]
        h = Enthalpy[c]
        S_v = var.tab_vap(p, h)
        S[1, c] = 1 - S_v
        S[2, c] = S_v
    end
    return S
end

struct EnthalpyFromPT{T} <: ScalarVariable
    tab::T
end

@jutul_secondary function update_enthalpy_from_pt!(H, var::EnthalpyFromPT, model, Pressure, Temperature, ix)
    for c in ix
        p = Pressure[c]
        T = Temperature[c]
        H[c] = var.tab(p, T)
    end
    return H
end


# struct LVTotalThermalEnergy <: ScalarVariable end

# @jutul_secondary function update_lv_total_thermal_energy!(E_total, var::LVTotalThermalEnergy, model,
#     Pressure, Enthalpy, TotalMasses, Saturations, PhaseMassDensities, RockDensity, RockInternalEnergy, BulkVolume, FluidVolume, ix)
#     U_r = RockInternalEnergy
#     ρ_r = RockDensity
#     M_f = TotalMasses
#     V_f = FluidVolume
#     V = BulkVolume
#     h = Enthalpy
#     ρ_f = PhaseMassDensities
#     S = Saturations
#     println("Updating total thermal energy for $(length(ix)) cells...")
#     for i in ix
#         ρ_fm = ρ_f[1, i]*S[1, i] + ρ_f[2, i]*S[2, i]
#         V_r = V[i] - V_f[i]
#         E_r = ρ_r[i] * U_r[i] * V_r
#         p = Pressure[i]
#         E_f = M_f[i] * (h[i] - p/ρ_fm)
#         E_total[i] = E_r + E_f
#     end
#     return E_total
# end


# const MSWellDomain = DiscretizedDomain{<:MultiSegmentWell}
# const MSWellFlowModel = SimulationModel{<:MSWellDomain, <:MultiPhaseSystem}
# @jutul_secondary function update_lv_total_thermal_energy!(E_total, te::LVTotalThermalEnergy, model::MSWellFlowModel,
#     Pressure, Enthalpy, TotalMasses, Saturations, PhaseMassDensities, MaterialDensities, MaterialInternalEnergy, BulkVolume, FluidVolume, ix)
#     update_lv_total_thermal_energy!(E_total, te::LVTotalThermalEnergy, nothing, Pressure, Enthalpy, TotalMasses, Saturations, PhaseMassDensities, MaterialDensities, MaterialInternalEnergy, BulkVolume, FluidVolume, ix)
# end