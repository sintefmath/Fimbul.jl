## Variable types for the pressure-enthalpy (P,H) thermal formulation

"""
    Enthalpy{T} <: ScalarVariable

Specific enthalpy [J/kg] as a primary variable, for use in the pressure-enthalpy
(P,H) thermal formulation. The default bounds cover liquid water up through
superheated steam.
"""
Base.@kwdef struct Enthalpy{T} <: ScalarVariable
    min::T = 0.0
    max::T = 4e6
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
const PCRIT_WATER_CP = 22.064e6  # Critical pressure of water in Pa

function water_phase(p, h, h_l, h_v)
    # 1 -> Subcooled liquid
    # 2 -> Superheated vapor
    # 3 -> Two-phase
    # 4 -> Supercritical or liquid
    p, h, h_l, h_v = value(p), value(h), value(h_l), value(h_v)
    if p < PCRIT_WATER_CP
        if h <= h_l
            return 1
        elseif h_l < h < h_v
            return 2
        elseif h >= h_v
            return 3
        end
    else
        return 4
    end
end

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

struct PressureEnthalpyDependentPhaseVariable{TLV, TM, N} <: JutulDarcy.PhaseVariables
    tab_l::TLV
    tab_v::TLV
    tab_m::TM
    h_l::TLV
    h_v::TLV
    
    function PressureEnthalpyDependentPhaseVariable(tab_l, tab_v, tab_m, h_l, h_v)
        N = length(tab_m(1e8, 125e3))
        new{typeof(tab_l), typeof(tab_m), N}(tab_l, tab_v, tab_m, h_l, h_v)
    end
end

@jutul_secondary function update_ph_dependent_phase!(result, var::PressureEnthalpyDependentPhaseVariable{TLV, TM, N}, model, Pressure, Enthalpy, FluidEnthalpy, ix) where {TLV, TM, N}
    for c in ix
        h = Enthalpy[c]
        p = Pressure[c]
        h_l = var.h_l(p)
        h_v = var.h_v(p)
        phase = water_phase(p, h, h_l, h_v)
        if phase == 2 # Two-phase
            # Subcooled liquid
            x_l = var.tab_l(p)
            x_v = var.tab_v(p)
            result[1, c] = x_l
            result[2, c] = x_v
        else # Supercooled liquid/superheated vapor/supercritical
            x = var.tab_m(p, h)
            result[1, c] = x
            result[2, c] = x
        end
    end
    return result
end

struct GeothermalLVFluidEnthalpy{T, N} <: JutulDarcy.PhaseVariables
    tab_l::T
    tab_v::T
    function GeothermalLVFluidEnthalpy(tab_l, tab_v)
        N = 2
        new{typeof(tab_l), N}(tab_l, tab_v) # N is fixed at 2 for liquid and vapor phases
    end
end
@jutul_secondary function update_geothermal_lv_fluid_enthalpy!(H, var::GeothermalLVFluidEnthalpy, model, Pressure, Enthalpy, ix)
    for c in ix
        h = Enthalpy[c]
        p = Pressure[c]
        h_l = var.tab_l(p)
        h_v = var.tab_v(p)
        phase = water_phase(p, h, h_l, h_v)
        if phase == 2 # Two-phase
            H[1, c] = h_l
            H[2, c] = h_v
        else
            # In single-phase regions, assign the available phase enthalpy to both rows for consistency
            H[1, c] = h
            H[2, c] = h
        end
    end
    return H
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

struct GeothermalLVSaturation{TLV, N} <: JutulDarcy.PhaseVariables
    h_l::TLV
    h_v::TLV
    function GeothermalLVSaturation(h_l, h_v)
        N = 2
        new{typeof(h_l), N}(h_l, h_v)
    end
end
@jutul_secondary function update_geothermal_liquid_vapor_saturation!(S, var::GeothermalLVSaturation, model, Pressure, Enthalpy, FluidEnthalpy, PhaseMassDensities, ix)
    for c in ix
        h = Enthalpy[c]
        p = Pressure[c]
        h_l = var.h_l(p)
        h_v = var.h_v(p)
        phase = water_phase(p, h, h_l, h_v)
        if phase == 1
            S[1, c] = 1.0  # Fully liquid
            S[2, c] = 0.0
        elseif phase == 2 # Two-phase
            ρ_l = PhaseMassDensities[1, c]
            ρ_v = PhaseMassDensities[2, c]
            x = (h - h_l) / (h_v - h_l)  # Vapor quality
            # Convert to volume saturation
            S_v = (x/ρ_v) / (x/ρ_v + (1 - x)/ρ_l)
            S[1, c] = 1 - S_v  # Liquid saturation
            S[2, c] = S_v      # Vapor saturation
        elseif phase == 3 || phase == 4
            S[1, c] = 0.0  # Fully vapor
            S[2, c] = 1.0
        else
            error("Invalid phase detected in saturation update")
        end
    end
    return S
end