## Fimbul overload of JutulDarcy.setup_reservoir_model_geothermal
## that adds support for the pressure-enthalpy (P,H) formulation.
##
## Because Julia method dispatch is positional-only, defining a keyword-extended
## version of the same function replaces JutulDarcy's original for the
## (DataDomain,) positional signature. The :temperature branch therefore
## replicates JutulDarcy's original body via the private helper below.

# ── Private helper: mirrors JutulDarcy's setup_reservoir_model_geothermal body ─

function _setup_geothermal_temperature(
        reservoir::DataDomain;
        thermal = true,
        extra_out = false,
        parameters = Dict{Symbol, Any}(),
        salt_mole_fractions = Float64[],
        salt_names = String[],
        table_arg = NamedTuple(),
        single_phase = true,
        update_reservoir = true,
        table_cache = Dict(),
        kwarg...
    )
    thermal || throw(ArgumentError("Cannot setup geothermal reservoir model with thermal = false"))
    tables = JutulDarcy.Geothermal.geothermal_setup_tables(
        table_cache, salt_names, salt_mole_fractions, table_arg
    )

    rhoWS = first(JutulDarcy.reference_densities(:co2brine))
    if single_phase
        sys = SinglePhaseSystem(AqueousPhase(), reference_density = rhoWS)
    else
        error("Multiphase geothermal is not implemented")
    end
    cond_water = tables[:phase_conductivity](1*si_unit(:atm), 273.15 + 20.0)
    if update_reservoir
        reservoir[:fluid_thermal_conductivity] .= cond_water
    end
    model = setup_reservoir_model(reservoir, sys; thermal = true, extra_out = false, kwarg...)

    rho = JutulDarcy.PressureTemperatureDependentVariable(tables[:density])
    c_p = JutulDarcy.PressureTemperatureDependentVariable(tables[:heat_capacity_constant_pressure])
    mu  = JutulDarcy.PTViscosities(tables[:viscosity])
    for (k, m) in pairs(model.models)
        if k == :Reservoir || JutulDarcy.model_or_domain_is_well(m)
            set_secondary_variables!(m;
                PhaseMassDensities   = rho,
                PhaseViscosities     = mu,
                ComponentHeatCapacity = c_p,
            )
        end
    end

    rmodel = reservoir_model(model)
    outvar = rmodel.output_variables
    push!(outvar, :PhaseMassDensities, :RockInternalEnergy, :FluidInternalEnergy, :TotalThermalEnergy)
    unique!(outvar)

    if extra_out
        parameters = setup_parameters(model, parameters)
        return (model, parameters)
    else
        return model
    end
end

# ── Enthalpy-formulation model modifier ──────────────────────────────────────

function _apply_enthalpy_formulation!(model, enthalpy_tables)
    # Switch Temperature from primary to secondary.
    # set_secondary_variables! internally calls delete_variable!, which removes
    # Temperature from primary_variables before adding it to secondary_variables.
    set_primary_variables!(model, Enthalpy = Enthalpy())
    set_secondary_variables!(model,
        Temperature = TemperatureFromEnthalpy(enthalpy_tables[:T]),
    )
    # Replace P,T-dependent fluid properties with P,H-dependent ones
    set_secondary_variables!(model,
        PhaseMassDensities    = PressureEnthalpyDependentVariable(enthalpy_tables[:rho]),
        PhaseViscosities      = PressureEnthalpyDependentVariable(enthalpy_tables[:mu]),
        ComponentHeatCapacity = PressureEnthalpyDependentVariable(enthalpy_tables[:c_p]),
        FluidInternalEnergy   = FluidInternalEnergyFromEnthalpy(),
    )
    if haskey(enthalpy_tables, :H_phases)
        set_secondary_variables!(model,
            FluidEnthalpy = PressureEnthalpyDependentVariable(enthalpy_tables[:H_phases]),
        )
    end
    out = model.output_variables
    push!(out, :Enthalpy)
    unique!(out)
    return model
end

# ── Public overload ───────────────────────────────────────────────────────────

"""
    JutulDarcy.setup_reservoir_model_geothermal(reservoir; formulation = :temperature, enthalpy_tables = nothing, kwarg...)

Fimbul extension of `JutulDarcy.setup_reservoir_model_geothermal` that adds
support for the pressure-enthalpy (P,H) formulation via `formulation = :enthalpy`.

When `formulation = :temperature` (default), behaviour is identical to
JutulDarcy's original implementation.

When `formulation = :enthalpy`, the model is first set up with the T-formulation
(using JutulDarcy's geothermal (P,T) tables for the model structure), then
reconfigured so that:
- `Enthalpy` replaces `Temperature` as the primary variable.
- `Temperature` becomes a secondary variable computed via a `(P,H) → T` table.
- `PhaseMassDensities`, `PhaseViscosities`, and `ComponentHeatCapacity` are
  replaced by `PressureEnthalpyDependentVariable` instances.

`enthalpy_tables` must be a `NamedTuple` with callable entries:
- `enthalpy_tables[:T]`   : `(P,H) → T [K]` (scalar)
- `enthalpy_tables[:rho]` : `(P,H) → [ρ]` (vector, kg/m³ per phase)
- `enthalpy_tables[:mu]`  : `(P,H) → [μ]` (vector, Pa·s per phase)
- `enthalpy_tables[:c_p]` : `(P,H) → [cₚ]` (vector, J/(kg·K) per phase)

All other keyword arguments are forwarded to JutulDarcy's geothermal setup.
"""
function JutulDarcy.Geothermal.setup_reservoir_model_geothermal(
        reservoir::DataDomain, formulation::Symbol;
        enthalpy_tables = nothing,
        kwarg...
    )
    if formulation == :temperature
        return JutulDarcy.Geothermal.setup_reservoir_model_geothermal(reservoir; kwarg...)
    elseif formulation == :enthalpy
        isnothing(enthalpy_tables) && throw(ArgumentError(
            "enthalpy_tables must be provided for formulation = :enthalpy"))
        out = JutulDarcy.Geothermal.setup_reservoir_model_geothermal(reservoir; kwarg...)
        # out is either a model or (model, parameters) depending on extra_out
        model = out isa Tuple ? first(out) : out
        for (k, m) in pairs(model.models)
            if k == :Reservoir || JutulDarcy.model_or_domain_is_well(m)
                _apply_enthalpy_formulation!(m, enthalpy_tables)
            end
        end
        return out
    else
        throw(ArgumentError("Unknown formulation :$formulation — use :temperature or :enthalpy"))
    end
end

# ── Convergence criterion override ────────────────────────────────────────────
#
# JutulDarcy's temperature_increment calls update_report[:Temperature], which
# throws a KeyError when Temperature is not a primary variable. Override it
# here to handle both formulations safely.

function JutulDarcy.temperature_increment(model, state, update_report)
    if haskey(update_report, :Temperature)
        t_report = update_report[:Temperature]
        return haskey(t_report, :max) ? t_report.max : 1.0
    elseif haskey(update_report, :Enthalpy)
        h_report = update_report[:Enthalpy]
        return haskey(h_report, :max) ? h_report.max : 1.0
    end
    return 1.0
end
