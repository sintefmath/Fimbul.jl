## Fimbul overload of JutulDarcy.setup_reservoir_model_geothermal
## that adds support for the pressure-enthalpy (P,H) formulation.
##
## Because Julia method dispatch is positional-only, defining a keyword-extended
## version of the same function replaces JutulDarcy's original for the
## (DataDomain,) positional signature. The :temperature branch therefore
## replicates JutulDarcy's original body via the private helper below.

# ── Private helper: mirrors JutulDarcy's setup_reservoir_model_geothermal body ─

# function _setup_geothermal_temperature___(
#         reservoir::DataDomain;
#         thermal = true,
#         extra_out = false,
#         parameters = Dict{Symbol, Any}(),
#         salt_mole_fractions = Float64[],
#         salt_names = String[],
#         table_arg = NamedTuple(),
#         single_phase = true,
#         update_reservoir = true,
#         table_cache = Dict(),
#         kwarg...
#     )
#     thermal || throw(ArgumentError("Cannot setup geothermal reservoir model with thermal = false"))
#     tables = JutulDarcy.Geothermal.geothermal_setup_tables(
#         table_cache, salt_names, salt_mole_fractions, table_arg
#     )

#     rhoWS = first(JutulDarcy.reference_densities(:co2brine))
#     if single_phase
#         sys = SinglePhaseSystem(AqueousPhase(), reference_density = rhoWS)
#     else
#         error("Multiphase geothermal is not implemented")
#     end
#     cond_water = tables[:phase_conductivity](1*si_unit(:atm), 273.15 + 20.0)
#     if update_reservoir
#         reservoir[:fluid_thermal_conductivity] .= cond_water
#     end
#     model = setup_reservoir_model(reservoir, sys; thermal = true, extra_out = false, kwarg...)

#     rho = JutulDarcy.PressureTemperatureDependentVariable(tables[:density])
#     c_p = JutulDarcy.PressureTemperatureDependentVariable(tables[:heat_capacity_constant_pressure])
#     mu  = JutulDarcy.PTViscosities(tables[:viscosity])
#     for (k, m) in pairs(model.models)
#         if k == :Reservoir || JutulDarcy.model_or_domain_is_well(m)
#             set_secondary_variables!(m;
#                 PhaseMassDensities   = rho,
#                 PhaseViscosities     = mu,
#                 ComponentHeatCapacity = c_p,
#             )
#         end
#     end

#     rmodel = reservoir_model(model)
#     outvar = rmodel.output_variables
#     push!(outvar, :PhaseMassDensities, :RockInternalEnergy, :FluidInternalEnergy, :TotalThermalEnergy)
#     unique!(outvar)

#     if extra_out
#         parameters = setup_parameters(model, parameters)
#         return (model, parameters)
#     else
#         return model
#     end
# end

function JutulDarcy.add_thermal_to_model!(model::SimulationModel{D, S, F, C}) where {D, S<:GeothermalTwoPhaseSystem, F, C}

    invoke(JutulDarcy.add_thermal_to_model!, Tuple{SimulationModel}, model)
    Jutul.delete_variable!(model, :Temperature)

    pvt = model.system.pvt_tables

    set_secondary_variables!(model,
        Temperature = TemperatureFromEnthalpy(pvt[:temperature]),
        FluidEnthalpy = PhaseEnthalpyH2O(pvt[:enthalpy_liquid_ph], pvt[:enthalpy_vapor_ph]),
        FluidInternalEnergy = FluidInternalEnergyFromEnthalpy(),
    )
    set_primary_variables!(model, Enthalpy = Enthalpy())
    model.extra[:enthalpy] = pvt[:enthalpy]
    return model
    # set_primary_variables!(model, Temperature = Temperature())
    # set_parameters!(model,
    #     RockHeatCapacity = RockHeatCapacity(),
    #     RockDensity = RockDensity(),
    #     BulkVolume = BulkVolume(),
    #     ComponentHeatCapacity = ComponentHeatCapacity(),
    # )
    # set_secondary_variables!(model,
    #     FluidInternalEnergy = FluidInternalEnergy(),
    #     FluidEnthalpy = FluidEnthalpy(),
    #     TotalThermalEnergy = TotalThermalEnergy(),
    # )
    # is_reservoir = !model_or_domain_is_well(model)
    # if is_reservoir
    #     set_parameters!(model,
    #         RockThermalConductivities = RockThermalConductivities(),
    #         FluidThermalConductivities = FluidThermalConductivities()
    #     )
    #     set_secondary_variables!(model,
    #         RockInternalEnergy = RockInternalEnergy()
    #     )
    # else
    #     if model_or_domain_is_well(model)
    #         w = physical_representation(model.domain)

    #         set_parameters!(model,
    #             WellIndicesThermal = WellIndicesThermal(),
    #         )
    #         if w isa MultiSegmentWell
    #             set_parameters!(model,
    #                 MaterialThermalConductivities = MaterialThermalConductivities(),
    #                 MaterialHeatCapacities = MaterialHeatCapacities(),
    #                 MaterialDensities = MaterialDensities()
    #             )
    #             set_secondary_variables!(model,
    #                 MaterialInternalEnergy = MaterialInternalEnergy()
    #             )

    #         else
    #             w::SimpleWell
    #             set_secondary_variables!(model,
    #                 RockInternalEnergy = RockInternalEnergy()
    #             )
    #         end
    #     end
    # end
    # disc = model.domain.discretizations.heat_flow
    # model.equations[:energy_conservation] = ConservationLaw(disc, :TotalThermalEnergy, 1)

    # out = model.output_variables
    # push!(out, :TotalThermalEnergy)
    # push!(out, :FluidEnthalpy)
    # push!(out, :Temperature)

    # unique!(out)
    # return model
end

# ── Enthalpy-formulation model modifier ──────────────────────────────────────

function _apply_enthalpy_formulation!(model, enthalpy_tables)

    out = model.output_variables
    push!(out, :Enthalpy)
    # Switch Temperature from primary to secondary.
    # set_secondary_variables! internally calls delete_variable!, which removes
    # Temperature from primary_variables before adding it to secondary_variables.
    set_primary_variables!(model, Enthalpy = Enthalpy())
    set_secondary_variables!(model,
        Temperature = TemperatureFromEnthalpy(enthalpy_tables[:temperature]),
    )
    # Replace P,T-dependent fluid properties with smooth per-phase (P,H) tables
    set_secondary_variables!(model,
        PhaseMassDensities = LVPhaseDensity(
            enthalpy_tables[:density_liquid_ph],
            enthalpy_tables[:density_vapor_ph]),
        PhaseViscosities = LVPhaseViscosity(
            enthalpy_tables[:viscosity_liquid_ph],
            enthalpy_tables[:viscosity_vapor_ph]),
        FluidEnthalpy = LVPhaseEnthalpy(
            enthalpy_tables[:enthalpy_liquid_ph],
            enthalpy_tables[:enthalpy_vapor_ph]),
        FluidInternalEnergy = FluidInternalEnergyFromEnthalpy(),
        # TotalThermalEnergy = LVTotalThermalEnergy()
    )
    Jutul.delete_variable!(model, :ComponentHeatCapacity)
    set_secondary_variables!(model,
        Saturations = LVPhaseSaturation(enthalpy_tables[:saturation_vapor_ph]),
    )
    push!(out, :Saturations)
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
            if k == :Reservoir
                _apply_enthalpy_formulation!(m, enthalpy_tables)
            elseif JutulDarcy.model_or_domain_is_well(m)
                # Skip phase-split (Saturations) for wells: the single-phase well
                # equations (saturation_mixed) only support 1-row Saturations.
                _apply_enthalpy_formulation!(m, enthalpy_tables)
            end
        end
        return out
    else
        throw(ArgumentError("Unknown formulation :$formulation — use :temperature or :enthalpy"))
    end
end

# ── Two-phase geothermal enthalpy setup ───────────────────────────────────────

"""
    setup_reservoir_model_geothermal_2ph(reservoir; enthalpy_tables, kwarg...)

Set up a two-phase (liquid + vapour) geothermal reservoir model using the
pressure-enthalpy (P, H) formulation backed by `GeothermalTwoPhaseSystem`.

The system has **1 component** (pure H₂O) and **2 phases** (liquid + vapour),
giving exactly 2 equations (1 mass + 1 energy) for 2 unknowns (P, H).
Phase saturations, densities, viscosities and enthalpies are all derived from
the supplied `enthalpy_tables` (built by `build_steam_tables_2ph()`).

Relative permeabilities default to `BrooksCoreyRelativePermeabilities` with
exponent = 1, residual = 0, endpoint = 1 (i.e. `kr_α = S_α`, linear).

`enthalpy_tables` must be a `Dict` with callable entries:
- `:T`        : `(P,H) → T [K]`
- `:rho`      : `(P,H) → (ρ_l, ρ_v) [kg/m³]`
- `:mu`       : `(P,H) → (μ_l, μ_v) [Pa·s]`
- `:S`        : `(P,H) → (S_l, S_v)` (volume saturations)
- `:H_phases` : `(P,H) → (H_l, H_v) [J/kg]`
- `:c_p`      : `(P,H) → (cp_l, cp_v) [J/(kg·K)]`
- `:H_pT`     : `(P,T) → H [J/kg]`  (for T-based initialisation)
"""
function setup_reservoir_model_geothermal_2ph(
        reservoir::DataDomain;
        enthalpy_tables,
        extra_out = false,
        update_reservoir = true,
        kwarg...
    )
    isnothing(enthalpy_tables) && throw(ArgumentError("enthalpy_tables must be provided"))

    # Reference densities probed from the tables at representative conditions
    rhoL_ref = first(enthalpy_tables[:density_liquid](1e6))   # liquid  at 10 bar, ~95 °C
    rhoV_ref = last(enthalpy_tables[:density_vapor](1e6))   # vapour  at 10 bar, ~220 °C

    sys = GeothermalTwoPhaseSystem(
        reference_densities = (rhoL_ref, rhoV_ref),
    )

    # Nominal fluid thermal conductivity (liquid water at ~20 °C)
    if update_reservoir
        reservoir[:fluid_thermal_conductivity] .= 0.6   # W/(m·K)
    end

    model = setup_reservoir_model(reservoir, sys; thermal = true, kwarg...)

    # Apply the (P, H) primary-variable formulation to every sub-model.
    # add_phase_split = true is required for both reservoir and wells:
    #   - Reservoir: 2-row SaturationsFromEnthalpy feeds RelativePermeabilities.
    #   - Wells:     2-row SaturationsFromEnthalpy is accessed by the
    #                GeothermalTwoPhaseSystem perforation-flux dispatch when
    #                computing injection fluxes (multisegment path).
    for (k, m) in pairs(model.models)
        if k == :Reservoir || JutulDarcy.model_or_domain_is_well(m)
            _apply_enthalpy_formulation!(m, enthalpy_tables)
        end
    end

    rmodel = reservoir_model(model)
    outvar = rmodel.output_variables
    push!(outvar, :PhaseMassDensities, :Saturations, :RockInternalEnergy,
                  :FluidInternalEnergy, :TotalThermalEnergy)
    unique!(outvar)

    if extra_out
        parameters = setup_parameters(model)
        return (model, parameters)
    else
        return model
    end
end

# ── setup_reservoir_state override for GeothermalTwoPhaseSystem ──────────────
#
# JutulDarcy's MultiModel setup_reservoir_state uses `model_is_thermal` to
# decide whether to read `:Temperature` from each well's init dict.  With
# GeothermalTwoPhaseSystem the thermal primary variable is `:Enthalpy`, not
# `:Temperature`, so the upstream code throws a KeyError.
#
# This overload replicates the upstream logic but:
#   - reads `:Enthalpy` from the well dicts, and
#   - converts it to temperature via the system's T(P,H) table to populate
#     the facility's `SurfaceTemperature`.

# function JutulDarcy.setup_reservoir_state___(
#         model::MultiModel,
#         equil::Union{Missing, Vector, JutulDarcy.EquilibriumRegion} = missing;
#         kwarg...,
#     )
#     rmodel = JutulDarcy.reservoir_model(model)

#     # Dispatch to our specialised version for GeothermalTwoPhaseSystem
#     if rmodel.system isa GeothermalTwoPhaseSystem
#         return _setup_reservoir_state_geothermal_2ph(model, equil; kwarg...)
#     end

#     # Fall back to JutulDarcy's original implementation for all other systems
#     return invoke(
#         JutulDarcy.setup_reservoir_state,
#         Tuple{MultiModel, Union{Missing, Vector, JutulDarcy.EquilibriumRegion}},
#         model, equil; kwarg...,
#     )
# end

function _setup_reservoir_state_geothermal_2ph(
        model::MultiModel,
        equil::Union{Missing, Vector, JutulDarcy.EquilibriumRegion} = missing;
        kwarg...,
    )
    rmodel = JutulDarcy.reservoir_model(model)
    pvars  = collect(keys(Jutul.get_primary_variables(rmodel)))

    # Initialise the reservoir submodel (handles P, H keyword args)
    res_state = JutulDarcy.setup_reservoir_state(rmodel, equil; kwarg...)

    init = Dict{Symbol, Any}(:Reservoir => res_state)

    perf_subset(v::AbstractVector, i) = v[i]
    perf_subset(v::AbstractMatrix, i) = v[:, i]
    perf_subset(v, i) = v

    # Retrieve the T(P,H) table from the secondary variable definition
    t_tab = rmodel.secondary_variables[:Temperature].tab

    # ── Individual well models ──────────────────────────────────────────────
    for k in keys(model.models)
        k == :Reservoir && continue
        W = model.models[k]
        W.domain isa JutulDarcy.WellGroup && continue  # handled in second pass

        init_w = Dict{Symbol, Any}()
        wg = Jutul.physical_representation(W.domain)
        if wg isa JutulDarcy.MultiSegmentWell
            init_w[:TotalMassFlux] = 0.0
        end
        c = JutulDarcy.map_well_nodes_to_reservoir_cells(wg, rmodel.data_domain)
        for pk in pvars
            pv = res_state[pk]
            init_w[pk] = perf_subset(pv, c)
        end
        init[k] = init_w
    end

    # ── Facility / WellGroup ───────────────────────────────────────────────
    T_float = Float64
    for (k, _) in JutulDarcy.get_model_wells(model)
        T_float = promote_type(T_float, eltype(init[k][:Pressure]))
    end

    for (k, W) in pairs(model.models)
        W.domain isa JutulDarcy.WellGroup || continue

        init_arg = Dict{Symbol, Any}()
        init_arg[:TotalSurfaceMassRate] = 0.0
        init_arg[:SurfacePhaseRates]    = 0.0

        own_wells = W.domain.well_symbols
        bh   = zeros(T_float, length(own_wells))
        temp = similar(bh)

        for (i, w) in enumerate(own_wells)
            bh[i] = init[w][:Pressure][1]
            H_w   = init[w][:Enthalpy][1]
            P_w   = init[w][:Pressure][1]
            temp[i] = t_tab(P_w, H_w)   # T [K] from (P,H)
        end

        init_arg[:BottomHolePressure] = bh
        init_arg[:SurfaceTemperature] = temp
        init[k] = Jutul.setup_state(W; pairs(init_arg)...)
    end

    return Jutul.setup_state(model, init)
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

JutulDarcy.system_uses_cnv_mb(system::GeothermalTwoPhaseSystem) = true
