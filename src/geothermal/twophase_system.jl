## Two-phase (liquid + vapour), single-component geothermal system.
##
## Contrast with `ImmiscibleSystem`: for an ImmiscibleSystem with N phases,
## `number_of_components = N`, yielding N mass equations.  For pure water that
## would give 2 mass equations + 1 energy equation for only 2 unknowns (P, H).
## `GeothermalTwoPhaseSystem` fixes this by having `number_of_components = 1`,
## so the equation count (1 mass + 1 energy) matches the primary variables
## (Pressure + Enthalpy).
##
## The two phases share one component (H₂O).  Phase split, densities,
## viscosities and enthalpies are derived from (P, H) look-up tables supplied
## by `build_steam_tables_2ph`.

const WATER_CRITICAL_PRESSURE = 22.064e6
const WATER_CRITICAL_ENTHALPY = 2085e3

# ── Struct ────────────────────────────────────────────────────────────────────

"""
    GeothermalTwoPhaseSystem

Two-phase (liquid + vapour), single-component fluid system for geothermal
(P, H) simulations.

| quantity               | value              |
|:-----------------------|:-------------------|
| `number_of_phases`     | 2                  |
| `number_of_components` | 1  (pure H₂O)      |

One mass conservation equation is closed by the `Enthalpy` primary variable
(alongside `Pressure`).  Phase split, densities, viscosities and enthalpies are
obtained via Steam-Tables look-ups from `build_steam_tables_2ph`.

Use `setup_reservoir_model_geothermal_2ph` to obtain a fully configured model.
"""
struct GeothermalTwoPhaseSystem{T <: Tuple, F <: NTuple} <: JutulDarcy.AbstractCompositionalSystemLV
    phases  :: T
    rho_ref :: F
    reference_phase_index::Int
    pvt_tables::Dict
end

"""
    GeothermalTwoPhaseSystem(; reference_densities = (1000.0, 1.0))

Construct a `GeothermalTwoPhaseSystem` for a `(AqueousPhase, VaporPhase)` pair.
`reference_densities` should be `(ρ_liquid_ref, ρ_vapour_ref)` in kg/m³.
"""
function GeothermalTwoPhaseSystem(pvt_tables::Dict;
        reference_densities :: NTuple{2, <:Real} = (1000.0, 1.0),
    )
    phases  = (AqueousPhase(), VaporPhase())
    rho_ref = tuple(Float64.(reference_densities)...)
    reference_phase_index = 1
    return GeothermalTwoPhaseSystem(phases, rho_ref, reference_phase_index, pvt_tables)
end

Base.show(io::IO, ::GeothermalTwoPhaseSystem) =
    print(io, "GeothermalTwoPhaseSystem (AqueousPhase + VaporPhase, 1 component)")

# ── Core MultiPhaseSystem interface ───────────────────────────────────────────

JutulDarcy.get_phases(sys::GeothermalTwoPhaseSystem)          = sys.phases
JutulDarcy.eachphase(sys::GeothermalTwoPhaseSystem)           = (1, 2)
JutulDarcy.number_of_phases(sys::GeothermalTwoPhaseSystem)    = 2
JutulDarcy.number_of_components(sys::GeothermalTwoPhaseSystem) = 1
JutulDarcy.reference_densities(sys::GeothermalTwoPhaseSystem) = sys.rho_ref
JutulDarcy.phase_indices(sys::GeothermalTwoPhaseSystem)       = (1, 2)
JutulDarcy.component_names(sys::GeothermalTwoPhaseSystem)       = (:H₂O,)

# ── Primary variables ─────────────────────────────────────────────────────────
# Only Pressure is selected here.  The Enthalpy primary variable (and all
# (P,H)-dependent secondary variables) are added by _apply_enthalpy_formulation!
# in enthalpy_setup.jl.

function JutulDarcy.select_primary_variables!(
        S, ::GeothermalTwoPhaseSystem, model::SimulationModel,
    )
    S[:Pressure] = Pressure()
    # S[:Enthalpy] = Enthalpy()
end

function Jutul.select_secondary_variables!(S, system::GeothermalTwoPhaseSystem, model)
    JutulDarcy.select_default_darcy_secondary_variables!(S, model.domain, system, model.formulation)
    set_secondary_variables!(model,
        PhaseViscosities = PHDependentPhaseVariableH2O(system.pvt_tables[:viscosity_liquid_ph], system.pvt_tables[:viscosity_vapor_ph]),
        PhaseMassDensities = PHDependentPhaseVariableH2O(system.pvt_tables[:density_liquid_ph], system.pvt_tables[:density_vapor_ph]),
        FluidEnthalpy = PHDependentPhaseVariableH2O(system.pvt_tables[:enthalpy_liquid_ph], system.pvt_tables[:enthalpy_vapor_ph]),
        Saturations = SaturationH2O(system.pvt_tables[:saturation_vapor_ph]),
        Enthalpy = EnthalpyFromPT(system.pvt_tables[:enthalpy]),
    )
end

function JutulDarcy.select_parameters!(prm, system::GeothermalTwoPhaseSystem, model)
    JutulDarcy.select_default_darcy_parameters!(prm, model.domain, system, model.formulation)
    prm[:LiquidMassFractions] = JutulDarcy.PhaseMassFractions(:liquid)
    prm[:VaporMassFractions] = JutulDarcy.PhaseMassFractions(:vapor)
    prm[:Temperature] = JutulDarcy.Temperature()
end

function JutulDarcy.set_reservoir_variable_defaults!(
    model::SimulationModel{O, S, F, C};
        p_min,
        p_max,
        dp_max_abs,
        dp_max_rel,
        ds_max,
        dz_max,
        dr_max,
        dT_max_rel = nothing,
        dT_max_abs = nothing,
        T_min = convert_to_si(0.0, :Celsius),
        flash_reuse_guess = false,
        flash_stability_bypass = flash_reuse_guess
    ) where {O, S<:GeothermalTwoPhaseSystem, F, C}
    # Replace various variables - if they are available
    # replace_variables!(model, OverallMoleFractions = OverallMoleFractions(dz_max = dz_max), throw = false)
    pvt = model.system.pvt_tables
    # Jutul.delete!(model.parameters, :Temperature)

    # replace_variables!(model, Saturations = SaturationsFromEnthalpy(pvt[:saturation_vapor_ph]), throw = false)
    # replace_variables!(model, Temperature = TemperatureFromEnthalpy(pvt[:temperature]), throw = false)
    
    p_def = Pressure(
        max_abs = dp_max_abs,
        max_rel = dp_max_rel,
        minimum = p_min,
        maximum = p_max
    )
    replace_variables!(model, Pressure = p_def, throw = false)
    return model
end

function Jutul.select_equations!(eqs, sys::GeothermalTwoPhaseSystem, model::SimulationModel)
    nc = JutulDarcy.number_of_components(sys)
    mdisc = model.domain.discretizations.mass_flow
    eqs[:mass_conservation] = ConservationLaw(mdisc, :TotalMasses, nc)
    # hdisc = model.domain.discretizations.heat_flow
    # eqs[:energy_conservation] = ConservationLaw(hdisc, :TotalThermalEnergy, nc)
end

# ── Total masses: 1 water component = liquid mass + vapour mass ───────────────
#
# For ImmiscibleSystem:  totmass[ph, i] = ρ_ph * S_ph * V   (2 rows, 2 phases)
# For GeothermalTwoPhaseSystem: totmass[1, i]  = (ρ_l*S_l + ρ_v*S_v) * V  (1 row)
#
# Saturations here come from SaturationsFromEnthalpy (registered as a secondary
# variable by _apply_enthalpy_formulation!).

@jutul_secondary function update_total_masses!(
        totmass,
        tv   :: TotalMasses,
        model :: SimulationModel{G, S},
        PhaseMassDensities,
        Saturations,
        FluidVolume,
        ix,
    ) where {G, S <: GeothermalTwoPhaseSystem}
    rho = PhaseMassDensities
    sat = Saturations
    @inbounds for i in ix
        V = FluidVolume[i]
        totmass[1, i] = (rho[1, i] * sat[1, i] + rho[2, i] * sat[2, i]) * V
    end
end

# ── TotalMasses allocation: 1 row (1 component), not 2 rows (2 phases) ────────
#
# The generic dispatch uses number_of_phases, which gives a 2-row TotalMasses.
# But the mass_conservation ConservationLaw is sized by number_of_components = 1,
# so WellFromFacilityFlowCT gets out[1] vs mix[1:2] → DimensionMismatch at
# sparsity detection.  Override to keep them consistent.

function JutulDarcy.degrees_of_freedom_per_entity(
        model::SimulationModel{G, <:GeothermalTwoPhaseSystem},
        ::TotalMasses,
    ) where {G}
    return 1
end

function JutulDarcy.apply_flow_bc!(
        acc, q,
        bc, model :: SimulationModel{<:Any, <:GeothermalTwoPhaseSystem},
        state, time,
    )
    acc[] += sum(q)
end

# @jutul_secondary function update_total_thermal_energy!(
#         E_total, te::TotalThermalEnergyH2O,
#         model::SimulationModel{<: Any, <: GeothermalTwoPhaseSystem, <: Any, <: Any},
#         Pressure, Enthalpy, Saturations, PhaseMassDensities, TotalMasses,
#         RockDensity, RockInternalEnergy, BulkVolume, FluidVolume, ix
#     )
#     # U_f = FluidInternalEnergy
#     # U_r = RockInternalEnergy
#     # ρ_f = PhaseMassDensities
#     # S = Saturations
#     # V = BulkVolume
#     println("update_total_thermal_energy! with ix = $ix")
#     for i in ix
#         E_i = compute_total_thermal_energy(
#             RockInternalEnergy[i], RockDensity[i], Saturations[:, i], BulkVolume[i], FluidVolume[i],
#             TotalMasses[1, i], Pressure[i], Enthalpy[i], PhaseMassDensities[:, i]
#         )
#         # V_f = FluidVolume[i]
#         # E_i = RockDensity[i]*U_r[i]*(V[i] - V_f)
#         # M = TotalMasses[1, i]
#         # p = Pressure[i]
#         # h = Enthalpy[i]
#         # ρ_l = ρ_f[1, i]
#         # ρ_v = ρ_f[2, i]
#         # ρ_mix = ρ_l*S[1, i] + ρ_v*S[2, i]
#         # E_i += M * (h - p/ρ_mix)
#         # # for ph in axes(U_f, 1)
#         # #     E_i += ρ_f[ph, i]*S[ph, i]*U_f[ph, i]*V_f
#         # # end
#         E_total[i] = E_i
#     end
# end

# const MSWellFlowModelGeothermal = SimulationModel{<:JutulDarcy.MSWellDomain, <:GeothermalTwoPhaseSystem}
# @jutul_secondary function update_total_thermal_energy!(E_total, te::TotalThermalEnergyH2O, model::MSWellFlowModelGeothermal,
#     Pressure, Enthalpy, Saturations, PhaseMassDensities, TotalMasses,
#     MaterialDensities, MaterialInternalEnergy, BulkVolume, FluidVolume, ix)

#     for i in ix
#         E_i = compute_total_thermal_energy(
#             MaterialInternalEnergy[i], MaterialDensities[i], Saturations[:, i], BulkVolume[i], FluidVolume[i],
#             TotalMasses[1, i], Pressure[i], Enthalpy[i], PhaseMassDensities[:, i]
#         )
#         # V_f = FluidVolume[i]
#         # E_i = RockDensity[i]*U_r[i]*(V[i] - V_f)
#         # M = TotalMasses[1, i]
#         # p = Pressure[i]
#         # h = Enthalpy[i]
#         # ρ_l = ρ_f[1, i]
#         # ρ_v = ρ_f[2, i]
#         # ρ_mix = ρ_l*S[1, i] + ρ_v*S[2, i]
#         # E_i += M * (h - p/ρ_mix)
#         # # for ph in axes(U_f, 1)
#         # #     E_i += ρ_f[ph, i]*S[ph, i]*U_f[ph, i]*V_f
#         # # end
#         E_total[i] = E_i
#     end
#     # update_total_thermal_energy!(E_total, te::TotalThermalEnergyH2O, nothing,
#     # Pressure, Enthalpy, Saturations, PhaseMassDensities, TotalMasses,
#     # MaterialDensities, MaterialInternalEnergy, BulkVolume, FluidVolume, ix)
# end

# function compute_total_thermal_energy(U_r, ρ_r, S, V, V_f, M_f, p, h, ρ_f)

#     ρ_mix = ρ_f[1]*S[1] + ρ_f[2]*S[2]
#     return ρ_r*U_r*(V - V_f) + M_f * (h - p/ρ_mix)

# end

# ── Multi-segment well perforation flux ───────────────────────────────────────
#
# For ImmiscibleSystem, out[ph] = phase_mass_flux(ph).
# For GeothermalTwoPhaseSystem, out has 1 entry: total water mass flux.
#
# Note: state_well.Saturations must have 2 rows, which requires
# _apply_enthalpy_formulation! to be called with add_phase_split = true for
# well models (done in Phase 3 / setup_reservoir_model_geothermal_2ph).

# Base.@propagate_inbounds function JutulDarcy.multisegment_well_perforation_flux!(
#         out,
#         sys       :: GeothermalTwoPhaseSystem,
#         state_res,
#         state_well,
#         rhoS,
#         conn,
#     )
#     λ_l, λ_v = JutulDarcy.perforation_reservoir_mobilities(
#         state_res, state_well, sys, conn.reservoir, conn.well,
#     )
    
#     # μ = state_well[:PhaseViscosities][:, conn.well]
#     # λ_l = 1/μ[1]
#     # λ_v = 1/μ[2]
#     λ_t = λ_l + λ_v
#     # λ_t = sum(JutulDarcy.perforation_reservoir_mobilities(
#     #     state_res, state_well, sys, conn.reservoir, conn.well,
#     # ))
#     q_total = zero(eltype(out))
#     for ph in 1:2
#         q_total += JutulDarcy.perforation_phase_mass_flux(
#             λ_t, conn, state_res, state_well, ph,
#         )
#     end
#     if q_total < 0
#         @info "multisegment_well_perforation_flux!: computed total mass flux = $(value(q_total))"
#     end
#     out[] = q_total
#     return out
# end

# ── Simple-well perforation flux ─────────────────────────────────────────────
#
# The SimpleWell (StandardWellFlowModel) model uses MassFractions and a
# simplified mass-balance.  For GeothermalTwoPhaseSystem (nc=1), out[1] is the
# total water mass flux summed over phases.
#
# Injection direction is determined from the liquid phase (phase 1) potential.
# For production both phase fluxes are accumulated.

# Base.@propagate_inbounds function JutulDarcy.simple_well_perforation_flux!(
#         out,
#         sys       :: GeothermalTwoPhaseSystem,
#         state_res,
#         state_well,
#         rhoS,
#         conn,
#     )
#     error()
#     rc  = conn.reservoir
#     ρ   = state_res.PhaseMassDensities
#     mob = state_res.PhaseMobilities

#     # Total mass mobility in the reservoir cell
#     ρλ_t = ρ[1, rc] * mob[1, rc] + ρ[2, rc] * mob[2, rc]

#     # Use liquid phase (1) to determine flow direction
#     ψ₁ = JutulDarcy.perforation_phase_potential_difference(
#         conn, state_res, state_well, 1,
#     )

#     q_total = zero(eltype(out))
#     if ψ₁ < 0
#         # Injection: single-component water is injected; use liquid-phase
#         # pressure drive and total reservoir mass mobility.
#         q_total = ψ₁ * ρλ_t
#     else
#         # Production: accumulate Darcy flux from each phase independently.
#         ψ₂ = JutulDarcy.perforation_phase_potential_difference(
#             conn, state_res, state_well, 2,
#         )
#         q_total = mob[1, rc] * ρ[1, rc] * ψ₁ +
#                   mob[2, rc] * ρ[2, rc] * ψ₂
#     end
#     out[1] = q_total
#     return out
# end

# ── SurfaceWellConditions: phase densities and saturations from property evaluators ──
#
# The generic flash_wellstream_at_surface dispatches on system type to compute
# surface densities (rhoS) and volume fractions. For GeothermalTwoPhaseSystem
# the (P,H)-based property evaluators (PhaseMassDensities, Saturations) already
# provide the correct per-phase values — no additional flash calculation is needed.
# We read them directly from the well state at the top node.

function JutulDarcy.flash_wellstream_at_surface(
        var,
        well_model,
        system  :: GeothermalTwoPhaseSystem,
        well_state,
        rhoS,
        cond    = JutulDarcy.default_surface_cond(),
    )
    wc  = JutulDarcy.well_top_node()
    rho = well_state.PhaseMassDensities[:, wc]   # SVector-like, length 2
    sat = well_state.Saturations[:, wc]           # volume fractions, length 2
    # Guard against all-zero saturation (e.g. during AD sparsity detection)
    s_total = sum(sat)
    volfrac = s_total > 0 ? sat ./ s_total : sat .* 0 .+ 0.5
    return (rho, volfrac)
end

# ── Convergence criterion for GeothermalTwoPhaseSystem ────────────────────────
#
# The generic multiphase convergence_criterion loops `ph = 1:number_of_phases`
# over the residual `r`, but for GeothermalTwoPhaseSystem `r` has only 1 row
# (nc = 1, single H₂O component). Reading r[2, c] is a silent memory read of
# the wrong cell, producing a spurious "Vapor" mass balance residual normalised
# by the (low) vapor density — which never converges in superheated cases (:c).
#
# Fix: use Val(1) and normalize by the mixture density ρ_l*S_l + ρ_v*S_v.

function JutulDarcy.convergence_criterion(
        model::SimulationModel{D, <:GeothermalTwoPhaseSystem},
        storage, eq::ConservationLaw{:TotalMasses}, eq_s, r;
        dt = 1.0, update_report = missing,
    ) where D
    M = global_map(model.domain)
    v = x -> as_value(Jutul.active_view(x, M, for_variables = false))
    Φ    = v(storage.state.FluidVolume)
    ρ_ph = v(storage.state.PhaseMassDensities)
    S_ph = v(storage.state.Saturations)
    nc   = length(Φ)
    ρ_mix = reshape([ρ_ph[1,c]*S_ph[1,c] + ρ_ph[2,c]*S_ph[2,c] for c in 1:nc], 1, nc)
    cnv, mb = JutulDarcy.cnv_mb_errors(r, Φ, ρ_mix, dt, Val(1))
    dp_abs, dp_rel = JutulDarcy.pressure_increments(model, storage.state, update_report)
    names = (:Water,)
    return (
        CNV = (errors = cnv, names = names),
        MB  = (errors = mb,  names = names),
        increment_dp_abs = (errors = (dp_abs/1e6,), names = (raw"Δp (abs, MPa)",)),
        increment_dp_rel = (errors = (dp_rel,),     names = (raw"Δp (rel)",)),
    )
end

function JutulDarcy.temperature_increment(model::SimulationModel{D, <:GeothermalTwoPhaseSystem}, state, update_report) where D
    return 0.0
end

# function JutulDarcy.temperature_increment(model::SimulationModel{D, <:GeothermalTwoPhaseSystem}, state, update_report) where D
#     return 0.0
# end

function JutulDarcy.temperature_increment(model::SimulationModel{D, <:GeothermalTwoPhaseSystem}, state, update_report::Missing) where D
    return 1.0
end