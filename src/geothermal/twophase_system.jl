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
struct GeothermalTwoPhaseSystem{T <: Tuple, F <: NTuple} <: MultiPhaseSystem
    phases  :: T
    rho_ref :: F
    reference_phase_index::Int
end

"""
    GeothermalTwoPhaseSystem(; reference_densities = (1000.0, 1.0))

Construct a `GeothermalTwoPhaseSystem` for a `(AqueousPhase, VaporPhase)` pair.
`reference_densities` should be `(ρ_liquid_ref, ρ_vapour_ref)` in kg/m³.
"""
function GeothermalTwoPhaseSystem(;
        reference_densities :: NTuple{2, <:Real} = (1000.0, 1.0),
    )
    phases  = (AqueousPhase(), VaporPhase())
    rho_ref = tuple(Float64.(reference_densities)...)
    reference_phase_index = 1
    return GeothermalTwoPhaseSystem(phases, rho_ref, reference_phase_index)
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

# ── Darcy face flux: 1-component water = sum of liquid and vapour fluxes ──────
#
# The generic component_mass_fluxes! (for ImmiscibleSystem) copies phase_fluxes[ph]
# into q[ph].  For GeothermalTwoPhaseSystem we sum both phase fluxes into q[1].

@inline function component_mass_fluxes!(
        q, face, state,
        model    :: SimulationModel{<:Any, <:GeothermalTwoPhaseSystem},
        flux_type, kgrad, upw,
    )
    phase_fluxes = JutulDarcy.darcy_phase_mass_fluxes(
        face, state, model, flux_type, kgrad, upw,
    )
    return setindex(q, phase_fluxes[1] + phase_fluxes[2], 1)
end

# ── Boundary-condition mass flux: sum phases into single component ────────────
#
# compute_bc_mass_fluxes returns a nph-length vector q (one entry per phase).
# For ImmiscibleSystem apply_flow_bc! does acc[ph] += q[ph].
# For GeothermalTwoPhaseSystem, acc has 1 row so we sum all phases.

function JutulDarcy.apply_flow_bc!(
        acc, q,
        bc, model :: SimulationModel{<:Any, <:GeothermalTwoPhaseSystem},
        state, time,
    )
    q_total = zero(eltype(q))
    for i in eachindex(q)
        q_total += q[i]
    end
    acc[1] += q_total
end

# ── Multi-segment well perforation flux ───────────────────────────────────────
#
# For ImmiscibleSystem, out[ph] = phase_mass_flux(ph).
# For GeothermalTwoPhaseSystem, out has 1 entry: total water mass flux.
#
# Note: state_well.Saturations must have 2 rows, which requires
# _apply_enthalpy_formulation! to be called with add_phase_split = true for
# well models (done in Phase 3 / setup_reservoir_model_geothermal_2ph).

Base.@propagate_inbounds function JutulDarcy.multisegment_well_perforation_flux!(
        out,
        sys       :: GeothermalTwoPhaseSystem,
        state_res,
        state_well,
        rhoS,
        conn,
    )
    λ_t = sum(JutulDarcy.perforation_reservoir_mobilities(
        state_res, state_well, sys, conn.reservoir, conn.well,
    ))
    q_total = zero(eltype(out))
    for ph in 1:2
        q_total += JutulDarcy.perforation_phase_mass_flux(
            λ_t, conn, state_res, state_well, ph,
        )
    end
    out[1] = q_total
    return out
end

# ── Simple-well perforation flux ─────────────────────────────────────────────
#
# The SimpleWell (StandardWellFlowModel) model uses MassFractions and a
# simplified mass-balance.  For GeothermalTwoPhaseSystem (nc=1), out[1] is the
# total water mass flux summed over phases.
#
# Injection direction is determined from the liquid phase (phase 1) potential.
# For production both phase fluxes are accumulated.

Base.@propagate_inbounds function JutulDarcy.simple_well_perforation_flux!(
        out,
        sys       :: GeothermalTwoPhaseSystem,
        state_res,
        state_well,
        rhoS,
        conn,
    )
    rc  = conn.reservoir
    ρ   = state_res.PhaseMassDensities
    mob = state_res.PhaseMobilities

    # Total mass mobility in the reservoir cell
    ρλ_t = ρ[1, rc] * mob[1, rc] + ρ[2, rc] * mob[2, rc]

    # Use liquid phase (1) to determine flow direction
    ψ₁ = JutulDarcy.perforation_phase_potential_difference(
        conn, state_res, state_well, 1,
    )

    q_total = zero(eltype(out))
    if ψ₁ < 0
        # Injection: single-component water is injected; use liquid-phase
        # pressure drive and total reservoir mass mobility.
        q_total = ψ₁ * ρλ_t
    else
        # Production: accumulate Darcy flux from each phase independently.
        ψ₂ = JutulDarcy.perforation_phase_potential_difference(
            conn, state_res, state_well, 2,
        )
        q_total = mob[1, rc] * ρ[1, rc] * ψ₁ +
                  mob[2, rc] * ρ[2, rc] * ψ₂
    end
    out[1] = q_total
    return out
end

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