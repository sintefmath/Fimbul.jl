const WATER_CRITICAL_PRESSURE = 22.064e6
const WATER_CRITICAL_TEMPERATURE = 647.096
const WATER_CRITICAL_ENTHALPY = 2085e3

"""
    H2OSystem

Two-phase (liquid + vapour), single-component fluid system for geothermal
(P, H) simulations.

| quantity               | value              |
|:-----------------------|:-------------------|
| `number_of_phases`     | 2                  |
| `number_of_components` | 1  (pure H₂O)      |

One mass conservation equation is closed by the `Enthalpy` primary variable
(alongside `Pressure`).  Phase split, densities, viscosities and enthalpies are
obtained via Steam-Tables look-ups from `build_steam_tables_2ph`.

Use `H2OSystem` together with `setup_reservoir_model` to obtain a fully
configured model.
"""
struct H2OSystem{T <: Tuple, F <: NTuple} <: JutulDarcy.AbstractCompositionalSystemLV
    phases  :: T
    rho_ref :: F
    reference_phase_index::Int
    pvt_tables::Dict
end

"""
    H2OSystem(; reference_densities = (1000.0, 1.0))

Construct a `H2OSystem` for a `(AqueousPhase, VaporPhase)` pair.
`reference_densities` should be `(ρ_liquid_ref, ρ_vapour_ref)` in kg/m³.
"""
function H2OSystem(pvt_tables::Dict;
        reference_densities :: NTuple{2, <:Real} = (1000.0, 1.0),
    )
    phases  = (AqueousPhase(), VaporPhase())
    rho_ref = tuple(Float64.(reference_densities)...)
    reference_phase_index = 1
    return H2OSystem(phases, rho_ref, reference_phase_index, pvt_tables)
end

Base.show(io::IO, ::H2OSystem) =
    print(io, "H2OSystem (AqueousPhase + VaporPhase, 1 component)")

# Core functionality
JutulDarcy.get_phases(sys::H2OSystem)           = sys.phases
JutulDarcy.eachphase(sys::H2OSystem)            = (1, 2)
JutulDarcy.number_of_phases(sys::H2OSystem)     = 2
JutulDarcy.number_of_components(sys::H2OSystem) = 1
JutulDarcy.reference_densities(sys::H2OSystem)  = sys.rho_ref
JutulDarcy.phase_indices(sys::H2OSystem)        = (1, 2)
JutulDarcy.component_names(sys::H2OSystem)      = (:H₂O,)

# Convenience const for easy dispatching on models with this H2OSystem.
const GeothermalModel = SimulationModel{<:Any, <:H2OSystem, <:Any, <:Any}
Jutul.default_value(model::GeothermalModel, ::JutulDarcy.PhaseMassFractions) = 1.0

# Only Pressure is selected here. The Enthalpy primary variable is added by add_thermal_to_model!.
function JutulDarcy.select_primary_variables!(S, ::H2OSystem, model::SimulationModel)
    S[:Pressure] = Pressure()
end

# Secondary variables, including the enthalpy from pressure/temperature
function Jutul.select_secondary_variables!(S, system::H2OSystem, model)
    JutulDarcy.select_default_darcy_secondary_variables!(S, model.domain, system, model.formulation)
    tab = system.pvt_tables
    set_secondary_variables!(model,
        PhaseViscosities = PHDependentPhaseVariableH2O(tab[:viscosity_liquid_ph], tab[:viscosity_vapor_ph]),
        PhaseMassDensities = PHDependentPhaseVariableH2O(tab[:density_liquid_ph], tab[:density_vapor_ph]),
        FluidEnthalpy = PHDependentPhaseVariableH2O(tab[:enthalpy_liquid_ph], tab[:enthalpy_vapor_ph]),
        Saturations = SaturationH2O(tab[:saturation_vapor_ph]),
        Enthalpy = EnthalpyFromPT(tab[:enthalpy]),
    )
end

# Liquid and vapor mass fractions are included to leverage the JutulDarcy
# functionality for compositional systems. 
function JutulDarcy.select_parameters!(prm, system::H2OSystem, model)
    JutulDarcy.select_default_darcy_parameters!(prm, model.domain, system, model.formulation)
    prm[:LiquidMassFractions] = JutulDarcy.PhaseMassFractions(:liquid)
    prm[:VaporMassFractions] = JutulDarcy.PhaseMassFractions(:vapor)
    prm[:Temperature] = JutulDarcy.Temperature()
end

# TODO: This overload is necessary to avoid replacing temperature and
# saturation, which are already registered as secondary variables during
# add_thermal_to_model!, and should be refactored to avoid the need for this
# special case.
function JutulDarcy.set_reservoir_variable_defaults!(model::GeothermalModel;
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
)
    p_def = Pressure(
        max_abs = dp_max_abs,
        max_rel = dp_max_rel,
        minimum = p_min,
        maximum = p_max
    )
    replace_variables!(model, Pressure = p_def, throw = false)
    return model
end

# TODO: This overload can likely be avoided by adding missing compositional
# logic to H2OSystem
@jutul_secondary function update_total_masses!(
        totmass,
        tv::TotalMasses,
        model::GeothermalModel,
        PhaseMassDensities,
        Saturations,
        FluidVolume,
        ix,
    )
    rho = PhaseMassDensities
    sat = Saturations
    @inbounds for i in ix
        V = FluidVolume[i]
        totmass[1, i] = (rho[1, i] * sat[1, i] + rho[2, i] * sat[2, i]) * V
    end
end

function JutulDarcy.degrees_of_freedom_per_entity(
    model::GeothermalModel, ::TotalMasses
    )
    return 1
end

# For H2OSystem the (P,H)-based property evaluators (PhaseMassDensities,
# Saturations) already provide the correct per-phase values — no additional
# flash calculation is needed. We read them directly from the well state at the
# top node.
function JutulDarcy.flash_wellstream_at_surface(
        var,
        well_model,
    system  :: H2OSystem,
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

"""
    JutulDarcy.set_default_cnv_mb_inner!(tol, model::GeothermalModel;
        inc_tol_dh_rel = Inf, inc_tol_dh_abs = Inf, kwargs...)

Extend the default convergence-tolerance dictionary for an `H2OSystem` model
with enthalpy-increment tolerances for the thermal energy equation.
"""
function JutulDarcy.set_default_cnv_mb_inner!(tol, model::GeothermalModel;
    inc_tol_dh_rel = Inf, inc_tol_dh_abs = Inf, kwargs...
    )
    invoke(JutulDarcy.set_default_cnv_mb_inner!, Tuple{Dict, SimulationModel}, tol, model; kwargs...)
    if haskey(model.equations, :energy_conservation)
        etol = tol[:energy_conservation]
        tol[:energy_conservation] = merge(etol, (
            increment_dh_abs = inc_tol_dh_abs,
            increment_dh_rel = inc_tol_dh_rel,
        ))
    end
end

"""
    Jutul.convergence_criterion(model::GeothermalModel, storage,
        eq::ConservationLaw{:TotalMasses}, eq_s, r; dt = 1.0,
        update_report = missing)

Compute mass-balance convergence metrics for an `H2OSystem` model.

The residual is normalized using the mixture density `ρ_l S_l + ρ_v S_v`
"""
function Jutul.convergence_criterion(
    model::SimulationModel{D, S}, storage,
    eq::ConservationLaw{:TotalMasses}, eq_s, r;
    dt = 1.0, update_report = missing,
    ) where {D, S <: H2OSystem}
    M = global_map(model.domain)
    v = x -> as_value(Jutul.active_view(x, M, for_variables = false))
    Φ    = v(storage.state.FluidVolume)
    ρ_ph = v(storage.state.PhaseMassDensities)
    S_ph = v(storage.state.Saturations)
    nc   = length(Φ)
    ρ_mix = reshape([ρ_ph[1,c]*S_ph[1,c] + ρ_ph[2,c]*S_ph[2,c] for c in 1:nc], 1, nc)
    cnv, mb = JutulDarcy.cnv_mb_errors(r, Φ, ρ_mix, dt, Val(1))
    dp_abs, dp_rel = JutulDarcy.pressure_increments(model, storage.state, update_report)
    names = JutulDarcy.component_names(model.system)
    return (
        CNV = (errors = cnv, names = names),
        MB  = (errors = mb,  names = names),
        increment_dp_abs = (errors = (dp_abs/1e6,), names = (raw"Δp (abs, MPa)",)),
        increment_dp_rel = (errors = (dp_rel,),     names = (raw"Δp (rel)",)),
    )
end

"""
    Jutul.convergence_criterion(model::GeothermalModel, storage,
        eq::ConservationLaw{:TotalThermalEnergy}, eq_s, r; dt = 1.0,
        update_report = missing)

Compute convergence metrics for the thermal energy equation in an `H2OSystem`
model, including absolute and relative enthalpy increments.
"""
function Jutul.convergence_criterion(
    model::GeothermalModel, storage,
    eq::ConservationLaw{:TotalThermalEnergy}, eq_s, r;
    dt = 1.0, update_report = missing
    )
    a = active_entities(model.domain, Cells())
    E0 = storage.state0.TotalThermalEnergy
    eb, cnv, Etot = 0.0, -Inf, 0.0
    Δh_abs, Δh_rel = enthalpy_increment(model, storage.state, update_report)
    for (i, c) in enumerate(a)
        eb += r[i]
        cnv = max(cnv, abs(r[i])*dt/value(E0[c]))
        Etot += value(E0[c])
    end
    eb = abs(eb)*dt/Etot

    return (
        CNV = (errors = (cnv, ), names = ("Max", )),
        EB = (errors = (eb, ), names = ("Energy balance", )), 
        increment_dh_abs = (errors = (Δh_abs, ), names = ("Δh", )),
        increment_dh_rel = (errors = (Δh_rel, ), names = ("Δh (rel)", )),
        )
end

# Enthalpy increment calculation for convergence criterion
function enthalpy_increment(model, state, update_report)
    max_h = maximum(value, state.Enthalpy)
    h_report = update_report[:Enthalpy]
    if haskey(h_report, :max)
        dh_abs = h_report.max
        dh_rel = dh_abs/max_h
    else
        dh_abs = 1.0
        dh_rel = 1.0
    end
    return (dh_abs, dh_rel)
end

function enthalpy_increment(model, state, update_report::Missing)
    return (1.0, 1.0)
end