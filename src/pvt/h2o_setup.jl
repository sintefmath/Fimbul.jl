"""
    JutulDarcy.add_thermal_to_model!(model::GeothermalModel)

Activate the pressure-enthalpy thermal formulation for an `H2OSystem` model.

This replaces `Temperature` as a primary variable with `Enthalpy`, registers
temperature and phase enthalpy as secondary variables derived from the
pressure-enthalpy tables, and stores the `(P, T) -> H` lookup table in
`model.extra[:enthalpy]` for later use.
"""
function JutulDarcy.add_thermal_to_model!(model::GeothermalModel)

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

    unique!(push!(model.output_variables, :Enthalpy))

    return model

end

JutulDarcy.system_uses_cnv_mb(system::H2OSystem) = true
