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

end

JutulDarcy.system_uses_cnv_mb(system::GeothermalTwoPhaseSystem) = true
