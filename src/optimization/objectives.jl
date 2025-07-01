function well_mismatch_thermal(
    model,
    wells,
    states_f,
    state_c,
    dt,
    step_info,
    forces;
    w_bhp = 1/(10.0si_unit(:bar)),
    w_energy = 1/(si_unit(:mega)*si_unit(:watt)),
    w_temp = 1/(si_unit(:Kelvin)),
    scale = 1.0
)

    step_no = step_info[:step]
    state_f = states_f[step_no]
    get_property = (state, well, prop) -> state[well][prop][1]
    mismatch = 0.0

    for well in wells

        if forces[:Facility].control[well] isa DisabledControl
            continue
        end
        
        if w_bhp > 0
            bhp_f = get_property(state_f, well, :Pressure)
            bhp_c = get_property(state_c, well, :Pressure)
            Δbhp = bhp_f - bhp_c
            mismatch += (w_bhp*Δbhp)^2*dt
        end

        if w_temp > 0
            temp_f = get_property(state_f, well, :Temperature)
            temp_c = get_property(state_c, well, :Temperature)
            Δtemp = temp_f - temp_c
            mismatch += (w_temp*Δtemp)^2*dt
        end

        if w_energy > 0
            rate_f = get_property(state_f, well, :TotalMassFlux)
            rate_c = get_property(state_c, well, :TotalMassFlux)
            h_f = get_property(state_f, well, :FluidEnthalpy)
            h_c = get_property(state_c, well, :FluidEnthalpy)
            energy_f = rate_f*h_f
            energy_c = rate_c*h_c
            Δenergy = energy_f - energy_c
            mismatch += (w_energy*Δenergy)^2*dt
        end

    end

    mismatch = mismatch./scale

    return mismatch

end

function reservoir_mismatch_thermal(
    model,
    states_f,
    state_c,
    dt,
    step_no,
    forces;
    w_pressure = 1/(1.0si_unit(:bar)),
    w_temperature = 1/(1.0si_unit(:Kelvin)),
    w_energy = /(1.0si_unit(:mega)*si_unit(:watt)),
    scale = 1.0
    )

    state_f = states_f[step_no][:Reservoir]
    mismatch = 0.0

    state_c = state_c[:Reservoir]

    weights = [w_pressure, w_temperature, w_energy]
    names = [:Pressure, :Temperature, :TotalThermalEnergy]
    for (n,ω) in zip(names, weights)
        if ω > 0
            Δ = state_f[n][1] - state_c[n][1]
            mismatch += (ω*Δ)^2*dt
        end
    end

    mismatch = mismatch./scale

    return mismatch

end