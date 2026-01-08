function injected_fluid_volume(rate, duration)
    Q, t = rate, duration
    return Q*t
end

function fluid_radius_aquifer(rate, duration, height, porosity)
    Q, t = rate, duration
    H, ϕ = height, porosity
    V = injected_fluid_volume(Q, t)
    return sqrt(V/(π*H*ϕ))
end

function aquifer_heat_capacity(porosity, fluid_heat_capacity, rock_heat_capacity)
    ϕ, Cf, Cr = porosity, fluid_heat_capacity, rock_heat_capacity
    return ϕ*Cf + (1 - ϕ)*Cr
end

function injected_thermal_volume(rate, duration, porosity, fluid_heat_capacity, rock_heat_capacity)
    Q, t = rate, duration
    ϕ, Cf, Cr = porosity, fluid_heat_capacity, rock_heat_capacity
    Ca = aquifer_heat_capacity(ϕ, Cf, Cr)
    Vf = injected_fluid_volume(Q, t)
    return Vf*Cf/Ca
end

function thermal_radius_aquifer(rate, duration, height, porosity, fluid_heat_capacity, rock_heat_capacity)
    Q, t = rate, duration
    H, ϕ, Cf, Cr = height, porosity, fluid_heat_capacity, rock_heat_capacity
    V = injected_thermal_volume(Q, t, ϕ, Cf, Cr)
    return sqrt(V/(π*H))
end

function injection_rate_from_thermal_radius(thermal_radius, duration, height, porosity, fluid_heat_capacity, rock_heat_capacity)
    r, H, ϕ, Cf, Cr = thermal_radius, height, porosity, fluid_heat_capacity, rock_heat_capacity
    VQ⁻¹ = injected_thermal_volume(1.0, duration, ϕ, Cf, Cr)
    return π*r^2*H/VQ⁻¹
end