using SpecialFunctions

function radial_analytical(r, t, Q, T_i, T_0, α_L, b, ϕ, D, R)

    A = Q/(2*π*b*ϕ*R)
    r_f = sqrt(2*A*t)

    den = 2*sqrt((4/3)*α_L*r_f^3 + (D/(A*R))*r_f^4)
    ξ = (r^2 - r_f^2)/den

    T = T_0 + (T_i - T_0)*0.5*erfc(ξ)

    return T

end

function radial_analtical(r, t, Q, T_i, T_0, ρ_f, α_L, b, domain::DataDomain)

    ϕ = domain[:porosity][1]
    ρ_r = domain[:rock_density][1]
    C_f = domain[:component_heat_capacity][1]
    C_r = domain[:rock_heat_capacity][1]
    λ_f = domain[:fluid_thermal_conductivity][1]
    λ_r = domain[:rock_thermal_conductivity][1]
    
    C_b = ϕ*C_f + (1-ϕ)*C_r
    ρ_b = ρ_f*ϕ + ρ_r*(1-ϕ)
    K_DT = C_b/(ρ_f*C_f)
    λ_b = ϕ*λ_f + (1-ϕ)*λ_r
    R = 1 + (ρ_b/ϕ)*K_DT
    D = λ_b/(ϕ*ρ_f*C_f)

    return radial_analytical.(r, t, Q, T_i, T_0, α_L, b, ϕ, D, R)

end