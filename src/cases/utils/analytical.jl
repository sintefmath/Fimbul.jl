using SpecialFunctions
using Statistics

function radial_analytical(r, t, Q, T_i, T_0, α_L, b, ϕ, D, R)

    A = Q/(2*π*b*ϕ*R)
    r_f = sqrt(2*A*t)

    den = 2*sqrt((4/3)*α_L*r_f^3 + (D/(A*R))*r_f^4)
    ξ = (r^2 - r_f^2)/den

    T = T_0 + (T_i - T_0)*0.5*erfc(ξ)

    return T

end

function representative_domain_value(domain::DataDomain, key::Symbol; rtol = 1e-6)
    raw = domain[key]
    vals = vec(collect(raw))
    isempty(vals) && error("Domain property $(key) is empty.")
    μ = mean(vals)
    scale = max(abs(μ), eps(Float64))
    maxdev = maximum(abs.(vals .- μ))
    if maxdev > rtol*scale
        @warn "Domain property $(key) is not constant; using mean value $(μ)." maxdev = maxdev mean = μ
    end
    return μ
end

function radial_analytical(r, t, Q, T_i, T_0, ρ_f, α_L, b, domain::DataDomain)
    ϕ = representative_domain_value(domain, :porosity)
    ρ_r = representative_domain_value(domain, :rock_density)
    C_f = representative_domain_value(domain, :component_heat_capacity)
    C_r = representative_domain_value(domain, :rock_heat_capacity)
    λ_f = representative_domain_value(domain, :fluid_thermal_conductivity)
    λ_r = representative_domain_value(domain, :rock_thermal_conductivity)

    ρ_bC_b = ϕ*ρ_f*C_f + (1 - ϕ)*ρ_r*C_r
    λ_b = ϕ*λ_f + (1 - ϕ)*λ_r
    R = ρ_bC_b/(ϕ*ρ_f*C_f)
    D = λ_b/(ϕ*ρ_f*C_f)

    return radial_analytical.(r, t, Q, T_i, T_0, α_L, b, ϕ, D, R)
end