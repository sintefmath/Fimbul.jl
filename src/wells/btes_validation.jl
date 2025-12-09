using Integrals

function temperature_pipe(u, Tin, Trock, ρ, Cp, A, L, R1Δ, R2Δ, R12Δ)

    f1, f2, f3, f4, f5 = f_functions(u, ρ, Cp, A, R1Δ, R2Δ, R12Δ)

    domain = z -> (0.0, z)
    
    integrand_supp = z -> (ξ,p) -> Trock*f4(z - ξ)
    I_supp = z -> IntegralProblem(integrand_supp(z), domain(z))
    
    integrand_ret = z -> (ξ,p) -> Trock*f5(z - ξ)
    I_ret = z -> IntegralProblem(integrand_ret(z), domain(z))

    # Find T_out from T_supply(L) = T_return(L)
    Tout = (
        Tin*(f1(L) .+ f2(L)) .+ 
        solve(I_supp(L), QuadGKJL()) .+ solve(I_ret(L), QuadGKJL())
        )/(f3(L) .- f2(L))

    T_supply = z -> Tin*f1(z) .+ Tout*f2(z) .+ solve(I_supp(z), QuadGKJL())
    T_return = z -> -Tin*f2(z) .+ Tout*f3(z) .- solve(I_ret(z), QuadGKJL())

    return T_supply, T_return

end

function f_functions(u, ρ, Cp, A, R1Δ, R2Δ, R12Δ)

    β1 = 1/(R1Δ*A*ρ*Cp*u)
    β2 = 1/(R2Δ*A*ρ*Cp*u)
    β12 = 1/(R12Δ*A*ρ*Cp*u)
    β = (β2 - β1)/2
    γ = sqrt((β1 + β2)^2/4 + β12*(β1 + β2))
    δ = 1/γ*(β12 + (β1 + β2)/2)

    f1 = z -> exp(β*z)*(cosh(γ*z) - δ*sinh(γ*z))
    f2 = z -> exp(β*z)*β12/γ*sinh(γ*z)
    f3 = z -> exp(β*z)*(cosh(γ*z) + δ*sinh(γ*z))
    f4 = z -> exp(β*z)*(β1*cosh(γ*z) - (δ*β1 + β2*β12/γ)*sinh(γ*z))
    f5 = z -> exp(β*z)*(β2*cosh(γ*z) + (δ*β2 + β1*β12/γ)*sinh(γ*z))
    return (f1, f2, f3, f4, f5)
end

function temperature_u1(Q, T_in, T_rock, ρ, Cp,
    length, radius_grout, radius_pipe, wall_thickness_pipe, pipe_spacing,
    thermal_conductivity_grout, thermal_conductivity_pipe)

    Rpg, Rgr, Rgg = Fimbul.closed_loop_thermal_resistance_u1(
        radius_grout, radius_pipe, wall_thickness_pipe, pipe_spacing,
        thermal_conductivity_grout, thermal_conductivity_pipe)
    A = π*(radius_pipe - wall_thickness_pipe)^2

    u = Q/A
    Tps, Tpr = temperature_pipe_u1(
        u, T_in, T_rock, ρ, Cp, A, length, Rpg, Rgr, Rgg)
    
    Tgs, Tgr = temperature_grout_u1(
        Tps, Tpr, T_rock, Rpg, Rgr, Rgg)
    
    sol = Dict()
    sol[:pipe_supply] = Tps
    sol[:pipe_return] = Tpr
    sol[:grout_supply] = Tgs
    sol[:grout_return] = Tgr

    return sol

end

function temperature_pipe_u1(u, T_in, T_rock, ρ, Cp, A, L, Rpg, Rgr, Rgg)

    u1 = 1/Rpg + 1/Rgr + 1/Rgg
    R1Δ = Rpg + Rgr
    R2Δ = Rpg + Rgr
    R12Δ = ((u1*Rpg*Rgg)^2 - Rpg^2)/Rgg

    return temperature_pipe(u, T_in, T_rock, ρ, Cp, A, L, R1Δ, R2Δ, R12Δ)

end

function temperature_grout_u1(T_supply, T_return, T_rock, Rpg, Rgr, Rgg)

    u1 = 1/Rpg + 1/Rgr + 1/Rgg
    T_grout_supply = z -> (
        T_rock/Rgr + 
        T_return(z)/Rpg + 
        (T_rock/Rgr + T_supply(z)/Rpg)*u1*Rgg
        )*Rgg/
        ((Rgg*u1)^2-1)

    T_grout_return = z -> (
        T_grout_supply(z)/Rgg + 
        T_return(z)/Rpg + 
        T_rock/Rgr
        )*1/u1

    return T_grout_supply, T_grout_return

end

Rpg = 0.15577 # pipe to grout
Rgg = 0.11516  # grout to grout (between pipes)
Rgr = 0.02574  # grout to rock
u1 = 1/Rpg + 1/Rgr + 1/Rgg
R1Δ = R2Δ = Rpg + Rgr
R12Δ = ((u1*Rpg*Rgg)^2 - Rpg^2)/Rgg

d = 32e-3
Ai = π*(d/2 - 2.9e-3)^2           # internal area of one pipe [m²]
ρf = 988.1
Cf = 4.1312e6/ρf               # J m⁻3 K⁻1

using Jutul
meter, day = si_units(:meter, :day)
Q = 21.86*meter^3/day          # m³ s⁻1
u = Q/(Ai)                 # m s⁻1

sol = temperature_u1(Q, 80.0 + 273.15, 10.0 + 273.15, ρf, Cf,
    100.0, 65e-3, 16e-3, 2.9e-3, 60e-3,
    2.3, 0.38)

# z = range(0, stop=100.0, length=10)
# f = f_functions(Ai, ρf, Cf, u, R1Δ, R2Δ, R12Δ)
# f_val = [fi.(z) for fi in f]

# T_supply, T_return = temperature_pipe(u, 80.0, 10.0,
#     ρf, Cf, Ai, 100.0, R1Δ, R2Δ, R12Δ)

# T_supply, T_return = temperature_pipe_u1(u, 80.0, 10.0,
#     ρf, Cf, Ai, 100.0, Rpg, Rgg, Rgr)

# T_supply_grout, T_return_grout = temperature_grout_u1(T_supply, T_return, 10.0, Rpg, Rgg, Rgr)



##
using GLMakie
fig = Figure(size = (800, 600))
ax = Axis(fig[1, 1]; yreversed=true)
z = collect(range(0, 100.0, step=0.1))
lines!(ax, convert_from_si.(sol[:pipe_supply].(z), :Celsius), z)
lines!(ax, convert_from_si.(sol[:pipe_return].(z), :Celsius), z)
lines!(ax, convert_from_si.(sol[:grout_supply].(z), :Celsius), z)
lines!(ax, convert_from_si.(sol[:grout_return].(z), :Celsius), z)
fig