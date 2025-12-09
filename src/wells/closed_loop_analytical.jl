function temperature_closed_loop_pipe(u, T_in, T_rock, ρ, Cp, A, L, R1Δ, R2Δ, R12Δ)

    f1, f2, f3, f4, f5 = f_functions(u, ρ, Cp, A, R1Δ, R2Δ, R12Δ)

    domain = z -> (0.0, z)
    
    integrand_supp = z -> (ξ,p) -> T_rock*f4(z - ξ)
    I_supp = z -> solve(IntegralProblem(integrand_supp(z), domain(z)), QuadGKJL())
    
    integrand_ret = z -> (ξ,p) -> T_rock*f5(z - ξ)
    I_ret = z -> solve(IntegralProblem(integrand_ret(z), domain(z)), QuadGKJL())

    # Find T_out from T_supply(L) = T_return(L)
    T_out = (T_in*(f1(L) .+ f2(L)) .+ I_supp(L) .+ I_ret(L))/(f3(L) .- f2(L))
    T_supply = z -> T_in*f1(z) .+ T_out*f2(z) .+ I_supp(z)
    T_return = z -> -T_in*f2(z) .+ T_out*f3(z) .- I_ret(z)

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

function analytical_closed_loop_u1(rate, temperature_in, temperature_rock,
    density_fluid, heat_capacity_fluid, well::DataDomain)

    # Identify pipe and grout cells
    section = well[:section]
    is_pipe_cell = last.(section) .== :pipe_left
    is_grout_cell = last.(section) .== :grout_left
    # Get length
    length = sum(well[:cell_length][is_pipe_cell])
    # Get pipe outer radius
    rp_i = well[:radius][is_pipe_cell][1]
    wall_thickness_pipe = well[:casing_thickness][is_pipe_cell][1]
    radius_pipe = rp_i + wall_thickness_pipe
    # Get grout radius
    L = well[:cell_length][is_grout_cell][1]
    vg = well[:volume_override_grouting][is_grout_cell][1]
    vp = radius_pipe^2*π*L
    radius_grout = sqrt((2*vg + 2*vp)/(π*L))
    # Get other parameters
    pipe_spacing = maximum(well[:pipe_spacing])
    thermal_conductivity_grout = maximum(well[:grouting_thermal_conductivity])
    thermal_conductivity_pipe = maximum(well[:casing_thermal_conductivity])

    return analytical_closed_loop_u1(
        rate, temperature_in, temperature_rock,
        density_fluid, heat_capacity_fluid,
        length, radius_grout, radius_pipe, wall_thickness_pipe, pipe_spacing,
        thermal_conductivity_grout, thermal_conductivity_pipe)

end

function analytical_closed_loop_u1(rate, temperature_in, temperature_rock,
    density_fluid, heat_capacity_fluid,
    length, radius_grout, radius_pipe, wall_thickness_pipe, pipe_spacing,
    thermal_conductivity_grout, thermal_conductivity_pipe)

    Rpg, Rgr, Rgg = Fimbul.closed_loop_thermal_resistance_u1(
        radius_grout, radius_pipe, wall_thickness_pipe, pipe_spacing,
        thermal_conductivity_grout, thermal_conductivity_pipe)
    A = π*(radius_pipe - wall_thickness_pipe)^2

    u = rate/A
    Tps, Tpr = temperature_pipe_u1(
        u, temperature_in, temperature_rock, density_fluid, heat_capacity_fluid, A, length, Rpg, Rgr, Rgg)
    
    Tgs, Tgr = temperature_grout_u1(
        Tps, Tpr, temperature_rock, Rpg, Rgr, Rgg)
    
    sol = Dict()
    sol[:pipe_left] = Tps
    sol[:pipe_right] = Tpr
    sol[:grout_left] = Tgs
    sol[:grout_right] = Tgr

    return sol

end

function temperature_pipe_u1(u, T_in, T_rock, ρ, Cp, A, L, Rpg, Rgr, Rgg)

    u1 = 1/Rpg + 1/Rgr + 1/Rgg
    R1Δ = Rpg + Rgr
    R2Δ = Rpg + Rgr
    R12Δ = ((u1*Rpg*Rgg)^2 - Rpg^2)/Rgg

    return temperature_closed_loop_pipe(u, T_in, T_rock, ρ, Cp, A, L, R1Δ, R2Δ, R12Δ)

end

function temperature_grout_u1(T_supply, T_return, T_rock, Rpg, Rgr, Rgg)

    u1 = 1/Rpg + 1/Rgr + 1/Rgg
    Tgs = z -> (
        T_rock/Rgr + 
        T_return(z)/Rpg + 
        (T_rock/Rgr + T_supply(z)/Rpg)*u1*Rgg
        )*Rgg/
        ((Rgg*u1)^2-1)

    Tgr = z -> (
        Tgs(z)/Rgg + 
        T_return(z)/Rpg + 
        T_rock/Rgr
        )*1/u1

    return Tgs, Tgr

end

function analytical_closed_loop_coaxial(rate, temperature_in, temperature_rock,
    density_fluid, heat_capacity_fluid, well::DataDomain; kwargs...)

    # Identify pipe and grout cells
    section = well[:section]
    is_inner_pipe_cell = last.(section) .== :pipe_inner
    is_outer_pipe_cell = last.(section) .== :pipe_outer
    is_grout_cell = last.(section) .== :grout
    # Get length
    length = sum(well[:cell_length][is_inner_pipe_cell])
    # Get pipe radii and wall thicknesses
    radius_inner_pipe = well[:radius_inner][is_outer_pipe_cell][1]
    wall_thickness_inner_pipe = well[:casing_thickness][is_inner_pipe_cell][1]
    rpo_i = well[:radius][is_outer_pipe_cell][1]
    wall_thickness_outer_pipe = well[:casing_thickness][is_outer_pipe_cell][1]
    radius_outer_pipe = rpo_i + wall_thickness_outer_pipe
    # Get grout radius
    L = well[:cell_length][is_grout_cell][1]
    vg = well[:volume_override_grouting][is_grout_cell][1]
    radius_grout = sqrt((vg + radius_outer_pipe^2*pi*L)/(π*L))
    # Get other parameters
    thermal_conductivity_grout = maximum(well[:grouting_thermal_conductivity])
    thermal_conductivity_inner_pipe = maximum(well[:casing_thermal_conductivity][is_inner_pipe_cell])
    thermal_conductivity_outer_pipe = maximum(well[:casing_thermal_conductivity][is_outer_pipe_cell])

    println("""
    radius_inner_pipe: $radius_inner_pipe
    wall_thickness_inner_pipe: $wall_thickness_inner_pipe
    radius_outer_pipe: $radius_outer_pipe
    wall_thickness_outer_pipe: $wall_thickness_outer_pipe
    radius_grout: $radius_grout
    length: $length
    thermal_conductivity_grout: $thermal_conductivity_grout
    thermal_conductivity_inner_pipe: $thermal_conductivity_inner_pipe
    thermal_conductivity_outer_pipe: $thermal_conductivity_outer_pipe
    """)

    # return nothing

    return analytical_closed_loop_coaxial(rate, temperature_in, temperature_rock,
    density_fluid, heat_capacity_fluid,
    length, radius_grout,
    radius_inner_pipe, wall_thickness_inner_pipe,
    radius_outer_pipe, wall_thickness_outer_pipe,
    thermal_conductivity_grout,
    thermal_conductivity_inner_pipe,
    thermal_conductivity_outer_pipe;
    kwargs...)

end

function analytical_closed_loop_coaxial(rate, temperature_in, temperature_rock,
    density_fluid, heat_capacity_fluid,
    length, radius_grout,
    radius_inner_pipe, wall_thickness_inner_pipe,
    radius_outer_pipe, wall_thickness_outer_pipe,
    thermal_conductivity_grout,
    thermal_conductivity_inner_pipe,
    thermal_conductivity_outer_pipe;
    inlet = :outer)

    Rpp, Rpg, Rgr = closed_loop_thermal_resistance_coaxial(
        radius_grout,
        radius_inner_pipe, wall_thickness_inner_pipe,
        radius_outer_pipe, wall_thickness_outer_pipe,
        thermal_conductivity_grout,
        thermal_conductivity_inner_pipe,
        thermal_conductivity_outer_pipe
    )
    if inlet == :outer
        A = π*((radius_outer_pipe - wall_thickness_outer_pipe)^2 - radius_inner_pipe^2)
    elseif inlet == :inner
        A = π*(radius_inner_pipe - wall_thickness_inner_pipe)^2
    else
        error("Inlet must be :outer or :inner")
    end
    u = rate/A
    Tps, Tpr = temperature_pipe_coaxial(
        u, temperature_in, temperature_rock,
        density_fluid, heat_capacity_fluid,
        A, length, Rpg, Rgr, Rpp, inlet
    )
    
    To = inlet == :inner ? Tpr : Tps
    Ti = inlet == :inner ? Tps : Tpr
    Tg = temperature_grout_coaxial(
        To, temperature_rock, Rpg, Rgr)
    
    sol = Dict()
    sol[:pipe_inner] = Ti
    sol[:pipe_outer] = To
    sol[:grout] = Tg

    return sol

end

function temperature_pipe_coaxial(u, T_in, T_rock, ρ, Cp, A, L, Rpg, Rgr, Rpp, inlet)

    if inlet == :outer
        R1Δ = Rpg + Rgr
        R2Δ = Inf
    elseif inlet == :inner
        R1Δ = Inf
        R2Δ = Rpg + Rgr
    else
        error("Inlet must be :outer or :inner")
    end
    R12Δ = Rpp

    return temperature_closed_loop_pipe(u, T_in, T_rock, ρ, Cp, A, L, R1Δ, R2Δ, R12Δ)

end

function temperature_grout_coaxial(T_annulus, T_rock, Rpg, Rgr)

    T_grout = z -> Rpg/(Rpg + Rgr)*(T_rock - T_annulus(z)) + T_annulus(z)
    return T_grout

end