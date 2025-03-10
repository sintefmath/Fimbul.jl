## # Heat equation in 1D
# This example goes through the the classical solution of the heat equation in
# 1D, and compares the analytical solution to the numerical solution obtained
# using JutulDary to verify the numerical scheme.
using Jutul, JutulDarcy
using Fimbul
using HYPRE
using GLMakie

to_kelvin = T -> convert_to_si(T, :Celsius)
to_celsius = T -> convert_from_si(T, :Celsius)

## ## The 1D heat equation
# We consider a piece of solid rock with length $L$ and constant thermal
# donductivity $\lambda$, and heat capacity $C_p$ and density $\rho$, for
# which conservation of energy takes the form of the archetypical heat equation:
# 
# ``\dot{T} = \alpha \Delta T``, \alpha = \frac{\lambda}{\rho C_p}``
#
# where $T$ is the temperature and $\alpha$ is called the thermal diffusivity.
# Imposing equal and fixed boundary conditions $T(0, t) = T(L, t) = T_b$ and
# initial conditions T(0, x) = T_0(x), the solution to this equation can be
# found in any textbook on partial differential equations, and is given by
#
# ``T(x, t) = T_b + \sum_{k = 1}^{\infty} C_k \exp\big(-\alpha \frac{(k\pi}{L}\big^2 t\right) \sin\big(\frac{k\pi}{L} x)``
#
# where the coefficients $C_k$ are determined by the initial condition and
# boundary conditions:xs
#
# ``C_k = \frac{2}{L} \int_0^L (T_0(x) - T_b) \sin\big(\frac{n\pi}{L} x) dx``
#

# ## Simple initial conditions
# We first consider a simple initial condition where the initial temperature
# profile takes the form of a sine curve,
#
# ``T_0(x) = T_\max\sin\big(\frac{2\pi}{L} x\big)``

L = 100.0
T_b = to_kelvin(10.0)
T_0 = x -> to_kelvin(90.0).*sin(π*x/L) .+ T_b;
setup_case = (nx, nt) -> analytical_1d(
    length_x = L, temperature_boundary = T_b, initial_condition = T_0,
    num_cells = nx, num_steps = nt);

## ### Simulate with 100 cells and 100 timesteps
out = setup_case(100, 100)
results = simulate_reservoir(out[1], info_level = 0)

## Plot the temperature profile
# We set up a simple function for plotting the numerical and analytical
# temperature profiles at a selected number of timesteps from the initial time
# up to a fraction of the total time.
function plot_temperature_1d(case, sol_n, sol_a, x_n, t_n, n, frac)
    fig = Figure(size = (600, 600), fontsize = 20)
    ax = Axis(fig[1, 1]; xlabel = "Distance (m)", ylabel = "Temperature (°C)")

    x_a = range(0, 100, length = 500)
    timesteps = Int.(round.(range(0,frac*length(t_n),n)))

    colors = cgrad(:Spectral, n, categorical = true)
    for (i,k) = enumerate(timesteps)
        if k == 0
            T_n = to_celsius.(case.state0[:Reservoir][:Temperature])
            lines!(ax, x_n, T_n, linestyle = (:dash, 1), linewidth = 6, color = colors[i], label = "Analytical")
            lines!(ax, x_n, T_n, linewidth = 2, color = colors[i], label = "Numerical")
        else
            T_a = to_celsius.(sol_a(x_a, t_n[k]))
            T_n = to_celsius.(sol_n.states[k][:Temperature])
            lines!(ax, x_a, T_a, linestyle = (:dash, 1), linewidth = 6, color = colors[i])
            lines!(ax, x_n, T_n, linewidth = 2, color = colors[i])
        end
    end
    Legend(fig[1,2], ax, loc = :northeast)
    fig
end

plot_temperature_1d(out[1], results, out[2], 10, 0.5)

## ## Piecewise constant initial conditions
# We now consider a more complex initial condition where the initial temperature
# profile is piecewise constant, with four different constant values.
initial_condition = x ->
    to_kelvin(100.0).*(x < 25) + 
    to_kelvin(20.0).*(25 <= x < 50) +
    to_kelvin(50.0).*(50 <= x < 75) +
    to_kelvin(75.0).*(75 <= x);

# We can still compute the analytical solution using the formula above, but the
# sum will now be an infinite series. The function `analytical_1d` handles this
# pragmatically by cutting off the series when the contribution of the next term
# is less than a 1e-6.
out = analytical_1d(
    length_x = L, temperature_boundary = T_b, initial_condition = initial_condition,
    num_cells = 500, num_steps = 500);

results = simulate_reservoir(out[1], info_level = 0)
plot_temperature_1d(out[1], results, out[2], 10, 0.175)

## ### Convergence study
# Next, we perform a convergence study by simulating the same problem with
# increasing number of cells and timesteps, and comparing the numerical and
# analytical solutions

Δx, err = [], []
for nx in 2 .^range(3, 12)
    for nt = 100
        out = setup_case(nx, nt)
        case, sol, x, t = out
        
        sim, cfg = setup_reservoir_simulator(case;
            relaxation=true,
            tol_cnv=1e-8,
            info_level=-1,
        );
        cfg[:tolerances][:Reservoir][:default] = 1e-8
       
        results = simulate_reservoir(out[1], info_level = 0, simulator = sim, config = cfg)
        nt_comp = nt÷2
        nt_comp = 5
        x = tpfv_geometry(physical_representation(reservoir_model(
            case.model).data_domain)).cell_centroids[1,:]
        dx = x[2] - x[1]
        T_n = results.states[nt_comp][:Temperature]
        T_a = sol(x, t[nt_comp])

        push!(Δx, dx)
        push!(err, norm(dx.*(T_n .- T_a)./L))
    end
end

fig = Figure(size = (800, 800), fontsize = 20)
ax = Axis(fig[1, 1]; xlabel = "Δx", ylabel = "Error", xscale = log2, yscale = log2)
lines!(ax, Δx, err, linewidth = 2, color = :black)
scatter!(ax, Δx, err, marker = :rect, markersize = 20, color = :black)
fig