using Jutul, JutulDarcy
using Fimbul
using HYPRE
using GLMakie

## Set up case

T0 = convert_to_si(10.0, :Celsius)
initial_condition = x -> convert_to_si(90.0, :Celsius) .+ 0.0*x - T0

initial_condition = x -> convert_to_si(90.0, :Celsius).*(x < 50) + convert_to_si(10.0, :Celsius).*(x >= 50) - T0

initial_condition = x -> convert_to_si(100.0, :Celsius).*(x < 25) + 
    convert_to_si(10.0, :Celsius).*(25 <= x < 50) +
    convert_to_si(50.0, :Celsius).*(50 <= x < 75) +
    convert_to_si(75.0, :Celsius).*(75 <= x) -
    T0


case, sol = analytical_1d(num_cells = 1000, num_steps = 500, initial_condition = initial_condition);
# case, sol = analytical_1d(num_cells = 8000, num_steps = 1500)

## Simulate
res = simulate_reservoir(case, info_level = 0)
plot_reservoir(case, res.states)

##
fig = Figure(size = (800, 400), fontsize = 20)
ax1 = Axis(fig[1, 1];) 
mesh = physical_representation(reservoir_model(case.model).data_domain)
geo = tpfv_geometry(mesh)
x = geo.cell_centroids[1,:]
xa = range(0, 100, length = 500)
t = cumsum(case.dt)   
n = 5
timesteps = Int.(round.(collect(range(1, stop = length(t), length = n))))
timesteps = 1:4:32
for n = timesteps
    T_sim = res.states[n][:Temperature]
    T_analytical = sol(xa, t[n])
    lines!(ax1, xa, T_analytical, linestyle = (:dash, 1), linewidth = 6, color = :black)
    lines!(ax1, x, T_sim, color = :black, linewidth = 2)
    # ΔT = norm((T_sim .- T_analytical)./T_analytical, Inf)
    # println("Time: $(t[n]), Relative error: $(ΔT)")
end
fig

## Compare to analytical solution
# t = cumsum(case.dt)
# ΔT = []
# for (k, tk) = enumerate(t)
#     T_sim = res.states[k][:Temperature]
#     T_analytical = sol(x, tk)
#     ΔTk = (T_sim .- T_analytical)./T_analytical
#     push!(ΔT, norm(ΔTk, Inf))
# end
