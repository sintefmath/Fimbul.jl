using Jutul, JutulDarcy
using Fimbul
using HYPRE
using GLMakie

## Set up case

# initial_condition = x -> convert_to_si(90.0, :Celsius) .+ 0.0*x
# case, sol = analytical_1d(num_cells = 4000, num_steps = 500, initial_condition = initial_condition);
case, sol = analytical_1d(num_cells = 4000, num_steps = 500)

## Simulate
res = simulate_reservoir(case, info_level = 0)
plot_reservoir(case, res.states)

##
fig = Figure(size = (800, 400), fontsize = 20)
ax1 = Axis(fig[1, 1];) 
mesh = physical_representation(reservoir_model(case.model).data_domain)
geo = tpfv_geometry(mesh)
x = geo.cell_centroids[1,:]
t = cumsum(case.dt)   
n = 5
timesteps = Int.(round.(collect(range(1, stop = length(t), length = n))))
for n = timesteps
    T_sim = res.states[n][:Temperature]
    T_analytical = sol(x, t[n])
    lines!(ax1, x, T_analytical, linestyle = (:dash, 1), linewidth = 6, color = :black)
    lines!(ax1, x, T_sim, color = :blue)
    ΔT = norm((T_sim .- T_analytical)./T_analytical, Inf)
    println("Time: $(t[n]), Relative error: $(ΔT)")
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
