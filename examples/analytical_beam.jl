using Jutul, JutulDarcy
using Fimbul
using HYPRE

##
case, sol = beam_thermal(num_cells = 1000000, num_steps = 200);

##

res = simulate_reservoir(case, info_level = 2)
plot_reservoir(case, res.states)

##

mesh = physical_representation(reservoir_model(case.model).data_domain)
geo = tpfv_geometry(mesh)
x = geo.cell_centroids[1,:]
t = cumsum(case.dt)
for (k, tk) = enumerate(t)
    T_sim = res.states[k][:Temperature]
    T_analytical = sol(x, tk)
    ΔT = (T_sim .- T_analytical).*geo.volumes./sum(geo.volumes)
    println(norm(ΔT))
end
