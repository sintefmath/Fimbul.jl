using Jutul, JutulDarcy
using Fimbul
using HYPRE

## Set up case
case, sol = beam_thermal(num_cells = 1000, num_steps = 200);

## Simulate
res = simulate_reservoir(case, info_level = 0)
plot_reservoir(case, res.states)

## Compare to analytical solution
mesh = physical_representation(reservoir_model(case.model).data_domain)
geo = tpfv_geometry(mesh)
x = geo.cell_centroids[1,:]
t = cumsum(case.dt)
ΔT = []
for (k, tk) = enumerate(t)
    T_sim = res.states[k][:Temperature]
    T_analytical = sol(x, tk)
    ΔTk = (T_sim .- T_analytical)./T_analytical
    push!(ΔT, norm(ΔTk, Inf))
end
