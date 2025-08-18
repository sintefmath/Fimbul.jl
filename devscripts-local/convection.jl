using Jutul, JutulDarcy, Fimbul
using HYPRE
using GLMakie

##
case, plot_args = Fimbul.fimbul_logo()

##
sim, cfg = setup_reservoir_simulator(case; max_timestep = Inf, info_level = 2)
res = simulate_reservoir(case; simulator = sim, config = cfg)

##
GLMakie.closeall()
Tf = res.states[end][:Temperature]
xy = geo.cell_centroids[1:2,:]

colors = cgrad(:seaborn_icefire_gradient, 8, categorical=true)
# colors = cgrad(:hot, 8, categorical=true)
msh = physical_representation(reservoir_model(case.model).data_domain)
plot_cell_data(msh, Tf; colormap = colors, shading = NoShading)

##
plot_reservoir(case, res.states;
colormap = :seaborn_icefire_gradient, shading = NoShading,
aspect = :data)

##

GLMakie.closeall()
fig = Figure()
steps = Int.(ceil.(range(1,length(case.dt), 6)))

# n_plot = 6
α = (length(case.dt)).^(1/(n_plot-1))
steps = Int.(ceil.(α.^(0:(n_plot-1))))
for (sno, step) in enumerate(steps)
    row, col = (sno-1)÷3+1, (sno-1)%3+1
    ax = Axis3(fig[row, col]; zreversed = true, plot_args...)
    T = res.states[step][:Temperature]
    plot_cell_data!(ax, msh, T; shading = NoShading, colormap = :seaborn_icefire_gradient)
    hidedecorations!(ax)
end
fig