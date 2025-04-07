
xw = fibonacci_pattern_2d(100; radius = missing, spacing = 5.0)
cell_constraints = map(x -> [x], xw)

depths = [0.0, 10.0, 20, 30]
mesh, layers, metrics = horizontal_fractured_mesh(cell_constraints, depths, 10; aperture = 1e-3)

##
using GLMakie
GLMakie.activate!()

fig = Figure(size = (1200, 800))
ax = Axis3(fig[1, 1], zreversed = true)
Jutul.plot_mesh_edges!(ax, mesh, alpha = 0.5)
fig

