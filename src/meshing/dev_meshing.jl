using Fimbul

##
<<<<<<< Updated upstream
xw = fibonacci_pattern_2d(100; radius = missing, spacing = 5.0)

##
mesh, m = extruded_mesh(xw, [0, 50.0, 75.0], hz = [0.5, 5.0])

##
using GLMakie
GLMakie.activate!()

fig = Figure(size = (1200, 800))
ax = Axis3(fig[1, 1], zreversed = true)
Jutul.plot_mesh_edges!(ax, mesh, alpha = 0.5)
fig

##
x, xc = Fimbul.get_convex_hull(xw)

v = map(x -> x .- xc, x)
v = map(v -> v./norm(v,2), v)

<<<<<<< Updated upstream
xb = Fimbul.offset_boundary(x, xc, 50)