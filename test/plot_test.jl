using GLMakie
using Fimbul

field = Fimbul.polygonal_pattern(200, 0.5,6; num_sectors = 6, depths = [0.0, 0.5, 100, 125])
xy = hcat([w[1:2, 1] for w in vcat(field...)]...)

fig = Figure()
ax = Axis(fig[1, 1], aspect = DataAspect(),limits = (-20, 20, -20, 20))
scatter!(ax, xy[1, :], xy[2, :])
fig
