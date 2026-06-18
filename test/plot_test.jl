using GLMakie
using Fimbul

field = Fimbul.circular_pattern(48, 5.0; num_sectors = 6, depths = [0.0, 0.5, 100, 125])
xy = hcat([w[1:2, 1] for w in vcat(field...)]...)

fig = Figure()
ax = Axis(fig[1, 1], aspect = DataAspect())
scatter!(ax, xy[1, :], xy[2, :])
fig
