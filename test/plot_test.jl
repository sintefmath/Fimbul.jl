using GLMakie
using Fimbul

field = Fimbul.polygonal_pattern(200, 0.5,6; num_sectors = 6, depths = [0.0, 0.5, 100, 125])

fig = Figure()
ax = Axis(fig[1, 1], aspect = DataAspect(),limits = (-20, 20, -20, 20))
colors = Makie.wong_colors()
for (sno, sector) in enumerate(field)
    xy = hcat([w[1:2, 1] for w in sector]...)
    scatter!(ax, xy[1, :], xy[2, :]; color = colors[mod1(sno, length(colors))], label = "Sector $sno")
end
axislegend(ax)
fig
