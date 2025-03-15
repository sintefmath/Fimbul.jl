using Fimbul

xw = fibonachi_pattern_2d(100; radius = missing, spacing = 5.0)

##
using GLMakie
fig = Figure()
ax = Axis(fig[1, 1], title = "BTES well coordinates", aspect = 1.0)
scatter!(ax, xw, markersize=10)
fig

##

for (i,x) in enumerate(xw)
    d = Inf
    for (j,y) in enumerate(xw)
        if i == j
            continue
        end
        d = min(d, norm(x .- y,2))
    end
    println("Minimum distance from well $i to other wells: $d")
end

##
mesh = extruded_mesh(xw, 100.0)