function Fimbul.plot_btes_field!(ax, case::JutulCase; kwargs...)
    data = case.input_data
    haskey(data, :field) || error("The case has no BTES well field. \
        Set it up with btes().")
    return Fimbul.plot_btes_field!(
        ax, data[:field];
        pattern = get(data, :pattern, missing),
        well_spacing = get(data, :well_spacing, missing),
        kwargs...
    )
end

function Fimbul.plot_btes_field!(
        ax, field::AbstractVector;
        pattern = missing,
        well_spacing = missing,
        structure = true,
        structure_color = (:gray, 0.4),
        sector_colors = Makie.to_colormap(:tab10),
        linewidth = 2,
        markersize = 10,
        first_markersize = 20,
    )

    xy_all = hcat([wc[1:2, 1] for sector in field for wc in sector]...)
    if structure && !ismissing(pattern)
        R = maximum(sqrt.(sum(xy_all .^ 2, dims = 1)))
        if pattern == :sunflower
            plot_sunflower_spines!(ax, size(xy_all, 2), R, structure_color)
        elseif pattern == :circular
            if ismissing(well_spacing)
                ## The innermost ring is one spacing from the center
                r = sqrt.(sum(xy_all .^ 2, dims = 1))
                well_spacing = minimum(r[r .> 1.0e-8 * R]; init = R)
            end
            plot_circular_rings!(ax, R, well_spacing, structure_color)
        elseif pattern == :rectangular
            plot_cartesian_grid!(ax, xy_all, structure_color)
        else
            error("Unknown pattern: $pattern. Supported patterns are \
                :sunflower, :circular and :rectangular.")
        end
    end
    for (sno, sector) in enumerate(field)
        color = sector_colors[mod1(sno, length(sector_colors))]
        xy = hcat([wc[1:2, 1] for wc in sector]...)
        lines!(ax, xy[1, :], xy[2, :]; color = color, linewidth = linewidth)
        scatter!(
            ax, xy[1, :], xy[2, :]; color = color, markersize = markersize,
            label = "Sector $sno"
        )
        scatter!(
            ax, xy[1, 1:1], xy[2, 1:1]; color = color,
            markersize = first_markersize
        )
    end
    return ax
end

function Fimbul.plot_btes_field(x; title = "", legend = true, kwargs...)
    fig = Figure()
    ax = Axis(
        fig[1, 1]; title = title, aspect = DataAspect(),
        xlabel = "x [m]", ylabel = "y [m]"
    )
    Fimbul.plot_btes_field!(ax, x; kwargs...)
    legend && Legend(fig[1, 2], ax)
    return fig
end

# Spines of a sunflower pattern of N wells with outer radius R (see
# fibonacci_pattern_2d): wells F apart in sunflower order, for a Fibonacci
# number F, lie on a spiral. The two largest families with at least two wells
# per spine are drawn, which include the family that spiral sectors are made of.
function plot_sunflower_spines!(ax, N, R, color; num_points = 200)
    N < 2 && return
    Δθ = 2π / Base.MathConstants.golden^2
    fib = [1, 2]
    while fib[end] + fib[end - 1] <= N / 2
        push!(fib, fib[end] + fib[end - 1])
    end
    for F in fib[(end - 1):end]
        ## Turning angle between wells F apart, taken as the smallest equivalent
        ## angle, so that the spine through them is a smooth spiral
        δ = mod(F * Δθ + π, 2π) - π
        for n in 1:min(F, N)
            t = range(0.0, (N - n) / F, length = num_points)
            r = R .* sqrt.((n .+ F .* t .- 0.5) ./ (N - 0.5))
            θ = n * Δθ .+ δ .* t
            lines!(ax, r .* cos.(θ), r .* sin.(θ); color = color, linewidth = 1)
        end
    end
    return
end

# Rings of a circular pattern, one well spacing apart
function plot_circular_rings!(ax, R, spacing, color)
    for k in 1:round(Int, R / spacing)
        arc!(ax, Point2f(0, 0), k * spacing, 0, 2π; color = color, linewidth = 1)
    end
    return
end

# Grid lines through the wells of a rectangular pattern
function plot_cartesian_grid!(ax, xy, color)
    xs, ys = sort(unique(xy[1, :])), sort(unique(xy[2, :]))
    for x in xs
        lines!(ax, [x, x], [ys[1], ys[end]]; color = color, linewidth = 1)
    end
    for y in ys
        lines!(ax, [xs[1], xs[end]], [y, y]; color = color, linewidth = 1)
    end
    return
end
