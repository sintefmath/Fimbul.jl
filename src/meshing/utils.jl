function fibonacci_pattern_2d(num_points; spacing = 5.0, radius = missing)
    
    # Define coordinate funciton using golden ratio
    N = num_points
    φ = (1 + sqrt(5))/2 # ≈ golden ratio
    Δθ = 2*π/φ^2
    p = (n, r) -> r.*sqrt(n/N).*(cos(n*Δθ), sin(n*Δθ))

    # Set radius from spacing if not provided
    if !ismissing(spacing)
        @assert ismissing(radius) "Please provide either radius or spacing"
        radius = 1.3513*spacing/norm(p(1,1.0) .- p(2,1.0),2)
    elseif !ismissing(radius)
        @assert ismissing(spacing) "Please provide either radius or spacing"
    end

    x = [p(n, radius) for n = 1:N];

    return x

end

using Meshes
function get_convex_hull(x)

    ptset = PointSet(x)
    ch = convexhull(ptset)
    xc = Tuple(centroid(ch).coords)
    x = map(p -> Tuple(p.coords), vertices(ch))
    return x, xc

end

function offset_boundary(x, xc, offset)

    @assert offset > 0 "offset must be larger than 0"

    v = map(x -> x .- xc, x)
    # v = map(v -> v./norm(v,2), v)
    xs = map((x,v) -> x .+ offset.*v, x, v)
    return xs

end

function chain_measure(x)

    x = vcat(x, x[1])
    dx = [x[i+1] .- x[i] for i = 1:length(x)-1]
    return sum(norm.(dx,2))

end

function min_max_distance(x)

    d_min, d_max = Inf, -Inf
    for (i,x1) in enumerate(x)
        for (j,x2) in enumerate(x)
            i == j ? continue : nothing
            d_ij = norm(x1 .- x2,2)
            d_min = min(d_min, d_ij)
            d_max = max(d_max, d_ij)
        end
    end
    return d_min, d_max

end

function min_distance(x)

    return min_max_distance(x)[1]

end

function max_distance(x)

    return min_max_distance(x)[2]
    
end