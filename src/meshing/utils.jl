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

function curve_measure(x)

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

function add_geometry_0d(points, h, tag = 0, embed = false)

    tags = Vector{Int}(undef, length(points))
    for (i, pt) = enumerate(points)
        x, y = pt
        tags[i] = gmsh.model.geo.addPoint(x, y, 0.0, h, tag + i)
    end
    if embed
        gmsh.model.geo.synchronize()
        gmsh.model.mesh.embed(0, tags, 2, 1)
    end

    return tags

end

function add_geometry_1d(points, h, tag_0d = 0, tag_1d = 0, embed = false)
    
    np = length(points)
    tags = Vector{Int}(undef, np)
    tags_0d = add_geometry_0d(points, h, tag_0d, embed)
    for i = 1:np
        tags[i] = gmsh.model.geo.addLine(tags_0d[i], tags_0d[mod(i, np) + 1], tag_1d + i)
    end
    if embed
        gmsh.model.geo.synchronize()
        gmsh.model.mesh.embed(1, tags, 2, 1)
    end

    return tags, tags_0d

end

function add_geometry_2d(points, h, tag_0d = 0, tag_1d = 0, tag_2d = 0, embed = false)

    np = length(points)
    tags = Vector{Int}(undef, np)
    tags_1d, tags_0d = add_geometry_1d(points, h, tag_0d, tag_1d, embed)
    tag = gmsh.model.geo.addCurveLoop(tags_1d, tag_1d + 1)
    gmsh.model.geo.addPlaneSurface([1], tag_2d + 1)

    return tag, tags_1d, tags_0d

end

# function make_domain_2d(x, h, tag_0d = 0, tag_1d = 0, tag_2d = 0)

#     n = length(x)
#     for i = 1:n
#         xi, yi = x[i]
#         gmsh.model.geo.addPoint(xi, yi, 0.0, h, tag_0d + i)
#     end
#     for i = 1:n
#         gmsh.model.geo.addLine(i, mod(i, n) + 1, tag_1d + i)
#     end
#     gmsh.model.geo.addCurveLoop(collect(1:n), tag_1d + 1)
#     gmsh.model.geo.addPlaneSurface([1], tag_2d + 1)

# end