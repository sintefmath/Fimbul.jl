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

function get_convex_hull(points, tol = eps(Float64))
    # Ensure there are at least 3 points
    @assert length(points) ≥ 3 "At least 3 points are required to compute a convex hull."

    # Helper function to compute the cross product of vectors (p1 -> p2) and (p1 -> p3)
    function cross_product(p1, p2, p3)
        (p2[1] - p1[1]) * (p3[2] - p1[2]) - (p2[2] - p1[2]) * (p3[1] - p1[1])
    end

    # Start with the leftmost point
    hull = Vector{Tuple{Float64, Float64}}()
    start = argmin(map(p -> p[1], points))  # Index of the leftmost point
    current = start

    while true
        push!(hull, points[current])  # Add the current point to the hull
        next_point = mod1(current + 1, length(points))  # Next candidate point

        for i in 1:length(points)
            if cross_product(points[current], points[next_point], points[i]) < 0
                next_point = i  # Update the next point if it is more counterclockwise
            end
        end

        current = next_point  # Move to the next point
        if current == start
            break  # Stop when we return to the starting point
        end
    end

    hull = first(ramer_douglas_peucker(hull; tolerance = tol))

    return hull
end

function point_to_line_distance_2d(point, line_start, line_end)
    """Calculate perpendicular distance from point to line segment in 2D"""
    x0, y0 = point
    x1, y1 = line_start
    x2, y2 = line_end
    
    # Line vector
    dx = x2 - x1
    dy = y2 - y1
    
    # If line segment has zero length, return distance to start point
    if dx == 0 && dy == 0
        return sqrt((x0 - x1)^2 + (y0 - y1)^2)
    end
    
    # Calculate perpendicular distance using cross product formula
    numerator = abs(dy * x0 - dx * y0 + x2 * y1 - y2 * x1)
    denominator = sqrt(dx^2 + dy^2)
    
    return numerator / denominator
end

function point_to_line_distance_3d(point, line_start, line_end)
    """Calculate perpendicular distance from point to line segment in 3D"""
    # Vector from line start to point
    v = point .- line_start
    # Line direction vector
    u = line_end .- line_start
    
    # If line segment has zero length, return distance to start point
    if norm(u) ≈ 0
        return norm(v)
    end
    
    # Cross product to get perpendicular component
    cross_prod = cross(v, u)
    
    # Distance is magnitude of cross product divided by line length
    return norm(cross_prod) / norm(u)
end

function ramer_douglas_peucker(points::Vector{T}, tolerance::Float64) where T <: Union{Tuple{Float64, Float64}, Tuple{Float64, Float64, Float64}}
    """
    Simplify a curve using the Ramer-Douglas-Peucker algorithm.
    
    Parameters:
    - points: Vector of tuples representing 2D (x,y) or 3D (x,y,z) points
    - tolerance: Maximum allowed deviation from original curve
    
    Returns:
    - simplified_points: Vector of simplified points
    - indices: Vector of indices into original points array
    """
    
    n = length(points)
    
    # Need at least 2 points
    if n <= 2
        return points, collect(1:n)
    end
    
    # Recursive helper function
    function rdp_recursive(start_idx::Int, end_idx::Int)
        # Find the point with maximum distance from line segment
        max_dist = 0.0
        max_idx = 0
        
        start_point = points[start_idx]
        end_point = points[end_idx]
        
        for i in (start_idx + 1):(end_idx - 1)
            point = points[i]
            
            # Calculate distance based on dimensionality
            if length(start_point) == 2
                dist = point_to_line_distance_2d(point, start_point, end_point)
            else  # 3D
                dist = point_to_line_distance_3d(point, start_point, end_point)
            end
            
            if dist > max_dist
                max_dist = dist
                max_idx = i
            end
        end
        
        # If max distance is greater than tolerance, recursively simplify
        if max_dist > tolerance && max_idx > 0
            # Recursively simplify the two segments
            left_indices = rdp_recursive(start_idx, max_idx)
            right_indices = rdp_recursive(max_idx, end_idx)
            
            # Combine results (remove duplicate point at junction)
            return vcat(left_indices, right_indices[2:end])
        else
            # All points are within tolerance, return only endpoints
            return [start_idx, end_idx]
        end
    end
    
    # Get indices of simplified points
    indices = rdp_recursive(1, n)
    
    # Extract simplified points
    simplified_points = points[indices]
    
    return simplified_points, indices
end

# Convenience function with default tolerance
function ramer_douglas_peucker(points::Vector{T}; tolerance::Float64 = 1.0) where T <: Union{Tuple{Float64, Float64}, Tuple{Float64, Float64, Float64}}
    return ramer_douglas_peucker(points, tolerance)
end

function offset_boundary(x, offset; h = missing, n = 6)

    @assert offset > 0 "Offset must be larger than 0"

    # Number of points to sample around each point
    if !ismissing(h)
        n = max(Int(ceil(2π*offset/h)) + 1, n)
    else
        h = 2π*offset/n
    end

    function sample_circle(x0, r)

        θ = LinRange(0, 2π, n+1)
        xc = [x0 .+ (r*cos(θ[i]), r*sin(θ[i])) for i = 1:n]
        return xc

    end

    xb = Vector{Tuple{Float64,Float64}}(undef, 0)
    for xi in x
        xi_c = sample_circle(xi, offset)
        xb = vcat(xb, xi_c)
    end

    xb = get_convex_hull(xb,h*0.1)

    return xb

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

"""
Determines whether a set of points is inside a polygon.

# Arguments
- `points::Vector{Tuple{Float64, Float64}}`: A vector of points (x, y) to check.
- `polygon::Vector{Tuple{Float64, Float64}}`: A vector of points (x, y) representing the vertices of the polygon.

# Returns
- `Vector{Bool}`: A vector of booleans where `true` indicates the corresponding point is inside the polygon.
"""
function point_in_polygon(point::Tuple{Float64, Float64}, polygon)
    x, y = point
    n = length(polygon)
    inside = false
    j = n
    for i in 1:n
        xi, yi = polygon[i]
        xj, yj = polygon[j]
        if (yi > y) != (yj > y) && (x < (xj - xi) * (y - yi) / (yj - yi) + xi)
            inside = !inside
        end
        j = i
    end
    return inside
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

function add_geometry_1d(points, h, tag_1d = 0, tag_0d = 0, embed = false)
    
    np = length(points)
    tags = Vector{Int}(undef, np-1)
    n = (points[1] == points[end]) ? np-1 : np

    tags_0d = add_geometry_0d(points[1:n], h, tag_0d)
    for i = 1:np-1
        start = i
        stop = (i + 1 - 1)%n + 1
        tags[i] = gmsh.model.geo.addLine(tags_0d[start], tags_0d[stop], tag_1d + i)
    end
    if embed
        gmsh.model.geo.synchronize()
        gmsh.model.mesh.embed(1, tags, 2, 1)
    end

    return tags, tags_0d

end

function add_geometry_2d(points, h, tag_2d = 0, tag_1d = 0, tag_0d = 0, embed = false)

    tags_1d, tags_0d = add_geometry_1d(points, h, tag_1d, tag_0d, embed)
    tag = gmsh.model.geo.addCurveLoop(tags_1d)
    gmsh.model.geo.addPlaneSurface([tag], tag_2d + 1)

    return tag, tags_1d, tags_0d

end