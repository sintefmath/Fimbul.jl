function fibonacci_pattern_2d(num_points; spacing = 5.0, radius = missing)
    
    # Define coordinate function using golden ratio
    N = num_points
    φ = (1 + sqrt(5))/2 # ≈ golden ratio
    Δθ = 2*π/φ^2
    
    # Generate points as matrix
    points = Matrix{Float64}(undef, N, 2)
    
    for n = 1:N
        r_n = sqrt(n/N)
        θ_n = n * Δθ
        points[n, 1] = r_n * cos(θ_n)
        points[n, 2] = r_n * sin(θ_n)
    end
    
    # Set radius from spacing if not provided
    if !ismissing(spacing)
        @assert ismissing(radius) "Please provide either radius or spacing"
        # Calculate spacing between first two points
        p1 = points[1, :]
        p2 = points[2, :]
        current_spacing = norm(p1 - p2)
        radius = 1.3513 * spacing / current_spacing
    elseif !ismissing(radius)
        @assert ismissing(spacing) "Please provide either radius or spacing"
    else
        radius = 1.0
    end
    
    # Scale points by radius
    points .*= radius
    
    return points
end

function get_convex_hull(points::AbstractMatrix, tol = eps(Float64))
    # Ensure there are at least 3 points
    @assert size(points, 1) >= 3 "At least 3 points are required to compute a convex hull."

    # Helper function to compute the cross product of vectors (p1 -> p2) and (p1 -> p3)
    function cross_product(p1, p2, p3)
        (p2[1] - p1[1]) * (p3[2] - p1[2]) - (p2[2] - p1[2]) * (p3[1] - p1[1])
    end

    n_points = size(points, 1)
    
    # Start with the leftmost point
    start = argmin(points[:, 1])  # Index of the leftmost point
    current = start

    hull_indices = Int[]
    
    while true
        push!(hull_indices, current)  # Add the current point to the hull
        next_point = mod1(current + 1, n_points)  # Next candidate point

        for i in 1:n_points
            if cross_product(points[current, :], points[next_point, :], points[i, :]) < 0
                next_point = i  # Update the next point if it is more counterclockwise
            end
        end

        current = next_point  # Move to the next point
        if current == start
            break  # Stop when we return to the starting point
        end
    end
    
    # Extract hull points as matrix
    hull_points = points[hull_indices, :]

    hull_simplified, _ = ramer_douglas_peucker(hull_points; tolerance = tol)

    return hull_simplified
end

function point_to_line_distance_2d(point::AbstractVector, line_start::AbstractVector, line_end::AbstractVector)
    """Calculate perpendicular distance from point to line segment in 2D"""
    x0, y0 = point[1], point[2]
    x1, y1 = line_start[1], line_start[2]
    x2, y2 = line_end[1], line_end[2]
    
    # Line vector
    dx = x2 - x1
    dy = y2 - y1
    
    # If line segment has zero length, return distance to start point
    if dx == 0 && dy == 0
        return norm(point - line_start)
    end
    
    # Calculate perpendicular distance using cross product formula
    numerator = abs(dy * x0 - dx * y0 + x2 * y1 - y2 * x1)
    denominator = sqrt(dx^2 + dy^2)
    
    return numerator / denominator
end

function point_to_line_distance_3d(point::AbstractVector, line_start::AbstractVector, line_end::AbstractVector)
    """Calculate perpendicular distance from point to line segment in 3D"""
    # Vector from line start to point
    v = point - line_start
    # Line direction vector
    u = line_end - line_start
    
    # If line segment has zero length, return distance to start point
    if norm(u) ≈ 0
        return norm(v)
    end
    
    # Cross product to get perpendicular component
    cross_prod = cross(v, u)
    
    # Distance is magnitude of cross product divided by line length
    return norm(cross_prod) / norm(u)
end

function ramer_douglas_peucker(points::AbstractMatrix, tolerance::Float64)
    """
    Simplify a curve using the Ramer-Douglas-Peucker algorithm.
    
    Parameters:
    - points: Matrix (Nx2 or Nx3) representing 2D (x,y) or 3D (x,y,z) points
    - tolerance: Maximum allowed deviation from original curve
    
    Returns:
    - simplified_points: Matrix of simplified points
    - indices: Vector of indices into original points array
    """
    
    n_points, n_dims = size(points)
    
    # Need at least 2 points
    if n_points <= 2
        return copy(points), collect(1:n_points)
    end
    
    # Recursive helper function
    function rdp_recursive(start_idx::Int, end_idx::Int)
        # Find the point with maximum distance from line segment
        max_dist = 0.0
        max_idx = 0
        
        start_point = points[start_idx, :]
        end_point = points[end_idx, :]
        
        for i in (start_idx + 1):(end_idx - 1)
            point = points[i, :]
            
            # Calculate distance based on dimensionality
            if n_dims == 2
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
    indices = rdp_recursive(1, n_points)
    
    # Extract simplified points
    simplified_points = points[indices, :]
    
    return simplified_points, indices
end

# Convenience function with default tolerance
function ramer_douglas_peucker(points::AbstractMatrix; tolerance::Float64 = 1.0)
    return ramer_douglas_peucker(points, tolerance)
end

function offset_boundary(x::AbstractMatrix, offset; h = missing, n = 6)

    @assert offset > 0 "Offset must be larger than 0"

    # Number of points to sample around each point
    if !ismissing(h)
        n = max(Int(ceil(2π*offset/h)) + 1, n)
    else
        h = 2π*offset/n
    end

    function sample_circle(x0::AbstractVector, r)
        θ = LinRange(0, 2π, n+1)
        # Create matrix for circle points
        circle_points = Matrix{Float64}(undef, n, 2)
        for i = 1:n
            circle_points[i, 1] = x0[1] + r * cos(θ[i])
            circle_points[i, 2] = x0[2] + r * sin(θ[i])
        end
        return circle_points
    end

    # Generate all circle points
    n_points = size(x, 1)
    total_circle_points = n * n_points
    xb = Matrix{Float64}(undef, total_circle_points, 2)
    
    idx = 1
    for i = 1:n_points
        circle = sample_circle(x[i, :], offset)
        xb[idx:idx+n-1, :] = circle
        idx += n
    end

    # Get convex hull
    hull = get_convex_hull(xb, h*0.1)

    return hull
end

function curve_measure(x::AbstractMatrix)
    n_points = size(x, 1)
    total_length = 0.0
    for i = 1:n_points-1
        total_length += norm(x[i+1, :] - x[i, :])
    end
    return total_length
end

function min_max_distance(x::AbstractMatrix)
    n_points = size(x, 1)
    d_min, d_max = Inf, -Inf
    
    for i in 1:n_points
        for j in 1:n_points
            i == j && continue
            d_ij = norm(x[i, :] - x[j, :])
            d_min = min(d_min, d_ij)
            d_max = max(d_max, d_ij)
        end
    end
    return d_min, d_max
end

function min_distance(x::AbstractMatrix)
    return min_max_distance(x)[1]
end

function max_distance(x::AbstractMatrix)
    return min_max_distance(x)[2]
end

"""
Determines whether a point is inside a polygon.

# Arguments
- `point::AbstractVector`: A point (x, y) to check.
- `polygon::AbstractMatrix`: Matrix (Nx2) representing the vertices of the polygon.

# Returns
- `Bool`: `true` if the point is inside the polygon.
"""
function point_in_polygon(point::AbstractVector, polygon::AbstractMatrix)
    x, y = point[1], point[2]
    n = size(polygon, 1)
    inside = false
    j = n
    for i in 1:n
        xi, yi = polygon[i, 1], polygon[i, 2]
        xj, yj = polygon[j, 1], polygon[j, 2]
        if (yi > y) != (yj > y) && (x < (xj - xi) * (y - yi) / (yj - yi) + xi)
            inside = !inside
        end
        j = i
    end
    return inside
end

"""
Determines whether multiple points are inside a polygon.

# Arguments
- `points::AbstractMatrix`: Matrix (Mx2) of points to check.
- `polygon::AbstractMatrix`: Matrix (Nx2) representing the vertices of the polygon.

# Returns
- `Vector{Bool}`: Vector where `true` indicates the point is inside the polygon.
"""
function points_in_polygon(points::AbstractMatrix, polygon::AbstractMatrix)
    return [point_in_polygon(points[i, :], polygon) for i in 1:size(points, 1)]
end

function add_geometry_0d(points::AbstractMatrix, h, tag = 0, embed = false)
    n_points = size(points, 1)
    tags = Vector{Int}(undef, n_points)
    for i = 1:n_points
        x, y = points[i, 1], points[i, 2]
        tags[i] = gmsh.model.geo.addPoint(x, y, 0.0, h, tag + i)
    end
    if embed
        gmsh.model.geo.synchronize()
        gmsh.model.mesh.embed(0, tags, 2, 1)
    end
    return tags
end

function add_geometry_1d(points::AbstractMatrix, h, tag_1d = 0, tag_0d = 0, embed = false)
    n_points = size(points, 1)
    tags = Vector{Int}(undef, n_points-1)
    
    # Check if curve is closed (first and last points are the same)
    is_closed = norm(points[1, :] - points[end, :]) < 1e-12
    n = is_closed ? n_points-1 : n_points

    tags_0d = add_geometry_0d(points[1:n, :], h, tag_0d)
    for i = 1:n_points-1
        start = i
        stop = (i + 1 - 1) % n + 1
        tags[i] = gmsh.model.geo.addLine(tags_0d[start], tags_0d[stop], tag_1d + i)
    end
    if embed
        gmsh.model.geo.synchronize()
        gmsh.model.mesh.embed(1, tags, 2, 1)
    end
    return tags, tags_0d
end

function add_geometry_2d(points::AbstractMatrix, h, tag_2d = 0, tag_1d = 0, tag_0d = 0, embed = false)
    tags_1d, tags_0d = add_geometry_1d(points, h, tag_1d, tag_0d, embed)
    tag = gmsh.model.geo.addCurveLoop(tags_1d)
    gmsh.model.geo.addPlaneSurface([tag], tag_2d + 1)
    return tag, tags_1d, tags_0d
end