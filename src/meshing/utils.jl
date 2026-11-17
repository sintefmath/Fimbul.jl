"""
    fibonacci_pattern_2d(num_points; spacing = 5.0, radius = missing)

Points in a sunflower (Vogel spiral) pattern, returned as a `2 x num_points`
matrix. Point `n` is placed at angle `n*Δθ`, with `Δθ = 2π/φ^2` and `φ` the
golden ratio, and at radius `radius*sqrt((n - 1/2)/(num_points - 1/2))`. The
half-step offset closes the gap that otherwise forms at the center of the
pattern. Points `F` apart in this order, for a Fibonacci number `F`, lie on a
spiral ("spine") through the pattern.

The pattern is scaled either by `spacing`, the median distance from a point to
its nearest neighbor, or by `radius`, the distance from the center to the
outermost point. A single point is placed at the origin.
"""
function fibonacci_pattern_2d(num_points; spacing = 5.0, radius = missing)
    N = num_points
    N >= 1 || error("num_points must be at least 1, got $N")
    Δθ = 2π/Base.MathConstants.golden^2
    points = Matrix{Float64}(undef, 2, N)
    for n = 1:N
        r_n = N == 1 ? 0.0 : sqrt((n - 0.5)/(N - 0.5))
        θ_n = n*Δθ
        points[1, n] = r_n*cos(θ_n)
        points[2, n] = r_n*sin(θ_n)
    end
    if !ismissing(spacing)
        ismissing(radius) || error("Please provide either radius or spacing")
        N == 1 && return points
        scale = spacing/median(nearest_neighbor_distances(points))
    elseif !ismissing(radius)
        scale = radius
    else
        scale = 1.0
    end
    return points.*scale
end

# Distance from each column of x to its nearest neighbor among the others
function nearest_neighbor_distances(x::AbstractMatrix)
    n = size(x, 2)
    return [minimum(norm(x[:, i] - x[:, j]) for j in 1:n if j != i) for i in 1:n]
end

function get_convex_hull(points::AbstractMatrix, tol = eps(Float64))
    # Ensure there are at least 3 points
    @assert size(points, 2) >= 3 "At least 3 points are required to compute a convex hull."

    # Helper function to compute the cross product of vectors (p1 -> p2) and (p1 -> p3)
    function cross_product(p1, p2, p3)
        (p2[1] - p1[1]) * (p3[2] - p1[2]) - (p2[2] - p1[2]) * (p3[1] - p1[1])
    end

    n_points = size(points, 2)
    
    # Start with the leftmost point
    start = argmin(points[1, :])  # Index of the leftmost point
    current = start

    hull_indices = Int[]
    while true
        push!(hull_indices, current)  # Add the current point to the hull
        next_point = mod1(current+1, n_points)  # Next candidate point

        for i in 1:n_points
            if cross_product(points[:, current], points[:, next_point], points[:, i]) < 0
                next_point = i  # Update the next point if it is more counterclockwise
            end
        end

        current = next_point  # Move to the next point
        if current == start || current ∈ hull_indices
            break  # Stop when we return to the starting point
        end
    end
    
    # Extract hull points as matrix
    hull_points = points[:, hull_indices]

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
    - points: Matrix (2×N or 3×N) representing 2D (x,y) or 3D (x,y,z) points
    - tolerance: Maximum allowed deviation from original curve
    
    Returns:
    - simplified_points: Matrix of simplified points
    - indices: Vector of indices into original points array
    """
    
    n_dims, n_points = size(points)
    
    # Need at least 2 points
    if n_points <= 2
        return copy(points), collect(1:n_points)
    end
    
    # Recursive helper function
    function rdp_recursive(start_idx::Int, end_idx::Int)
        # Find the point with maximum distance from line segment
        max_dist = 0.0
        max_idx = 0
        
        start_point = points[:, start_idx]
        end_point = points[:, end_idx]
        
        for i in (start_idx + 1):(end_idx - 1)
            point = points[:, i]
            
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
    simplified_points = points[:, indices]
    
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
        # Create matrix for circle points (2×n)
        circle_points = Matrix{Float64}(undef, 2, n)
        for i = 1:n
            circle_points[1, i] = x0[1] + r * cos(θ[i])
            circle_points[2, i] = x0[2] + r * sin(θ[i])
        end
        return circle_points
    end

    # Generate all circle points
    n_points = size(x, 2)
    total_circle_points = n * n_points
    xb = Matrix{Float64}(undef, 2, total_circle_points)
    
    idx = 1
    for i = 1:n_points
        circle = sample_circle(x[:, i], offset)
        xb[:, idx:idx+n-1] = circle
        idx += n
    end

    # Get convex hull
    hull = get_convex_hull(xb, h*0.1)

    return hull
end

function curve_measure(x::AbstractMatrix)
    n_points = size(x, 2)
    total_length = 0.0
    for i = 1:n_points-1
        total_length += norm(x[:, i+1] - x[:, i])
    end
    return total_length
end

function min_max_distance(x::AbstractMatrix)
    n_points = size(x, 2)
    d_min, d_max = Inf, -Inf
    
    for i in 1:n_points
        for j in 1:n_points
            i == j && continue
            d_ij = norm(x[:, i] - x[:, j])
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

function field_from_points(sector_division::Symbol,xy::AbstractMatrix, num_sectors::Int, depths::AbstractVector;
    block_size = missing)

    num_wells = size(xy, 2)
    well_coords = Vector{Matrix{Float64}}(undef, num_wells)
    for i in 1:num_wells
        x_top, y_top, z_top = xy[1, i], xy[2, i], 0.0 + 1e-3
        x_bottom, y_bottom, z_bottom = xy[1, i], xy[2, i], depths[end] - 1e-3
        well_coords[i] = permutedims([x_top y_top z_top; x_bottom y_bottom z_bottom])
    end
    if sector_division == :angular
        sector_indices = group_into_sectors_angular(xy, num_sectors)
    elseif sector_division == :cartesian
        sector_indices = group_into_sectors_cartesian(xy, num_sectors; block_size = block_size)
    elseif sector_division == :sequential
        # Consecutive groups of wells, in the order given
        sector_indices = [collect(a:b) for (a, b) in equal_index_ranges(num_wells, num_sectors)]
    elseif sector_division == :spiral
        # Bundles of sunflower spines, with the wells in sunflower order
        sector_indices = group_into_sectors_spiral(num_wells, num_sectors)
    else
        error("Unknown sector division: $sector_division")
    end

    return [well_coords[idx] for idx in sector_indices]

end

"""
    group_into_sectors_spiral(num_wells, num_sectors)

Group the points of a sunflower pattern (see [`fibonacci_pattern_2d`](@ref)),
given in sunflower order, into sectors made of adjacent spines. Points `F`
apart lie on the same spine, for a Fibonacci number `F`, so that a family of
`F` spines holds all points. The spines are ordered counter-clockwise and split
into `num_sectors` bundles of adjacent spines with as equal a number of points
as possible. Within a bundle, the points are ordered from the center
outwards: the first point of each spine, counter-clockwise, then the second
point of each spine, and so on.

The spine family is chosen among the Fibonacci numbers `num_sectors ≤ F ≤
num_wells/2` to give the most even sectors, and among those the family whose
spines connect the nearest points.
"""
function group_into_sectors_spiral(num_wells::Int, num_sectors::Int)
    N = num_wells
    num_sectors <= N || error("Cannot divide $N wells into $num_sectors sectors.")
    num_sectors == 1 && return [collect(1:N)]
    Δθ = 2π/Base.MathConstants.golden^2
    unit_points = fibonacci_pattern_2d(N; spacing = missing, radius = 1.0)
    # Candidate spine families
    fib = [1, 2]
    while fib[end] < N
        push!(fib, fib[end] + fib[end-1])
    end
    candidates = filter(F -> num_sectors <= F <= max(num_sectors, N ÷ 2), unique(fib))
    isempty(candidates) && (candidates = [num_sectors])
    best = nothing
    for F in candidates
        # Spines j = 1, ..., F hold the points j, j + F, j + 2F, ...
        spines = [collect(j:F:N) for j in 1:F]
        # Counter-clockwise order of the spines, compared at a common position
        # along the spiral: spine j turns by δ for every step of F points
        δ = mod(F*Δθ + π, 2π) - π
        φ = [mod(j*Δθ - δ*j/F, 2π) for j in 1:F]
        order = sortperm(φ)
        counts = length.(spines[order])
        # Split the spines into contiguous bundles, trying every starting spine
        for start in 0:F-1
            c = circshift(counts, -start)
            cuts = spine_bundle_cuts(c, num_sectors)
            isnothing(cuts) && continue
            sizes = [sum(c[a:b]) for (a, b) in cuts]
            imbalance = maximum(sizes) - minimum(sizes)
            # Mean distance between consecutive points on a spine
            step = sum(norm(unit_points[:, n + F] - unit_points[:, n]) for n in 1:N-F; init = 0.0)/max(N - F, 1)
            score = (imbalance, step, start)
            if isnothing(best) || score < best[1]
                best = (score, spines, circshift(order, -start), cuts)
            end
        end
    end
    _, spines, order, cuts = best
    # Each well has a level (its position along its spine, from the center)
    # and a rank (the counter-clockwise position of its spine)
    sectors = [[(level, rank, spines[order[rank]][level])
        for rank in a:b for level in eachindex(spines[order[rank]])]
        for (a, b) in cuts]
    rebalance_spiral_sectors!(sectors, length(order))
    return [last.(sort(sector)) for sector in sectors]
end

# Even out the number of wells in the sectors, to within one, by passing the
# outermost well on the boundary spine of a sector to its neighbor, hop by hop
# around the field from the largest to the smallest sector. Wells passed across
# the start of the counter-clockwise spine order get their rank shifted by the
# number of spines, so that the ranks within a sector stay contiguous.
function rebalance_spiral_sectors!(sectors, num_spines)
    K = length(sectors)
    for _ in 1:sum(length, sectors)
        sizes = length.(sectors)
        i, j = argmax(sizes), argmin(sizes)
        sizes[i] - sizes[j] <= 1 && break
        # Shortest way around the ring of sectors from i to j
        forward = mod(j - i, K) <= mod(i - j, K)
        step = forward ? 1 : -1
        k = i
        while k != j
            next = mod1(k + step, K)
            sector = sectors[k]
            # Boundary spine towards the neighbor, and its outermost well
            boundary = forward ? maximum(w[2] for w in sector) : minimum(w[2] for w in sector)
            on_boundary = findall(w -> w[2] == boundary, sector)
            pick = on_boundary[argmax([sector[m][1] for m in on_boundary])]
            level, rank, well = sector[pick]
            deleteat!(sector, pick)
            if forward && next < k
                rank -= num_spines
            elseif !forward && next > k
                rank += num_spines
            end
            push!(sectors[next], (level, rank, well))
            k = next
        end
    end
    return sectors
end

# Split items with the given counts into num_groups contiguous, non-empty
# groups, minimizing the sum of squared deviations of the group totals from an
# equal share (dynamic programming over the cut positions). Returns the index
# ranges, or nothing if there are fewer items than groups.
function spine_bundle_cuts(counts, num_groups)
    n = length(counts)
    n >= num_groups || return nothing
    target = sum(counts)/num_groups
    cum = [0; cumsum(counts)]
    cost(i, j) = (cum[j+1] - cum[i] - target)^2  # group of items i..j
    # best[k, j]: minimal cost of splitting items 1..j into k groups
    best = fill(Inf, num_groups, n)
    prev = zeros(Int, num_groups, n)
    for j in 1:n
        best[1, j] = cost(1, j)
    end
    for k in 2:num_groups, j in k:n, i in k:j
        c = best[k-1, i-1] + cost(i, j)
        if c < best[k, j]
            best[k, j] = c
            prev[k, j] = i
        end
    end
    ranges = Vector{Tuple{Int, Int}}(undef, num_groups)
    j = n
    for k in num_groups:-1:1
        i = k == 1 ? 1 : prev[k, j]
        ranges[k] = (i, j)
        j = i - 1
    end
    return ranges
end

function group_into_sectors_angular(xy::AbstractMatrix, num_sectors::Int)
    # Divide the well field into roughly equal angle sectors

    n = size(xy, 2)
    r = sqrt.(xy[1,:].^2 .+ xy[2,:].^2)
    θ = atan.(xy[2,:], xy[1,:]) .+ π
    order_θ = sortperm(θ)
    wells_per_sector = div(n, num_sectors)
    rem = n - wells_per_sector*num_sectors
    wells_per_sector = fill(wells_per_sector, num_sectors)
    wells_per_sector[1:rem] .+= 1

    sector_indices = Vector{Vector{Int}}()
    wtot = 0
    for sno in 1:num_sectors
        idx = order_θ[(1:wells_per_sector[sno]) .+ wtot]
        wtot += wells_per_sector[sno]
        order_r = sortperm(r[idx])
        push!(sector_indices, idx[order_r])
    end

    return sector_indices
end

function equal_index_ranges(total::Int, n_groups::Int)
    # Split `total` items into n_groups contiguous, as-equal-as-possible chunks
    base = div(total, n_groups)
    rem = total - base*n_groups
    sizes = fill(base, n_groups)
    sizes[1:rem] .+= 1

    ranges = Vector{Tuple{Int,Int}}()
    s = 0
    for sz in sizes
        push!(ranges, (s+1, s+sz))
        s += sz
    end
    return ranges
end

function group_into_sectors_cartesian(xy::AbstractMatrix, num_sectors::Int;
    block_size = missing)
    # Divide the well field into a grid of roughly equal Cartesian blocks
    # The order of the wells within a sector is done in a zigzag order, 
    # where the first column is populated from the top, then the next column is populated 
    # from the bottom etc. 

    n = size(xy, 2)
    cols = sort(unique(xy[1, :]))
    rows = sort(unique(xy[2, :]))
    n_cols_total, n_rows_total = length(cols), length(rows)

    # Split the rows and columns into bands, n_row_bands*n_col_bands =
    # num_sectors, choosing the split that gives the most even number of wells
    # per sector, and among those the one closest to square blocks
    function band_split(n_row_bands, n_col_bands)
        col_bands = equal_index_ranges(n_cols_total, n_col_bands)
        row_bands = equal_index_ranges(n_rows_total, n_row_bands)
        counts = [count(i -> cols[c0] <= xy[1, i] <= cols[c1] &&
            rows[r0] <= xy[2, i] <= rows[r1], 1:n)
            for (c0, c1) in col_bands for (r0, r1) in row_bands]
        aspect = abs(log((n_rows_total/n_row_bands)/(n_cols_total/n_col_bands)))
        return (maximum(counts) - minimum(counts), aspect), row_bands, col_bands
    end
    if ismissing(block_size)
        splits = [band_split(d, num_sectors ÷ d) for d in 1:num_sectors
            if num_sectors % d == 0 && d <= n_rows_total && num_sectors ÷ d <= n_cols_total]
        isempty(splits) && error("Cannot divide a $n_cols_total x $n_rows_total \
            grid of wells into $num_sectors sectors.")
        _, row_bands, col_bands = splits[argmin(first.(splits))]
    else
        # Blocks of block_size[1] columns and block_size[2] rows
        length(block_size) == 2 && all(>(0), block_size) ||
            error("block_size must be (nx, ny) with nx, ny > 0, got $block_size.")
        bx, by = block_size
        col_bands = [(c, min(c + bx - 1, n_cols_total)) for c in 1:bx:n_cols_total]
        row_bands = [(r, min(r + by - 1, n_rows_total)) for r in 1:by:n_rows_total]
        num_blocks = length(col_bands)*length(row_bands)
        num_blocks == num_sectors || error("Blocks of $bx x $by wells divide \
            the $n_cols_total x $n_rows_total grid of wells into $num_blocks \
            sectors, but num_sectors = $num_sectors.")
    end

    sector_indices = Vector{Vector{Int}}()
    for (c0, c1) in col_bands
        band_cols = cols[c0:c1]
        col_rank = Dict(xc => j for (j, xc) in enumerate(band_cols))
        for (r0, r1) in row_bands
            band_rows = rows[r0:r1]
            idx = [i for i in 1:n if xy[1, i] in band_cols && xy[2, i] in band_rows]
            idx = sort(idx, by = i -> (xy[1, i], iseven(col_rank[xy[1, i]]) ? -xy[2, i] : xy[2, i]))
            push!(sector_indices, idx)
        end
    end
    return sector_indices
end

"""
Determines whether a point is inside a polygon. Uses the "even-odd" rule to find whether a point is 
inside or outside a polygon; if a horizontal ray from the point to infinity crosses an edge of the polygon 
an odd number of times, then a point is inside; else, it is outside. 

# Arguments
- `point::AbstractVector`: A point (x, y) to check.
- `polygon::AbstractMatrix`: Matrix (2×N) representing the vertices of the polygon.

# Returns
- `Bool`: `true` if the point is inside the polygon.
"""
function point_in_polygon(point::AbstractVector, polygon::AbstractMatrix)
    x, y = point[1], point[2]
    n = size(polygon, 2)
    inside = false
    j = n
    for i in 1:n
        xi, yi = polygon[1, i], polygon[2, i]
        xj, yj = polygon[1, j], polygon[2, j]
        if (yi > y) != (yj > y) && (x < (xj - xi) * (y - yi) / (yj - yi) + xi)
            inside = !inside
        end
        j = i
    end
    return inside
end

function points_in_polygon(points::AbstractMatrix, polygon::AbstractMatrix)
    return [point_in_polygon(points[:, i], polygon) for i in 1:size(points, 2)]
end

function add_geometry_0d(points::AbstractMatrix, h, tag = 0, embed = false)
    n_points = size(points, 1)  # For N×2 matrix (transposed from 2×N)
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
    n_points = size(points, 1)  # For N×2 matrix (transposed from 2×N)
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