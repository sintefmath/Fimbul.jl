function extruded_mesh(well_coordinates, depths;
    radius_outer = missing,
    radius_outer_factor = 5,
    hxy_min = missing, hxy_max = missing, hz = missing,
    dist_min_factor = 1.1, dist_max_factor = 1.5
    )

    # Initialize Gmsh
    gmsh.initialize()
    gmsh.clear()
    gmsh.model.add("extruded_mesh")

    # ## Get and set domain metrics
    # Compute radius and center of mass of well region
    xw = well_coordinates
    xc = (0.0, 0.0)
    diameter, well_distance = -Inf, Inf
    for (i,x) in enumerate(xw)
        xc = xc .+ x
        for (j,y) in enumerate(xw)
            i == j ? continue : nothing
            d_ij = norm(x .- y,2)
            well_distance = min(well_distance, d_ij)
            diameter = max(diameter, d_ij)
        end
    end
    radius = diameter/2
    xc = xc./length(xw)
    # Set outer radius
    if ismissing(radius_outer)
        @assert !ismissing(radius_outer_factor) 
        "Please provide either radius_outer or radius_outer_factor"
        radius_outer = radius*radius_outer_factor
    end
    @assert radius_outer > radius 
    "radius_outer must be larger than the well region diameter"
    # Construct outer boundary
    perimeter = 2*π*radius_outer
    # Set mesh sizes
    if ismissing(hxy_min)
        hxy_min = well_distance/5
    end
    if ismissing(hxy_max)
        hxy_max = perimeter/25
    end
    @assert 0.0 < hxy_min < hxy_max < perimeter/4
    "Please ensure that 0 < hxy_min < hxy_max < perimeter/4"
    if ismissing(hz)
        hz = sum(depths)/50
    end

    # ## Create 2D mesh
    n = Int(ceil(perimeter/hxy_max))
    Δθ = 2*π/n
    for i = 1:n
        x, y = xc[1] + radius_outer*cos(i*Δθ), xc[2] + radius_outer*sin(i*Δθ)
        gmsh.model.geo.addPoint(x, y, 0.0, hxy_max, i)
    end
    for i = 1:n
        gmsh.model.geo.addLine(i, mod(i, n) + 1, i)
    end
    gmsh.model.geo.addCurveLoop(collect(1:n), 1)
    gmsh.model.geo.addPlaneSurface([1], 1)
    # Embed well coordinates
    for (i, x) in enumerate(xw)
        gmsh.model.geo.addPoint(x[1], x[2], 0.0, hxy_min, n + i)
    end
    nw = length(xw)
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.embed(0, collect(n+1:n+length(xw)), 2, 1)
    gmsh.model.geo.addPoint(xc[1], xc[2], 0.0, hxy_min, n + nw + 1)
    # Set mesh size
    gmsh.model.mesh.field.add("Distance", 1)
    gmsh.model.mesh.field.setNumbers(1, "PointsList", [n + nw + 1])
    gmsh.model.mesh.field.add("Threshold", 2)
    gmsh.model.mesh.field.setNumber(2, "InField", 1)
    gmsh.model.mesh.field.setNumber(2, "SizeMin", hxy_min)
    gmsh.model.mesh.field.setNumber(2, "SizeMax", hxy_max)
    dist_min = dist_min_factor*radius
    dist_max = min(dist_max_factor*radius, radius_outer*0.9)
    @assert dist_min < dist_max
    "dist_min must be smaller than dist_max"
    if dist_max > radius_outer
        @warn "Warning: dist_max is larger than outer radius"
    end
    gmsh.model.mesh.field.setNumber(2, "DistMin", dist_min)
    gmsh.model.mesh.field.setNumber(2, "DistMax", dist_max)
    gmsh.model.mesh.field.setAsBackgroundMesh(2)
    # Recombine mesh into quadrilaterals
    gmsh.model.geo.mesh.setRecombine(2, 1)

    println("dist_min: $dist_min, dist_max: $dist_max, hxy_min: $hxy_min, hxy_max: $hxy_max, radius: $radius, radius_outer: $radius_outer") 

    # ## Extrude to 3D and generate

    z, layer = interpolate_z(depths, hz)
    z = z[2:end]
    depth = sum(depths)
    height = z./depth
    println()
    num_elements = ones(Int, length(z))
    println("z: $z")
    println("height: $height")

    gmsh.model.geo.extrude([(2, 1)], 0, 0, depth, num_elements, height, true)
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(3)

    # ## Create mesh object
    mesh = Jutul.mesh_from_gmsh(z_is_depth=true)
    gmsh.finalize()

    metrics = (radius = radius, radius_outer = radius_outer,
        hxy_min = hxy_min, hxy_max = hxy_max, hz = hz,
        dist_min = dist_min, dist_max = dist_max)

    return mesh, metrics

end

function interpolate_z(depths, hz;
    interpolation = :default, transition = 0.33)

    n = length(depths)
    @assert n >= 2
    
    hz = length(hz) == 1 ? fill(hz, n-1) : hz
    @assert length(hz) == n-1

    interpolation = (interpolation isa Symbol) ? [interpolation] : interpolation
    interpolation = length(interpolation) == 1 ?
        fill(interpolation, n-1) : interpolation
    @assert length(hz) == n-1

    z = Vector{Float64}(undef, 0)
    layer = Vector{Int}(undef, 0)

    thickness = depths[2:end] - depths[1:end-1]
    hz = min.(hz, 0.5*thickness)

    for i = 1:n-1
        z_0, z_1 = depths[i], depths[i+1]

        hz_prev = hz[max(i-1,1)]
        hz_curr = hz[i]
        hz_next = hz[min(i+1, n-1)]

        z_0t, z_1t = z_0, z_1
        hz_0t, hz_1t = hz_curr, hz_curr
        hz_0, hz_1 = hz_prev, hz_next
        thickness = depths[i+1] - depths[i]

        if interpolation[i] == :default

            dz = transition*thickness
            if hz_curr > 1.5*hz_prev
                frac = i < n-1 ? transition : 1.0
                dz = frac*thickness
                if hz_curr < 1.0*dz
                    z_0t = z_0 + dz
                else
                    @warn "Transition from layer $(i-1) to layer $(i) not feasible"
                end
            end
            if hz_curr > 1.5*hz_next
                frac = i > 1 ? transition : 1.0
                dz = frac*thickness
                if hz_curr < 1.0*dz
                    z_1t = z_1 - dz
                else
                    @warn "Transition from layer $i to layer $(i+1) not feasible"
                end
            end
            
        else
            @error "Unknown interpolation method $interpolation"
        end

        z_top = interpolate_z(z_0, z_0t, hz_0, hz_0t)
        z_mid = interpolate_z(z_0t, z_1t, hz_curr, hz_curr)
        z_bottom = interpolate_z(z_1t, z_1, hz_1t, hz_1)
        z_i = unique(vcat(z_top, z_mid, z_bottom))

        push!(z, z_i[1:end-1]...)
        unique!(z)

        push!(layer, fill(i, length(z_i)-1)...)

    end

    push!(z, depths[end])

    return z, layer

end

function interpolate_z(z_a::Float64, z_b::Float64, dz_a::Float64, dz_b::Float64)

    L = z_b - z_a
    println()
    println("z_a: $z_a, z_b: $z_b, L: $L, dz_a: $dz_a, dz_b: $dz_b")
    if isapprox(L, 0.0)
        z = [z_a, z_b]

    elseif isapprox(dz_a, dz_b)
        n = max(Int(ceil(L/dz_a))+1,2)
        println("n: $n")
        z = collect(range(z_a, z_b, length=n))

    else
        α = (L-dz_a)/(L-dz_b)
        K = Int(ceil(log(dz_b/dz_a)/log(α)))
        α = (dz_b/dz_a)^(1/K)
        dz = dz_a*α.^(0:K)
        rem = sum(dz) - L
        dz .-= rem.*dz./sum(dz)
        z = z_a .+ cumsum(vcat(0, dz))

    end

    @assert isapprox(z[1], z_a) && isapprox(z[end], z_b)
    println("z: $z")

    return z

end