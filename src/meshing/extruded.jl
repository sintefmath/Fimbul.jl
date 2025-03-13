function extruded_mesh(well_coordinates, depth;
    radius_outer = missing,
    radius_outer_factor = 5,
    hxy_min = missing, hxy_max = missing, hz = depth/50,
    dist_min_factor = 1.1, dist_max_factor = 1.5)

    # Initialize Gmsh
    gmsh.initialize()
    gmsh.clear()
    gmsh.model.add("extruded_mesh")

    # ## Get and set doain metrics
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
    num_elements = Int(ceil(depth/hz))
    gmsh.model.geo.extrude([(2, 1)], 0, 0, depth, [num_elements], [1.0], true)
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