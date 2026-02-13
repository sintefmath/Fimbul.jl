function extruded_mesh(cell_constraints, depths;
    boundary = missing, offset = missing, offset_rel = 5,
    hxy_min = missing, hxy_max = missing, hz = missing,
    interpolation = :default,
    dist_min_factor = 1.1, dist_max_factor = 0.75,
    recombine_to_quads = true,
    info_level = -1
    )

    depths[1] == 0.0 || error("First depth must be 0.0")
    issorted(depths) || error("Depths must be sorted in ascending order")

    # Initialize Gmsh
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", info_level < 0 ? 0 : 1)
    gmsh.option.setNumber("General.Verbosity", max(info_level, 0))
    gmsh.clear()
    gmsh.model.add("extruded_mesh")

    # ## Get and set domain metrics
    # Compute radius and center of mass of well region
    x_cc = vcat(cell_constraints...)

    min_cc_distance, max_cc_distance = min_max_distance(x_cc)

    if ismissing(boundary)
        if ismissing(offset)
            @assert !ismissing(offset_rel)
            offset = max_cc_distance/2*offset_rel
        else
            @assert ismissing(offset_rel)
        end
        xb_outer = offset_boundary(x_cc, offset; n = 12)
    else
        xb_outer = boundary
        offset = min_distance(vcat(xb_outer, x_cc...))
    end
    if norm(xb_outer[1] .- xb_outer[end], 2) > 0.0
        xb_outer = push!(xb_outer, xb_outer[1])
    end
    perimeter = curve_measure(xb_outer)

    # Set mesh sizes
    if ismissing(hxy_min)
        hxy_min = min_cc_distance/5
    end
    if ismissing(hxy_max)
        hxy_max = perimeter/6
    end    
    radius_outer = max_distance(xb_outer)/2
    depth = depths[end] - depths[1]

    @assert 0.0 < hxy_min <= hxy_max < perimeter/4 "Please ensure that "*
        "0 < hxy_min = $hxy_min < hxy_max = $hxy_max < perimeter/4 = $(perimeter/4)"

    layer_thicknesss = diff(depths)
    num_layers = length(layer_thicknesss)
    if ismissing(hz)
        hz = layer_thicknesss./2.0
    elseif length(hz) == 1
        hz = fill(hz[1], num_layers)
    end
    length(hz) == num_layers || error("Length of hz ($(length(hz))) must match number of layers ($num_layers)")
    for (i, ht) in enumerate(hz)
        ht <= layer_thicknesss[i] || error("hz[$i] = $ht exceeds layer thickness $(layer_thicknesss[i])")
    end
    
    # ## Create 2D mesh
    # Make 2d domain
    tag2d_brd, tag1d_brd, tag0d_bdr = add_geometry_2d(xb_outer, hxy_max)
    # Embed well coordinates
    tag0d_cc, tag1d_cc = Vector{Int}(undef, 0), Vector{Int}(undef, 0)
    i0, i1 = tag0d_bdr[end], tag1d_brd[end]
    for (i, cc) in enumerate(cell_constraints)
        if length(cc) == 1
            t0d_i = add_geometry_0d(cc, hxy_min, i0, true)
            push!(tag0d_cc, t0d_i...)
            i0 += 1
        else
            t1d_i, t0d_i = add_geometry_1d(cc, hxy_min, i1, i0, true)
            push!(tag1d_cc, t1d_i...)
            i1 += length(t1d_i)
            i0 += length(t0d_i)
        end
    end
    # Set mesh size
    gmsh.model.mesh.field.add("Distance", 1)
    if !isempty(tag0d_cc)
        gmsh.model.mesh.field.setNumbers(1, "PointsList", tag0d_cc)
    end
    if !isempty(tag1d_cc)
        gmsh.model.mesh.field.setNumbers(1, "CurvesList", tag1d_cc)
    end
    gmsh.model.mesh.field.add("Threshold", 2)
    gmsh.model.mesh.field.setNumber(2, "InField", 1)
    gmsh.model.mesh.field.setNumber(2, "SizeMin", hxy_min)
    gmsh.model.mesh.field.setNumber(2, "SizeMax", hxy_max)
    # Cell size transition
    dist_min = min_cc_distance*dist_min_factor
    dist_max = (max_cc_distance/2 + offset)*dist_max_factor
    @assert dist_min < dist_max
    "dist_min must be smaller than dist_max"
    if dist_max > radius_outer
        @warn "Warning: dist_max is larger than outer radius"
    end
    gmsh.model.mesh.field.setNumber(2, "DistMin", dist_min)
    gmsh.model.mesh.field.setNumber(2, "DistMax", dist_max)
    gmsh.model.mesh.field.setAsBackgroundMesh(2)
    # Recombine mesh into quadrilaterals
    if recombine_to_quads
        gmsh.model.geo.mesh.setRecombine(2, 1)
    end

    # ## Extrude to 3D and generate
    z, layers = interpolate_z(depths, hz; interpolation = interpolation)
    z = z[2:end]
    height = z./depth
    num_elements = ones(Int, length(z))

    gmsh.model.geo.extrude([(2, 1)], 0, 0, depth, num_elements, height, true)
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(3)

    # ## Create mesh object
    mesh = Jutul.mesh_from_gmsh(z_is_depth=true)
    gmsh.finalize()

    nl = length(layers)
    nc_2d = Int(number_of_cells(mesh)/nl)
    layers = repeat(layers, nc_2d)

    metrics = (nc_2d = nc_2d, hxy_min = hxy_min, hxy_max = hxy_max, hz = hz, 
        offset = offset, dist_min = dist_min, dist_max = dist_max, 
        boundary = xb_outer)

    return mesh, layers, metrics

end

function interpolate_z(depths, hz;
    interpolation = :default, transition = 1/3, force_interpolation = true)

    n = length(depths)
    @assert n >= 2
    
    hz = length(hz) == 1 ? fill(hz, n-1) : hz
    @assert length(hz) == n-1

    interpolation = (interpolation isa Symbol) ? [interpolation] : interpolation
    interpolation = length(interpolation) == 1 ?
        repeat(interpolation, n-1) : interpolation
    for i = 1:n-1
        if interpolation[i] == :default
            hz_prev = hz[max(i-1,1)]
            hz_curr = hz[i]
            hz_next = hz[min(i+1, n-1)]
            top = hz_curr > 1.5*hz_prev
            bottom = hz_curr > 1.5*hz_next
            if top && bottom
                interpolation[i] = :both
            elseif top
                interpolation[i] = :top
            elseif bottom
                interpolation[i] = :bottom
            else
                interpolation[i] = :nothing
            end
        end
    end
    @assert all([itp ∈ [:top, :bottom, :both, :nothing] for itp in interpolation]) "
        interpolation must be either :top, :bottom, :both or :nothing"
    @assert length(interpolation) == n-1

    @assert 0 < transition <= 0.5

    z = Vector{Float64}(undef, 0)
    layer = Vector{Int}(undef, 0)

    for i = 1:n-1
        z_0, z_1 = depths[i], depths[i+1]
        z_0t, z_1t = z_0, z_1
        thickness = z_1 - z_0

        hz_prev = hz[max(i-1,1)]
        hz_curr = hz[i]
        hz_next = hz[min(i+1, n-1)]
        
        frac = (1 < i < n-1) ? transition : 2/3
        dz = frac*thickness
        do_interpolation = interpolation[i] != :nothing
        if do_interpolation && hz_curr > dz
            msg = "Transition in layer $(i) not feasible (hz = $hz_curr > dz = $dz)"
            if force_interpolation
                hz_curr = dz/2
                @warn msg * ". Forcing interpolation by reducing hz to $hz_curr"
            else
                do_interpolation = false
                @warn msg * ". No interpolation will be done."
            end
        end
        if do_interpolation
            if interpolation[i] ∈ [:top, :both]
                z_0t = z_0 + dz
            end
            if interpolation[i] ∈ [:bottom, :both]
                z_1t = z_1 - dz
            end
        end

        z_top = interpolate_z(z_0, z_0t, hz_prev, hz_curr)
        z_mid = interpolate_z(z_0t, z_1t, hz_curr, hz_curr)
        z_bottom = interpolate_z(z_1t, z_1, hz_curr, hz_next)
        z_i = unique(vcat(z_top[1:end-1], z_mid[1:end-1], z_bottom[1:end]))

        push!(z, z_i[1:end-1]...)
        unique!(z)

        push!(layer, fill(i, length(z_i)-1)...)

    end

    push!(z, depths[end])

    return z, layer

end

function interpolate_z(z_a::Float64, z_b::Float64, dz_a::Float64, dz_b::Float64)

    L = z_b - z_a
    if isapprox(L, 0.0)
        z = [z_a, z_b]

    elseif isapprox(dz_a, dz_b)
        n = max(Int(ceil(L/dz_a))+1,2)
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

    return z

end