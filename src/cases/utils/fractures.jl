function add_fractures(mesh, centers, normals, matrix_properties, fracture_properties)

    fracture_faces = Int[]
    fracture_cells = Int[]
    fracture_properties_exp = Dict(keys(fracture_properties) .=> [Float64[] for _ in keys(fracture_properties)])

    for (i, (normal, center)) in enumerate(zip(normals, centers))
        is_frac = falses(number_of_faces(mesh))
        is_frac[fracture_faces] .= true
        plane = PlaneCut(center, normal)
        mesh, info = cut_mesh(mesh, plane; min_cut_fraction = 0.0, extra_out=true)
        # Update fracture face vector
        face_index = filter(x->x>0, info["face_index"])
        is_frac = is_frac[face_index]
        new_faces = findall(info["face_index"] .== 0)
        if any(is_frac)
            old_to_new_faces = face_index[is_frac]
            face_to_cell = Dict((f, c) for (f, c) in zip(face_index, 1:length(face_index)))
            old_to_new_cells = face_to_cell[old_to_new_faces]
        else
            old_to_new_cells = Int[]
        end
        # Expand fracture porperties
        for (key, values) in fracture_properties
            value_exp = fracture_properties_exp[key]
            value_exp = value_exp[old_to_new_cells]
            append!(fracture_properties_exp[key], fill(values[i], length(new_faces)))
        end
        fracture_faces = findall(is_frac)
        append!(fracture_faces, new_faces)
        append!(fracture_cells, collect(1:length(new_faces)))
        println("Fracture $i: $(length(fracture_faces)) faces, new faces: $(length(new_faces))")
    end
    # Generate embedded mesh for fractures
    fracture_mesh = Jutul.EmbeddedMeshes.EmbeddedMesh(mesh, fracture_faces)

    return mesh, fracture_mesh, fracture_properties_exp

end

function strike_dip_to_normal(strike::Float64, dip::Float64)
    strike_rad = deg2rad(strike)
    dip_rad = deg2rad(dip)
    normal_x = -sin(dip_rad) * sin(strike_rad)
    normal_y = -sin(dip_rad) * cos(strike_rad)
    normal_z = cos(dip_rad)
    return [normal_x, normal_y, normal_z]
end