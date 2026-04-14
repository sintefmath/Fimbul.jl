function add_fractures(mesh::UnstructuredMesh, centers::Vector, normals::Vector;
    tags=collect(1:length(centers)))

    length(centers) == length(normals) || error("Number of fracture centers and normals must match. Got $(length(centers)) centers and $(length(normals)) normals.")
    length(tags) == length(centers) || error("Not enough fracture tags provided. Expected $(length(centers)), got $(length(tags)).")
    
    fracture_faces = Int[]
    fracture_tag = []
    for (i, (normal, center)) in enumerate(zip(normals, centers))
        is_frac = falses(number_of_faces(mesh))
        is_frac[fracture_faces] .= true
        plane = PlaneCut(center, normal)
        mesh, info = cut_mesh(mesh, plane; min_cut_fraction = 0.0, extra_out=true)
        # Update fracture face vector
        face_index = filter(x->x>0, info["face_index"])
        is_frac = is_frac[face_index]
        new_faces = findall(info["face_index"] .== 0)
        # Utils for mapping old to new faces and cells
        old_to_new_faces = face_index[is_frac]
        face_to_cell = zeros(Int, number_of_faces(mesh))
        face_to_cell[fracture_faces] .= 1:length(fracture_faces)
        old_to_new_cells = face_to_cell[old_to_new_faces]
        # Update existing fracture tags
        fracture_tag = fracture_tag[old_to_new_cells]
        # Add new fracture tags
        append!(fracture_tag, fill(tags[i], length(new_faces)))
        # Add new fracture faces to fracture face vector
        fracture_faces = findall(is_frac)
        append!(fracture_faces, new_faces)
    end
    # Generate embedded mesh for fractures
    fracture_mesh = Jutul.EmbeddedMeshes.EmbeddedMesh(mesh, fracture_faces)

    return mesh, fracture_mesh, fracture_tag

end

function strike_dip_to_normal(strike::Float64, dip::Float64)
    strike_rad = deg2rad(strike)
    dip_rad = deg2rad(dip)
    normal_x = -sin(dip_rad) * sin(strike_rad)
    normal_y = -sin(dip_rad) * cos(strike_rad)
    normal_z = cos(dip_rad)
    return [normal_x, normal_y, normal_z]
end
