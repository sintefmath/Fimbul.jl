function horizontal_fractured_mesh(cell_constraints, depths, num_fractures;
        aperture = 1e-3,
        hz = missing,
        kwargs...
    )

    if length(num_fractures) == 1
        @assert length(depths) == 4
        num_fractures = [0, num_fractures, 0]
    end
    if !ismissing(hz)
        @assert length(hz) == length(depths)-1
    end

    z = [depths[1]]
    hz_all = []
    layer_map = Int64[]
    fracture_map = []
    for i = 1:length(depths)-1
        dz = depths[i+1] - depths[i]
        if num_fractures[i] > 0
            @assert num_fractures[i] > 1
            nf = num_fractures[i]
            dz_mat = (dz - nf*aperture)/(nf-1)
            dz = repeat([aperture, dz_mat], nf)[1:end-1]
            zi = cumsum(dz) .+ depths[i]
            if ismissing(hz)
                hzi = dz_mat/5
            else
                hzi = hz[i]
            end
            hzi = fill(hzi, 2*nf-1)
            @assert isapprox(zi[end], depths[i+1])
            fmap = repeat([true, false], nf)[1:end-1]
        else
            zi = depths[i+1]
            if ismissing(hz)
                hzi = dz/3
            else
                hzi = hz[i]
            end
            fmap = false
        end
        nl = length(zi)
        @assert length(fmap) == nl
        push!(layer_map, fill(i, nl)...)
        push!(fracture_map, fmap...)
        push!(z, zi...)
        push!(hz_all, hzi...)
    end

    mesh, layers, metrics = extruded_mesh(cell_constraints, z; hz = hz_all, kwargs...)
    
    fractures = fracture_map[layers]
    layers = layer_map[layers]
    
    return mesh, layers, fractures, metrics

end