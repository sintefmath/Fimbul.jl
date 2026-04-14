function get_well_neighborship(mesh, coordinates_or_cells, connectivity::Matrix{Int64}, geometry=missing;
    top_node = false, output_directions=false, kwargs...)
    
    # If geometry is not provided, use the mesh geometry
    if ismissing(geometry)
        geometry = tpfv_geometry(mesh)
    end
    num_sections = length(coordinates_or_cells)
    # Allocate arrays
    reservoir_cells = Vector{Vector{Int64}}(undef, 0)
    well_cells = Vector{Vector{Int64}}(undef, 0)
    neighborship = Vector{Matrix{Int64}}(undef, 0)
    sections = Vector{Vector{Int64}}(undef, 0)
    if output_directions
        directions = Vector{Vector{Float64}}(undef, 0)
    end
    # For each well section, find the corresponding reservoir and set up neighborship
    if top_node
        push!(well_cells, [1])
        push!(neighborship, zeros(2,0))
        wc_max = 1
    else
        wc_max = 0
    end
    size(connectivity) == (length(coordinates_or_cells)+top_node, 2) || 
    error("Connectivity matrix must have size (number of well sections, 2) \
    or (number of well sections + 1, 2) if top_node is true.")
    for (sno, x) in enumerate(coordinates_or_cells)
        sno = sno + top_node
        # Find reservoir cells corresponding to the well section
        if x isa Matrix{Float64}
            out = Jutul.find_enclosing_cells(mesh, x;
                geometry=geometry, extra_out=output_directions, kwargs...)
            if output_directions
                rc, extra = out
            else
                rc = out
            end
        else
            rc = x
        end
        push!(reservoir_cells, rc)
        # Create well cells
        wc = collect(1:length(rc)) .+ wc_max
        wc_max += length(wc)
        push!(well_cells, wc)
        # Section neighborship
        n = vcat(wc[1:end-1]', wc[2:end]')
        push!(neighborship, n)
        push!(sections, fill(sno, length(wc)))

        # Direction vectors for each well cell
        if output_directions
            dir = Vector.(extra[:direction].*extra[:lengths])
            append!(directions, dir)
        end

    end

    for (sno, (wc, rc)) in enumerate(zip(well_cells, reservoir_cells))
        # Neighborship from previous section
        n = neighborship[sno]
        from_section = connectivity[sno, 1]
        if from_section > 0
            ix = findfirst(reservoir_cells[from_section] .== rc[1])
            if isnothing(ix)
                @warn "First reservoir cell of section $sno does not match any \
                reservoir cell in the from_section $from_section. Connecting \
                to last cell."
                ix = length(reservoir_cells[from_section])
            end
            wc_from = well_cells[from_section][ix]
            n = hcat([wc_from; wc[1]], n)
        end
        # Neighborship to next section
        to_section = connectivity[sno, 2]
        if to_section > 0
            ix = findfirst(reservoir_cells[to_section] .== rc[end])
            if isnothing(ix)
                @warn "Last reservoir cell of section $sno does not match any \
                reservoir cell in the to_section $to_section. Connecting \
                to first cell."
                ix = 1
            end
            wc_to = well_cells[to_section][ix]
            n = hcat(n, [wc[end]; wc_to])
        end
        neighborship[sno] = n
    end

    reservoir_cells = vcat(reservoir_cells...)
    well_cells = vcat(well_cells...)
    neighborship = hcat(neighborship...)
    sections = vcat(sections...)

    if output_directions
        return reservoir_cells, well_cells, neighborship, sections, directions
    else
        return reservoir_cells, well_cells, neighborship, sections
    end

end
