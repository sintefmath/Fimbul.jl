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
        push!(sections, [1])
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

    tn = Int(top_node)
    for sno in eachindex(well_cells)
        wc = well_cells[sno]
        # Correctly map well section index to reservoir_cells index,
        # accounting for the extra top-node entry in well_cells.
        rc_idx = sno - tn
        rc = (rc_idx >= 1 && rc_idx <= length(reservoir_cells)) ? reservoir_cells[rc_idx] : Int[]
        # Neighborship from previous section
        n = neighborship[sno]
        from_section = connectivity[sno, 1]
        if from_section > 0
            from_rc_idx = from_section - tn
            if from_rc_idx >= 1
                # Normal section: find matching reservoir cell
                ix = isempty(rc) ? nothing : findfirst(reservoir_cells[from_rc_idx] .== rc[1])
                if isnothing(ix)
                    @warn "First reservoir cell of section $sno does not match any \
                    reservoir cell in the from_section $from_section. Connecting \
                    to last cell."
                    ix = length(reservoir_cells[from_rc_idx])
                end
            else
                # from_section is the top node (no reservoir cells): connect to its only well cell
                ix = 1
            end
            wc_from = well_cells[from_section][ix]
            n = hcat([wc_from; wc[1]], n)
        end
        # Neighborship to next section
        to_section = connectivity[sno, 2]
        if to_section > 0
            to_rc_idx = to_section - tn
            if to_rc_idx >= 1
                ix = isempty(rc) ? nothing : findfirst(reservoir_cells[to_rc_idx] .== rc[end])
                if isnothing(ix)
                    @warn "Last reservoir cell of section $sno does not match any \
                    reservoir cell in the to_section $to_section. Connecting \
                    to first cell."
                    ix = 1
                end
            else
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

"""
    adjust_well_indices!(well, well_name, fractures=false)

Adjust the drainage radius of well perforations where the equivalent Peaceman
radius would exceed the cell dimensions (i.e. `re ≤ 3*r_perf`). In such cases
the drainage radius is set to `r_perf * (1 + 1e-3)` to satisfy the assertion
in `compute_peaceman_index`.

Set `fractures=true` to operate on the fracture well indices instead of the
matrix well indices.
"""
function adjust_well_indices!(well, well_name, fractures=false)
    if !fractures
        cd = :cell_dims
        pr = :perforation_radius
        dr = :drainage_radius
        pd = :perforation_direction
        e  = Perforations()
    else
        cd = :cell_dims_frac
        pr = :perforation_radius_frac
        dr = :drainage_radius_frac
        pd = :perforation_direction_frac
        e  = JutulDarcy.FracturePerforations()
    end
    Δ      = well[cd, e]
    radius = well[pr, e]
    dir    = well[pd, e]
    num_violations = 0
    for (k, (Δk, rk, dk)) in enumerate(zip(Δ, radius, dir))
        if dk == :x
            i, j = 1, 3
        elseif dk == :y
            i, j = 2, 3
        elseif dk == :z
            i, j = 1, 2
        else
            error("Unknown perforation direction $dk for well $well_name")
        end
        re = 0.14 * sqrt(Δk[i]^2 + Δk[j]^2)
        if re <= 3 * rk
            well[dr, e][k] = rk * (1.0 + 1e-3)
            num_violations += 1
        end
    end
    if num_violations > 0
        kind = fractures ? "fracture " : ""
        @warn "Adjusted $num_violations $(kind)drainage radii for well \
        $well_name due to cell dimensions smaller than perforation radius."
    end
    return well
end
