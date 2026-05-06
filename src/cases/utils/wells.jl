function get_well_neighborship(mesh, coordinates_or_cells, connectivity::Matrix{Int64}, geometry=missing;
    type = :standard, top_node = false, output_directions=false, kwargs...)

    if type == :closed_loop_u1
        return _get_cl_u1_neighborship(mesh, coordinates_or_cells, geometry; kwargs...)
    elseif type == :closed_loop_coaxial
        return _get_cl_coaxial_neighborship(mesh, coordinates_or_cells, geometry; kwargs...)
    elseif type != :standard
        error("Unknown well type: $type. Valid types are :standard, :closed_loop_u1, :closed_loop_coaxial")
    end

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

"""
    _get_cl_u1_neighborship(mesh, coordinates_or_cells, geometry; kwargs...)

Internal helper for `get_well_neighborship` with `type = :closed_loop_u1`.

Sets up the neighborship for one or more U-tube (U1) borehole heat exchangers.
Each entry in `coordinates_or_cells` represents one borehole: either a matrix of
3-D coordinates (one row per depth level) or a vector of pre-computed reservoir
cell indices (going downward).

The resulting well has two manifold nodes — node 1 (supply inlet) and node 2
(return outlet) — followed by pipe and grout cell groups for each borehole:

  * **pipe cells** (2n per borehole): supply leg (n cells going down) followed
    by return leg (n cells going back up).
  * **grout cells** (2n per borehole): thermally coupled to the pipe cells,
    with supply-side grout connected serially and return-side grout connected
    serially in reverse so that the grout path crosses at the bottom.

Returns `(reservoir_cells, well_cells, neighborship, sections, pipe_cells, grout_cells)`.
"""
function _get_cl_u1_neighborship(mesh, coordinates_or_cells, geometry; kwargs...)
    if ismissing(geometry)
        geometry = tpfv_geometry(mesh)
    end

    # --- find reservoir cells for each borehole ---
    reservoir_cells_per_bh = Vector{Vector{Int64}}(undef, 0)
    for x in coordinates_or_cells
        if x isa Matrix{Float64}
            rc = Jutul.find_enclosing_cells(mesh, x; geometry=geometry, kwargs...)
        else
            rc = collect(Int64, x)
        end
        push!(reservoir_cells_per_bh, rc)
    end

    # Manifold nodes: 1 = supply top, 2 = return top
    offset = 2
    all_pipe_cells  = Vector{Int64}(undef, 0)
    all_grout_cells = Vector{Int64}(undef, 0)
    all_neighborship = Vector{Matrix{Int64}}(undef, 0)
    pipe_sections   = Vector{Int64}(undef, 0)
    grout_sections  = Vector{Int64}(undef, 0)
    all_reservoir_cells = Vector{Int64}(undef, 0)

    for (bno, rc) in enumerate(reservoir_cells_per_bh)
        n = length(rc)
        pipe_cells  = collect(1:2n) .+ offset
        grout_cells = collect(1:2n) .+ (offset + 2n)

        # Reservoir cells: supply leg (down) + return leg (up)
        append!(all_reservoir_cells, rc, reverse(rc))

        # Neighborship connections
        s2p  = reshape([1, pipe_cells[1]], 2, 1)                            # supply → first pipe
        e2p  = reshape([pipe_cells[end], 2], 2, 1)                          # last pipe → return
        p2p  = vcat(pipe_cells[1:end-1]', pipe_cells[2:end]')               # pipe serial
        p2g  = vcat(pipe_cells', grout_cells')                              # pipe → grout (radial)
        ncg  = div(length(grout_cells), 2)
        g2g  = vcat(grout_cells[1:ncg]', grout_cells[end:-1:ncg+1]')       # grout crossing at bottom
        push!(all_neighborship, hcat(s2p, e2p, p2p, p2g, g2g))

        append!(all_pipe_cells,  pipe_cells)
        append!(all_grout_cells, grout_cells)
        append!(pipe_sections,   fill(bno, 2n))
        append!(grout_sections,  fill(bno, 2n))

        offset += 4n
    end

    well_cells   = vcat(all_pipe_cells, all_grout_cells)
    neighborship = hcat(all_neighborship...)
    sections     = vcat(pipe_sections, grout_sections)

    return all_reservoir_cells, well_cells, neighborship, sections, all_pipe_cells, all_grout_cells
end

"""
    _get_cl_coaxial_neighborship(mesh, coordinates_or_cells, geometry; kwargs...)

Internal helper for `get_well_neighborship` with `type = :closed_loop_coaxial`.

Sets up the neighborship for one or more coaxial borehole heat exchangers.
Each entry in `coordinates_or_cells` represents one coaxial unit: either a
matrix of 3-D coordinates (one row per depth level) or a vector of pre-computed
reservoir cell indices (going downward).

The resulting well has two manifold nodes — node 1 (supply inlet) and node 2
(return outlet) — followed by three cell groups for each unit:

  * **inner pipe cells** (n per unit): fluid flows downward.
  * **outer pipe cells** (n per unit): fluid returns upward (thermally coupled
    to the inner pipe and to the grout annulus).
  * **grout cells** (n per unit): thermally coupled to the outer pipe and to
    the reservoir.

Connections:
  - supply (node 1) → inner pipe top
  - inner pipe: serial downward
  - inner pipe bottom → outer pipe bottom  (hydraulic U-turn)
  - outer pipe: serial upward (array order bottom→top, flow top←bottom)
  - outer pipe top → return (node 2)
  - inner ↔ outer coupling at each depth (excluding the bottom node)
  - outer → grout coupling at each depth

Returns `(reservoir_cells, well_cells, neighborship, sections, pipe_cells_inner, pipe_cells_outer, grout_cells)`.
"""
function _get_cl_coaxial_neighborship(mesh, coordinates_or_cells, geometry; kwargs...)
    if ismissing(geometry)
        geometry = tpfv_geometry(mesh)
    end

    # --- find reservoir cells for each coaxial unit ---
    reservoir_cells_per_unit = Vector{Vector{Int64}}(undef, 0)
    for x in coordinates_or_cells
        if x isa Matrix{Float64}
            rc = Jutul.find_enclosing_cells(mesh, x; geometry=geometry, kwargs...)
        else
            rc = collect(Int64, x)
        end
        push!(reservoir_cells_per_unit, rc)
    end

    # Manifold nodes: 1 = supply top, 2 = return top
    offset = 2
    all_pipe_cells_inner = Vector{Int64}(undef, 0)
    all_pipe_cells_outer = Vector{Int64}(undef, 0)
    all_grout_cells      = Vector{Int64}(undef, 0)
    all_neighborship     = Vector{Matrix{Int64}}(undef, 0)
    inner_sections       = Vector{Int64}(undef, 0)
    outer_sections       = Vector{Int64}(undef, 0)
    grout_sections       = Vector{Int64}(undef, 0)
    all_reservoir_cells  = Vector{Int64}(undef, 0)

    for (uno, rc) in enumerate(reservoir_cells_per_unit)
        n = length(rc)
        pipe_cells_inner = collect(1:n) .+ offset
        pipe_cells_outer = collect(1:n) .+ (offset + n)
        grout_cells      = collect(1:n) .+ (offset + 2n)

        # Reservoir cells: one entry per depth level (grout perforations)
        append!(all_reservoir_cells, rc)

        # Neighborship connections
        s2p          = reshape([1, pipe_cells_inner[1]], 2, 1)                              # supply → inner top
        e2p          = reshape([pipe_cells_outer[1], 2], 2, 1)                              # outer top → return
        pi2pi        = vcat(pipe_cells_inner[1:end-1]', pipe_cells_inner[2:end]')           # inner serial (down)
        po2po        = vcat(pipe_cells_outer[1:end-1]', pipe_cells_outer[2:end]')           # outer serial
        pi2po_bottom = reshape([pipe_cells_inner[end], pipe_cells_outer[end]], 2, 1)        # bottom U-turn
        pi2po        = vcat(pipe_cells_inner[1:end-1]', pipe_cells_outer[1:end-1]')         # inner↔outer coupling
        po2g         = vcat(pipe_cells_outer', grout_cells')                                # outer → grout
        push!(all_neighborship, hcat(s2p, e2p, pi2pi, po2po, pi2po_bottom, pi2po, po2g))

        append!(all_pipe_cells_inner, pipe_cells_inner)
        append!(all_pipe_cells_outer, pipe_cells_outer)
        append!(all_grout_cells,      grout_cells)
        append!(inner_sections,       fill(uno, n))
        append!(outer_sections,       fill(uno, n))
        append!(grout_sections,       fill(uno, n))

        offset += 3n
    end

    well_cells   = vcat(all_pipe_cells_inner, all_pipe_cells_outer, all_grout_cells)
    neighborship = hcat(all_neighborship...)
    sections     = vcat(inner_sections, outer_sections, grout_sections)

    return all_reservoir_cells, well_cells, neighborship, sections,
           all_pipe_cells_inner, all_pipe_cells_outer, all_grout_cells
end
