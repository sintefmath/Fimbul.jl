"""
    setup_closed_loop_supply(type::ClosedLoopCoaxialWell, D, reservoir_cells; kwarg...)

Supply part of a coaxial closed-loop well. By default, the inner pipe goes
down through `reservoir_cells` and the outer annulus back up again, embedded
in grout (or another fill material). Nodes are tagged `:pipe_inner`,
`:pipe_outer` and `:grout`. If `grouting_thickness` is zero, there are no grout
nodes, and the outer pipe exchanges heat directly with the reservoir.

# Keyword arguments
- `grouting_thickness = 30e-3`: Thickness of the borehole fill outside the
  outer pipe [m]. Can be zero.
- `radius_pipe_inner = 50e-3`: Outer radius of the inner pipe [m].
- `radius_pipe_outer = sqrt(2)*radius_pipe_inner`: Outer radius of the outer
  pipe [m].
- `wall_thickness_pipe_inner = 3e-3`, `wall_thickness_pipe_outer = 3e-3`:
  Pipe wall thicknesses [m].
- `pipe_thermal_conductivity_inner = 0.38`,
  `pipe_thermal_conductivity_outer = 0.38`: Pipe wall thermal conductivities
  [W/(m K)].
- `fill = :grout`: Borehole fill material, see
  [`Fimbul.borehole_fill_properties`](@ref).
- `fill_thermal_conductivity`, `fill_density`, `fill_heat_capacity`: Override
  the thermal conductivity [W/(m K)], density [kg/m³] and heat capacity
  [J/(kg K)] of `fill`.
- `dir = :z`: Direction of the well in each reservoir cell.
- `neighborship`, `pipe_cells_inner`, `pipe_cells_outer`, `grout_cells`,
  `well_cell_centers`, `end_nodes`: General well topology. Must be given
  together (`grout_cells` must be omitted or empty if `grouting_thickness` is
  zero).
- `segment_models`, `tag`: Segment models and node tags for a general topology.
"""
function setup_closed_loop_supply(type::ClosedLoopCoaxialWell, D::DataDomain, reservoir_cells;
    name = :CL_inner,
    cell_centers = D[:cell_centroids],
    neighborship = missing,
    pipe_cells_inner = missing,
    pipe_cells_outer = missing,
    grout_cells = missing,
    tag = missing,
    well_cell_centers = missing,
    dir = :z,
    segment_models = missing,
    end_nodes = missing,
    grouting_thickness = missing,
    radius_pipe_inner = 50e-3,
    radius_pipe_outer = sqrt(2)*radius_pipe_inner,
    wall_thickness_pipe_inner = 3e-3,
    wall_thickness_pipe_outer = 3e-3,
    pipe_thermal_conductivity_inner = 0.38,
    pipe_thermal_conductivity_outer = 0.38,
    fill = :grout,
    fill_thermal_conductivity = missing,
    fill_density = missing,
    fill_heat_capacity = missing,
    # Deprecated keyword arguments
    radius_grout = missing,
    grouting_thermal_conductivity = missing,
    grouting_heat_capacity = missing,
    grouting_density = missing,
    kwarg...)

    if !ismissing(radius_grout)
        Base.depwarn("Keyword argument `radius_grout` is deprecated for \
            coaxial closed-loop wells, use `grouting_thickness` instead.",
            :setup_closed_loop_well)
        ismissing(grouting_thickness) || error("Cannot give both \
            `radius_grout` and `grouting_thickness`.")
        grouting_thickness = radius_grout - radius_pipe_outer
    end
    if ismissing(grouting_thickness)
        grouting_thickness = 30e-3
    end
    grouting_thickness >= 0.0 || error("Grouting thickness must be non-negative.")
    has_grout = grouting_thickness > 0.0
    fill_props = closed_loop_fill_properties(fill, fill_thermal_conductivity,
        fill_density, fill_heat_capacity;
        grouting_thermal_conductivity = grouting_thermal_conductivity,
        grouting_density = grouting_density,
        grouting_heat_capacity = grouting_heat_capacity)

    if ismissing(neighborship)
        check_closed_loop_topology_args(pipe_cells_inner, pipe_cells_outer,
            grout_cells, well_cell_centers, end_nodes)
        # Set inner/outer pipe and grout cells
        nc_r = length(reservoir_cells)
        pipe_cells_inner = collect(1:nc_r)
        pipe_cells_outer = pipe_cells_inner .+ nc_r
        grout_cells = has_grout ? pipe_cells_outer .+ nc_r : Int[]
        num_cells = 2*nc_r + length(grout_cells)
        # Set connectivity
        pi2pi = vcat(pipe_cells_inner[1:end-1]', pipe_cells_inner[2:end]')
        po2po = vcat(pipe_cells_outer[1:end-1]', pipe_cells_outer[2:end]')
        pi2po_bottom = [pipe_cells_inner[end], pipe_cells_outer[end]]
        pi2po = vcat(pipe_cells_inner[1:end-1]', pipe_cells_outer[1:end-1]')
        neighborship = hcat(pi2pi, po2po, pi2po_bottom, pi2po)
        if has_grout
            po2g = vcat(pipe_cells_outer', grout_cells')
            neighborship = hcat(neighborship, po2g)
        end
        # Set centers and end nodes
        well_cell_centers = repeat(cell_centers[:, reservoir_cells], 1,
            div(num_cells, nc_r))
        end_nodes = [nc_r+1]
        if !ismissing(tag)
            @warn "tag argument is ignored when neighborship is not \
                provided. Tags will be created automatically."
        end
        # Set tag for identifying cell groups (e.g., inner/outer pipes, grout)
        tag = Vector{Symbol}(undef, num_cells)
        tag[pipe_cells_inner] .= :pipe_inner
        tag[pipe_cells_outer] .= :pipe_outer
        tag[grout_cells] .= :grout
    else
        if !has_grout && ismissing(grout_cells)
            grout_cells = Int[]
        end
        if ismissing(pipe_cells_inner) || ismissing(pipe_cells_outer) ||
            ismissing(grout_cells) || ismissing(well_cell_centers) ||
            ismissing(end_nodes)
            error("If neighborship is provided, pipe_cells_inner, \
                pipe_cells_outer, grout_cells, well_cell_centers, and \
                end_nodes must also be provided.")
        end
        has_grout || isempty(grout_cells) || error("grout_cells must be \
            empty when grouting_thickness is zero.")
        num_cells = length(pipe_cells_inner) + length(pipe_cells_outer) + length(grout_cells)
    end

    if ismissing(segment_models)
        # Set up segment flow models
        segment_models = Vector{Any}()
        pipe_connection = [pipe_cells_inner[end], pipe_cells_outer[end]]
        # Check that segment is either from inner to inner, outer to outer, or inner to outer at bottom
        is_flow(seg) =
            all([n ∈ pipe_cells_inner for n in seg]) ||
            all([n ∈ pipe_cells_outer for n in seg]) ||
            all([n ∈ pipe_connection for n in seg])
        for seg in eachcol(neighborship)
            if is_flow(seg)
                # Pipe segments use standard wellbore friction model
                seg_model = SegmentWellBoreFrictionHB()
            else
                seg_model = JutulDarcy.ClosedSegment()
            end
            push!(segment_models, seg_model)
        end
    end

    # Set cell radii
    cell_radius = zeros(num_cells)
    cell_radius_inner = zeros(num_cells)
    cell_radius[pipe_cells_inner] .= radius_pipe_inner - wall_thickness_pipe_inner
    cell_radius[pipe_cells_outer] .= radius_pipe_outer - wall_thickness_pipe_outer
    cell_radius_inner[pipe_cells_outer] .= radius_pipe_inner
    # Set casing thicknesses
    wall_thickness = zeros(num_cells)
    wall_thickness[pipe_cells_inner] .= wall_thickness_pipe_inner
    wall_thickness[pipe_cells_outer] .= wall_thickness_pipe_outer
    # Set pipe wall thermal conductivities
    λ_pipe = zeros(num_cells)
    λ_pipe[pipe_cells_inner] .= pipe_thermal_conductivity_inner
    λ_pipe[pipe_cells_outer] .= pipe_thermal_conductivity_outer
    # Set grouting thermal conductivities
    λ_grout = zeros(num_cells)
    λ_grout[grout_cells] .= fill_props.thermal_conductivity

    # Common well arguments
    args = (
        type = type,
        simple_well = false,
    )
    if has_grout
        perforation_cells_well = collect(grout_cells)
    else
        # Without grout, the outer pipe is perforated. It must not exchange
        # mass with the reservoir.
        perforation_cells_well = collect(pipe_cells_outer)
        args = merge(args, (WI = 0.0,))
    end

    # Setup supply well
    supply_well = setup_well(D, reservoir_cells;
        name = name,
        neighborship = neighborship,
        perforation_cells_well = perforation_cells_well,
        well_cell_centers = well_cell_centers,
        cell_radius = cell_radius,
        cell_radius_inner = cell_radius_inner,
        casing_thickness = wall_thickness,
        casing_thermal_conductivity = λ_pipe,
        grouting_thermal_conductivity = λ_grout,
        segment_models = segment_models,
        end_nodes = end_nodes,
        dir = dir,
        args..., kwarg...)
    # Augment supply well with closed loop specific data
    augment_closed_loop_domain_coaxial!(supply_well,
        radius_pipe_outer .+ grouting_thickness,
        fill_props.heat_capacity,
        fill_props.density,
        tag = tag
    )
    # Set default thermal indices
    set_default_closed_loop_thermal_indices_coaxial!(supply_well,
        pipe_cells_inner, pipe_cells_outer, grout_cells)

    return supply_well

end

function augment_closed_loop_domain_coaxial!(well::DataDomain,
    radius_grout,
    grouting_heat_capacity,
    grouting_density;
    pipe_pipe_thermal_index = missing,
    pipe_grout_thermal_index = missing,
    tag = missing
)

    c = Cells()
    p = Perforations()
    f = Faces()

    well[:radius_grout, c] = radius_grout
    well[:material_heat_capacity, c] = grouting_heat_capacity
    well[:material_density, c] = grouting_density

    λpp = treat_defaulted(pipe_pipe_thermal_index)
    λpg = treat_defaulted(pipe_grout_thermal_index)
    well[:pipe_pipe_thermal_index, p] = λpp
    well[:pipe_grout_thermal_index, p] = λpg

    if !ismissing(tag)
        well[:tag, c] = tag
    end

end

function set_default_closed_loop_thermal_indices_coaxial!(well::DataDomain,
    pipe_cells_inner, pipe_cells_outer, grout_cells)

    cell = Cells()
    face = Faces()
    perf = Perforations()

    num_nodes = well.representation.num_nodes

    hole_volumes = zeros(Float64, num_nodes)
    casing_volumes = zeros(Float64, num_nodes)
    grout_volumes = zeros(Float64, num_nodes)

    N = well.representation.neighborship
    rep = physical_representation(well)
    for (pno, perforated_cell) in enumerate(rep.perforations.self)

        if perforated_cell ∈ grout_cells
            # Grout cell connected to the outer pipe
            gc = perforated_cell
            seg_pg = findall(N[2, :] .== gc)
            @assert length(seg_pg) == 1
            seg_pg = seg_pg[1]
            pc_out = N[1, seg_pg]
        else
            # No grout: the outer pipe is perforated directly
            gc = nothing
            seg_pg = nothing
            pc_out = perforated_cell
        end
        @assert pc_out ∈ pipe_cells_outer
        seg_pp = findall(j -> N[2, j] == pc_out && N[1, j] ∈ pipe_cells_inner,
            axes(N, 2))
        @assert 1 <= length(seg_pp) <= 2
        seg_pp = seg_pp[end]
        pc_in = N[1, seg_pp]

        r_grout = well[:radius_grout, cell][pc_out]
        r_inner_pipe = well[:radius, cell][pc_in]
        wall_thickness_inner = well[:casing_thickness, cell][pc_in]
        r_outer_pipe = well[:radius, cell][pc_out]
        wall_thickness_outer = well[:casing_thickness, cell][pc_out]
        L = well[:cell_length, cell][pc_in]
        vol_ip, vol_iw, vol_op, vol_ow, vol_g, L = closed_loop_volume_coaxial(
            L, r_grout,
            r_inner_pipe + wall_thickness_inner, wall_thickness_inner,
            r_outer_pipe + wall_thickness_outer, wall_thickness_outer
        )

        hole_volumes[pc_in] = vol_ip
        hole_volumes[pc_out] = vol_op
        casing_volumes[pc_in] = vol_iw
        casing_volumes[pc_out] = vol_ow
        grout_volumes[pc_in] = 0.0
        grout_volumes[pc_out] = 0.0
        if !isnothing(gc)
            # Small fluid volume in the grout to avoid singular mass balance
            # (consistent with the U1 closed-loop well)
            hole_volumes[gc] = 1e-3*vol_g
            casing_volumes[gc] = 0.0
            grout_volumes[gc] = vol_g
        end

        λg = isnothing(gc) ? 0.0 : well[:grouting_thermal_conductivity, cell][gc]
        λpi = well[:casing_thermal_conductivity, cell][pc_in]
        λpo = well[:casing_thermal_conductivity, cell][pc_out]
        Rpp, Rpg, Rgr = closed_loop_thermal_resistance_coaxial(
            r_grout,
            r_inner_pipe + wall_thickness_inner, wall_thickness_inner,
            r_outer_pipe + wall_thickness_outer, wall_thickness_outer,
            λg, λpi, λpo)
        λpp = L/Rpp
        λpg = L/Rpg
        λgr = L/Rgr

        # The thermal well index only accounts for the resistance inside the
        # borehole, and contains no term for conduction from the borehole wall
        # to the center of the reservoir cell. The cell temperature is thus
        # assumed to approximate the borehole wall temperature, as in the
        # thermal resistance and capacity model of Diersch et al. (2011),
        # which requires a reservoir grid that is refined around the well.
        # Without grout, the outer pipe wall is the only resistance.
        WIth = isnothing(gc) ? L/(Rpg + Rgr) : λgr
        isfinite(WIth) || error("Thermal well index is infinite. Without \
            grout, the outer pipe wall must have a positive thickness.")
        if isnan(well[:thermal_well_index, perf][pno])
            well[:thermal_well_index, perf][pno] = WIth
        end

        if !isnothing(seg_pg) && well[:material_thermal_conductivity, face][seg_pg] == 0.0
            well[:material_thermal_conductivity, face][seg_pg] = λpg
        end

        if well[:material_thermal_conductivity, face][seg_pp] == 0.0
            well[:material_thermal_conductivity, face][seg_pp] = λpp
        end
    end

    well[:volume_override_hole, cell] = hole_volumes
    well[:volume_override_casing, cell] = casing_volumes
    well[:volume_override_grouting, cell] = grout_volumes
end

function closed_loop_volume_coaxial(
    length,
    radius_grout,
    radius_inner_pipe, wall_thickness_inner_pipe,
    radius_outer_pipe, wall_thickness_outer_pipe
    )

    # Compute pipe and grout volume
    r_g, r_ipi, r_ipo, r_opi, r_opo, L =
        radius_grout,
        radius_inner_pipe - wall_thickness_inner_pipe,
        radius_inner_pipe,
        radius_outer_pipe - wall_thickness_outer_pipe,
        radius_outer_pipe,
        length

    vol_hole_inner = π*r_ipi^2*L
    vol_wall_inner = π*(r_ipo^2 - r_ipi^2)*L
    vol_hole_outer = π*(r_opi^2 - r_ipo^2)*L
    vol_wall_outer = π*(r_opo^2 - r_opi^2)*L
    vol_grout = π*(r_g^2 - r_opo^2)*L

    return vol_hole_inner, vol_wall_inner, vol_hole_outer, vol_wall_outer, vol_grout, L

end

function closed_loop_thermal_resistance_coaxial(
    radius_grout,
    radius_inner_pipe, wall_thickness_inner_pipe,
    radius_outer_pipe, wall_thickness_outer_pipe,
    thermal_conductivity_grout,
    thermal_conductivity_inner_pipe,
    thermal_conductivity_outer_pipe
    )

    r_g, r_ipi, r_ipo, r_opi, r_opo =
        radius_grout,
        radius_inner_pipe - wall_thickness_inner_pipe,
        radius_inner_pipe,
        radius_outer_pipe - wall_thickness_outer_pipe,
        radius_outer_pipe

    λ_g, λ_pi, λ_po =
    thermal_conductivity_grout,
    thermal_conductivity_inner_pipe,
    thermal_conductivity_outer_pipe

    R_ai = 0.0 # TODO: Implement this
    R_ao = 0.0

    R_ci = log(r_ipo/r_ipi)/(2*π*λ_pi)
    R_co = log(r_opo/r_opi)/(2*π*λ_po)

    d_g, d_opo = 2*r_g, 2*r_opo
    if d_g > d_opo
        x = log(sqrt((d_g^2 + d_opo^2)/(2*d_opo^2)))/log(d_g/d_opo)
        R_g = log(d_g/d_opo)/(2*π*λ_g)
    else
        # No grout
        x = 0.0
        R_g = 0.0
    end
    R_cg = x*R_g

    Rpp = R_ai + R_ci
    Rpg = R_ao + R_co + R_cg
    Rgr = (1-x)*R_g

    return Rpp, Rpg, Rgr

end