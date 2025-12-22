function setup_closed_loop_well_coaxial(D::DataDomain, reservoir_cells;
    name = :CL,
    return_reservoir_cell = reservoir_cells[1],
    cell_centers = D[:cell_centroids],
    neighborship = missing,
    pipe_cells_inner = missing,
    pipe_cells_outer = missing,
    grout_cells = missing,
    section = missing,
    well_cell_centers = missing,
    segment_models = missing,
    end_nodes = missing,
    radius_grout = 65e-3,
    radius_pipe_inner = 15e-3,
    radius_pipe_outer = 20e-3,
    wall_thickness_pipe_inner = 3e-3,
    wall_thickness_pipe_outer = 3e-3,
    grouting_thermal_conductivity = 2.3,
    inner_pipe_thermal_conductivity = 0.38,
    outer_pipe_thermal_conductivity = 0.38,
    grouting_heat_capacity = 1460,
    grouting_density = 1500.0,
    kwarg...)


     if ismissing(neighborship)
        if !ismissing(pipe_cells_inner) && !ismissing(pipe_cells_outer) && !ismissing(grout_cells)
            error("""Provide either neighborship, pipe_cells_inner,
                pipe_cells_outer, grout_cells, "well_cell_centers, end_nodes, or
                none of them.""")
        end
        # Set inner/outer pipe and grout cells
        nc_pipe = length(reservoir_cells)*2
        nc_grout = length(reservoir_cells)
        nc_r = length(reservoir_cells)
        pipe_cells_inner = collect(1:nc_r)
        pipe_cells_outer = pipe_cells_inner .+ nc_r
        grout_cells = pipe_cells_outer .+ nc_r
        # Set connectivity
        num_cells = nc_pipe + nc_grout
        pi2pi = vcat(pipe_cells_inner[1:end-1]', pipe_cells_inner[2:end]')
        po2po = vcat(pipe_cells_outer[1:end-1]', pipe_cells_outer[2:end]')
        pi2po_bottom = [pipe_cells_inner[end]', pipe_cells_outer[end]']
        pi2po = vcat(pipe_cells_inner[1:end-1]', pipe_cells_outer[1:end-1]')
        po2g = vcat(pipe_cells_outer', grout_cells')
        neighborship = hcat(pi2pi, po2po, pi2po_bottom, pi2po, po2g)
        # Set centers and end nodes
        well_cell_centers = repeat(cell_centers[:, reservoir_cells], 1, 3)
        end_nodes = [nc_r+1]
        if !ismissing(section)
            warn(["section argument is ignored when neighborship is not provided. ",
                "Sections will be created automatically."])
        end
        section = Vector{Any}(undef, nc_pipe + nc_grout)
        section[pipe_cells_inner] .= [(1, :pipe_inner)]
        section[pipe_cells_outer] .= [(1, :pipe_outer)]
        section[grout_cells] .= [(1, :grout)]
    end

    if ismissing(segment_models)
        # Set up segment flow models
        segment_models = Vector{Any}()
        flow_segs = size(pi2pi, 2) + size(po2po, 2) + size(pi2po_bottom, 2)
        nseg = size(neighborship,2)
        for s in 1:nseg
            if s <= flow_segs
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
    λ_pipe[pipe_cells_inner] .= inner_pipe_thermal_conductivity
    λ_pipe[pipe_cells_outer] .= outer_pipe_thermal_conductivity
    # Set grouting thermal conductivities
    λ_grout = zeros(num_cells)
    λ_grout[grout_cells] .= grouting_thermal_conductivity

    # Common well arguments
    args = (
        type = :closed_loop,
        simple_well = false,
        WI = 0.0
    )

    # Setup supply well
    supply_well = setup_well(D, reservoir_cells;
        name = Symbol(name, "_supply"),
        neighborship = neighborship,
        perforation_cells_well = collect(grout_cells),
        well_cell_centers = well_cell_centers,
        cell_radius = cell_radius,
        cell_radius_inner = cell_radius_inner,
        casing_thickness = wall_thickness,
        casing_thermal_conductivity = λ_pipe,
        grouting_thermal_conductivity = λ_grout,
        segment_models = segment_models,
        end_nodes = end_nodes,
        args..., kwarg...)
    # Augment supply well with closed loop specific data
    augment_closed_loop_domain_coaxial!(supply_well,
        radius_grout,
        grouting_heat_capacity,
        grouting_density,
        section = section
    )
    # Set default thermal indices
    set_default_closed_loop_thermal_indices_coaxial!(supply_well)

    # Setup return well
    return_well = setup_well(D, return_reservoir_cell;
        name = Symbol(name, "_return"),
        WIth = 0.0,
        args...)

    return [supply_well, return_well]

end

function augment_closed_loop_domain_coaxial!(well::DataDomain,
    radius_grout,
    grouting_heat_capacity,
    grouting_density;
    pipe_pipe_thermal_index = missing,
    pipe_grout_thermal_index = missing,
    section = missing
)

    c = Cells()
    p = Perforations()
    f = Faces()

    well[:radius_grout, c] = radius_grout
    well[:material_heat_capacity, c] = grouting_heat_capacity
    well[:material_density, c] = grouting_density

    treat_defaulted(x) = x
    treat_defaulted(::Missing) = NaN
    treat_defaulted(::Nothing) = NaN

    λpp = treat_defaulted(pipe_pipe_thermal_index)
    λpg = treat_defaulted(pipe_grout_thermal_index)
    well[:pipe_pipe_thermal_index, p] = λpp
    well[:pipe_grout_thermal_index, p] = λpg

    if !ismissing(section)
        well[:section, c] = section
    end

end

function set_default_closed_loop_thermal_indices_coaxial!(well::DataDomain)

    cell = Cells()
    face = Faces()
    perf = Perforations()

    num_nodes = well.representation.num_nodes

    hole_volumes = zeros(Float64, num_nodes)
    casing_volumes = zeros(Float64, num_nodes)
    grout_volumes = zeros(Float64, num_nodes)

    nc = Int.(num_nodes/3)
    pipe_cells_inner = Int.(1:nc)
    pipe_cells_outer = Int.(1:nc) .+ nc
    grout_cells = Int.(1:nc) .+ 2*nc
    N = well.representation.neighborship
    rep = physical_representation(well)
    for (pno, gc) in enumerate(rep.perforations.self)

        @assert gc ∈ grout_cells
        seg_pg = findall(N[2, :] .== gc)
        @assert length(seg_pg) == 1
        seg_pg = seg_pg[1]
        pc_out = N[1, seg_pg]
        @assert pc_out ∈ pipe_cells_outer
        seg_pp = findall(N[2, :] .== pc_out)
        @assert length(seg_pp) <= 2
        seg_pp = seg_pp[end]
        pc_in = N[1, seg_pp]
        @assert pc_in ∈ pipe_cells_inner

        r_grout = well[:radius_grout, cell][gc]
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
        hole_volumes[gc] = 0.0
        casing_volumes[pc_in] = vol_iw
        casing_volumes[pc_out] = vol_ow
        casing_volumes[gc] = 0.0
        grout_volumes[pc_in] = 0.0
        grout_volumes[pc_out] = 0.0
        grout_volumes[gc] = vol_g

        λg = well[:grouting_thermal_conductivity, cell][gc]
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

        if isnan(well[:thermal_well_index, perf][pno])
            well[:thermal_well_index, perf][pno] = λgr
        end

        if well[:material_thermal_conductivity, face][seg_pg] == 0.0
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
    x = log(sqrt((d_g^2 + d_opo^2)/(2*d_opo^2)))/log(d_g/d_opo)
    R_g = log(d_g/d_opo)/(2*π*λ_g)
    R_cg = x*R_g

    Rpp = R_ai + R_ci
    Rpg = R_ao + R_co + R_cg
    Rgr = (1-x)*R_g

    return Rpp, Rpg, Rgr

end