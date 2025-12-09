function setup_closed_loop_well_u1(D::DataDomain, reservoir_cells;
    name = :CL,
    return_reservoir_cell = reservoir_cells[1],
    cell_centers = D[:cell_centroids],
    neighborship = missing,
    pipe_cells = missing,
    grout_cells = missing,
    well_cell_centers = missing,
    segment_models = missing,
    end_nodes = missing,
    radius_grout = 65e-3,
    radius_pipe = 15e-3,
    wall_thickness = 3e-3,
    pipe_spacing = 60e-3,
    grouting_thermal_conductivity = 2.3,
    pipe_thermal_conductivity = 0.38,
    heat_capacity_grout = 1460,
    density_grout = 1500.0,
    kwarg...)

    if ismissing(neighborship)
        if !ismissing(pipe_cells) && !ismissing(grout_cells)
            error(["Provide either neighborship, pipe_cells, grout_cells, ",
                "well_cell_centers, end_nodes, or none of them."])
        end
        # Create grout and pipe cells
        reservoir_cells = vcat(reservoir_cells, reverse(reservoir_cells))
        nc_pipe = length(reservoir_cells)
        nc_grout = length(reservoir_cells)
        nc_mid = div(nc_pipe, 2)
        down_ix = Int.(1:nc_mid)
        up_ix = Int.(1:nc_mid) .+ nc_mid
        pipe_cells = collect(1:nc_pipe)
        grout_cells = collect(1:nc_grout) .+ nc_pipe
        # Set up connectivity
        pipe_to_pipe = vcat(pipe_cells[1:end-1]', pipe_cells[2:end]')
        pipe_to_grout = vcat(pipe_cells', grout_cells')
        grout_to_grout = vcat(grout_cells[down_ix]', reverse(grout_cells[up_ix]'))
        neighborship = hcat(pipe_to_pipe, pipe_to_grout, grout_to_grout)
        # Set up well cell centers and end nodes
        well_cell_centers = repeat(cell_centers[:, reservoir_cells], 1, 2)
        end_nodes = [pipe_cells[end]]
    else
        if ismissing(pipe_cells) || ismissing(grout_cells) || 
            ismissing(well_cell_centers) || ismissing(end_nodes) 
            error(["If neighborship is provided, pipe_cells, grout_cells, ",
                "well_cell_centers, and end_nodes must also be provided."])
        end
    end

    # Set up segment flow models
    if ismissing(segment_models)
        segment_models = Vector{Any}()
        for seg in eachcol(neighborship)
            if all([c ∈ pipe_cells for c in seg])
                # Pipe segments use standard wellbore friction model
                seg_model = SegmentWellBoreFrictionHB()
            else
                seg_model = JutulDarcy.ClosedSegment()
            end
            push!(segment_models, seg_model)
        end
    end

    # Setup pipe radius
    nc = maximum(neighborship)
    if radius_pipe isa Number
        radius_pipe = fill(radius_pipe, nc)
        radius_pipe[grout_cells] .= 0.0
    end
    # Setup wall thickness
    if wall_thickness isa Number
        wall_thickness = fill(wall_thickness, nc)
        wall_thickness[grout_cells] .= 0.0
    end
    # Setup pipe wall thermal conductivities
    if pipe_thermal_conductivity isa Number
        pipe_thermal_conductivity = fill(pipe_thermal_conductivity, nc)
        pipe_thermal_conductivity[grout_cells] .= 0.0
    end
    # Setup grouting thermal conductivities
    if grouting_thermal_conductivity isa Number
        grouting_thermal_conductivity = fill(grouting_thermal_conductivity, nc)
        grouting_thermal_conductivity[pipe_cells] .= 0.0
    end

    # Set up common well properties
    args = (
        type = :closed_loop,
        simple_well = false,
        WI = 0.0
    )
    # Set up supply well
    supply_well = setup_well(D, reservoir_cells;
        name = Symbol(name, "_supply"),
        neighborship = neighborship,
        perforation_cells_well = grout_cells,
        well_cell_centers = well_cell_centers,
        cell_radius = radius_pipe - wall_thickness,
        casing_thickness = wall_thickness,
        casing_thermal_conductivity = pipe_thermal_conductivity,
        grouting_thermal_conductivity = grouting_thermal_conductivity,
        segment_models = segment_models,
        end_nodes = end_nodes,
        args..., kwarg...)
    # Augment supply well with closed_loop-specific properties
    augment_closed_loop_domain_u1!(supply_well,
        radius_grout,
        pipe_spacing,
        heat_capacity_grout,
        density_grout
    )
    # Set default thermal indices
    set_default_closed_loop_thermal_indices_u1!(
    supply_well, pipe_cells, grout_cells)
    # Set up return well
    return_well = setup_well(D, return_reservoir_cell;
        name = Symbol(name, "_return"),
        WIth = 0.0,
        args...)

    return [supply_well, return_well]

end

function augment_closed_loop_domain_u1!(well::DataDomain,
    radius_grout,
    pipe_spacing,
    heat_capacity_grout,
    density_grout;
    pipe_grout_thermal_index = missing,
    grout_grout_thermal_index = missing
)

    c = Cells()
    p = Perforations()
    f = Faces()

    well[:radius_grout, c] = radius_grout
    well[:pipe_spacing, p] = pipe_spacing
    well[:material_heat_capacity, c] = heat_capacity_grout
    well[:material_density, c] = density_grout

    treat_defaulted(x) = x
    treat_defaulted(::Missing) = NaN
    treat_defaulted(::Nothing) = NaN

    λpg = treat_defaulted(pipe_grout_thermal_index)
    λgg = treat_defaulted(grout_grout_thermal_index)
    well[:pipe_grout_thermal_index, p] = λpg
    well[:grout_grout_thermal_index, p] = λgg

end

function set_default_closed_loop_thermal_indices_u1!(well::DataDomain, pipe_cells, grout_cells)
    
    cell = Cells()
    face = Faces()
    perf = Perforations()

    num_nodes = well.representation.num_nodes
    hole_volumes = zeros(Float64, num_nodes)
    casing_volumes = zeros(Float64, num_nodes)
    grout_volumes = zeros(Float64, num_nodes)
    N = well.representation.neighborship

    rep = physical_representation(well)
    for (pno, grout_cell) in enumerate(rep.perforations.self)

        @assert grout_cell ∈ grout_cells
        segs = findall(N[2, :] .== grout_cell)
        @assert length(segs) <= 2
        seg_pg = segs[1]
        pipe_cell = N[1, seg_pg]
        @assert pipe_cell ∈ pipe_cells
        seg_gg = nothing
        if length(segs) > 1
            seg_gg = segs[2]
            other_grout_cell = N[1, seg_gg]
        end

        r_grout = well[:radius_grout, cell][grout_cell]
        r_pipe = well[:radius, cell][pipe_cell]
        wall_thickness = well[:casing_thickness, cell][pipe_cell]
        L = well[:cell_length, cell][pipe_cell]
        vol_p, vol_w, vol_g = closed_loop_volume_u1(
            L, r_grout, r_pipe, wall_thickness
        )
       
        hole_volumes[pipe_cell] = vol_p
        hole_volumes[grout_cell] = 0.0
        casing_volumes[pipe_cell] = 0.0
        casing_volumes[grout_cell] = 0.0
        grout_volumes[pipe_cell] = 0.0
        grout_volumes[grout_cell] = vol_g
        
        pipe_spacing = well[:pipe_spacing, perf][pno]
        λg = well[:grouting_thermal_conductivity, cell][grout_cell]
        λp = well[:casing_thermal_conductivity, cell][pipe_cell]

        Rpg, Rgr, Rgg = closed_loop_thermal_resistance_u1(
            r_grout, r_pipe + wall_thickness, wall_thickness, pipe_spacing, λg, λp)
        λpg = L/Rpg
        λgr = L/Rgr
        λgg = L/Rgg

        if isnan(well[:thermal_well_index, perf][pno])
            well[:thermal_well_index, perf][pno] = λgr
        end

        if well[:material_thermal_conductivity, face][seg_pg] == 0.0
            well[:material_thermal_conductivity, face][seg_pg] = λpg
        end

        if seg_gg !== nothing && well[:material_thermal_conductivity, face][seg_gg] == 0.0
            well[:material_thermal_conductivity, face][seg_gg] = λgg
        end
    end

    for pipe_cell in pipe_cells
        if hole_volumes[pipe_cell] == 0.0
            # This pipe cell was not assigned a volume, likely because it is not
            # connected to a grout cell. Assign a default volume based on geometry.
            r_pipe = well[:radius, cell][pipe_cell]
            wall_thickness = well[:casing_thickness, cell][pipe_cell]
            L = well[:cell_length, cell][pipe_cell]
            vol_p, vol_w, vol_g = closed_loop_volume_u1(
                L, 0.0, r_pipe, wall_thickness
            )
            hole_volumes[pipe_cell] = vol_p
            casing_volumes[pipe_cell] = 0.0
        end
    end

    well[:volume_override_hole, cell] = hole_volumes
    well[:volume_override_casing, cell] = casing_volumes
    well[:volume_override_grouting, cell] = grout_volumes
end

function closed_loop_volume_u1(length, radius_grout, radius_pipe, wall_thickness)

    # Compute pipe and grout volume
    L = length
    rg, rp_in, rp_out = radius_grout, radius_pipe - wall_thickness, radius_pipe
    vol_hole = π*rp_in^2*L
    vol_pipe = π*rp_out^2*L
    vol_wall = vol_pipe - vol_hole
    vol_grout = π*rg^2*L/2 - vol_pipe

    return vol_hole, vol_wall, vol_grout, L

end

function closed_loop_thermal_resistance_u1(
    radius_grout, radius_pipe, wall_thickness_pipe, pipe_spacing,
    thermal_conductivity_grout, thermal_conductivity_pipe)
    # Conveient short-hand notation
    radius_pipe_outer = radius_pipe
    radius_pipe_inner = radius_pipe - wall_thickness_pipe
    rg, rp_in, rp_out, w = 
        radius_grout, radius_pipe_inner, radius_pipe_outer, pipe_spacing
    λg, λp = thermal_conductivity_grout, thermal_conductivity_pipe

    ## Compute thermal resistances
    # Advection-dependent pipe thermal resistance
    Ra = 0.0 # TODO: Implement this
    # Conduction-dependent pipe thermal resistance
    Rc_a = log(rp_out/rp_in)/(2*π*λp)
    dg, dp_out = 2*rg, 2*rp_out
    x = log(sqrt(dg^2 + 2*dp_out^2)/(2*dp_out))/log(dg/(sqrt(2)*dp_out))
    # Grout thermal resistance
    Rg = acosh((dg^2 + dp_out^2 - w^2)/(2*dg*dp_out))/(2*π*λg)*(1.601 - 0.888*w/dg)
    # Conduction-dependent grout thermal resistance
    Rc_b = x*Rg
    # Combined thermal resistance of pipe and grout
    Rpg = Ra + Rc_a + Rc_b
    Rgr = (1-x)*Rg
    # Grout-to-grout thermal resistance
    Rar = acosh((2*w^2 - dp_out^2)/dp_out^2)/(2*π*λg)
    Rgg = 2*Rgr*(Rar - 2*x*Rg)/(2*Rgr - Rar + 2*x*Rg)

    return Rpg, Rgr, Rgg

end