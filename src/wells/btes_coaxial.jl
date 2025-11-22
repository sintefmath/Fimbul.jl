struct BTESTypeCoaxial <: AbstractBTESType end

function setup_btes_well_coaxial(D::DataDomain, reservoir_cells;
    cell_centers = D[:cell_centroids],
    name = :BTES,
    radius_grout = 65e-3,
    radius_inner_pipe = 15e-3,
    wall_thickness_inner = 3e-3,
    radius_outer_pipe = 20e-3,
    wall_thickness_outer = 3e-3,
    thermal_conductivity_grout = 2.3,
    thermal_conductivity_inner_pipe = 0.38,
    thermal_conductivity_outer_pipe = 0.38,
    heat_capacity_grout = 1460,
    density_grout = 1500.0,
    kwarg...)

    ## Set up connectivity
    reservoir_cells = vcat(reservoir_cells, reverse(reservoir_cells))
    nc_pipe = length(reservoir_cells)
    nc_grout = length(reservoir_cells)÷2
    nc_mid = div(nc_pipe, 2)
    inner_ix = Int.(1:nc_mid)
    outer_ix = Int.(1:nc_mid) .+ nc_mid

    pipe_cells = (1:nc_pipe)
    grout_cells = (1:nc_grout) .+ nc_pipe
    num_cells = nc_pipe + nc_grout

    pipe_to_pipe = vcat(pipe_cells[1:end-1]', pipe_cells[2:end]')
    pipe_to_grout = vcat(pipe_cells[outer_ix]', grout_cells')
    
    N = hcat(pipe_to_pipe, pipe_to_grout)

    ## Set up segment flow models
    segment_models = Vector{Any}()

    # Set centers and depths
    well_cell_centers = repeat(cell_centers[:, reservoir_cells], 1, 2)
    
    # Add top node
    nseg = size(N,2)

    # Set material thermal conducivities
    nr = length(reservoir_cells)

    for seg in eachcol(N)
        if all([c ∈ pipe_cells for c in seg])
            # Pipe segments use standard wellbore friction model
            seg_model = SegmentWellBoreFrictionHB()
        else
            seg_model = JutulDarcy.ClosedSegment()
        end
        push!(segment_models, seg_model)
    end

    cell_radius = zeros(num_cells)
    cell_radius_inner = zeros(num_cells)
    cell_radius[inner_ix] .= radius_inner_pipe - wall_thickness_inner
    cell_radius[outer_ix] .= radius_outer_pipe - wall_thickness_outer
    cell_radius_inner[outer_ix] .= radius_inner_pipe
    
    wall_thickness = zeros(num_cells)
    wall_thickness[inner_ix] .= wall_thickness_inner
    wall_thickness[outer_ix] .= wall_thickness_outer

    thermal_conductivity_pipe = zeros(num_cells)
    thermal_conductivity_pipe[inner_ix] .= thermal_conductivity_inner_pipe
    thermal_conductivity_pipe[outer_ix] .= thermal_conductivity_outer_pipe

    thermal_conductivity_grout0 = thermal_conductivity_grout
    thermal_conductivity_grout = zeros(num_cells)
    thermal_conductivity_grout[grout_cells] .= thermal_conductivity_grout0

    ## Set up supply and return wells
    args = (
        type = :closed_loop,
        simple_well = false,
        WI = 0.0
    )

    supply_well = setup_well(D, reservoir_cells;
        name = Symbol(name, "_supply"),
        neighborship = N,
        perforation_cells_well = collect(grout_cells),
        well_cell_centers = well_cell_centers,
        cell_radius = cell_radius,
        casing_thickness = wall_thickness,
        thermal_conductivity_casing = thermal_conductivity_pipe,
        thermal_conductivity_grout = thermal_conductivity_grout,
        segment_models = segment_models,
        end_nodes = [nc_pipe],
        args..., kwarg...)

    augment_btes_domain!(supply_well,
        radius_grout,
        pipe_spacing,
        heat_capacity_grout,
        density_grout
    )
    
    return_well = setup_well(D, reservoir_cells[end];
        name = Symbol(name, "_return"),
        args...)

    set_default_btes_thermal_indices!(BTESTypeCoaxial(), supply_well)
    # set_default_btes_thermal_indices!(return_well)

    return [supply_well, return_well]

end

function augment_btes_domain!(type::BTESTypeCoaxial, well::DataDomain,
    radius_grout,
    heat_capacity_grout,
    density_grout;
    pipe_pipe_thermal_index = missing,
    pipe_grout_thermal_index = missing
)

    c = Cells()
    p = Perforations()
    f = Faces()

    well[:radius_grout, c] = radius_grout
    # well[:heat_capacity_grout, p] = heat_capacity_grout
    well[:casing_heat_capacity, c] = heat_capacity_grout
    well[:casing_density, c] = density_grout

    treat_defaulted(x) = x
    treat_defaulted(::Missing) = NaN
    treat_defaulted(::Nothing) = NaN

    λpp = treat_defaulted(pipe_pipe_thermal_index)
    λpg = treat_defaulted(pipe_grout_thermal_index)
    well[:pipe_pipe_thermal_index, p] = λpp
    well[:pipe_grout_thermal_index, p] = λpg

end

function set_default_btes_thermal_indices!(type::BTESTypeCoaxial, well::DataDomain)
    
    cell = Cells()
    face = Faces()
    perf = Perforations()

    num_nodes = well.representation.num_nodes

    hole_volumes = zeros(Float64, num_nodes)
    casing_volumes = zeros(Float64, num_nodes)
    grout_volumes = zeros(Float64, num_nodes)

    nc = Int.(num_nodes/3)
    inner_pipe_cells = Int.(1:nc)
    outer_pipe_cells = Int.(1:nc) .+ nc
    grout_cells = Int.(1:nc) .+ 2*nc
    N = well.representation.neighborship

    rep = physical_representation(well)
    for (pno, grout_cell) in enumerate(rep.perforations.self)

        @assert grout_cell ∈ grout_cells
        seg_pg = findall(N[2, :] .== grout_cell)
        @assert length(seg_pg) == 1
        outer_pipe_cell = N[1, seg_pg[1]]
        @assert outer_pipe_cell ∈ outer_pipe_cells
        seg_pp = findall(N[2, :] .== outer_pipe_cell)
        @assert length(seg_pp) == 1
        inner_pipe_cell = N[1, seg_pp[1]]
        @assert inner_pipe_cell ∈ inner_pipe_cells

        r_grout = well[:radius_grout, cell][grout_cell]
        r_inner_pipe = well[:radius, cell][inner_pipe_cell]
        wall_thickness_inner = well[:casing_thickness, cell][inner_pipe_cell]
        r_outer_pipe = well[:radius, cell][outer_pipe_cell]
        wall_thickness_outer = well[:casing_thickness, cell][outer_pipe_cell]
        L = well[:cell_length, cell][inner_pipe_cell]
        vol_ip, vol_iw, vol_op, vol_ow, vol_g, L = btes_volume(
            BTESTypeCoaxial(),
            L, r_grout,
            r_inner_pipe + wall_thickness_inner, wall_thickness_inner,
            r_outer_pipe + wall_thickness_outer, wall_thickness_outer
        )
       
        hole_volumes[inner_pipe_cell] = vol_ip
        hole_volumes[outer_pipe_cell] = vol_op
        hole_volumes[grout_cell] = 0.0
        casing_volumes[inner_pipe_cell] = vol_iw
        casing_volumes[outer_pipe_cell] = vol_ow
        casing_volumes[grout_cell] = 0.0
        grout_volumes[inner_pipe_cell] = 0.0
        grout_volumes[outer_pipe_cell] = 0.0
        grout_volumes[grout_cell] = vol_g
        
        λg = well[:thermal_conductivity_grout, cell][grout_cell]
        λp = well[:thermal_conductivity_casing, cell][pipe_cell]
        λpp, λpg, λgr = btes_thermal_conductivity(
            BTESTypeCoaxial(),
            r_grout, r_pipe + wall_thickness, wall_thickness, pipe_spacing, L, λg, λp)

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

function btes_volume(type::BTESTypeCoaxial,
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

function btes_thermal_conductivity(type::BTESTypeCoaxial,
    radius_grout,
    radius_inner_pipe, wall_thickness_inner_pipe, 
    radius_outer_pipe, wall_thickness_outer_pipe,
    length,
    thermal_conductivity_grout,
    thermal_conductivity_inner_pipe,
    thermal_conductivity_outer_pipe
    )

    r_g, r_ipi, r_ipo, r_opi, r_opo, L =
        radius_grout,
        radius_inner_pipe - wall_thickness_inner_pipe,
        radius_inner_pipe,
        radius_outer_pipe - wall_thickness_outer_pipe,
        radius_outer_pipe,
        length

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
    println("x = $x")
    R_g = log(d_g/d_opo)/(2*π*λ_g)
    R_cg = x*R_g

    R_pp = R_ai + R_ci
    R_pg = R_ao + R_co + R_cg
    R_gr = (1-x)*R_g
    
    λpp = L/R_pp
    λpg = L/R_pg
    λgr = L/R_gr

    return λpp, λpg, λgr

end