function setup_closedloop_well_coaxial(D::DataDomain, reservoir_cells;
    cell_centers = D[:cell_centroids],
    name = :ClosedLoop,
    radius_grout = 65e-3,
    radius_inner_pipe = 15e-3,
    radius_outer_pipe = 20e-3,
    pipe_thickness_outer = 3e-3,
    thermal_conductivity_grout = 2.3,
    thermal_conductivity_pipe = 0.38,
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

    radius_pipe = zeros(num_cells)
    radius_pipe[inner_ix] = radius_inner_pipe
    radius_outer[outer_ix] = radius_outer_pipe

    pipe_thickness = zeros(num_cells)
    pipe_thickness[inner_ix] = pipe_thickness_inner
    pipe_thickness[outer_ix] = pipe_thickness_outer

    thermal_conductivity_pipe = zeros(num_cells)
    thermal_conductivity_pipe[inner_ix] = thermal_conductivity_inner_pipe
    thermal_conductivity_pipe[outer_ix] = thermal_conductivity_outer_pipe

    thermal_conductivity_grout0 = thermal_conductivity_grout
    thermal_conductivity_grout = zeros(num_cells)
    thermal_conductivity_grout[grout_cells] = thermal_conductivity_grout0

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
        radius = radius,
        casing_thickness = pipe_thickness,
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

    set_default_btes_thermal_indices!(BTESTypeU1(), supply_well)
    # set_default_btes_thermal_indices!(return_well)

    return [supply_well, return_well]

end

function augment_btes_domain!(type::BTESTypeCoaxial, well::DataDomain,
    radius_grout,
    heat_capacity_grout,
    density_grout;
    pipe_grout_thermal_index = missing,
)

    c = Cells()
    p = Perforations()
    f = Faces()

    well[:radius_grout, c] = radius_grout
    well[:pipe_spacing, p] = pipe_spacing
    # well[:heat_capacity_grout, p] = heat_capacity_grout
    well[:casing_heat_capacity, c] = heat_capacity_grout
    well[:casing_density, c] = density_grout

    treat_defaulted(x) = x
    treat_defaulted(::Missing) = NaN
    treat_defaulted(::Nothing) = NaN

    λpg = treat_defaulted(pipe_grout_thermal_index)
    well[:pipe_grout_thermal_index, p] = λpg

end


struct BTESTypeCoaxial <: AbstractBTESType end

function btes_volume(type::BTESTypeCoaxial, length, radius_grout, radius_pipe_inner, radius_pipe_outer)

    # Compute pipe and grout volume
    L = length
    rg, rp_in, rp_out = radius_grout, radius_pipe_inner, radius_pipe_outer
    vol_hole = π*rp_in^2*L
    vol_pipe = π*rp_out^2*L
    vol_wall = vol_pipe - vol_hole
    vol_grout = π*rg^2*L/2 - vol_pipe

    return vol_hole, vol_wall, vol_grout, L

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
        radius_inner_pipe,
        radius_inner_pipe + wall_thickness_inner_pipe,
        radius_outer_pipe,
        radius_outer_pipe + wall_thickness_outer_pipe,
        length

    λ_g, λ_pi, λ_po =
    thermal_conductivity_grout,
    thermal_conductivity_inner_pipe,
    thermal_conductivity_outer_pipe

    R_ai = 0.0 # TODO: Implement this
    R_ao = 0.0

    R_ci = log(r_ipo/r_ipi)/(2*π*λ_pi)
    R_co = log(r_opo/r_opi)/(2*π*λ_po)

    d_g, d_ipo, d_opo = 2*r_g, 2*r_ipo, 2*r_opo
    x = log(sqrt(d_g^2 + d_opo^2)/(2*d_opo^2))/log(d_g/d_opo)
    R_cg = (1-x)*log(d_g/d_opo)/(2*π*λ_g)

    R_pp = R_ai + R_ci
    R_pg = R_ao + R_co + R_cg

    λ_pp = L/R_pp
    λ_pg = L/R_pg

    return λ_pp, λ_pg

end