using LinearAlgebra

# Utility functions for setting up BTES wells
function setup_closed_loop_well(D::DataDomain, reservoir_cells;
    closed_loop_type = :simple,
    name = :CL, kwarg...)

    if closed_loop_type == :simple
        # Simple closed loop well with direct pipe/reservoir thermal communication
        return setup_closed_loop_well_simple(D, reservoir_cells; name = name, kwarg...)
    elseif closed_loop_type == :u1
        # U-type closed loop well with grout annulus
        return setup_closed_loop_well_u1(D, reservoir_cells; name = name, kwarg...)
    elseif closed_loop_type == :coaxial
        # Coaxial closed loop well with grout annulus
        return setup_closed_loop_well_coaxial(D, reservoir_cells; name = name, kwarg...)
    else
        # Unknown closed loop type
        error("Unknown closed loop type: $closed_loop_type")
    end

end

function setup_vertical_closed_loop_well(D::DataDomain, i, j;
    heel = 1, toe = missing, kwarg...)

    # Get reservoir cells from ijk-indices
    g = physical_representation(D)
    if ismissing(toe)
        toe = grid_dims_ijk(g)[3]
    end
    @assert heel <= toe
    @assert heel > 0
    @assert toe > 0
    k_range = heel:toe
    n = length(k_range)
    @assert n > 0
    reservoir_cells = zeros(Int64, n)
    for (ix, k) in enumerate(k_range)
        reservoir_cells[ix] = cell_index(g, (i, j, k))
    end
    # Set up BTES well
    return setup_closed_loop_well(D, reservoir_cells; kwarg...)
end

function setup_btes_well(D::DataDomain, reservoir_cells;
    name = :BTES, kwarg...)

    return setup_closed_loop_well(D, reservoir_cells;
        name = name, kwarg...)

end

function setup_vertical_btes_well(D::DataDomain, i, j;
    name = :BTES, heel = 1, toe = missing, kwarg...)

    return setup_vertical_closed_loop_well(D, i, j;
        name = name, heel = heel, toe = toe, kwarg...)

end

function setup_closed_loop_well_simple(D::DataDomain, reservoir_cells;
    return_reservoir_cell = missing,
    cell_centers = D[:cell_centroids],
    neighborship = missing,
    well_cell_centers = missing,
    segment_models = missing,
    end_nodes = missing,
    section = missing,
    name = :CL,
    radius_pipe = 20e-3,
    wall_thickness = 2.5e-3,
    grouting_thickness = 50e-3,
    pipe_thermal_conductivity = 0.35,
    grouting_thermal_conductivity = 2.3,
    kwarg...)

    if ismissing(neighborship)
        # Create pipe cells
        reservoir_cells = vcat(reservoir_cells, reverse(reservoir_cells))
        nc_pipe = length(reservoir_cells)
        pipe_cells = collect(1:nc_pipe)
        # Set up connectivity
        pipe_to_pipe = vcat(pipe_cells[1:end-1]', pipe_cells[2:end]')
        neighborship = hcat(pipe_to_pipe)
        # Set up well cell centers and end nodes
        nc_mid = div(nc_pipe, 2)
        left_ix = Int.(1:nc_mid)
        right_ix = Int.(1:nc_mid) .+ nc_mid
        pipe_spacing = 60e-3
        wc_left = cell_centers[:, reservoir_cells[1:nc_mid]] .- [pipe_spacing/2, 0.0, 0.0]
        wc_right = cell_centers[:, reservoir_cells[nc_mid+1:end]] .+ [pipe_spacing/2, 0.0, 0.0]
        well_cell_centers = hcat(wc_left, wc_right)
        end_nodes = [pipe_cells[end]]
        # Set section for easy lookup
        if !ismissing(section)
            @warn(["section argument is ignored when neighborship is not provided. ",
                "Sections will be created automatically."])
        end
        section = Vector{Any}(undef, nc_pipe)
        section[pipe_cells[left_ix]] .= [(1, :pipe_left)]
        section[pipe_cells[right_ix]] .= [(1, :pipe_right)]
    else
        if ismissing(well_cell_centers) || ismissing(end_nodes) 
            error("""If neighborship is provided, well_cell_centers, and 
            end_nodes must also be provided.""")
        end
    end
    # Set return reservoir cell if missing
    if ismissing(return_reservoir_cell)
        return_reservoir_cell = reservoir_cells[end]
    end

    # Common properties
    args = (
        WI = 0.0,
        radius = radius_pipe,
        simple_well = false,
        type = :closed_loop
    )

    # Set up supply and return wells
    supply_well = setup_well(D::DataDomain, reservoir_cells;
        name = Symbol(name, "_supply"),
        casing_thickness = wall_thickness,
        grouting_thickness = grouting_thickness,
        casing_thermal_conductivity = pipe_thermal_conductivity,
        grouting_thermal_conductivity = grouting_thermal_conductivity,
        end_nodes = end_nodes,
        args..., kwarg...)

    return_well = setup_well(D::DataDomain, return_reservoir_cell;
        name = Symbol(name, "_return"),
        WIth = 0.0,
        casing_thickness = wall_thickness,
        args..., kwarg...)

    return [supply_well, return_well]

end