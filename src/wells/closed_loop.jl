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
    name = :CL,
    radius_pipe = 20e-3,
    wall_thickness = 2.5e-3,
    grouting_thickness = 50e-3,
    pipe_thermal_conductivity = 0.35,
    grouting_thermal_conductivity = 2.3,
    kwarg...)

    # Common properties
    args = (
        WI = 0.0,
        radius = radius_pipe,
        casing_thickness = wall_thickness,
        grouting_thickness = grouting_thickness,
        casing_thermal_conductivity = pipe_thermal_conductivity,
        grouting_thermal_conductivity = grouting_thermal_conductivity,
        end_nodes = [length(reservoir_cells)],
        simple_well = false,
        type = :closed_loop
    )

    # Set up supply and return wells
    supply_well = setup_well(D::DataDomain, reservoir_cells;
        name = Symbol(name, "_supply"),
        args..., kwarg...)
    return_well = setup_well(D::DataDomain, reservoir_cells;
        name = Symbol(name, "_return"),
        args..., kwarg...)

    return supply_well, return_well

end