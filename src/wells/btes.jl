using LinearAlgebra

## Utility functions for setting up BTES wells
function setup_btes_well(D::DataDomain, reservoir_cells;
    btes_type = :simple,
    name = :BTES, kwarg...)

    if btes_type == :simple
        # Simple BTES well with direct pipe/reservoir thermal communication
        return setup_btes_well_simple(D, reservoir_cells; name = name, kwarg...)
    else
        # Unknown BTES type
        error("Unknown BTES type: $btes_type")
    end

end

function setup_vertical_btes_well(D::DataDomain, i, j;
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
    return setup_btes_well(D, reservoir_cells; kwarg...)
end

function setup_btes_well_simple(D::DataDomain, reservoir_cells;
    name = :BTES,
    radius_pipe = 20e-3,
    pipe_thickness = 2.5e-3,
    grouting_thickness = 50e-3,
    thermal_conductivity_pipe = 0.35,
    kwarg...)

    # Common properties
    args = (
        WI = 0.0,
        radius = radius_pipe,
        casing_thickness = pipe_thickness,
        grouting_thickness = grouting_thickness,
        thermal_conductivity_casing = thermal_conductivity_pipe,
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