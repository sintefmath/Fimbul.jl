using LinearAlgebra

# Convert missing/nothing to NaN for defaulted thermal indices
treat_defaulted(x) = x
treat_defaulted(::Missing) = NaN
treat_defaulted(::Nothing) = NaN

## Closed-loop well types

"""
    AbstractClosedLoopWell

Abstract type for closed-loop wells. A closed-loop well is set up as a single
multisegment well with two controllable parts: a supply well, which holds the
full closed-loop topology (pipes and grout), and a return well, which is a
single outlet node connected to the end node of the supply well by a regular
well segment. The two parts are merged into one well domain using
`JutulDarcy.WellMerging.merge_well_domains`, so that each part can be
controlled independently without the need for a cross term between them. The
names of the two parts are given by [`closed_loop_well_names`](@ref).

The closed-loop type is stored as the `type` of the `MultiSegmentWell`.
Concrete types implement [`setup_closed_loop_supply`](@ref).
"""
abstract type AbstractClosedLoopWell end

"""
    GenericClosedLoopWell()

Closed-loop well with direct thermal communication between the pipe and the
reservoir. If no topology is given, the well is set up as a simple U-shaped
pipe. A general topology (e.g., for advanced geothermal systems with multiple
laterals) can be given using `neighborship`, `perforation_cells_well`,
`well_cell_centers` and `end_nodes`.
"""
struct GenericClosedLoopWell <: AbstractClosedLoopWell end

"""
    U1ClosedLoopWell()

Single U-tube closed-loop well with two pipes embedded in grout.
"""
struct U1ClosedLoopWell <: AbstractClosedLoopWell end

"""
    CoaxialClosedLoopWell()

Coaxial closed-loop well with an inner pipe and an outer annulus embedded in
grout.
"""
struct CoaxialClosedLoopWell <: AbstractClosedLoopWell end

function closed_loop_well_type(type::Symbol)
    if type == :generic || type == :simple
        return GenericClosedLoopWell()
    elseif type == :u1
        return U1ClosedLoopWell()
    elseif type == :coaxial
        return CoaxialClosedLoopWell()
    else
        error("Unknown closed loop type: $type")
    end
end

closed_loop_well_type(type::AbstractClosedLoopWell) = type

"""
    is_closed_loop_well(well)

Check if a well (domain or representation) is a closed-loop well.
"""
is_closed_loop_well(well::DataDomain) = is_closed_loop_well(physical_representation(well))
is_closed_loop_well(well::JutulDarcy.MultiSegmentWell) = well.type isa AbstractClosedLoopWell
is_closed_loop_well(well) = false

"""
    supply_name, return_name = closed_loop_well_names(type, name)
    supply_name, return_name = closed_loop_well_names(well::DataDomain)

Names of the supply and return parts of a closed-loop well, used for controls
and well results. The names are `name` with a suffix that depends on the
closed-loop type:

- [`GenericClosedLoopWell`](@ref): `_a` and `_b`
- [`U1ClosedLoopWell`](@ref): `_supply` and `_return`
- [`CoaxialClosedLoopWell`](@ref): `_inner` and `_outer`

For a merged closed-loop well domain, the names are taken from the domain.
"""
function closed_loop_well_names(type, name)
    suffix_supply, suffix_return = closed_loop_well_suffixes(closed_loop_well_type(type))
    return (Symbol(name, suffix_supply), Symbol(name, suffix_return))
end

function closed_loop_well_names(well::DataDomain)
    is_closed_loop_well(well) || error("Well is not a closed-loop well.")
    names = physical_representation(well).multiwell.names
    return (names[1], names[2])
end

closed_loop_well_suffixes(::AbstractClosedLoopWell) = ("_a", "_b")
closed_loop_well_suffixes(::U1ClosedLoopWell) = ("_supply", "_return")
closed_loop_well_suffixes(::CoaxialClosedLoopWell) = ("_inner", "_outer")

## Setup functions

"""
    setup_closed_loop_well(D::DataDomain, reservoir_cells;
        closed_loop_type = :generic, name = :CL, kwarg...)

Set up a closed-loop well in `reservoir_cells` of the reservoir domain `D`.
`closed_loop_type` is either an [`AbstractClosedLoopWell`](@ref) or one of
`:generic` (alias `:simple`), `:u1` and `:coaxial`. The result is a single
merged well domain named `name`, with supply and return parts named by
[`closed_loop_well_names`](@ref). The return part is connected to the
reservoir in `return_reservoir_cell` (defaults to the first reservoir cell)
with zero well indices.
"""
function setup_closed_loop_well(D::DataDomain, reservoir_cells;
    closed_loop_type = :generic,
    kwarg...)

    type = closed_loop_well_type(closed_loop_type)
    return setup_closed_loop_well(type, D, reservoir_cells; kwarg...)

end

function setup_closed_loop_well(type::AbstractClosedLoopWell, D::DataDomain, reservoir_cells;
    name = :CL,
    return_reservoir_cell = reservoir_cells[1],
    kwarg...)

    supply_name, return_name = closed_loop_well_names(type, name)
    supply_well = setup_closed_loop_supply(type, D, reservoir_cells;
        name = supply_name, kwarg...)
    return merge_closed_loop_well(supply_well, return_reservoir_cell;
        name = name, return_name = return_name)

end

"""
    setup_closed_loop_supply(type::AbstractClosedLoopWell, D, reservoir_cells; name, kwarg...)

Set up the supply part of a closed-loop well of the given type. This is a
multisegment well holding the full closed-loop topology, with exactly one end
node where the fluid leaves the loop.
"""
function setup_closed_loop_supply end

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

function setup_closed_loop_supply(type::GenericClosedLoopWell, D::DataDomain, reservoir_cells;
    cell_centers = D[:cell_centroids],
    neighborship = missing,
    perforation_cells_well = missing,
    well_cell_centers = missing,
    dir = :z,
    segment_models = missing,
    end_nodes = missing,
    tag = missing,
    name = :CL_a,
    radius_pipe = 20e-3,
    wall_thickness = 2.5e-3,
    grouting_thickness = 50e-3,
    pipe_thermal_conductivity = 0.35,
    grouting_thermal_conductivity = 2.3,
    pipe_spacing = 60e-3,
    kwarg...)

    if ismissing(neighborship)
        if !ismissing(perforation_cells_well) || !ismissing(end_nodes)
            error("Provide either neighborship, perforation_cells_well, \
                well_cell_centers, end_nodes, or none of them.")
        end
        # Create pipe cells
        reservoir_cells = vcat(reservoir_cells, reverse(reservoir_cells))
        if dir isa AbstractVector
            dir = vcat(dir, reverse(dir))
        end
        nc_pipe = length(reservoir_cells)
        pipe_cells = collect(1:nc_pipe)
        perforation_cells_well = pipe_cells
        # Set up connectivity
        neighborship = vcat(pipe_cells[1:end-1]', pipe_cells[2:end]')
        # Set up well cell centers and end nodes
        nc_mid = div(nc_pipe, 2)
        left_ix = Int.(1:nc_mid)
        right_ix = Int.(1:nc_mid) .+ nc_mid
        wc_left = cell_centers[:, reservoir_cells[left_ix]] .- [pipe_spacing/2, 0.0, 0.0]
        wc_right = cell_centers[:, reservoir_cells[right_ix]] .+ [pipe_spacing/2, 0.0, 0.0]
        well_cell_centers = hcat(wc_left, wc_right)
        end_nodes = [pipe_cells[end]]
        # Set tag for identifying cell groups (e.g., pipe legs)
        if !ismissing(tag)
            @warn "tag argument is ignored when neighborship is not \
                provided. Tags will be created automatically."
        end
        tag = Vector{Symbol}(undef, nc_pipe)
        tag[pipe_cells[left_ix]] .= :pipe_left
        tag[pipe_cells[right_ix]] .= :pipe_right
    else
        if ismissing(well_cell_centers) || ismissing(end_nodes)
            error("If neighborship is provided, well_cell_centers and \
                end_nodes must also be provided.")
        end
        if ismissing(perforation_cells_well)
            perforation_cells_well = collect(eachindex(reservoir_cells))
        end
    end
    if ismissing(segment_models)
        segment_models = nothing
    end

    supply_well = setup_well(D, reservoir_cells;
        name = name,
        type = type,
        simple_well = false,
        neighborship = neighborship,
        perforation_cells_well = perforation_cells_well,
        well_cell_centers = well_cell_centers,
        segment_models = segment_models,
        end_nodes = end_nodes,
        dir = dir,
        WI = 0.0,
        radius = radius_pipe,
        casing_thickness = wall_thickness,
        grouting_thickness = grouting_thickness,
        casing_thermal_conductivity = pipe_thermal_conductivity,
        grouting_thermal_conductivity = grouting_thermal_conductivity,
        kwarg...)
    if !ismissing(tag)
        supply_well[:tag, Cells()] = tag
    end

    return supply_well

end

## Merging of supply and return wells

"""
    merge_closed_loop_well(supply_well::DataDomain, return_reservoir_cell;
        name, return_name = Symbol(name, "_b"))

Merge the supply part of a closed-loop well with a return well consisting of a
single node at the end node of the supply well. The return node gets a copy of
the data in the supply end node, and a zero-length segment connects the two.
The return well is perforated in `return_reservoir_cell` with zero well
indices.
"""
function merge_closed_loop_well(supply_well::DataDomain, return_reservoir_cell;
    name, return_name = Symbol(name, "_b"))

    w = physical_representation(supply_well)
    length(w.end_nodes) == 1 || error("Closed-loop supply well must have exactly one end node.")
    end_node = only(w.end_nodes)
    return_well = setup_closed_loop_return_well(
        supply_well, return_reservoir_cell; name = return_name)
    merged = JutulDarcy.WellMerging.merge_well_domains(
        [supply_well, return_well], name)
    # Connect the supply end node to the (only) return node
    return_node = number_of_cells(w) + 1
    return add_closed_loop_connection(merged, end_node, return_node)

end

function setup_closed_loop_return_well(supply_well::DataDomain, return_reservoir_cell;
    name)

    w = physical_representation(supply_well)
    end_node = only(w.end_nodes)
    perf = findfirst(isequal(return_reservoir_cell), w.perforations.reservoir)
    if isnothing(perf)
        perf = 1
    end
    return_well = JutulDarcy.MultiSegmentWell(zeros(Int, 2, 0), [return_reservoir_cell], [1];
        type = w.type,
        name = name,
        end_nodes = [1],
        segment_models = eltype(w.segment_models)[],
        include_inertia = w.include_inertia,
        surface_conditions = w.surface
    )
    return_well = DataDomain(return_well)
    # Copy data from the supply end node (cells) and from the perforation in
    # the return reservoir cell (perforations)
    index = Dict(
        Cells() => [end_node],
        Perforations() => [perf],
        Faces() => Int[]
    )
    for key in keys(supply_well)
        entity = Jutul.associated_entity(supply_well, key)
        value = supply_well[key]
        if entity == NoEntity()
            return_well[key, entity] = value
        else
            haskey(index, entity) || error("Cannot set up return well with data $key on $entity.")
            return_well[key, entity] = select_entity_values(value, index[entity])
        end
    end
    # No flow or heat exchange with the reservoir
    return_well[:well_index, Perforations()] = 0.0
    return_well[:thermal_well_index, Perforations()] = 0.0
    # The return node holds fluid only
    for key in (:volume_override_casing, :volume_override_grouting)
        if haskey(return_well, key, Cells())
            return_well[key, Cells()] = 0.0
        end
    end
    if haskey(return_well, :tag, Cells())
        return_well[:tag, Cells()] = [:return]
    end

    return return_well

end

function add_closed_loop_connection(well::DataDomain, from_node, to_node;
    segment_model = SegmentWellBoreFrictionHB())

    w = physical_representation(well)
    info = w.multiwell
    # Add the segment to the last well (i.e., the return well)
    neighborship = hcat(w.neighborship, [from_node, to_node])
    segment_models = vcat(w.segment_models, [segment_model])
    face_counts = diff(info.faces.pos)
    face_counts[end] += 1
    faces = JutulDarcy.WellMerging.entity_indirection_map(face_counts)
    info = JutulDarcy.MultiWellInfo(info.names, info.domains, info.nodes,
        faces, info.perforations, info.top_nodes)
    connected_well = JutulDarcy.MultiSegmentWell(w.type, w.num_nodes,
        w.num_segments + 1, w.num_perforations, w.perforations, neighborship,
        w.end_nodes, w.surface, w.name, segment_models, w.include_inertia, info)

    connected = DataDomain(connected_well)
    # Face data for the new segment is copied from a segment connected to the
    # from node, with no heat conduction in the well material
    face = findfirst(seg -> from_node in seg, collect(eachcol(w.neighborship)))
    if isnothing(face)
        face = 1
    end
    for key in keys(well)
        entity = Jutul.associated_entity(well, key)
        value = well[key]
        if entity == Faces()
            value = cat(value, select_entity_values(value, [face]); dims = ndims(value))
            if key == :material_thermal_conductivity
                value[end] = 0.0
            end
        end
        connected[key, entity] = value
    end

    return connected

end

select_entity_values(v::AbstractVector, ix) = v[ix]
select_entity_values(v::AbstractMatrix, ix) = v[:, ix]
