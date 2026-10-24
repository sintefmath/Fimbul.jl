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
    ClosedLoopGenericWell()

Closed-loop well with direct thermal communication between the pipe and the
reservoir. If no topology is given, the well is a single pipe that starts in
the first and ends in the last of the given reservoir cells. A general
topology (e.g., for advanced geothermal systems with multiple laterals) can be
given using `neighborship`, `perforation_cells_well`, `well_cell_centers` and
`end_nodes`.
"""
struct ClosedLoopGenericWell <: AbstractClosedLoopWell end

"""
    ClosedLoop1UWell()

Single U-tube closed-loop well with two pipes embedded in grout.
"""
struct ClosedLoop1UWell <: AbstractClosedLoopWell end

"""
    ClosedLoopCoaxialWell()

Coaxial closed-loop well with an inner pipe and an outer annulus embedded in
grout.
"""
struct ClosedLoopCoaxialWell <: AbstractClosedLoopWell end

function closed_loop_well_type(type::Symbol)
    if type == :generic
        return ClosedLoopGenericWell()
    elseif type == :u1
        return ClosedLoop1UWell()
    elseif type == :coaxial
        return ClosedLoopCoaxialWell()
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

- [`ClosedLoopGenericWell`](@ref) and [`ClosedLoop1UWell`](@ref): `_supply`
  and `_return`
- [`ClosedLoopCoaxialWell`](@ref): `_inner` and `_outer`

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

closed_loop_well_suffixes(::AbstractClosedLoopWell) = ("_supply", "_return")
closed_loop_well_suffixes(::ClosedLoopCoaxialWell) = ("_inner", "_outer")

## Setup functions

"""
    setup_closed_loop_well(D::DataDomain, reservoir_cells;
        closed_loop_type = :generic, name = :CL, kwarg...)

Set up a closed-loop well in `reservoir_cells` of the reservoir domain `D`.
`closed_loop_type` is either an [`AbstractClosedLoopWell`](@ref) or one of
`:generic`, `:u1` and `:coaxial`. The result is a single
merged well domain named `name`, with supply and return parts named by
[`closed_loop_well_names`](@ref). The return part is connected to the
reservoir in `return_reservoir_cell` with zero well indices. By default, this
is the reservoir cell perforated by the end node of the supply part, or the
first perforated reservoir cell if the end node is not perforated.

# Pressure in the loop
The fluid in a closed loop is a single liquid phase, and the pressure must stay
above the minimum pressure of the model everywhere in the loop. If one part is
controlled by bottom-hole pressure (BHP) at a node below the highest point of
the loop, the BHP must exceed the hydrostatic pressure of the fluid column
above that node, ρ g Δz, plus the minimum pressure. Otherwise the simulation
fails to converge (physically, the siphon would break). This is typically the
case for a generic well that flows downwards to an outlet at depth. Use
[`check_closed_loop_pressure`](@ref) to check the controls.
"""
function setup_closed_loop_well(D::DataDomain, reservoir_cells;
    closed_loop_type = :generic,
    kwarg...)

    type = closed_loop_well_type(closed_loop_type)
    return setup_closed_loop_well(type, D, reservoir_cells; kwarg...)

end

function setup_closed_loop_well(type::AbstractClosedLoopWell, D::DataDomain, reservoir_cells;
    name = :CL,
    return_reservoir_cell = missing,
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

"""
    setup_vertical_closed_loop_well(D::DataDomain, i, j; heel = 1, toe = missing, kwarg...)

Set up a vertical closed-loop well in column `(i, j)` of a logically
structured mesh, from layer `heel` to layer `toe` (defaults to the bottom
layer). See [`setup_closed_loop_well`](@ref) for the keyword arguments.
"""
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
    return setup_closed_loop_well(D, reservoir_cells; kwarg...)
end

"""
    setup_btes_well(D::DataDomain, reservoir_cells; name = :BTES, closed_loop_type = :u1, kwarg...)

Set up a borehole thermal energy storage (BTES) well. This is a closed-loop
well (see [`setup_closed_loop_well`](@ref)) that defaults to a 1U well.
"""
function setup_btes_well(D::DataDomain, reservoir_cells;
    name = :BTES, closed_loop_type = :u1, kwarg...)

    return setup_closed_loop_well(D, reservoir_cells;
        name = name, closed_loop_type = closed_loop_type, kwarg...)

end

"""
    setup_vertical_btes_well(D::DataDomain, i, j; name = :BTES, closed_loop_type = :u1, heel = 1, toe = missing, kwarg...)

Set up a vertical BTES well in column `(i, j)` of a logically structured mesh.
See [`setup_btes_well`](@ref) and [`setup_vertical_closed_loop_well`](@ref).
"""
function setup_vertical_btes_well(D::DataDomain, i, j;
    name = :BTES, closed_loop_type = :u1, heel = 1, toe = missing, kwarg...)

    return setup_vertical_closed_loop_well(D, i, j;
        name = name, closed_loop_type = closed_loop_type,
        heel = heel, toe = toe, kwarg...)

end

"""
    setup_closed_loop_supply(type::ClosedLoopGenericWell, D, reservoir_cells; kwarg...)

Supply part of a generic closed-loop well. By default, this is a single pipe
that starts in the first and ends in the last of `reservoir_cells`. All nodes
are tagged `:pipe`.

# Keyword arguments
- `radius_pipe = 20e-3`: Outer pipe radius [m].
- `wall_thickness_pipe = 2.5e-3`: Pipe wall thickness [m].
- `grouting_thickness = 50e-3`: Grout thickness outside the pipe [m].
- `pipe_thermal_conductivity = 0.35`: Pipe wall thermal conductivity [W/(m K)].
- `grouting_thermal_conductivity = 2.3`: Grout thermal conductivity [W/(m K)].
- `dir = :z`: Direction of the well in each reservoir cell.
- `neighborship`, `perforation_cells_well`, `well_cell_centers`, `end_nodes`:
  General well topology. `neighborship`, `well_cell_centers` and `end_nodes`
  must be given together, with exactly one end node. `perforation_cells_well`
  defaults to one node per reservoir cell.
- `segment_models`, `tag`: Segment models and node tags for a general topology.
"""
function setup_closed_loop_supply(type::ClosedLoopGenericWell, D::DataDomain, reservoir_cells;
    cell_centers = D[:cell_centroids],
    neighborship = missing,
    perforation_cells_well = missing,
    well_cell_centers = missing,
    dir = :z,
    segment_models = missing,
    end_nodes = missing,
    tag = missing,
    name = :CL_supply,
    radius_pipe = 20e-3,
    wall_thickness_pipe = 2.5e-3,
    grouting_thickness = 50e-3,
    pipe_thermal_conductivity = 0.35,
    grouting_thermal_conductivity = 2.3,
    kwarg...)

    nc = length(reservoir_cells)
    if ismissing(neighborship)
        check_closed_loop_topology_args(perforation_cells_well,
            well_cell_centers, end_nodes)
        # Single pipe from the first to the last reservoir cell
        pipe_cells = collect(1:nc)
        perforation_cells_well = pipe_cells
        neighborship = vcat(pipe_cells[1:end-1]', pipe_cells[2:end]')
        end_nodes = [pipe_cells[end]]
        well_cell_centers = cell_centers[:, reservoir_cells]
    else
        if ismissing(well_cell_centers) || ismissing(end_nodes)
            error("If neighborship is provided, well_cell_centers and \
                end_nodes must also be provided.")
        end
        if ismissing(perforation_cells_well)
            perforation_cells_well = collect(1:nc)
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
        casing_thickness = wall_thickness_pipe,
        grouting_thickness = grouting_thickness,
        casing_thermal_conductivity = pipe_thermal_conductivity,
        grouting_thermal_conductivity = grouting_thermal_conductivity,
        kwarg...)
    if ismissing(tag)
        tag = fill(:pipe, number_of_cells(physical_representation(supply_well)))
    end
    supply_well[:tag, Cells()] = tag

    return supply_well

end

# Topology arguments are only used together with neighborship
function check_closed_loop_topology_args(args...)
    if any(!ismissing, args)
        error("Topology arguments (e.g., cell groups, well_cell_centers and \
            end_nodes) can only be given together with neighborship.")
    end
end

## Pressure check

"""
    check_closed_loop_pressure(model, controls; density = missing, p_min = missing)

Check that bottom-hole pressure (BHP) controls of closed-loop wells keep the
pressure in the loop above the minimum pressure of the model. `controls` is a
`Dict` of well controls (as passed to `setup_reservoir_forces`), forces with a
`:Facility` entry, or a vector of forces.

For each part of a closed-loop well with a BHP target, the hydrostatic
pressure at the highest node of the loop is estimated as
`bhp - density*g*(z_bhp - z_top)`, where `z_bhp` is the depth of the node where
the BHP applies. A warning is given if this is below `p_min`. `density`
defaults to the largest reference density of the well system, and `p_min` to
the minimum pressure of the well model. Friction is not accounted for.

Returns `true` if all checks pass, `false` otherwise.
"""
function check_closed_loop_pressure(model::MultiModel, controls::AbstractVector; kwarg...)
    ok = true
    for forces in unique(objectid, controls)
        ok = check_closed_loop_pressure(model, forces; kwarg...) && ok
    end
    return ok
end

function check_closed_loop_pressure(model::MultiModel, controls; kwarg...)
    if haskey(controls, :Facility)
        controls = controls[:Facility].control
    end
    ok = true
    for (key, submodel) in pairs(model.models)
        well = submodel.data_domain
        is_closed_loop_well(well) || continue
        ok = check_closed_loop_pressure(submodel, controls; kwarg...) && ok
    end
    return ok
end

function check_closed_loop_pressure(well_model::SimulationModel, controls;
    density = missing,
    p_min = missing)

    well = well_model.data_domain
    W = physical_representation(well)
    if ismissing(density)
        density = maximum(JutulDarcy.reference_densities(well_model.system))
    end
    if ismissing(p_min)
        p_min = JutulDarcy.DEFAULT_MINIMUM_PRESSURE
        pvar = get(well_model.primary_variables, :Pressure, nothing)
        if pvar isa JutulDarcy.Pressure
            p_min = pvar.minimum_pressure
        end
    end
    # Depth of the highest node in the loop
    z = well[:cell_centroids][3, :]
    z_top = minimum(z)
    ok = true
    for (ix, name) in enumerate(W.multiwell.names)
        ctrl = get(controls, name, nothing)
        if !(ctrl isa Union{InjectorControl, ProducerControl}) ||
            !(ctrl.target isa BottomHolePressureTarget)
            continue
        end
        bhp = ctrl.target.value
        Δp = density*Jutul.gravity_constant*(z[W.multiwell.top_nodes[ix]] - z_top)
        if bhp - Δp < p_min
            ok = false
            @warn "BHP of closed-loop well $name is too low to keep the pressure \
                in the loop above the minimum pressure \
                ($(round(p_min/1e5, digits = 3)) bar). The hydrostatic pressure \
                at the top of the loop is estimated to \
                $(round((bhp - Δp)/1e5, digits = 3)) bar, and the BHP should be \
                at least $(round((p_min + Δp)/1e5, digits = 3)) bar."
        end
    end
    return ok
end

## Merging of supply and return wells

"""
    merge_closed_loop_well(supply_well::DataDomain, return_reservoir_cell = missing;
        name, return_name = Symbol(name, "_return"))

Merge the supply part of a closed-loop well with a return well consisting of a
single node at the end node of the supply well. The return node gets a copy of
the data in the supply end node, and a zero-length segment connects the two.
The return well is perforated in `return_reservoir_cell` with zero well
indices (see [`setup_closed_loop_well`](@ref) for the default).
"""
function merge_closed_loop_well(supply_well::DataDomain, return_reservoir_cell = missing;
    name, return_name = Symbol(name, "_return"))

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
    if ismissing(return_reservoir_cell)
        # Use the reservoir cell at the end node if it is perforated
        perf = findfirst(isequal(end_node), w.perforations.self)
        if isnothing(perf)
            perf = 1
        end
        return_reservoir_cell = w.perforations.reservoir[perf]
    else
        perf = findfirst(isequal(return_reservoir_cell), w.perforations.reservoir)
        if isnothing(perf)
            perf = 1
        end
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
