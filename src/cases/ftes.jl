import Jutul.CutCellMeshes: PlaneCut, PolygonalSurface, cut_mesh

"""
    ftes(well_coordinates, fractures; <keyword arguments>)

Create a Fractured Thermal Energy Storage (FTES) simulation case from explicit
well trajectories and fracture geometry.

The FTES system consists of one injector and one multi-lateral producer well
connected through a discrete fracture network (DFM) embedded in a low-permeability
matrix. During charging, hot water is injected; during discharging the flow
direction is reversed.

# Arguments
- `well_coordinates`: Vector of `3×n` trajectory matrices (one per well). The
  first entry is the injector; the remaining entries are producer laterals. Each
  column is an `(x, y, z)` waypoint.
- `fractures`: `Dict{Symbol, Any}` with fracture geometry. Required keys:
  - `:normal` – vector of normal vectors, one per fracture
  - `:centers` – vector of center point vectors `[x, y, z]`
  - `:radius` – vector of radii (`Inf` for unbounded plane cuts)
  - `:aperture` – vector of aperture values [m]
  - `:porosity` – vector of fracture porosity values [-]

# Keyword Arguments

## Geometry
- `depths = nothing`: Layer interface depths [m]. If `nothing`, layers are
  derived automatically from the well depths. Additional required depths are
  always inserted.
- `mesh_args = NamedTuple()`: Extra keyword arguments forwarded to `extruded_mesh`.

## Rock Properties
- `matrix_properties = NamedTuple()`: Named tuple of matrix rock properties
  (e.g. `permeability`, `porosity`). Defaults to
  `(permeability = 1e-4 darcy, porosity = 0.01)`.

## Operational Parameters
- `rate_charge = missing`: Total injection rate during charging [m³/s]. If
  missing, estimated automatically from the fracture domain geometry.
- `rate_discharge = rate_charge`: Total injection rate during discharging [m³/s].
- `temperature_charge = 95°C`: Temperature of injected fluid during charging [K].
- `temperature_discharge = 20°C`: Temperature of injected fluid during discharging [K].
- `charge_period = ["April", "November"]`: Start and end month of the charging period.
- `discharge_period = ["December", "March"]`: Start and end month of the discharging period.
- `utes_schedule_args = NamedTuple()`: Extra keyword arguments forwarded to
  `make_utes_schedule` (e.g. `num_years`, `report_interval`).

## Diagnostics
- `info_level = 0`: Verbosity level. Set to 1 to print progress messages.

# Returns
A `JutulCase` for the FTES system.
"""
function ftes(well_coordinates::Vector{Matrix{Float64}}, fractures::Dict{Symbol, Any};
    depths = nothing,
    matrix_properties = NamedTuple(),
    rate_charge = missing,
    rate_discharge = missing,
    temperature_charge = convert_to_si(95.0, :Celsius),
    temperature_discharge = convert_to_si(20.0, :Celsius),
    charge_period = ["April", "November"],
    discharge_period = ["December", "March"],
    utes_schedule_args = NamedTuple(),
    mesh_args = NamedTuple(),
    info_level = 0,
    )

    # ## Build mesh constraints from well collar positions
    info_level > 0 && @info "Setting up wells and making mesh"
    collars = hcat([x[1:2, 1] for x in well_coordinates]...)
    Δx_min, Δx_max = Fimbul.min_max_distance(collars)
    hxy_min = Δx_min/4
    r_given = filter(!isinf, fractures[:radius])
    r_max = ifelse(isempty(r_given), 0.0, maximum(r_given))
    well_offset = max(Δx_max/2, r_max - Δx_max/2)
    well_outline = Fimbul.offset_boundary(collars, well_offset; n=24)
    well_outline = hcat(well_outline, well_outline[:, 1]) # Close the loop
    collars = [permutedims([x[1] x[2]]) .+ hxy_min/2 for x in eachcol(collars)]
    cell_constraints = [x for x in collars]
    push!(cell_constraints, well_outline)

    # ## Determine layer depths from well geometry
    well_depth = maximum(maximum(x[3, :]) for x in well_coordinates)
    if isnothing(depths)
        depths = [0.0, well_depth + 1e-2, well_depth*1.25]
    else
        extra_depths = [0.0, well_depth + 1e-2, well_depth*1.25]
        for d in extra_depths
            if all(!isapprox.(depths, d, atol=0.5))
                push!(depths, d)
            end
        end
        depths = sort(depths)
    end

    # ## Generate matrix mesh
    num_fractures = length(fractures[:normal])
    hz = diff(depths) ./ [num_fractures*3, 2]
    matrix_mesh, layers, _ = extruded_mesh(cell_constraints, depths;
        hxy_min=hxy_min, hz=hz, offset=well_offset*4, offset_rel=missing, mesh_args...)

    # ## Embed fractures as polygonal disk cuts
    info_level > 0 && @info "Adding fractures to the mesh"
    n_poly = 16
    θ_poly = range(0.0, 2π, length = n_poly + 1)[1:n_poly]
    function fracture_plane_basis(n_hat)
        ref = abs(n_hat[1]) < 0.9 ? [1.0, 0.0, 0.0] : [0.0, 1.0, 0.0]
        u = ref .- dot(ref, n_hat) .* n_hat; u ./= norm(u)
        v = cross(n_hat, u); v ./= norm(v)
        return u, v
    end
    cuts = []
    geo = tpfv_geometry(matrix_mesh)
    for (normal, center, r) in zip(fractures[:normal], fractures[:centers], fractures[:radius])
        if isinf(r)
            push!(cuts, PlaneCut(center, normal))
        else
            n_hat = normal ./ norm(normal)
            u, v = fracture_plane_basis(n_hat)
            polygon = [Jutul.SVector{3, Float64}(
                center[1] + r * (cos(θ) * u[1] + sin(θ) * v[1]),
                center[2] + r * (cos(θ) * u[2] + sin(θ) * v[2]),
                center[3] + r * (cos(θ) * u[3] + sin(θ) * v[3])) for θ in θ_poly]
            push!(cuts, PolygonalSurface([polygon]))
        end
    end
    matrix_mesh, cut_info = cut_mesh(matrix_mesh, cuts; extra_out = true, min_cut_fraction = 0.0)
    fracture_faces = findall(cut_info[:face_index] .== 0)
    layers = layers[cut_info[:cell_index]]
    fracture_mesh = Jutul.EmbeddedMeshes.EmbeddedMesh(matrix_mesh, fracture_faces)

    # ## Set up matrix and fracture domains
    geo = tpfv_geometry(matrix_mesh)
    if isempty(matrix_properties)
        matrix_properties = (permeability=1e-4si_unit(:darcy), porosity=0.01)
    end
    matrix_domain = layered_reservoir_domain(matrix_mesh, layers, matrix_properties)
    fracture_domain = JutulDarcy.fracture_domain(fracture_mesh, matrix_domain;
        aperture=fractures[:aperture][1],
        porosity=fractures[:porosity][1])

    # ## Set up injector well
    cells = Jutul.find_enclosing_cells(matrix_mesh, permutedims(well_coordinates[1]), n=1_000_000)
    well_inj = setup_well(matrix_domain, cells;
        name=:Injector, radius=75e-3, simple_well=false)

    # ## Set up multi-lateral producer well
    x_prod = [permutedims(x) for x in well_coordinates[2:end]]
    connectivity = zeros(Int, length(x_prod) + 1, 2)
    connectivity[2:end, 1] .= 1
    cells, wcells, neighborship = Fimbul.get_well_neighborship(
        matrix_mesh, x_prod, connectivity, geo; top_node=true, n=1_000_000)
    well_cell_centers = hcat([0; 0; 0], geo.cell_centroids[:, cells])
    well_prod = setup_well(matrix_domain, cells;
        name=:Producer,
        radius=75e-3,
        neighborship=neighborship,
        perforation_cells_well=wcells[2:end],
        well_cell_centers=well_cell_centers,
        use_top_node=true,
        simple_well=false)
    wells = [well_inj, well_prod]

    # ## Set up DFM geothermal model
    info_level > 0 && @info "Setting up DFM model"
    model = JutulDarcy.setup_fractured_reservoir_model(matrix_domain, fracture_domain, :geothermal;
        wells=wells, block_backend=true)

    # ## Set up initial state and boundary conditions
    info_level > 0 && @info "Setting up initial state, boundary conditions, and well controls"
    ρ = reservoir_model(model).system.rho_ref[1]
    p0(z) = 20si_unit(:atm)
    T0(z) = convert_to_si(10.0, :Celsius)
    z = geo.cell_centroids[3, :]
    state0 = setup_reservoir_state(model; Pressure=p0(z), Temperature=T0(z))
    bc_cells = geo.boundary_neighbors
    bc_temperature = state0[:Reservoir][:Temperature][bc_cells]
    bc_pressure = state0[:Reservoir][:Pressure][bc_cells]
    bc = flow_boundary_condition(bc_cells, matrix_domain, bc_pressure, bc_temperature)

    # ## Set up well controls for charging and discharging phases
    if ismissing(rate_charge)
        rate_charge = 20_000*Fimbul.scaled_rate(
            fracture_domain, wells, charge_period; mean_well_coordinate=true)
    end
    if ismissing(rate_discharge)
        rate_discharge = rate_charge
    end
    ctrl_prod = ProducerControl(BottomHolePressureTarget(p0(0.0)*0.1))
    ctrl_inj = InjectorControl(TotalRateTarget(rate_charge), [1.0];
        density=ρ, temperature=temperature_charge)
    control_charge = Dict()
    for w in wells
        name = w.representation.name
        control_charge[name] = (name == :Injector) ? ctrl_inj : ctrl_prod
    end
    forces_charge = setup_reservoir_forces(model, bc=bc, control=control_charge)

    ctrl_inj = InjectorControl(TotalRateTarget(rate_discharge), [1.0];
        density=ρ, temperature=temperature_discharge)
    control_discharge = Dict()
    for w in wells
        name = w.representation.name
        control_discharge[name] = (name == :Injector) ? ctrl_inj : ctrl_prod
    end
    forces_discharge = setup_reservoir_forces(model, bc=bc, control=control_discharge)
    forces_rest = setup_reservoir_forces(model, bc=bc)

    # ## Build operational schedule
    dt, forces, timestamps = make_utes_schedule(forces_charge, forces_discharge, forces_rest;
        charge_period=charge_period, discharge_period=discharge_period,
        utes_schedule_args...)

    # ## Assemble and return simulation case
    info = Dict{Symbol, Any}(
        :well_coordinates => well_coordinates,
        :fractures        => fractures,
        :timestamps       => timestamps,
    )
    case = JutulCase(model, dt, forces; state0=state0, input_data=info)
    return case

end

"""
    ftes(wells, fractures; <keyword arguments>)

Convenience constructor for a Fractured Thermal Energy Storage (FTES) simulation
case. Accepts high-level descriptions of the well layout and fracture network
and delegates to `ftes(well_coordinates, fractures_dict; kwargs...)`.

# Arguments
- `wells`: Well layout, specified as one of:
  - `Vector{Matrix{Float64}}` – explicit trajectory matrices (passed through unchanged)
  - `NamedTuple` with fields `:num_producers`, `:radius`, and `:depth` – generates a
    symmetric layout with `num_producers` producer wells arranged in a circle of
    `radius` [m] around a central injector, all reaching `depth` [m].
- `fractures`: Fracture network, specified as one of:
  - `Dict{Symbol, Any}` – explicit fracture geometry (passed through unchanged).
    Must contain the keys `:normal`, `:centers`, `:radius`, `:aperture`, and `:porosity`.
  - `Int` – number of fractures; depth range is derived automatically from the well depths.
  - `NamedTuple` with at least the fields `:num`, `:z_min`, and `:z_max`. Any
    additional fields are forwarded as keyword arguments to `setup_ftes_fractures`.

# Keyword Arguments
All keyword arguments are forwarded to `ftes(well_coordinates, fractures_dict; kwargs...)`.
See that method for full documentation.

# Returns
A `JutulCase` for the FTES system.

# Example
```julia
case = ftes(
    (num_producers = 6, radius = 30.0, depth = 250.0),
    (num = 20, z_min = 50.0, z_max = 230.0, radius = 60.0);
    rate_charge = 30si"litre/second",
    utes_schedule_args = (num_years = 3,),
)
```
"""
function ftes(wells, fractures; kwargs...)

    # ## Resolve well coordinates
    well_coordinates = wells
    if well_coordinates isa NamedTuple
        if !haskey(well_coordinates, :num_producers) || !haskey(well_coordinates, :radius) || !haskey(well_coordinates, :depth)
            error("Named tuple defining wells must contain the keys: :num_producers, :radius, and :depth")
        end
        well_coordinates = setup_ftes_well_coordinates(wells.num_producers, wells.radius, wells.depth)
    end
    well_coordinates isa Vector{Matrix{Float64}} || error("well_coordinates must be a Vector of Matrix{Float64}")

    # ## Resolve fracture geometry
    if fractures isa Int
        z_min = minimum(minimum(x[3, :]) for x in well_coordinates)
        z_max = maximum(maximum(x[3, :]) for x in well_coordinates)
        Δz = z_max - z_min
        z_min += Δz/8
        z_max -= Δz/8
        fractures = (num=fractures, z_min=z_min, z_max=z_max)
    end
    if fractures isa NamedTuple
        if !haskey(fractures, :num) || !haskey(fractures, :z_min) || !haskey(fractures, :z_max)
            error("Named tuple defining fractures must contain the keys: :num, :z_min, and :z_max")
        end
        num_fractures = fractures.num
        z_min = fractures.z_min
        z_max = fractures.z_max
        # Strip the bookkeeping keys before forwarding remaining options
        fractures = (; (k => v for (k, v) in pairs(fractures) if k ∉ (:num, :z_min, :z_max))...)
        fractures = setup_ftes_fractures(num_fractures, z_min, z_max; fractures...)
    end
    fractures isa Dict{Symbol, Any} || error("fractures must be a Dict{Symbol, Any}")
    required_keys = [:normal, :centers, :radius, :aperture, :porosity]
    for key in required_keys
        haskey(fractures, key) || error("fractures must contain the key: $key")
    end

    return ftes(well_coordinates, fractures; kwargs...)

end

"""
    setup_ftes_well_coordinates(num_producers, radius, depth)

Generate well trajectory matrices for a symmetric FTES well layout.

Places `num_producers` producer wells evenly on a circle of `radius` [m] around
a central injector. All wells are vertical and extend from the surface (z = 0)
to `depth` [m].

# Arguments
- `num_producers`: Number of producer wells.
- `radius`: Radial distance of the producers from the central injector [m].
- `depth`: Well depth [m].

# Returns
`Vector{Matrix{Float64}}` of `3×2` trajectory matrices. The first entry is the
injector; the remaining entries are the producers. Each matrix has columns
`[x, y, z]` for the wellhead (z = 0) and the toe (z = depth).
"""
function setup_ftes_well_coordinates(num_producers::Int, radius::Float64, depth::Float64)
    # Place producers evenly on a circle; first entry is the central injector
    Δθ = 2π / num_producers
    producer_coordinates = [[radius*cos(i*Δθ), radius*sin(i*Δθ), 0.0] for i in 0:num_producers-1]
    wc = pushfirst!(producer_coordinates, [0.0, 0.0, 0.0])

    well_coordinates = Vector{Matrix{Float64}}(undef, length(wc))
    for (wno, x) in enumerate(wc)
        x = repeat(x, 1, 2)
        x[3, 2] = depth
        well_coordinates[wno] = x
    end
    well_coordinates = [wc for wc in well_coordinates]

    return well_coordinates
end

"""
    setup_ftes_fractures(num, z_min, z_max; <keyword arguments>)

Generate a randomised fracture geometry dictionary for use in an FTES simulation.

Fracture centers are distributed uniformly in depth between `z_min` and `z_max`.
Strike, dip, radius, aperture, and porosity can be specified as fixed values or
as `(mean, std)` tuples, in which case they are drawn from a normal distribution.

# Arguments
- `num`: Number of fractures to generate.
- `z_min`: Minimum depth of fracture centers [m].
- `z_max`: Maximum depth of fracture centers [m].

# Keyword Arguments
- `strike = (0.0, 5.0)`: Strike angle [°], or `(mean, std)` tuple.
- `dip = (0.0, 5.0)`: Dip angle [°], or `(mean, std)` tuple.
- `radius = Inf`: Fracture radius [m] (`Inf` for unbounded plane cuts), or
  `(mean, std)` tuple.
- `aperture = 1e-4`: Fracture aperture [m], or `(mean, std)` tuple.
- `porosity = 0.5`: Fracture porosity [-], or `(mean, std)` tuple.
- `boundary_or_center = [0.0, 0.0]`: Horizontal position of fracture centers.
  Either a 2-element `[x, y]` vector (same center for all fractures) or a
  `2×n` polygon matrix defining a region from which centers are sampled
  uniformly at random.

# Returns
`Dict{Symbol, Any}` with keys `:normal`, `:centers`, `:radius`, `:aperture`,
and `:porosity`.
"""
function setup_ftes_fractures(num::Int, z_min::Float64, z_max::Float64;
    strike::Union{Float64, Tuple{Float64, Float64}}=(0.0, 5.0),
    dip::Union{Float64, Tuple{Float64, Float64}}=(0.0, 5.0),
    radius::Union{Float64, Tuple{Float64, Float64}}=Inf,
    aperture::Union{Float64, Tuple{Float64, Float64}}=1.0*1e-4,
    porosity::Union{Float64, Tuple{Float64, Float64}}=0.5,
    boundary_or_center=[0.0, 0.0],
    )

    # Promote scalar values to (mean, std) tuples with zero spread
    strike   = strike   isa Tuple{Float64, Float64} ? strike   : (strike,   0.0)
    dip      = dip      isa Tuple{Float64, Float64} ? dip      : (dip,      0.0)
    radius   = radius   isa Tuple{Float64, Float64} ? radius   : (radius,   0.0)
    aperture = aperture isa Tuple{Float64, Float64} ? aperture : (aperture, 0.0)
    porosity = porosity isa Tuple{Float64, Float64} ? porosity : (porosity, 0.0)

    # Sample strike, dip, and derived normals
    strike_vals = strike[1] .+ randn(num) .* strike[2]
    dip_vals    = dip[1]    .+ randn(num) .* dip[2]
    normal = [strike_dip_to_normal(s, d) for (s, d) in zip(strike_vals, dip_vals)]

    # Sample remaining fracture properties
    radius_vals   = radius[1]   .+ randn(num) .* radius[2]
    aperture_vals = aperture[1] .+ randn(num) .* aperture[2]
    porosity_vals = porosity[1] .+ randn(num) .* porosity[2]

    # Sample fracture center depths uniformly in [z_min, z_max]
    z = z_min .+ rand(num) .* (z_max - z_min)

    # Determine horizontal positions of fracture centers
    if boundary_or_center isa Vector
        centers = [vcat(boundary_or_center, z[i]) for i in 1:num]
    else
        x_min = minimum(boundary_or_center[1, :])
        x_max = maximum(boundary_or_center[1, :])
        y_min = minimum(boundary_or_center[2, :])
        y_max = maximum(boundary_or_center[2, :])
        centers = Vector{Vector{Float64}}()
        while length(centers) < num
            x = x_min + rand() * (x_max - x_min)
            y = y_min + rand() * (y_max - y_min)
            if Fimbul.point_in_polygon([x, y], boundary_or_center)
                push!(centers, [x, y, z[length(centers) + 1]])
            end
        end
    end

    fractures = Dict{Symbol, Any}(
        :normal   => normal,
        :centers  => centers,
        :radius   => radius_vals,
        :aperture => aperture_vals,
        :porosity => porosity_vals,
    )
    return fractures

end
