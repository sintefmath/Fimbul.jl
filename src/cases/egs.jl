import Jutul.CutCellMeshes: PlaneCut, PolygonalSurface, cut_mesh

"""
    egs_well_coordinates(; kwargs...)

Generate default EGS well trajectories with a smooth curved bend from the
vertical to the horizontal section. Returns `(injector_coords, producer_coords)`
as vectors of `n×3` trajectory matrices (each row is an `(x, y, z)` waypoint).

# Keyword arguments
- `well_depth = 2500.0`: Vertical depth of the injector horizontal section [m]
- `well_spacing_x = 100.0`: Horizontal separation between injector and
  producers in x [m]
- `well_lateral = 1000.0`: Length of horizontal well section [m]
- `bend_radius = 200.0`: Radius of curvature at the vertical-to-horizontal
  bend [m]
- `n_bend = 8`: Number of waypoints along the curved bend
- `producer_depth_offset = nothing`: Vertical offset of producer horizontal
  sections above the injector [m]. Defaults to the height of an equilateral
  triangle with side `well_spacing_x`.
"""
function egs_well_coordinates(;
    well_depth = 2500.0,
    well_spacing_x = 100.0,
    well_lateral = 1000.0,
    bend_radius = 200.0,
    n_bend = 8,
    producer_depth_offset = nothing,
)
    if isnothing(producer_depth_offset)
        producer_depth_offset = sqrt(3.0/4.0) * well_spacing_x
    end
    r = bend_radius
    function make_trajectory(x, wd)
        # Quarter-arc bend in the y-z plane:
        #   θ = 0   → top of arc: position (x, 0,   wd),     going downward (+z)
        #   θ = π/2 → end of arc: position (x, r,   wd),     going horizontal (+y)
        # Arc center at (x, r, wd) so the well curves from vertical to horizontal
        θ_arr = range(0.0, π/2, length = n_bend + 2)   # 0 → π/2
        y_bend = r .- r .* cos.(θ_arr)                  # 0 → r
        z_bend = wd .+ r .* (sin.(θ_arr) .- 1)                 # wd-r → wd
        traj = [x 0.0 0.0;
                x 0.0 (wd - r)]                         # surface → start of bend
        for i in 2:length(θ_arr)-1
            traj = vcat(traj, [x y_bend[i] z_bend[i]])  # interior bend waypoints
        end
        traj = vcat(traj, [x well_lateral wd])           # end of horizontal section
        return traj
    end
    injector_coords = [make_trajectory(0.0, well_depth)]
    producer_coords = [
        make_trajectory(-well_spacing_x/2, well_depth - producer_depth_offset),
        make_trajectory( well_spacing_x/2, well_depth - producer_depth_offset),
    ]
    return injector_coords, producer_coords
end

"""
    egs(injector_coords, producer_coords, fracture_radius, fracture_spacing; kwargs...)

Set up an Enhanced Geothermal System (EGS) simulation case using a Discrete
Fracture Model (DFM) for the fracture network.

Each well group (injectors, producers) is represented as a single multi-segment
well with branches coupled at a shared top node.

# Arguments
- `injector_coords`: Vector of `n×3` trajectory matrices for injector legs.
- `producer_coords`: Vector of `n×3` trajectory matrices for producer legs.
- `fracture_radius`: Radius of the stimulated fracture disks (in the x-z plane) [m].
- `fracture_spacing`: Spacing between fractures along the horizontal well [m].

# Keyword arguments
- `fracture_aperture = 0.5e-3`: Physical fracture aperture [m]
- `fracture_porosity = 0.5`: Fracture porosity [-]
- `fracture_permeability = missing`: Fracture permeability (defaults to `a²/12`)
- `permeability = 1e-4*darcy`: Matrix permeability
- `porosity = 0.01`: Matrix porosity
- `rock_thermal_conductivity = 2.5*watt/(meter*Kelvin)`: Rock thermal conductivity
- `rock_heat_capacity = 900*joule/(kilogram*Kelvin)`: Rock heat capacity
- `temperature_inj = 25°C`: Injection temperature
- `rate`: Total injection rate
- `num_years = 20`: Number of years to simulate
- `fracture_offset = 50.0`: Minimum distance along the horizontal section before
  the first fracture [m], ensuring the well is fully horizontal before fracturing
- `fracture_theta = 0.0`: Tilt angle(s) between the fracture normal and the well
  axis (y-direction). A scalar applies the same angle to all fractures; a vector
  of length `n_frac` sets one angle per fracture; a 2-tuple `(theta, theta_std)`
  samples angles from `N(theta, theta_std)` independently per fracture [rad]
- `mesh_args = NamedTuple()`: Extra keyword arguments forwarded to `extruded_mesh`
- `schedule_args = NamedTuple()`: Extra keyword arguments forwarded to
  `make_schedule` (e.g. `report_interval`)
"""
function egs(
    injector_coords::Vector{<:Matrix{Float64}},
    producer_coords::Vector{<:Matrix{Float64}},
    fracture_radius::Real,
    fracture_spacing::Real;
    fracture_aperture = 0.5e-3,
    fracture_porosity = 0.5,
    fracture_permeability = missing,
    permeability = 1e-4 * darcy,
    porosity = 0.01,
    rock_thermal_conductivity = 2.5 * watt/(meter*Kelvin),
    rock_heat_capacity = 900.0 * joule/(kilogram*Kelvin),
    temperature_inj = convert_to_si(25, :Celsius),
    rate = 100kilogram/second/(1000kilogram/meter^3),
    num_years = 20,
    fracture_offset = 50.0,
    fracture_theta = 0.0,
    mesh_args = NamedTuple(),
    schedule_args = NamedTuple(),
)
    all_coords = vcat(injector_coords, producer_coords)

    # ── Key depths ────────────────────────────────────────────────────────────
    z_wells  = [maximum(traj[:, 3]) for traj in all_coords]
    wd_min, wd_max = extrema(z_wells)

    # Identify horizontal-section y-range (waypoints within fracture_radius of max depth)
    y_deep = vcat([traj[traj[:, 3] .>= wd_max - fracture_radius, 2]
                   for traj in all_coords]...)
    isempty(y_deep) && (y_deep = vcat([traj[:, 2] for traj in all_coords]...))
    y_horiz_min = minimum(y_deep)
    y_horiz_max = maximum(y_deep)

    # Fracture y-positions, evenly distributed along the horizontal section.
    # fracture_offset ensures the well is fully horizontal before the first fracture.
    frac_y_start = y_horiz_min + fracture_offset + fracture_spacing / 2
    frac_y_end   = y_horiz_max - fracture_spacing / 2
    n_frac = max(1, Int(round((frac_y_end - frac_y_start) / fracture_spacing)) + 1)
    y_fracs = collect(range(frac_y_start, frac_y_end, length = n_frac))

    # Depth layers: overburden | reservoir zone | buffer
    wd_top = max(0.0, wd_min - max(fracture_radius * 2.0, wd_min * 0.15))
    depths = [0.0,
              wd_top,
              wd_max + fracture_radius * 0.5,
              wd_max + fracture_radius * 2.0]
    hz = diff(depths) ./ [2.0, n_frac * 3.0, 2.0]

    # ── 2D cell constraints (well x-y footprints) ─────────────────────────────
    cell_constraints = Matrix{Float64}[]
    for traj in all_coords
        push!(cell_constraints, traj[:, 1:2]')  # 2×n_points
    end

    x_wells = [mean(traj[:, 1]) for traj in all_coords]
    Δx_min  = length(x_wells) > 1 ? minimum(abs.(diff(sort(x_wells)))) : 2.0 * fracture_radius
    hxy_min = min(Δx_min / 3.0, fracture_radius / 4.0)

    # Shift cell constraints by hxy_min/2 so well coordinates don't land on cell boundaries
    cell_constraints = [c .+ hxy_min/2 for c in cell_constraints]

    # ── Build extruded matrix mesh ─────────────────────────────────────────────
    msh, layers, _ = extruded_mesh(cell_constraints, depths;
        hxy_min  = hxy_min,
        hz       = hz,
        offset_rel = 1.0,
        mesh_args...
    )

    # ── Add DFM fractures ────────────────────────────────────────────────────
    # Each fracture is a disk centered on the well at the fracture y-position.
    # fracture_theta controls the angle between the fracture normal and the well
    # direction (y-axis): theta=0 means normal parallel to well (fracture ⊥ well).
    # fracture_theta may be:
    #   - a scalar (same angle for all fractures)
    #   - a vector of length n_frac (one angle per fracture)
    #   - a 2-tuple (theta_mean, theta_std) → sample N(theta_mean, theta_std)
    x_c    = mean(x_wells)
    z_c    = mean(z_wells)
    n_poly = 16
    θ_poly = range(0.0, 2π, length = n_poly + 1)[1:n_poly]

    # Resolve per-fracture tilt angles
    _n_frac = length(y_fracs)
    if fracture_theta isa Tuple{<:Real, <:Real}
        θ_mean, θ_std = fracture_theta
        tilt_angles = θ_mean .+ θ_std .* randn(_n_frac)
    elseif fracture_theta isa AbstractVector
        length(fracture_theta) == _n_frac ||
            error("fracture_theta vector must have length n_frac = $_n_frac")
        tilt_angles = Float64.(fracture_theta)
    else
        tilt_angles = fill(Float64(fracture_theta), _n_frac)
    end

    cuts = PolygonalSurface[]
    cuts = PlaneCut[]
    for (fno, y_frac) in enumerate(y_fracs)
        α = tilt_angles[fno]   # tilt around x-axis: 0 → disk in x-z plane, normal=[0,1,0]
        # Disk parameterisation: local u = x-axis, v = z-axis rotated by α around x
        u_vec = Jutul.SVector(1.0, 0.0, 0.0)
        v_vec = Jutul.SVector(0.0, -sin(α), cos(α))  # tilt in y-z plane
        # polygon = [Jutul.SVector{3, Float64}(
        #     x_c + fracture_radius * (cos(θ) * u_vec[1] + sin(θ) * v_vec[1]),
        #     y_frac + fracture_radius * (cos(θ) * u_vec[2] + sin(θ) * v_vec[2]),
        #     z_c + fracture_radius * (cos(θ) * u_vec[3] + sin(θ) * v_vec[3])) for θ in θ_poly]
        # push!(cuts, PolygonalSurface([polygon]))
        push!(cuts, PlaneCut((x_c, y_frac, z_c), (0.0, cos(α), sin(α))))
    end

    msh, info = cut_mesh(msh, cuts; extra_out = true, min_cut_fraction = 0.0)
    fracture_faces = findall(info[:face_index] .== 0)
    layers = layers[info[:cell_index]]
    # Filter out fracture faces outside the fracture radius
    geo = tpfv_geometry(msh)
    frac_centroids = geo.face_centroids[:, fracture_faces]
    dist_to_well = sqrt.((frac_centroids[1, :] .- x_c).^2 .+ (frac_centroids[3, :] .- z_c).^2)
    valid_fracture_faces = dist_to_well .<= fracture_radius * 1.5
    fracture_faces = fracture_faces[valid_fracture_faces]

    fracture_mesh = Jutul.EmbeddedMeshes.EmbeddedMesh(msh, fracture_faces)
    geo = tpfv_geometry(msh)

    # return fracture_mesh, msh

    # ── Matrix domain ──────────────────────────────────────────────────────────
    matrix_domain = layered_reservoir_domain(msh, layers,
        (
            permeability              = [permeability],
            porosity                  = [porosity],
            rock_thermal_conductivity = [rock_thermal_conductivity],
            rock_heat_capacity        = [rock_heat_capacity],
        )
    )

    # ── Fracture domain ────────────────────────────────────────────────────────
    if ismissing(fracture_permeability)
        fracture_permeability = fracture_aperture^2 / 12
    end
    frac_domain = JutulDarcy.fracture_domain(fracture_mesh, matrix_domain;
        aperture     = fracture_aperture,
        porosity     = fracture_porosity,
        permeability = fracture_permeability,
    )

    # ── Injector well (multi-branch MSW, all legs from a shared top node) ─────
    inj_connectivity = zeros(Int, length(injector_coords) + 1, 2)
    inj_connectivity[2:end, 1] .= 1
    inj_cells, inj_wcells, inj_N, _, inj_dir = Fimbul.get_well_neighborship(
        msh, injector_coords, inj_connectivity, geo; top_node = true, output_directions=true, n = 1_000_000)
    # Change to :x/:y/:z as vector directions are not supported for fractures yet
    inj_dir_sym = Vector{Symbol}(undef, length(inj_dir))
    for (i, d) in enumerate(inj_dir)
        inj_dir_sym[i] = [:x, :y, :z][last(findmax(abs.(d)))]
    end
    # TODO: This is a hack to cater for how fractures are added to wells, and
    # should be removed once this is fixed in setup_fractured_reservoir_model 
    inj_dir_sym[1] = :y
    inj_cell_centers = hcat(zeros(3), geo.cell_centroids[:, inj_cells])
    well_inj = setup_well(matrix_domain, inj_cells;
        name                   = :Injector,
        neighborship           = inj_N,
        perforation_cells_well = inj_wcells[2:end],
        well_cell_centers      = inj_cell_centers,
        dir                    = inj_dir_sym,
        use_top_node           = true,
        simple_well            = false,
    )

    # ── Producer well (multi-branch MSW, all legs from a shared top node) ─────
    prod_connectivity = zeros(Int, length(producer_coords) + 1, 2)
    prod_connectivity[2:end, 1] .= 1
    prod_cells, prod_wcells, prod_N, _, prod_dir = Fimbul.get_well_neighborship(
        msh, producer_coords, prod_connectivity, geo; top_node = true, output_directions=true, n = 1_000_000)
    prod_dir_sym = Vector{Symbol}(undef, length(prod_dir))
    for (i, d) in enumerate(prod_dir)
        prod_dir_sym[i] = [:x, :y, :z][last(findmax(abs.(d)))]
    end
    # TODO: This is a hack to cater for how fractures are added to wells, and
    # should be removed once this is fixed in setup_fractured_reservoir_model
    prod_dir_sym[1] = :y
    prod_cell_centers = hcat(zeros(3), geo.cell_centroids[:, prod_cells])
    well_prod = setup_well(matrix_domain, prod_cells;
        name                   = :Producer,
        neighborship           = prod_N,
        perforation_cells_well = prod_wcells[2:end],
        well_cell_centers      = prod_cell_centers,
        dir                    = prod_dir_sym,
        use_top_node           = true,
        simple_well            = false,
    )

    wells = [well_inj, well_prod]

    # ── DFM model ─────────────────────────────────────────────────────────────
    model = JutulDarcy.setup_fractured_reservoir_model(matrix_domain, frac_domain, :geothermal;
        wells = wells, block_backend = true)

    for (well_name, well_model) in get_model_wells(model; data_domain = true)
        adjust_well_indices!(well_model, well_name, true)
    end

    bc, p, T = set_dirichlet_bcs(model; pressure_surface = 100atm, output_state=false)
    z = geo.cell_centroids[3, :]
    z_min = minimum(z)
    z_hat = z .- z_min
    state0 = setup_reservoir_state(model; Pressure = p(0.0), Temperature = T(0.0))
    state0[:Reservoir][:Pressure] .= p(z)
    state0[:Reservoir][:Temperature] .= T(z)

    geo_frac = tpfv_geometry(fracture_mesh)
    z = geo_frac.cell_centroids[3, :]
    z_hat = z .- z_min
    state0[:Fractures][:Pressure] .= p(z)
    state0[:Fractures][:Temperature] .= T(z)

    z = model.models[:Injector].data_domain[:cell_centroids][3, :]
    z_hat = z .- z_min
    println(size(state0[:Injector][:Pressure]), " injector cells")
    println(size(z_hat), " injector cell centroids")
    state0[:Injector][:Pressure] .= p(z)
    state0[:Injector][:Temperature] .= T(z)

    z = model.models[:Producer].data_domain[:cell_centroids][3, :]
    z_hat = z .- z_min
    state0[:Producer][:Pressure] .= p(z)
    state0[:Producer][:Temperature] .= T(z)

    rho = reservoir_model(model).system.rho_ref[1]
    ctrl_inj = InjectorControl(TotalRateTarget(rate), [1.0];
        density = rho, temperature = temperature_inj, check = false)
    ctrl_prod = ProducerControl(BottomHolePressureTarget(50si"atm"))
    control = Dict(:Injector => ctrl_inj, :Producer => ctrl_prod)

    forces = setup_reservoir_forces(model, control = control, bc = bc)
    dt, forces = make_schedule([forces], [(1, 1), (1, 1)];
        num_years = num_years, schedule_args...)

    input_data = Dict{Symbol, Any}(
        :injector_coords => injector_coords,
        :producer_coords => producer_coords,
        :fracture_radius => fracture_radius,
        :fracture_spacing => fracture_spacing,
        :fracture_offset  => fracture_offset,
        :fracture_theta   => fracture_theta,
        :tilt_angles      => tilt_angles,
        :y_fracs => y_fracs,
        :cut_no  => info[:cut_no],
    )

    return JutulCase(model, dt, forces; state0 = state0, input_data = input_data)
end

"""
    egs(well_coords, fracture_radius, fracture_spacing; kwargs...)

Backward-compatible dispatch: `well_coords[1]` is the injector trajectory,
`well_coords[2:end]` are the producer trajectories.
"""
function egs(well_coords::Vector{<:Matrix}, fracture_radius, fracture_spacing; kwargs...)
    egs([well_coords[1]], well_coords[2:end], fracture_radius, fracture_spacing; kwargs...)
end

"""
    get_egs_fracture_data(states, model)

Extract per-fracture temperature and thermal-energy data from EGS DFM
simulation results. Fracture cells are grouped by y-coordinate (fracture
position) and mean temperature / total thermal energy are returned for every
timestep.

Returns a `Dict` with:
- `:y`:                 y-coordinates of each fracture [m]
- `:Temperature`:       `(n_steps × n_frac)` matrix of mean fracture temperature [K]
- `:TotalThermalEnergy`: `(n_steps × n_frac)` matrix of total thermal energy [J]
"""
function get_egs_fracture_data(states, model)
    frac_domain = model.models[:Fractures].data_domain
    frac_mesh   = physical_representation(frac_domain)
    frac_geo    = tpfv_geometry(frac_mesh)
    y_all = frac_geo.cell_centroids[2, :]

    # Group fracture cells by y-position (round to nearest metre)
    y_rounded = round.(y_all; digits = 0)
    y_unique  = sort(unique(y_rounded))
    n_frac    = length(y_unique)
    n_step    = length(states)

    T_mat = zeros(n_step, n_frac)
    E_mat = zeros(n_step, n_frac)

    for (sno, state) in enumerate(states)
        T_frac = state[:Fractures][:Temperature]
        E_frac = state[:Fractures][:TotalThermalEnergy]
        for (fno, y) in enumerate(y_unique)
            ix = y_rounded .== y
            T_mat[sno, fno] = mean(T_frac[ix])
            E_mat[sno, fno] = sum(E_frac[ix])
        end
    end

    return Dict{Symbol, Any}(
        :y                   => y_unique,
        :Temperature         => T_mat,
        :TotalThermalEnergy  => E_mat,
    )
end