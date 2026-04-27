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
- `fracture_aperture = 0.5e-3`: Physical fracture aperture [m]. A scalar applies
  the same aperture to all fractures; a vector of length `n_frac` sets one
  aperture per fracture; a 2-tuple `(aperture_mean, aperture_std)` samples from
  `N(aperture_mean, aperture_std)` independently per fracture
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
- `fracture_theta = 0.0`: Tilt angle(s) of the fracture around the z-axis (vertical).
  At θ=0 the fracture normal points along y (perpendicular to the well). A scalar
  applies the same angle to all fractures; a vector of length `n_frac` sets one angle
  per fracture; a 2-tuple `(theta, theta_std)` samples angles from
  `N(theta, theta_std)` independently per fracture [rad]
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
    x_wells = [mean(traj[:, 1]) for traj in all_coords]
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

    x_c    = mean(x_wells)
    z_c    = mean(z_wells)    

    # Depth layers: overburden | reservoir zone | buffer
    depths = [0.0,
              z_c - fracture_radius * 0.75,
              z_c + fracture_radius * 0.75,
              z_c + fracture_radius * 2.0]
    hz = diff(depths) ./ [5.0, n_frac * 3.0, 5.0]

    # println("Well depths: ", z_wells)
    # println("Depth layers: ", depths)
    # println("Layer thicknesses: ", hz)
    # return nothing

    # ── 2D cell constraints (well x-y footprints) ─────────────────────────────
    cell_constraints = Matrix{Float64}[]
    for traj in all_coords
        push!(cell_constraints, traj[:, 1:2]')  # 2×n_points
    end

    Δx_min  = length(x_wells) > 1 ? minimum(abs.(diff(sort(x_wells)))) : 2.0 * fracture_radius
    hxy_min = min(Δx_min / 3.0, fracture_radius / 4.0)

    # Add fracture footprint lines to enforce mesh refinement across the full
    # fracture diameter, not just near the wellbore.
    x_c_pre = mean(x_wells)
    for y_f in y_fracs
        push!(cell_constraints, [x_c_pre - fracture_radius  x_c_pre + fracture_radius;
                                  y_f                        y_f                       ])
    end

    # Shift cell constraints by hxy_min/2 so well coordinates don't land on cell boundaries
    cell_constraints = [c .+ hxy_min/2 for c in cell_constraints]

    @info "Building extruded mesh with"
    # ── Build extruded matrix mesh ─────────────────────────────────────────────
    msh, layers, _ = extruded_mesh(cell_constraints, depths;
        hxy_min  = hxy_min,
        hz       = hz,
        offset_rel = 1.0,
        mesh_args...
    )
    @info "Building extruded mesh completed"

    # ── Add DFM fractures ────────────────────────────────────────────────────
    @info "Generating fractures and cutting mesh..."
    # Each fracture is a disk centered on the well at the fracture y-position.
    # fracture_theta controls rotation around the vertical z-axis: theta=0 means
    # the fracture lies in the x-z plane (normal along y, perpendicular to well).
    # fracture_theta may be:
    #   - a scalar (same angle for all fractures)
    #   - a vector of length n_frac (one angle per fracture)
    #   - a 2-tuple (theta_mean, theta_std) → sample N(theta_mean, theta_std)
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

    # Resolve per-fracture apertures
    if fracture_aperture isa Tuple{<:Real, <:Real}
        a_mean, a_std = fracture_aperture
        aperture_per_frac = a_mean .+ a_std .* randn(_n_frac)
    elseif fracture_aperture isa AbstractVector
        length(fracture_aperture) == _n_frac ||
            error("fracture_aperture vector must have length n_frac = $_n_frac")
        aperture_per_frac = Float64.(fracture_aperture)
    else
        aperture_per_frac = fill(Float64(fracture_aperture), _n_frac)
    end

    cuts = PolygonalSurface[]
    # cuts = PlaneCut[]
    for (fno, y_frac) in enumerate(y_fracs)
        α = tilt_angles[fno]   # tilt around z-axis: 0 → disk in x-z plane, normal=[0,1,0]
        # Disk parameterisation: local u = x-axis rotated by α around z, v = z-axis
        u_vec = Jutul.SVector(cos(α), sin(α), 0.0)   # rotation in x-y plane
        v_vec = Jutul.SVector(0.0, 0.0, 1.0)          # z-axis stays fixed
        polygon = [Jutul.SVector{3, Float64}(
            x_c + fracture_radius * (cos(θ) * u_vec[1] + sin(θ) * v_vec[1]),
            y_frac + fracture_radius * (cos(θ) * u_vec[2] + sin(θ) * v_vec[2]),
            z_c + fracture_radius * (cos(θ) * u_vec[3] + sin(θ) * v_vec[3])) for θ in θ_poly]
        push!(cuts, PolygonalSurface([polygon]))
        # push!(cuts, PlaneCut((x_c, y_frac, z_c), (0.0, cos(α), sin(α))))
    end

    msh, info = cut_mesh(msh, cuts; extra_out = true, min_cut_fraction = 0.0)
    fracture_faces = findall(info[:face_index] .== 0)
    layers = layers[info[:cell_index]]
    @info "Generating fractures and cutting mesh completed"
    # Filter out fracture faces outside the fracture radius
    # geo = tpfv_geometry(msh)
    # frac_centroids = geo.face_centroids[:, fracture_faces]
    # dist_to_well = sqrt.((frac_centroids[1, :] .- x_c).^2 .+ (frac_centroids[3, :] .- z_c).^2)
    # valid_fracture_faces = dist_to_well .<= fracture_radius * 1.5
    # fracture_faces = fracture_faces[valid_fracture_faces]

    @info "Generating embedded mesh"
    fracture_mesh = Jutul.EmbeddedMeshes.EmbeddedMesh(msh, fracture_faces)
    geo = tpfv_geometry(msh)
    @info "Generating embedded mesh completed"

    # Map per-fracture apertures to per-cell values using cut_no from cut_mesh.
    # Intersection junction cells (beyond n_frac_regular) get the mean aperture.
    n_frac_regular = length(fracture_mesh.parent_faces)
    cut_no_regular = info[:cut_no][fracture_mesh.parent_faces]
    cell_aperture  = fill(mean(aperture_per_frac), number_of_cells(fracture_mesh))
    cell_aperture[1:n_frac_regular] .= aperture_per_frac[cut_no_regular]

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
        fracture_permeability = cell_aperture .^ 2 ./ 12
    end
    frac_domain = JutulDarcy.fracture_domain(fracture_mesh, matrix_domain;
        aperture     = cell_aperture,
        porosity     = fracture_porosity,
        permeability = fracture_permeability,
    )

    # ── Injector well (multi-branch MSW, all legs from a shared top node) ─────
    inj_connectivity = zeros(Int, length(injector_coords) + 1, 2)
    inj_connectivity[2:end, 1] .= 1
    inj_cells, inj_wcells, inj_N, _, inj_dir = Fimbul.get_well_neighborship(
        msh, injector_coords, inj_connectivity, geo; top_node = true, output_directions=true, n = 10_000)
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
        msh, producer_coords, prod_connectivity, geo; top_node = true, output_directions=true, n = 10_000)
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
    @info "Setting up DFM model..."
    model = JutulDarcy.setup_fractured_reservoir_model(matrix_domain, frac_domain, :geothermal;
        wells = wells, block_backend = true)

    for (well_name, well_model) in get_model_wells(model; data_domain = true)
        adjust_well_indices!(well_model, well_name, true)
    end
    @info "Setting up DFM model completed"

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
        :fracture_theta    => fracture_theta,
        :tilt_angles       => tilt_angles,
        :fracture_aperture => fracture_aperture,
        :aperture_per_frac => aperture_per_frac,
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
    _well_fracture_energy_exchange(states, model, well_name, cut_no_frac, n_frac)

Compute the net wellbore-to-fracture thermal energy exchange for every fracture
and every timestep, using the MSW face-based energy balance.

For each fracture-connected well cell `wc`, the net wellbore energy entering
that cell equals the energy diverted to (or received from) the connected
fracture, under the quasi-steady assumption that the well-cell stored energy
changes negligibly over a timestep.

Energy flux at face `f`: `Φ(f) = TotalMassFlux(f) × H(upstream_cell)`.  
Convention: `TotalMassFlux > 0` means flow from `N[1,f]` to `N[2,f]`.

Returns an `(n_steps × n_frac)` matrix where:
- For the injector: positive entry = energy carried by cold fluid INTO the fracture.
- For the producer: negative entry = energy flowing OUT of the fracture into the
  well (negate to obtain `q_out`).
"""
function _well_fracture_energy_exchange(states, model, well_name, cut_no_frac, n_frac)
    well_domain = physical_representation(model.models[well_name].data_domain)
    haskey(well_domain.perforations, :fracture) || return zeros(length(states), n_frac)

    fc  = vec(well_domain.perforations.fracture)       # fracture cell indices
    wc  = vec(well_domain.perforations.self_fracture)  # well cell indices

    # Fracture number for each perforation (0 if cell is an intersection node)
    fno = [k <= length(cut_no_frac) ? cut_no_frac[k] : 0 for k in fc]

    N            = well_domain.neighborship        # 2 × n_faces
    n_well_cells = well_domain.num_nodes
    n_faces      = size(N, 2)

    # Pre-build cell → faces mapping with sign:
    #   sign = +1: face enters cell  (cell == N[2,f])
    #   sign = -1: face leaves cell  (cell == N[1,f])
    cell_faces = [Int[] for _ in 1:n_well_cells]
    cell_sign  = [Int[] for _ in 1:n_well_cells]
    for f in 1:n_faces
        l, r = N[1, f], N[2, f]
        push!(cell_faces[l], f); push!(cell_sign[l], -1)
        push!(cell_faces[r], f); push!(cell_sign[r], +1)
    end

    n_step = length(states)
    q_mat  = zeros(n_step, n_frac)

    for (sno, state) in enumerate(states)
        sw  = state[well_name]
        tmf = sw[:TotalMassFlux]  # mass flux at each face [kg/s], positive = N[1]→N[2]
        H   = sw[:FluidEnthalpy]  # specific enthalpy [J/kg], size (n_phases, n_cells)

        # Thermal energy flux at each face [W]: positive = from N[1] to N[2]
        Φ = Vector{Float64}(undef, n_faces)
        for f in 1:n_faces
            upstream = tmf[f] >= 0 ? N[1, f] : N[2, f]
            Φ[f] = tmf[f] * sum(H[:, upstream])
        end

        # Accumulate net wellbore energy entering each fracture-connected cell
        for (wc_i, k) in zip(wc, fno)
            k == 0 && continue
            q_net = zero(eltype(Φ))
            for (f, s) in zip(cell_faces[wc_i], cell_sign[wc_i])
                q_net += s * Φ[f]
            end
            q_mat[sno, k] += q_net
        end
    end

    return q_mat
end

"""
    get_egs_fracture_data(states, case::JutulCase)

Extract per-fracture temperature and thermal-energy data from EGS DFM
simulation results. Fracture cells are assigned to their originating fracture
plane using the `:cut_no` field stored in `case.input_data`, which is produced
by `cut_mesh` when given a vector of cuts. This is more robust than grouping by
rounded y-coordinate.

Well-fracture energy exchange is also computed from the MSW wellbore energy
balance. For each fracture, `q_in` is the thermal energy carried by the
injected fluid into the fracture and `q_out` is the thermal energy carried by
the produced fluid out of the fracture. The matrix-to-fracture heat extraction
rate (geothermal power) is then `power = dE/dt + q_out - q_in`.

Returns a `Dict` with:
- `:y`:                  y-coordinates of each fracture plane [m]
- `:Temperature`:        `(n_steps × n_frac)` matrix of mean fracture temperature [K]
- `:TotalThermalEnergy`: `(n_steps × n_frac)` matrix of total thermal energy [J]
- `:q_in`:               `(n_steps × n_frac)` matrix of injector→fracture energy flux [W]
- `:q_out`:              `(n_steps × n_frac)` matrix of fracture→producer energy flux [W]
"""
function get_egs_fracture_data(states, case::JutulCase)
    model       = case.model
    cut_no_all  = case.input_data[:cut_no]
    y_fracs     = case.input_data[:y_fracs]

    frac_domain = model.models[:Fractures].data_domain
    frac_mesh   = physical_representation(frac_domain)

    # Each regular fracture cell (1:n_regular) maps to a parent face in the cut
    # mesh. Look up which cut (fracture plane) produced that face.
    n_regular = length(frac_mesh.parent_faces)
    cut_no    = cut_no_all[frac_mesh.parent_faces]

    # Full mapping from any fracture cell index → fracture plane number
    # (intersection junction cells, if any, get 0 and are ignored)
    n_total      = number_of_cells(frac_mesh)
    cut_no_frac  = zeros(Int, n_total)
    cut_no_frac[1:n_regular] .= cut_no

    n_frac = length(y_fracs)
    n_step = length(states)
    T_mat  = zeros(n_step, n_frac)
    E_mat  = zeros(n_step, n_frac)

    for (sno, state) in enumerate(states)
        # Index only regular cells; any intersection junction cells (numbered
        # after n_regular) are excluded as their volume is negligible.
        T_frac = state[:Fractures][:Temperature][1:n_regular]
        E_frac = state[:Fractures][:TotalThermalEnergy][1:n_regular]
        for fno in 1:n_frac
            ix = cut_no .== fno
            any(ix) || continue
            T_mat[sno, fno] = mean(T_frac[ix])
            E_mat[sno, fno] = sum(E_frac[ix])
        end
    end

    # q_in:  net wellbore energy entering fracture k from the injector [W]
    #        (positive = cold fluid injected into fracture)
    # q_out: net wellbore energy leaving fracture k to the producer [W]
    #        (positive = warm fluid extracted from fracture)
    #        = minus the "net entering" value for the producer cell, because
    #          the producer cell receives fluid from the fracture but sends
    #          more up the well than it receives from below.
    q_in_mat  = _well_fracture_energy_exchange(states, model, :Injector, cut_no_frac, n_frac)
    q_out_mat = -_well_fracture_energy_exchange(states, model, :Producer, cut_no_frac, n_frac)

    return Dict{Symbol, Any}(
        :y                   => y_fracs,
        :Temperature         => T_mat,
        :TotalThermalEnergy  => E_mat,
        :q_in                => q_in_mat,
        :q_out               => q_out_mat,
    )
end

"""
    get_egs_fracture_data(states, model)

Fallback dispatch. Groups fracture cells by rounded y-coordinate.
Prefer passing the full `JutulCase` to use the exact `:cut_no` mapping.
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