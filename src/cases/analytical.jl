"""
    analytical_radial(; <keyword arguments>)

Setup function for radial thermal advection in a homogeneous aquifer, with
piston-like displacement analytical solution.

A single injection well is placed at the center of a 3D Cartesian domain with
logarithmic radial grid refinement. Hot fluid is injected at constant rate,
creating a radially expanding thermal plume. Under the assumption of advection-
dominated transport and piston-like displacement, the thermal front radius is
given analytically by

    r_th(t) = sqrt(Cf * Q * t / (Caq * π * H))

where `Cf = fluid_density * fluid_heat_capacity` is the volumetric heat
capacity of the injected fluid, `Q` is the volumetric injection rate,
`Caq = Cf*ϕ + Cr*(1 - ϕ)` is the volumetric heat capacity of the aquifer
(fluid + rock), `ϕ` is the porosity, and `H` is the aquifer thickness.

# Keyword arguments
- `domain_radius = 500.0`: Radius of the simulation domain (m). Should be at
  least twice the expected thermal front radius.
- `aquifer_depth_top = 850.0`: Depth to the top of the aquifer (m).
- `aquifer_depth_bottom = 950.0`: Depth to the bottom of the aquifer (m).
- `porosity = 0.3`: Aquifer porosity (-).
- `permeability = 500.0e-3*darcy`: Aquifer permeability (m²).
- `rock_density = 2600.0*kilogram/meter^3`: Rock density (kg/m³).
- `rock_heat_capacity = 900.0*joule/(kilogram*Kelvin)`: Rock specific heat
  capacity (J/(kg K)).
- `rock_thermal_conductivity = 2.0*watt/(meter*Kelvin)`: Rock thermal
  conductivity (W/(m K)).
- `fluid_density = 1000.0*kilogram/meter^3`: Fluid reference density (kg/m³).
- `fluid_heat_capacity = 4184.0*joule/(kilogram*Kelvin)`: Fluid specific heat
  capacity (J/(kg K)).
- `temperature_initial = 10°C`: Initial reservoir temperature (K).
- `temperature_injection = 80°C`: Injection temperature (K).
- `rate = 100.0*meter^3/day`: Volumetric injection rate (m³/s).
- `num_cells_radial = 40`: Number of cells from the center to the domain edge
  in each horizontal direction. The mesh has `2*num_cells_radial` cells in
  both x and y with logarithmic spacing.
- `num_steps = 10`: Number of simulation time steps.
- `thermal_radius_target = missing`: Target thermal front radius at end of
  simulation (m). Determines the simulation duration. Defaults to
  `0.4*domain_radius` if not provided.

# Returns
- `case`: `JutulCase` ready for simulation.
- `r_th`: Analytical thermal radius function `r_th(t)` → radius in metres.
- `T_analytical`: Analytical temperature function `T_analytical(r, t)` → temperature in K.
- `t`: Vector of cumulative simulation times (s).
"""
function analytical_radial(;
    domain_radius = 500.0meter,
    aquifer_depth_top = 850.0meter,
    aquifer_depth_bottom = 950.0meter,
    porosity = 0.3,
    permeability = 500.0e-3*darcy,
    rock_density = 2600.0*kilogram/meter^3,
    rock_heat_capacity = 900.0*joule/(kilogram*Kelvin),
    rock_thermal_conductivity = 2.0*watt/(meter*Kelvin),
    fluid_density = 1000.0*kilogram/meter^3,
    fluid_heat_capacity = 4184.0*joule/(kilogram*Kelvin),
    temperature_initial = convert_to_si(10.0, :Celsius),
    temperature_injection = convert_to_si(80.0, :Celsius),
    rate = 100.0*meter^3/day,
    num_cells_radial = 40,
    num_steps = 10,
    thermal_radius_target = missing,
)

    # ## Derived quantities
    aquifer_thickness = aquifer_depth_bottom - aquifer_depth_top
    Cf = fluid_density * fluid_heat_capacity  # Volumetric heat capacity of fluid
    Cr = rock_density * rock_heat_capacity    # Volumetric heat capacity of rock
    Caq = Cf * porosity + Cr * (1.0 - porosity) # Effective aquifer heat capacity

    # ## Analytical solution functions
    # Thermal radius as a function of time (piston-like displacement)
    r_th = t -> sqrt(Cf * rate * t / (Caq * π * aquifer_thickness))
    # Temperature as a step function: injected T inside the front, initial T outside
    T_analytical = (r, t) -> r <= r_th(t) ? temperature_injection : temperature_initial

    # ## Determine simulation time from target thermal radius
    if ismissing(thermal_radius_target)
        thermal_radius_target = 0.4 * domain_radius
    end
    simulation_time = Caq * π * aquifer_thickness * thermal_radius_target^2 / (Cf * rate)

    # ## Build a radially-graded Cartesian mesh
    # Logarithmic spacing from a small inner radius to the domain edge,
    # mirrored symmetrically about zero, gives fine cells near the well.
    r_min = domain_radius / (10.0 * num_cells_radial)
    r_edges = vcat(0.0, exp.(range(log(r_min), log(domain_radius), length = num_cells_radial)))
    dr = diff(r_edges)
    dx = vcat(reverse(dr), dr)  # Symmetric: fine near center, coarse at boundary
    dy = dx
    dz = [aquifer_thickness]
    nx = 2 * num_cells_radial
    ny = 2 * num_cells_radial
    nz = 1
    msh = CartesianMesh((nx, ny, nz), (dx, dy, dz),
        origin = [-domain_radius, -domain_radius, aquifer_depth_top])

    # ## Set up reservoir domain
    domain = reservoir_domain(msh;
        porosity = porosity,
        permeability = permeability,
        rock_thermal_conductivity = rock_thermal_conductivity,
        rock_heat_capacity = rock_heat_capacity,
        rock_density = rock_density,
        fluid_thermal_conductivity = rock_thermal_conductivity,
        component_heat_capacity = fluid_heat_capacity,
    )

    # ## Set up injection well at the domain center
    # The center indices are (num_cells_radial, num_cells_radial)
    ix_c = num_cells_radial
    iy_c = num_cells_radial
    well = setup_vertical_well(domain, ix_c, iy_c; name = :Injector, simple_well = true)

    # ## Set up reservoir model
    sys = SinglePhaseSystem(AqueousPhase(); reference_density = fluid_density)
    model = setup_reservoir_model(domain, sys; wells = [well], thermal = true)

    # ## Set up initial state (uniform pressure and temperature)
    p0 = 10atm
    state0 = setup_reservoir_state(model;
        Pressure = p0,
        Temperature = temperature_initial,
    )

    # ## Set up boundary conditions (fixed pressure and temperature at outer boundary)
    geo = tpfv_geometry(msh)
    bc_cells = geo.boundary_neighbors
    bc = flow_boundary_condition(bc_cells, domain, p0, temperature_initial)

    # ## Set up injection control
    ctrl_inj = InjectorControl(
        TotalRateTarget(rate), [1.0];
        density = fluid_density,
        temperature = temperature_injection,
    )
    control = Dict(:Injector => ctrl_inj)
    forces = setup_reservoir_forces(model; control = control, bc = bc)

    # ## Set up time steps
    dt = fill(simulation_time / num_steps, num_steps)
    t = cumsum(dt)

    # ## Assemble simulation case
    case = JutulCase(model, dt, forces; state0 = state0)

    return case, r_th, T_analytical, t

end

"""
    analytical_1d(; <keyword arguments>)

Setup function for conductive heat transfer in 1D, with analytical solution

# Keyword arguments
- `L = 100.0`: Length of the domain (m)
- `thermal_conductivity = 2.0`: Thermal conductivity of the rock (W/(m K))
- `heat_capacity = 900.0`: Heat capacity of the rock (J/(kg K))
- `density = 2600`: Density of the rock (kg/m^3)
- `temperature_boundary = 283.15`: Temperature at the boundary (K)
- `initial_condition = missing`: Initial temperature profile. Set to sine curve if not provided
- `num_cells = 100`: Number of cells in the mesh
- `num_steps = 100`: Number of time steps

"""
function analytical_1d(;
    L = 100.0,
    thermal_conductivity = 2.0watt/(meter*Kelvin),
    heat_capacity = 900.0joule/(kilogram*Kelvin),
    density = 2600kilogram/meter^3,
    temperature_boundary = convert_to_si(10.0, :Celsius),
    initial_condition = missing,
    num_cells = 100,
    num_steps = 100,
    )

    # ## Set initial conditions if not provided
    if ismissing(initial_condition)
        T_b = temperature_boundary
        T_max = convert_to_si(90.0, :Celsius)
        initial_condition = x -> (T_max - T_b)*sin(π*x/L) .+ T_b
    end

    # ## Make 1D mesh
    mesh = CartesianMesh((num_cells, 1, 1), (L, 1.0, 1.0))
    domain = reservoir_domain(mesh;
        porosity = 1e-10,
        permeability = 1e-6si_unit(:darcy),
        rock_thermal_conductivity = thermal_conductivity,
        rock_heat_capacity = heat_capacity,
        rock_density = density,
        fluid_thermal_conductivity = thermal_conductivity,
        component_heat_capacity = heat_capacity,
    )

    # ## Make resservoir model
    sys = SinglePhaseSystem(AqueousPhase(); reference_density = density)
    model = setup_reservoir_model(
        domain, sys,
        thermal = true,
    );

    # ## Define analytical solution
    geo = tpfv_geometry(mesh)
    x = geo.cell_centroids[1,:]
    sol = (x,t) -> analytical_solution_1d(x, t,
        L, thermal_conductivity, heat_capacity, density,
        initial_condition, temperature_boundary
    )

    # ## Set up initial state
    state0 = setup_reservoir_state(model,
    Pressure = 1atm,
    Temperature = initial_condition.(x)
    )

    # ## Set up boundary conditions
    bc = flow_boundary_condition(
        [1, num_cells], domain, 1atm, temperature_boundary)
    forces = setup_reservoir_forces(model, bc = bc)

    # ## Compute time at which the temperature has decreased by 99%
    λ, Cₚ, ρ = thermal_conductivity, heat_capacity, density
    α = λ/(ρ*Cₚ)
    time = -log(0.01)/(α*(π/L)^2)
    dt = fill(time/num_steps, num_steps)
    t = cumsum(dt)

    # ## Set up case
    case = JutulCase(model, dt, forces, state0 = state0)

    return case, sol, x, t

end

function analytical_solution_1d(x, t,
        L, 
        thermal_conductivity,
        heat_capacity,
        density,
        initial_condition,
        temperature_boundary,
        k_max = 500
    )

    λ, Cₚ, ρ = thermal_conductivity, heat_capacity, density
    α = λ/(ρ*Cₚ) # Thermal diffusivity

    T_b = temperature_boundary
    T = fill(T_b, length(x))
    f = initial_condition
    domain = (0.0, L)

    for k = 1:k_max
        fk = (x,p) -> (f(x) - T_b)*sin((k*π/L)*x)
        prob = IntegralProblem(fk, domain)
        Bk = 2/L*solve(prob, QuadGKJL())
        ΔT = Bk*exp((-α*(k*π/L)^2).*t).*sin.((k*π/L).*x)
        norm(ΔT) < 1e-6 ? break : nothing
        T .+= ΔT
    end

    return T

end
