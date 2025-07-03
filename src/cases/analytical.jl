# Useful SI units
meter, = si_units(:meter)
atm, = si_units(:atm)
kilogram = si_unit(:kilogram)
second, hour, year = si_units(:second, :hour, :year)
Kelvin, joule, watt, = si_units(:Kelvin, :joule, :watt)

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
    model, parameters = setup_reservoir_model(
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
    case = JutulCase(model, dt, forces, state0 = state0, parameters = parameters)

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