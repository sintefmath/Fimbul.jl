meter, = si_units(:meter)
atm, = si_units(:atm)
kilogram = si_unit(:kilogram)
second, hour, year = si_units(:second, :hour, :year)
Kelvin, joule, watt, = si_units(:Kelvin, :joule, :watt)

function analytical_1d(;
    num_cells = 1000,
    length = 100.0,
    thermal_conductivity = 25.0watt/(meter*Kelvin),
    heat_capacity = 800.0joule/(kilogram*Kelvin),
    temperature_boundary = convert_to_si(10.0, :Celsius),
    initial_condition = missing,
    density = 2600kilogram/meter^3,
    time = 10.0year,
    num_steps = 100,
    )

    if ismissing(initial_condition)
        Tb = temperature_boundary
        T0 = convert_to_si(90.0, :Celsius)
        L = length
        initial_condition = x -> (T0 - Tb)*sin(π*x/L)
    end

    dx = length/num_cells
    mesh = CartesianMesh((num_cells, 1, 1), (length + dx, 1.0, 1.0), origin = (-dx/2, 0.0, 0.0))
    domain = reservoir_domain(mesh;
        porosity = 1e-4,
        rock_thermal_conductivity = thermal_conductivity,
        rock_heat_capacity = heat_capacity,
        rock_density = density,
        fluid_thermal_conductivity = thermal_conductivity,
        fluid_heat_capacity = heat_capacity,
    )

    sys = SinglePhaseSystem(AqueousPhase(), 1000.0)
    model, parameters = setup_reservoir_model(
        domain, sys,
        thermal = true,
    );
    
    geo = tpfv_geometry(mesh)
   
    x = geo.cell_centroids[1,:]

    sol = (x,t) -> analytical_solution_1d(x, t,
    length, thermal_conductivity, heat_capacity, density,
    initial_condition, temperature_boundary
    )
    # ## Set up initial state
    state0 = setup_reservoir_state(model,
    Pressure = 1atm,
    Temperature = sol(x, 0),
    )

    bc = flow_boundary_condition(
        [1, num_cells], domain, 1atm, temperature_boundary)

    forces = setup_reservoir_forces(model, bc = bc)

    dt = fill(time/num_steps, num_steps)

    case = JutulCase(model, dt, forces, state0 = state0, parameters = parameters)

    return case, sol

end

function analytical_solution_1d(x, t,
        beam_length, 
        thermal_conductivity,
        heat_capacity,
        density,
        initial_condition,
        temperature_boundary,
        k_cutoff = 100
    )

    L, λ, Cₚ, ρ = beam_length, thermal_conductivity, heat_capacity, density
    α = λ/(ρ*Cₚ)

    T = fill(temperature_boundary, length(x))
    f = initial_condition
    domain = (0.0, L)

    for k = 1:k_cutoff
        fk = (x,p) -> f(x)*sin((k*π/L)*x)
        prob = IntegralProblem(fk, domain)
        Bk = 2/L*solve(prob, QuadGKJL())
        T .+= Bk*exp((-α*(k*π/L)^2).*t).*sin.((k*π/L).*x)
    end

    return T

end