meter, = si_units(:meter)
atm, = si_units(:atm)
kilogram = si_unit(:kilogram)
second, hour, year = si_units(:second, :hour, :year)
Kelvin, joule, watt, = si_units(:Kelvin, :joule, :watt)

function analytical_1d(;
    num_cells = 1000,
    length_x = 100.0,
    thermal_conductivity = 2.0watt/(meter*Kelvin),
    heat_capacity = 900.0joule/(kilogram*Kelvin),
    density = 2600kilogram/meter^3,
    temperature_boundary = convert_to_si(10.0, :Celsius),
    initial_condition = missing,
    num_steps = 100,
    )

    if ismissing(initial_condition)
        Tb = temperature_boundary
        T0 = convert_to_si(90.0, :Celsius)
        L = length_x
        initial_condition = x -> (T0 - Tb)*sin(π*x/L)
    end

    dx = length_x/num_cells
    dx = 0.0
    mesh = CartesianMesh((num_cells, 1, 1), (length_x + dx, 1.0, 1.0), origin = (-dx/2, 0.0, 0.0))
    domain = reservoir_domain(mesh;
        porosity = 1e-10,
        permeability = 1e-6si_unit(:darcy),
        rock_thermal_conductivity = thermal_conductivity,
        rock_heat_capacity = heat_capacity,
        rock_density = density,
        fluid_thermal_conductivity = thermal_conductivity,
        component_heat_capacity = heat_capacity,
    )

    sys = SinglePhaseSystem(AqueousPhase(); reference_density = density)
    model, parameters = setup_reservoir_model(
        domain, sys,
        thermal = true,
    );
    
    geo = tpfv_geometry(mesh)
   
    x = geo.cell_centroids[1,:]

    sol = (x,t) -> analytical_solution_1d(x, t,
    length_x, thermal_conductivity, heat_capacity, density,
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


    L, λ, Cₚ, ρ = length_x, thermal_conductivity, heat_capacity, density
    α = λ/(ρ*Cₚ)
    println("α = $α")
    time = -log(0.1)/(α*(π/L)^2)
    dt = fill(time/num_steps, num_steps)

    case = JutulCase(model, dt, forces, state0 = state0, parameters = parameters)

    return case, sol

end

function analytical_solution_1d(x, t,
        length_x, 
        thermal_conductivity,
        heat_capacity,
        density,
        initial_condition,
        temperature_boundary,
        k_cutoff = 200
    )

    L, λ, Cₚ, ρ = length_x, thermal_conductivity, heat_capacity, density
    α = λ/(ρ*Cₚ) # Thermal diffusivity

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