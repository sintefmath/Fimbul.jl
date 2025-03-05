meter, = si_units(:meter)
atm, = si_units(:atm)
kilogram = si_unit(:kilogram)
second, hour, year = si_units(:second, :hour, :year)
Kelvin, joule, watt, = si_units(:Kelvin, :joule, :watt)

function beam_thermal(;
    num_cells = 1000,
    length = 100.0,
    thermal_conductivity = 25.0watt/(meter*Kelvin),
    heat_capacity = 800.0joule/(kilogram*Kelvin),
    temperature_boundary = convert_to_si(10.0, :Celsius),
    temperature_max = convert_to_si(90.0, :Celsius),
    density = 2600kilogram/meter^3,
    time = 10.0year,
    num_steps = 100,
    )

    mesh = CartesianMesh((num_cells, 1, 1), (length, 1.0, 1.0))
    domain = reservoir_domain(mesh;
        porosity = 1e-3,
        rock_thermal_conductivity = thermal_conductivity,
        rock_heat_capacity = heat_capacity,
        rock_density = density,
    )
    model, parameters = setup_reservoir_model(
        domain, :geothermal,
        thermal = true,
    );
    
    geo = tpfv_geometry(mesh)
    # ϵ = 1e-3
    # parameters[:Reservoir][:FluidVolume] = ϵ*geo.volumes
    # parameters[:Reservoir][:BulkVolume] = (1-ϵ)*geo.volumes

    x = geo.cell_centroids[1,:]

    sol = (x,t) -> analytical_solution(x, t,
    length, thermal_conductivity, heat_capacity, density,
    temperature_max, temperature_boundary
    )
    # ## Set up initial state
    state0 = setup_reservoir_state(model,
    Pressure = 1atm,
    Temperature = sol(x,0)
    )

    bc = flow_boundary_condition(
        [1, num_cells], domain, 1atm, sol(x[[1,end]],0))

    forces = setup_reservoir_forces(model, bc = bc)

    dt = fill(time/num_steps, num_steps)

    case = JutulCase(model, dt, forces, state0 = state0, parameters = parameters)

    return case, sol

end

function analytical_solution(x, t,
        beam_length, 
        thermal_conductivity,
        heat_capacity,
        density,
        temperature_max,
        temperature_boundary;
        k_cutoff = 1
    )

    L, λ, Cₚ, ρ = beam_length, thermal_conductivity, heat_capacity, density

    Κ = λ/(ρ*Cₚ)

    T_max = temperature_max
    T_b = temperature_boundary
    T = fill(T_b, length(x))
    for k = 1:k_cutoff
        T .+= (T_max - T_b)*exp.((-Κ*(k*π/L)^2).*t).*sin.((k*π/L).*x)
    end

    return T

end