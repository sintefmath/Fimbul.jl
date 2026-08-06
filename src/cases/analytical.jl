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

"""
    analytical_ates(; <keyword arguments>)

Create the radial ATES benchmark case used for comparison against the Gelhar &
Collins analytical heat-transport solution.

# Keyword arguments
- `nang = 32`: Number of angular sectors in the radial mesh.
- `radii = vcat(collect(1.0:100.0), 100.0 .+ cumsum([1.0 * 1.3^(i - 1) for i in 1:14]))`:
  Radial cell-edge positions in meters.
- `layer_depths = [0.0, 5.0, 10.0, 15.0, 20.0]`: Vertical layer interfaces in meters.
- `permeability_horizontal = 2.31e-11*meter^2`: Horizontal permeability.
- `permeability_vertical = permeability_horizontal/3`: Vertical permeability.
- `porosity = 0.35`: Rock porosity.
- `rock_thermal_conductivity = 2.0*watt/(meter*Kelvin)`: Rock thermal conductivity.
- `fluid_thermal_conductivity = 0.59*watt/(meter*Kelvin)`: Fluid thermal conductivity.
- `rock_heat_capacity = 800.0*joule/(kilogram*Kelvin)`: Rock heat capacity.
- `rock_density = 2650.0*kilogram/meter^3`: Rock density.
- `component_heat_capacity = 4184.0*joule/(kilogram*Kelvin)`: Fluid heat capacity.
- `boundary_pressure = 1.0*atm`: Outer boundary pressure.
- `boundary_temperature = convert_to_si(10.0, :Celsius)`: Outer boundary temperature.
- `injection_rate = 400.0*meter^3/day`: Total injection rate.
- `injection_temperature = convert_to_si(20.0, :Celsius)`: Injection temperature.
- `well_radius = 0.1`: Well radius.
- `num_steps = 90`: Number of report steps.
- `step_size = day`: Length of each report step.

# Returns
A `JutulCase` for the analytical radial ATES benchmark.
"""
function analytical_ates(;
    nang::Int = 32,
    nrad::Int = 100,
    # radii = vcat(collect(0.5:99.5), 100.0 .+ cumsum([1.0 * 1.3^(i - 1) for i in 1:14])),
    layer_depths = [0.0, 5.0, 10.0, 15.0, 20.0],
    permeability_horizontal = 2.31e-11*meter^2,
    permeability_vertical = permeability_horizontal/3,
    porosity = 0.35,
    rock_thermal_conductivity = 2.0*watt/(meter*Kelvin),
    thermal_dispersivity = 0.0*watt/(meter*Kelvin),
    fluid_thermal_conductivity = 0.59*watt/(meter*Kelvin),
    rock_heat_capacity = 800.0*joule/(kilogram*Kelvin),
    rock_density = 2650.0*kilogram/meter^3,
    component_heat_capacity = 4184.0*joule/(kilogram*Kelvin),
    boundary_pressure = 1.0*atm,
    boundary_temperature = convert_to_si(10.0, :Celsius),
    injection_rate = 400.0*meter^3/day,
    injection_temperature = convert_to_si(20.0, :Celsius),
    well_radius = 0.1,
    num_steps::Int = 90,
    step_size = day,
)

    radii = vcat(range(0.5, 99.5, length = nrad), 100.0 .+ cumsum([1.0 * 1.3^(i - 1) for i in 1:14]))
    # @assert length(radii) == 250 + 14 "Expected 250 + 14 radial cells for the benchmark discretization."
    # @assert length(layer_depths) == 5 "Expected 4 vertical layers for the benchmark discretization."

    mesh_2d = Jutul.RadialMeshes.radial_mesh(nang, radii; centerpoint = false)
    mesh = Jutul.extrude_mesh(mesh_2d, layer_depths)

    domain = reservoir_domain(mesh;
        permeability = [permeability_horizontal, permeability_horizontal, permeability_vertical],
        porosity = porosity,
        # rock_thermal_conductivity = rock_thermal_conductivity .+ thermal_dispersivity*porosity/(1-porosity),
        rock_thermal_conductivity = rock_thermal_conductivity,
        fluid_thermal_conductivity = fluid_thermal_conductivity,
        rock_heat_capacity = rock_heat_capacity,
        rock_density = rock_density,
        component_heat_capacity = component_heat_capacity,
    )

    geometry = tpfv_geometry(mesh)
    well_cells = isapprox.(geometry.cell_centroids[1,:], -0.5) .&& isapprox.(geometry.cell_centroids[2,:], -0.5)
    well_cells = findall(well_cells)
    # return mesh, geometry, well_cells

    ρ = 998.0

    well = setup_well(domain, well_cells;
        name = :Well,
        simple_well = false,
        radius = well_radius,
        WIth = 10
    )
    control = Dict(
        :Well => InjectorControl(TotalRateTarget(injection_rate), [1.0];
            density = ρ,
            temperature = injection_temperature,
            check = false,
        )
    )
    system = SinglePhaseSystem(AqueousPhase(); reference_density = ρ)
    # system = :geothermal

    model = setup_reservoir_model(domain, system; wells = [well], thermal = true)
    λ = fluid_thermal_conductivity' .+ thermal_dispersivity
    # model.models[:Reservoir].data_domain[:fluid_thermal_conductivity] = fluid_thermal_conductivity

    # ϕ = model.models[:Reservoir].data_domain[:porosity]
    # λ = model.models[:Reservoir].data_domain[:fluid_thermal_conductivity]
    T = compute_face_trans(domain, porosity.*λ)
    T = repeat(T', 1, 1)
    model.models[:Reservoir].data_domain[:fluid_thermal_conductivities, Faces()] = T


    density = ConstantCompressibilityDensities(system, si_unit(:atm), [ρ], [1e-10/si_unit(:bar)]) # Replace density with a lighter pair
    replace_variables!(model, PhaseMassDensities = density);
    
    # z_boundary = geometry.boundary_centroids[3, :]
    # top = isapprox.(z_boundary, minimum(z_boundary))
    # bottom = isapprox.(z_boundary, maximum(z_boundary))
    # sides = .!(top .| bottom)

    bc, state0, _, _ = set_dirichlet_bcs(model, :sides;
        pressure_surface = boundary_pressure,
        temperature_surface = boundary_temperature,
        geothermal_gradient = 0.0,
        output_state = true,
    )

    # bc = flow_boundary_condition(
    #     geometry.boundary_neighbors[sides],
    #     domain,
    #     boundary_pressure,
    #     boundary_temperature,
    # )

    # state0 = setup_reservoir_state(model;
    #     Pressure = boundary_pressure,
    #     Temperature = boundary_temperature,
    # )

    force = setup_reservoir_forces(model; bc = bc, control = control)
    dt = fill(step_size, num_steps)
    forces = fill(force, num_steps)

    info = Dict{Symbol, Any}(
        :description => "Analytical radial ATES benchmark case generated using Fimbul.analytical_ates()",
        :grid => (nang = nang, nr = length(radii), nz = length(layer_depths) - 1),
        :radii => radii,
        :layer_depths => layer_depths,
        :boundary_pressure => boundary_pressure,
        :boundary_temperature => boundary_temperature,
        :injection_rate => injection_rate,
        :injection_temperature => injection_temperature,
        :well_name => :Well,
    )

    # parameters = setup_parameters(model)
    # parameters[:Reservoir][:TwoPointGravityDifference] .= 0.0
    parameters = nothing
    return JutulCase(model, dt, forces; state0 = state0, input_data = info, parameters = parameters)
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
