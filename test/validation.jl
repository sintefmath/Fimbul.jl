using Jutul, JutulDarcy, Fimbul, Test, LinearAlgebra

@testset "Validation 1D" begin

    function compute_error(res, sol, x, t)
        dx = x[2] - x[1]
        dt = t[2] - t[1]
        ϵ = 0.0
        for (k, tk) in enumerate(t)
            ϵ += dt*norm(dx*(res.states[k][:Temperature] .- sol(x, tk))/L, 2)
        end
        ϵ /= sum(t)
        return ϵ
    end

    L = 100.0
    case, sol, x, t = analytical_1d(L = L, num_cells = 100, num_steps = 100)
    res = simulate_reservoir(case, info_level = -1)
    ϵ = compute_error(res, sol, x, t)
    @test ϵ < 1e-3

    T_b = convert_to_si(10.0, :Celsius)
    T_0 = x -> convert_to_si(100.0, :Celsius).*(40 <= x < 60) .+ T_b
    case, sol, x, t = analytical_1d(L = L, num_cells = 100, num_steps = 100,
        temperature_boundary = T_b, initial_condition = T_0)
    res = simulate_reservoir(case, info_level = -1)
    ϵ = compute_error(res, sol, x, t)
    @test ϵ < 1e-2

end

@testset "Validation closed-loop" begin

    meter = si_unit(:meter)
    atm = si_unit(:atm)
    kilogram = si_unit(:kilogram)
    Kelvin, Joule = si_units(:Kelvin, :joule)
    darcy = si_unit(:darcy)
    # Operational conditions
    T_in = convert_to_si(80.0, :Celsius) # Inlet temperature
    T_rock = convert_to_si(10.0, :Celsius) # Ground temperature (assumed constant)
    Q = 20.0*meter^3/day # Volumetric flow rate
    # Fluid properties
    ρf = 988.1*kilogram/meter^3 # Fluid density
    Cpf = 4184.0*Joule/(kilogram*Kelvin) # Fluid heat capacity
    λf = 0.6405*Watt/(meter*Kelvin) # Fluid thermal conductivity
    
    L = 100.0*meter # BTES length
    geo_u1 = ( # U1 geometry
        closed_loop_type = :u1,
        radius_grout = 65e-3*meter, # Grout radius
        radius_pipe = 16e-3*meter, # Pipe outer radius
        wall_thickness_pipe = 2.9e-3*meter, # Pipe wall thickness
        pipe_spacing = 60e-3*meter, # Pipe spacing
        pipe_thermal_conductivity = λp, # Pipe thermal conductivity
    )
    geo_coax = ( # Coaxial geometry
        closed_loop_type = :coaxial,
        radius_grout = 50e-3*meter, # Grout radius
        radius_pipe_inner = 12e-3*meter, # Inner pipe outer radius
        wall_thickness_pipe_inner = 3e-3*meter, # Inner pipe wall thickness
        radius_pipe_outer = 25e-3*meter, # Outer pipe outer radius
        wall_thickness_pipe_outer = 4e-3*meter, # Outer pipe wall thickness
        inner_pipe_thermal_conductivity = λp, # Pipe thermal conductivity
        outer_pipe_thermal_conductivity = λp # Pipe thermal conductivity
    );

    function setup_btes_single( # Utility function to set up a BTES simulation case
        type; # Type of closed-loop geometry (:u1 or :coaxial)
        nz = 110, # Number of cells in vertical direction
        n_step = 1, # Number of time steps
        inlet = :outer # Inlet pipe for coaxial geometry (:outer or :inner
        )

        ## Select well arguments based on geometry type
        if type == :u1
            well_args = geo_u1
        elseif type == :coaxial
            well_args = geo_coax
        else
            error("Unknown closed loop type: $type")
        end
        ## Create reservoir domain
        dims = (1,1,nz) 
        Δ = (100000.0*dims[3]/L, 100000.0*dims[3]/L, L)
        msh = CartesianMesh(dims, Δ)
        reservoir = reservoir_domain(msh;
            rock_density = 2650.0,
            rock_heat_capacity = 900.0,
            rock_thermal_conductivity = 2.5,
            fluid_thermal_conductivity = λf,
            component_heat_capacity = Cpf,
            porosity = 0.01,
            permeability = 1e-3*darcy
        )
        ## Set up BTES well
        wells = Fimbul.setup_vertical_btes_well(reservoir, 1, 1;
            name = :BTES,
            well_args...,
            grouting_density = 2000.0,
            grouting_heat_capacity = 1500.0,
            grouting_thermal_conductivity = 2.3
        )
        ## Set up reservoir model
        sys = SinglePhaseSystem(AqueousPhase(), reference_density = 1000.0)
        model = setup_reservoir_model(reservoir, sys; wells = wells, thermal = true)
        ## Set up controls and boundary conditions
        ctrl_inj = InjectorControl( # Injection control with fixed temperature
            TotalRateTarget(Q), [1.0], 
            density=1000.0, temperature=T_in)
        ctrl_prod = ProducerControl( # Production control with fixed BHP
            BottomHolePressureTarget(5.0*atm))
        geo = tpfv_geometry(msh)
        bc_cells = geo.boundary_neighbors
        domain = reservoir_model(model).data_domain
        bc = flow_boundary_condition( # Fixed pressure and temperature BCs
            bc_cells, domain, 5.0*atm, T_rock)
        ## Set up forces and initial state
        if type == :coaxial && inlet == :outer # Coax with injection in outer pipe
            ctrl_supply = ctrl_prod
            ctrl_return = ctrl_inj
        else # Coax with injection in inner pipe or U1
            ctrl_supply = ctrl_inj
            ctrl_return = ctrl_prod
        end
        controls = Dict(
            :BTES_supply => ctrl_supply,
            :BTES_return => ctrl_return
        )
        forces = setup_reservoir_forces(model; control=controls, bc = bc)
        state0 = setup_reservoir_state( # Initialize with uniform pressure and temperature
            model; Pressure = 5*atm, Temperature = T_rock)
        ## Set time large enough to ensure steady-state
        rg = well_args.radius_grout
        time = 5/4*2*rg*(ϕ*ρf*Cpf + (1 - ϕ)*ρr*Cpr)/(ϕ*λf + (1 - ϕ)*λr)*10
        dt = fill(time/n_step, n_step)
        ## Create Jutul case
        case = JutulCase(model, dt, forces; state0 = state0)
        ## Set up simulator with temperature-based timestep selector
        sim, cfg = setup_reservoir_simulator(
            case; initial_dt = 1.0, info_level = -1);
        sel = VariableChangeTimestepSelector(:Temperature, 2.5; 
            relative = false, model = :BTES_supply)
        push!(cfg[:timestep_selectors], sel);

        return case, sim, cfg

    end

    function compute_error(case, simulated, analytical)

        well = case.model.models[:BTES_supply].data_domain
        sections = last.(well[:section]) |> unique |> collect
        sim_temp = simulated.result.states[end][:BTES_supply][:Temperature]

        ϵ = 0.0
        for section in sections
            cells = findall(last.(well[:section]) .== section)
            z = well[:cell_centroids][3, cells]
            dz = well[:cell_length][cells]
            L = sum(dz)
            Tn = sim_temp[cells]
            Ta = analytical[section].(z)
            ϵ += sqrt(sum((Tn .- Ta).^2.0.*dz))/L
        end

        return ϵ

    end

    tolerance = 1e-1

    # U1 validation
    u1, sim_u1, cfg_u1 = setup_btes_single(:u1);
    res_u1 = simulate_reservoir(u1; simulator=sim_u1, config=cfg_u1)
    analytical_u1 = Fimbul.analytical_closed_loop_u1(Q, T_in, T_rock,
    ρf, Cpf, u1.model.models[:BTES_supply].data_domain)
    ϵ_u1 = compute_error(u1, res_u1, analytical_u1)
    @test ϵ_u1 < tolerance

    # Coaxial validation with outer inlet
    coax_out, sim_coax_out, cfg_coax_out = setup_btes_single(:coaxial, inlet = :outer);
    res_coax_out = simulate_reservoir(coax_out; simulator=sim_coax_out, config=cfg_coax_out)
    analytical_coax_out = Fimbul.analytical_closed_loop_coaxial(Q, T_in, T_rock,
    ρf, Cpf, coax_out.model.models[:BTES_supply].data_domain;
    inlet = :outer)
    ϵ_coax_out = compute_error(coax_out, res_coax_out, analytical_coax_out)
    @test ϵ_coax_out < tolerance

    # Coaxial validation with inner inlet
    coax_in, sim_coax_in, cfg_coax_in = setup_btes_single(:coaxial, inlet = :inner);
    res_coax_in = simulate_reservoir(coax_in; simulator=sim_coax_in, config=cfg_coax_in)
    analytical_coax_in = Fimbul.analytical_closed_loop_coaxial(Q, T_in, T_rock,
    ρf, Cpf, coax_in.model.models[:BTES_supply].data_domain;
    inlet = :inner)
    ϵ_coax_in = compute_error(coax_in, res_coax_in, analytical_coax_in)
    @test ϵ_coax_in < tolerance

end