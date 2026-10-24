using Jutul, JutulDarcy, Fimbul, Test

@testset "Closed-loop wells" begin

    # Vertical column of cells for the closed-loop wells
    nz = 10
    msh = CartesianMesh((1, 1, nz), (1000.0, 1000.0, 100.0))
    reservoir = reservoir_domain(msh;
        rock_density = 2650.0,
        rock_heat_capacity = 900.0,
        rock_thermal_conductivity = 2.5,
        fluid_thermal_conductivity = 0.6,
        component_heat_capacity = 4184.0,
        porosity = 0.01,
        permeability = 1e-15
    )
    cells = collect(1:nz)
    T_rock = convert_to_si(10.0, :Celsius)
    T_in = convert_to_si(80.0, :Celsius)
    rate = 20.0si_unit(:meter)^3/si_unit(:day)

    setup(type; name = :CL, cells = cells, kwarg...) = setup_closed_loop_well(
        reservoir, cells; closed_loop_type = type, name = name, kwarg...)

    # Expected names, number of supply nodes and supply end node for each type
    types = [
        (:generic, ClosedLoopGenericWell(), (:CL_supply, :CL_return), nz, nz),
        (:u1, ClosedLoop1UWell(), (:CL_supply, :CL_return), 4nz, 2nz),
        (:coaxial, ClosedLoopCoaxialWell(), (:CL_inner, :CL_outer), 3nz, nz + 1),
    ]

    @testset "Types and names" begin
        for (symbol, type, names, _, _) in types
            @test Fimbul.closed_loop_well_type(symbol) == type
            @test Fimbul.closed_loop_well_type(type) == type
            @test closed_loop_well_names(symbol, :CL) == names
            @test closed_loop_well_names(type, :CL) == names
            @test type isa AbstractClosedLoopWell
        end
        @test_throws ErrorException Fimbul.closed_loop_well_type(:simple)
        @test_throws ErrorException Fimbul.closed_loop_well_type(:unknown)
        # Wells that are not closed-loop wells
        w = setup_vertical_well(reservoir, 1, 1; name = :Inj)
        @test !Fimbul.is_closed_loop_well(w)
        @test_throws ErrorException closed_loop_well_names(w)
    end

    @testset "Merged well domain: $symbol" for (symbol, type, names, n_supply, end_node) in types
        for t in (symbol, type)
            well = setup(t)
            W = physical_representation(well)
            info = W.multiwell
            # One merged well with the closed-loop type and two named parts
            @test Fimbul.is_closed_loop_well(well)
            @test W.type == type
            @test W.name == :CL
            @test closed_loop_well_names(well) == names
            @test info.names == collect(names)
            # The return node comes after the supply nodes
            @test number_of_cells(W) == n_supply + 1
            @test info.top_nodes == [1, n_supply + 1]
            @test W.end_nodes == [end_node, n_supply + 1]
            # The last segment connects the supply end node to the return node
            @test W.neighborship[:, end] == [end_node, n_supply + 1]
            @test W.segment_models[end] isa SegmentWellBoreFrictionHB
            @test Jutul.indirection_range(info.faces, 2) == W.num_segments:W.num_segments
            # Zero-length connection without heat conduction
            x = well[:cell_centroids]
            @test x[:, end] == x[:, end_node]
            @test well[:material_thermal_conductivity][end] == 0.0
            # Return node copies the end node, but holds fluid only
            @test well[:radius][end] == well[:radius][end_node]
            if haskey(well, :volume_override_hole)
                @test well[:volume_override_hole][end] == well[:volume_override_hole][end_node]
                @test well[:volume_override_casing][end] == 0.0
                @test well[:volume_override_grouting][end] == 0.0
            end
            # Return perforation exchanges neither mass nor heat
            @test W.perforations.self[end] == n_supply + 1
            @test well[:well_index][end] == 0.0
            @test well[:thermal_well_index][end] == 0.0
            if haskey(well, :tag)
                @test well[:tag][end] == :return
                @test count(isequal(:return), well[:tag]) == 1
            end
        end
    end

    @testset "Generic well topology" begin
        # Default: single pipe from the first to the last cell
        well = setup(:generic)
        W = physical_representation(well)
        @test W.neighborship[:, 1:end-1] == vcat((1:nz-1)', (2:nz)')
        @test W.perforations.reservoir == vcat(cells, cells[end])
        @test well[:cell_centroids][:, 1:nz] == reservoir[:cell_centroids][:, cells]
        @test well[:tag] == vcat(fill(:pipe, nz), :return)
        # Return reservoir cell can be given explicitly
        well = setup(:generic; return_reservoir_cell = cells[3])
        @test physical_representation(well).perforations.reservoir[end] == cells[3]
        # Branching topology: two parallel branches from node 2 that join in node 6
        N = [1 2 3 2 4 5 6;
             2 3 6 4 5 6 7]
        rcells = [1, 2, 3, 2, 3, 4, 5]
        well = setup(:generic;
            cells = rcells,
            neighborship = N,
            perforation_cells_well = collect(1:7),
            well_cell_centers = reservoir[:cell_centroids][:, rcells],
            end_nodes = [7],
            name = :AGS
        )
        W = physical_representation(well)
        @test W.neighborship[:, 1:end-1] == N
        @test W.neighborship[:, end] == [7, 8]
        @test W.perforations.reservoir[end] == rcells[end]
        @test closed_loop_well_names(well) == (:AGS_supply, :AGS_return)
        # Topology arguments must be given together
        @test_throws ErrorException setup(:generic; end_nodes = [nz])
        @test_throws ErrorException setup(:generic;
            well_cell_centers = reservoir[:cell_centroids][:, cells])
        @test_throws ErrorException setup(:generic; cells = rcells, neighborship = N)
        # The supply part must have a single end node
        @test_throws ErrorException setup(:generic;
            cells = rcells,
            neighborship = N,
            perforation_cells_well = collect(1:7),
            well_cell_centers = reservoir[:cell_centroids][:, rcells],
            end_nodes = [6, 7]
        )
    end

    @testset "1U and coaxial return cell" begin
        # The end node is not perforated: return defaults to the first cell
        for type in (:u1, :coaxial)
            well = setup(type)
            @test physical_representation(well).perforations.reservoir[end] == cells[1]
        end
    end

    @testset "Topology arguments: $symbol" for (symbol, args) in [
            (:u1, (:pipe_cells, :grout_cells, :well_cell_centers, :end_nodes)),
            (:coaxial, (:pipe_cells_inner, :pipe_cells_outer, :grout_cells,
                :well_cell_centers, :end_nodes)),
        ]
        # Topology arguments are only allowed together with neighborship
        for arg in args
            @test_throws ErrorException setup(symbol; NamedTuple{(arg,)}(([1],))...)
        end
    end

    @testset "Well direction" begin
        # A direction per reservoir cell is accepted by all types
        dir = fill([0.0, 0.0, 10.0], nz)
        for type in (:generic, :u1, :coaxial)
            well = setup(type; dir = dir)
            @test Fimbul.is_closed_loop_well(well)
            @test all(isapprox.(well[:cell_length][1:end-1][well[:cell_length][1:end-1] .> 0], 10.0))
        end
    end

    @testset "BTES wells" begin
        # BTES wells default to 1U wells
        well = setup_vertical_btes_well(reservoir, 1, 1)
        @test physical_representation(well).type == ClosedLoop1UWell()
        @test physical_representation(well).name == :BTES
        @test closed_loop_well_names(well) == (:BTES_supply, :BTES_return)
        well = setup_btes_well(reservoir, cells; closed_loop_type = :coaxial)
        @test physical_representation(well).type == ClosedLoopCoaxialWell()
        # Vertical closed-loop wells default to generic wells
        well = setup_vertical_closed_loop_well(reservoir, 1, 1; heel = 2, toe = 5)
        @test physical_representation(well).type == ClosedLoopGenericWell()
        @test physical_representation(well).perforations.reservoir == [2, 3, 4, 5, 5]
    end

    @testset "Pressure check" begin
        sys = SinglePhaseSystem(AqueousPhase(), reference_density = 1000.0)
        function setup_model(type)
            well = setup(type)
            model = setup_reservoir_model(reservoir, sys;
                wells = [well], thermal = true, can_shut_wells = false)
            return model, well
        end
        ctrl_inj = InjectorControl(TotalRateTarget(rate), [1.0],
            density = 1000.0, temperature = T_in)
        ctrl_bhp(p) = ProducerControl(BottomHolePressureTarget(p))
        # Generic pipe flowing downwards: the BHP applies at the bottom node
        model, well = setup_model(:generic)
        inj, prod = closed_loop_well_names(well)
        z = well[:cell_centroids][3, :]
        p_min = JutulDarcy.DEFAULT_MINIMUM_PRESSURE
        p_required = p_min + 1000.0*Jutul.gravity_constant*(maximum(z) - minimum(z))
        low = Dict(inj => ctrl_inj, prod => ctrl_bhp(0.99*p_required))
        high = Dict(inj => ctrl_inj, prod => ctrl_bhp(1.01*p_required))
        @test (@test_logs (:warn, r"too low") check_closed_loop_pressure(model, low)) == false
        @test (@test_logs check_closed_loop_pressure(model, high)) == true
        # Forces and vectors of forces
        forces_low = setup_reservoir_forces(model, control = low)
        forces_high = setup_reservoir_forces(model, control = high)
        @test (@test_logs (:warn, r"too low") check_closed_loop_pressure(model, forces_low)) == false
        @test (@test_logs check_closed_loop_pressure(model, [forces_high, forces_high])) == true
        @test (@test_logs (:warn, r"too low") check_closed_loop_pressure(model,
            [forces_high, forces_low, forces_low])) == false
        # Custom minimum pressure
        @test (@test_logs check_closed_loop_pressure(model, low; p_min = 0.0)) == true
        # Rate-controlled parts are not checked
        rate_only = Dict(inj => ctrl_inj,
            prod => ProducerControl(TotalRateTarget(-rate)))
        @test (@test_logs check_closed_loop_pressure(model, rate_only)) == true
        # 1U and coaxial wells have both ends at the top
        for type in (:u1, :coaxial)
            model, well = setup_model(type)
            sup, ret = closed_loop_well_names(well)
            for ctrl in (Dict(sup => ctrl_inj, ret => ctrl_bhp(1.0si_unit(:atm))),
                    Dict(sup => ctrl_bhp(1.0si_unit(:atm)), ret => ctrl_inj))
                @test (@test_logs check_closed_loop_pressure(model, ctrl)) == true
            end
        end
    end

    function simulate_closed_loop(type; n_step = 5)
        # The generic pipe flows downwards from the top to the bottom cell. The
        # outlet pressure must exceed the hydrostatic pressure of the pipe, so
        # that the pressure stays positive at the top.
        well = setup(type)
        p_loop = 15.0si_unit(:atm)
        inj, prod = closed_loop_well_names(well)
        sys = SinglePhaseSystem(AqueousPhase(), reference_density = 1000.0)
        model = setup_reservoir_model(reservoir, sys;
            wells = [well], thermal = true, can_shut_wells = false)
        control = Dict(
            inj => InjectorControl(TotalRateTarget(rate), [1.0],
                density = 1000.0, temperature = T_in),
            prod => ProducerControl(BottomHolePressureTarget(p_loop))
        )
        forces = setup_reservoir_forces(model, control = control)
        state0 = setup_reservoir_state(model,
            Pressure = p_loop, Temperature = T_rock)
        dt = fill(si_unit(:day), n_step)
        res = simulate_reservoir(state0, model, dt;
            forces = forces, info_level = -1)
        return model, res, inj, prod
    end

    @testset "Model and simulation: $symbol" for (symbol, _, names, n_supply, _) in types
        model, res, inj, prod = simulate_closed_loop(symbol)
        # One well model, with both parts controlled by the facility
        @test haskey(model.models, :CL)
        @test !haskey(model.models, inj)
        @test model.models[:Facility].domain.well_symbols == collect(names)
        @test sort(well_symbols(model)) == sort(collect(names))
        # Only standard well cross terms: the supply and return parts are
        # connected inside the merged well
        well_cts = filter(ct -> :CL in (ct.target, ct.source), model.cross_terms)
        @test !isempty(well_cts)
        @test all(ct -> JutulDarcy.WellMerging.is_standard_well_cross_term(ct.cross_term), well_cts)
        # Mass in equals mass out, and the fluid is cooled by the rock
        ws = res.wells
        q_in = ws[inj][:mass_rate][end]
        q_out = ws[prod][:mass_rate][end]
        @test q_in > 0
        @test q_in ≈ -q_out rtol = 1e-6
        T_out = ws[prod][:temperature][end]
        @test T_rock < T_out < T_in
        # States hold supply and return nodes of the merged well
        state = res.result.states[end][:CL]
        @test length(state[:Temperature]) == n_supply + 1
        @test all(T_rock - 0.1 .<= state[:Temperature] .<= T_in + 0.1)
    end

    @testset "Generic pipe temperature" begin
        _, res, _, _ = simulate_closed_loop(:generic)
        T = res.result.states[end][:CL][:Temperature]
        # Fluid cools monotonically from inlet to outlet, and the return node
        # has the temperature of the outlet
        @test issorted(T[1:nz], rev = true)
        @test T[end] ≈ T[nz] rtol = 1e-3
    end

end
