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

    @testset "Borehole fill" begin
        # Predefined material with overrides
        grout = Fimbul.borehole_fill_properties(:grout)
        @test grout == Fimbul.BOREHOLE_FILL_MATERIALS[:grout]
        props = Fimbul.borehole_fill_properties(:grout; density = 1800.0)
        @test props.density == 1800.0
        @test props.thermal_conductivity == grout.thermal_conductivity
        @test props.heat_capacity == grout.heat_capacity
        # Custom material
        custom = (thermal_conductivity = 1.0, density = 1000.0, heat_capacity = 800.0)
        @test Fimbul.borehole_fill_properties(custom) == custom
        @test_throws ErrorException Fimbul.borehole_fill_properties(:unknown)
        @test_throws ErrorException Fimbul.borehole_fill_properties((density = 1.0,))
        # Fill properties end up in the well domain
        for type in (:u1, :coaxial)
            well = setup(type; fill_material = custom, fill_heat_capacity = 900.0)
            gc = findall(isequal(:grout_left), well[:tag])
            gc = isempty(gc) ? findall(isequal(:grout), well[:tag]) : gc
            @test all(well[:grouting_thermal_conductivity][gc] .== 1.0)
            @test all(well[:material_density] .== 1000.0)
            @test all(well[:material_heat_capacity] .== 900.0)
        end
        well = setup(:generic; fill_thermal_conductivity = 1.0)
        @test all(well[:grouting_thermal_conductivity] .== 1.0)
    end

    @testset "Removed grout keyword arguments" begin
        # Removed keyword arguments give an error pointing to the new ones
        for type in (:u1, :coaxial)
            for (k, v) in ((:radius_grout, 70e-3),
                    (:grouting_thermal_conductivity, 2.0),
                    (:grouting_density, 2000.0),
                    (:grouting_heat_capacity, 1500.0))
                @test_throws ErrorException setup(type; NamedTuple{(k,)}((v,))...)
            end
        end
        @test_throws ErrorException setup(:generic; grouting_thermal_conductivity = 2.0)
    end

    @testset "Coaxial grouting thickness" begin
        r_po = sqrt(2)*50e-3
        # Default: grout nodes with borehole radius r_po + 30 mm
        well = setup(:coaxial)
        @test count(isequal(:grout), well[:tag]) == nz
        @test all(well[:radius_grout] .≈ r_po + 30e-3)
        # Zero thickness: no grout nodes, the outer pipe is perforated and
        # exchanges heat (but not mass) with the reservoir through its wall
        well = setup(:coaxial; grouting_thickness = 0.0)
        W = physical_representation(well)
        @test number_of_cells(W) == 2nz + 1
        @test count(isequal(:grout), well[:tag]) == 0
        @test all(well[:tag][W.perforations.self[1:end-1]] .== :pipe_outer)
        @test all(well[:well_index] .== 0.0)
        @test all(well[:volume_override_grouting] .== 0.0)
        L = well[:cell_length][1]
        R_co = log(r_po/(r_po - 3e-3))/(2π*0.38)
        @test all(well[:thermal_well_index][1:end-1] .≈ L/R_co)
        @test_throws ErrorException setup(:coaxial; grouting_thickness = -1e-3)
    end

    function simulate_closed_loop(type; n_step = 5, kwarg...)
        # The generic pipe flows downwards from the top to the bottom cell. The
        # outlet pressure must exceed the hydrostatic pressure of the pipe, so
        # that the pressure stays positive at the top.
        well = setup(type; kwarg...)
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

    @testset "Coaxial well without grout" begin
        _, res, inj, prod = simulate_closed_loop(:coaxial; grouting_thickness = 0.0)
        ws = res.wells
        @test ws[inj][:rate][end] ≈ -ws[prod][:rate][end] rtol = 1e-6
        T_out = ws[prod][:temperature][end]
        @test T_rock < T_out < T_in
        # Less thermal resistance than with grout: more heat is extracted
        _, res_grout, _, _ = simulate_closed_loop(:coaxial)
        @test T_out < res_grout.wells[prod][:temperature][end]
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
