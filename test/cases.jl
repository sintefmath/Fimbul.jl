using Fimbul, Jutul, JutulDarcy, Test

@testset "Egg" begin
    case = egg_geothermal()
    res = simulate_reservoir(case[1:1])
    @test true
    case = egg_geothermal_doublet()
    res = simulate_reservoir(case[1:1])
    @test true
    case = egg_ates()
    res = simulate_reservoir(case[1:1])
    @test true
end
@testset "Schedule from durations" begin

    forces = [:a, :b, :c]
    durations = [10.0, 20.0, 30.0]

    function do_test(dt, f, ts, num_cycles)
        @test length(f) == length(dt)
        @test length(ts) == length(dt) + 1
        @test sum(dt) ≈ sum(durations)*num_cycles
        start = 1
        for (n, t) in enumerate(cumsum(repeat(durations, num_cycles)))
            stop = findall(isapprox.(ts, t))
            @test length(stop) == 1
            f[start:stop[1]-1] .= forces[mod(n-1, length(forces)) + 1]
            start = stop[1]
        end
    end

    for num_cycles in (1, 5)

        for num_reports in (missing, 2, [2, 4, 6])
            dt, f, ts = make_schedule(forces, durations;
            num_reports = num_reports, num_cycles = num_cycles)
            do_test(dt, f, ts, num_cycles)
        end

        for report_interval in (5.0, [5.0, 10.0, 15.0])
            dt, forces_out, timestamps = make_schedule(forces, durations;
            report_interval = report_interval, num_cycles = num_cycles)
            do_test(dt, forces_out, timestamps, num_cycles)
        end
    end

end

using Dates
@testset "Schedule from time periods" begin

    start_year = Dates.year(now())
    forces = [:a, :b, :c]
    function do_test(periods, dt, f, ts, num_years)
        @test length(f) == length(dt)
        @test length(ts) == length(dt) + 1
        start = 1
        for year in (0:(num_years-1)) .+ start_year
            periods_year = Fimbul.process_periods(year, periods)
            for (n, t) in enumerate(periods_year)
                stop = findall(ts .== t)
                @test length(stop) == 1
                f[start:stop[1]-1] .= forces[mod(n-1, length(forces)) + 1]
            end
        end
        @test Dates.value(ts[end]-ts[1])*1e-3 ≈ sum(dt)
    end

    for num_years in (1, 5)

        for periods in (
            [3,5,10,3],
            ["March", "May", "October", "March"],
            [(3,1), (5,1), (10,1), (3,1)],
            [(3,1,12), (5,1,12), (10,1,12), (3,1,12)]
        )
            for num_reports in (missing, 2, [2, 4, 6])
                dt, f, ts = make_schedule(forces, periods;
                num_reports = num_reports, num_years = num_years, start_year)
                do_test(periods, dt, f, ts, num_years)
            end

            for report_interval in (7si_unit(:day), [3, 7, 14].*si_unit(:day))
                dt, f, ts = make_schedule(forces, periods;
                report_interval = report_interval, num_years = num_years, start_year = start_year)
                do_test(periods, dt, f, ts, num_years)
            end
        end
    end

end

@testset "H2O pressure-enthalpy formulation" begin
    case = benchmark_ht_1d(benchmark_case = :a, nx = 10)
    rmodel = reservoir_model(case.model)
    @test haskey(rmodel.primary_variables, :Enthalpy)
    @test !haskey(rmodel.primary_variables, :Temperature)
end

@testset "Magmatic intrusion" begin
    # Tensor grid: total cell count, uniform refined cells and cap-rock layers
    dz = Fimbul.magmatic_intrusion_cell_sizes(5000.0, 40, 1750.0, 3250.0;
        top = (500.0, 3))
    @test length(dz) == 40
    @test sum(dz) ≈ 5000.0
    @test dz[1:3] ≈ fill(500.0/3, 3)
    @test maximum(max.(dz[2:end]./dz[1:end-1], dz[1:end-1]./dz[2:end])) < 1.3
    # Automatic number of cap-rock layers, refined all the way to the cap
    dz = Fimbul.magmatic_intrusion_cell_sizes(5000.0, 40, 500.0, 3000.0;
        top = (500.0, nothing))
    @test length(dz) == 40
    @test sum(dz) ≈ 5000.0
    @test any(isapprox.(cumsum(dz), 500.0))
    @test maximum(max.(dz[2:end]./dz[1:end-1], dz[1:end-1]./dz[2:end])) < 1.3
    case = magmatic_intrusion(dims = (12, 12, 16), num_years = 0)
    @test number_of_cells(physical_representation(
        reservoir_model(case.model).data_domain)) == 12*12*16
    # Simulations on a small uniform grid (with a moderate intrusion
    # temperature, since the grid is very coarse)
    args = (dims = (6, 6, 10), injector_position = (-100.0, 250.0),
        cell_growth = 1.0, temperature_intrusion = convert_to_si(500.0, :Celsius))
    # Natural state without wells
    case_ns = magmatic_intrusion(; args..., num_years = 0,
        natural_state_time = 300*si_unit(:year))
    rmodel = reservoir_model(case_ns.model)
    @test haskey(rmodel.primary_variables, :Enthalpy)
    @test isempty(get_model_wells(case_ns.model))
    res_ns = simulate_reservoir(case_ns, info_level = -1)
    @test length(res_ns.states) == length(case_ns.dt)
    # Production with wells, restarted from the natural state
    case = magmatic_intrusion(; args..., natural_state_time = 0.0,
        num_years = 1, report_interval = si_unit(:year)/4,
        initial_state = res_ns.states[end])
    @test Set(keys(get_model_wells(case.model))) == Set([:Producer, :Injector])
    @test case.state0[:Reservoir][:Pressure] ≈ res_ns.states[end][:Pressure]
    res = simulate_reservoir(case[1:2], info_level = -1)
    @test length(res.states) == 2
    # The injector is not included if it is never active
    case_noinj = magmatic_intrusion(; args..., natural_state_time = 0.0,
        num_years = 1, injection_years = 0)
    @test collect(keys(get_model_wells(case_noinj.model))) == [:Producer]
end

##
@testset "BTES topology" begin

    kw = (num_wells = 6, num_sectors = 2,
        depths = [0.0, 0.5, 50.0, 65.0], n_z = [1, 2, 1], n_xy = 2, num_years = 1)

    is_closed_loop(case) = [k for (k, m) in pairs(case.model.models)
        if Fimbul.is_closed_loop_well(m.data_domain)]
    part_names(case, well) = closed_loop_well_names(case.model.models[well].data_domain)

    # ## Sectors in parallel, boreholes within each sector in series
    # Default: one closed-loop well per sector, with the boreholes in series
    case = btes(:sunflower; kw...)
    @test case.input_data[:topology] == :sectors_parallel
    @test length(is_closed_loop(case)) == kw.num_sectors
    control = case.forces[1][:Facility].control
    geo = tpfv_geometry(physical_representation(reservoir_model(case.model).data_domain))
    for sno in 1:kw.num_sectors
        name = Symbol("S$sno")
        @test case.input_data[:sectors][name] == [name]
        sup, ret = part_names(case, name)
        @test control[sup].target isa TotalRateTarget
        @test control[ret] isa ProducerControl
        well_model = case.model.models[name]
        well = well_model.domain.representation
        N = well.neighborship
        boreholes = well_model.data_domain[:borehole, Cells()]
        nb = maximum(boreholes)
        @test sort(unique(boreholes[1:end-1])) == collect(1:nb)
        # Fluid enters the first borehole and leaves the last one
        @test boreholes[1] == 1
        @test boreholes[well.end_nodes[1]] == nb
        # Consecutive boreholes are connected by a resistance-free surface
        # pipe that carries flow. The last segment connects the return part.
        connectors = [j for j in 1:size(N, 2) - 1
            if boreholes[N[1, j]] != boreholes[N[2, j]]]
        @test length(connectors) == nb - 1
        @test all(boreholes[N[2, j]] == boreholes[N[1, j]] + 1 for j in connectors)
        @test all(!JutulDarcy.well_segment_is_closed(well.segment_models[j]) for j in connectors)
        @test all(case.parameters[name][:SegmentLength][connectors] .== 0.0)
        # Perforations pair well nodes with reservoir cells at the same depth
        perfs = 1:well.num_perforations - 1
        z_well = well_model.data_domain[:cell_centroids, Cells()][3, well.perforations.self[perfs]]
        z_reservoir = geo.cell_centroids[3, well.perforations.reservoir[perfs]]
        @test z_well ≈ z_reservoir
    end
    # The mesh is refined around the wells: the perforated cells have the
    # horizontal size set by well_spacing/n_xy
    hxy = 5.0/kw.n_xy
    rcells = unique(vcat([physical_representation(case.model.models[w]).perforations.reservoir
        for w in keys(case.input_data[:sectors])]...))
    msh = physical_representation(reservoir_model(case.model).data_domain)
    area = [geo.volumes[c]/cell_dims(msh, c)[3] for c in rcells]
    @test all(sqrt.(area) .< 1.5*hxy)
    simulate_reservoir(case[1:1])
    # Reversed discharge: the sector well is fed at its return part, and
    # produced from its supply part
    case = btes(:sunflower; reversed_discharge = true, kw...)
    T_dch = convert_to_si(10.0, :Celsius)
    is_dch(c) = c isa InjectorControl && c.temperature == T_dch
    dch = findfirst(i -> any(is_dch, values(case.forces[i][:Facility].control)),
        eachindex(case.forces))
    control = case.forces[dch][:Facility].control
    for sno in 1:kw.num_sectors
        sup, ret = part_names(case, Symbol("S$sno"))
        @test is_dch(control[ret])
        @test control[sup] isa ProducerControl
        @test case.forces[1][:Facility].control[sup].target isa TotalRateTarget
    end
    res = simulate_reservoir(case[1:dch]; error_on_incomplete = true, info_level = -1)
    q = res.result.states[end][:Facility][:TotalSurfaceMassRate]
    ws = case.model.models[:Facility].domain.well_symbols
    sup, ret = part_names(case, :S1)
    @test q[findfirst(==(ret), ws)] > 0
    @test q[findfirst(==(sup), ws)] < 0
    @test_throws ErrorException btes(:sunflower; well_grouping = :nonsense, kw...)
    @test_throws ErrorException btes(:sunflower; topology = :sectors_series,
        well_grouping = :borehole, kw...)

    # One closed-loop well per borehole, each sector fed at its first well
    case = btes(:sunflower; well_grouping = :borehole, kw...)
    @test case.input_data[:topology] == :sectors_parallel
    @test length(is_closed_loop(case)) == kw.num_wells
    @test length(well_symbols(case.model)) == 2*kw.num_wells
    control = case.forces[1][:Facility].control
    sectors = case.input_data[:sectors]
    @test length(sectors) == kw.num_sectors
    # Sectors list the closed-loop wells
    @test sort(vcat(values(sectors)...)) == sort(is_closed_loop(case))
    for wells in values(sectors)
        supply = [first(part_names(case, w)) for w in wells]
        ret = [last(part_names(case, w)) for w in wells]
        @test control[supply[1]].target isa TotalRateTarget
        for k in 2:length(supply)
            target = control[supply[k]].target
            @test target isa JutulDarcy.ReinjectionTarget
            @test target.wells == [ret[k-1]]
            @test target.fraction == 1.0
        end
        for well in ret
            @test control[well] isa ProducerControl
        end
    end
    simulate_reservoir(case[1:1])

    # ## Sectors in series, boreholes within each sector in parallel
    # One closed-loop well per sector, with the boreholes of the sector hung
    # between a shared inlet and a shared outlet node
    case = btes(:sunflower; topology = :sectors_series, kw...)
    @test case.input_data[:topology] == :sectors_series
    @test length(is_closed_loop(case)) == kw.num_sectors
    @test length(well_symbols(case.model)) == 2*kw.num_sectors
    control = case.forces[1][:Facility].control
    @test control[:S1_supply].target isa TotalRateTarget
    for sno in 2:kw.num_sectors
        target = control[Symbol("S$sno", "_supply")].target
        @test target isa JutulDarcy.ReinjectionTarget
        @test target.wells == [Symbol("S$(sno-1)", "_return")]
    end
    geo = tpfv_geometry(physical_representation(reservoir_model(case.model).data_domain))
    for sno in 1:kw.num_sectors
        name = Symbol("S$sno")
        @test part_names(case, name) == (Symbol("S$sno", "_supply"), Symbol("S$sno", "_return"))
        @test control[Symbol("S$sno", "_return")] isa ProducerControl
        @test case.input_data[:sectors][name] == [name]
        well_model = case.model.models[name]
        well = well_model.domain.representation
        N = well.neighborship
        # Supply part: two manifold nodes, plus a pipe and a grout node per
        # perforation. The return part is a single node (and perforation)
        # connected to the shared outlet.
        n_supply = well.num_nodes - 1
        @test n_supply == 2 + 2*(well.num_perforations - 1)
        @test well.end_nodes == [2, n_supply + 1]
        @test N[:, end] == [2, n_supply + 1]
        # Every borehole is fed from the shared inlet and returns to the shared
        # outlet, and the segments that do so carry flow
        manifold = findall(j -> N[1, j] == 1 || N[2, j] == 2, axes(N, 2))
        @test count(N[1, :] .== 1) == count(N[2, :] .== 2)
        @test length(manifold) == 2*count(N[1, :] .== 1)
        @test all(!JutulDarcy.well_segment_is_closed(well.segment_models[j])
            for j in manifold)
        # The manifold is resistance-free, so the split between the boreholes
        # cannot depend on where the shared nodes were placed
        @test all(case.parameters[name][:SegmentLength][manifold] .== 0.0)
        # Perforations pair well nodes with reservoir cells at the same depth
        perfs = 1:well.num_perforations - 1
        z_well = well_model.data_domain[:cell_centroids, Cells()][3, well.perforations.self[perfs]]
        z_reservoir = geo.cell_centroids[3, well.perforations.reservoir[perfs]]
        @test z_well ≈ z_reservoir
        # Boreholes are identifiable after the fact, manifold nodes tagged 0
        boreholes = well_model.data_domain[:borehole, Cells()]
        @test sort(unique(boreholes)) == collect(0:count(N[1, :] .== 1))
    end
    simulate_reservoir(case[1:1])

    # ## Rest periods
    # The loops stop circulating, but are not sealed: the supply side holds the
    # wellhead pressure so that the loop fluid can expand and contract with
    # temperature. A sealed, cooling loop does not converge.
    forces = case.forces
    rest = findfirst(i -> forces[i][:Facility].control[:S1_supply] isa ClosedLoopRestControl, eachindex(forces))
    @test !isnothing(rest)
    @test rest > 1
    control = forces[rest][:Facility].control
    for sno in 1:kw.num_sectors
        ctrl = control[Symbol("S$sno", "_supply")]
        @test ctrl isa ClosedLoopRestControl
        @test ctrl.target isa BottomHolePressureTarget
        @test ctrl.target.value ≈ 1.0si_unit(:atm)
        @test control[Symbol("S$sno", "_return")] isa DisabledControl
    end
    # Charge, then take the first rest step from the charged state
    simulate_reservoir(case[1:rest]; error_on_incomplete = true, info_level = -1)

    @test_throws ErrorException btes(:sunflower; topology = :nonsense, kw...)

end

##
@testset "BTES well patterns" begin
    field_wells(field) = vcat(field...)
    xy(field) = hcat([w[1:2, 1] for w in field_wells(field)]...)
    r_max(field) = maximum(sqrt.(sum(xy(field).^2, dims = 1)))
    for pattern in (:sunflower, :circular, :rectangular)
        size_kw = pattern == :rectangular ? (extent = (20.0, 10.0),) : (radius = 20.0,)
        # The number of wells is honored exactly, with and without a field size
        for N in (1, 7, 17, 37)
            for kw in ((num_wells = N,), (num_wells = N, size_kw...))
                field, s = Fimbul.btes_well_field(pattern; num_sectors = 3, kw...)
                @test length(field_wells(field)) == N
                @test s > 0
            end
        end
        # A field size and a spacing give the number of wells, which fit
        # within the field
        field, s = Fimbul.btes_well_field(pattern; well_spacing = 4.0, size_kw...)
        @test s == 4.0
        if pattern == :rectangular
            span = maximum(xy(field), dims = 2) - minimum(xy(field), dims = 2)
            @test all(span .<= [20.0, 10.0] .+ 1e-8)
        else
            @test r_max(field) <= 20.0 + 1e-8
        end
        # Defaults: 48 wells at 5 m spacing
        field, s = Fimbul.btes_well_field(pattern)
        @test length(field_wells(field)) == 48
        @test s == 5.0
        # The number of wells and the spacing cannot both be given
        @test_throws ErrorException Fimbul.btes_well_field(pattern;
            num_wells = 10, well_spacing = 5.0)
    end
    # Circular sectors are consecutive groups of wells from the center
    # outwards, counter-clockwise on each ring
    field, s = Fimbul.btes_well_field(:circular; num_wells = 48, num_sectors = 6,
        sector_grouping = :sequential)
    @test all(length.(field) .== 8)
    r(w) = sqrt(w[1, 1]^2 + w[2, 1]^2)
    θ(w) = mod(atan(w[2, 1], w[1, 1]), 2π)
    wells = vcat(field...)
    rings = round.(Int, r.(wells)./s)
    @test issorted(rings)
    for k in unique(rings)
        @test issorted(θ.(wells[rings .== k]))
    end
    @test maximum(r.(field[1])) <= minimum(r.(field[2])) + 1e-8
    # By default, circular sectors are circle sectors, each spanning the rings
    field, s = Fimbul.btes_well_field(:circular; num_wells = 48, num_sectors = 6)
    @test all(length.(field) .== 8)
    @test all(maximum(r.(sector)) > 2*s for sector in field)
    @test_throws ErrorException Fimbul.btes_well_field(:rectangular;
        sector_grouping = :angular)
    # Sunflower sectors are bundles of spiral spines, ordered from the center
    # outwards, with the wells of each sector in sunflower order level by level
    field, s = Fimbul.btes_well_field(:sunflower; num_wells = 30, num_sectors = 5)
    @test all(length.(field) .== 6)
    sectors = Fimbul.group_into_sectors_spiral(30, 5)
    @test sort(vcat(sectors...)) == 1:30
    # With one spine per sector, the wells of a sector are F apart
    @test all(length(unique(diff(sector))) == 1 for sector in sectors)
    @test all(issorted(r.(sector)) for sector in field)
    for (N, ns) in ((12, 3), (36, 6), (48, 6), (100, 8))
        sectors = Fimbul.group_into_sectors_spiral(N, ns)
        @test sort(vcat(sectors...)) == 1:N
        sizes = length.(sectors)
        @test maximum(sizes) - minimum(sizes) <= 1
    end
    # Rectangular sectors as blocks of a given size
    field, s = Fimbul.btes_well_field(:rectangular; num_wells = 36,
        num_sectors = 6, sector_dims = (3, 2))
    @test all(length.(field) .== 6)
    for sector in field
        xy = hcat([w[1:2, 1] for w in sector]...)
        @test length(unique(xy[1, :])) == 3
        @test length(unique(xy[2, :])) == 2
    end
    # Without a block size, the blocks are as square as possible, and the
    # wells left at the top form the last sector
    field, s = Fimbul.btes_well_field(:rectangular; num_wells = 45, num_sectors = 5)
    @test all(length.(field) .== 9)
    nx_ny(sector) = (length(unique([w[1, 1] for w in sector])),
        length(unique([w[2, 1] for w in sector])))
    @test nx_ny.(field) == [(3, 3), (3, 3), (3, 3), (3, 3), (6, 2)]
    for N in 1:40, ns in 1:min(N, 6)
        field, s = Fimbul.btes_well_field(:rectangular; num_wells = N, num_sectors = ns)
        @test length(field) == ns
        @test sum(length, field) == N
    end
    @test_throws ErrorException Fimbul.btes_well_field(:sunflower;
        num_wells = 36, num_sectors = 6, sector_dims = (3, 2))
    @test_throws ErrorException Fimbul.btes_well_field(:sunflower;
        sector_grouping = :angular)
    # The sunflower pattern has no gap at the center: the largest empty circle
    # around the center is smaller than the spacing
    x = fibonacci_pattern_2d(30; spacing = 5.0)
    @test minimum(sqrt.(sum(x.^2, dims = 1))) < 5.0/2
    @test Fimbul.median(Fimbul.nearest_neighbor_distances(x)) ≈ 5.0
    # The field size matches the pattern
    @test_throws ErrorException Fimbul.btes_well_field(:rectangular; radius = 10.0)
    @test_throws ErrorException Fimbul.btes_well_field(:sunflower; extent = (10.0, 10.0))
    @test_throws ErrorException Fimbul.btes_well_field(:polygonal)
end
