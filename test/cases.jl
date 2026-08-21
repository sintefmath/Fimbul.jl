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

##
@testset "BTES topology" begin

    kw = (num_wells = 6, num_sectors = 2, well_spacing = 5.0,
        depths = [0.0, 0.5, 50.0, 65.0], n_z = [1, 2, 1], n_xy = 2, num_years = 1)

    # ## Sectors in parallel, boreholes within each sector in series
    # One pair of wells per borehole, each sector fed at its first well
    case = btes(:sunflower; kw...)
    @test case.input_data[:topology] == :sectors_parallel
    @test length(well_symbols(case.model)) == 2*kw.num_wells
    control = case.forces[1][:Facility].control
    sectors = case.input_data[:sectors]
    @test length(sectors) == kw.num_sectors
    for wells in values(sectors)
        supply = filter(w -> endswith(String(w), "_supply"), wells)
        @test control[supply[1]].target isa TotalRateTarget
        for k in 2:length(supply)
            target = control[supply[k]].target
            @test target isa JutulDarcy.ReinjectionTarget
            @test target.wells == [Symbol(replace(String(supply[k-1]), "_supply" => "_return"))]
            @test target.fraction == 1.0
        end
        for well in filter(w -> endswith(String(w), "_return"), wells)
            @test control[well] isa ProducerControl
        end
    end
    simulate_reservoir(case[1:1])

    # ## Sectors in series, boreholes within each sector in parallel
    # One pair of wells per sector, with the boreholes of the sector hung
    # between a shared inlet and a shared outlet node
    case = btes(:sunflower; topology = :sectors_series, kw...)
    @test case.input_data[:topology] == :sectors_series
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
        name = Symbol("S$sno", "_supply")
        @test control[Symbol("S$sno", "_return")] isa ProducerControl
        @test case.input_data[:sectors][Symbol("S$sno")] ==
            [name, Symbol("S$sno", "_return")]
        well_model = case.model.models[name]
        well = well_model.domain.representation
        N = well.neighborship
        # Two manifold nodes, plus a pipe and a grout node per perforation
        @test well.num_nodes == 2 + 2*well.num_perforations
        @test well.end_nodes == [2]
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
        z_well = well_model.data_domain[:cell_centroids, Cells()][3, well.perforations.self]
        z_reservoir = geo.cell_centroids[3, well.perforations.reservoir]
        @test z_well ≈ z_reservoir
        # Boreholes are identifiable after the fact, manifold nodes tagged 0
        boreholes = well_model.data_domain[:borehole, Cells()]
        @test sort(unique(boreholes)) == collect(0:count(N[1, :] .== 1))
    end
    simulate_reservoir(case[1:1])

    @test_throws ErrorException btes(:sunflower; topology = :nonsense, kw...)

end
