using Fimbul, JutulDarcy, Test

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