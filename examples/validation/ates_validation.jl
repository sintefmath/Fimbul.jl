# # Analytical ATES validation against reference outputs
# <tags: Validation, ATES>
# This example sets up and runs the analytical ATES benchmark using
# `Fimbul.analytical_ates()`, then compares the numerical outputs against
# reference CSV files. The reference path currently points to the
# `jutul-agent-test` repository; you can later replace this with packaged
# artifacts or analytical expressions.

using Jutul, JutulDarcy, Fimbul
using Statistics
using DataFrames
using CSV
using CairoMakie, GLMakie

function state_temperature_c(state)
    return convert_from_si.(state[:Temperature], :Celsius)
end

function radial_index_cell_map(g)
    ncell = number_of_cells(g)
    mapping = Dict{Int, Vector{Int}}()
    # for rix in 0:100
    #     r = 
    #     cells = findall(isapprox.(radii, r; atol = 0.1))
    #     mapping[rix + 1] = cells
    # end
    # return radii
    for c in 1:ncell
        ijk = cell_ijk(g, c)
        j = ijk[2]
        push!(get!(mapping, j, Int[]), c)
    end
    return mapping
end

function layer_index_cell_map(g)
    ncell = number_of_cells(g)
    mapping = Dict{Int, Vector{Int}}()
    for c in 1:ncell
        ijk = cell_ijk(g, c)
        k = ijk[3]
        push!(get!(mapping, k, Int[]), c)
    end
    return mapping
end

# function layer_index_cell_map(g)
#     ncell = number_of_cells(g)
#     mapping = Dict{Int, Vector{Int}}()
#     geo = tpfv_geometry(g)
#     zmid = [2.5, 7.5, 12.5, 17.5]
#     for k in 1:4
#         cells = findall(isapprox.(geo.cell_centroids[3, :], zmid[k]; atol = 0.1))
#         mapping[k] = cells
#     end
#     return mapping
# end

function extract_layer_averaged_radial_profile(states, g; max_radius = 100)
    rmap = radial_index_cell_map(g)
    out = DataFrame(time_d = Float64[])
    for r in 0:max_radius
        out[!, Symbol("T_r$(lpad(string(r), 3, '0'))")] = Float64[]
    end
    state_offset = length(states) == 91 ? 1 : 0
    for d in 1:90
        s = states[d + state_offset]
        T = state_temperature_c(s)
        row = Float64[d]
        for r in 0:max_radius
            cells = rmap[r + 1]
            push!(row, mean(T[cells]))
        end
        push!(out, Tuple(row))
    end
    return out
end

function extract_cross_section(states, g, target_day; max_radius = 100)
    rmap = radial_index_cell_map(g)
    lmap = layer_index_cell_map(g)
    zmid = [-2.5, -7.5, -12.5, -17.5]
    out = DataFrame(elevation_m = zmid)
    state_offset = length(states) == 91 ? 1 : 0
    s = states[target_day + state_offset]
    T = state_temperature_c(s)
    for r in 0:max_radius
        col = Symbol("T_r$(lpad(string(r), 3, '0'))")
        vals = Float64[]
        for k in 1:4
            cells = intersect(rmap[r + 1], lmap[k])
            push!(vals, mean(T[cells]))
        end
        out[!, col] = vals
    end
    return out
end

function default_reference_dir()
    return get(
        ENV,
        "FIMBUL_ATES_REFERENCE_DIR",
        joinpath(homedir(), "Documents", "repositories", "papers", "jutul-agent-test", "Comparison_Analytical_Solution", "Results_analytical"),
    )
end

function compare_to_reference(daily, outdir, reference_dir)
    ref_daily_path = joinpath(reference_dir, "mf6_results_analytical_T_radial.csv")
    isfile(ref_daily_path) || error("Missing reference file: $(ref_daily_path)")

    ref_daily = CSV.read(ref_daily_path, DataFrame)
    common_cols = [c for c in names(ref_daily) if c != "time_d" && c in names(daily)]
    A = Matrix(select(ref_daily, common_cols))
    B = Matrix(select(daily, common_cols))
    daily_max = maximum(abs.(A .- B))
    daily_mean = mean(abs.(A .- B))

    rows = DataFrame(
        comparison_label = String[],
        report_day = String[],
        reference_file = String[],
        numerical_file = String[],
        compared_columns = String[],
        max_abs_diff_C = Float64[],
        mean_abs_diff_C = Float64[],
    )

    for d in [10, 30, 50, 70, 90]
        ref_path = joinpath(reference_dir, "mf6_cross_sections", "mf6_cross_section_t0$(lpad(string(d), 2, '0'))d.csv")
        out_path = joinpath(outdir, "cross_section_$(d)d.csv")
        isfile(ref_path) || error("Missing reference file: $(ref_path)")
        isfile(out_path) || error("Missing output file: $(out_path)")

        ref = CSV.read(ref_path, DataFrame)
        out = CSV.read(out_path, DataFrame)
        common_cols = [c for c in names(ref) if c != "elevation_m" && c in names(out)]
        A = Matrix(select(ref, common_cols))
        B = Matrix(select(out, common_cols))
        push!(rows, (
            "Cross-section at day $d",
            string(d),
            basename(ref_path),
            basename(out_path),
            "T_r000:T_r100",
            maximum(abs.(A .- B)),
            mean(abs.(A .- B)),
        ))
    end

    push!(rows, (
        "Overall radial profile over full simulation",
        "overall",
        basename(ref_daily_path),
        basename(joinpath(outdir, "daily_layer_averaged_radial_temperatures.csv")),
        "temperature_vs_radius_over_full_simulation",
        daily_max,
        daily_mean,
    ))

    CSV.write(joinpath(outdir, "comparison_summary.csv"), rows)

    numeric_rows = DataFrame(
        comparison_label = String[],
        report_day = String[],
        max_abs_diff_C = Float64[],
        mean_abs_diff_C = Float64[],
    )
    for r in eachrow(rows)
        push!(numeric_rows, (r.comparison_label, r.report_day, r.max_abs_diff_C, r.mean_abs_diff_C))
    end
    CSV.write(joinpath(outdir, "comparison_metrics.csv"), numeric_rows)

    return daily_max, daily_mean, rows
end

function radial_column_name(r)
    return "T_r$(lpad(string(r), 3, '0'))"
end

function write_timeseries_plot(outdir, reference_dir; distances = 0:10:100, report_days = [10, 30, 50, 70, 90])
    ref_dir = joinpath(reference_dir, "mf6_cross_sections")
    sample_cross = CSV.read(joinpath(outdir, "cross_section_$(first(report_days))d.csv"), DataFrame)
    sample_ref = CSV.read(joinpath(ref_dir, "mf6_cross_section_t0$(lpad(string(first(report_days)), 2, '0'))d.csv"), DataFrame)

    columns = [
        radial_column_name(r)
        for r in distances
        if radial_column_name(r) in names(sample_cross) && radial_column_name(r) in names(sample_ref)
    ]
    isempty(columns) && error("No overlapping radial columns found in numerical and reference cross-section files.")

    colors = [:blue, :orange, :green, :red, :purple, :brown, :magenta, :gray40, :olive, :cyan, :black]
    fig = Figure(size = (1000, 550))
    ax = Axis(fig[1, 1], xlabel = "time [d]", ylabel = "temperature [°C]", title = "Analytical ATES validation: temperature vs time")

    for (i, col) in enumerate(columns)
        distance = parse(Int, col[end-2:end])
        color = colors[mod1(i, length(colors))]
        numerical_values = Float64[]
        reference_values = Float64[]
        for d in report_days
            cross = CSV.read(joinpath(outdir, "cross_section_$(d)d.csv"), DataFrame)
            ref = CSV.read(joinpath(ref_dir, "mf6_cross_section_t0$(lpad(string(d), 2, '0'))d.csv"), DataFrame)
            push!(numerical_values, mean(skipmissing(cross[!, col])))
            push!(reference_values, mean(skipmissing(ref[!, col])))
        end
        lines!(ax, report_days, numerical_values, color = color, linewidth = 2, label = "r = $(distance) m")
        lines!(ax, report_days, reference_values, color = color, linewidth = 4, linestyle = :dash, label = i == 1 ? "reference" : nothing)
    end

    axislegend(ax, position = :rb)
    save(joinpath(outdir, "ates_timeseries_comparison.png"), fig)
    return fig
end

function write_vertical_average_time_plot(daily, outdir, reference_dir; distances = 0:10:100)
    ref_daily_path = joinpath(reference_dir, "mf6_results_analytical_T_radial.csv")
    isfile(ref_daily_path) || error("Missing reference file: $(ref_daily_path)")
    ref_daily = CSV.read(ref_daily_path, DataFrame)

    daily_cols = names(daily)
    radial_cols = sort(filter(c -> startswith(c, "T_r"), daily_cols))
    ref_radial_cols = Set(filter(c -> startswith(c, "T_r"), names(ref_daily)))
    columns = [
        radial_column_name(r)
        for r in distances
        if radial_column_name(r) in radial_cols && radial_column_name(r) in ref_radial_cols
    ]
    if isempty(columns)
        columns = [c for c in radial_cols if c in ref_radial_cols]
    end
    isempty(columns) && error("No overlapping radial temperature columns found in daily profile and reference file.")

    colors = [:blue, :orange, :green, :red, :purple, :brown, :magenta, :gray40, :olive, :cyan, :black]
    fig = Figure(size = (1000, 550))
    ax = Axis(
        fig[1, 1],
        xlabel = "time [d]",
        ylabel = "vertically averaged temperature [°C]",
        title = "Vertically averaged radial temperatures vs time",
    )

    for (i, col) in enumerate(columns)
        distance = parse(Int, col[end-2:end])
        color = colors[mod1(i, length(colors))]
        lines!(ax, daily.time_d, daily[!, col], color = color, linewidth = 2, label = "r = $(distance) m")
        lines!(ax, ref_daily.time_d, ref_daily[!, col], color = color, linewidth = 4, linestyle = :dash, label = i == 1 ? "reference" : nothing)
    end

    axislegend(ax, position = :rb)
    save(joinpath(outdir, "ates_vertical_average_temperature_vs_time.png"), fig)
    return fig
end

function run_ates_validation(; outdir = joinpath(pwd(), "ates_validation_outputs"), reference_dir = default_reference_dir(), nang = 32)
    mkpath(outdir)

    case = Fimbul.analytical_ates(nang = nang)
    mesh = physical_representation(reservoir_model(case.model).data_domain)

    sim, cfg = setup_reservoir_simulator(case; info_level = 2, transport_scheme = :tvd)
    result = simulate_reservoir(case, simulator = sim, config = cfg)

    daily = extract_layer_averaged_radial_profile(result.states, mesh)
    CSV.write(joinpath(outdir, "daily_layer_averaged_radial_temperatures.csv"), daily)

    for d in [10, 30, 50, 70, 90]
        cross = extract_cross_section(result.states, mesh, d)
        CSV.write(joinpath(outdir, "cross_section_$(d)d.csv"), cross)
    end

    summary = DataFrame(
        metric = ["n_cells", "n_steps", "n_states", "final_T_r000_C", "final_T_r010_C", "final_T_r050_C"],
        value = [number_of_cells(mesh), length(case.dt), length(result.states), daily.T_r000[end], daily.T_r010[end], daily.T_r050[end]],
    )
    CSV.write(joinpath(outdir, "summary.csv"), summary)

    daily_max, daily_mean, comparison = compare_to_reference(daily, outdir, reference_dir)
    fig = write_timeseries_plot(outdir, reference_dir)
    fig_vertical_avg = write_vertical_average_time_plot(daily, outdir, reference_dir)

    println("Validation outputs written to: $(outdir)")
    println("Reference directory: $(reference_dir)")
    println("Overall max abs diff [C]: $(daily_max)")
    println("Overall mean abs diff [C]: $(daily_mean)")

    return (
        case = case,
        result = result,
        daily = daily,
        comparison = comparison,
        summary = summary,
        figure = fig,
        figure_vertical_average = fig_vertical_avg,
        outdir = outdir,
        reference_dir = reference_dir,
        daily_max = daily_max,
        daily_mean = daily_mean,
    )
end

# Run validation with defaults. Override reference path via:
#   ENV["FIMBUL_ATES_REFERENCE_DIR"] = "<path to Results_analytical>"

GLMakie.closeall()
validation = run_ates_validation(nang = 64)
# validation.figure
validation.figure_vertical_average

##





##
# using GLMakie
# mesh, geometry, well_cells = analytical_ates(nang = 32)

# ##
mesh = physical_representation(reservoir_model(validation.case.model).data_domain)
geo = tpfv_geometry(mesh)
radii = vec(sqrt.(sum((geo.cell_centroids[1:2, :] .+ 0.5).^2, dims = 1)))

##

radii = radial_index_cell_map(mesh)
bc_cells = [bc.cell for bc in validation.case.forces[1][:Reservoir].bc]

mesh = physical_representation(reservoir_model(validation.case.model).data_domain)
plot_mesh(mesh, cells = bc_cells)

##
state = Jutul.evaluate_all_secondary_variables(validation.case.model, validation.case.state0)

##
mesh = physical_representation(reservoir_model(validation.case.model).data_domain)

rix = radial_index_cell_map(mesh)

rix_vec = zeros(number_of_cells(mesh))
for (r, cells) in rix
    for c in cells
        rix_vec[c] = r
    end
end