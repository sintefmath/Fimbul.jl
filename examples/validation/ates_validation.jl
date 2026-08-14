# # Analytical radial benchmark with and without dispersion
# <tags: Validation>
# This example validates the radial thermal benchmark against the analytical
# solution in two configurations:
# 1. Pure advection-conduction without thermal dispersion.
# 2. Advection-conduction with thermal dispersion.
#
# For each configuration, we run a reference simulation, compare numerical and
# analytical temperatures at selected radii, and perform a simple spatial
# convergence study.

using Jutul, JutulDarcy, Fimbul
using HYPRE
using Statistics
using DataFrames
using GLMakie

const DAYS = collect(1:90)

to_celsius(T) = convert_from_si.(T, :Celsius)

# ## Utilities
function radial_shell_lookup(case)
    mesh = physical_representation(reservoir_model(case.model).data_domain)
    geo = tpfv_geometry(mesh)
    x = geo.cell_centroids[1:2, :]
    r_all = vec(sqrt.(sum((x .+ 0.5).^2, dims = 1)))
    r_unique = sort(unique(round.(r_all; digits = 6)))
    cells = [findfirst(isapprox.(r_all, ri; atol = 1e-6)) for ri in r_unique]
    r_map = Dict{Int, Float64}(round(Int, r_unique[i]) => r_unique[i] for i in eachindex(r_unique))
    return cells, r_unique, r_map
end

function numerical_temperature_timeseries(results, cells)
    return [
        [to_celsius(state[:Temperature][cell]) for state in results.states]
        for cell in cells
    ]
end

function analytical_temperature_timeseries(r_real, α_L)
    ϕ = 0.35
    b = 20.0
    ρ_f = 998.0
    ρ_r = 2650.0
    C_r = 800.0
    C_f = 4184.0

    ρ_bC_b = ρ_f*C_f*ϕ + ρ_r*C_r*(1 - ϕ)
    R = ρ_bC_b/(ϕ*ρ_f*C_f)

    λ_f = 0.59
    λ_r = 2.0
    λ_b = ϕ*λ_f + (1 - ϕ)*λ_r
    D = λ_b/(ϕ*ρ_f*C_f)

    Q = 400.0*si"meter^3/day"
    T_0 = 10.0
    T_i = 20.0
    t = DAYS .* si"day"

    return Fimbul.radial_analytical.(r_real, t, Q, T_i, T_0, α_L, b, ϕ, D, R)
end

function plot_radial_time_comparison(case, results; α_L = 0.0, plot_radii = 1:5:40, title = "Radial benchmark")
    cells, r_unique, r_map = radial_shell_lookup(case)
    T_sim_all = numerical_temperature_timeseries(results, cells)
    T_sim = Dict{Int, Vector{Float64}}(round(Int, r_unique[i]) => T_sim_all[i] for i in eachindex(r_unique))

    fig = Figure(size = (1100, 600))
    ax = Axis(
        fig[1, 1],
        xlabel = "time [d]",
        ylabel = "temperature [°C]",
        title = title,
    )

    colors = cgrad(:BrBG_4, length(plot_radii), categorical = true)
    for (k, r_key) in enumerate(plot_radii)
        r = round(Int, r_key)
        if !haskey(T_sim, r) || !haskey(r_map, r)
            continue
        end
        r_real = r_map[r]
        T_a = analytical_temperature_timeseries(r_real, α_L)

        lines!(ax, DAYS, T_sim[r], color = colors[k], linewidth = 2, label = "r=$(round(r_real, digits = 1)) m")
        lines!(ax, DAYS, T_a, color = colors[k], linewidth = 6, linestyle = :dash)
    end
    Legend(fig[1, 2], ax)
    return fig
end

function setup_radial_case(; with_dispersion = false, nrad = 1600, nang = 32, α_L = 1.0, α_T = 0.1)
    thermal_dispersivity = if with_dispersion
        permutedims([α_L α_T])
    else
        nothing
    end

    return Fimbul.analytical_radial(
        nang = nang,
        nrad = nrad,
        thermal_dispersivity = thermal_dispersivity,
        layer_depths = [0.0, 20.0],
    )
end

function simulate_radial_case(case)
    sim, cfg = setup_reservoir_simulator(
        case;
        initial_dt = 1.0,
        info_level = 0,
        max_dt = 0.1si"day",
    )
    tss = VariableChangeTimestepSelector(:Temperature, 0.25, model = :Reservoir, relative = false)
    push!(cfg[:timestep_selectors], tss)
    return simulate_reservoir(case; simulator = sim, config = cfg)
end

function compute_case_error(case, results; α_L = 0.0, plot_radii = 5:5:40)
    cells, r_unique, r_map = radial_shell_lookup(case)
    T_sim_all = numerical_temperature_timeseries(results, cells)
    T_sim = Dict{Int, Vector{Float64}}(round(Int, r_unique[i]) => T_sim_all[i] for i in eachindex(r_unique))

    errs = Float64[]
    for r_key in plot_radii
        r = round(Int, r_key)
        if !haskey(T_sim, r) || !haskey(r_map, r)
            continue
        end
        T_a = analytical_temperature_timeseries(r_map[r], α_L)
        push!(errs, sqrt(mean((T_sim[r] .- T_a).^2)))
    end
    return mean(errs)
end

function convergence_study(; with_dispersion = false, nrad_values = [400, 800, 1600], nang = 32, α_L = 1.0, α_T = 0.1)
    rows = DataFrame(
        nrad = Int[],
        n_cells = Int[],
        h = Float64[],
        rms_error_C = Float64[],
    )

    for nrad in nrad_values
        case = setup_radial_case(
            with_dispersion = with_dispersion,
            nrad = nrad,
            nang = nang,
            α_L = α_L,
            α_T = α_T,
        )
        results = simulate_radial_case(case)
        err = compute_case_error(case, results; α_L = with_dispersion ? α_L : 0.0)
        mesh = physical_representation(reservoir_model(case.model).data_domain)

        push!(rows, (
            nrad,
            number_of_cells(mesh),
            1.0/nrad,
            err,
        ))
    end

    return rows
end

function plot_convergence(rows; title = "Convergence")
    h0 = rows.h ./ rows.h[1]
    e0 = rows.rms_error_C ./ rows.rms_error_C[1]

    fig = Figure(size = (800, 550))
    ax = Axis(
        fig[1, 1],
        xlabel = "h/h0",
        ylabel = "RMS error / error0",
        xscale = log2,
        yscale = log2,
        title = title,
    )
    lines!(ax, h0, e0, linewidth = 2, color = :black, label = "Numerical")
    scatter!(ax, h0, e0, marker = :rect, markersize = 14, color = :black)
    lines!(ax, h0, 0.5*h0.^0.5, linewidth = 6, linestyle = :dash, color = :black, label = "1st-order reference")
    axislegend(ax, position = :rb)
    return fig
end

# ## Part 1: No dispersion
# Baseline radial benchmark with thermal dispersion disabled.
case_no_disp = setup_radial_case(with_dispersion = false, nrad = 800, nang = 32)
results_no_disp = simulate_radial_case(case_no_disp)
fig_no_disp = plot_radial_time_comparison(
    case_no_disp,
    results_no_disp;
    α_L = 0.0,
    title = "Analytical radial benchmark (no dispersion)",
)
fig_no_disp

# ### Simple convergence study (no dispersion)
conv_no_disp = convergence_study(
    with_dispersion = false,
    nrad_values = [50, 100, 200, 400, 800, 1600],
    nang = 32,
)
fig_conv_no_disp = plot_convergence(conv_no_disp; title = "Convergence (no dispersion)")
fig_conv_no_disp
conv_no_disp

# ## Part 2: With dispersion
# Radial benchmark with longitudinal and transverse thermal dispersivity.
αL_disp = 1.0
αT_disp = 0.1
case_disp = setup_radial_case(
    with_dispersion = true,
    nrad = 1600,
    nang = 32,
    α_L = αL_disp,
    α_T = αT_disp,
)
results_disp = simulate_radial_case(case_disp)
fig_disp = plot_radial_time_comparison(
    case_disp,
    results_disp;
    α_L = αL_disp,
    title = "Analytical radial benchmark (with dispersion)",
)
fig_disp

# ### Simple convergence study (with dispersion)
conv_disp = convergence_study(
    with_dispersion = true,
    nrad_values = [400, 800, 1600],
    nang = 32,
    α_L = αL_disp,
    α_T = αT_disp,
)
fig_conv_disp = plot_convergence(conv_disp; title = "Convergence (with dispersion)")
fig_conv_disp
conv_disp