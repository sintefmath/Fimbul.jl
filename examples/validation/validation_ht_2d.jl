# # 2D high-temperature geothermal benchmark
# This example demonstrates the 2D high-temperature geothermal benchmark from
# [weis_hydrothermal_2014](@cite) using the pressure-enthalpy formulation in
# Fimbul.
#
# We demonstrate the currently implemented magmatic-fluid-source variants: a
# moderate-enthalpy source that remains single-phase near the injector, and a
# high-enthalpy source that produces a two-phase plume during ascent. The model
# represents a 9 km by 3 km vertical crustal section with an open top boundary
# at atmospheric pressure and 10 °C, and a hot H2O source at the center of the
# bottom boundary. The alternative bottom heat-flux benchmark is not included
# here because heat-flux boundary conditions are not yet supported in Fimbul.

using Jutul, JutulDarcy, Fimbul, HYPRE, GLMakie

to_celsius(T) = convert_from_si.(T, :Celsius)
to_megapascal(p) = convert_from_si.(p, "megapascal")

nx = 91
nz = 30

const HYDROTHERM_RESULTS_ROOT = normpath(joinpath(
    @__DIR__, "..", "..", "..", "..", "..", "misc", "hydrotherm-dev", "validation_ht_2d", "profiles"
))
const X_LIMITS_KM = (-2.0, 2.0)
const SNAPSHOT_YEARS = (500.0, 2000.0, 5000.0)
const SINGLE_PHASE_TEMPERATURE_LEVELS = collect(0.0:25.0:125.0)
const TWO_PHASE_TEMPERATURE_LEVELS = collect(0.0:50.0:350.0)
const SINGLE_PHASE_PRESSURE_LEVELS = collect(0.0:5.0:50.0)
const TWO_PHASE_PRESSURE_LEVELS = collect(0.0:5.0:35.0)
const VAPOR_SATURATION_LEVELS = [0.0, 1e-6, 0.2, 0.4, 0.6, 0.8, 1.0]
const HYDROTHERM_CONTOUR_LINEWIDTH = 5
const HYDROTHERM_CONTOUR_COLOR = :black

##
tables = Fimbul.build_steam_tables_2ph()

# ## Set up and simulate the benchmark cases
# We simulate both fluid-source variants from Weis et al. (2014): the
# moderate-enthalpy single-phase plume and the hotter two-phase plume.
function simulate_benchmark_case(benchmark_case, tables; nx = 120, nz = 60)
    case = benchmark_ht_2d(
        benchmark_case = benchmark_case,
        nx = nx,
        nz = nz,
        enthalpy_tables = tables,
    )
    case = replace_case_timesteps(case, load_hydrotherm_timesteps(case))

    simulator, config = setup_reservoir_simulator(
        case;
        tol_cnv = 1e-3,
        tol_mb = 1e-7,
        max_timestep = Inf,
        timesteps = :none,
        relaxation = true,
    )
    results = simulate_reservoir(case; simulator = simulator, config = config)
    return (case = case, results = results)
end

single_phase = simulate_benchmark_case(:single_phase_source, tables; nx = nx, nz = nz)
two_phase = simulate_benchmark_case(:two_phase_source, tables; nx = nx, nz = nz)

##
function section_axes(case)
    domain = case.model.models[:Reservoir].data_domain
    centroids = tpfv_geometry(physical_representation(domain)).cell_centroids
    x0 = case.input_data[:x_coordinate_origin]
    x = sort(unique(vec(centroids[1, :] .- x0))) ./ 1e3
    depth = sort(unique(vec(centroids[3, :]))) ./ 1e3
    return (x_km = x, depth_km = depth)
end

function section_data(case, values)
    axes = section_axes(case)
    nx = length(axes.x_km)
    nz = length(axes.depth_km)
    return reshape(vec(values), nx, nz)
end

function hydrotherm_property_path(case_name, property_name; root = HYDROTHERM_RESULTS_ROOT)
    path = joinpath(root, "$(case_name)_$(property_name).txt")
    isfile(path) || error("Missing HYDROTHERM export $(path). Run export_profiles.jl in validation_ht_2d first.")
    return path
end

function hydrotherm_timestep_path(case_name; root = HYDROTHERM_RESULTS_ROOT)
    path = joinpath(root, "$(case_name)_timesteps.txt")
    isfile(path) || error("Missing HYDROTHERM timestep export $(path). Run export_profiles.jl in validation_ht_2d first.")
    return path
end

function load_hydrotherm_vector(path)
    values = Float64[]
    for line in eachline(path)
        stripped = strip(line)
        (isempty(stripped) || startswith(stripped, "#")) && continue
        push!(values, parse(Float64, stripped))
    end
    return values
end

function load_hydrotherm_timesteps(case)
    timesteps = Float64[]
    for line in eachline(hydrotherm_timestep_path(hydrotherm_case_name(case)))
        stripped = strip(line)
        (isempty(stripped) || startswith(stripped, "#")) && continue
        tokens = split(stripped)
        length(tokens) >= 2 || error("Expected at least two columns in HYDROTHERM timestep line: $(line)")
        push!(timesteps, parse(Float64, tokens[2]))
    end
    return timesteps .* si_unit(:year)
end

function replace_case_timesteps(case, timesteps)
    isapprox(sum(timesteps), sum(case.dt); rtol = 1e-3) ||
        error("HYDROTHERM timestep sum $(sum(timesteps)) does not match case duration $(sum(case.dt))")

    (; model, forces, state0, parameters, input_data) = case
    return JutulCase(model, timesteps, forces;
        state0 = state0,
        parameters = parameters,
        input_data = input_data,
    )
end

function load_hydrotherm_matrix(path)
    rows = Vector{Vector{Float64}}()
    for line in eachline(path)
        stripped = strip(line)
        (isempty(stripped) || startswith(stripped, "#")) && continue
        push!(rows, parse.(Float64, split(stripped)))
    end
    isempty(rows) && error("Empty HYDROTHERM matrix file: $(path)")
    matrix = reduce(vcat, permutedims.(rows))
    return permutedims(matrix)
end

function align_hydrotherm_depth_to_fimbul(z_m, values...)
    length(z_m) >= 2 || error("Need at least two HYDROTHERM depth coordinates to align the reference grid")
    dz_m = z_m[2] - z_m[1]
    z_aligned_m = z_m[2:end] .- dz_m
    aligned_values = map(values) do value
        value[:, 2:end]
    end
    return (z_aligned_m, aligned_values...)
end

hydrotherm_case_name(case) = String(case.input_data[:benchmark_case])

function load_hydrotherm_reference(case)
    case_name = hydrotherm_case_name(case)
    x_m = load_hydrotherm_vector(hydrotherm_property_path(case_name, "x_m"))
    z_m = load_hydrotherm_vector(hydrotherm_property_path(case_name, "z_m"))
    pressure_mpa = load_hydrotherm_matrix(hydrotherm_property_path(case_name, "pressure_mpa"))
    temperature_c = load_hydrotherm_matrix(hydrotherm_property_path(case_name, "temperature_c"))
    liquid_saturation = load_hydrotherm_matrix(hydrotherm_property_path(case_name, "liquid_saturation"))
    vapor_saturation = map(liquid_saturation) do value
        isnan(value) ? NaN : 1.0 - value
    end
    # z_m = z_m .+ (z_m[2] - z_m[1])/2  # shift from HYDROTHERM cell centers to Fimbul cell centers
    # z_m, pressure_mpa, temperature_c, vapor_saturation = align_hydrotherm_depth_to_fimbul(
    #     z_m,
    #     pressure_mpa,
    #     temperature_c,
    #     vapor_saturation,
    # )

    x0 = x_m[cld(length(x_m), 2)]
    return (
        x_km = (x_m .- x0) ./ 1e3,
        depth_km = z_m ./ 1e3,
        pressure = pressure_mpa,
        temperature = temperature_c,
        vapor_saturation = vapor_saturation,
    )
end

function snapshot_index(case, years)
    times_years = convert_from_si.(cumsum(case.dt), :year)
    return argmin(abs.(times_years .- years))
end

to_vapor_saturation(S) = vec(S[2, :])

final_time_years(case) = round(Int, convert_from_si(sum(case.dt), :year))

function source_location_km(case)
    domain = case.model.models[:Reservoir].data_domain
    centroids = tpfv_geometry(physical_representation(domain)).cell_centroids
    source_cell = case.input_data[:source_cell]
    x0 = case.input_data[:x_coordinate_origin]
    return (
        x_km = (centroids[1, source_cell] - x0) / 1e3,
        depth_km = centroids[3, source_cell] / 1e3,
    )
end


function contour_colorbar_ticks(levels)
    ticks = collect(levels)
    labels = map(ticks) do value
        rounded = round(value; digits = 1)
        if isapprox(rounded, round(rounded))
            string(Int(round(rounded)))
        else
            string(rounded)
        end
    end
    return (ticks, labels)
end

function final_state_field_specs(case, results)
    axes = section_axes(case)
    source = source_location_km(case)
    state = results.states[end]
    final_years = final_time_years(case)
    source_regime = get(case.input_data, :source_regime, :single_phase)
    hydrotherm = load_hydrotherm_reference(case)

    temperature = to_celsius(section_data(case, state[:Temperature]))
    pressure = to_megapascal(section_data(case, state[:Pressure]))

    temperature_levels = source_regime == :two_phase ? TWO_PHASE_TEMPERATURE_LEVELS : SINGLE_PHASE_TEMPERATURE_LEVELS
    pressure_levels = source_regime == :two_phase ? TWO_PHASE_PRESSURE_LEVELS : SINGLE_PHASE_PRESSURE_LEVELS

    specs = Any[
        (
            title = "Temperature after $(final_years) years",
            values = temperature,
            hydrotherm_values = hydrotherm.temperature,
            levels = temperature_levels,
            colormap = :seaborn_icefire_gradient,
            colorbar_label = "Temperature [°C]",
        ),
        (
            title = "Pressure after $(final_years) years",
            values = pressure,
            hydrotherm_values = hydrotherm.pressure,
            levels = pressure_levels,
            colormap = :vik,
            colorbar_label = "Pressure [MPa]",
        ),
    ]

    if source_regime == :two_phase
        vapor_saturation = clamp.(section_data(case, to_vapor_saturation(state[:Saturations])), 0.0, 1.0)
        push!(specs,
            (
                title = "Vapor saturation after $(final_years) years",
                values = vapor_saturation,
                hydrotherm_values = hydrotherm.vapor_saturation,
                levels = VAPOR_SATURATION_LEVELS,
                colormap = :dense,
                colorbar_label = "Vapor saturation [-]",
            ),
        )
    end

    return (axes = axes, source = source, specs = specs, hydrotherm = hydrotherm)
end

function plot_final_state(case, results)
    state = final_state_field_specs(case, results)
    axes = state.axes
    source = state.source
    specs = state.specs
    hydrotherm = state.hydrotherm

    fig = Figure(size = (520*length(specs), 520))
    for (i, spec) in enumerate(specs)
        ax = Axis(
            fig[1, i];
            limits = ((-2.0, 2.0), (0.0, 3.0)),
            title = spec.title,
            xlabel = "Distance [km]",
            ylabel = ifelse(i == 1, "Depth [km]", ""),
            yreversed = true,
            aspect = AxisAspect(4/3),
        )
        plt = contourf!(ax, axes.x_km, axes.depth_km, spec.values;
            colormap = spec.colormap,
            levels = spec.levels,
        )
        contour!(ax, axes.x_km, axes.depth_km, spec.values;
            color = :white,
            levels = spec.levels,
        )
        contour!(ax, hydrotherm.x_km, hydrotherm.depth_km, spec.hydrotherm_values;
            levels = spec.levels,
            color = :white,
            linewidth = HYDROTHERM_CONTOUR_LINEWIDTH,
            linestyle = :dash,
        )
        scatter!(ax, [source.x_km], [source.depth_km]; color = :black, marker = :star5, markersize = 14)
        if i > 1
            hideydecorations!(ax, ticks = false)
        end
        Colorbar(fig[2, i], plt;
            vertical = false,
            label = spec.colorbar_label,
            ticks = contour_colorbar_ticks(spec.levels),
            flipaxis = false,
        )
    end
    return fig
end

# ## Single-phase source case
# The moderate-enthalpy source generates an upward-rising thermal plume that
# stays in the liquid regime. At late time, the deepest part of the plume also
# carries the highest absolute pressures near the source.
plot_reservoir(single_phase.case, single_phase.results.states;
    key = :Temperature, colormap = :seaborn_icefire_gradient, aspect = (9,0.1,3))

# ### Validate against reference data
fig_single_phase_final = plot_final_state(single_phase.case, single_phase.results)
fig_single_phase_final

# ## Two-phase source case
# Raising the source enthalpy to 1.5 MJ/kg produces a hotter plume that enters
# the two-phase field during ascent. At late time, the plume maintains high
# deep pressures near the source and develops a liquid-saturation deficit where
# vapor forms in the rising core.
plot_reservoir(two_phase.case, two_phase.results.states;
    key = :Temperature, colormap = :seaborn_icefire_gradient, aspect = (9,0.1,3))

# ### Validate against reference data
fig_two_phase = plot_final_state(two_phase.case, two_phase.results)
fig_two_phase
