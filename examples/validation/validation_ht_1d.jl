# # 1D high-temperature geothermal benchmark
# <tags: Validation>
# This example reproduces the 1D high-temperature geothermal benchmark cases
# from [weis_hydrothermal_2014](@cite) using the pressure-enthalpy formulation
# in Fimbul, and validates the results against HYDROTHERM
# [kipp_hydrotherm_2008](@cite) reference solutions.

using Jutul, JutulDarcy, Fimbul, GLMakie

to_celsius(T) = convert_from_si.(T, :Celsius)
to_megapascal(p) = convert_from_si.(p, "megapascal")
to_kj_per_kg(h) = h ./ 1e3

const SINGLE_PHASE_CASES = (:a, :b, :c)
const TWO_PHASE_CASES = (:d, :e)
const CASE_COLORS = Dict(
    :a => :blue,
    :b => :red,
    :c => :lightgreen,
    :d => :cyan,
    :e => :magenta,
)
const PROFILE_SPECS = (
    (name = :Pressure, label = "Pressure [MPa]", transform = to_megapascal,
        hydrotherm_column = "pressure_mpa"),
    (name = :Temperature, label = "Temperature [°C]", transform = to_celsius,
        hydrotherm_column = "temperature_c"),
    (name = :Saturations, label = "Liquid saturation [-]", transform = x -> vec(x[1,:]),
        hydrotherm_column = "liquid_saturation"),
)

nx = 200
cell_size = 10.0

# ## Utilities
function simulate_benchmark_case(case_symbol; vertical = false, nx = 100, cell_size = 10.0)

    case = benchmark_ht_1d(
        benchmark_case = case_symbol,
        nx = nx,
        cell_size = cell_size,
        vertical = vertical,
    )
    case = Fimbul.replace_case_timesteps(
        case,
        Fimbul.load_hydrotherm_1d_timesteps(case_symbol, vertical);
        check_sum = true,
    )
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

function simulate_case_family(case_symbols; nx = 100, cell_size = 10.0)

    outputs = Dict{Tuple{Symbol, Bool}, Any}()
    for case_symbol in case_symbols
        for vertical in Fimbul.available_vertical_modes_ht_1d(case_symbol)
            @info "Simulating case $(case_symbol) with vertical = $(vertical)"
            outputs[(case_symbol, vertical)] = simulate_benchmark_case(
                case_symbol;
                vertical = vertical,
                nx = nx,
                cell_size = cell_size,
            )
        end
    end
    return outputs
end

function reservoir_coordinate(out)
    domain = out.case.model.models[:Reservoir].data_domain
    centroids = tpfv_geometry(physical_representation(domain)).cell_centroids
    if out.case.input_data[:vertical]
        return vec(centroids[3, :])
    else
        return vec(centroids[1, :])
    end
end

function ordered_values(values, vertical)
    values = vec(values)
    return vertical ? reverse(values) : values
end

function plot_property_maps(tables)
    specs = (
        (
            variable = :temperature,
            label = "Temperature [°C]",
        ),
        (
            variable = :density_mix,
            label = "Density [kg/m³]",
        ),
        (
            variable = :saturation_vapor_ph,
            label = "Vapor saturation [-]",
        ),
    )

    fig = Figure(size = (1200, 500))
    for (i, spec) in enumerate(specs)
        yaxisposition = ifelse(i < length(specs), :left, :right)
        ax = Axis(fig[2, i];
            xticksmirrored = true, yticksmirrored = true,
            yaxisposition = yaxisposition, aspect = AxisAspect(1))
        handles = Fimbul.plot_phase_diagram_contours!(
            ax,
            tables;
            variable = spec.variable,
            pressure_limits = (1e5, 52.5e6),
            enthalpy_limits = (500e3, 3500e3),
            n_pressure = 500,
            n_enthalpy = 500,
            levels = 18,
            lines = true,
            contourf_kwargs = (; colormap = :seaborn_icefire_gradient),
        )       
        Colorbar(fig[1, i], handles.filled, vertical = false, flipaxis = true, width = ax.width,
            label = spec.label, labelsize = 20)
        if i > 1
            hideydecorations!(ax, ticks = false)
        end
    end
    is_ax = [f isa Axis for f in fig.content]
    linkaxes!(fig.content[is_ax]...)
    return fig
end

function plot_case_profiles(case_symbol, results)
    vertical_modes = Fimbul.available_vertical_modes_ht_1d(case_symbol)
    fig_height = 400 * length(vertical_modes)
    fig = Figure(size = (800 + 400*(case_symbol ∈ TWO_PHASE_CASES), fig_height))

    for (k, vertical) in enumerate(vertical_modes)
        out = results[(case_symbol, vertical)]
        x = reservoir_coordinate(out)
        state = out.results.states[end]

        row = 2*(k-1)
        vertical_label = vertical ? "vertical" : "horizontal"
        for (col, spec) in enumerate(PROFILE_SPECS)
            if spec.name == :Saturations && case_symbol ∉ TWO_PHASE_CASES
                continue
            end
            values = spec.transform(state[spec.name])
            hydrotherm = Fimbul.load_hydrotherm_1d_property(case_symbol, vertical, spec.hydrotherm_column)
            if vertical
                ax = Axis(
                    fig[row+1, col];
                    xlabel = spec.label,
                    ylabel = "Depth [m]",
                    yreversed = true,
                )
                if hydrotherm !== nothing
                    lines!(ax, hydrotherm.values, hydrotherm.coordinate_m .- cell_size/2;
                        linewidth = 8, linestyle = :dash, color = :black, label = "HYDROTHERM")
                end
                lines!(ax, values, x, color = CASE_COLORS[case_symbol], linewidth = 3, label = "Fimbul")
            else
                ax = Axis(
                    fig[row+1, col];
                    xlabel = "Distance [m]",
                    ylabel = spec.label,
                )
                y = ordered_values(values, vertical)
                if hydrotherm !== nothing
                    lines!(ax, hydrotherm.coordinate_m .- cell_size/2, hydrotherm.values;
                        linewidth = 8, linestyle = :dash, color = :black, label = "HYDROTHERM")
                end
                lines!(ax, x, y, color = CASE_COLORS[case_symbol], linewidth = 3, label = "Fimbul")
                if col == 1
                    axislegend(ax; position = :rt)
                end
            end
        end
        fig[row, :] = Label(fig, "Case $(case_symbol) ($(vertical_label))"; fontsize = 20)
    end
    return fig
end

# ## Simulate all cases
all_results = simulate_case_family((SINGLE_PHASE_CASES..., TWO_PHASE_CASES...); nx = nx, cell_size = cell_size)

# ## H2O properties in pressure-enthalpy space
# The steam tables have been generated using the `CoolProp` library
# [bell_pure_2014](@cite)--see also `FimbulCoolPropExt``). They are available
# in the model's fluid system, but can also be loaded directly using Artifacts.
# To generate your own steam tables, see the `build_steam_tables_h2o` function
# in `FimbulCoolPropExt`.

# We first inspect the steam tables directly. The figure below shows three key
# properties in $(p, h)$-space: temperature, density, and vapor saturation. The
# two-phase envelope is drawn on each subplot.
tables = Fimbul.steam_tables_h2o()
fig_properties = plot_property_maps(tables)
fig_properties

# ## Single-phase benchmark cases
# Cases `:a` to `:c` remain in the single-phase region for the plotted states.
# We therefore use them to compare Fimbul and HYDROTHERM pressure and
# temperature profiles along the 1D column, both with and without gravity.

# ### Case a
# Case `:a` spans 50 to 25 MPa and 350 to 150 °C. The pressure stays high
# enough that the fluid remains in the compressed-liquid region throughout the
# column, so the solution is single-phase liquid with smooth pressure and
# temperature variations.
fig_case_a = plot_case_profiles(:a, all_results)
fig_case_a

# ### Case b
# Case `:b` spans 40 to 20 MPa and 450 to 300 °C. These conditions are hotter
# than case `:a`, but the pressure is still high enough to avoid flashing, so
# the response remains single-phase while showing stronger thermal contrasts.
fig_case_b = plot_case_profiles(:b, all_results)
fig_case_b

# ### Case c
# Case `:c` spans 15 to 1 MPa and 500 to 350 °C. Here the fluid is much hotter
# and less compressed, placing it in the vapor-dominated single-phase region,
# so the profiles represent hot steam rather than liquid water.
fig_case_c = plot_case_profiles(:c, all_results)
fig_case_c

# ## Two-phase benchmark cases
# Cases `:d` and `:e` traverse the two-phase region. Here we compare Fimbul and
# HYDROTHERM pressure, temperature, and liquid-saturation profiles, and then
# compare the no-gravity paths in pressure-enthalpy space.

# ### Case d
# Case `:d` spans 20 to 1 MPa and 400 to 150 °C. The inlet starts as hot,
# pressurized water, but the strong pressure and temperature drop drives the
# state path across the saturation envelope, producing a genuine two-phase
# liquid-vapor transition along the column.
fig_case_d = plot_case_profiles(:d, all_results)
fig_case_d

# ### Case e
# Case `:e` spans 4 to 1 MPa and 300 to 150 °C. Because the pressures are low,
# boiling is easier to trigger than in case `:d`, so the system develops a broad
# two-phase region with more pronounced vapor formation.
fig_case_e = plot_case_profiles(:e, all_results)
fig_case_e

# ## Phase diagram comparison
# We compare the horizontal cases directly in pressure-enthalpy space,
# overlaying Fimbul state paths and HYDROTHERM reference paths on temperature
# contours.
fig = Figure(size = (700, 640))
Label(fig[0, 1:2], "Phase diagram comparison"; fontsize = 22)
ax = Axis(fig[1, 1]; xticksmirrored = true, yticksmirrored = true, aspect = AxisAspect(1))
handles = Fimbul.plot_phase_diagram_contours!(
    ax,
    tables;
    variable = :temperature,
    pressure_limits = (1e5, 52.5e6),
    enthalpy_limits = (500e3, 3500e3),
    levels = 20,
    lines = true,
)

for case_symbol in (SINGLE_PHASE_CASES..., TWO_PHASE_CASES...)
    hydrotherm = Fimbul.load_hydrotherm_1d_phase_path(case_symbol)
    lines!(ax, hydrotherm.enthalpy_kj_per_kg, hydrotherm.pressure_mpa; linewidth = 6, linestyle = :dash, color = CASE_COLORS[case_symbol])
    out = all_results[(case_symbol, false)]
    state = out.results.states[end]
    Fimbul.plot_reservoir_state_ph!(
        ax,
        state;
        color = CASE_COLORS[case_symbol],
        linewidth = 3,
        label = "Case $(case_symbol)",
    )
end

axislegend(ax; position = :rt)
Colorbar(fig[1, 2], handles.filled, label = "Temperature [°C]")
fig
