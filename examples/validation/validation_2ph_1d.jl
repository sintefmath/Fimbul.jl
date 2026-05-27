# # 1D two-phase geothermal benchmark
# This example reproduces the 1D geothermal benchmark cases from
# [weis_hydrothermal_2014](@cite) using the pressure-enthalpy formulation in
# Fimbul.

using Jutul, JutulDarcy, Fimbul, CoolProp, GLMakie

to_celsius(T) = convert_from_si.(T, :Celsius)
to_megapascal(p) = convert_from_si.(p, "megapascal")
to_kj_per_kg(h) = h ./ 1e3

const TABLE_PRESSURE_LIMITS = (1e5, 52.5e6)
const TABLE_ENTHALPY_LIMITS = (500e3, 3500e3)
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

# ## Governing equations for two-phase flow
# In the two-phase H2O benchmark, the conserved quantities are total mass and
# total energy. Using pressure $p$ and specific enthalpy $h$ as primary
# variables, the governing equations can be written as
#
# ``\frac{\partial}{\partial t}\left(\phi \sum_{\alpha \in \{l, v\}} \rho_\alpha S_\alpha\right)
# + \nabla \cdot \left(\sum_{\alpha \in \{l, v\}} \rho_\alpha \mathbf{u}_\alpha\right) = q_m``
#
# ``\frac{\partial}{\partial t}\left(\phi \sum_{\alpha \in \{l, v\}} \rho_\alpha S_\alpha U_\alpha + (1 - \phi) \rho_r c_r T\right)
# + \nabla \cdot \left(\sum_{\alpha \in \{l, v\}} \rho_\alpha h_\alpha \mathbf{u}_\alpha - \lambda \nabla T\right) = q_e``
#
# where phase velocities follow Darcy's law,
#
# ``\mathbf{u}_\alpha = -\frac{k k_{r,\alpha}}{\mu_\alpha}\left(\nabla p - \rho_\alpha g \nabla z\right)``
#
# The nonlinear closure is provided by steam-table lookups that map $(p, h)$ to
# temperature, phase densities, phase viscosities, phase enthalpies, and vapor
# saturation.

# Shared setup used throughout the example.
table_resolution = 50
nx = 200
cell_size = 10.0

##
tables = Fimbul.build_steam_tables_2ph(
    n_pressure = table_resolution,
    n_enthalpy = table_resolution,
    info_level = 1,
)

##
available_vertical_modes(case_symbol) = case_symbol == :e ? (false,) : (false, true)

hydrotherm_profile_path(case_symbol, vertical; root = HYDROTHERM_RESULTS_ROOT) =
    joinpath(root, "case_$(case_symbol)_$(vertical ? "vertical" : "horizontal").txt")

function load_hydrotherm_table(case_symbol, vertical; root = HYDROTHERM_RESULTS_ROOT)
    path = hydrotherm_profile_path(case_symbol, vertical; root)
    isfile(path) || return nothing

    lines = readlines(path)
    isempty(lines) && error("Empty HYDROTHERM profile file: $(path)")
    header = split(strip(replace(first(lines), "#" => "")))
    columns = Dict(name => Float64[] for name in header)

    for line in Iterators.drop(lines, 1)
        stripped = strip(line)
        isempty(stripped) && continue
        tokens = split(stripped)
        length(tokens) == length(header) || error("Expected $(length(header)) columns in $(path), got $(length(tokens))")
        for (name, token) in zip(header, tokens)
            push!(columns[name], token == "NaN" ? NaN : parse(Float64, token))
        end
    end

    return columns
end

function load_hydrotherm_property(case_symbol, vertical, spec)
    table = load_hydrotherm_table(case_symbol, vertical)
    table === nothing && return nothing
    haskey(table, spec.hydrotherm_column) || return nothing

    coordinate_column = vertical ? "depth_m" : "distance_m"
    return (coordinate_m = table[coordinate_column], values = table[spec.hydrotherm_column])
end

function load_hydrotherm_phase_path(case_symbol)
    table = load_hydrotherm_table(case_symbol, false)
    table === nothing && return nothing
    return (pressure_mpa = table["pressure_mpa"], enthalpy_kj_per_kg = table["enthalpy_kj_per_kg"])
end

function simulate_benchmark_case(case_symbol, tables; vertical = false, nx = 100, cell_size = 10.0)
    case = benchmark_2ph_1d(
        benchmark_case = case_symbol,
        nx = nx,
        cell_size = cell_size,
        enthalpy_tables = tables,
        vertical = vertical,
    )
    simulator, config = setup_reservoir_simulator(
        case;
        tol_cnv = 1e-3,
        tol_mb = 1e-7,
        # info_level = 1,
        # max_timestep = maximum(case.dt),
        max_timestep = 0.5si"year",
        relaxation = true,
    )
    results = simulate_reservoir(case; simulator = simulator, config = config)
    return (case = case, results = results)
end

function simulate_case_family(case_symbols, tables; nx = 100, cell_size = 10.0)
    outputs = Dict{Tuple{Symbol, Bool}, Any}()
    for case_symbol in case_symbols
        for vertical in available_vertical_modes(case_symbol)
            @info "Simulating case $(case_symbol) with vertical = $(vertical)"
            outputs[(case_symbol, vertical)] = simulate_benchmark_case(
                case_symbol,
                tables;
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
            pressure_limits = TABLE_PRESSURE_LIMITS,
            enthalpy_limits = TABLE_ENTHALPY_LIMITS,
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
    vertical_modes = available_vertical_modes(case_symbol)
    hydrotherm_present = true
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
            hydrotherm = load_hydrotherm_property(case_symbol, vertical, spec)
            hydrotherm_present = hydrotherm_present || hydrotherm !== nothing
            if vertical
                ax = Axis(
                    fig[row+1, col];
                    xlabel = spec.label,
                    ylabel = "Depth [m]",
                    yreversed = true,
                )
                if hydrotherm !== nothing
                    lines!(ax, hydrotherm.values, hydrotherm.coordinate_m;
                        linewidth = 6, linestyle = :dash, color = :black)
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
                    lines!(ax, hydrotherm.coordinate_m, hydrotherm.values;
                        linewidth = 6, linestyle = :dash, color = :black, label = "HYDROTHERM")
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

##
# all_results = simulate_case_family((SINGLE_PHASE_CASES..., TWO_PHASE_CASES...), tables; nx = nx, cell_size = cell_size)
# all_results = simulate_case_family((SINGLE_PHASE_CASES..., TWO_PHASE_CASES...), tables; nx = nx, cell_size = cell_size)
all_results = simulate_case_family((SINGLE_PHASE_CASES..., TWO_PHASE_CASES...), tables; nx = nx, cell_size = cell_size)

# ## H2O properties in pressure-enthalpy space
# We first inspect the steam tables directly. The figure below shows three key
# properties in $(p, h)$-space: temperature, density, and vapor saturation. The
# two-phase envelope is drawn on each subplot.
fig_properties = plot_property_maps(tables)
fig_properties

# ## Single-phase benchmark cases
# Cases `:a` to `:c` remain in the single-phase region for the plotted states.
# We therefore use them to inspect how pressure, enthalpy, and temperature vary
# along the 1D column both with and without gravity.

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
# Cases `:d` and `:e` traverse the two-phase region. As above, we first inspect
# the 1D pressure, enthalpy, and temperature profiles, and then compare the
# no-gravity paths in pressure-enthalpy space.

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
# We compare the cases directly in pressure-enthalpy space, overlaying the state
# paths on temperature contours.
fig = Figure(size = (700, 640))
Label(fig[0, 1:2], "Phase diagram comparison"; fontsize = 22)
ax = Axis(fig[1, 1]; xticksmirrored = true, yticksmirrored = true, aspect = AxisAspect(1))
handles = Fimbul.plot_phase_diagram_contours!(
    ax,
    tables;
    variable = :temperature,
    pressure_limits = TABLE_PRESSURE_LIMITS,
    enthalpy_limits = TABLE_ENTHALPY_LIMITS,
    levels = 20,
    lines = true,
)

for case_symbol in (SINGLE_PHASE_CASES..., TWO_PHASE_CASES...)
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
