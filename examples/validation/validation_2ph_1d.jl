# # 1D two-phase geothermal benchmark
# This example reproduces the 1D geothermal benchmark cases from
# [weis_hydrothermal_2014](@cite) using the pressure-enthalpy formulation in
# Fimbul.

using Jutul, JutulDarcy, Fimbul, CoolProp, GLMakie

to_celsius(T) = convert_from_si.(T, :Celsius)
to_megapascal(p) = convert_from_si.(p, "megapascal")
to_kj_per_kg(h) = h ./ 1e3

const TABLE_PRESSURE_LIMITS = (1e5, 50e6)
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
    (name = :Pressure, label = "Pressure [MPa]", transform = to_megapascal),
    (name = :Temperature, label = "Temperature [°C]", transform = to_celsius),
    (name = :Saturations, label = "Liquid saturation [-]", transform = x -> vec(x[1,:])),
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
table_resolution = 100
nx = 100
cell_size = 10.0

tables = Fimbul.build_steam_tables_2ph(
    n_pressure = table_resolution,
    n_enthalpy = table_resolution,
)

##
available_vertical_modes(case_symbol) = case_symbol == :e ? (false,) : (false, true)

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
        info_level = 0,
        max_timestep = maximum(case.dt),
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

function plotting_timestep(out; time = nothing)
    case = out.case
    results = out.results
    plot_time = isnothing(time) ? get(case.input_data, :plot_time, last(results.time)) : time
    timestep = findfirst(t -> t >= plot_time, results.time)
    return isnothing(timestep) ? length(results.states) : timestep
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
            title = "Temperature",
            label = "Temperature [°C]",
            transform = nothing,
        ),
        (
            variable = :density_mix,
            title = "Density",
            label = "Density [kg/m³]",
            transform = nothing,
        ),
        (
            variable = :saturation_vapor_ph,
            title = "Vapor saturation",
            label = "Vapor saturation [-]",
            transform = nothing,
        ),
    )

    fig = Figure(size = (1200, 500))
    for (i, spec) in enumerate(specs)
        ax = Axis(fig[2, i])
        handles = Fimbul.plot_phase_diagram_contours!(
            ax,
            tables;
            variable = spec.variable,
            pressure_limits = TABLE_PRESSURE_LIMITS,
            enthalpy_limits = TABLE_ENTHALPY_LIMITS,
            levels = 18,
            lines = true,
            value_transform = spec.transform,
            contourf_kwargs = (; colormap = :seaborn_icefire_gradient),
        )
        Colorbar(fig[1, i], handles.filled, vertical = false, flipaxis = true, label = spec.label, labelsize = 20)
        if i > 1
            hideydecorations!(ax, ticks = false)
        end
    end
    return fig
end

function plot_case_profiles(case_symbol, results)
    vertical_modes = available_vertical_modes(case_symbol)
    fig = Figure(size = (800 + 400*(case_symbol ∈ TWO_PHASE_CASES), 400 * length(vertical_modes)))

    for (k, vertical) in enumerate(vertical_modes)
        out = results[(case_symbol, vertical)]
        x = reservoir_coordinate(out)
        state = out.results.states[plotting_timestep(out)]

        row = 2*(k-1)
        vertical_label = vertical ? "vertical" : "horizontal"
        for (col, spec) in enumerate(PROFILE_SPECS)
            if spec.name == :Saturations && case_symbol ∉ TWO_PHASE_CASES
                continue
            end
            values = spec.transform(state[spec.name])
            if vertical
                ax = Axis(
                    fig[row+1, col];
                    xlabel = spec.label,
                    ylabel = "Depth [m]",
                    yreversed = true,
                )
                lines!(ax, values, x, color = CASE_COLORS[case_symbol], linewidth = 3)
            else
                ax = Axis(
                    fig[row+1, col];
                    xlabel = "Distance [m]",
                    ylabel = spec.label,
                )
                y = ordered_values(values, vertical)
                lines!(ax, x, y, color = CASE_COLORS[case_symbol], linewidth = 3)
            end
        end
        fig[row, :] = Label(fig, "Case $(case_symbol) ($(vertical_label))"; fontsize = 20)
    end
    return fig
end

function plot_case_family_diagram(case_symbols, results, tables; title = "Phase diagram comparison")
    fig = Figure(size = (700, 640))
    Label(fig[0, 1:2], title; fontsize = 22)
    ax = Axis(fig[1, 1])
    handles = Fimbul.plot_phase_diagram_contours!(
        ax,
        tables;
        variable = :temperature,
        pressure_limits = TABLE_PRESSURE_LIMITS,
        enthalpy_limits = TABLE_ENTHALPY_LIMITS,
        levels = 20,
        lines = true,
    )

    for case_symbol in case_symbols
        Fimbul.plot_reservoir_state_ph!(
            ax,
            results[(case_symbol, false)];
            color = CASE_COLORS[case_symbol],
            linewidth = 3,
            label = "Case $(case_symbol)",
        )
    end

    # Legend(fig[1, 2], framevisible = false)
    axislegend(ax; position = :rt)
    Colorbar(fig[1, 2], handles.filled, label = "Temperature [°C]")
    return fig
end
##
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
fig_two_phase_diagram = plot_case_family_diagram(
    vcat(SINGLE_PHASE_CASES..., TWO_PHASE_CASES...),
    all_results,
    tables;
)
fig_two_phase_diagram