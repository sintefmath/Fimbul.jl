# # 1D two-phase geothermal benchmark
# This example reproduces the 1D geothermal benchmark cases from Weis et al.
# (2014) using the pressure-enthalpy formulation in Fimbul.

using Jutul, JutulDarcy, Fimbul, CoolProp, GLMakie

to_celsius(T) = convert_from_si(T, :Celsius)
to_megapascal(p) = convert_from_si(p, "megapascal")
to_kj_per_kg(h) = h ./ 1e3

const TABLE_PRESSURE_LIMITS = (1e5, 50e6)
const TABLE_ENTHALPY_LIMITS = (500e3, 3500e3)
const SINGLE_PHASE_CASES = (:a, :b, :c)
const TWO_PHASE_CASES = (:d, :e)
const CASE_COLORS = Dict(
    :a => :blue,
    :b => :red,
    :c => :yellow,
    :d => :cyan,
    :e => :magenta,
)
const PROFILE_SPECS = (
    (name = :Pressure, label = "Pressure [MPa]", transform = to_megapascal),
    (name = :Enthalpy, label = "Enthalpy [kJ/kg]", transform = to_kj_per_kg),
    (name = :Temperature, label = "Temperature [C]", transform = to_celsius),
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
table_resolution = 500
nx = 100
cell_size = 10.0

tables = Fimbul.build_steam_tables_2ph(
    n_pressure = table_resolution,
    n_enthalpy = table_resolution,
)

##
available_gravity_modes(case_symbol) = case_symbol == :e ? (false,) : (false, true)

function trim_case_to_plot_time(case)
    plot_time = get(case.input_data, :plot_time, sum(case.dt))
    cumulative_time = cumsum(case.dt)
    nstep = something(findfirst(t -> t >= plot_time, cumulative_time), length(case.dt))
    return case[1:nstep]
end

function simulate_benchmark_case(case_symbol, tables; gravity = false, nx = 100, cell_size = 10.0)
    case = benchmark_2ph_1d(
        benchmark_case = case_symbol,
        nx = nx,
        cell_size = cell_size,
        enthalpy_tables = tables,
        gravity = gravity,
    )
    case = trim_case_to_plot_time(case)
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
        for gravity in available_gravity_modes(case_symbol)
            @info "Simulating case $(case_symbol) with gravity = $(gravity)"
            outputs[(case_symbol, gravity)] = simulate_benchmark_case(
                case_symbol,
                tables;
                gravity = gravity,
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
    if out.case.input_data[:gravity]
        return vec(centroids[3, :])
    else
        return vec(centroids[1, :])
    end
end

function ordered_values(values, gravity)
    values = vec(values)
    return gravity ? reverse(values) : values
end

function plot_property_maps(tables)
    specs = (
        (
            variable = :temperature,
            title = "Temperature",
            label = "Temperature [C]",
            transform = nothing,
            colormap = :seaborn_icefire_gradient,
        ),
        (
            variable = :density_mix,
            title = "Density",
            label = "Density [kg/m^3]",
            transform = nothing,
            colormap = :haline,
        ),
        (
            variable = :saturation_vapor_ph,
            title = "Vapor saturation",
            label = "Vapor saturation [-]",
            transform = nothing,
            colormap = :viridis,
        ),
    )

    fig = Figure(size = (1800, 700), fontsize = 18)
    for (i, spec) in enumerate(specs)
        ax = Axis(fig[1, i]; title = spec.title)
        _, hcf, hcl = Fimbul.plot_phase_diagram_contours!(
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
        Colorbar(fig[2, i], hcf, vertical = false, flipaxis = false, label = spec.label)
        if i > 1
            hideydecorations!(ax, ticks = false)
        end
    end
    return fig
end

function plot_case_profiles(case_symbol, results)
    gravity_modes = available_gravity_modes(case_symbol)
    fig = Figure(size = (1200, 350 * length(gravity_modes)), fontsize = 18)

    for (row, gravity) in enumerate(gravity_modes)
        out = results[(case_symbol, gravity)]
        x = reservoir_coordinate(out)
        state = out.results.states[plotting_timestep(out)]
        x_label = gravity ? "Depth [m]" : "Distance [m]"
        gravity_label = gravity ? "with gravity" : "without gravity"

        for (col, spec) in enumerate(PROFILE_SPECS)
            ax = Axis(
                fig[row, col];
                title = "$(String(spec.name)) ($(gravity_label))",
                xlabel = x_label,
                ylabel = spec.label,
            )
            y = ordered_values(spec.transform.(state[spec.name]), gravity)
            lines!(ax, x, y, color = CASE_COLORS[case_symbol], linewidth = 3)
        end
    end
    return fig
end

function plot_case_family_diagram(case_symbols, results, tables; title)
    fig = Figure(size = (700, 600), fontsize = 18)
    ax = Axis(fig[1, 1]; title = title)
    _, hcf, hcl = Fimbul.plot_phase_diagram_contours!(
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
    axislegend(ax; position = :rt, framevisible = false)
    Colorbar(fig[1, 2], hcf, label = "Temperature [C]")
    return fig
end
##
all_results = simulate_case_family((SINGLE_PHASE_CASES..., TWO_PHASE_CASES...), tables; nx = nx, cell_size = cell_size)

# ## H2O properties in pressure-enthalpy space
# We first inspect the steam tables directly. The figure below shows four key
# properties in $(p, h)$-space: temperature, density, a saturation-weighted
# effective viscosity, and vapor saturation. The two-phase envelope is drawn on
# each subplot.
fig_properties = plot_property_maps(tables)
fig_properties

# ## Single-phase benchmark cases
# Cases `:a` to `:c` remain in the single-phase region for the plotted states.
# We therefore use them to inspect how pressure, enthalpy, and temperature vary
# along the 1D column both with and without gravity.

# ### Case a
fig_case_a = plot_case_profiles(:a, all_results)
fig_case_a

# ### Case b
fig_case_b = plot_case_profiles(:b, all_results)
fig_case_b

# ### Case c
fig_case_c = plot_case_profiles(:c, all_results)
fig_case_c

# ###
# The three single-phase cases without gravity can also be compared directly in
# pressure-enthalpy space. Here their state paths are overlaid on temperature
# contours.
fig_single_phase_diagram = plot_case_family_diagram(
    SINGLE_PHASE_CASES,
    all_results,
    tables;
    title = "Single-phase cases without gravity",
)
fig_single_phase_diagram

# ## Two-phase benchmark cases
# Cases `:d` and `:e` traverse the two-phase region. As above, we first inspect
# the 1D pressure, enthalpy, and temperature profiles, and then compare the
# no-gravity paths in pressure-enthalpy space.

# ### Case d
fig_case_d = plot_case_profiles(:d, all_results)
fig_case_d

# ### Case e
fig_case_e = plot_case_profiles(:e, all_results)
fig_case_e

fig_two_phase_diagram = plot_case_family_diagram(
    TWO_PHASE_CASES,
    all_results,
    tables;
    title = "Two-phase cases without gravity",
)
fig_two_phase_diagram