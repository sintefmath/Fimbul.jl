# # Deep coaxial geothermal well
# This example demonstrates simulation and analysis of geothermal energy
# production from a deep coaxial closed-loop well. The well is defined by a
# general trajectory (m×3 matrix) and uses coaxial heat exchange with the
# surrounding rock formation.
#
# We explore three scenarios:
# 1. **Homogeneous reservoir**: A single uniform layer to establish a baseline.
# 2. **Layered reservoir**: Four distinct geological layers (styrofoam, clay,
#    sandstone, granite) with exaggerated thermal conductivities to highlight
#    the effect of heterogeneity.
# 3. **Inner pipe conductivity study**: Homogeneous reservoir with varying
#    inner pipe wall thermal conductivity (from 0 to 2× default) to
#    demonstrate heat loss between supply and return pipes.

# Add required modules to namespace
using Jutul, JutulDarcy, Fimbul
using HYPRE
using GLMakie

# Useful SI units
meter, hour, watt, Kelvin, joule, kilogram = si_units(
    :meter, :hour, :watt, :Kelvin, :joule, :kilogram);

# ## Common parameters
# All scenarios share the same well trajectory, flow rate, injection
# temperature, simulation duration, and mesh settings.
well_trajectory = [
    0.0  0.0     0.0;
    0.0  0.0  2500.0;
]
base_args = (
    well_trajectory = well_trajectory,
    rate = 25meter^3/hour,
    temperature_inj = convert_to_si(25.0, :Celsius),
    num_years = 10,
    report_interval = si_unit(:year),
    hxy_min = 2.5,
    mesh_args = (offset = 250.0, offset_rel = missing),
);

# ## Simulator helper
# A reusable function to run any case with consistent solver settings.
function run_case(case)
    sim, cfg = setup_reservoir_simulator(case;
        output_substates = true,
        info_level = 0,
        initial_dt = 120.0,
        presolve_wells = true,
        relaxation = true)
    sel = VariableChangeTimestepSelector(:Temperature, 5.0;
        relative = false, model = :Reservoir)
    push!(cfg[:timestep_selectors], sel)
    sel = VariableChangeTimestepSelector(:Temperature, 5.0;
        relative = false, model = :CoaxialWell_supply)
    push!(cfg[:timestep_selectors], sel)
    return simulate_reservoir(case; simulator = sim, config = cfg)
end

# ---
# ## Scenario 1 – Homogeneous reservoir
# A single layer from the surface to 3000 m with uniform rock properties.
# This serves as a baseline for comparison.
homogeneous_args = (;
    base_args...,
    depths = [0.0, 2550.0, 3000.0],
    permeability = [1e-2, 1e-2]*si_unit(:darcy),
    porosity = [0.01, 0.01],
    rock_thermal_conductivity = [2.5, 2.5]*watt/(meter*Kelvin),
    rock_heat_capacity = [900, 900]*joule/(kilogram*Kelvin),
    rock_density = [2600, 2600]*kilogram/meter^3,
);

case_hom_inner = coaxial_well_branches(; inject_into = :inner, homogeneous_args...);
case_hom_outer = coaxial_well_branches(; inject_into = :outer, homogeneous_args...);

# ### Inspect mesh and well
msh_homo = physical_representation(reservoir_model(case_hom_inner.model).data_domain)
fig = Figure(size = (800, 800))
ax = Axis3(fig[1, 1]; zreversed = true, perspectiveness = 0.5, aspect = (1, 1, 4),
    title = "Homogeneous reservoir – mesh")
Jutul.plot_mesh_edges!(ax, msh_homo, alpha = 0.2)
wells = get_model_wells(case_hom_inner.model)
for (name, well) in wells
    plot_well!(ax, msh_homo, well)
end
fig

# ### Simulate
results_hom_inner = run_case(case_hom_inner);
results_hom_outer = run_case(case_hom_outer);

# ### Well temperature profile
fig_hom = Figure(size = (600, 600))
colors = cgrad(:BrBG_4, 4, categorical = true)[[1, 2, 4]]

for (i, (c, r, l)) in enumerate(zip([results_hom_inner, results_hom_outer], [case_hom_inner, case_hom_outer], ["Injection in Inner pipe", "Injection in Outer pipe"]))
    ax = Axis(fig_hom[1, i];
        title = l,
        xlabel = "Temperature (°C)",
        ylabel = "Depth (m)",
        yreversed = true)

        well = r.model.models[:CoaxialWell_supply].data_domain
        T_well = convert_from_si.(
            c.result.states[end][:CoaxialWell_supply][:Temperature], :Celsius)
        tags = well[:tag] |> unique |> collect

    for (j, tag) in enumerate(tags)
        name = titlecase(replace(string(tag), "_" => " "))
        cells = well[:tag] .== tag
        Tn = T_well[cells]
        zn = well[:cell_centroids][3, cells]
        lines!(ax, Tn, zn; color = colors[j], linewidth = 4, label = name)
    end
    reservoir_cells = well.representation.perforations.reservoir
    T_reservoir = convert_from_si.(
        c.result.states[end][:Reservoir][:Temperature], :Celsius)
    T_reservoir_perf = T_reservoir[reservoir_cells]
    zn_reservoir = r.model.models[:Reservoir].data_domain[:cell_centroids][3, reservoir_cells]
    lines!(ax, T_reservoir_perf, zn_reservoir; color = :black, linewidth = 2, label = "Reservoir")
    axislegend(ax; position = :rt)
end
fig_hom

# ---
# ## Scenario 2 – Layered reservoir
# Four geological layers span the same 0–3000 m interval, each with distinct
# thermal and hydraulic properties. The layers mimic (from top to bottom):
# - **Clay** (0–250 m): moderate conductivity, low permeability
# - **Sandstone** (250–750 m): higher conductivity, higher permeability
# - **Granite** (750–1250 m): high conductivity, very low permeability
# - **Sandstone** (1250–3000 m): very low conductivity, low permeability (representing a shallow cap rock or insulating layer)
#
# Thermal conductivities are somewhat exaggerated to clearly show the impact
# of layering on the temperature profiles.
layered_args = (;
    base_args...,
    depths = [0.0, 250.0, 750.0, 1250.0, 3000.0],
    permeability  = [1e-3, 1e-1, 1e-4, 1e-1]*si_unit(:darcy),
    porosity      = [0.15, 0.20, 0.01, 0.20],
    rock_thermal_conductivity = [1.0, 2.0, 4.0, 2.0]*watt/(meter*Kelvin),
    rock_heat_capacity = [800, 900, 950, 900]*joule/(kilogram*Kelvin),
    rock_density       = [2650, 2600, 2580, 2600]*kilogram/meter^3,
);

case_layered = coaxial_well_branches(; inject_into = :outer, layered_args...);

# ### Plot reservoir properties
# The interactive viewer shows how conductivity and other properties vary with
# depth across the four layers.
plot_reservoir(case_layered.model)

# ### Simulate
results_layered = run_case(case_layered);

# ### Compare homogeneous vs layered temperature profiles
fig_layers = Figure(size = (1200, 600))

for (i, (results, case, label)) in enumerate(zip(
    [results_hom_outer, results_layered],
    [case_hom_outer, case_layered],
    ["Injection in Outer pipe", "Layered (styrofoam / clay / sandstone / granite)"]))

    ax = Axis(fig_layers[1, i];
        title = label,
        xlabel = "Temperature (°C)",
        ylabel = "Depth (m)",
        yreversed = true)

    well = case.model.models[:CoaxialWell_supply].data_domain
    T_well = convert_from_si.(
        results.result.states[end][:CoaxialWell_supply][:Temperature], :Celsius)
    tags = well[:tag] |> unique |> collect

    for (j, tag) in enumerate(tags)
        name = titlecase(replace(string(tag), "_" => " "))
        cells = well[:tag] .== tag
        Tn = T_well[cells]
        zn = well[:cell_centroids][3, cells]
        lines!(ax, Tn, zn; color = colors[j], linewidth = 4, label = name)
    end
    res_cells = well.representation.perforations.reservoir
    T_res = convert_from_si.(
        results.result.states[end][:Reservoir][:Temperature], :Celsius)
    T_res_perf = T_res[res_cells]
    zn_res = case.model.models[:Reservoir].data_domain[:cell_centroids][3, res_cells]
    scatter!(ax, T_res_perf, zn_res; color = :black, markersize = 4, label = "Reservoir")
    axislegend(ax; position = :rt)
end
fig_layers

# ---
# ## Scenario 3 – Effect of inner pipe wall thermal conductivity
# Using the homogeneous reservoir, we vary the inner pipe wall thermal
# conductivity from 0.0 (perfectly insulating) to 2× the default value
# (0.76 W/(m·K)). A higher conductivity increases heat transfer between the
# supply and return pipes, causing more thermal short-circuiting and lower
# production temperatures.
λ_default = 0.38  # default inner pipe wall thermal conductivity [W/(m·K)]
λ_values = [0.0*λ_default, λ_default, 4.0*λ_default]
labels_λ = ["λ = 0.00", "λ = 0.38 (default)", "λ = 1.52"]

cases_λ = []
results_λ = []
for λ in λ_values
    case_i = coaxial_well_branches(;
        inject_into = :outer,
        homogeneous_args...,
        well_args = (inner_pipe_thermal_conductivity = λ,
        outer_pipe_thermal_conductivity = λ_default,  # keep outer pipe conductivity constant),
    ))
    push!(cases_λ, case_i)
    push!(results_λ, run_case(case_i))
end

# ### Well temperature profiles for different inner pipe conductivities
fig_lambda = Figure(size = (1600, 600))
colors_λ = cgrad(:viridis, length(λ_values), categorical = true)

for (i, (results, case, label)) in enumerate(
    zip(results_λ, cases_λ, labels_λ))

    ax = Axis(fig_lambda[1, i];
        title = label,
        xlabel = "Temperature (°C)",
        ylabel = "Depth (m)",
        yreversed = true)

    well = case.model.models[:CoaxialWell_supply].data_domain
    T_well = convert_from_si.(
        results.result.states[end][:CoaxialWell_supply][:Temperature], :Celsius)
    tags = well[:tag] |> unique |> collect

    for (j, tag) in enumerate(tags)
        name = titlecase(replace(string(tag), "_" => " "))
        cells = well[:tag] .== tag
        Tn = T_well[cells]
        zn = well[:cell_centroids][3, cells]
        lines!(ax, Tn, zn; color = colors[j], linewidth = 4, label = name)
    end
    res_cells = well.representation.perforations.reservoir
    T_res = convert_from_si.(
        results.result.states[end][:Reservoir][:Temperature], :Celsius)
    T_res_perf = T_res[res_cells]
    zn_res = case.model.models[:Reservoir].data_domain[:cell_centroids][3, res_cells]
    scatter!(ax, T_res_perf, zn_res; color = :black, markersize = 4, label = "Reservoir")
    if i == length(λ_values)
        axislegend(ax; position = :rt)
    end
end
fig_lambda

# ### Production temperature comparison
# Overlay the outlet temperature over time for each inner pipe conductivity
# to show the effect on thermal short-circuiting.
fig_prod = Figure(size = (800, 500))
ax = Axis(fig_prod[1, 1];
    title = "Effect of inner pipe thermal conductivity",
    xlabel = "Time (years)",
    ylabel = "Production temperature (°C)")
for (i, (results, case, label)) in enumerate(zip(results_λ, cases_λ, labels_λ))
    T_prod = convert_from_si.(results.wells[:CoaxialWell_return][:temperature], :Celsius)
    t_years = convert_from_si.(cumsum(case.dt), :year)
    lines!(ax, t_years, T_prod; color = colors_λ[i], linewidth = 2, label = label)
end
axislegend(ax; position = :rb)
fig_prod
