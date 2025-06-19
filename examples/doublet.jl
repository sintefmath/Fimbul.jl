using Jutul, JutulDarcy, Fimbul
using HYPRE
using GLMakie

# ## Set up case
case, plot_args = Fimbul.geothermal_doublet()

# ## Inspect model
# We first plot the computational mesh and wells. The mesh is refined around the
# wells in the horizontal plane and vertically in and near the target aquifer.
msh = physical_representation(reservoir_model(case.model).data_domain)
fig = Figure(size = (1200, 800))
ax = Axis3(fig[1, 1], zreversed = true, aspect = plot_args.aspect)
Jutul.plot_mesh_edges!(ax, msh, alpha = 1.0)
wells = get_model_wells(case.model)
for (k, w) in wells
    plot_well!(ax, msh, w)
end
display(GLMakie.Screen(), fig)

# ### Plot reservoir properties
# Next, we Visualize the reservoir interactively.
plot_reservoir(case.model; plot_args...)

# ## Simulate geothermal energy production
# We simulate the geothermal doublet for 200 years
results = simulate_reservoir(case; info_level = 0)

# ## Visualize results
# We first plot the reservoir state interactively. You can notice how the
# cold front propagates from the injector well by filtering out high values.
plot_reservoir(case.model, results.states;
colormap = :seaborn_icefire_gradient, plot_args...)

# ##
plot_well_results(results.wells)

# ## Analyze well temperature evolution
geo = tpfv_geometry(msh)
well = case.model.models[:Producer]
cells = well.data_domain.representation.perforations.reservoir
z = geo.cell_centroids[3, cells]
states = results.result.states
times = convert_from_si.(cumsum(case.dt), :year)

fig = Figure(size = (1200, 800))
ax = Axis(fig[1, 1], xlabel = "Temperature (C)", ylabel = "Depth (m)", yreversed = true)
colors = cgrad(:seaborn_icefire_gradient, length(states), categorical = true)
for (n, state) in enumerate(states)
    T = state[:Producer][:Temperature][2:end]
    T = convert_from_si.(T, :Celsius)
    lines!(ax, T, z; linewidth = 4, color = colors[n], alpha = 0.25)
end

timesteps = Int64.(ceil.(range(1, length(states), length = 6)))
timesteps = [7, 21, 65, 200]
for n in timesteps
    T = states[n][:Producer][:Temperature][2:end]
    T = convert_from_si.(T, :Celsius)
    lines!(ax, T, z; linewidth = 4, color = colors[n], label = "$(times[n]) years")
end
axislegend(ax; position = :lc, fontsize = 20)
display(GLMakie.Screen(), fig)

##
Δstates = []
for state in results.states
    Δstate = Dict{Symbol, Any}()
    for (k, v) in state
        v0 = case.state0[:Reservoir][k]
        Δstate[k] = v .- v0
    end
    push!(Δstates, Δstate)
end
plot_reservoir(case.model, Δstates; colormap = :seaborn_icefire_gradient, plot_args...)

##
T_min = minimum(Δstates[end][:Temperature])
T_max = maximum(Δstates[1][:Temperature])
fig = Figure(size = (1500, 600))
for (i, n) in enumerate(timesteps)
    ax = Axis3(fig[1, i];
    title = "$(times[n]) years",
    zreversed = true,
    elevation = plot_args.elevation,
    aspect = :data)
    ΔT = Δstates[n][:Temperature]
    cells = geo.cell_centroids[1, :] .> -1000.0/2
    cells = cells .&& geo.cell_centroids[2, :] .> 50.0

    plot_cell_data!(ax, msh, ΔT; 
        cells = cells, colormap = :seaborn_icefire_gradient, colorrange = (T_min, T_max))
    plot_well!(ax, msh, wells[:Injector]; color = :black, linewidth = 4, top_factor = 0.4, markersize = 0.0)
    plot_well!(ax, msh, wells[:Producer]; color = :black, linewidth = 4)
    hidedecorations!(ax)
end
Colorbar(fig[2, 1:length(timesteps)]; 
    colormap = :seaborn_icefire_gradient, colorrange = (T_min, T_max), 
    label = "ΔT (C)", vertical = false)
display(GLMakie.Screen(), fig)