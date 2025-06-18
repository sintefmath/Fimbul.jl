using Jutul, JutulDarcy, Fimbul
using HYPRE
using GLMakie

# ## Set up case
case = Fimbul.geothermal_doublet()

# ## Inspect model
# We first plot the computational mesh and wells. The mesh is refined around the
# wells in the horizontal plane and vertically in and near the target aquifer.
msh = physical_representation(reservoir_model(case.model).data_domain)
fig = Figure(size = (1200, 800))
ax = Axis3(fig[1, 1], zreversed = true)
Jutul.plot_mesh_edges!(ax, msh, alpha = 1.0)
wells = get_model_wells(case.model)
for (k, w) in wells
    plot_well!(ax, msh, w)
end
display(GLMakie.Screen(), fig)

plot_reservoir(case.model)

# ## Simulate
results = simulate_reservoir(case; info_level = 2)

# ## Visualize results
plot_reservoir(case.model, results.states)

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
timesteps = [1, 10, 50, 100, 200]
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
plot_reservoir(case.model, Δstates)

for n in timesteps
    fig = Figure(size = (1200, 800))
    ax = Axis3(fig[1, 1], zreversed = true)
    ΔT = Δstates[n][:Temperature]
    cells = ΔT .< -15.0 .&& geo.cell_centroids[2, :] .< 0.0
    plot_cell_data!(ax, msh, ΔT; cells = cells, colormap = :seaborn_icefire_gradient)
    plot_well!(ax, msh, wells[:Injector]; color = :black, linewidth = 4)
    plot_well!(ax, msh, wells[:Producer]; color = :black, linewidth = 4)
    display(GLMakie.Screen(), fig)
end

##
msh = physical_representation(reservoir_model(case.model).data_domain)
fig = Figure(size = (1200, 800))
ax = Axis3(fig[1, 1], zreversed = true)
Jutul.plot_mesh_edges!(ax, msh, alpha = 1.0)
wells = get_model_wells(case.model)
for (k, w) in wells
    plot_well!(ax, msh, w)
end
display(GLMakie.Screen(), fig)

##

cc = map(bc -> bc.cell, case.forces[:Reservoir].bc)