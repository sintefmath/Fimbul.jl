using Fimbul
using Jutul, JutulDarcy
using HYPRE

# ## Set up case
case = Fimbul.geothermal_doublet(num_years=50)

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