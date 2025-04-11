using Fimbul
using Jutul, JutulDarcy
using HYPRE

##
case = Fimbul.geothermal_doublet(num_years=50)

##
results = simulate_reservoir(case)

##
plot_reservoir(case.model, results.states)

##
plot_well_results(results.wells)

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