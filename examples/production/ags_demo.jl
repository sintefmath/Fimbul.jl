using Jutul, JutulDarcy
using Fimbul
using HYPRE
using GLMakie

##

case = Fimbul.ags();

##
msh = physical_representation(reservoir_model(case.model).data_domain)
fig = Figure(size = (800, 600))
ax = Axis3(fig[1, 1], title = "Extruded mesh for AGS demo case",
    xlabel = "X (m)", ylabel = "Y (m)", zlabel = "Z (m)", aspect = :data)
plot_mesh_edges!(ax, msh)

for (wn, well) in get_model_wells(case.model)
    plot_well!(ax, msh, well)
end

fig

##
results = simulate_reservoir(case; info_level = 2)

##


well_coords0 = Fimbul.get_ags_trajectory()

Δ = 10.0
well_coords = []
for wc in well_coords0
    wc_left = copy(wc)
    wc_left[:, 2] .-= Δ/2
    push!(well_coords, wc_left)
    wc_right = copy(wc)
    wc_right[:, 2] .+= Δ/2
    push!(well_coords, wc_right)
end

cell_constraints = []
for wc in well_coords
    xw = [Tuple(x) for x in eachrow(unique(wc[:, 1:2], dims=1))]
    for cc in cell_constraints
        cor = [norm(collect(x) .- collect(y), 2) for x in xw, y in cc]
        println("Minimum distance between wells: ", cor)
        xw = [xw[i] for i in eachindex(xw) if all(cor[i, :] .>= Δ)]
    end
    push!(cell_constraints, xw)
end

##
hxy_min = 20.0
msh = Fimbul.extruded_mesh(
    cell_constraints,
    [0.0, 1500.0, 2300, 2400, 2800],
    hz = [250.0, 100.0, 10.0, 100.0],
    hxy_min = hxy_min,
    offset_rel = 1.0,
    dist_min_factor = 50.0,
);

##



##
xc = vcat(cell_constraints...)

