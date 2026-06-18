using GLMakie
using Jutul
using Fimbul
using JutulDarcy

println(pathof(Fimbul))

field = [
    [ [8.0 8.0; 4.0 4.0; 65.0 65.0], [10.0 10.0; 20.0 20.0; 80.0 80.0] ],   # sector 1: 2 wells
    [ [0.0 5.0; 0.0 15.0; 0.0 30.0], [5.0 5.0; 5.0 5.0; 65.0 65.0] ]    # sector 2: 2 wells
] 

pattern = :rectangular

case = btes(pattern = pattern,num_wells = 49, depths = [0.0, 0.5, 100, 125],
    charge_period = ["April", "September"],
    discharge_period = ["October", "March"],
    num_years = 4,
);

# Nominal grid positions, straight from the pattern function (no mesh involved)
nominal_field = Fimbul.rectangular_pattern(48, 5.0; num_sectors = 6, depths = [0.0, 0.5, 100, 125])
nominal_xy = hcat([w[1:2, 1] for w in vcat(nominal_field...)]...)
println("First 5 nominal (x,y): ", [nominal_xy[:, i] for i in 1:5])

# Actual positions used by the built wells (cell centroids the well ended up at)
geo = tpfv_geometry(physical_representation(reservoir_model(case.model).data_domain))
wells = filter(w -> contains(String(w), "_supply"), well_symbols(case.model))
actual_xy = [geo.cell_centroids[1:2, case.model.models[w].domain.representation.perforations.reservoir[1]] for w in wells]
println("First 5 actual (x,y): ", actual_xy[1:5])

# ------- Important: well pattern

msh = physical_representation(reservoir_model(case.model).data_domain)
fig = Figure(size = (800, 800))
ax = Axis3(fig[1, 1]; zreversed = true, aspect = :data,
azimuth = 0, elevation = π/2)
Jutul.plot_mesh_edges!(ax, msh, alpha = 0.2)
colors = Makie.wong_colors()
lns, labels = [], String[]
sector_keys = sort(collect(keys(case.input_data[:sectors])),
    by = k -> parse(Int, replace(String(k), "S" => "")))
for (sno, sector_key) in enumerate(sector_keys)
    wells = case.input_data[:sectors][sector_key]
    for (wno, wname) in enumerate(wells)
        well = case.model.models[wname].domain.representation
        l = plot_well!(ax, msh, well; color = colors[sno], fontsize = 0)
        if wno == 1
            push!(lns, l)
            push!(labels, "Sector $sno")
        end
    end
end
Legend(fig[1, 2], lns, labels);
fig