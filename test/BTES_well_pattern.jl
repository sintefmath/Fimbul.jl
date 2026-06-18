using GLMakie
using Jutul
using Fimbul
using JutulDarcy

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