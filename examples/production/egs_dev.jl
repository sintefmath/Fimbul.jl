# Add required modules to namespace 
using Jutul, JutulDarcy, Fimbul
using HYPRE
using Statistics
using GLMakie
# Useful SI units

##
Kelvin, joule, watt = si_units(:Kelvin, :joule, :watt)
kilogram = si_unit(:kilogram)
meter = si_unit(:meter)
darcy = si_unit(:darcy);

well_spacing = 200.0 .* meter
fracture_radius = 300.0 .* meter
well_depth = 2500.0 .* meter
well_lateral = 2000.0 .* meter
fracture_aperture = 1e-3 .* meter
hx_min = well_spacing/5
ws, wd, wl = well_spacing, well_depth, well_lateral

well_coords = [
    [-ws/2 0.0  0.0; -ws/2 0.0 wd; -ws/2 wl wd], 
    [ ws/2 0.0  0.0;  ws/2 0.0 wd;  ws/2 wl wd]
]

day, hour = si_units(:day, :hour)
case = Fimbul.egs(well_coords, fracture_radius, 200.0;
    num_years = 1,
    mesh_args = (hx_min = hx_min, ),
    schedule_args = (report_interval = 1day,)
);

##
msh = physical_representation(reservoir_model(case.model).data_domain)
plot_mesh_edges(msh)

##


##
output_path = jutul_output_path("egs_demo")
sim, cfg = setup_reservoir_simulator(case;
    method = :nldd,
    nldd_partition = p,
    solve_tol_temperature = 0.5,
    info_level = 2,
    output_path = output_path,
    tol_cnv = 1e-2,
    tol_mb = 1e-6,
    timesteps = :auto,
    max_timestep = 3.0si_unit(:hour),
    initial_dt = 5.0,
    relaxation = true);
cfg[:tolerances][:Reservoir][:default] = 1e-2
cfg[:linear_solver].config.max_iterations = 100
# The transition from charging to discharging creates a thermal shock that is
# numerically challenging for the nonlinear solver. We use a specialized
# timestep selector that reduces the timestep to 5 seconds during control
# changes to maintain numerical stability and convergence.
sel = JutulDarcy.ControlChangeTimestepSelector(
    case.model, 0.0, convert_to_si(5.0, :second))
push!(cfg[:timestep_selectors], sel)
cfg[:timestep_max_decrease] = 1e-6
for ws in well_symbols(case.model)
    cfg[:tolerances][ws][:default] = 1e-2
end

##
res = simulate_reservoir(case, simulator = sim, config = cfg, restart = true)

##
plot_reservoir(case, res.states; aspect = :data)

##
r_domain = reservoir_model(case.model).data_domain
msh = physical_representation(r_domain)

block_size = 1000
num_blocks = Int.(ceil(number_of_cells(msh)./block_size))
p = JutulDarcy.partition_reservoir(case.model, num_blocks; wells_in_single_block = true)

is_frac = r_domain[:porosity] .== 0.5
is_well = falses(number_of_cells(msh))

for (name, well) in get_model_wells(case.model)
   cells = well.perforations.reservoir
   is_well[cells] .= true
end

# Pad domain
using SparseArrays
ii = [msh.faces.neighbors[f][1] for f in 1:number_of_faces(msh)]
jj = [msh.faces.neighbors[f][2] for f in 1:number_of_faces(msh)]
M = sparse(vcat(ii, jj), vcat(jj, ii), 1.0)
for k = 1:1
   global is_well = is_well .|| M*is_well .> 0
   global is_frac = is_frac .|| M*is_frac .> 0
end
pmax = maximum(p)
p[is_frac] .= pmax + 1
p[is_well] .= pmax + 1

p = Jutul.compress_partition(p)


