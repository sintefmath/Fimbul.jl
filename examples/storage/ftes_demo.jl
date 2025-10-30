# # Fractured Thermal Energy Storage (FTES)
# This example demonstrates the simulation of a fractured thermal energy storage
# system in Fimbul. The simulation includes setting up the reservoir model,
# configuring the simulator, running the simulation, and visualizing the
# results.

# Add required modules to namespace 
using Jutul, JutulDarcy, Fimbul # Core reservoir simulation framework
using HYPRE # High-performance linear solvers
using GLMakie # 3D visualization and plotting capabilities

# ## Set up simulation case
# We use the ftes setup function to set up a simulation case. The function is
# parametrized, and supports
# * `matrix_permeability`: Permeability of the reservoir matrix [m^2]
# * `fracture_permeability`: Permeability of the fractures [m^2]
# * `depths`: Depths of the reservoir layers [m]
# * `num_fractures`: Number of fractures per layer [-]

case = Fimbul.ftes()

##
plot_reservoir(case.model; aspect = :data, key = :permeability, colormap = :viridis)

##
sim, cfg = setup_reservoir_simulator(case;
    info_level = 2, # 0=progress bar, 1=basic, 2=detailed
    output_substates = true,
    tol_cnv = Inf, # Disable CNV tolerance - use other criteria
    inc_tol_dT = 1e-2, # Temperature increment tolerance [K]
    inc_tol_dp_rel = 1e-3, # Relative pressure increment tolerance [-]
    initial_dt = 5.0, # Initial timestep [s]
    relaxation = true); # Enable relaxation in Newton solver

sel = JutulDarcy.ControlChangeTimestepSelector(case.model, 0.1, 5.0)
push!(cfg[:timestep_selectors], sel)
sel = VariableChangeTimestepSelector(:Temperature, 5.0; model = :Reservoir, relative = false)
push!(cfg[:timestep_selectors], sel)
cfg[:timestep_max_decrease] = 1e-6; # Prevent excessive timestep reduction
results = simulate_reservoir(case; simulator = sim, config = cfg)

##
plot_reservoir(case.model, results.states;
aspect = :data, key = :Temperature, colormap = :seaborn_icefire_gradient)

##
plot_well_results(results.wells)