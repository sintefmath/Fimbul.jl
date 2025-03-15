## # Borehole Thermal Energy Storage (BTES)
# NOTE: This example currently relies on a timestepper that is not yet part of
# the main branch of Jutul. To run the example, check out the ´timestepping´
# branch on Jutul. The example will be updated once the branch is merged.
using Jutul, JutulDarcy
using Fimbul
using HYPRE

## ## Set up simulation case
# We consider a domain with 50 BTES wells that all reach 100 m depth. The BTES
# wells are charged during the summer months and discharged during the winter
# months. The simulation is run for 10 years.
case = btes(num_wells = 50, depth = 100;
    charge_months = ["April", "May", "June", "July", "August", "September"],
    discharge_months = ["October", "November", "December", "January", "February", "March"],
    num_years = 10,
)

# ## Set up reservoir simulator
Simulator, config = setup_reservoir_simulator(case;
    tol_cnv = 1e-2,
    tol_mb = 1e-6,
    timesteps = :auto,
    initial_dt = 5.0,
    target_its = 5,
    max_nonlinear_iterations = 8,
    relaxation = true
);
# Changing from charging to discharging results in a thermal shock that is
# challenging to resolve for the nonlinear solver. We therefore use a timestep
# selector that reduces the timestep to 5 seconds when the control changes.
thresholds = Dict(:B1_supply => 0.0)
sel = setup_control_change_timestep_selector(
    case.model, 0.0, convert_to_si(5.0, :second))
push!(config[:timestep_selectors], sel)
config[:timestep_max_decrease] = 1e-6
for ws in well_symbols(case.model)
    config[:tolerances][ws][:default] = 1e-2
end

## ## Simulate the case
results = simulate_reservoir(case, simulator=simulator, config=config, info_level=2);

## ## Visualize the results

## ### Interactive visualization of the reservoir state
plot_reservoir(case.model, results.states)

## ### Interactive plot of the well output
plot_well_results(results.wells)
