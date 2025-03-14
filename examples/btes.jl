using Jutul, JutulDarcy
using Fimbul
using HYPRE

##
case = btes(num_wells = 50, depth = 100, meshing_args = (hz = 10.0,), 
    temperature_charge = convert_to_si(50.0, :Celsius),
    charge_months = ["April", "May", "June", "July", "August", "September"],
    discharge_months = ["October", "November", "December", "January", "February", "March"],
    rate_charge = 0.5si_unit(:litre)/si_unit(:second),
    rate_discharge = 0.5si_unit(:litre)/si_unit(:second),
    num_years = 20,
    report_interval = Inf,
)

simulator, config = setup_reservoir_simulator(case;
    tol_cnv = 1e-2,
    tol_mb = 1e-6,
    timesteps = :auto,
    initial_dt = 5.0,
    target_its = 5,
    max_nonlinear_iterations = 8,
    relaxation = true
);
thresholds = Dict(:B1_supply => 0.0)
sel = Fimbul.ControlChangeTimestepSelector(thresholds, convert_to_si(5.0, :second))
push!(config[:timestep_selectors], sel)
config[:timestep_max_decrease] = 1e-6
for ws in well_symbols(case.model)
    config[:tolerances][ws][:default] = 1e-2
end

##
results = simulate_reservoir(case, simulator=simulator, config=config, info_level=2);

##
plot_reservoir(case.model, results.states)

##
plot_well_results(results.wells)

##
c = case.model.models[:B1_return].domain.representation.perforations.reservoir
mesh = physical_representation(reservoir_model(case.model).data_domain);
geo = tpfv_geometry(mesh)
z = geo.cell_centroids[3,c]
println(z)