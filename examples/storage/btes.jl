# # Borehole Thermal Energy Storage (BTES)
using Jutul, JutulDarcy
using Fimbul
using HYPRE

# ## Set up simulation case
# We consider a domain with 50 BTES wells that all reach 100 m depth. The BTES
# wells are charged during the summer months and discharged during the winter
# months. The simulation is run for 10 years.
case = btes(num_wells = 50, depths = [0.0, 0.5, 100, 125],
    charge_months = ["April", "May", "June", "July", "August", "September"],
    discharge_months = ["October", "November", "December", "January", "February", "March"],
    num_years = 10,
);

# ## Set up reservoir simulator
simulator, config = setup_reservoir_simulator(case;
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
sel = JutulDarcy.ControlChangeTimestepSelector(
    case.model, 0.0, convert_to_si(5.0, :second))
push!(config[:timestep_selectors], sel)
config[:timestep_max_decrease] = 1e-6
for ws in well_symbols(case.model)
    config[:tolerances][ws][:default] = 1e-2
end

# ## Simulate the case
results = simulate_reservoir(case[1:19], simulator=simulator, config=config, info_level=2, restart = true);

# ## Visualize the results

# ### Interactive visualization of the reservoir state
plot_reservoir(case.model, results.states)

# ### Interactive plot of the well output
plot_well_results(results.wells)


##

number_of_section = 8
rate_charge = 25.0e-3
temperature_charge = 90.0 + 273.15 # K
temperature_discharge = 10.0 + 273.15 # K
rate_discharge = rate_charge
rho = reservoir_model(model).system.rho_ref[1]

model = case.model
rate_target = TotalRateTarget(rate_charge)
ctrl_charge = InjectorControl(rate_target, [1.0], 
    density=rho, temperature=temperature_charge)
rate_target = TotalRateTarget(rate_discharge)
ctrl_discharge = InjectorControl(rate_target, [1.0],
    density=rho, temperature=temperature_discharge);
# BHP control for return side
bhp_target = BottomHolePressureTarget(5.0si_unit(:atm))
ctrl_prod = ProducerControl(bhp_target);
# Set up forces
control_charge = Dict()
control_discharge = Dict()

number_of_section = 8

model = case.model
msh = physical_representation(reservoir_model(model).data_domain)
geo = tpfv_geometry(msh)
xy = geo.cell_centroids[1:2,:]
# map from (x,y) to polar coordiates
r = sqrt.(xy[1,:].^2 .+ xy[2,:].^2)
θ = atan.(xy[2,:], xy[1,:]) .+ π
assigned = []
wells = well_symbols(model)
filter!(well -> contains(String(well), "_supply"), wells)
get_return = (well) -> Symbol(replace(String(well), "_supply" => "_return"))
for sno in 1:number_of_section
    θ_min = (sno - 1)*2π/number_of_section
    θ_max = sno*2π/number_of_section
    in_section = θ_min .<= θ .<= θ_max
    section_wells, well_radii = Symbol[], Float64[]
    for well in wells
        well in assigned ? continue : nothing
        wmodel = model.models[well]
        wc = wmodel.domain.representation.perforations.reservoir[1]
        !in_section[wc] ? continue : nothing
        push!(section_wells, well)
        push!(well_radii, r[wc])
    end
    order = sortperm(well_radii)
    println("Section $sno: ", section_wells)
    for wno in order
        well_sup = section_wells[wno]
        well_ret = get_return(well_sup)
        println("Assigned wells: $assigned")
        @assert well_sup ∉ assigned
        @assert well_ret ∉ assigned
        if wno == 1 
            control_charge[well_sup] = ctrl_charge
            control_discharge[well_sup] = ctrl_discharge
        else
            well_prev = get_return(section_wells[wno-1])
            target = JutulDarcy.ReinjectionTarget(NaN, [well_prev])
            ctrl = InjectorControl(target, [1.0],
                density=rho, temperature=NaN; check=false)
            control_charge[well_sup] = ctrl
            control_discharge[well_sup] = ctrl
        end
        control_charge[well_ret] = ctrl_prod
        control_discharge[well_ret] = ctrl_prod
        push!(assigned, well_sup, well_ret)
    end
end

@assert sort(assigned) == sort(well_symbols(model))