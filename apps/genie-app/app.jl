module App
using GenieFramework
using JutulDarcy
using Jutul
using Fimbul
using DataFrames
using CairoMakie
using Dates
using FileIO
using Base64

@genietools

@app begin
    @in temperature_charge = 95.0
    @in temperature_discharge = 25.0
    @in rate_charge_l_s = 25.0
    @in permeability_md = 1000.0
    @in porosity = 0.35
    @in well_distance = 500.0
    @in aquifer_thickness = 100.0
    @in depth = 1000.0
    @in rock_thermal_conductivity = 2.0
    @in rock_heat_capacity = 900.0
    
    # Cap rock / Basement properties
    @in permeability_cap_md = 1.0
    @in porosity_cap = 0.01
    @in rock_thermal_conductivity_cap = 2.0
    @in rock_heat_capacity_cap = 900.0
    
    @in thermal_gradient = 30.0
    @in temperature_surface = 10.0
    @in selected_tab = "parameters"

    @in run_simulation = false
    
    @out simulation_status = "Ready"
    @out plot_data = PlotData()
    @out reservoir_img = ""
    @in plot_step = 1
    @out max_step = 1
    
    @onchange run_simulation begin
        if run_simulation
            simulation_status = "Running..."
            try
                # Run simulation
                results, case = run_ates_simulation(
                    temperature_charge, 
                    temperature_discharge, 
                    rate_charge_l_s, 
                    permeability_md, 
                    porosity,
                    well_distance,
                    aquifer_thickness,
                    depth,
                    rock_thermal_conductivity,
                    rock_heat_capacity,
                    permeability_cap_md,
                    porosity_cap,
                    rock_thermal_conductivity_cap,
                    rock_heat_capacity_cap,
                    thermal_gradient,
                    temperature_surface
                )
                simulation_status = "Completed"
                
                # Store results globally (simple hack for single-user local app)
                global current_results = results
                global current_case = case
                
                # Calculate global min/max for plotting
                T_all = reduce(vcat, [s[:Temperature] for s in results.states]) .- 273.15
                global current_limits = (minimum(T_all), maximum(T_all))

                max_step = length(results.states)
                plot_step = 1
                reservoir_img = generate_plot_image(plot_step)
                
            catch e
                simulation_status = "Error: $e"
                @error e
            end
            run_simulation = false
        end
    end
    
    @onchange plot_step begin
        reservoir_img = generate_plot_image(plot_step)
    end
end

# Global storage for results
current_results = nothing
current_case = nothing
current_limits = nothing

function generate_plot_image(step_index)
    global current_results, current_case, current_limits
    if isnothing(current_results)
        return ""
    end
    
    # Generate reservoir image using CairoMakie
    # We will plot the temperature field
    
    # Extract data
    state = current_results.states[step_index]
    T = state[:Temperature] .- 273.15 # Convert to Celsius
    
    # Get mesh
    model = current_case.model
    msh = physical_representation(reservoir_model(model).data_domain)
    
    # Create plot
    fig = Figure(size = (800, 600))
    ax = Axis3(fig[1, 1], title = "Temperature at step $(step_index) (°C)", aspect = :data, zreversed=true)
    
    # Use Jutul's plotting if possible, but we need to capture it.
    # Jutul.plot_mesh_edges!(ax, msh, alpha = 0.2)
    # We can use plot_cell_data! from Jutul if we import it or use Fimbul's re-exports.
    # Fimbul exports plot_well_data! but maybe not plot_cell_data! directly?
    # JutulDarcy exports plot_cell_data!
    
    # Let's try to use JutulDarcy.plot_cell_data!
    JutulDarcy.plot_cell_data!(ax, msh, T, colormap = :seaborn_icefire_gradient, colorrange = current_limits)
    
    # Add wells
    # wells = get_model_wells(model)
    # colors = [:red, :blue]
    # for (i, (k, w)) in enumerate(wells)
    #     display(w)
    #     JutulDarcy.plot_well!(ax, msh, w;
    #     # color = colors[i], linewidth = 6
    #     )
    # end
    
    # Save to buffer
    io = IOBuffer()
    show(io, MIME("image/png"), fig)
    img_data = take!(io)
    base64_img = base64encode(img_data)
    return "data:image/png;base64,$base64_img"
end

function run_ates_simulation(t_charge, t_discharge, rate_l_s, perm_md, phi, well_dist, aq_thick, d, rock_cond, rock_cp, perm_cap_md, phi_cap, rock_cond_cap, rock_cp_cap, grad, t_surf)
    darcy, litre, second, watt, meter, Kelvin, joule, kilogram = si_units(:darcy, :litre, :second, :watt, :meter, :Kelvin, :joule, :kilogram)
    
    # Simple properties: [aquifer, cap/basement]
    # Cap rock permeability
    perms = [perm_md * 1e-3 * darcy, perm_cap_md * 1e-3 * darcy]
    
    # Cap rock porosity
    porosities = [phi, phi_cap]

    # Thermal properties
    rock_conds = [rock_cond, rock_cond_cap] .* watt / (meter * Kelvin)
    rock_cps = [rock_cp, rock_cp_cap] .* joule / (kilogram * Kelvin)
    
    # Coarser mesh for speed
    # hxy_min default is nearwell_radius/6. nearwell_radius ~ 125m/2 = 62.5m. hxy_min ~ 10m.
    # Let's set hxy_min to 20m.
    
    case = Fimbul.ates_simple(
        temperature_charge = convert_to_si(t_charge, :Celsius),
        temperature_discharge = convert_to_si(t_discharge, :Celsius),
        rate_charge = rate_l_s * litre / second,
        permeability = perms,
        porosity = porosities,
        well_distance = well_dist,
        aquifer_thickness = aq_thick,
        depth = d,
        rock_thermal_conductivity = rock_conds,
        rock_heat_capacity = rock_cps,
        thermal_gradient = (grad / 1000.0) * Kelvin / meter,
        temperature_surface = convert_to_si(t_surf, :Celsius),
        use_2d = true,
        mesh_args = (hxy_min = 20.0, hxy_max = 100.0, hz_min = 10.0, hz_max = 50.0),
        utes_schedule_args = (num_years = 1,)
    )
    
    sim, cfg = setup_reservoir_simulator(case; info_level = 0)
    sel = JutulDarcy.ControlChangeTimestepSelector(case.model)
    push!(cfg[:timestep_selectors], sel)
    cfg[:timestep_max_decrease] = 1e-3
    
    results = simulate_reservoir(case, simulator = sim, config = cfg)
    
    return results, case
end

ui() = [
    h1("Fimbul ATES Simulator"),
    p("Simulate Aquifer Thermal Energy Storage (ATES) systems using Fimbul and JutulDarcy."),
    tabgroup(:selected_tab, class="bg-primary text-white shadow-2", [
        tab(name="parameters", label="Parameters"),
        tab(name="results", label="Results")
    ]),
    tabpanels(:selected_tab, animated=true, [
            tabpanel(name="parameters", [
                # h4("Parameters"),
                
                # h5("Geology"),

                row([
                    cell(p("Aquifer Thickness (m)"), size=3),
                    cell(slider(1.0:1.0:500.0, :aquifer_thickness; label=true), size=6),
                    cell(textfield("", :aquifer_thickness, type="number", inputclass="text-right"), size=3, class="q-pl-sm")
                ], class="items-center"),

                row([
                    cell(p("Aquifer Depth (m)"), size=3),
                    cell(slider(10.0:1.0:3000.0, :depth; label=true), size=6),
                    cell(textfield("", :depth, type="number", inputclass="text-right"), size=3, class="q-pl-sm")
                ], class="items-center"),

                row([
                    cell(p("Thermal Gradient (°C/km)"), size=3),
                    cell(slider(0.0:1.0:100.0, :thermal_gradient; label=true), size=6),
                    cell(textfield("", :thermal_gradient, type="number", inputclass="text-right"), size=3, class="q-pl-sm")
                ], class="items-center"),

                row([
                    cell(p("Surface Temperature (°C)"), size=3),
                    cell(slider(1.0:1.0:30.0, :temperature_surface; label=true), size=6),
                    cell(textfield("", :temperature_surface, type="number", inputclass="text-right"), size=3, class="q-pl-sm")
                ], class="items-center"),

                row([
                    cell(size=2),
                    cell(h6("Aquifer"), size=5),
                    cell(h6("Cap/Basement"), size=5)
                ]),

                row([
                    cell(p("Permeability (mD)"), size=2),
                    cell(slider(1.0:1.0:2000.0, :permeability_md; label=true), size=3),
                    cell(textfield("", :permeability_md, type="number", inputclass="text-right"), size=2, class="q-pl-sm"),
                    cell(slider(0.001:0.001:2000.0, :permeability_cap_md; label=true), size=3, class="q-pl-md"),
                    cell(textfield("", :permeability_cap_md, type="number", inputclass="text-right"), size=2, class="q-pl-sm")
                ], class="items-center"),

                row([
                    cell(p("Porosity"), size=2),
                    cell(slider(0.01:0.01:1.0, :porosity; label=true), size=3),
                    cell(textfield("", :porosity, type="number", inputclass="text-right"), size=2, class="q-pl-sm"),
                    cell(slider(0.001:0.001:1.0, :porosity_cap; label=true), size=3, class="q-pl-md"),
                    cell(textfield("", :porosity_cap, type="number", inputclass="text-right"), size=2, class="q-pl-sm")
                ], class="items-center"),

                row([
                    cell(p("Rock Thermal Conductivity (W/m/K)"), size=2),
                    cell(slider(0.5:0.1:10.0, :rock_thermal_conductivity; label=true), size=3),
                    cell(textfield("", :rock_thermal_conductivity, type="number", inputclass="text-right"), size=2, class="q-pl-sm"),
                    cell(slider(0.5:0.1:10.0, :rock_thermal_conductivity_cap; label=true), size=3, class="q-pl-md"),
                    cell(textfield("", :rock_thermal_conductivity_cap, type="number", inputclass="text-right"), size=2, class="q-pl-sm")
                ], class="items-center"),

                row([
                    cell(p("Rock Heat Capacity (J/kg/K)"), size=2),
                    cell(slider(500.0:50.0:2000.0, :rock_heat_capacity; label=true), size=3),
                    cell(textfield("", :rock_heat_capacity, type="number", inputclass="text-right"), size=2, class="q-pl-sm"),
                    cell(slider(500.0:50.0:2000.0, :rock_heat_capacity_cap; label=true), size=3, class="q-pl-md"),
                    cell(textfield("", :rock_heat_capacity_cap, type="number", inputclass="text-right"), size=2, class="q-pl-sm")
                ], class="items-center"),

                h5("Operational Parameters"),
                row([
                    cell(p("Well Distance (m)"), size=3),
                    cell(slider(100.0:50.0:2000.0, :well_distance; label=true), size=6),
                    cell(textfield("", :well_distance, type="number", inputclass="text-right"), size=3, class="q-pl-sm")
                ], class="items-center"),

                row([
                    cell(p("Charge Temp (°C)"), size=3),
                    cell(slider(1:1:150, :temperature_charge; label=true), size=6),
                    cell(textfield("", :temperature_charge, type="number", inputclass="text-right"), size=3, class="q-pl-sm")
                ], class="items-center"),

                row([
                    cell(p("Discharge Temp (°C)"), size=3),
                    cell(slider(1:1:150, :temperature_discharge; label=true), size=6),
                    cell(textfield("", :temperature_discharge, type="number", inputclass="text-right"), size=3, class="q-pl-sm")
                ], class="items-center"),

                row([
                    cell(p("Charge Rate (L/s)"), size=3),
                    cell(slider(0.1:0.1:100, :rate_charge_l_s; label=true), size=6),
                    cell(textfield("", :rate_charge_l_s, type="number", inputclass="text-right"), size=3, class="q-pl-sm")
                ], class="items-center"),

                btn("Run Simulation", @click("run_simulation = true"), loading=:run_simulation)
            ]),
            tabpanel(name="results", [
                h4("Status: {{simulation_status}}"),
                h5("Reservoir Temperature"),
                imageview(src = :reservoir_img, width="100%"),
                slider(1:1:10, :plot_step; min=1, max=:max_step, step=1, label=true)
            ])
        ])
]

@page("/", ui)

end
