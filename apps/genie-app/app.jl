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
    @in permeability_darcy = 1.0
    @in porosity = 0.35
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
                    permeability_darcy, 
                    porosity
                )
                simulation_status = "Completed"
                
                # Store results globally (simple hack for single-user local app)
                global current_results = results
                global current_case = case
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

function generate_plot_image(step_index)
    global current_results, current_case
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
    JutulDarcy.plot_cell_data!(ax, msh, T, colormap = :seaborn_icefire_gradient)
    
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

function run_ates_simulation(t_charge, t_discharge, rate_l_s, perm_d, phi)
    darcy, litre, second = si_units(:darcy, :litre, :second)
    
    # Default properties
    perms = [1.0, 5.0, 1000.0, 5.0, 1.0] .* 1e-3 .* darcy
    # Update aquifer layer (index 3)
    perms[3] = perm_d * darcy
    
    porosities = [0.01, 0.05, 0.35, 0.05, 0.01]
    porosities[3] = phi
    
    # Coarser mesh for speed
    # hxy_min default is nearwell_radius/6. nearwell_radius ~ 125m/2 = 62.5m. hxy_min ~ 10m.
    # Let's set hxy_min to 20m.
    
    case = ates(
        temperature_charge = convert_to_si(t_charge, :Celsius),
        temperature_discharge = convert_to_si(t_discharge, :Celsius),
        rate_charge = rate_l_s * litre / second,
        permeability = perms,
        porosity = porosities,
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
    row([
        cell([
            h4("Parameters"),
            textfield("Charge Temp (°C)", :temperature_charge),
            textfield("Discharge Temp (°C)", :temperature_discharge),
            textfield("Charge Rate (L/s)", :rate_charge_l_s),
            textfield("Permeability (Darcy)", :permeability_darcy),
            textfield("Porosity", :porosity),
            btn("Run Simulation", @click("run_simulation = true"), loading=:run_simulation)
        ]),
        cell([
            h4("Status: {{simulation_status}}"),
            h5("Reservoir Temperature"),
            imageview(src = :reservoir_img, width="100%"),
            slider(1:1:10, :plot_step; min=1, max=:max_step, step=1, label=true)
        ])
    ])
]

@page("/", ui)

end
