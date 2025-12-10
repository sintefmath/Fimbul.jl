# # Analytical solutions for Borehole Thermal Energy Storage (BTES)
# This example validates the numerical simulation of a BTES system against
# analytical solutions for both U1 and coaxial closed-loop configurations
# assuming fixed ground temperature.

# Add modules to namespace
using Jutul, JutulDarcy, Fimbul
using GLMakie

# Useful SI units
meter = si_unit(:meter)
kilogram = si_unit(:kilogram)
atm = si_unit(:atm)
second, day = si_units(:second, :day)
Kelvin, Joule, Watt = si_units(:kelvin, :joule, :watt)
darcy = si_unit(:darcy);

# ## Problem parameters
# We set up common parameters for both BTES configurations.

# Operational conditions
T_in = convert_to_si(80.0, :Celsius) # Inlet temperature
T_rock = convert_to_si(10.0, :Celsius) # Ground temperature (assumed constant)
Q = 20.0*meter^3/day # Volumetric flow rate
# Fluid properties
ρf = 988.1*kilogram/meter^3 # Fluid density
Cpf = 4184.0*Joule/(kilogram*Kelvin) # Fluid heat capacity
λf = 0.6405*Watt/(meter*Kelvin) # Fluid thermal conductivity
# Rock properties
ϕ = 0.01 # Porosity
K = 1e-3*darcy # Permeability
ρr = 2650.0*kilogram/meter^3 # Rock density
Cpr = 900.0*Joule/(kilogram*Kelvin) # Rock heat capacity
λr = 2.5*Watt/(meter*Kelvin) # Rock thermal conductivity
# BTES material properties
ρg = 2000.0*kilogram/meter^3 # Grout density
Cpg = 1500.0*Joule/(kilogram*Kelvin) # Grout heat capacity
λg = 2.3*Watt/(meter*Kelvin) # Grout thermal conductivity
λp = 0.38*Watt/(meter*Kelvin); # Pipe wall thermal conductivity

# ### BTES geometries
# We will consider two different closed-loop geometries: U1 and coaxial. In a U1
# configuration, the fluid flows down one pipe, makes a U-turn and returns to
# the surface through a parallel pipe. In a coaxial configuration, the fluid
# flows down an outer pipe and returns through an inner pipe concentrically
# located inside the outer pipe, or vice versa. In both configurations, the
# pipes are surrounded by grout material, which in turn is embedded in the rock
# formation.

L = 100.0*meter # BTES length
geo_u1 = ( # U1 geometry
    closed_loop_type = :u1,
    radius_grout = 65e-3*meter, # Grout radius
    radius_pipe = 16e-3*meter, # Pipe outer radius
    wall_thickness_pipe = 2.9e-3*meter, # Pipe wall thickness
    pipe_spacing = 60e-3*meter, # Pipe spacing
    pipe_thermal_conductivity = λp, # Pipe thermal conductivity
)
geo_coax = ( # Coaxial geometry
    closed_loop_type = :coaxial,
    radius_grout = 50e-3*meter, # Grout radius
    radius_pipe_inner = 12e-3*meter, # Inner pipe outer radius
    wall_thickness_pipe_inner = 3e-3*meter, # Inner pipe wall thickness
    radius_pipe_outer = 25e-3*meter, # Outer pipe outer radius
    wall_thickness_pipe_outer = 4e-3*meter, # Outer pipe wall thickness
    inner_pipe_thermal_conductivity = λp, # Pipe thermal conductivity
    outer_pipe_thermal_conductivity = λp # Pipe thermal conductivity
);

# ## Utility functions
# We set up a convenience function to create a BTES simulation case based on the
# specified geometry type and discretization parameters, and a function to plot
# the results of a BTES simulation against the analytical solution.

function setup_btes_single( # Utility function to set up a BTES simulation case
    type; # Type of closed-loop geometry (:u1 or :coaxial)
    nz = 100, # Number of cells in vertical direction
    n_step = 1, # Number of time steps
    inlet = :outer # Inlet pipe for coaxial geometry (:outer or :inner
    )

    ## Select well arguments based on geometry type
    if type == :u1
        well_args = geo_u1
    elseif type == :coaxial
        well_args = geo_coax
    else
        error("Unknown closed loop type: $type")
    end
    ## Create reservoir domain
    dims = (1,1,nz) 
    Δ = (100000.0*dims[3]/L, 100000.0*dims[3]/L, L)
    msh = CartesianMesh(dims, Δ)
    reservoir = reservoir_domain(msh;
        rock_density = ρr,
        rock_heat_capacity = Cpr,
        rock_thermal_conductivity = λr,
        fluid_thermal_conductivity = λf,
        component_heat_capacity = Cpf,
        porosity = ϕ,
        permeability = K
    )
    ## Set up BTES well
    wells = Fimbul.setup_vertical_btes_well(reservoir, 1, 1;
        well_args...,
        grouting_density = ρg,
        grouting_heat_capacity = Cpg,
        grouting_thermal_conductivity = λg
    )
    ## Set up reservoir model
    sys = SinglePhaseSystem(AqueousPhase(), reference_density = ρf)
    model = setup_reservoir_model(reservoir, sys; wells = wells, thermal = true)
    ## Set up controls and boundary conditions
    ctrl_inj = InjectorControl( # Injection control with fixed temperature
        TotalRateTarget(Q), [1.0], 
        density=ρf, temperature=T_in)
    ctrl_prod = ProducerControl( # Production control with fixed BHP
        BottomHolePressureTarget(5.0*atm))
    geo = tpfv_geometry(msh)
    bc_cells = geo.boundary_neighbors
    domain = reservoir_model(model).data_domain
    bc = flow_boundary_condition( # Fixed pressure and temperature BCs
        bc_cells, domain, 5.0*atm, T_rock)
    ## Set up forces and initial state
    if type == :coaxial && inlet == :outer # Coax with injection in outer pipe
        ctrl_supply = ctrl_prod
        ctrl_return = ctrl_inj
    else # Coax with injection in inner pipe or U1
        ctrl_supply = ctrl_inj
        ctrl_return = ctrl_prod
    end
    controls = Dict(
        :BTES_supply => ctrl_supply,
        :BTES_return => ctrl_return
    )
    forces = setup_reservoir_forces(model; control=controls, bc = bc)
    state0 = setup_reservoir_state( # Initialize with uniform pressure and temperature
        model; Pressure = 5*atm, Temperature = T_rock)
    ## Set time large enough to ensure steady-state
    rg = well_args.radius_grout
    time = 5/4*2*rg*(ϕ*ρf*Cpf + (1 - ϕ)*ρr*Cpr)/(ϕ*λf + (1 - ϕ)*λr)*10
    dt = fill(time/n_step, n_step)
    ## Create Jutul case
    case = JutulCase(model, dt, forces; state0 = state0)
    ## Set up simulator with temperature-based timestep selector
    sim, cfg = setup_reservoir_simulator(
        case; initial_dt = 1.0);
    sel = VariableChangeTimestepSelector(:Temperature, 2.5; 
        relative = false, model = :BTES_supply)
    push!(cfg[:timestep_selectors], sel);

    return case, sim, cfg

end

function plot_btes( # Utility function to plot BTES simulation results
    case, simulated, analytical)

    ## Create figure
    fig = Figure(size = (800, 600))
    ax = Axis(fig[1, 1], 
    title = "BTES Temperature", xlabel = "Temperature (°C)", ylabel = "Depth (m)",
    yreversed = true)
    ## Extract numerical solution
    well = case.model.models[:BTES_supply].data_domain
    sim_temp = convert_from_si.(
        simulated.result.states[end][:BTES_supply][:Temperature], :Celsius)
    za = collect(range(0, L, step=0.1)) # Analytical depth points
    section_handles = []
    solution_handles = []
    section_names = String[]
    solution_names = String[]
    ## Plot results for each well section
    sections = last.(well[:section]) |> unique |> collect
    colors = cgrad(:BrBG_4, 4, categorical=true)[[1,2,4,3]]
    err_inf = -Inf
    for (i, section) in enumerate(sections)
        ## Section name for printing
        name = titlecase(replace(string(section), "_" => " "))
        ## Plot analytical solution
        Ta = convert_from_si.(analytical[section].(za), :Celsius)
        la = lines!(ax, Ta, za; color=colors[i], linewidth = 8, linestyle = :dash, label=name)
        ## Plot numerical solution
        cells = last.(well[:section]) .== section
        Tn = sim_temp[cells]
        zn = well[:cell_centroids][3, cells]
        ln = lines!(ax, Tn, zn; color=colors[i], linewidth = 2)
        ## Store handles and names for legend
        push!(section_handles, la)
        push!(section_names, name)
        if i == 1
           push!(solution_handles, la, ln)
           push!(solution_names, "Analytical", "Numerical")
        end
        ## Compute maximum error in section
        Ta_n = convert_from_si.(analytical[section].(zn), :Celsius)
        e_inf = maximum(abs.(Tn .- Ta_n))
        err_inf = max(err_inf, e_inf)
    end
    ## Add legend
    Legend(fig[1, 2], [section_handles, solution_handles], [section_names, solution_names],
        ["Sections", "Solutions"]; orientation = :vertical)

    return fig, err_inf
    
end;

# ## Validate U1 configuration
# We set up and simulate a BTES system with U1 closed-loop geometry, and
# compare the numerical results against the analytical solution.
case_u1, sim, cfg = setup_btes_single(:u1; nz=125)
res_u1 = simulate_reservoir(case_u1; simulator = sim, config = cfg, info_level = 0)

# ### Analytical solution
# The function `analytical_closed_loop_u1` computes the analytical steady-state
# solution for a vertical U1 closed-loop system using [CITE]. The function can
# be called using all geometric and material parameters as defined above, but
# also comes in a convenience-style version that extracts the necessary
# parameters directly from the DataDomain of the BTES well model.
analytical_u1 = Fimbul.analytical_closed_loop_u1(Q, T_in, T_rock,
    ρf, Cpf, case_u1.model.models[:BTES_supply].data_domain)

# ### Plot comparison
# The numerical solution agrees very well with the analytical solution.
fig, err_inf = plot_btes(case_u1, res_u1, analytical_u1)
println("Maximum error in U1 configuration: $(round(err_inf, digits=3)) °C")
fig

# ## Validate coaxial configuration with outer pipe inlet
# Next, we set up and simulate a BTES system with coaxial closed-loop geometry,
# injecting fluid through the outer pipe, and compare the numerical results
# against the analytical solution.
case_coax_outer, sim, cfg = setup_btes_single(:coaxial; nz=125, inlet = :outer)
res_coax_outer = simulate_reservoir(case_coax_outer; simulator = sim, config = cfg, info_level = 0)

# ### Analytical solution
analytical_coax_outer = Fimbul.analytical_closed_loop_coaxial(Q, T_in, T_rock,
    ρf, Cpf, case_coax_outer.model.models[:BTES_supply].data_domain; inlet = :outer)

# ### Plot comparison
fig, err_inf = plot_btes(case_coax_outer, res_coax_outer, analytical_coax_outer)
println("Maximum error in coaxial outer configuration: $(round(err_inf, digits=3)) °C")
fig

# ## Validate coaxial configuration with inner pipe inlet
# Finally, we set up and simulate a BTES system with coaxial closed-loop
# geometry, injecting fluid through the inner pipe, and compare the numerical
# results.
case_coax_inner, sim, cfg = setup_btes_single(:coaxial; nz=125, inlet = :inner)
res_coax_inner = simulate_reservoir(case_coax_inner; simulator = sim, config = cfg, info_level = 0)

# ### Analytical solution
analytical_coax_inner = Fimbul.analytical_closed_loop_coaxial(Q, T_in, T_rock,
    ρf, Cpf, case_coax_inner.model.models[:BTES_supply].data_domain; inlet = :inner)

# ### Plot comparison
fig, err_inf = plot_btes(case_coax_inner, res_coax_inner, analytical_coax_inner)
println("Maximum error in coaxial inner configuration: $(round(err_inf, digits=3)) °C")
fig