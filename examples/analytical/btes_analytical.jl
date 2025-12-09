using Jutul, JutulDarcy, Fimbul

meter = si_unit(:meter)
kilogram = si_unit(:kilogram)
atm = si_unit(:atm)
second, day = si_units(:second, :day)
Kelvin, Joule, Watt = si_units(:kelvin, :joule, :watt)
darcy = si_unit(:darcy)

## Operating conditions
T_in = convert_to_si(80.0, :Celsius)
T_rock = convert_to_si(10.0, :Celsius)
Q = 21.86*meter^3/day

# Fluid
ρf = 988.1*kilogram/meter^3
Cpf = 4184.0*Joule/(kilogram*Kelvin)
λf = 0.6405*Watt/(meter*Kelvin)

# Rock
ϕ = 0.01
K = 1e-3*darcy
ρr = 2650.0*kilogram/meter^3
Cpr = 900.0*Joule/(kilogram*Kelvin)
λr = 2.5*Watt/(meter*Kelvin)

# Materials
ρg = 2000.0*kilogram/meter^3 # Grout
Cpg = 1500.0*Joule/(kilogram*Kelvin) # Grout
λg = 2.3*Watt/(meter*Kelvin) # Grout
λp = 0.38*Watt/(meter*Kelvin) # Pipe

# BTES geometries
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
)

##
function setup_btes_single(type; nz = 100, n_step = 1, inlet = :outer)

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

    if type == :u1
        well_args = geo_u1
    elseif type == :coaxial
        well_args = geo_coax
    else
        error("Unknown closed loop type: $type")
    end
    wells = Fimbul.setup_vertical_btes_well(reservoir, 1, 1;
        well_args...,
        grouting_density = ρg,
        grouting_heat_capacity = Cpg,
        grouting_thermal_conductivity = λg
    )

    #
    sys = SinglePhaseSystem(AqueousPhase(), reference_density = ρf)
    model = setup_reservoir_model(reservoir, sys; wells = wells, thermal = true)

    ## Injection control
    rate_target = TotalRateTarget(Q)
    ctrl_inj = InjectorControl(rate_target, [1.0], 
        density=ρf, temperature=T_in)
    ## BHP control for return side
    bhp_target = BottomHolePressureTarget(5.0*atm)
    ctrl_prod = ProducerControl(bhp_target);
    geo = tpfv_geometry(msh)
    bc_cells = geo.boundary_neighbors
    domain = reservoir_model(model).data_domain
    bc = flow_boundary_condition(bc_cells, domain, 5.0*atm, T_rock)

    if type == :coaxial && inlet == :outer
        ctrl_supply = ctrl_prod
        ctrl_return = ctrl_inj
    else
        ctrl_supply = ctrl_inj
        ctrl_return = ctrl_prod
    end
    controls = Dict(
        :BTES_supply => ctrl_supply,
        :BTES_return => ctrl_return
    )

    forces = setup_reservoir_forces(model; control=controls, bc = bc)
    rg = well_args.radius_grout
    time = 5/4*2*rg*(ϕ*ρf*Cpf + (1 - ϕ)*ρr*Cpr)/(ϕ*λf + (1 - ϕ)*λr)*10
    dt = fill(time/n_step, n_step)
    state0 = setup_reservoir_state(model; Pressure = 5*atm, Temperature = T_rock)
    case = JutulCase(model, dt, forces; state0 = state0)

    return case

end

##
function plot_btes(case, simulated, analytical)

    fig = Figure(size = (800, 600))
    ax = Axis(fig[1, 1], 
    title = "BTES Temperature", xlabel = "Temperature", ylabel = "Depth",
    yreversed = true)

    well = case.model.models[:BTES_supply].data_domain
    za = analytical[:z]
    sections = last.(well[:section]) |> unique |> collect
    # colors = cgrad(:Paired_10, 10, categorical=true)
    # colors = cgrad(:BrBg, 5, categorical=true)[[1,2,end-1,end]]
    colors = cgrad(:BrBG_4, 4, categorical=true)[[1,2,4,3]]
    sim_temp = convert_from_si.(
        simulated.result.states[end][:BTES_supply][:Temperature], :Celsius)
    section_handles = []
    solution_handles = []
    section_names = String[]
    solution_names = String[]
    for (i, section) in enumerate(sections)
        # Section name for printing
        name = titlecase(replace(string(section), "_" => " "))
        # Plot analytical solution
        Ta = convert_from_si.(analytical[section].(za), :Celsius)
        la = lines!(ax, Ta, za; color=colors[i], linewidth = 8, linestyle = :dash, label=name)
        # Plot numerical solution
        cells = last.(well[:section]) .== section
        Tn = sim_temp[cells]
        zn = well[:cell_centroids][3, cells]
        ln = lines!(ax, Tn, zn; color=colors[i], linewidth = 2)
        push!(section_handles, la)
        push!(section_names, name)
        if i == 1
           push!(solution_handles, la, ln)
           push!(solution_names, "Analytical", "Numerical")
        end
    end
    Legend(fig[1, 2], [section_handles, solution_handles], [section_names, solution_names],
        ["Sections", "Solutions"]; orientation = :vertical)


    return fig
    
end

##
case_u1 = setup_btes_single(:u1; nz=100)
sim, cfg = setup_reservoir_simulator(case;
initial_dt = 1.0);
sel = VariableChangeTimestepSelector(:Temperature, 2.5; 
    relative = false, model = :BTES_supply)
push!(cfg[:timestep_selectors], sel);
res_u1 = simulate_reservoir(case; simulator = sim, config = cfg, info_level = 2)

##
analytical_u1 = Fimbul.analytical_closed_loop_u1(Q, T_in, T_rock,
    ρf, Cpf, L, geo_u1.radius_grout, geo_u1.radius_pipe,
    geo_u1.wall_thickness_pipe, geo_u1.pipe_spacing, λg, geo_u1.pipe_thermal_conductivity)
analytical_u1[:z] = collect(range(0, L, step=0.1));

##
fig = plot_btes(case_u1, res_u1, analytical_u1)

##
case_coax_outer = setup_btes_single(:coaxial; nz=100, inlet = :outer)
sim, cfg = setup_reservoir_simulator(case_coax_outer;
initial_dt = 1.0);
sel = VariableChangeTimestepSelector(:Temperature, 2.5; 
    relative = false, model = :BTES_supply)
push!(cfg[:timestep_selectors], sel);
res_coax_outer = simulate_reservoir(case_coax_outer; simulator = sim, config = cfg, info_level = 2)

##
analytical_coax_outer = Fimbul.analytical_closed_loop_coaxial(Q, T_in, T_rock,
    ρf, Cpf, L, geo_coax.radius_grout,
    geo_coax.radius_pipe_inner, geo_coax.wall_thickness_pipe_inner,
    geo_coax.radius_pipe_outer, geo_coax.wall_thickness_pipe_outer,
    λg, geo_coax.inner_pipe_thermal_conductivity, geo_coax.outer_pipe_thermal_conductivity,
    inlet = :outer)
analytical_coax_outer[:z] = collect(range(0, L, step=0.1));

##
plot_btes(case_coax_outer, res_coax_outer, analytical_coax_outer)

##
case_coax_inner = setup_btes_single(:coaxial; nz=100, inlet = :inner)
sim, cfg = setup_reservoir_simulator(case_coax_inner;
initial_dt = 1.0);
sel = VariableChangeTimestepSelector(:Temperature, 2.5; 
    relative = false, model = :BTES_supply)
push!(cfg[:timestep_selectors], sel);
res_coax_inner = simulate_reservoir(case_coax_inner; simulator = sim, config = cfg, info_level = 2)

##
analytical_coax_inner = Fimbul.analytical_closed_loop_coaxial(Q, T_in, T_rock,
    ρf, Cpf, L, geo_coax.radius_grout,
    geo_coax.radius_pipe_inner, geo_coax.wall_thickness_pipe_inner,
    geo_coax.radius_pipe_outer, geo_coax.wall_thickness_pipe_outer,
    λg, geo_coax.inner_pipe_thermal_conductivity, geo_coax.outer_pipe_thermal_conductivity,
    inlet = :inner)
analytical_coax_inner[:z] = collect(range(0, L, step=0.1));

##
plot_btes(case_coax_inner, res_coax_inner, analytical_coax_inner)