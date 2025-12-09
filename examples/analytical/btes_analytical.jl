using Jutul, JutulDarcy, Fimbul

meter = si_unit(:meter)
kilogram = si_unit(:kilogram)
atm = si_unit(:atm)
second, day = si_units(:second, :day)
Kelvin, Joule, Watt = si_units(:kelvin, :joule, :watt)
darcy = si_unit(:darcy)

## Operation conditions
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
    time = 5/4*2*rg*(ϕ*ρf*Cpf + (1 - ϕ)*ρr*Cpr)/(ϕ*λf + (1 - ϕ)*λr)*10
    dt = fill(time/n_step, n_step)
    state0 = setup_reservoir_state(model; Pressure = 5*atm, Temperature = T_rock)
    case = JutulCase(model, dt, forces; state0 = state0)

    return case

end

##
function plot_btes!(ax, case, simulated, analytical)

    well = case.model.models[:BTES_supply].data_domain
    za = analytical[:z]
    sections = last.(well[:section]) |> unique |> collect

    sprops = (color = :darkred, linewidth = 2)
    aprops = (color = :black, linewidth = 8, linestyle = :dash)

    sim_temp = convert_from_si.(
        simulated.result.states[end][:BTES_supply][:Temperature], :Celsius)
    for section in sections
        
        Ta = convert_from_si.(analytical[section].(za), :Celsius)
        lines!(ax, Ta, za; aprops...)
        
        cells = last.(well[:section]) .== section
        Tn = sim_temp[cells]
        zn = well[:cell_centroids][3, cells]
        lines!(ax, Tn, zn; sprops...)
    end

end

##
case = setup_btes_single(:u1; nz=100)
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
fig = Figure(size = (800, 600))
ax = Axis(fig[1, 1], title = "BTES U1 Temperature", xlabel = "Temperature", ylabel = "Depth", yreversed = true)
plot_btes!(ax, case, res, analytical_u1)
fig

##
case = setup_btes_single(:coaxial; nz=100, inlet = :outer)
sim, cfg = setup_reservoir_simulator(case;
initial_dt = 1.0);
sel = VariableChangeTimestepSelector(:Temperature, 2.5; 
    relative = false, model = :BTES_supply)
push!(cfg[:timestep_selectors], sel);
res_coax_outer = simulate_reservoir(case; simulator = sim, config = cfg, info_level = 2)

##
analytical_coax_outer = Fimbul.analytical_closed_loop_coaxial(Q, T_in, T_rock,
    ρf, Cpf, L, geo_coax.radius_grout,
    geo_coax.radius_pipe_inner, geo_coax.wall_thickness_pipe_inner,
    geo_coax.radius_pipe_outer, geo_coax.wall_thickness_pipe_outer,
    λg, geo_coax.inner_pipe_thermal_conductivity, geo_coax.outer_pipe_thermal_conductivity,
    inlet = :outer)
analytical_coax_outer[:z] = collect(range(0, L, step=0.1));

##
fig = Figure(size = (800, 600))
ax = Axis(fig[1, 1], title = "BTES Coaxial Temperature", xlabel = "Temperature", ylabel = "Depth", yreversed = true)
plot_btes!(ax, case, res_coax_outer, analytical_coax_outer)
fig

##
case = setup_btes_single(:coaxial; nz=100, inlet = :inner)
sim, cfg = setup_reservoir_simulator(case;
initial_dt = 1.0);
sel = VariableChangeTimestepSelector(:Temperature, 2.5; 
    relative = false, model = :BTES_supply)
push!(cfg[:timestep_selectors], sel);
res_coax_inner = simulate_reservoir(case; simulator = sim, config = cfg, info_level = 2)

##
analytical_coax_inner = Fimbul.analytical_closed_loop_coaxial(Q, T_in, T_rock,
    ρf, Cpf, L, geo_coax.radius_grout,
    geo_coax.radius_pipe_inner, geo_coax.wall_thickness_pipe_inner,
    geo_coax.radius_pipe_outer, geo_coax.wall_thickness_pipe_outer,
    λg, geo_coax.inner_pipe_thermal_conductivity, geo_coax.outer_pipe_thermal_conductivity,
    inlet = :inner)
analytical_coax_inner[:z] = collect(range(0, L, step=0.1));

##
fig = Figure(size = (800, 600))
ax = Axis(fig[1, 1], title = "BTES Coaxial Temperature", xlabel = "Temperature", ylabel = "Depth", yreversed = true)
plot_btes!(ax, case, res_coax_inner, analytical_coax_inner)
fig









##
sol = Fimbul.temperature_u1(Q, T_in, T_rock,
    ρf, Cpf, L, rg, rp, wtp, ps, λg, λp)
    

function compute_error(res, case, Tsol)
    well = case.model.models[:BTES_supply].data_domain
    Tnum = convert_from_si.(res.result.states[end][:BTES_supply][:Temperature], :Celsius)
    znum = well[:cell_centroids][3,:]
    dz = well[:cell_length]

    err_l2 = 0.0
    err_inf = -Inf
    for section in (:pipe_supply, :pipe_return, :grout_supply, :grout_return)
        cells = last.(well[:section]) .== section
        # dz = diff(vcat(0.0, zs))
        Ts = convert_from_si.(Tsol[section].(znum[cells]), :Celsius)
        Tn = Tnum[cells]
        dz_s = dz[cells]
        e_l2 = sqrt(sum((Tn .- Ts).^2.0.*dz_s)./sum(dz_s))
        e_inf = maximum(abs.(Tn .- Ts))
        err_l2 += e_l2
        err_inf = max(err_inf, e_inf)
        println("Max error in section $section: $e_inf °C")
    end
    return err_l2, err_inf
end

errors = Dict(:l2 => Float64[], :inf => Float64[])
cell_size = Float64[]
num_cells = 5*2.0.^(0:5)
for nz = num_cells
    println("Running with nz = $nz")
    case = setup_btes_single(Int(nz))
    sim, cfg = setup_reservoir_simulator(case;
    initial_dt = 1.0);
    sel = VariableChangeTimestepSelector(:Temperature, 5.0; 
        relative = false, model = :BTES_supply)
    push!(cfg[:timestep_selectors], sel);
    res = simulate_reservoir(case; simulator = sim, config = cfg, info_level = 1)
    
    err_l2, err_inf = compute_error(res, case, Tsol)
    push!(errors[:l2], err_l2)
    push!(errors[:inf], err_inf)
    push!(cell_size, L/Int(nz))
end

##
fig = Figure(size = (800, 600))
ax = Axis(fig[1, 1]; xreversed=true, xscale = log2, yscale = log2, aspect = AxisAspect(1),
title = "L² Error", xlabel = "Cell size (m)")
lines!(ax, cell_size, errors[:l2]; linewidth = 2, color = :black)
scatter!(ax, cell_size, errors[:l2]; markersize = 10, color = :black)

ax = Axis(fig[1, 2]; xreversed=true, xscale = log2, yscale = log2, aspect = AxisAspect(1),
title = "L∞ Error", xlabel = "Cell size (m)")
lines!(ax, cell_size, errors[:inf]; linewidth = 2, color = :black)
scatter!(ax, cell_size, errors[:inf]; markersize = 10, color = :black)
fig

##
rg = 50e-3
rp_out = 25e-3
rp_in = 12e-3
wtp_out = 4.0e-3
wtp_in = 3e-3
sol = Fimbul.analytical_closed_loop_coaxial(Q, Tin, Trock,
    ρf, Cpf, L, rg, rp_in, wtp_in, rp_out, wtp_out, λg, λp, λp; inlet = :outer)


fig = Figure(size = (800, 600))
ax = Axis(fig[1, 1]; yreversed=true)
aprops = (color = :black, linewidth = 8, linestyle = :dash)

plot_analytical! = (ax, section) -> begin
    lines!(ax, convert_from_si.(sol[section].(zsol), :Celsius), zsol; aprops...)
end

for section in (:pipe_inner, :pipe_outer, :grout)
    plot_analytical!(ax, section)
end
fig