using Jutul, JutulDarcy, Fimbul

meter = si_unit(:meter)
kilogram = si_unit(:kilogram)
atm = si_unit(:atm)
second, day = si_units(:second, :day)
Kelvin, Joule, Watt = si_units(:kelvin, :joule, :watt)
darcy = si_unit(:darcy)

## Operation conditions
Tin = convert_to_si(80.0, :Celsius)
Trock = convert_to_si(10.0, :Celsius)
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

# Geometry
L = 100.0*meter # BTES length
rg = 65e-3*meter # Grout radius
rp = 16e-3*meter # Pipe outer radius
wtp = 2.9e-3*meter # Pipe wall thickness
ps = 60e-3*meter # Pipe spacing

##
function setup_btes_single(nz = 100, n_step = 1)

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

    #
    wells = Fimbul.setup_vertical_btes_well(reservoir, 1, 1;
        closed_loop_type = :u1,
        radius_grout = rg, 
        pipe_spacing = ps,
        radius_pipe = rp,
        wall_thickness = wtp,
        grouting_density = ρg,
        grouting_heat_capacity = Cpg,
        pipe_thermal_conductivity = λp,
        grouting_thermal_conductivity = λg
    )

    #
    sys = SinglePhaseSystem(AqueousPhase(), reference_density = ρf)
    model = setup_reservoir_model(reservoir, sys; wells = wells, thermal = true)

    #
    rate_target = TotalRateTarget(Q)
    ctrl_charge = InjectorControl(rate_target, [1.0], 
        density=ρf, temperature=Tin)

    # BHP control for return side
    bhp_target = BottomHolePressureTarget(5.0*atm)
    ctrl_ret = ProducerControl(bhp_target);
    geo = tpfv_geometry(msh)
    bc_cells = geo.boundary_neighbors
    domain = reservoir_model(model).data_domain
    bc = flow_boundary_condition(bc_cells, domain, 5.0*atm, Trock)
    controls = Dict(
        :BTES_supply => ctrl_charge,
        :BTES_return => ctrl_ret
    )

    forces = setup_reservoir_forces(model; control=controls, bc = bc)
    time = 5/4*2*rg*(ϕ*ρf*Cpf + (1 - ϕ)*ρr*Cpr)/(ϕ*λf + (1 - ϕ)*λr)*10
    dt = fill(time/n_step, n_step)
    state0 = setup_reservoir_state(model; Pressure = 5*atm, Temperature = Trock)
    case = JutulCase(model, dt, forces; state0 = state0)

    return case

end

##
case = setup_btes_single(100)
sim, cfg = setup_reservoir_simulator(case;
initial_dt = 1.0);
sel = VariableChangeTimestepSelector(:Temperature, 2.5; 
    relative = false, model = :BTES_supply)
push!(cfg[:timestep_selectors], sel);
res = simulate_reservoir(case; simulator = sim, config = cfg, info_level = 2)

##
Tsol = Fimbul.analytical_closed_loop_u1(Q, Tin, Trock,
ρf, Cpf, L, rg, rp, wtp, ps, λg, λp)
zsol = collect(range(0, L, step=0.1));

##
using GLMakie

well = case.model.models[:BTES_supply].data_domain

fig = Figure(size = (800, 600))
ax = Axis(fig[1, 1]; yreversed=true)
z = collect(range(0, 100.0, step=0.1))

Tnum = convert_from_si.(res.result.states[end][:BTES_supply][:Temperature], :Celsius)
znum = well[:cell_centroids][3,:]

sprops = (color = :darkred, linewidth = 2)
aprops = (color = :black, linewidth = 8, linestyle = :dash)

plot_analytical! = (ax, section) -> begin
    lines!(ax, convert_from_si.(Tsol[section].(zsol), :Celsius), zsol; aprops...)
end

plot_numerical! = (ax, section) -> begin 
    cells = last.(well[:section]) .== section
    lines!(ax, Tnum[cells], znum[cells]; sprops...)
end

for section in (:pipe_supply, :pipe_return, :grout_supply, :grout_return)
    plot_analytical!(ax, section)
    plot_numerical!(ax, section)
end

fig

##
sol = Fimbul.temperature_u1(Q, Tin, Trock,
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