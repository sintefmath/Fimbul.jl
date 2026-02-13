"""
    ags(; kwargs...)

Create an Advanced Geothermal System (AGS) simulation case.

# Keyword Arguments

## Well trajectory parameters
- `well_coords`: A vector of matrices specifying the (x,y,z) coordinates of each
  well section. Default is the AGS well trajectory given by
  `get_ags_trajectory()`.
- `well_connectivity`: A matrix specifying the connectivity between well sections.
  Each row corresponds to a well section, with the first column indicating the
  index of the section from which it receives flow (0 if none), and the second
  column indicating the index of the section to which it sends flow (0 if none).
  Default is the connectivity for the AGS well trajectory given by
  `get_ags_trajectory()`. Note that each connection should only be specified once.

## Geological parameters
- `depths`: Vector of depths defining layer boundaries [m]
- `porosity`: Porosity as scalar or one value for each layer [-]
- `permeability`: Permeability as scalar or one value for each layer [m²]
- `rock_thermal_conductivity`: Thermal conductivity as scalar or one value for each layer [W/(m·K)]
- `rock_heat_capacity`: Rock heat capacity as scalar or one value for each layer [J/(kg·K)]

## Thermal Parameters
- `temperature_surface = 10°C`: Surface temperature [K]
- `thermal_gradient = 0.03`: Geothermal gradient [K/m]

## Operational Parameters
- `rate = 25 m³/h`: Injection/production rate [m³/s]
- `temperature_inj = 25°C`: Injection temperature [K]
- `num_years = 50`: Total simulation time [years]
- `report_interval = year/4`: Reporting interval [s]
- `schedule_args`: Additional arguments passed to make_schedule

## Mesh Parameters
- `hz`: Vector of vertical cell sizes for each layer [m]. If `missing`, a default
  sizing is used.
- `hxy_min = 25.0`: Minimum horizontal cell size [m]
- `hxy_max = 250.0`: Maximum horizontal cell size [m]
- `mesh_args`: Additional arguments passed to Fimbul.extruded_mesh

# Returns
A `JutulCase` for AGS
```
"""
function ags(;
    well_coords = first(get_ags_trajectory()),
    well_connectivity = last(get_ags_trajectory()),
    depths = [0.0, 1500.0, 2300, 2400, 2800],
    porosity = 0.01,
    permeability = 1e-3 .* darcy,
    rock_thermal_conductivity = 2.5 .* watt/(meter*Kelvin),
    rock_heat_capacity = 900.0*joule/(kilogram*Kelvin),
    temperature_surface = convert_to_si(10, :Celsius),
    thermal_gradient = 0.03Kelvin/meter,
    rate = 25meter^3/hour,
    temperature_inj = convert_to_si(25, :Celsius),
    num_years = 50,
    report_interval = year/4,
    schedule_args = NamedTuple(),
    hz = missing,
    hxy_min = 25.0,
    hxy_max = 250.0,
    mesh_args = NamedTuple()
)

    # ## Set up mesh
    if ismissing(hz)
        hz = fill(250.0, length(depths)-1)
        hz[end-1] = 10.0
        dz = diff(depths)
        n = dz./hz
        ok = n .>= 5
        hz[.!ok] = dz[.!ok]./5
    end
    
    constraints = get_ags_constraints(well_coords; hxy_min = hxy_min)

    domain, layers, metrics = layered_reservoir_domain(constraints, depths;
        mesh_args = (;
            hz = hz, 
            hxy_min = hxy_min, 
            hxy_max = hxy_max, 
            offset_rel = 1.0, 
            dist_min_factor = 50.0, 
            mesh_args...
        ),
        layer_properties = (
            porosity = porosity,
            permeability = permeability,
            rock_thermal_conductivity = rock_thermal_conductivity,
            rock_heat_capacity = rock_heat_capacity
        )
    )

    # ## Set up model
    wells, section_info = setup_ags_wells(domain, well_coords, well_connectivity)

    model = setup_reservoir_model(
        domain, :geothermal; wells = wells)
    bc, state0 = set_dirichlet_bcs(model;
        pressure_surface = 10atm,
        temperature_surface = temperature_surface,
        geothermal_gradient = thermal_gradient,
    )

    rho = reservoir_model(model).system.rho_ref[1]
    # Injector: inject cooled water at a rate equal to production
    rate_target = TotalRateTarget(rate)
    ctrl_supply = InjectorControl(
        rate_target, [1.0], density = rho, temperature = temperature_inj)
    # Producer
    # TODO: support multiple periods with different rates
    bhp_target = BottomHolePressureTarget(10atm)
    ctrl_return = ProducerControl(bhp_target)
    control = Dict(:AGS_supply => ctrl_supply, :AGS_return => ctrl_return)

    forces = setup_reservoir_forces(model, control = control, bc = bc)#, limits = limits)

    dt, forces = make_schedule(
        [forces], [(1,1), (1,1)];
        num_years = num_years,
        report_interval = report_interval,
        schedule_args...)

    # ## Additional case info
    info = Dict(
        :description => "AGS closed-loop geothermal system",
        :wells_coords => well_coords,
        :sections => section_info,
    )

    # ## Create and return the complete simulation case
    case = JutulCase(model, dt, forces; state0 = state0, input_data = info)

    return case

end

function get_ags_constraints(well_coords; hxy_min)

    Δ = hxy_min/2
    well_coords_2x = []
    for wc in well_coords
        wc_left = copy(wc)
        wc_left[:, 2] .-= Δ/2
        push!(well_coords_2x, wc_left)
        wc_right = copy(wc)
        wc_right[:, 2] .+= Δ/2
        push!(well_coords_2x, wc_right)
    end

    cell_constraints = []
    for wc in well_coords_2x
        xw = [Tuple(x) for x in eachrow(unique(wc[:, 1:2], dims=1))]
        for cc in cell_constraints
            cor = [norm(collect(x) .- collect(y), 2) for x in xw, y in cc]
            xw = [xw[i] for i in eachindex(xw) if all(cor[i, :] .>= Δ)]
        end
        push!(cell_constraints, xw)
    end

    return cell_constraints

end

function setup_ags_wells(domain, well_coords, well_connectivity)

    msh = physical_representation(domain)
    geo = tpfv_geometry(msh)
    wells = []
    directions = []
    wcells = []
    rcells = []
    neighbors = []
    cell_to_section = []
    wc0 = 0
    for (k, wc) in enumerate(well_coords)
        println("Processing well section $k/$(length(well_coords))")

        rc, extra = Jutul.find_enclosing_cells(
            msh, wc, n = 1000; geometry=geo, extra_out = true)
        dir = Vector.(extra[:direction].*extra[:lengths])

        wc = collect(1:length(rc)) .+ wc0
        wc0 += length(wc)
        push!(rcells, rc)
        push!(wcells, wc)
        push!(directions, dir)
        push!(neighbors, vcat(wc[1:end-1]', wc[2:end]'))
        push!(cell_to_section, fill(k, length(rc)))
        println("  Number of cells in section $k: ", length(rc))
    end

    connection_segments = []

    seg_no = 0
    for k = 1:length(well_coords)
        from = well_connectivity[k,1]
        if from != 0
            ix = findfirst(rcells[from] .== rcells[k][1])
            if isnothing(ix)
                error("Could not find connection from well section $from to well $k")
            end
            in_seg = [wcells[from][ix], wcells[k][1]]
            neighbors[k] = hcat(in_seg, neighbors[k])
        end
        to = well_connectivity[k,2]
        if to != 0
            ix = findfirst(rcells[to] .== rcells[k][end])
            if isnothing(ix)
                error("Could not find connection from well section $k to well $to")
            end
            out_seg = [wcells[k][end], wcells[to][ix]]
            neighbors[k] = hcat(neighbors[k], out_seg)
        end
        push!(connection_segments, [1, size(neighbors[k], 2)].+seg_no)
        seg_no += size(neighbors[k], 2)
    end

    rcells = vcat(rcells...)
    wcells = vcat(wcells...)
    directions = vcat(directions...)
    neighbors = hcat(neighbors...)
    cell_to_section = vcat(cell_to_section...)

    # return rcells, wcells, directions, neighbors

    common_args = (
        simple_well = false,
        type = :closed_loop,
        radius = 150e-3,
        WI = 0.0
    )

    w_supply = setup_well(domain, rcells;
        neighborship = neighbors,
        perforation_cells_well = wcells,
        dir = directions,
        well_cell_centers = geo.cell_centroids[:, rcells],
        end_nodes = [length(wcells)],
        name = :AGS_supply,
        common_args...
    )

    w_return = setup_well(domain, rcells[end];
        name = :AGS_return,
        WIth = 0.0,
        common_args...
    )

    wells = [w_supply, w_return]

    section_info = (
        cell_to_section = cell_to_section,
        connection_segments = connection_segments
    )

    return wells, section_info

end

function get_section_data_ags(case, states, well=missing)

    model = case.model
    if ismissing(well)
        model_names = keys(model.models)
        supply_ix = findfirst(map(k -> contains(string(k), "_supply"), model_names))
        well = model_names[supply_ix]
    end

    section_info = case.input_data[:sections]

    section_data = Dict()

    N = model.models[well].data_domain.representation.neighborship
    cell_to_section = section_info.cell_to_section

    num_sections = maximum(cell_to_section)

    for section = 1:num_sections
        cells = findall(section .== cell_to_section)
        segs = section_info.connection_segments[section]
        @assert any(N[:,segs[1]] .== cells[1])
        @assert any(N[:,segs[2]] .== cells[end])

        sd = Dict(
            :Temperature => [],
            :TotalMassFlux => [],
            :Power => [],
        )

        for state in states

            add_section_data_ags!(sd, state, well, cells, segs[1], segs[2])
            if haskey(state, :substates)
                map(substate -> add_section_data_ags!(
                    sd, substate, well, cells, segs[1], segs[2]),
                    state[:substates])
            end

        end
        
        for (key, _) in sd
            if !haskey(section_data, key)
                section_data[key] = Matrix{Float64}(
                    undef, length(sd[:Temperature]), num_sections)
            end
            section_data[key][:, section] = sd[key]
        end
        
    end

    return section_data

end

function add_section_data_ags!(data, state, well, cells, seg_in, seg_out)
    T = state[well][:Temperature][cells][end]
    h_in = state[well][:FluidEnthalpy][cells][1]
    h_out = state[well][:FluidEnthalpy][cells][end]

    q_in = state[well][:TotalMassFlux][seg_in]
    q_out = state[well][:TotalMassFlux][seg_out]
    power = q_out*h_out-q_in*h_in

    push!(data[:Temperature], T)
    push!(data[:TotalMassFlux], q_out)
    push!(data[:Power], power)
end

function get_ags_trajectory()
    # AGS well setup provided by Alexander Rath (OMV)
    injector = [
        0	0	0;
        0	0	-25;
        0	0	-50;
        0	0	-75;
        0	0	-100;
        0	0	-125;
        0	0	-150;
        0	0	-175;
        0	0	-200;
        0	0	-225;
        0	0	-250;
        0	0	-275;
        0	0	-300;
        0	0	-325;
        0	0	-350;
        0	0	-375;
        0	0	-400;
        0	0	-425;
        0	0	-450;
        0	0	-475;
        0	0	-500;
        0	0	-525;
        0	0	-550;
        0	0	-575;
        0	0	-600;
        0	0	-625;
        0	0	-650;
        0	0	-675;
        0	0	-700;
        0	0	-725;
        0	0	-750;
        0	0	-775;
        0	0	-800;
        0	0	-825;
        0	0	-850;
        0	0	-875;
        0	0	-900;
        0	0	-925;
        0	0	-950;
        0	0	-975;
        0	0	-1000;
        0	0	-1025;
        0	0	-1050;
        0	0	-1075;
        0	0	-1100;
        0	0	-1125;
        0	0	-1150;
        0	0	-1175;
        0	0	-1200;
        0	0	-1225;
        0	0	-1250;
        0	0	-1275;
        0	0	-1300;
        0	0	-1325;
        0	0	-1350;
        0	0	-1375;
        0	0	-1400;
        0	0	-1425;
        0	0	-1450;
        0	0	-1475;
        0	0	-1500;
        0	0	-1525;
        0	0	-1550;
        0	0	-1575;
        0	0	-1600;
        0	0	-1625;
        0	0	-1650;
        0	0	-1675;
        0	0	-1700;
        0	0	-1725;
        0	0	-1750;
        0	0	-1775;
        0	0	-1800;
        0	0	-1825;
        0	0	-1850;
        0	0	-1875;
        0	0	-1900;
        0	0	-1925;
        0	0	-1950;
        0	0	-1975;
        0	0	-2000;
        0	0	-2025;
        0	0	-2050;
        0	0	-2075;
        0	0	-2100;
        0	0	-2125;
        0	0	-2150;
        0	0	-2175;
        0	0	-2200;
        0	0	-2225;
        0	0	-2250;
        0	0	-2275;
        0	0	-2300;
        0	0	-2325;
        0	0	-2350;
        0	0	-2375;
        2	0	-2398;
        27	0	-2398;
        52	0	-2398;
        77	0	-2398;
        100.410000000149	0	-2398;
    ]

    lateral1 = [
        100.410000000149	0	-2398;
        120.785909930244	-14.4852440259419	-2398;
        144.313397579826	-22.4170780730201	-2398;
        168.350060170516	-29.2901306319982	-2398;
        192.386722760275	-36.1631831899285	-2398;
        216.938165090047	-40.5950822109589	-2398;
        241.612873909995	-44.5655155599816	-2398;
        266.173622090369	-49.2313111869153	-2398;
        290.947639769875	-52.2590283619938	-2398;
        315.915907430463	-53.1500000000233	-2398;
        340.915907430463	-53.1500000000233	-2398;
        365.915907430463	-53.1500000000233	-2398;
        390.915907430463	-53.1500000000233	-2398;
        415.915907430463	-53.1500000000233	-2398;
        440.915907430463	-53.1500000000233	-2398;
        465.915907430463	-53.1500000000233	-2398;
        490.915907430463	-53.1500000000233	-2398;
        515.915907430463	-53.1500000000233	-2398;
        540.915907430463	-53.1500000000233	-2398;
        565.915907430463	-53.1500000000233	-2398;
        590.915907430463	-53.1500000000233	-2398;
        615.915907430463	-53.1500000000233	-2398;
        640.915907430463	-53.1500000000233	-2398;
        665.915907430463	-53.1500000000233	-2398;
        690.915907430463	-53.1500000000233	-2398;
        715.915907430463	-53.1500000000233	-2398;
        740.915907430463	-53.1500000000233	-2398;
        765.915907430463	-53.1500000000233	-2398;
        790.915907430463	-53.1500000000233	-2398;
        815.915907430463	-53.1500000000233	-2398;
        840.915907430463	-53.1500000000233	-2398;
        865.915907430463	-53.1500000000233	-2398;
        876.839999999851	-53.1500000000233	-2398;
        901.839999999851	-53.1500000000233	-2398;
        926.839999999851	-53.1500000000233	-2398;
        951.839999999851	-53.1500000000233	-2398;
        976.839999999851	-53.1500000000233	-2398;
        1001.83999999985	-53.1500000000233	-2398;
        1026.83999999985	-53.1500000000233	-2398;
        1051.83999999985	-53.1500000000233	-2398;
        1076.83999999985	-53.1500000000233	-2398;
        1101.83999999985	-53.1500000000233	-2398;
        1126.83999999985	-53.1500000000233	-2398;
        1151.83999999985	-53.1500000000233	-2398;
        1176.83999999985	-53.1500000000233	-2398;
        1201.83999999985	-53.1500000000233	-2398;
        1226.83999999985	-53.1500000000233	-2398;
        1251.83999999985	-53.1500000000233	-2398;
        1276.83999999985	-53.1500000000233	-2398;
        1301.83999999985	-53.1500000000233	-2398;
        1326.83999999985	-53.1500000000233	-2398;
        1351.83999999985	-53.1500000000233	-2398;
        1376.83999999985	-53.1500000000233	-2398;
        1401.83999999985	-53.1500000000233	-2398;
        1426.41214019991	-49.0466844260227	-2398;
        1450.87435765006	-43.8891736490186	-2398;
        1475.3365751002	-38.7316628720146	-2398;
        1499.79879255034	-33.5741520950105	-2398;
        1522.40766253043	-23.0685325369705	-2398;
        1544.81573797017	-11.9834905809257	-2398;
        1567.22381340992	-0.898448625928722	-2398;
        1569.04000000004	0	-2398;
    ]

    lateral2 = [
        100.410000000149	0	-2398;
        127	0	-2398;
        152	0	-2398;
        177	0	-2398;
        202	0	-2398;
        227	0	-2398;
        252	0	-2398;
        277	0	-2398;
        302	0	-2398;
        327	0	-2398;
        352	0	-2398;
        377	0	-2398;
        402	0	-2398;
        427	0	-2398;
        452	0	-2398;
        477	0	-2398;
        502	0	-2398;
        527	0	-2398;
        552	0	-2398;
        577	0	-2398;
        602	0	-2398;
        627	0	-2398;
        652	0	-2398;
        677	0	-2398;
        702	0	-2398;
        727	0	-2398;
        752	0	-2398;
        777	0	-2398;
        802	0	-2398;
        827	0	-2398;
        852	0	-2398;
        875.839999999851	0	-2398;
        901.839999999851	0	-2398;
        926.839999999851	0	-2398;
        951.839999999851	0	-2398;
        976.839999999851	0	-2398;
        1001.83999999985	0	-2398;
        1026.83999999985	0	-2398;
        1051.83999999985	0	-2398;
        1076.83999999985	0	-2398;
        1101.83999999985	0	-2398;
        1126.83999999985	0	-2398;
        1151.83999999985	0	-2398;
        1176.83999999985	0	-2398;
        1201.83999999985	0	-2398;
        1226.83999999985	0	-2398;
        1251.83999999985	0	-2398;
        1276.83999999985	0	-2398;
        1301.83999999985	0	-2398;
        1326.83999999985	0	-2398;
        1351.83999999985	0	-2398;
        1376.83999999985	0	-2398;
        1401.83999999985	0	-2398;
        1426.83999999985	0	-2398;
        1451.83999999985	0	-2398;
        1476.83999999985	0	-2398;
        1501.83999999985	0	-2398;
        1526.83999999985	0	-2398;
        1569.04000000004	0	-2398;
    ]

    producer = [
        1569.04000000004	0	-2398;
        1576.83999999985	0	-2398;
        1601.83999999985	0	-2398;
        1626.83999999985	0	-2398;
        1651.83999999985	0	-2398;
        1676.83999999985	0	-2398;
        1701.83999999985	0	-2398;
        1703.49000000022	0	-2374.65000000037;
        1703.49000000022	0	-2349.65000000037;
        1703.49000000022	0	-2324.65000000037;
        1703.49000000022	0	-2299.65000000037;
        1703.49000000022	0	-2274.65000000037;
        1703.49000000022	0	-2249.65000000037;
        1703.49000000022	0	-2224.65000000037;
        1703.49000000022	0	-2199.65000000037;
        1703.49000000022	0	-2174.65000000037;
        1703.49000000022	0	-2149.65000000037;
        1703.49000000022	0	-2124.65000000037;
        1703.49000000022	0	-2099.65000000037;
        1703.49000000022	0	-2074.65000000037;
        1703.49000000022	0	-2049.65000000037;
        1703.49000000022	0	-2024.65000000037;
        1703.49000000022	0	-1999.65000000037;
        1703.49000000022	0	-1974.65000000037;
        1703.49000000022	0	-1949.65000000037;
        1703.49000000022	0	-1924.65000000037;
        1703.49000000022	0	-1899.65000000037;
        1703.49000000022	0	-1874.65000000037;
        1703.49000000022	0	-1849.65000000037;
        1703.49000000022	0	-1824.65000000037;
        1703.49000000022	0	-1799.65000000037;
        1703.49000000022	0	-1774.65000000037;
        1703.49000000022	0	-1749.65000000037;
        1703.49000000022	0	-1724.65000000037;
        1703.49000000022	0	-1699.65000000037;
        1703.49000000022	0	-1674.65000000037;
        1703.49000000022	0	-1649.65000000037;
        1703.49000000022	0	-1624.65000000037;
        1703.49000000022	0	-1599.65000000037;
        1703.49000000022	0	-1574.65000000037;
        1703.49000000022	0	-1549.65000000037;
        1703.49000000022	0	-1524.65000000037;
        1703.49000000022	0	-1499.65000000037;
        1703.49000000022	0	-1474.65000000037;
        1703.49000000022	0	-1449.65000000037;
        1703.49000000022	0	-1424.65000000037;
        1703.49000000022	0	-1399.65000000037;
        1703.49000000022	0	-1374.65000000037;
        1703.49000000022	0	-1349.65000000037;
        1703.49000000022	0	-1324.65000000037;
        1703.49000000022	0	-1299.65000000037;
        1703.49000000022	0	-1274.65000000037;
        1703.49000000022	0	-1249.65000000037;
        1703.49000000022	0	-1224.65000000037;
        1703.49000000022	0	-1199.65000000037;
        1703.49000000022	0	-1174.65000000037;
        1703.49000000022	0	-1149.65000000037;
        1703.49000000022	0	-1124.65000000037;
        1703.49000000022	0	-1099.65000000037;
        1703.49000000022	0	-1074.65000000037;
        1703.49000000022	0	-1049.65000000037;
        1703.49000000022	0	-1024.65000000037;
        1703.49000000022	0	-999.650000000372;
        1703.49000000022	0	-974.65000000037;
        1703.49000000022	0	-949.650000000373;
        1703.49000000022	0	-924.650000000372;
        1703.49000000022	0	-899.65000000037;
        1703.49000000022	0	-874.650000000373;
        1703.49000000022	0	-849.650000000372;
        1703.49000000022	0	-824.650000000371;
        1703.49000000022	0	-799.650000000374;
        1703.49000000022	0	-774.650000000372;
        1703.49000000022	0	-749.650000000371;
        1703.49000000022	0	-724.650000000374;
        1703.49000000022	0	-699.650000000372;
        1703.49000000022	0	-674.650000000371;
        1703.49000000022	0	-649.650000000374;
        1703.49000000022	0	-624.650000000373;
        1703.49000000022	0	-599.650000000371;
        1703.49000000022	0	-574.650000000374;
        1703.49000000022	0	-549.650000000373;
        1703.49000000022	0	-524.650000000371;
        1703.49000000022	0	-499.650000000374;
        1703.49000000022	0	-474.650000000373;
        1703.49000000022	0	-449.650000000371;
        1703.49000000022	0	-424.650000000374;
        1703.49000000022	0	-399.650000000373;
        1703.49000000022	0	-374.650000000372;
        1703.49000000022	0	-349.650000000374;
        1703.49000000022	0	-324.650000000373;
        1703.49000000022	0	-299.650000000372;
        1703.49000000022	0	-274.650000000374;
        1703.49000000022	0	-249.650000000373;
        1703.49000000022	0	-224.650000000372;
        1703.49000000022	0	-199.650000000375;
        1703.49000000022	0	-174.650000000373;
        1703.49000000022	0	-149.650000000372;
        1703.49000000022	0	-124.650000000371;
        1703.49000000022	0	-99.6500000003734;
        1703.49000000022	0	-74.6500000003721;
        1703.49000000022	0	-49.6500000003707;
        1703.49000000022	0	-24.6500000003734;
        1703.49000000022	0	0;
    ]

    well_coords = [injector, lateral1, lateral2, producer];
    well_coords = [[wc[:, 1] wc[:, 2] .-wc[:, 3]] for wc in well_coords];

    well_connectivity = [
        0 0;
        1 4;
        1 4;
        0 0;
    ]

    return well_coords, well_connectivity

end

##
function get_ags_temperature_depth_profile()
    temp_vs_depth = [
        0   10;
        -200    16;
        -400    22;
        -600    28;
        -800	34;
        -1000	40;
        -1200	46;
        -1400	52;
        -1600	58;
        -1800	64;
        -2000	70;
        -2200	76;
        -2400	82;
        -2600	88;
        -2800	94;
        -3000	100;
        -3200	106;
        -3400	112;
        -3600	118;
        -3800	124;
        -4000	130;
        -4200	136;
        -4400	142;
        -4600	148;
        -4800	154;
        -5000	160;
        -5200	166;
        -5400	172;
        -5600	178;
        -5800	184;
        -6000	190;
        -6200	196;
        -6400	202;
        -6600	208;
        -6800	214;
        -7000	220;
        -7200	226;
        -7400	232;
        -7600	238;
    ]

    temp_vs_depth[:,3] .*= -1

    return temp_vs_depth
end