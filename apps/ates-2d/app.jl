module App

    using GenieFramework
    using JutulDarcy
    using Jutul
    using HYPRE
    using DataFrames
    using Base64
    import CairoMakie

    @genietools

    DEFAULT_CYCLES = 3
    DEFAULT_DELTA_HOT = 50.0
    DEFAULT_DELTA_COLD = 30.0
    DEFAULT_BASE_TEMP = 40.0
    DEFAULT_CHARGE_PERIOD = 2:5
    DEFAULT_DISCHARGE_PERIOD = 8:12

    function simulate_hates(;
            t_res = DEFAULT_BASE_TEMP,
            delta_hot = DEFAULT_DELTA_HOT,
            delta_cold = DEFAULT_DELTA_COLD,
            charge_period = DEFAULT_CHARGE_PERIOD,
            discharge_period = DEFAULT_DISCHARGE_PERIOD,
            ncycles = DEFAULT_CYCLES
        )
        darcy, litre, year, second = si_units(:darcy, :litre, :year, :second)
        
        nx = 100
        nz = 100

        nx = 20
        ny = 20
        
        temperature_top = convert_to_si(40.0, :Celsius)
        pressure_top = convert_to_si(120.0, :bar)
        delta_charge = delta_hot
        delta_discharge = delta_cold
        
        grad_p = 1000*9.81
        grad_T = 0.3
        
        # ## Set up the reservoir
        g = CartesianMesh((nx, 1, nz), (250.0, 250.0, 75.0))
        reservoir = reservoir_domain(g,
            permeability = [0.3, 0.3, 0.1].*darcy,
            porosity = 0.3,
            rock_thermal_conductivity = 2.0,
            fluid_thermal_conductivity = 0.6
        )
        
        depth = reservoir[:cell_centroids][3, :]
        # ## Define wells and model
        di = Int(ceil(nx/4))
        k = Int(ceil(nz/2))
        Whot = setup_vertical_well(reservoir, 0+di   , 1, toe = k, name = :Hot)
        Wcold = setup_vertical_well(reservoir, nx-di+1, 1, toe = k, name = :Cold)
        
        model, parameters = setup_reservoir_model(reservoir, :geothermal, wells = [Whot, Wcold]);
        # ## Set up boundary and initial conditions
        bcells = Int[]
        pressure_res = Float64[]
        temperature_res = Float64[]
        for cell in 1:number_of_cells(g)
            d = depth[cell]
            push!(pressure_res, pressure_top + grad_p*d)
            push!(temperature_res, temperature_top + grad_T*d)
        
            I, J, K = cell_ijk(g, cell)
            if I == 1 || I == nx
                push!(bcells, cell)
            end
        end
        
        bc = flow_boundary_condition(bcells, reservoir, pressure_res[bcells], temperature_res[bcells])
        # ## Set up the schedule
        
        # ### Set up forces
        
        # We assume we have a supply amounting to 90 C at 25 l/s for storage. During the
        # rest period, we assume the same discharge rate and a temperature of 10 C.
        charge_rate = 25litre/second
        discharge_rate = charge_rate
        temperature_charge = temperature_top + delta_charge
        temperature_discharge = temperature_top - delta_discharge
        
        # Set up forces for charging
        rate_target = TotalRateTarget(charge_rate)
        ctrl_hot  = InjectorControl(rate_target, [1.0], density = 1000.0, temperature = temperature_charge)
        rate_target = TotalRateTarget(-charge_rate)
        ctrl_cold = ProducerControl(rate_target)
        forces_charge = setup_reservoir_forces(model, control = Dict(:Hot => ctrl_hot, :Cold => ctrl_cold), bc = bc)
        
        # Set up forces for discharging
        rate_target = TotalRateTarget(discharge_rate)
        ctrl_cold = InjectorControl(rate_target, [1.0], density = 1000.0, temperature = temperature_discharge)
        rate_target = TotalRateTarget(-discharge_rate)
        ctrl_hot = ProducerControl(rate_target)
        forces_discharge = setup_reservoir_forces(model, control = Dict(:Hot => ctrl_hot, :Cold => ctrl_cold), bc = bc)
        
        # ### Set up forces for rest period
        forces_rest = setup_reservoir_forces(model, bc = bc)
        
        # ### Set up timesteps and assign forces to each timestep
        dt = Float64[]
        forces = []
        month = year/12
        num_charge = 0
        num_discharge = 0
        num_rest = 0
        for year in 1:ncycles
            for mno in vcat(6:12, 1:5)
                if mno in charge_period
                    push!(dt, month)
                    push!(forces, forces_charge)
                    num_charge += 1
                elseif mno in discharge_period
                    push!(dt, month)
                    push!(forces, forces_discharge)
                    num_discharge += 1
                else
                    push!(dt, month)
                    push!(forces, forces_rest)
                    num_rest += 1
                end
            end
        end
        @info "Set up schedule" num_charge num_discharge num_rest
        # ## Set up initial state
        state0 = setup_reservoir_state(model, Pressure = pressure_res, Temperature = temperature_res)
        # ## Simulate the case
        ws, states, time_s = simulate_reservoir(state0, model, dt,
            forces = forces,
            parameters = parameters,
            info_level = 1
        )
        @info "Simulation done"
        # ## Plot energy recovery factor
        # The energy recovery factor η is defined as the amount of stored to produced
        # energy. We plot this both cumulatively and for each of the 25 yearly cycles
        wd = ws.wells[:Hot]
        c_p_water = 4.186 # kJ/kgK
        well_temp = wd[:temperature]
        well_temp_cold = ws.wells[:Cold][:temperature]

        q = wd[:mass_rate]
        storage = q .> 0
        q_store = q.*storage
        q_prod = q.*(.!storage)
        stored_energy = well_temp.*q_store.*c_p_water.*dt
        produced_energy = -well_temp.*q_prod.*c_p_water.*dt
        η_cumulative = cumsum(produced_energy)./cumsum(stored_energy)
        t = cumsum(dt)./si_unit(:day)

        @info "Cycles"
        num_years = ncycles
        eta, T = zeros(num_years), zeros(num_years)
        for i = 1:num_years
            ix = (1:12) .+ 12*(i-1)
            se = sum(stored_energy[ix])
            pe = sum(produced_energy[ix])
            eta[i] = pe/se
            T[i] = t[ix[end]]
        end
        @info "Cycles done."

        return Dict(
            :wtemp => well_temp .- 273.15,
            :ctemp => well_temp_cold .- 273.15,
            :time => time_s./si_unit(:day),
            :T => T,
            :eta => eta,
            :T_spatial => map(x -> reshape(x[:Temperature] .- 273.15, nx, nz), states)
        )
    end

    function figure_to_html(fig)
        # buffer = Base.IOBuffer()
        fname = "tmpfile.png"
        CairoMakie.save(fname, fig)
        buffer = open(fname, "r")
        data = base64encode(buffer)
        close(buffer)
        rm(fname)
        # return html("""<img src="data:image/png;base64,$(data)">""")
        return "data:image/png;base64,$(data)"
    end

    @app begin
        @in ncycles = DEFAULT_CYCLES
        @in delta_hot = DEFAULT_DELTA_HOT
        @in delta_cold = DEFAULT_DELTA_COLD
        @in base_temp = DEFAULT_BASE_TEMP
        @in name = "Genie"
        @in start = false
        @in running = false
        @in ButtonProgress_process = false
        @in ButtonProgress_progress = 0.0
        @in ChargePeriod = RangeData(DEFAULT_CHARGE_PERIOD)
        @in DisChargePeriod = RangeData(DEFAULT_DISCHARGE_PERIOD)
        @in tab_selected = "hot_temp"
        @in tplot_stepno = 0.5
        @out hotplot = PlotData()
        @out coldplot = PlotData()
        @out etaplot = PlotData()
        @out sim_result = simulate_hates()
        @out imgstr = ""
        @private u_x = []
        @private u_y = []
        @onchange tplot_stepno begin
            temperature_spatial = sim_result[:T_spatial]
            nstep = length(temperature_spatial)
            ix = clamp(Int(round(tplot_stepno*nstep)), 1, nstep)
            data = temperature_spatial[ix]
            data = data[:, end:-1:1]
            fig = CairoMakie.Figure(size = (800, 300))
            ax = CairoMakie.Axis(fig[1, 1], title = "Step $ix/$nstep")
            plt = CairoMakie.heatmap!(ax, data, colormap = :hot, colorrange = (base_temp - delta_cold, base_temp + delta_hot))
            CairoMakie.Colorbar(fig[1, 2], plt, label = "Temperature / °C")
            imgstr = figure_to_html(fig)
        end
        @onchange ChargePeriod begin
            # Should truncate the discharge period here.
        end
        @onbutton ButtonProgress_process begin
            @info "Hello button clicked" running
            running = false
            # u_x = []
            # u_y = []
            empty!(u_x)
            empty!(u_y)
            t = 0.0
            if running == false
                running = true
                res = simulate_hates(
                    t_res = base_temp,
                    delta_hot = delta_hot,
                    delta_cold = delta_cold,
                    charge_period = ChargePeriod.range,
                    discharge_period = DisChargePeriod.range,
                    ncycles = ncycles
                )
                @info "Simulation ok? "
                ButtonProgress_progress = 0.0
                wtemp = res[:wtemp]
                t = res[:time]
                @info "" res

                for i in eachindex(wtemp, t)
                    push!(u_x, t[i])
                    push!(u_y, wtemp[i])
                end
                @show ButtonProgress_progress
                hotplot = PlotData(
                    x = u_x,
                    y = u_y,
                    plot = StipplePlotly.Charts.PLOT_TYPE_LINE
                )
                coldplot = PlotData(
                    x = u_x,
                    y = res[:ctemp],
                    plot = StipplePlotly.Charts.PLOT_TYPE_LINE
                )
                yr = res[:T]
                eta = res[:eta]
                etaplot = PlotData(
                    x = yr,
                    y = eta,
                    plot = StipplePlotly.Charts.PLOT_TYPE_LINE
                )
            end
        end
    end


    function ui()
    [
        h1("High-temperature aquifer thermal energy storage (HT-ATES)")
        p("Fast simulation of energy storage with Fimbul+JutulDarcy.jl - View the results in the tabs below")
        [
            tabgroup(
                :tab_selected,
                inlinelabel = true,
                class = "bg-primary text-white shadow-2",
                [
                    tab(name = "hot_temp", icon = "local_fire_department", label = "Hot well"),
                    tab(name = "cold_temp", icon = "bolt", label = "Cold well"),
                    tab(name = "energy", icon = "bolt", label = "Energy recovered"),
                    tab(name = "reservoir", icon = "volcano", label = "Reservoir"),
                ],
            ),
            tabpanels(
                :tab_selected,
                animated = true,
                var"transition-prev" = "scale",
                var"transition-next" = "scale",
                [
                    tabpanel(name = "hot_temp", [
                        plot(:hotplot, layout = PlotLayout(
                                xaxis = [
                                        PlotLayoutAxis(xy="x", title="Time / days")
                                    ],
                                yaxis = [
                                        PlotLayoutAxis(xy="y", title="Temperature / °C")
                                    ],
                            )
                        )
                    ]),
                    tabpanel(name = "cold_temp", [
                        plot(:coldplot, layout = PlotLayout(
                            xaxis = [
                                    PlotLayoutAxis(xy="x", title="Time / years")
                                ],
                            yaxis = [
                                    PlotLayoutAxis(xy="y", title="Temperature / °C")
                                ],
                        )
                    )
                    ]),
                    tabpanel(name = "energy", [
                        plot(:etaplot, layout = PlotLayout(
                            xaxis = [
                                    PlotLayoutAxis(xy="x", title="Time / years")
                                ],
                            yaxis = [
                                    PlotLayoutAxis(xy="y", title="Produced energy / kJ")
                                ],
                        )
                    )
                    ]),
                    tabpanel(name = "reservoir", [
                        p("Reservoir temperature"),
                        # html("{{imgstr}}"),
                        imageview(src = :imgstr),
                        itemsection(slider(0.0:0.001:1.0, :tplot_stepno, label = "", color = "red")),
                    ]),
                ],
            ),
        ]

        p("Charging and discharging temperature difference")
        item([
            itemsection(avatar = "", icon("local_fire_department", color = "red")),
            itemsection(slider(0:1:50, :delta_hot, label = "", color = "red")),
            itemsection(avatar = "", icon("ac_unit", color = "blue")),
            itemsection(slider(0:1:50, :delta_cold, label = "", color = "blue")),
        ])
        p("Number of yearly cycles to simulate")
        item(
            [
                itemsection(avatar = "", icon("keyboard_double_arrow_up", color = "teal")),
                itemsection(slider(1:1:50, :ncycles, label = "", color = "teal"))
            ]
        )
        p("Simulating {{ncycles}} cycles, charge ΔT of {{delta_hot}}°C and discharge ΔT of {{delta_cold}}°C", class="st-module")
        p("Charging period")
        item(
            [
                itemsection(avatar = "", icon("keyboard_double_arrow_down", color = "red")),
                range(1:1:12, :ChargePeriod, markers = true, label = true, color = "red"),
            ]
        )
        p("Discharge period")
        item(
            [
                itemsection(avatar = "", icon("keyboard_double_arrow_up", color = "blue")),
                range(1:1:12, :DisChargePeriod, markers = true, label = true, color = "blue")
            ]
        )

        btn(
            "Run simulation",
            @click(:ButtonProgress_process),
            loading = :ButtonProgress_process,
            percentage = :ButtonProgress_progress,
            color = "primary",
            icon = "rocket_launch",
            class = "q-mr-sm",
        )
        separator()
    ]
    end

    @page("/", ui)
end
