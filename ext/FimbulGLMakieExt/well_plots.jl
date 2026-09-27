function Fimbul.plot_well_data!(ax, time, reference, other; 
    wells = :all, 
    names = missing, 
    field = :energy, 
    nan_ix = missing,
    legend = true,
    kwargs...)

    if wells == :all
        wells = setdiff(keys(reference[1]), [:Reservoir, :Facility])
    end

    for well in wells
        vr = get_field(reference, well, field)
        lines!(ax, time, vr; label=names[1], linestyle = (:dash, 1), linewidth = 6, color = :black)
        for (i, proxy) in enumerate(other)
            vp = get_field(proxy, well, field)
            if !ismissing(nan_ix) && !ismissing(nan_ix[i])
                vp[nan_ix[i]] .= NaN
            end
            lines!(ax, time, vp, label=names[i+1], 
                linewidth = 3, kwargs...)
        end
    end

    if legend
        axislegend(ax, loc = :best)
    end
    
end
function get_field(data, well, field)

    getter = (data, field) -> [d[well][field][1] for d in data]
    if field == :Energy
        T0 = convert_to_si(20.0, :Celsius)
        Cp = 4.1864
        p = 1.0si_unit(:atm)
        ρ = 1000.0
        h0 = Cp*T0 + p/ρ
        q = getter(data, :TotalMassFlux)
        h = getter(data, :FluidEnthalpy) .- h0
        v = abs.(q.*h)
    else
        v = getter(data, field)
    end

    if field == :Temperature
        v = convert_from_si.(v, :Celsius)
    end

    return v

end

function Fimbul.plot_mswell_values!(ax, model, well, values; nodes = missing, label = nothing, geo = missing, kwargs...)

    if ismissing(geo)
        rmodel = reservoir_model(model)
        msh = physical_representation(rmodel.data_domain)
        geo = tpfv_geometry(msh)
    end
        
    well_representation = model.models[well].domain.representation
    
    N = well_representation.neighborship
    if !ismissing(nodes)
        keep = [n[1] ∈ nodes && n[2] ∈ nodes for n in eachcol(N)]
        N = N[:, keep]
    end

    branches = dfs_branches(N)
    rcells = well_representation.perforations.reservoir
    
    z = model.models[well].data_domain[:cell_centroids][3, :]
    # Plot edges
    for (i, branch) in enumerate(branches)
        zb = z[branch]
        lbl = (i == 1) ? label : nothing
        if size(values, 2) == 1
            vb = values[branch]
            lines!(ax, vb, zb; color = :blue, label = lbl, kwargs...)
        else
            vb = values[branch,:]
            lines!(ax, vb[:, 1], vb[:,2], zb; color = :blue, label = lbl, kwargs...)
        end
    end
    
end

