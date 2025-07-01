function Fimbul.plot_well_data(time, reference, other; 
    wells = :all, 
    names = missing, 
    field = :energy, 
    nan_ix = missing,
    separate_legend = false,
    kwargs...)

    if wells == :all
        wells = setdiff(keys(reference[1]), [:Reservoir, :Facility])
    end

    fig = Figure(size = (800, 400), fontsize = 20)
    ax = Axis(fig[1, 1], xlabel = "Time (years)", ylabel = "Temperature (°C)")

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

    if separate_legend
        fig_legend = Figure(size = (800, 400), fontsize = 20)
        Legend(fig_legend[1, 1], ax, loc = :best)
        return fig, fig_legend
    else
        Legend(fig[1, 2], ax, loc = :best)
        return fig
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

function Fimbul.plot_mswell_values!(ax, model, well, values; label = nothing, geo = missing, kwargs...)

    if ismissing(geo)
        rmodel = reservoir_model(model)
        msh = physical_representation(rmodel.data_domain)
        geo = tpfv_geometry(msh)
    end
        
    well_representation = model.models[well].domain.representation
    
    N = well_representation.neighborship

    branches = dfs_branches(N)
    rcells = well_representation.perforations.reservoir
    
    z = vcat(well_representation.top.reference_depth, geo.cell_centroids[3, rcells])
    # Plot edges
    for (i, branch) in enumerate(branches)

        zb = z[branch]
        lbl = (i == 1) ? label : nothing
        if size(values, 2) == 1
            vb = values[branch]
            lines!(ax, vb, zb; color = :blue, label = lbl, kwargs...)
        else
            vb = values[branch,:]
            lines!(ax, vb[:, 1], vb[:,2], z_i; color = :blue, label = lbl, kwargs...)
        end
    end
    
end

function dfs_branches(edges)
    # Build adjacency list
    from_nodes = edges[1,:]
    to_nodes = edges[2,:]
    nodes = unique(vcat(from_nodes, to_nodes))
    adj = Dict(n => Int[] for n in nodes)
    for (f, t) in zip(from_nodes, to_nodes)
        push!(adj[f], t)
    end

    # Find root nodes (nodes that never appear as a "to" node)
    roots = setdiff(from_nodes, to_nodes)
    branches = Vector{Vector{Int}}()

    function dfs(node, path)
        push!(path, node)
        if isempty(adj[node])
            push!(branches, copy(path))
        else
            for child in adj[node]
                dfs(child, path)
            end
        end
        pop!(path)
    end

    for root in roots
        dfs(root, Int[])
    end

    return branches
end