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
    reached = Set{Int}()

    function dfs(node, path, on_path)
        push!(path, node)
        push!(on_path, node)
        push!(reached, node)
        if isempty(adj[node])
            push!(branches, copy(path))
        else
            for child in adj[node]
                if child in on_path
                    # Include the closing edge, then stop this branch.
                    push!(branches, [path; child])
                else
                    dfs(child, path, on_path)
                end
            end
        end
        delete!(on_path, node)
        pop!(path)
    end

    for root in roots
        dfs(root, Int[], Set{Int}())
    end

    # A component made entirely of cycles has no root.
    for node in nodes
        if !(node in reached)
            dfs(node, Int[], Set{Int}())
        end
    end

    return branches
end
