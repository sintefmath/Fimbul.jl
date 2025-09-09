# function egs(well_spacing, well_depth, well_lateral;
#     fracture_radius=well_spacing*1.25,
#     fracture_spacing=well_lateral/10,
#     kwargs...
# )

#     ws, wd, wl = well_spacing, well_depth, well_lateral
#     well_coords = [
#         [(-ws/2.0, 0.0, 0.0),(-ws/2.0, 0.0, wd), (-ws/2.0, wl, wd)],
#         [(ws/2.0, 0.0, 0.0), (ws/2.0, 0.0, wd), (ws/2.0, wl, wd)],
#     ]

#     return egs(well_coords, fracture_radius, fracture_spacing; kwargs...)

# end

function egs(well_coords, fracture_radius, fracture_spacing;
    well_names = missing,
    fracture_aperture=1e-3,
    porosity = 0.1,
    permeability = 1.0 .* darcy,
    rock_thermal_conductivity = 2.5 .* watt/(meter*Kelvin),
    rock_heat_capacity = 900.0*joule/(kilogram*Kelvin),
    temperature_inj = convert_to_si(25, :Celsius),
    msh_args = NamedTuple()
)

    msh = Fimbul.make_egs_cart_mesh(
        well_coords, fracture_radius, fracture_spacing;
        msh_args...
    )

    Δy = [cell_dims(msh, c)[2] for c in 1:number_of_cells(msh)]

    is_frac = isapprox.(Δy, fracture_aperture)
    geo = tpfv_geometry(msh)
    xz = geo.cell_centroids[[1,3],:]
    xc = well_coords[1][2, [1,3]]
    println(xc)
    r = vec(sum((xc .- xz).^2, dims=1).^0.5)
    is_frac = is_frac .&& (r .<= fracture_radius)

    nc = number_of_cells(msh)
    function process_prop(v, v_frac)
        v = fill(v, nc)
        v[is_frac] .= v_frac
        return v
    end

    permeability_frac = fracture_aperture^2/12
    porosity = process_prop(porosity, 0.5)
    permeability = process_prop(permeability, permeability_frac)
    domain = reservoir_domain(msh;
        porosity = porosity,
        permeability = permeability,
        rock_thermal_conductivity = rock_thermal_conductivity,
        rock_heat_capacity = rock_heat_capacity
    )

    if ismissing(well_names)
        well_names = ["inj"]
        for i in 2:length(well_coords)
            push!(well_names, "Prod$(i-1)")
        end
    end

    return domain, is_frac, msh

    for well in well_coords
        cells = Jutul.find_enclosing_cells(msh, well)
        println(length(cells))
        println(sum(is_frac[cells]))
    end

    return domain

end

function make_egs_cart_mesh(well_coords, fracture_radius, fracture_spacing;
    hxz_min = missing,
    hxz_max = missing,
    hy_min = missing,
    hy_max = missing,
    offset_rel = 2.0,
)
    # Calculate domain offset to minimize boundary effects
    well_distance = abs(well_coords[1][1][1] - well_coords[2][1][1])
    y_min = minimum(well_coords[1][:, 2])
    y_max = maximum(well_coords[1][:, 2])
    well_lateral = abs(y_max - y_min)

    offset = offset_rel*(well_distance + fracture_radius)
    # Well coordinates
    hxz_min = ismissing(hxz_min) ? well_distance/10 : hxz_min
    hxz_max = ismissing(hxz_max) ? offset/5 : hxz_max

    hy_min = ismissing(hy_min) ? fracture_spacing/3 : hy_min
    hy_max = ismissing(hy_max) ? offset/5 : hy_max

    function make_coords_xz(x0, xw)

        x_min, x_max = extrema(xw)
        x_mid = (x_min + x_max)/2

        xw .-= hxz_min/2
        xw = vcat(x_mid - fracture_radius, xw, x_mid + fracture_radius)
        xd = sort(vcat(x0, xw))
        nx = length(xd)-1
        hxz = fill(hxz_max, nx)
        x_mid = diff(xd)./2 .+ xd[1:end-1]
        is_frac = x_min - fracture_radius .<= x_mid .<= x_max + fracture_radius
        hxz[is_frac] .= hxz_min

        xd = first(Fimbul.interpolate_z(xd, hxz))

        return xd

    end

    xw = []
    for wc in well_coords
        x = unique!(wc[:, 1])
        push!(xw, x...)
    end
    x_min, x_max = extrema(xw)
    x_mid = (x_min + x_max)/2
    x0 = [
        x_mid - offset,
        x_mid + offset
    ]
    x = make_coords_xz(x0, xw)
    println(x)

    zw = []
    for wc in well_coords
        z = wc[2,3]
        push!(zw, z...)
    end
    z0 = [
        0.0,
        4000.0
    ]
    z = make_coords_xz(z0, zw)

    yw = []
    for wc in well_coords
        y = unique!(wc[:, 2])
        push!(yw, y...)
    end
    unique!(yw)

    y_min, y_max = extrema(yw)

    fracture_start = y_min + fracture_spacing/2
    fracture_end = y_max - fracture_spacing/2
    aperture = 1e-3

    fracture_y = fracture_start:fracture_spacing:fracture_end
    fracture_y = sort(vcat(fracture_y, fracture_y .+ aperture))

    y = [
        y_min - offset,
        y_min - hy_min/2,
        fracture_y...,
        y_max - hy_min/2,
        y_max + offset
    ]

    println(y)

    hy = fill(hy_max, length(y)-1)
    y_mid = diff(y)./2 .+ y[1:end-1]
    is_frac = y_min .< y_mid .< y_max
    hy[is_frac] .= hy_min

    y = first(Fimbul.interpolate_z(y, hy))

    xyz = (x, y, z)
    sizes = map(x->diff(x), xyz)
    dims = Tuple([length(s) for s in sizes])
    msh = CartesianMesh(dims, sizes, origin = minimum.(xyz))

    return msh

end