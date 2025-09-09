# Add required modules to namespace 
using Jutul, JutulDarcy, Fimbul
using HYPRE
using Statistics
using GLMakie
# Useful SI units

function topo_sort_well(cells, msh, geo = missing)

    sorted_cells = Int[]
    N = get_neighborship(UnstructuredMesh(msh))

    geo = ismissing(geo) ? tpfv_geometry(msh) : geo
    z = geo.cell_centroids[3, cells]

    top = last(findmin(z))
    println(top)
    push!(sorted_cells, popat!(cells, top))

    current_cell = sorted_cells[end]
    while !isempty(cells)
        for r = 1:2
            neighbors = N[r,:] .== current_cell
            for n in N[mod(r,2)+1, neighbors]
                if n in cells
                    push!(sorted_cells, n)
                    popat!(cells, findfirst(isequal(n), cells))
                    current_cell = n
                    continue
                end
            end
        end
    end

    return sorted_cells

end

##
Kelvin, joule, watt = si_units(:Kelvin, :joule, :watt)
kilogram = si_unit(:kilogram)
meter = si_unit(:meter)
darcy = si_unit(:darcy);

well_spacing = 200.0 .* meter
fracture_radius = 200.0 .* meter
well_depth = 3000.0 .* meter
well_lateral = 2000.0 .* meter
fracture_aperture = 1e-3 .* meter
ws, wd, wl = well_spacing, well_depth, well_lateral
# well_coords = [
#     [(-ws/2.0, 0.0, 0.0),(-ws/2.0, 0.0, wd), (-ws/2.0, wl, wd)],
#     [(ws/2.0, 0.0, 0.0), (ws/2.0, 0.0, wd), (ws/2.0, wl, wd)]
# ]

well_coords = [
    [-ws/2 0.0  0.0; -ws/2 0.0 wd; -ws/2 wl wd], 
    [ ws/2 0.0  0.0;  ws/2 0.0 wd;  ws/2 wl wd]
]

out = Fimbul.egs(well_coords, fracture_radius, 100.0;
);

##
domain, is_frac, msh = out


mm = UnstructuredMesh(msh)
N = get_neighborship(mm)
# for well in well_coords
well = well_coords[1]
    cells = Jutul.find_enclosing_cells(msh, well, n = 1000)
    println(length(cells))
    println(sum(is_frac[cells]))
    all_cells = Int[]
    for cell in cells
        faces = mm.faces.cells_to_faces[cell]
        c = vcat(N[:, faces]...)
        keep = is_frac[c]
        push!(all_cells, c[keep]...)
    end
    push!(all_cells, cells...)
    all_cells = unique!(all_cells)

# end
##

depths = [0.0, 2000.0, 2750.0, 3500.0] .* meter
well_distance = 200.0 .* meter
heel_yz = [0.0, 2250.0] .* meter
toe_yz = [2000.0, 2350.0] .* meter
wx_inj = [(-well_distance/2.0, 0.0)]
xw_prod = [(well_distance/2.0, 0.0)]

##

fracture_radius = 200meter
hxy_min = missing
hxy_max = missing
hz_min = missing
hz_max = missing
offset_rel = 2.0
# Calculate domain offset to minimize boundary effects
offset = offset_rel*(well_distance + fracture_radius)
# Well coordinates
wd = well_distance
ld = 2000meter
hd = 3000meter

well_coords = [
    [(-wd/2.0, 0.0, 0.0), (-wd/2.0, 0.0, hd), (-wd/2.0, ld, hd)],
    [(wd/2.0, 0.0, 0.0), (wd/2.0, 0.0, hd), (wd/2.0, ld, hd)],
]

# Define layer depths
depths = [0.0, 2000.0, 2750.0, 4000.0] .* meter
fracture_radius = 200meter

hxz_min = well_distance/10
hxz_max = offset/5
# set x coordinates based on well x coordinates

##


function make_coords_xz(x0, xw, h_min, h_max)

    x_min = minimum(xw)
    x_max = maximum(xw)
    x_mid = (x_min + x_max)/2

    xw .-= hxz_min/2

    xw = vcat(x_mid - fracture_radius, xw, x_mid + fracture_radius)
    
    x = sort(vcat(x0, xw))
    nx = length(x)-1

    hxz = fill(h_max, nx)
    x_mid = diff(x)./2 .+ x[1:end-1]
    println(x_mid)
    is_frac = x_min - fracture_radius .<= x_mid .<= x_max + fracture_radius
    println(is_frac)
    hxz[is_frac] .= h_min

    x = first(Fimbul.interpolate_z(x, hxz))

    return x

end

##
x0 = [
    -800.0,
    800.0
]
xw = []
for wc in well_coords
    x = unique!([w[1] for w in wc])
    push!(xw, x...)
end
x = make_coords_xz(x0, xw, 10.0, 200.0)

z0 = [
    0.0,
    4000.0
]
zw = []
for wc in well_coords
    z = wc[2][3]
    push!(zw, z...)
end
z = make_coords_xz(z0, zw, 10.0, 500.0)

yw = []
for wc in well_coords
    y = unique!([w[2] for w in wc])
    push!(yw, y...)
end
unique!(yw)

y_min, y_max = extrema(yw)

fracture_start = 50.0
fracture_end = ld - 50.0
fracture_spacing = 100.0
aperture = 1e-3

fracture_y = fracture_start:fracture_spacing:fracture_end
fracture_y = sort(vcat(fracture_y, fracture_y .+ aperture))

offset_y = 200
y = [
    y_min - offset_y,
    y_min,
    fracture_y...,
    y_max,
    y_max + offset_y
]

hy_min = fracture_spacing/3
hy_max = 200meter
hy = fill(hy_max, length(y)-1)
y_mid = diff(y)./2 .+ y[1:end-1]
is_frac = y_min .< y_mid .< y_max
hy[is_frac] .= hy_min

y = first(Fimbul.interpolate_z(y, hy))

xyz = (x, y, z)
sizes = map(x->diff(x), xyz)
dims = Tuple([length(s) for s in sizes])
msh = CartesianMesh(dims, sizes, origin = minimum.(xyz))

