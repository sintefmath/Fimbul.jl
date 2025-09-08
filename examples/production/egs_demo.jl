# Add required modules to namespace 
using Jutul, JutulDarcy, Fimbul
using HYPRE
using Statistics
using GLMakie
# Useful SI units
Kelvin, joule, watt = si_units(:Kelvin, :joule, :watt)
kilogram = si_unit(:kilogram)
meter = si_unit(:meter)
darcy = si_unit(:darcy);

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

hy_min = fracture_spacing/5
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
