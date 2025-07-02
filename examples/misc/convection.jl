using Jutul, JutulDarcy, Fimbul
using HYPRE

##

L = 1000
circle = (x0, r, n) -> [r.*(cos(θ), sin(θ)) .+ x0 for θ in range(π/2, 2π+π/2; length = n+1)]

x0 = (0.0, 0.0)
xb = circle((0.0, 0.0), L/3, 3)[1:end-1]

cc = []
for x0 in xb
    cc_k = circle(x0, L/5, 20)
    push!(cc, cc_k)
end

boundary = [(-L, -L), (L, -L), (L, L), (-L, L)]

hxy = L/50
msh, _ = Fimbul.extruded_mesh(cc, [0.0, 1.0]; boundary = boundary, hxy_min = hxy, hxy_max = hxy, hz = 1.0, recombine_to_quads = true)
# np = [[x[1], x[3], x[2]] for x in msh.node_points]

# g = UnstructuredMesh(
#         msh.structure,
#         msh.faces,
#         msh.boundary_faces,
#         np,
#         msh.cell_map,
#         msh.face_map,
#         msh.boundary_map,
#         msh.node_map,
#         msh.tags,
#         msh.z_is_depth
#     )

fig = Figure()
ax = Axis3(fig[1,1])
plot_mesh_edges!(ax, msh)
fig

##
domain = reservoir_domain(msh)
model, parameters = setup_reservoir_model(domain, :geothermal);

##
Kelvin = si_unit(:Kelvin)
meter = si_unit(:meter)
year = si_unit(:year)

geo = tpfv_geometry(msh)

xy = geo.cell_centroids[1:2,:]
xy = [(x[1], x[2]) for x in eachcol(xy)]

surface_temperature = convert_to_si(10.0, :Celsius)
surface_pressure = convert_to_si(5.0, :atm)

dTdz = 0.03Kelvin/meter
rho = reservoir_model(model).system.rho_ref[1]
dpdz = rho*gravity_constant
y = geo.cell_centroids[2,:]
z = maximum(y) .- y
T = surface_temperature .+ dTdz.*z
p = surface_pressure .+ dpdz.*z
T_blobs = convert_to_si.([10.0, 95.0, 95.0], :Celsius)

bc = []

for (bno, xc) in enumerate(cc)
    inside = map(p -> Fimbul.point_in_polygon(p, xc), xy)
    cells = findall(inside)
    T[cells] .= T_blobs[bno]

    # faces = [xor(inside[n[1]], inside[n[2]]) for n in eachcol(geo.neighbors)]
    # faces = findall(faces)
    # cells = []
    # for face in faces
    #     c1, c2 = geo.neighbors[:, face]
    #     if inside[c1]
    #         push!(cells, c1)
    #     else
    #         push!(cells, c2)
    #     end
    # end

    # JutulDarcy.flow_boundary_condition!(bc, domain, cells, p[cells], T_blobs[bno])

end

##
plot_cell_data(msh, T)

##
state0 = setup_reservoir_state(model; Pressure = p, Temperature = T)

forces = setup_reservoir_forces(model; bc = bc)
dt = fill(year, 100)

case = JutulCase(model, dt, forces; state0 = state0)

##
res = simulate_reservoir(case; info_level = 2)

##
Tf = res.states[end][:Temperature]
xy = geo.cell_centroids[1:2,:]

colors = cgrad(:seaborn_icefire_gradient, 8, categorical=true)
# colors = cgrad(:hot, 8, categorical=true)
plot_cell_data(msh, Tf; colormap = colors, shading = NoShading)

##
plot_reservoir(case, res.states; colormap = colors)