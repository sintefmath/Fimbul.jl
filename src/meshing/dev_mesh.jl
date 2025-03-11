import Gmsh: gmsh

gmsh.initialize()

gmsh.model.add("test")

lcar1 = .1
lcar2 = .0005
lcar3 = .055

wells = [
    [
    (0.0, -10.0, 0.0),
    (0.0, -10.0, 2000.0),
    (2000.0, -10.0, 2000.0)
    ],
    [
    (0.0, 10.0, 0.0),
    (0.0, 10.0, 2000.0),
    (2000.0, 10.0, 2000.0)
    ]
]

xw = Array{Float64}(undef, 0, 3)
for well in wells
    xw_k = reinterpret(reshape, Float64, well)
    xw = vcat(xw, xw_k)
end

xw_max = maximum(xw, dims=1)
xw_min = minimum(xw, dims=1)

pno, lno = 0, 0
for (wno, well) in enumerate(wells)

    for (k, x) in enumerate(well)
        pno += 1
        gmsh.model.geo.addPoint(x[1], x[2], x[3], lcar2, pno)
        if k > 1
            lno += 1
            gmsh.model.geo.addLine(pno - 1, pno, lno)
        end
    end
    
end


##
