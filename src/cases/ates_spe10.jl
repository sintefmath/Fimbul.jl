const SPE10_FORMATIONS = Dict(:tarbert => 1:35, :upper_ness => 36:85)
const SPE10_DATA = Ref{Any}(nothing)

"""
    spe10_data()

Permeabilities `kx`, `ky`, `kz` [m²] and porosity `poro` of the full SPE10
model 2 (Christie & Blunt, 2001) as 60 × 220 × 85 arrays (layer 1 at the top),
from `JutulDarcy.SPE10.load_data`, read on first use and kept.
"""
function spe10_data()
    isnothing(SPE10_DATA[]) && (SPE10_DATA[] = JutulDarcy.SPE10.load_data())
    return SPE10_DATA[]
end

"""
    spe10_block(layers; ix = 1:60, iy = 1:220, coarsen = (1, 1, 1), minporo = 0.01)

The SPE10 model 2 for the given `layers` (and optionally a lateral subset
`ix × iy` of its columns): `kx`, `ky`, `kz` [m²] and `poro` as `nx × ny × nl`
arrays (layer 1 at the top), and the cell sizes `dx`, `dy`, `dz` [m]. With
`coarsen = (cx, cy, cz)`, blocks of `cx × cy × cz` cells are merged:
arithmetic mean of the porosity and, for each permeability direction, the
harmonic mean along that direction of the arithmetic means across it.
"""
function spe10_block(layers; ix = 1:60, iy = 1:220, coarsen = (1, 1, 1), minporo = 0.01)
    d = spe10_data()
    kx, ky, kz = d.kx[ix, iy, layers], d.ky[ix, iy, layers], d.kz[ix, iy, layers]
    phi = max.(d.poro[ix, iy, layers], minporo)
    n = size(phi)
    all(n .% coarsen .== 0) || throw(ArgumentError("coarsen = $coarsen must divide the block size $n"))
    m = n .÷ coarsen
    blocks(I) = ntuple(a -> ((I[a] - 1) * coarsen[a] + 1):(I[a] * coarsen[a]), 3)
    harmonic(x) = length(x) / sum(1 ./ max.(x, 1.0e-11 * darcy))
    function upscale(k, dir)
        out = zeros(m)
        for I in CartesianIndices(m)
            b = k[blocks(Tuple(I))...]
            out[I] = harmonic(vec(mean(b, dims = filter(!=(dir), (1, 2, 3)))))
        end
        return out
    end
    poro = [mean(phi[blocks(Tuple(I))...]) for I in CartesianIndices(m)]
    h = JutulDarcy.SPE10.CELL_DIMS .* coarsen
    return (
        kx = upscale(kx, 1), ky = upscale(ky, 2), kz = upscale(kz, 3), poro = poro,
        dx = h[1], dy = h[2], dz = h[3],
    )
end

"""
    ates_spe10(; <keyword arguments>)

Aquifer thermal energy storage (ATES) in the SPE10 model 2 (Christie & Blunt,
2001): a hot and a cold well in the Upper Ness (fluvial channels, layers
36–85) or Tarbert (shallow marine, layers 1–35) formation, hot water injected
in the charge period and produced in the discharge period, with balanced
rates (the injector reinjects what the other well produces).

The SPE10 block (60 × 220 cells of 20 × 10 ft and 2 ft thick layers, or a
lateral subset `ix × iy`, optionally coarsened) is the aquifer, in the 3D
model or, with `slice`, a vertical slice along the wells. Low-permeability cap
and base rock with layers growing away from the aquifer carry heat by
conduction, and padding cells with the layer-averaged properties around the
block take the hydrostatic / geothermal boundary conditions away from the
wells. The wells are on the long (y) axis, at `well_positions` of the block
length, perforated over the full aquifer thickness.

# Keyword arguments

## Geometry
- `formation = :upper_ness`: `:upper_ness` or `:tarbert`
- `layers`: SPE10 layers of the aquifer (default those of `formation`)
- `ix = 1:60`, `iy = 1:220`: SPE10 columns of the aquifer block
- `slice = missing`: column `i` of a vertical slice along the wells (sets
  `ix = i:i` without padding in x, and no boundary conditions on the sides of
  the slice)
- `coarsen = (2, 2, 2)`: SPE10 cells merged in each block cell (see
  `spe10_block`)
- `depth_top = 800.0`: depth of the top of the aquifer [m]
- `cap`, `base`: thicknesses of the cap and base rock layers, from the
  aquifer outwards [m]
- `n_pad = 4`, `pad_growth = 1.6`: padding cells on each lateral side and
  their geometric growth; `pad_x = n_pad` in x
- `well_positions = (0.25, 0.75)`: positions of the hot and cold wells as
  fractions of the block length

## Rock properties
- `cap_permeability = 1e-3 mD`, `cap_porosity = 0.05`: cap and base rock
- `rock_thermal_conductivity = 2.5`, `cap_thermal_conductivity = 2.0` [W/(m K)]
- `rock_heat_capacity = 900.0` [J/(kg K)], `rock_density = 2600.0` [kg/m³]
- `minporo = 0.01`: lower bound on the SPE10 porosity

## Thermal and operational parameters
- `temperature_charge = 85 °C`, `temperature_discharge = 20 °C`: injection
  temperatures [K]
- `temperature_surface = 10 °C` [K], `thermal_gradient = 0.03` [K/m]
- `thermal_radius = missing`: thermal radius of the hot plume [m] (default a
  quarter of the well distance)
- `rate_charge = missing`: rate of the cold well during charging [m³/s]
  (default the rate that fills a cylinder of the thermal radius over the
  aquifer thickness during the charge period, or a slab of width twice the
  thermal radius in a slice)
- `rate_discharge = rate_charge`: rate of the hot well during discharging
- `charge_period = ["June", "September"]`,
  `discharge_period = ["December", "March"]`, `num_years = 1`,
  `start_year = 2025`, `utes_schedule_args = NamedTuple()`: the schedule (see
  `make_utes_schedule`)

# Returns
A `JutulCase` with the grid indices and setup in `input_data`: `:dims`,
`:dx`, `:dy`, `:dz`, `:aquifer_layers`, `:block_columns`, `:wells` (the `(i, j)`
of the hot and cold wells), `:aquifer_cells`, `:block_cells`, `:rate`,
`:well_distance`, `:thermal_radius`, `:slice`, `:coarsen` and `:timestamps`.
"""
function ates_spe10(;
        formation = :upper_ness,
        layers = SPE10_FORMATIONS[formation],
        ix = 1:60, iy = 1:220,
        slice = missing,
        coarsen = (2, 2, 2),
        minporo = 0.01,
        depth_top = 800.0meter,
        cap = [2.0, 4.0, 8.0, 16.0, 32.0, 64.0] .* meter,
        base = [2.0, 4.0, 8.0, 16.0, 32.0, 64.0] .* meter,
        n_pad = 4, pad_growth = 1.6,
        pad_x = ismissing(slice) ? n_pad : 0,
        cap_permeability = 1.0e-6 * darcy, cap_porosity = 0.05,
        rock_thermal_conductivity = 2.5watt / (meter * Kelvin),
        cap_thermal_conductivity = 2.0watt / (meter * Kelvin),
        rock_heat_capacity = 900.0joule / (kilogram * Kelvin),
        rock_density = 2600.0kilogram / meter^3,
        well_positions = (0.25, 0.75),
        temperature_charge = convert_to_si(85.0, :Celsius),
        temperature_discharge = convert_to_si(20.0, :Celsius),
        temperature_surface = convert_to_si(10.0, :Celsius),
        thermal_gradient = 0.03Kelvin / meter,
        thermal_radius = missing,
        rate_charge = missing,
        rate_discharge = rate_charge,
        charge_period = ["June", "September"],
        discharge_period = ["December", "March"],
        num_years = 1, start_year = 2025,
        utes_schedule_args = NamedTuple()
    )
    ismissing(slice) || (ix = slice:slice)
    blk = spe10_block(layers; ix = ix, iy = iy, coarsen = coarsen, minporo = minporo)
    mx, my, ml = size(blk.poro)

    # ## Grid: the block with padding around it, cap rock above and base rock below
    padded(n, h, np = n_pad) = [reverse([h * pad_growth^i for i in 1:np]); fill(h, n); [h * pad_growth^i for i in 1:np]]
    dx, dy = padded(mx, blk.dx, pad_x), padded(my, blk.dy)
    dz = [reverse(cap); fill(blk.dz, ml); base]
    nx, ny, nz = length(dx), length(dy), length(dz)
    i_blk, j_blk = (pad_x + 1):(pad_x + mx), (n_pad + 1):(n_pad + my)
    k_aq = (length(cap) + 1):(length(cap) + ml)
    mesh = CartesianMesh((nx, ny, nz), (dx, dy, dz), origin = [0.0, 0.0, depth_top - sum(cap)])
    cell(i, j, k) = i + (j - 1) * nx + (k - 1) * nx * ny
    nc = nx * ny * nz

    # ## Properties: SPE10 in the block, layer means in the padding, cap and base rock elsewhere
    perm = fill(cap_permeability, 3, nc)
    phi = fill(cap_porosity, nc)
    lambda = fill(cap_thermal_conductivity, nc)
    gmean(x) = exp(mean(log.(x)))
    for (l, k) in enumerate(k_aq)
        layer_mean = (gmean(blk.kx[:, :, l]), gmean(blk.ky[:, :, l]), gmean(blk.kz[:, :, l]), mean(blk.poro[:, :, l]))
        for j in 1:ny, i in 1:nx
            c = cell(i, j, k)
            if i in i_blk && j in j_blk
                a, b = i - pad_x, j - n_pad
                perm[:, c] .= (blk.kx[a, b, l], blk.ky[a, b, l], blk.kz[a, b, l])
                phi[c] = blk.poro[a, b, l]
            else
                perm[:, c] .= layer_mean[1:3]
                phi[c] = layer_mean[4]
            end
            lambda[c] = rock_thermal_conductivity
        end
    end
    domain = reservoir_domain(
        mesh;
        porosity = phi,
        permeability = perm,
        rock_thermal_conductivity = lambda,
        rock_heat_capacity = rock_heat_capacity,
        rock_density = rock_density
    )

    # ## Wells over the full aquifer thickness, on the long axis
    i_w = first(i_blk) + mx ÷ 2
    j_hot, j_cold = [first(j_blk) + round(Int, f * (my - 1)) for f in well_positions]
    hot = setup_well(domain, [cell(i_w, j_hot, k) for k in k_aq], name = :Hot, simple_well = false)
    cold = setup_well(domain, [cell(i_w, j_cold, k) for k in k_aq], name = :Cold, simple_well = false)
    model = setup_reservoir_model(domain, :geothermal, wells = [hot, cold])
    rmodel = reservoir_model(model)
    rho = rmodel.system.rho_ref[1]

    # ## Hydrostatic pressure and geothermal temperature on the outer faces
    # (in a slice, not on the faces across the slice)
    p(z) = 1atm .+ rho * gravity_constant .* z
    T(z) = temperature_surface .+ thermal_gradient .* z
    geo = tpfv_geometry(mesh)
    keep = nx == 1 ? abs.(geo.boundary_normals[1, :]) .< 0.5 : trues(length(geo.boundary_neighbors))
    zb = geo.boundary_centroids[3, keep]
    bc = flow_boundary_condition(geo.boundary_neighbors[keep], rmodel.data_domain, p(zb), T(zb))
    zc = geo.cell_centroids[3, :]
    state0 = setup_reservoir_state(model, Pressure = p(zc), Temperature = T(zc))

    # ## Rate: the injected heat fills a cylinder (slab in a slice) of the thermal radius
    ch_start = process_time(start_year, charge_period[1])
    ch_end = process_time(start_year, charge_period[2], true)
    ch_end < ch_start && (ch_end += Dates.Year(1))
    charge_duration = (ch_end - ch_start).value * 1.0e-3
    well_distance = (j_cold - j_hot) * blk.dy
    r_th = ismissing(thermal_radius) ? well_distance / 4 : thermal_radius
    Cf = 4184.0 * 1000.0
    phi_aq = mean(blk.poro)
    Caq = Cf * phi_aq + rock_heat_capacity * rock_density * (1 - phi_aq)
    H = ml * blk.dz
    is_slice = nx == 1
    if ismissing(rate_charge)
        rate_charge = (is_slice ? 2r_th * sum(dx) : π * r_th^2) * H * Caq / Cf / charge_duration
    end
    ismissing(rate_discharge) && (rate_discharge = rate_charge)

    # ## Controls: the injector reinjects what the other well produces
    inj(name, other) = InjectorControl(
        JutulDarcy.ReinjectionTarget([other]), [1.0], density = rho,
        temperature = name == :Hot ? temperature_charge : temperature_discharge; check = false
    )
    f_charge = setup_reservoir_forces(
        model, bc = bc,
        control = Dict(:Hot => inj(:Hot, :Cold), :Cold => ProducerControl(TotalRateTarget(-rate_charge)))
    )
    f_discharge = setup_reservoir_forces(
        model, bc = bc,
        control = Dict(:Hot => ProducerControl(TotalRateTarget(-rate_discharge)), :Cold => inj(:Cold, :Hot))
    )
    f_rest = setup_reservoir_forces(model, bc = bc)
    dt, forces, timestamps = make_utes_schedule(
        f_charge, f_discharge, f_rest;
        charge_period = charge_period, discharge_period = discharge_period, num_years = num_years,
        start_year = start_year, utes_schedule_args...
    )

    aquifer_cells = [cell(i, j, k) for k in k_aq for j in 1:ny for i in 1:nx]
    info = Dict(
        :dims => (nx, ny, nz), :dx => dx, :dy => dy, :dz => dz, :aquifer_layers => k_aq,
        :block_columns => (i_blk, j_blk), :wells => ((i_w, j_hot), (i_w, j_cold)), :rate => rate_charge,
        :slice => is_slice, :well_distance => well_distance, :thermal_radius => r_th, :timestamps => timestamps,
        :aquifer_cells => aquifer_cells, :coarsen => coarsen,
        :block_cells => [cell(i, j, k) for k in k_aq for j in j_blk for i in i_blk]
    )
    return JutulCase(model, dt, forces; state0 = state0, input_data = info)
end
