# # Radial thermal advection with piston-like analytical solution
#
# This example validates Fimbul's thermal transport for radial single-well
# injection by comparing the numerical simulation against the well-known
# piston-like (sharp-front) thermal-radius formula.
#
# A single injection well is placed at the center of a square aquifer domain.
# Hot fluid is injected at a constant rate, and the thermal plume expands
# radially outward. Under the assumption of advection-dominated transport and
# piston-like displacement, the thermal front radius satisfies:
#
#     r_th(t) = sqrt( Cf * Q * t / (Caq * π * H) )
#
# where Cf = ρ_f Cp_f is the volumetric heat capacity of the injected fluid,
# Caq = Cf ϕ + Cr (1−ϕ) is the effective aquifer heat capacity, Q is the
# volumetric injection rate, and H is the aquifer thickness.
#
# Parameters are read from the accompanying Excel workbook
# `analytical_radial_parameters.xlsx` located in the same directory.

# ## Required packages
# XLSX.jl must be installed:  using Pkg; Pkg.add("XLSX")
using Jutul, JutulDarcy, Fimbul
using XLSX
using GLMakie

# ## Read parameters from xlsx
xlsx_path = joinpath(@__DIR__, "analytical_radial_parameters.xlsx")
@assert isfile(xlsx_path) "Parameters file not found: $xlsx_path"

xf = XLSX.readxlsx(xlsx_path)
ws = xf["Parameters"]

# Build a lookup dict:  parameter_name => value
params = Dict{String, Float64}()
for row in XLSX.eachrow(ws)
    XLSX.row_number(row) == 1 && continue   # skip header
    name = row[1]
    val  = row[2]
    (name isa String && val isa Number) || continue
    params[name] = Float64(val)
end

# ## Convenience unit factors
meter    = si_unit(:meter)
kilogram = si_unit(:kilogram)
second   = si_unit(:second)
day      = si_unit(:day)
Kelvin   = si_unit(:Kelvin)
joule    = si_unit(:joule)
watt     = si_unit(:watt)
darcy    = si_unit(:darcy)
atm      = si_unit(:atm)

to_kelvin = T -> convert_to_si(T, :Celsius)
to_celsius = T -> convert_from_si(T, :Celsius)

# ## Build the case
# Translate xlsx values (SI-friendly units) into the unit system expected by
# analytical_radial.
case, r_th, T_analytical, t = analytical_radial(
    domain_radius           = params["domain_radius"] * meter,
    aquifer_depth_top       = params["aquifer_depth_top"] * meter,
    aquifer_depth_bottom    = params["aquifer_depth_bottom"] * meter,
    porosity                = params["porosity"],
    permeability            = params["permeability"] * 1e-3 * darcy,
    rock_density            = params["rock_density"] * kilogram / meter^3,
    rock_heat_capacity      = params["rock_heat_capacity"] * joule / (kilogram * Kelvin),
    rock_thermal_conductivity = params["rock_thermal_conductivity"] * watt / (meter * Kelvin),
    fluid_density           = params["fluid_density"] * kilogram / meter^3,
    fluid_heat_capacity     = params["fluid_heat_capacity"] * joule / (kilogram * Kelvin),
    temperature_initial     = to_kelvin(params["temperature_initial"]),
    temperature_injection   = to_kelvin(params["temperature_injection"]),
    rate                    = params["rate"] * meter^3 / day,
    num_cells_radial        = Int(params["num_cells_radial"]),
    num_steps               = Int(params["num_steps"]),
    thermal_radius_target   = params["thermal_radius_target"] * meter,
)

# ## Run the simulation
result = simulate_reservoir(case, info_level = 0)

# ## Postprocessing
# Extract the reservoir mesh and geometry
msh  = physical_representation(reservoir_model(case.model).data_domain)
geo  = tpfv_geometry(msh)
xc   = geo.cell_centroids[1, :]   # x-coordinates of cell centers
yc   = geo.cell_centroids[2, :]   # y-coordinates of cell centers

# Identify cells that lie along the positive x-axis (y ≈ 0, x > 0).
# Tolerance: half the smallest cell size.
num_cells_radial = Int(params["num_cells_radial"])
domain_radius    = params["domain_radius"] * meter
r_min_cell       = domain_radius / (10.0 * num_cells_radial) / 2
y_tol            = r_min_cell * 50

radial_idx = findall(yc .> -y_tol .&& yc .< y_tol .&& xc .> 0.0)
sort!(radial_idx; by = i -> xc[i])
r_num = xc[radial_idx]   # numerical radial coordinates

T_inj  = to_kelvin(params["temperature_injection"])
T_init = to_kelvin(params["temperature_initial"])

# ── Helper: find numerical thermal front radius ──────────────────────────────
# Define the front as the location where T crosses the midpoint between
# initial and injection temperature.
T_mid = (T_inj + T_init) / 2.0
function numerical_r_th(state)
    T = state[:Reservoir][:Temperature]
    T_radial = T[radial_idx]
    # Front is between the last cell above T_mid and the first cell below
    above = findlast(T_radial .> T_mid)
    isnothing(above) && return 0.0
    above == length(r_num) && return r_num[end]
    # Linear interpolation
    r1, r2 = r_num[above], r_num[above + 1]
    T1, T2 = T_radial[above], T_radial[above + 1]
    return r1 + (r2 - r1) * (T_mid - T1) / (T2 - T1)
end

# ── Compute thermal radii ────────────────────────────────────────────────────
r_th_numerical   = [numerical_r_th(s) for s in result.states]
r_th_analytical  = r_th.(t)

# ── Plot 1: Temperature profiles at selected time steps ─────────────────────
n_plot   = min(5, length(t))
step_idx = round.(Int, range(1, length(t), length = n_plot))

fig1 = Figure(size = (900, 550), fontsize = 18)
ax1  = Axis(fig1[1, 1];
    xlabel = "Radial distance r (m)",
    ylabel = "Temperature (°C)",
    title  = "Radial temperature profiles: numerical vs analytical")

r_a = range(0.0, domain_radius * 0.9, length = 500)
colors = cgrad(:Spectral, n_plot + 1, categorical = true)

for (ci, k) in enumerate(step_idx)
    T_num   = to_celsius.(result.states[k][:Reservoir][:Temperature][radial_idx])
    T_ana   = to_celsius.(T_analytical.(r_a, t[k]))
    lines!(ax1, r_a, T_ana; linewidth = 6, linestyle = :dash, color = colors[ci],
           label = ci == 1 ? "Analytical" : nothing)
    lines!(ax1, r_num, T_num; linewidth = 2, color = colors[ci],
           label = ci == 1 ? "Numerical" : nothing)
end
axislegend(ax1; position = :rt)
display(fig1)

# ── Plot 2: Thermal radius vs time ───────────────────────────────────────────
t_days = t ./ day

fig2 = Figure(size = (800, 500), fontsize = 18)
ax2  = Axis(fig2[1, 1];
    xlabel = "Time (days)",
    ylabel = "Thermal front radius (m)",
    title  = "Thermal front radius: numerical vs analytical")

lines!(ax2, t_days, r_th_analytical; linewidth = 4, linestyle = :dash,
       color = :black, label = "Analytical  r_th = sqrt(Cf Q t / (Caq π H))")
scatter!(ax2, t_days, r_th_numerical; markersize = 14, color = :steelblue,
         label = "Numerical (Fimbul)")
lines!(ax2, t_days, r_th_numerical; linewidth = 2, color = :steelblue)
axislegend(ax2; position = :lt)
display(fig2)

# ── Print summary ────────────────────────────────────────────────────────────
println("\n=== Radial Thermal Advection – Summary ===")
println("Final time          : $(round(t[end]/day, digits=1)) days")
println("Analytical r_th     : $(round(r_th_analytical[end], digits=1)) m")
println("Numerical  r_th     : $(round(r_th_numerical[end], digits=1)) m")
rel_err = abs(r_th_numerical[end] - r_th_analytical[end]) / r_th_analytical[end] * 100
println("Relative error      : $(round(rel_err, digits=2)) %")
