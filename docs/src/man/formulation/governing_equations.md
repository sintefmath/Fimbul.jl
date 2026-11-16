# Governing equations in geothermal energy

Fimbul currently supports pure H₂O as a single-component, two-phase
fluid system. Water may be present in a liquid phase ($L$) and/or a vapor phase
($V$).

## Conservation of mass
Mass conservation is written as

```math
    \partial_t\left([1-\phi]\rho_R + \phi[\rho_L s_L + \rho_V s_V]\right)
    + \nabla \cdot \left(\rho_L \vec{v}_L + \rho_V \vec{v}_V\right)
    = Q_M,
```

where $\phi$ is porosity, $\rho_\alpha$ is density for the rock ($\alpha = R$)
or fluid phase ($\alpha \in \{L, V\}$), and $s_\alpha$ is phase saturation
(the fraction of pore volume occupied by phase $\alpha$, with $s_L + s_V = 1$).
The term $Q_M$ denotes mass sources and sinks, including wells and boundary
fluxes. Phase velocities are given by Darcy's law,

```math
    \vec{v}_\alpha = - \lambda_\alpha \mathbf{K} \left(\nabla p - \rho_\alpha g \nabla z\right),
```

where $\lambda_\alpha = k_{r,\alpha}/\mu_\alpha$ is phase mobility (ratio of
relative permeability to viscosity), $\mathbf{K}$ is absolute permeability, $p$
is pressure, $g$ is gravitational acceleration, and $z$ is elevation.

> NOTE: In practice, the interface between liquid and vapor phases in an
> H₂O system may exhibit non-negligible capillary pressure effects.
> This is currently not modeled in Fimbul; support is planned for a future
> release.

## Conservation of thermal energy
Thermal energy conservation is expressed as

```math
    \partial_t\left([1-\phi]\rho_R u_R + \phi[\rho_L s_L u_L + \rho_V s_V u_V]\right)
    + \nabla \cdot \left(\rho_L h_L \vec{v}_L + \rho_V h_V \vec{v}_V\right)
    - \nabla \cdot \left(\mathbf{\Lambda}\nabla T\right)
    = Q_H,
```

where $u_\alpha$ is specific internal energy, $h_\alpha = u_\alpha +
p/\rho_\alpha$ is specific enthalpy, $T$ is temperature, and $\mathbf{\Lambda}$
is the effective thermal conductivity tensor of the rock–fluid mixture. The
second term on the left represents advective heat transport by the flowing
phases; the third term is conductive heat transport via Fourier's law,

```math
    \vec{q}_c = - \mathbf{\Lambda}\nabla T.
```

The term $Q_H$ denotes thermal energy sources and sinks.

## Fluid properties
The fluid properties of H₂O—such as density, viscosity, enthalpy,
and phase state—depend on both pressure and enthalpy and are provided through
precomputed lookup tables.

For the high-enthalpy (two-phase) regime, the steam tables are generated using
the [`CoolProp`](https://www.coolprop.org/) library [bell_pure_2014](@cite) via
the `FimbulCoolPropExt` extension. They are stored as `Artifacts` that are
loaded automatically when using the high-enthalpy formulation. The tables can
also be loaded directly using the `Artifacts` API. To generate your own tables,
see `build_steam_tables_h2o` in `FimbulCoolPropExt`.

For the low-enthalpy (single-phase liquid) regime, Fimbul uses PVT tables from
JutulDarcy generated using the NIST database.

The figure below shows the temperature, density, and vapor saturation of
H₂O as functions of pressure and specific enthalpy, covering the
full range of both single-phase and two-phase states. The saturation line
separates the single-phase liquid region (left), two-phase region (center), and
single-phase vapor region (right). These phase diagram contours can be
visualized using `plot_phase_diagram_contours`.

## Single- and two-phase systems
Most Fimbul examples currently operate in the low-enthalpy regime, i.e.,
single-phase liquid conditions. The high-enthalpy examples use the two-phase
steam-table formulation described above.

## Multi-continuum and fractures
The formulation above assumes local thermal equilibrium between rock and fluid,
i.e., a single shared temperature field $T$. This is appropriate for many
porous-medium geothermal applications, but may not be adequate when flow is
dominated by fractures or other high-contrast heterogeneities.

To model fractured systems, Fimbul can be combined with the discrete fracture
modeling framework in JutulDarcy—see e.g., 
[the FTES setup function](@ref ftes).

## Wells
Flow in wells can be modelled in two different ways:
- As a simple well where flow in and out of the wellbore is computed based on an
  explicitly computed pressure drawdown from the top to each perforated cell
  -- see JutulDarcy's [`SimpleWell`](https://sintefmath.github.io/JutulDarcy.jl/stable/man/basics/wells#JutulDarcy.SimpleWell). In this
  formulation, the entire wellbore gets a single temperature, which in many
  cases is a poor approximation for geothermal wells.
- As a multisegment well where each wellbore is discretized as a 1D object with
  its own conservation equations for mass and energy, with accompanying pressure
  and heat variables (either temperature or enthalpy depending on the
  formulation) -- see JutulDarcy's
  [`MultiSegmentWell`](https://sintefmath.github.io/JutulDarcy.jl/stable/man/basics/wells#JutulDarcy.MultiSegmentWell).

In geothermal energy, the main product (heat) is in principle leaking from or
into the formation along the entire wellbore, so a multisegment well
formulation is the preferred modelling choice. For closed-loop geothermal wells
(e.g., AGS, BHE and BTES), a multisegment well formulation is necessary, since
temperature changes along the wellbore are crucial to model. All examples in
Fimbul use the multisegment well formulation.

Both simple and multisegment wells (with the exception of closed-loop wells)
compute mass and conductive heat flow using Peaceman-type well indices. Mass
flow between the well and the reservoir is handled by JutulDarcy, using a well
index $\mathrm{WI}$ computed with Peaceman's formula -- see
[`compute_peaceman_index`](https://sintefmath.github.io/JutulDarcy.jl/stable/man/basics/wells#JutulDarcy.compute_peaceman_index) and the
[JutulDarcy documentation on wells](https://sintefmath.github.io/JutulDarcy.jl/stable/man/basics/wells). For a perforation connecting well
node $w$ to reservoir cell $c$, the heat flow from the reservoir into the well
is the sum of advective and conductive heat flow,

```math
    q_H = \sum_\alpha h_\alpha q_\alpha + \mathrm{WI}_{th}\left(T_c - T_w\right),
```

where $q_\alpha$ is the mass flow of phase $\alpha$ and the phase enthalpy
$h_\alpha$ is taken from the reservoir cell for production and from the well
for injection.

### Thermal well index
The thermal well index $\mathrm{WI}_{th}$ accounts for conduction through the
well casing, the grout surrounding the casing, and the rock. These are treated
as thermal resistances in series,

```math
    \mathrm{WI}_{th} = \frac{2\pi L}{R_c + R_g + R_r},
```

where the casing, grout and rock terms are

```math
    R_c = \frac{\ln(r_w/r_{c})}{\lambda_c}, \qquad
    R_g = \frac{\ln(r_g/r_w)}{\lambda_g}, \qquad
    R_r = \frac{\ln(r_e^{th}/r_g)}{\Lambda_e}.
```

Here, $r_{c} = r_w - t_c$ is the inner radius of a casing of thickness
$t_c$, $r_g = r_w + t_g$ is the outer radius of a grout layer of thickness
$t_g$, and $\lambda_c$ and $\lambda_g$ are the casing and grout thermal
conductivities. A term is omitted if the corresponding thickness is zero.
The rock term uses Peaceman's formula with thermal conductivity in place of
permeability: $\Lambda_e = \sqrt{\Lambda_1\Lambda_2}$ is computed from the
effective rock–fluid conductivity
$\Lambda = \phi\Lambda_f + (1-\phi)\Lambda_R$, and $r_e^{th}$ is the
Peaceman equivalent radius computed with $\Lambda$ in place of $\mathbf{K}$
and the Peaceman constant doubled ($0.56$ instead of $0.28$). Both indices
can also be set directly when setting up
the well (`WI` and `WIth` in `setup_well`).

### Closed-loop wells
Closed-loop wells (see [`setup_closed_loop_well`](@ref)) do not exchange mass with the
reservoir, so $\mathrm{WI} = 0$, and heat is exchanged by conduction only.
Fimbul represents the cross-section of the borehole as a network of nodes
(pipes and grout) connected by thermal resistances, following the thermal
resistance and capacity model described in [feflow_btes](@cite). The heat
flow between two connected nodes $a$ and $b$ (or between a grout node and the
reservoir) is

```math
    q_{H} = \frac{L}{R_{ab}}\left(T_a - T_b\right),
```

where $R_{ab}$ is the thermal resistance per unit length. In the expressions
below, $d_g = 2r_g$ is the borehole diameter, $\lambda_g$ is the grout
conductivity, and subscripts $i$ and $o$ denote the inner and outer radius of a
pipe wall. Advective resistances inside the pipes are currently neglected.

**Generic closed-loop well** (`ClosedLoopGenericWell`): the pipe exchanges heat
directly with the reservoir through $\mathrm{WI}_{th}$, as computed above with
the pipe wall as casing and the surrounding grout.

**Single U-tube** (`ClosedLoop1UWell`): two pipes with outer diameter
$d_p = 2r_{p,o}$ and shank spacing $w$, each connected to its own grout node.
With

```math
    x = \frac{\ln\left(\sqrt{d_g^2 + 2d_p^2}/(2d_p)\right)}{\ln\left(d_g/(\sqrt{2}d_p)\right)}, \qquad
    R_g = \left(1.601 - 0.888\frac{w}{d_g}\right)
    \frac{\cosh^{-1}\left(\frac{d_g^2 + d_p^2 - w^2}{2d_g d_p}\right)}{2\pi\lambda_g},
```

the resistances are

```math
\begin{aligned}
    R_{pg} &= \frac{\ln(r_{p,o}/r_{p,i})}{2\pi\lambda_p} + xR_g
    &&\text{(pipe to grout)},\\
    R_{gr} &= (1-x)R_g &&\text{(grout to reservoir)},\\
    R_{gg} &= \frac{2R_{gr}\left(R_{ar} - 2xR_g\right)}{2R_{gr} - R_{ar} + 2xR_g}
    &&\text{(grout to grout)},
\end{aligned}
```

where
$R_{ar} = \cosh^{-1}\left((2w^2 - d_p^2)/d_p^2\right)/(2\pi\lambda_g)$.

**Coaxial** (`ClosedLoopCoaxialWell`): an inner pipe inside an outer pipe
(the annulus), surrounded by grout. With $d_p = 2r_{op,o}$ the outer
diameter of the outer pipe,

```math
    x = \frac{\ln\sqrt{(d_g^2 + d_p^2)/(2d_p^2)}}{\ln(d_g/d_p)}, \qquad
    R_g = \frac{\ln(d_g/d_p)}{2\pi\lambda_g},
```

the resistances are

```math
\begin{aligned}
    R_{pp} &= \frac{\ln(r_{ip,o}/r_{ip,i})}{2\pi\lambda_{ip}}
    &&\text{(inner pipe to annulus)},\\
    R_{pg} &= \frac{\ln(r_{op,o}/r_{op,i})}{2\pi\lambda_{op}} + xR_g
    &&\text{(annulus to grout)},\\
    R_{gr} &= (1-x)R_g &&\text{(grout to reservoir)}.
\end{aligned}
```

For the U-tube and coaxial wells, the grout-to-reservoir index $L/R_{gr}$ is
used as the thermal well index $\mathrm{WI}_{th}$ of the perforation, and all
indices can be overridden when setting up the well.

> NOTE: For the U-tube and coaxial wells, $\mathrm{WI}_{th}$ only accounts for
> the thermal resistance inside the borehole. Unlike the Peaceman-type index
> above (used by the generic closed-loop well), it contains no term for
> conduction from the borehole wall to the center of the reservoir cell, i.e.,
> the temperature of the perforated cell is assumed to approximate the borehole
> wall temperature. This is the same assumption as in [feflow_btes](@cite), and
> requires the reservoir grid to be refined around the well.
