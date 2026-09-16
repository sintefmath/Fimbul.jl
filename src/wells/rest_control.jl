"""
    ClosedLoopRestControl(pressure)
    ClosedLoopRestControl(target::BottomHolePressureTarget)

Well control for a closed-loop well that is at rest, i.e. one that is not
circulating. The top node of the well is held at the given wellhead pressure,
and mass is free to flow either into or out of the well as the fluid in the
loop contracts or expands with temperature. Fluid that enters the loop carries
the enthalpy of the top node, so the exchange is energy-neutral.

This mimics the expansion vessel of a real closed-loop installation. Without it,
a shut-in loop is a sealed, rigid volume, and a temperature change of the fluid
it contains must be balanced by pressure alone (of order 10 bar/K for water),
which cannot be resolved for cooling loops in particular.

The control is meant to replace `DisabledControl` for the supply well of a
closed loop during rest periods. The return well should be disabled, so that
the loop has a single connection to the surface.

See also [`DisabledControl`](@ref), [`InjectorControl`](@ref),
[`ProducerControl`](@ref).
"""
struct ClosedLoopRestControl{T<:JutulDarcy.WellTarget} <: JutulDarcy.WellControlForce
    target::T
end

function ClosedLoopRestControl(pressure::Real)
    return ClosedLoopRestControl(BottomHolePressureTarget(pressure))
end

# The surface rate may take either sign, so no adjustment is made when the
# control is activated (the generic update leaves the sign free as well).
JutulDarcy.valid_surface_rate_for_control(q_t, ::ClosedLoopRestControl) = q_t

# No limits: the rate is whatever the loop needs to hold the pressure.
JutulDarcy.default_limits(::ClosedLoopRestControl) = nothing
JutulDarcy.check_well_limit(name::Symbol, limit_value, cond, ::ClosedLoopRestControl) = missing

# Only used by JutulDarcy's well output, which evaluates the well against every
# output target type in turn. The control itself always holds a pressure target.
JutulDarcy.replace_target(::ClosedLoopRestControl, target::JutulDarcy.WellTarget) = ClosedLoopRestControl(target)

JutulDarcy.effective_surface_rate(qts, ::ClosedLoopRestControl) = qts

function Base.show(io::IO, c::ClosedLoopRestControl)
    print(io, "ClosedLoopRestControl($(c.target))")
end
