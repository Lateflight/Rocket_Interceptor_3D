import numpy as np
import linalg_utils as la

G0 = 9.81          # m/s^2, standard gravity
MIN_SPEED = 1e-2   # m/s, below this the airspeed is too small for the drag and
                   # lift direction unit vectors to be numerically meaningful


def g(m):
    """Returns force of gravity in the world frame
    m: Mass
    Return -> F_g"""

    return G0*np.array([0.0,0.0,-1.0])*m


def drag(rho,v,S, drag_coeff):

    if np.linalg.norm(v)>MIN_SPEED:
           # tune this
        drag = -(1/2*rho*S* drag_coeff) * v @ v*la.quick_unit(v)  # quadratic drag
    else:
        drag = np.zeros(3)

    return drag

def s_lift(rho, v, S, lift_coeff):
    """
    Returns lift in local (body-frame) coordinates.

    rho: air density
    v: velocity, expressed in the BODY frame (caller must pass R.T @ vel_world)
    S: surface area of lifting body
    lift_coeff: lift coefficient

    s_ is for simple, see note below
    NOTE: Lift acts perpendicular to the relative wind, in the plane
    containing v and the body's +Z (nose) axis -- the standard
    axisymmetric-body "normal force" direction. It's the component of +Z
    left over after removing its part along v_hat, i.e. z rejected from v.
    This is zero at zero angle of attack (v parallel to +Z) and has no
    dependence on roll (no fin model).

    MAGNITUDE SCALES WITH sin(AoA): `perp` is deliberately used
    UN-normalized, because |z - (z.v_hat)v_hat| is exactly sin(angle
    between the nose axis and the relative wind). That makes this
    L = 1/2 rho S Cl V^2 sin(AoA), the usual linear normal-force model.
    Normalizing perp (as this did previously) gave a constant-magnitude
    force that merely flipped direction through zero AoA -- 1272 N at
    0.01 deg and the same 1272 N at 20 deg. That is survivable while the
    force acts at the CoM and generates no moment, but the moment it puts
    on a center-of-pressure lever arm is then bang-bang, which limit-cycles
    instead of damping.

    The caller decides where this acts: pass the airflow at the center of
    pressure (v + w x r_cp) and apply the result at r_cp to get restoring
    and damping moments; pass the CoM velocity to get a pure force."""

    speed = np.linalg.norm(v)
    if speed > MIN_SPEED:
        v_hat = v / speed
        z = np.array([0.0, 0.0, 1.0])
        perp = z - np.dot(z, v_hat) * v_hat   # |perp| == sin(AoA)
        lift = 1/2*rho*S*lift_coeff*(v@v) * perp
    else:
        lift = np.zeros(3)
    return lift
