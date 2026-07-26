import numpy as np

def apply_rot(tau,w,I_inv,I_principal):

    return I_inv @ (tau - np.cross(w, I_principal @ w))



def step(tau,w,I_inv,I_principal,dt):
    """Classic RK4 step for angular velocity w under torque tau.
    Returns the updated w. Does not mutate the w passed in (previously
    this used w += ..., which relied on NumPy's in-place-add semantics
    to update the caller's array as a side effect and returned None;
    callers that assign the return value now get the same value either
    way, but callers that ignored the return value and expected the
    passed-in array to be updated in place no longer work silently)."""

    k1= apply_rot(tau,w,I_inv,I_principal)
    k2 = apply_rot(tau,w+0.5*dt*k1,I_inv,I_principal)
    k3 = apply_rot(tau,w+0.5*k2*dt,I_inv,I_principal)
    k4 = apply_rot(tau,w+dt*k3,I_inv,I_principal)

    w_new = w + (dt/6.0) * (k1 + 2*k2 + 2*k3 + k4)
    return w_new