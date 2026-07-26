import numpy as np

def mmoinertia(sec1, sec2, m):
    """
    Moment of inertia for a two-segment rod (tip and tail vectors from COM).
    Ix/Iy: transverse axes — use axial (z) length components only.
    Iz: spin axis — use radial (x,y) components only.
    """
    # Axial (z) distances from COM to each end
    L1z = sec1[2]
    L2z = sec2[2]

    # Radial distances from spin axis for each end
    r1_sq = sec1[0]**2 + sec1[1]**2
    r2_sq = sec2[0]**2 + sec2[1]**2

    Ix = (1/3) * m * (L1z**2 + L2z**2)
    Iy = Ix
    Iz = (1/2) * m * (r1_sq + r2_sq)

    # Avoid degenerate zero spin inertia (perfectly thin rod on axis).
    #
    # NUMERICAL GUARD, NOT PHYSICS: this geometry is a pure line (both ends
    # have zero radial offset), so Iz is exactly 0 and has no physically
    # meaningful value to fall back on. The floor is taken relative to the
    # transverse moment rather than as an absolute constant, because an
    # absolute 1e-6 gave I_inv[2,2] = 1e6 -- a millionfold gain on anything
    # landing in the spin channel. That was harmless only by coincidence:
    # the gyroscopic z-term is wx*wy*(Iy-Ix), which is exactly zero solely
    # because Iy is assigned from Ix verbatim. Any future change that makes
    # the transverse moments differ even slightly would have fed that 1e6.
    #
    # For real spin inertia, model the finite body radius (0.15 m is used as
    # the aero reference radius elsewhere) instead of relying on this floor.
    Iz_floor = 1e-3 * Ix
    if Iz < Iz_floor:
        Iz = Iz_floor

    return np.array([[Ix, 0.0, 0.0],
                     [0.0, Iy, 0.0],
                     [0.0, 0.0, Iz]])

