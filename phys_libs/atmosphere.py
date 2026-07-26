import numpy as np

# ─────────────────────────────────────────
# ICAO / ISA standard atmosphere constants
# ─────────────────────────────────────────
T0   = 288.15      # K, sea-level standard temperature
P0   = 1e5         # Pa, sea-level pressure (ISA is 101325; 1e5 used here)
RHO0 = 1.25        # kg/m^3, sea-level density

M_AIR = 0.0289644  # kg/mol, molar mass of dry air
R_GAS = 8.3144     # J/(mol K), universal gas constant
G0    = 9.81       # m/s^2, standard gravity

# Layer boundaries (m) and lapse rates (K/m). Sign is part of the value: the
# troposphere cools with altitude, the stratosphere warms.
Z_TROPOPAUSE  = 11000
Z_STRAT1      = 20000
Z_STRAT2      = 32000
Z_STRATOPAUSE = 47000     # model ceiling; held constant above this
LAPSE_TROPO   = -0.0065
LAPSE_STRAT1  = +0.0010
LAPSE_STRAT2  = +0.0028


def T(z):
    """International Standard Atmosphere temperature model, in Kelvin.
    z: geopotential altitude in meters.

    Layer boundaries and lapse rates (ICAO standard atmosphere):
      0 - 11,000 m:  T0=288.15K, lapse -6.5 K/km   (troposphere)
      11,000 - 20,000 m:  isothermal at 216.65K     (tropopause)
      20,000 - 32,000 m:  lapse +1.0 K/km           (stratosphere)
      32,000 - 47,000 m:  lapse +2.8 K/km           (stratosphere)
      above 47,000 m:     isothermal at 270.65K     (stratopause, model ceiling)

    Each branch is anchored to the temperature at the START of its own
    layer plus the lapse rate times distance INTO that layer (z - boundary),
    so T(z) is continuous across every boundary. Previously the 20-30km and
    30-50km branches used the lapse rate times the raw altitude z instead of
    distance-into-layer, which caused the temperature to climb steeply and
    become physically wrong (60C+ too warm) through the stratosphere, then
    jump by ~270K at the old 50,000m seam because the last branch was an
    unrelated formula that didn't connect to the one before it.
    """
    T11 = T0 + LAPSE_TROPO*(Z_TROPOPAUSE - 0)          # temperature at 11,000 m
    T20 = T11                                          # isothermal, so T20 == T11
    T32 = T20 + LAPSE_STRAT1*(Z_STRAT2 - Z_STRAT1)     # temperature at 32,000 m

    if z <= Z_TROPOPAUSE:
        return T0 + LAPSE_TROPO*z
    elif z <= Z_STRAT1:
        return T11
    elif z <= Z_STRAT2:
        return T20 + LAPSE_STRAT1*(z - Z_STRAT1)
    elif z <= Z_STRATOPAUSE:
        return T32 + LAPSE_STRAT2*(z - Z_STRAT2)
    else:
        # Model ceiling: hold the 47,000 m value rather than extrapolating
        # indefinitely (stratopause is approximately isothermal in reality).
        return T32 + LAPSE_STRAT2*(Z_STRATOPAUSE - Z_STRAT2)

def return_atmo_state(z):
    """Returns (pressure Pa, density kg/m^3, temperature K) at altitude z.

    NOTE: no floor at z = 0. Called with negative altitude this extrapolates
    the exponential downward and returns density INCREASING with depth, which
    is not physical. Callers are expected to stop at the ground themselves.
    physics_engine.py does this via GROUND_ALT."""
    Tem = T(z)
    scale = np.e**(-M_AIR*G0*z/(R_GAS*Tem))   # barometric factor, shared by both
    return (P0*scale, RHO0*scale, Tem)