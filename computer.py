import numpy as np
from scipy.linalg import solve_continuous_are
import linalg_utils as la

# Numerical guards. These are thresholds, not tuning knobs -- each marks the
# point below which a quantity stops being meaningful rather than merely small.
MIN_FORWARD_RANGE = 1e-5   # m, forward component below which the target is
                           # level with or behind the focal plane and cannot
                           # be projected onto it at all
MIN_PERP          = 1e-8   # below this the perpendicular component carries no
                           # usable direction, so the azimuth is arbitrary

DEFAULT_CALIB_SAMPLES = 200  # stationary gyro samples averaged on the pad;
                             # noise falls as sqrt(n)


def sight(tar, Center_world, F_height, F_len, R,N):
    """Pinhole-camera seeker: projects the target onto a flat focal plane
    at standoff F_height, gated by a circular FOV cone of half-angle
    atan(F_len/F_height). Returns [alpha, beta, valid] -- valid is 1.0 for
    a fresh in-cone detection, 0.0 otherwise (target behind the focal
    plane or outside the cone). alpha/beta are only meaningful when
    valid==1.0; callers must not infer "no detection" from alpha==0/beta==0,
    since a genuine on-boresight detection also reads exactly zero."""
    rel = R.T @ (tar-Center_world)
    if rel[2] <= MIN_FORWARD_RANGE:
        return np.array([0.0, 0.0, 0.0])

    t = F_height/rel[2]
    alpha = rel[0]*t
    beta = rel[1]*t

    if rel[0]**2+rel[1]**2 <= (F_len/F_height)**2*rel[2]**2:
        return np.array([alpha, beta, 1.0])
    else:
        #print("Lock broken, attempting to re-lock...")
        return np.array([0.0, 0.0, 0.0])


# ─────────────────────────────────────────
# Strapdown IMU (rate gyro)
# ─────────────────────────────────────────

def imu(w_true, bias=None, noise_sigma=0.0, rng=None):
    """Strapdown 3-axis rate gyro. Returns the MEASURED body angular rate.

    w_true:      true body-frame angular velocity (rad/s) -- ground truth,
                 available here only because this is a simulation.
    bias:        constant per-axis turn-on bias (rad/s). Real MEMS parts have
                 a bias that is repeatable within a run but changes between
                 power cycles, which is why calibrate_bias() exists.
    noise_sigma: per-sample white noise standard deviation (rad/s).
    rng:         np.random.Generator, so runs stay reproducible. If None, the
                 output is noise-free and fully deterministic.

    The gyro is body-FIXED, exactly like the seeker, so it measures the same
    body-frame rates the equations of motion integrate -- no frame conversion
    is needed by the caller.

    Note on why this is worth having: the seeker alone cannot tell "the target
    is drifting across my view" from "I am rotating". It measures
    d(theta)/dt = lambda_dot - omega, one equation in two unknowns. The gyro
    supplies omega, which makes the inertial LOS rate lambda_dot observable.
    See los_rate() below."""
    m = np.asarray(w_true, dtype=float).copy()
    if bias is not None:
        m = m + np.asarray(bias, dtype=float)
    if rng is not None and noise_sigma > 0:
        m = m + rng.normal(0.0, noise_sigma, size=3)
    return m


def calibrate_bias(rng=None, bias=None, noise_sigma=0.0, n=DEFAULT_CALIB_SAMPLES):
    """Estimate gyro bias by averaging samples while the vehicle is stationary.

    This is what a real system does on the pad before launch: with the true
    rate known to be zero, whatever the gyro reads IS the bias. Averaging n
    samples cuts the noise contribution by sqrt(n).

    Bias matters here in a specific, bounded way. Guidance consumes the
    instantaneous RATE, never an integrated attitude, so bias does not
    accumulate without limit the way it does in dead reckoning. It shows up
    instead as a fixed offset in the estimated LOS rate -- a steady lead-angle
    error, i.e. a systematic miss offset rather than a divergence."""
    acc = np.zeros(3)
    for _ in range(n):
        acc += imu(np.zeros(3), bias=bias, noise_sigma=noise_sigma, rng=rng)
    return acc/n


def los_rate(dtheta_a, dtheta_b, dphi, dt_sample):
    """Inertial line-of-sight rate from a body-fixed seeker plus a gyro.

    dtheta_a/dtheta_b: change in the TRUE bearing angles (radians) over the
                       sample interval, i.e. arctan(alpha/F_height) differenced
                       between two consecutive fresh seeker frames.
    dphi:              integral of the measured body rate over that SAME
                       interval (radians) -- the attitude change, 3-vector.
    dt_sample:         length of the interval (s).

    Returns (lambda_dot_a, lambda_dot_b): the rotation rate of the line of
    sight in inertial space.

    SIGN CONVENTION. A rotation of the body by +phi_y swings the boresight
    toward +x, which makes the target's measured x-bearing DECREASE:
        dtheta_a = dlambda_a - dphi_y   ->   dlambda_a = dtheta_a + dphi_y
    A rotation by +phi_x swings the boresight toward -y, so the y-bearing
    INCREASES:
        dtheta_b = dlambda_b + dphi_x   ->   dlambda_b = dtheta_b - dphi_x

    IMPORTANT -- combine BEFORE filtering. Both terms must be measured over
    the same interval and summed raw. Filtering the seeker term but not the
    gyro term (or vice versa) leaves a residual of omega - filtered(omega),
    which does not vanish and injects the vehicle's own rotation straight into
    the guidance command. Filter the SUM afterwards if smoothing is wanted."""
    return ((dtheta_a + dphi[1]) / dt_sample,
            (dtheta_b - dphi[0]) / dt_sample)



def los_change(sight_old,sight_new,dt):

    v1 = sight_old
    v2 = sight_new

    cross = np.cross(v1, v2)
    dot = np.dot(v1, v2)

    angle = np.atan2(np.linalg.norm(cross), dot)

    # axis-angle vector (what you actually want)
    if np.linalg.norm(cross) != 0:
        axis = cross / np.linalg.norm(cross)
        angle_vec = axis * angle
    else:
        angle_vec = np.zeros(3)  # vectors are parallel or one is zero
    
    
    los_roc = angle_vec/dt
    
    return los_roc
    


def clamp_gimbal(F_desired, thrust_mag, theta_max_rad):
    # Normalize desired direction
    norm = np.linalg.norm(F_desired)
    if norm == 0:
        return np.array([0.0, 0.0, thrust_mag])  # default: straight up

    d = F_desired / norm

    # Engine axis (+Z)
    a = np.array([0.0, 0.0, 1.0])

    cos_theta = np.dot(d, a)

    # If already within cone, just use it
    if cos_theta >= np.cos(theta_max_rad):
        v = d
    else:
        # Decompose into parallel and perpendicular components
        d_parallel = cos_theta * a
        d_perp = d - d_parallel

        perp_norm = np.linalg.norm(d_perp)
        if perp_norm < MIN_PERP:
            # Direction is basically opposite or aligned; just clamp to edge
            v = np.array([
                np.sin(theta_max_rad),
                0.0,
                np.cos(theta_max_rad)
            ])
        else:
            u = d_perp / perp_norm
            v = np.cos(theta_max_rad) * a + np.sin(theta_max_rad) * u

    return thrust_mag * v