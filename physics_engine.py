from phys_libs import grapher as grp
import numpy as np
from phys_libs import PiForce as PF
from phys_libs import Rotator_module as rt
from phys_libs import linalg_utils as la
from phys_libs import Rocket as rock
from phys_libs import computer as cmp
from phys_libs import atmosphere as atm

"""
Reminders
state:
x,y,z of Center
x,y,z of top
x,y,z of bottom
"""

# ─────────────────────────────────────────
# Simulation Parameters
# ─────────────────────────────────────────
dt = 0.01                          # s, integration step (100 Hz)

# ─────────────────────────────────────────
# Thrust (in principal/body frame, from COM)
# ─────────────────────────────────────────
thrust_mag = 4000                  # N, constant while fuel remains
thrust_direction = la.quick_unit(np.array([0.0, 0.0, 1.0]))   # body +Z, the nose
thrust = thrust_mag * thrust_direction                        # N

# ─────────────────────────────────────────
# Rocket Definition
# ─────────────────────────────────────────
Top    = np.array([0, 0, 4.5])     # nose, body frame
Bottom = np.array([0, 0, 0])       # nozzle, body frame
m_fuel = 7/2                       # kg, TOTAL fuel across both sections
m = 5 + m_fuel/2                   # kg PER SECTION; total vehicle mass is 2*m
Center = (Top+Bottom)/2            # centre of mass, body frame

# ─────────────────────────────────────────
# Aerodynamic reference values
# ─────────────────────────────────────────
BODY_RADIUS   = 0.15                    # m, airframe reference radius
S_REF         = np.pi*BODY_RADIUS**2    # m^2, reference area for drag AND lift
DRAG_COEFF    = 0.5
LIFT_COEFF    = 0.2
VISCOUS_DAMP  = 0.2                     # N m per rad/s, artificial body-rate damping

TARGET_RADIUS     = 0.5                 # m
S_REF_TARGET      = np.pi*TARGET_RADIUS**2
DRAG_COEFF_TARGET = 0.09

# ─────────────────────────────────────────
# World-Frame Geometry
# ─────────────────────────────────────────
Center_world = np.array([0,0,2.25])                  # m, CoM position in world
Top_world    = Center_world+np.array([0,0,2.25])     # m, nose position in world
Bottom_world = Center_world-np.array([0,0,2.25])     # m, nozzle position in world
Rk = rock.Rocket(m, Bottom-Center, Top-Center, thrust,BODY_RADIUS)   # holds I_principal, I_inv

# ─────────────────────────────────────────
# Initial Rotation
# ─────────────────────────────────────────
R = np.identity(3)                 # body -> world rotation matrix

# ─────────────────────────────────────────
# Initial State
# ─────────────────────────────────────────
w   = np.array([0.0, 0.0, 0.0])    # rad/s, BODY-frame angular velocity
vel = np.zeros(3)                  # m/s, world-frame velocity
tar = np.array([-6000, 00, 10000.0]) # m, target position (world)
vel_tar = np.array([500,00.0,0])   # m/s, target velocity (world)
tar_static = tar.copy()            # target start point, for the plot marker

# ─────────────────────────────────────────
# Dummy lists for debugging, if needed or if you wish to export to CSV
# ─────────────────────────────────────────
U = []          # Top_world per step
D = []          # Bottom_world per step
X = []          # range to target per step
Q = []          # world-frame thrust vector per step (collected, not plotted)

# ─────────────────────────────────────────
# Data Collection Buffers
# ─────────────────────────────────────────
TAR_log  = []               # target position per step
LOS_view = []               # world-frame LOS unit vector per step (quiver)
LOS_log  = []               # [alpha, beta] seeker-plane bearing per step
closest  = float('inf')     # m, best range achieved so far

# ─────────────────────────────────────────
# Simulation parameters
# ─────────────────────────────────────────
SIM_LEN     = 45000    # max iterations (25 s at dt = 0.01)
HIT_RADIUS  = 5.0     # m, range below which the engagement counts as an intercept
GROUND_ALT  = 0.0     # m, both bodies stop here (atmosphere.py has no floor)
ROT_EPS     = 1e-18   # rad, below this a Rodrigues step is treated as no rotation

mg     = PF.g(m)      # N, gravity on one section (recomputed in loop as m depletes)
m_tar  = 1000         # kg, target mass (constant)
mg_tar = PF.g(m_tar)  # N, gravity on the target
c_av = 0              # previous frame's raw alpha (for the bearing difference)
c_bv = 0              # previous frame's raw beta

# ─────────────────────────────────────────
# Seeker
# ─────────────────────────────────────────
F_height = 4          # focal-plane standoff
F_len    = 1          # focal-plane radius; FOV half-angle = atan(F_len/F_height) = 14.04 deg
                      #   -- must stay above the lead angle the geometry needs (~12 deg)
SEEKER_LAG_STEPS = 2  # seeker refreshes every Nth step; also sets dt_sample below

theta_max = np.radians(10)   # rad, mechanical gimbal limit
fuel_depletion_rate = 0.04   # kg/s, total
d_alpha_f = 0                # filtered bearing rate, x (pursuit only)
d_beta_f  = 0                # filtered bearing rate, y (pursuit only)
deflect = np.zeros(2)        # rad, commanded gimbal deflection, held between frames

# ─────────────────────────────────────────
# Onboard IMU (strapdown rate gyro)
# ─────────────────────────────────────────
GYRO_SEED          = 12345   # fixed so runs stay reproducible
GYRO_CALIB_SAMPLES = 200     # stationary samples averaged for the pad calibration
gyro_rng   = np.random.default_rng(GYRO_SEED)
gyro_bias  = np.array([0.004, -0.003, 0.002])   # rad/s, turn-on bias (~0.2 deg/s)
gyro_noise = 0.0015                             # rad/s, per-sample white noise
gyro_bias_est = cmp.calibrate_bias(rng=gyro_rng, bias=gyro_bias,
                                   noise_sigma=gyro_noise, n=GYRO_CALIB_SAMPLES)
gyro_accum = np.zeros(3)     # rad, integrated rate since the last seeker frame
lam_a = 0.0                  # rad/s, inertial LOS rate, x channel
lam_b = 0.0                  # rad/s, inertial LOS rate, y channel
LOS_rate_log = []            # [t, lam_a, lam_b, seeker-only rate a, seeker-only rate b]

# ─────────────────────────────────────────
# Guidance mode
# ─────────────────────────────────────────
GUIDANCE = 'pn'       # 'pn' = proportional navigation | 'pursuit' = old bearing PD law

NAV = 5.0             # navigation constant N (PN only)
KW  = 0.10            # inner rate-loop gain (PN only)
KP_PURSUIT  = 0.1     # pursuit only
KD_PURSUIT  = 0.25    # pursuit only
TAU_SAMPLES = 2       # derivative filter time constant, in seeker sample periods

wdes_x = 0.0          # rad/s, commanded body rate about x
wdes_y = 0.0          # rad/s, commanded body rate about y

# ─────────────────────────────────────────
# Launch state
# ─────────────────────────────────────────

launched = False      #Is the vehicle launched?
prev_valid = False    # was the PREVIOUS seeker frame a fresh detection? A bearing
                      
valid = 0.0          #NOTE: Consider replacing with a proper boolean?
                     

print("gyro bias  true %s" % np.round(gyro_bias, 5))
print("     estimated %s   residual %s rad/s"
      % (np.round(gyro_bias_est, 5), np.round(gyro_bias - gyro_bias_est, 5)))

# ─────────────────────────────────────────
# Loop
# ─────────────────────────────────────────
for i in range(SIM_LEN):

    if np.linalg.norm(tar-Center_world)<closest:
        closest=np.linalg.norm(tar-Center_world)
    if np.linalg.norm(tar-Center_world)<HIT_RADIUS:
        print(f"Simulation terminated at {i} steps ({i*dt}s)")
        break

    if tar[2] <= GROUND_ALT:
        print("target hit the ground at t = %.2f s" % (i*dt))
        break
    if Center_world[2] <= GROUND_ALT:
        print("interceptor hit the ground at t = %.2f s" % (i*dt))
        break

    # ─────────────────────────────────────────
    # Read the gyro ONCE per step
    # ─────────────────────────────────────────
    w_meas = cmp.imu(w, bias=gyro_bias, noise_sigma=gyro_noise, rng=gyro_rng) - gyro_bias_est

    # ─────────────────────────────────────────
    # INNER RATE LOOP -- every step, off the gyro
    # ─────────────────────────────────────────
    if GUIDANCE == 'pn':
        deflect = np.array([-KW*(wdes_y - w_meas[1]),
                             KW*(wdes_x - w_meas[0])])

    # ─────────────────────────────────────────
    # Lever arm definition
    # ─────────────────────────────────────────
    r  = Bottom - Center      # CoM -> nozzle, where thrust acts

    # ─────────────────────────────────────────
    # Computer return
    # ─────────────────────────────────────────
    if ((i) % SEEKER_LAG_STEPS == 0):
        alpha, beta, valid = cmp.sight(tar,Center_world,F_height,F_len,R,1)
        if not valid:
            alpha = c_av      # no fresh detection: hold the last known bearing
            beta = c_bv

        # ─────────────────────────────────────────
        # Bearing recovery (used by BOTH guidance modes)
        # ─────────────────────────────────────────
        theta_a = np.arctan(alpha/F_height)     # rad, true off-boresight angle, x
        theta_b = np.arctan(beta/F_height)      # rad, true off-boresight angle, y

        Kp = KP_PURSUIT
        Kd = KD_PURSUIT
        dt_sample = SEEKER_LAG_STEPS*dt         # s, must match the gate above
        tau = TAU_SAMPLES*dt_sample             # s, derivative filter time constant

        if valid and prev_valid:
            d_alpha = (theta_a - np.arctan(c_av/F_height)) / dt_sample
            d_beta  = (theta_b - np.arctan(c_bv/F_height)) / dt_sample

            d_alpha_f = d_alpha_f + dt_sample/tau * (d_alpha - d_alpha_f)
            d_beta_f  = d_beta_f  + dt_sample/tau * (d_beta  - d_beta_f)

        if GUIDANCE == 'pursuit':
            deflect = -(Kp*np.array([theta_a, theta_b]) + Kd*np.array([d_alpha_f, d_beta_f]))

        # ─────────────────────────────────────────
        # OUTER GUIDANCE LOOP -- gyro-compensated inertial LOS rate
        # ─────────────────────────────────────────
        # Sum the bearing delta and the gyro integral RAW, before any filtering.
        if valid and prev_valid:
            d_th_a = theta_a - np.arctan(c_av/F_height)
            d_th_b = theta_b - np.arctan(c_bv/F_height)
            lam_a, lam_b = cmp.los_rate(d_th_a, d_th_b, gyro_accum, dt_sample)
            LOS_rate_log.append((i*dt, lam_a, lam_b, d_th_a/dt_sample, d_th_b/dt_sample))
            wdes_y = NAV*lam_a
            wdes_x = -NAV*lam_b
        else:
            wdes_y = 0
            wdes_x = 0
        prev_valid = bool(valid)
        gyro_accum = np.zeros(3)

    # ─────────────────────────────────────────
    # Position step
    # ─────────────────────────────────────────
    P,rho,T_air = atm.return_atmo_state(Center_world[2])
    drag = PF.drag(rho,vel,S_REF,DRAG_COEFF)
    lift = PF.s_lift(rho,R.T@vel,S_REF,LIFT_COEFF)
    P_tar,rho_tar,T_tar = atm.return_atmo_state(tar[2])
    drag_tar = PF.drag(rho_tar,vel_tar,S_REF_TARGET,DRAG_COEFF_TARGET)

    # ─────────────────────────────────────────
    # Acquisition latch (the launch condition)
    # ─────────────────────────────────────────
    
    if valid:
        launched = True

    if launched:
        angle_cmd = np.linalg.norm(deflect)                 # rad, commanded gimbal angle
        angle = min(angle_cmd, theta_max)                   # rad, after the mechanical clamp
        u = deflect / angle_cmd if angle_cmd > 0 else np.zeros(2)   # gimbal azimuth
        F_ret = thrust_mag * (np.cos(angle)*np.array([0.0, 0.0, 1.0]) + np.sin(angle)*np.array([u[0], u[1], 0.0]))
    else:
        F_ret = np.zeros(3)   # engine unlit: no thrust, and hence no gimbal torque below

    if launched:
        mg = PF.g(m)
        a = (2*mg + R@F_ret + drag + R@lift) / (2*m)
        vel += a*dt
        Center_world += vel*dt
        m_fuel -= fuel_depletion_rate*dt
        m -= fuel_depletion_rate/2 * dt
        Rk = rock.Rocket(m, Bottom-Center, Top-Center, thrust,BODY_RADIUS)

    if m_fuel <= 0:
        print("no fuel")
        print(f"Simulation terminated at {i} steps ({i*dt}s)")
        break

    a_tar = (mg_tar+drag_tar)/m_tar
    vel_tar += a_tar*dt
    tar+=vel_tar*dt

    gyro_accum = gyro_accum + w_meas * dt   # the gyro runs whether or not we have launched

    # ─────────────────────────────────────────
    # Attitude step -- only once off the rail
    # ─────────────────────────────────────────
    # This must stay inside the launch gate. Left running while held, the rate
    # loop chases gyro noise, deflects the gimbal and torques the airframe, so
    # the vehicle rotates on the pad and lifts off with an attitude built from
    # noise. It also used to be what opened the old launch gate: any tilt makes
    # Bottom_world[2] = 2.25*(1 - cos(tilt)) > 0, which tripped at 1.3e-10 m on
    # step 1 of every run and defeated the hold entirely.
    if launched:
        T_control = np.cross(r,F_ret) - w*VISCOUS_DAMP      # N m, body frame
        w = rt.step(T_control,w,Rk.I_inv,Rk.I_principal,dt)
        omega = w*dt                                        # rad, rotation vector this step
        angle = np.linalg.norm(omega)                       # NOTE: reused, now the Rodrigues angle

        # ─────────────────────────────────────────
        # Rodrigues
        # ─────────────────────────────────────────
        if angle>ROT_EPS:
            axis = omega / angle
            K = la.skew(axis)
            R_delta=(np.eye(3)+np.sin(angle)*K+(1-np.cos(angle))*K@K)
            R=R@R_delta

        # Re-orthogonalise R; the det check rejects reflections.
        U_svd, _, Vt_svd = np.linalg.svd(R)
        if np.linalg.det(U_svd @ Vt_svd) < 0:
            U_svd[:, -1] = -U_svd[:, -1]
        R = U_svd @ Vt_svd

    # ─────────────────────────────────────────
    # Apply rotation
    # ─────────────────────────────────────────
    Top_world = Center_world + R @ (Top - Center)
    Bottom_world = Center_world + R @ (Bottom - Center)
    c_av = alpha
    c_bv = beta
    U.append((Top_world))
    D.append((Bottom_world))
    Q.append((R@(F_ret)))
    X.append(np.linalg.norm(tar-Center_world))
    TAR_log.append(tar.copy())
    LOS_view.append(R@la.quick_unit(np.array([alpha,beta,F_height])))
    LOS_log.append(np.array([alpha, beta]))     # raw, to match the unit-circle FOV plot


print("Closest approach: ", closest)
print("Final approach: ", np.linalg.norm(tar-Center_world))
print("Overshoot ratio: ", (closest-np.linalg.norm(tar-Center_world))/np.linalg.norm(tar-Center_world))

ani3d, ani2d = grp.animate(U, D, TAR_log, LOS_view, LOS_log, X, tar_static)
