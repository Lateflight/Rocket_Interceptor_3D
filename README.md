# Bearing-Only Interceptor Simulation

A 6-DOF flight simulation of a thrust-vectored interceptor guided onto a
maneuvering target using **bearing-only** measurements, a body-fixed seeker,
and a strapdown rate gyro. No range, closing velocity, or target state is
available to the guidance law.

The vehicle has no fins and no reaction control. Every steering input is a
torque produced by gimballing the engine about the center of mass: the
airframe rotates first and translates as a consequence.

---

## Quick start

```bash
python physics_engine.py
```

Requires `numpy`, `scipy`, and `matplotlib`. Run from the repository root. Modules are imported as `from phys_libs import ...`.

The run prints a short summary and then opens two live animations (see [Output](#output)). Close both windows to exit.

---

## Repository layout

```
physics_engine.py       entry point: parameters, the integration loop, logging
phys_libs/
  computer.py           seeker, IMU, LOS-rate estimator  (the "flight computer")
  PiForce.py            gravity, quadratic drag, sin(AoA) normal force
  Rocket.py             vehicle object; holds I_principal and its inverse
  Constantgen.py        moment of inertia for the two-segment rod
  Rotator_module.py     RK4 step of Euler's rigid-body equations
  atmosphere.py         ICAO/ISA temperature, pressure, density to 47 km
  linalg_utils.py       skew, Rodrigues helpers, unit vectors, LQR solve
  grapher.py            matplotlib animation of a completed engagement
  PiTorque.py          quick torque generator, kept for simplicity in case of many forces
```

## The scenario conditions\*

- range below `HIT_RADIUS` = 5 m — intercept
- either body reaching `GROUND_ALT` = 0 m
- interceptor fuel exhaustion

*\*Note: scenario paramters can be changed, for more information, see [Key parameters](#key-parameters)]*

---

## How it works

### Physics

- **Attitude** — Euler's equations integrated with RK4 (`Rotator_module.step`),
  then the rotation vector `w·dt` applied to `R` via Rodrigues. `R` is
  re-orthogonalised each step by SVD, with a determinant check that rejects
  reflections.
- **Translation** — semi-implicit Euler on thrust + gravity + drag + normal
  force, all resolved into the world frame.
- **Mass** — the vehicle is modelled as two sections of mass `m` each (total
  `2m`). Fuel burns at 0.04 kg/s and the inertia tensor is rebuilt each step
  from the current mass, so the vehicle gets more responsive as it empties.
- **Aerodynamics** — quadratic drag, plus a normal force whose magnitude scales
  with `sin(AoA)`. Both use the ISA atmosphere at each body's own altitude.

### Guidance

The seeker is bolted to the airframe, so a target drifting across the field of
view is ambiguous — the target may be moving, or the vehicle may be rotating
underneath the seeker:

```
dθ/dt  =  λ̇  −  ω
```

One equation with two unknowns. The gyro supplies `ω`, which makes the **inertial
LOS rate `λ̇`** observable (`computer.los_rate`). That single measurement is
what the whole design rests on: driving `λ̇` to zero is a collision course,
and it needs no range.

Two nested loops:

| Loop | Rate | Input | Output |
|---|---|---|---|
| Outer (guidance) | 50 Hz | bearing delta + gyro integral over the same interval | commanded body rate `w_des = N·λ̇` |
| Inner (rate) | 100 Hz | gyro | gimbal deflection `= K_w·(w_des − w_meas)` |

The commanded deflection is clamped to the mechanical gimbal limit, then the thrust vector is rebuilt geometrically from that angle, so deflection stays proportional below the limit instead of saturating. 

**Critical detail:** the seeker term and the gyro term are summed raw, before any filtering. Filtering one but not the other leaves a residual of `ω − filtered(ω)`, which injects the vehicle's own rotation straight into the guidance command.

For the tick of measurements, 20 ms (depending on the scenario, generalizes to `2·dt` [see [Known Limitations](#known-limitations)]) are spent with no guidance in the PN law. Removing this break causes the rocket to use oversaturated values in the derivative term, which can command up to 61 rad/s for trajectories at the edge of the seeker's cone.

Set `GUIDANCE = 'pursuit'` for the earlier bearing-PD law, kept for comparison.

---

## Key parameters

All in `physics_engine.py`:

| Symbol | Value | Meaning |
|---|---|---|
| `dt` | 0.01 s | integration step (100 Hz) |
| `SEEKER_LAG_STEPS` | 2 | seeker refreshes every 2nd step (50 Hz) |
| `F_height`, `F_len` | 4, 1 | focal-plane standoff and radius → FOV half-angle 14.0° |
| `theta_max` | 10° | gimbal mechanical limit |
| `thrust_mag` | 4000 N | constant while fuel remains |
| `VISCOUS_DAMP` | 0.2 | artificial body-rate damping in the torque sum |
| `NAV` | 5.0 | navigation constant `N` |
| `KW` | 0.10 | inner rate-loop gain |
| `gyro_bias` | ~0.2 °/s | turn-on bias, calibrated on the pad |
| `gyro_noise` | 0.0015 rad/s | per-sample white noise |

The FOV half-angle must stay above the lead angle the geometry demands (~12°
for this scenario), or the seeker loses the target while turning into the
intercept.

Gyro noise is drawn from a seeded generator (`GYRO_SEED = 12345`), so runs are
reproducible.

---

## Output

Console:

```
gyro bias  true [...]
     estimated [...]   residual [...] rad/s
Closest approach:  ...
Final approach:  ...
Overshoot ratio:  ...
```

Two matplotlib animations:

1. **3D world view** — interceptor nose and tail traces, target track, LOS
   line, and a boresight quiver.
2. **Seeker view** — the raw `[alpha, beta]` bearing inside the FOV circle,
   with a live range readout. With `F_len = 1` the FOV gate reduces to the unit
   disc, so bearings are plotted unscaled.

Per-step buffers (`U`, `D`, `X`, `TAR_log`, `LOS_log`, `LOS_rate_log`, …) are
kept in memory after the loop and can be dumped to CSV for analysis.
`LOS_rate_log` in particular holds both the gyro-compensated LOS rate and the
raw seeker-only rate, which is the quickest way to see the compensation working.

---

## Known limitations

- No `Vc` gain scheduling — fixed `N` means the effective loop gain drifts over
  the engagement.
- No target-acceleration compensation; augmented PN would need time-to-go,
  which needs range.
- Lock loss is handled by *holding* the last bearing, not by reacquiring.
- One timestep of attitude lag in the translational step.
- The airframe has no aerodynamic stability or damping model; `VISCOUS_DAMP` is
  a stand-in.
- `atmosphere.return_atmo_state` has no floor at `z = 0`; callers must stop at
  the ground themselves.
