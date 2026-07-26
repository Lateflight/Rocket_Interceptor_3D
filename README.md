# Rocket Interceptor in 6 DOF Simulator
NOTE: PROTOTYPE, SPLIT INTO GOOD .md LATER

## Specifics of Design

Goal: interception within <= 5m on some trajectory

### On board sensors:
 - Circular seeker - Tracks the target within a specified POV as described by two variables alpha and beta (see variable definition for specifics) 
 - On-board gyroscope - Tracks the only current rotational velocity of the rocket interceptor with noise (see physics_engine.py to seed)

### Control authority (see physics.py):

 - Fixed thrust engine
 - Engine gimbal limited at 10 degrees 

## Simulation Basics

### The simulation rests on a few assumptions:

 - Simple drag: 1/2\* rho\* v^2\* S\*|v| \*C_d
 - Simple lift: 1/2\* rho\* v^2\* S\* cross(|v|,z) \*C_l applying force at CoM
 - ISA Standard Atmosphere: https://www.engineeringtoolbox.com/international-standard-atmosphere-d_985.html
 - Terrain curvature does not matter (due to all intercepts happening due to on-board computation)
 - Thrust only torque control
 - Target is point mass
 - Rigid body dynamics for the rocket
 - Linear viscous torque from atmosphere
 - For moments of inertia: Ix = Iy

 ### There are also two coordinate frames:
 - Body frame: tied directly to the rocket
 - World frame: defined initially
 Conversion happens through Rodrigues's rotation formula (see physics_engine.py , planned to move into linalg_utils.py)

 ## Planned Features
 ### Acutators
 - Time dependent response
 - Thrust change due to height  \*
 ### Sensors and Noise
 - Seeker noise
 - Kalman filter on Alpha/Beta instead of derivative filter (outer loop) \*
 ### Atmosphere and Flight
 - Lift torque
 - Wind  \*

  *NOTE: none of these are guaranteed, as the scope of the project mostly was to demonstrate control law understanding. Everything marked with a \* is a consideration\blind guesses.*
