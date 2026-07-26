import numpy as np
import linalg_utils as la
def torque(forces,radius,frames,return_frame=np.eye(3)):

    """Returns torques given forces and the frame they are in. Default return frame is the identiy
    forces: A tuple of force values as vectors
    radius: a tuple of radii
    frames: A tuple of np.array()s representing the frame->world
    return_frame: the frame in which the torque acts"""


    standard_frames = []
    conversion_frame=return_frame.T
    for i in frames:
        
        standard_frames.append(np.dot(conversion_frame,i))

    f_frame=[]

    for i in range(len(standard_frames)):

        f_frame.append(standard_frames[i]@forces[i])

    torques = np.array([0.0,0.0,0.0])
    radius_conv = []

    for i in radius:

        radius_conv.append(i)

    for i in range(len(f_frame)):
        torques+=np.cross(radius_conv[i],f_frame[i])

    return torques






