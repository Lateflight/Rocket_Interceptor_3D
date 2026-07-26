import numpy as np
def quick_unit(V, guard = 1e-9):
    """Returns a unit vector of vector v.
    v: a vector.
    guard: smallest accepted value of norm(v)
    
    NOTE: If guard fails, will return the 0 vector"""
    norm = np.linalg.norm(V)
    if norm > guard:
        return V/norm
    else:
        return np.zeros(3)

def project(v, u, guard = 1e-9):
    """Returns the projection of vector v onto vector u.
    v: vector (projected)
    u: vector (projected onto)
    guard: smallest accepted value of |u|^2.
    """
    denom = np.dot(u, u)

    if denom < guard:
        raise ValueError("Cannot project onto the zero vector.")

    return (np.dot(v, u) / denom) * u
def skew(w):

    """Returns the 3 dimensional skew-symmetric matrix of a given vector.
    w: a vector"""

    s = np.array([[0, -w[2],w[1]],
              [w[2], 0,-w[0]],
              [-w[1],w[0],0]])
    
    return s

def rotation_matrix(axis, theta):

    """Returns a rotation matrix given an axis and an angle (0-2pi)
    Axis: a vector about which rotation happens
    Theta: angle of rotation (radians)"""
    axis = axis / np.linalg.norm(axis)
    x, y, z = axis
    c = np.cos(theta)
    s = np.sin(theta)
    C = 1 - c

    return np.array([
        [c + x*x*C,     x*y*C - z*s, x*z*C + y*s],
        [y*x*C + z*s,   c + y*y*C,   y*z*C - x*s],
        [z*x*C - y*s,   z*y*C + x*s, c + z*z*C]
    ])

def angle_between(v1, v2):

    """Returns the angle between two vectors v1 and v2, in radians.
    v1: Vector defined as 0 radians.
    v2: Vector to which the measurement goes."""
    dot = v1[0]*v2[0] + v1[1]*v2[1]      # Dot product
    det = v1[0]*v2[1] - v1[1]*v2[0]      # 2D Cross product / Determinant
    return np.arctan2(det, dot)          # Returns angle in radians

def normalize(x):

    """Returns the part of a order-of-magnitude-independent number.
    x: a float or integer
    """
    if x == 0:
        return 0
    while abs(x) < 1:
        x *= 10
    while abs(x) >= 10:
        x /= 10
    return x