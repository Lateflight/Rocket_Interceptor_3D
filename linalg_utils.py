import numpy as np
from scipy.linalg import solve_continuous_are
def quick_unit(V):
    return V/np.linalg.norm(V)


def project(v, u):
    return (np.dot(v, u) / np.dot(u, u)) * u

def skew(w):

    s = np.array([[0, -w[2],w[1]],
              [w[2], 0,-w[0]],
              [-w[1],w[0],0]])
    
    return s

def rotation_matrix(axis, theta):


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

def dLQR(A,B,Q,R):

    P = solve_continuous_are(A,B,Q,R)

    k = np.linalg.inv(B.T@P @ B + R) @ (B.T @ P @ A)
    return k,P

def angle_between(v1, v2):
    # v1 = (x1, y1), v2 = (x2, y2)
    dot = v1[0]*v2[0] + v1[1]*v2[1]      # Dot product
    det = v1[0]*v2[1] - v1[1]*v2[0]      # 2D Cross product / Determinant
    return np.arctan2(det, dot)          # Returns angle in radians

def normalize(x):
    if x == 0:
        return 0
    while abs(x) < 1:
        x *= 10
    while abs(x) >= 10:
        x /= 10
    return x