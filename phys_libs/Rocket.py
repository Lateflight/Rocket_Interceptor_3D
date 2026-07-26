import numpy as np
from phys_libs import Constantgen as cnst


class Rocket:
    
    """Creates a Rocket. Whoosh!
    m: mass
    Bottom: Bottom vector in body frame
    Top: Top vector in body frame
    """
    def __init__(self,m,Bottom,Top,thrust):

        
        self.m=m #mass of each section
        self.Bottom = Bottom
        self.Top = Top
        self.Center = (Bottom+Top)/2
        self.I_principal = cnst.mmoinertia(self.Top,self.Bottom,m) 
        self.I_inv=np.linalg.inv(self.I_principal)
        self.thrust = thrust
    def print(self):
        for k,v in self.__dict__.items():
            print(f"{str(k)} = {str(v)}")

