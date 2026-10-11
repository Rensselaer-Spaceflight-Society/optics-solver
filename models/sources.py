# file for the light source classes -> they determine behavior through lenses, etc
import numpy as np
# what's neat here is you can use ray matrices to compute optical propagation
# shout out shout out rpi library

# generic light source class - we can also just define power as A U
# multi wavelength source by default -> for now might truncate in a basic class
class LightSource:
    def __init__(self):

        pass

    # helpers for the other element sources
    @staticmethod
    def translation(d, n):
        return np.array([[1.0, d / n], [0.0, 1.0]])

    @staticmethod
    def surface(R, n1, n2):
        P = 0.0 if np.isinf(R) else (n2 - n1) / R
        return np.array([[1.0, 0.0], [-P, 1.0]])

    def __init__(self, wavelength, power=1.0, n=1.0):
        self.wavelength = wavelength   # wavelength (might add multi-wavelength source later)
        self.power = power
        self.n = n                     # index of the medium it is currently in
        self.system = np.eye(2)        # cumulative M = M_k -> M_1 (set to np.eye(2) to reset)
 
    def apply(self, M, n_out=None):
        self.system = M @ self.system
        if n_out is not None:
            self.n = n_out
        return self
    
# okay first child class: paraxialsource -> transmits like a basic geom_optics source
# we can define apply to be a matrix multiplication
# we can define clip to be like based on radius/divergence angles as well
# might need something special for non-standard divergence angles (optic fiber acceptance angle)

class ParaxialSource(LightSource):
    profile = "flat"
 
    def __init__(self, wavelength, y=(0.0, 0.0), u=(0.0, 0.0), w=0.0, divergence=0.0, power=1.0, n=1.0):
        super().__init__(wavelength, power, n)
        # one row per transverse axis: [y, n*u]
        self.state = np.array([[y[0], n * u[0]], [y[1], n * u[1]]], dtype=float)
        self.marginal = np.array([w, n * divergence], dtype=float)   # [y_m, n*u_m]

    def apply(self, M, n_out=None):
        self.state = self.state @ M.T
        self.marginal = M @ self.marginal
        return super().apply(M, n_out)

    # all this below is matrix parameters from the kloos book
    @property
    def y(self):
        return self.state[:, 0]
 
    @property
    def nu(self):
        return self.state[:, 1]
 
    @property
    def u(self):
        return self.nu / self.n
    
    # spot radius
    @property
    def w(self):
        return abs(self.marginal[0])
    
    # add like tilted optics end up being a change of basis so they're aligned with a non-lens element