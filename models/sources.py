# file for the light source classes -> they determine behavior through lenses, etc
import copy
import numpy as np
# what's neat here is you can use ray matrices to compute optical propagation
# shout out shout out rpi library
from uuid import uuid4
# generic light source class - we can also define power as A U
# multi wavelength source by default -> for now might truncate in a basic class
from models import *


class LightSource:
    # helpers for the other element sources
    @staticmethod
    def translation(d, n):
        return np.array([[1.0, d / n], [0.0, 1.0]])

    @staticmethod
    def surface(R, n1, n2):
        P = 0.0 if np.isinf(R) else (n2 - n1) / R
        return np.array([[1.0, 0.0], [-P, 1.0]])

    def __init__(self, wavelength, power=1.0, n=1.0, options={
        "rotation": CoordinateVector(0, 0, 0)
    }):
        self.wavelength = wavelength
        self.power = power
        self.n = n
        self.system = np.eye(2)  # cumulative M = M_k -> M_1 (set to np.eye(2) to reset)
        self.id = uuid4()
        self.options = options

    def apply(self, M, n_out=None):
        self.system = M @ self.system
        if n_out is not None:
            self.n = n_out
        return self

    def transform_transverse(self, basis):
        # move frames of the transverse beam
        return self

    def fork(self, fraction, basis=None):
        # create a new branch
        child = copy.deepcopy(self)
        child.id = uuid4()
        child.power = self.power * fraction
        if basis is not None:
            child.transform_transverse(basis)
        return child


# okay first child class: paraxialsource -> transmits like a basic geom_optics source
# we can define apply to be a matrix multiplication
# we can define clip to be like based on radius/divergence angles as well
# might need something special for non-standard divergence angles (optic fiber acceptance angle)
class ParaxialSource(LightSource):
    profile = "flat"

    def __init__(self, wavelength, y=(0.0, 0.0), u=(0.0, 0.0), w=0.0, divergence=0.0, power=1.0, n=1.0, options={}):
        super().__init__(wavelength, power, n)
        # one row per transverse axis: [y, n*u]
        self.state = np.array([[y[0], n * u[0]], [y[1], n * u[1]]], dtype=float)
        self.marginal = np.array([w, n * divergence], dtype=float)  # [y_m, n*u_m]

    def apply(self, M, n_out=None):
        self.state = self.state @ M.T
        self.marginal = M @ self.marginal
        return super().apply(M, n_out)

    def transform_transverse(self, basis):
        basis = np.asarray(basis, dtype=float)
        if basis.shape != (2, 2):
            raise ValueError(f"transverse transform must be 2x2, got {basis.shape}")
        # The rows of state are transverse axes; rotate both position and angle.
        self.state = basis @ self.state
        return self

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



# okay essentially extends all of the properties of a paraxial source to
# a keyed dictionary where an element is a wavelength
class ParaxialMultiwavelengthSource(LightSource):
    def __init__(self, wavelengths=None, powers=None, y=(0.0, 0.0), u=(0.0, 0.0),
                 w=0.0, divergence=0.0, n=1.0, options=None):
        wavelengths = [540e-9] if wavelengths is None else list(wavelengths)
        powers = [1.0] * len(wavelengths) if powers is None else list(powers)
        if len(powers) != len(wavelengths):
            raise ValueError("wavelengths and powers must have the same length")
        if len(set(wavelengths)) != len(wavelengths):
            raise ValueError("duplicate wavelengths")
        if options is None:
            options = {"rotation": CoordinateVector(0, 0, 0)}

        super().__init__(tuple(wavelengths), sum(powers), n, options)

        self.children = {
            wv: ParaxialSource(wavelength=wv, y=y, u=u, w=w, divergence=divergence,
                               power=p, n=n, options=options)
            for wv, p in zip(wavelengths, powers)
        }

    def apply(self, M, n_out=None):
        for wv, child in self.children.items():
            M_wv = M[wv] if isinstance(M, dict) else M
            n_wv = n_out[wv] if isinstance(n_out, dict) else n_out
            child.apply(M_wv, n_wv)
        if n_out is not None and not isinstance(n_out, dict):
            self.n = n_out
        return self

    def transform_transverse(self, basis):
        for child in self.children.values():
            child.transform_transverse(basis)
        return self

    def fork(self, fraction, basis=None):
        branch = copy.deepcopy(self)
        branch.id = uuid4()
        branch.power = self.power * fraction
        for src in branch.children.values():
            src.id = uuid4()
            src.power *= fraction
        if basis is not None:
            branch.transform_transverse(basis)
        return branch

    @property
    def powers(self):
        return {wv: c.power for wv, c in self.children.items()}

    @property
    def y(self):
        return {wv: c.y for wv, c in self.children.items()}

    @property
    def nu(self):
        return {wv: c.nu for wv, c in self.children.items()}

    @property
    def u(self):
        return {wv: c.u for wv, c in self.children.items()}

    @property
    def w(self):
        return {wv: c.w for wv, c in self.children.items()}