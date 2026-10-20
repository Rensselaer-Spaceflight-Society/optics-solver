from models.space import CoordinateVector, UnOrderedRotation
from models.sources import LightSource
from typing import Callable, Union
import numpy as np
from uuid import uuid4
#from archive.hull import LensHull

#ELEMENT CLASSES
class Element:
    ends_arm = False# True -> no "through" beam (a mirror)
    def __init__(
            self,
            position: CoordinateVector,
            rotation: UnOrderedRotation,
            # param_map,
            *params
        ):
        self.id = uuid4()
        self.position = position
        self.rotation = rotation
        # self.hull = param_map(*params)

    def matrix(self, src: LightSource):# subclasses define this that's the name of the game
        raise NotImplementedError

    def transmission(self, src:LightSource):# generally no loss
        return 1.0


    def transfer(self, src: LightSource):
        src.power *= self.transmission(src)
        src.apply(*self.matrix(src))
        return src

# params are : center, rotation, thickness, radius, R1, R2 (where R1,2 = radii of curvature)
class CircularLens(Element):
    def __init__(self, position, rotation, thickness, radius, R1, R2, index: Union[float, Callable[[float], float]] = 1.0): # default lens doesn't do anything
        super().__init__(position, rotation)
        for R in (R1, R2):
            if abs(R) < radius:
                raise ValueError(f"|R|={abs(R)} must be >= lens radius {radius}")
        self.thickness, self.radius, self.R1, self.R2 = thickness, radius, R1, R2
        self.set_index(index)

    def set_index(self, index):
        self.index = index if callable(index) else (lambda wavelength: index)

    def matrix(self, src: LightSource):
        n_g = self.index(src.wavelength)
        M = (src.surface(self.R2, n_g, src.n)
             @src.translation(self.thickness, n_g)
             @src.surface(self.R1, src.n, n_g))
        return M, src.n

    def transmission(self, src: LightSource):
        return self._circular_transmission(src) 

    def _circular_transmission(self, src: LightSource):
        a = self.radius
        if src.w == 0:
            return 1.0 if np.hypot(*src.y) <= a else 0.0
        if src.profile == "gaussian":
            return 1.0 - np.exp(-2 * a**2 / src.w**2)
        return min(1.0, a**2 / src.w**2)

# params are center, [lengths in each axis]
# where the order of the lengths in each axis corresponds to the order of the coords in the position
# in models/elements.py
class MirrorElement(Element):
    def __init__(self, position, rotation, reflectivity=1.0):
        super().__init__(position, rotation)
        self.reflectivity = reflectivity

    def matrix(self, src):
        return np.eye(2), src.n

    def transmission(self, src):
        return src.power * self.reflectivity

    def normal(self):
        return np.array(self.rotation.transform(CoordinateVector(3, 0.0, 0.0, 1.0)).coords)