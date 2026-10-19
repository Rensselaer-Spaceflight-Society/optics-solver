# file meant to carry your sensors, etc. stuff
import numpy as np
from models.elements import Element
from uuid import uuid4

class Diagnostic(Element):
    def __init__(self, position, rotation, name):
        super().__init__(position, rotation)
        self.data = []
        self.name = name
        self.id = uuid4()

    def matrix(self, src):
        return np.eye(2), src.n

    def transfer(self, src):
        self.data.append(self.record(src))
        return super().transfer(src)

    def record(self, src):
        raise NotImplementedError

class SpotDiagnostic(Diagnostic):
    def __init__(self, position, rotation, name):
        super().__init__(position, rotation, name)

    def record(self, src):
        return {"w": float(src.w), "y": src.y.copy(), "power": src.power}