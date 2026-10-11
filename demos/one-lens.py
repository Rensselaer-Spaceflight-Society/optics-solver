import numpy as np
from models import *
import matplotlib.pyplot as plt

def at(z):
    return CoordinateVector(3, 0.0, 0.0, z), UnOrderedRotation(3, 0.0, 0.0, 0.0)

# 600 nm paraxial beam, 5 mm radius (u=0 -> parallel rays)
src = ParaxialSource(wavelength=600e-9, w=5e-3, divergence=0.0)

stage = Stage(arms=1)
stage.set_arm(0, {"length": 0.5, "source": src})

# biconvex lens, R = 10 cm each side, n = 1.5 -> f ~ 10 cm; front vertex at z = 0.10 m
stage.add_obj(0, CircularLens(*at(0.10), thickness=5e-3, radius=1e-2, R1=0.1, R2=-0.1, index=1.5))

# plus spot size reader sweeping through the focus
zs = np.linspace(0.12, 0.40, 15)
probes = [SpotDiagnostic(*at(z), name=f"z={z:.3f}") for z in zs]
for p in probes:
    stage.add_obj(0, p)

results = stage.run()[0]

print(results)

plt.plot(zs * 1e3, [results[p]["w"] * 1e3 for p in probes], "o-")
plt.xlabel("z (mm)"); plt.ylabel("spot radius (mm)"); plt.title("One-lens focus")
plt.show()