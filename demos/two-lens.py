import numpy as np
from models import *
import matplotlib.pyplot as plt


def at(z):
    return (
        CoordinateVector(3, 0.0, 0.0, z),
        UnOrderedRotation(3, 0.0, 0.0, 0.0)
    )


# 600 nm paraxial beam, 5 mm radius (u=0 -> parallel rays)
src = ParaxialSource(
    wavelength=600e-9,
    w=5e-3,
    divergence=0.0
)

# Create a stage with one arm.
# The arm's optical axis is +z; its transverse frame is generated automatically.
stage = Stage()
stage.add_arm(
    arm_id=0,
    optical_axis=CoordinateVector(3, 0.0, 0.0, 1.0),
    source=src,
    start=0.0
)

# Biconvex lens, R = 10 cm each side, n = 1.5 -> f ~ 10 cm
# Front vertex at z = 0.10 m
stage.add_elem(
    0,
    CircularLens(
        position=CoordinateVector(3, 0.0, 0.0, 0.10),
        thickness=5e-3,
        radius=1e-2,
        R1=0.1,
        R2=-0.1,
        index=1.5,
        rotation=UnOrderedRotation(3, 0.0, 0.0, 0.0)
    )
)

stage.add_elem(
    0,
    CircularLens(
        position=CoordinateVector(3, 0.0, 0.0, 0.50),
        thickness=5e-3,
        radius=1e-2,
        R1=0.5,
        R2=-0.5,
        index=1.5,
        rotation=UnOrderedRotation(3, 0.0, 0.0, 0.0)
    )
)

# Spot-size readers sweeping through the focus
zs = np.linspace(0.0, 1.0, 100)
probes = [
    SpotDiagnostic(*at(z), name=f"z={z:.3f}")
    for z in zs
]

for probe in probes:
    stage.add_elem(0, probe)

# Run the stage and collect diagnostic readings.
# Each diagnostic may be hit by a beam once; readings are keyed by diagnostic.
results = stage.run()

spot_radii = np.array([results[p][0]["w"] for p in probes])

print(results)

plt.plot(zs * 1e3, spot_radii * 1e3, "o-")
plt.xlabel("z (mm)")
plt.ylabel("spot radius (mm)")
plt.title("Two-lens focii")
plt.grid(True)
plt.show()