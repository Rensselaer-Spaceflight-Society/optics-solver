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
    w=1e-3,
    divergence=0.0
)

# Create a stage with one arm.
# The arm's optical axis is +z; its transverse frame is generated automatically
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
    CircularLens( # symmetric biconvex lens (LB1406, idk if its same model as ours)
        position=CoordinateVector(3, 0.0, 0.0, 50e-3 - 24e-3), # make them 2x focal length apart
        thickness=2.3e-3,
        radius=3e-3,
        R1=12e-3,
        R2=-12e-3,
        index=1.5, # give or take, need to add a index curve later maybe
        rotation=UnOrderedRotation(3, 0.0, 0.0, 0.0)
    )
)

stage.add_elem(
    0,
    CircularLens(
        position=CoordinateVector(3, 0.0, 0.0, 50e-3),
        thickness=2.3e-3,
        radius=3e-3,
        R1=12e-3,
        R2=-12e-3,
        index=1.5,
        rotation=UnOrderedRotation(3, 0.0, 0.0, 0.0)
    )
)

# Spot-size readers sweeping through the focus
zs = np.linspace(0.0, 100e-3, 100)
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


# annotation for lens positions
plt.axvline(x=50-24, color='black', linestyle='--', alpha=0.7) 
plt.annotate(
    text='Lens 1',
    textcoords="offset points",
    xy=(50-24, 0.8),
    xytext=(-40, 30),
    fontsize=10,
    color='black'
)
plt.axvline(x=50, color='black', linestyle='--', alpha=0.7) 
plt.annotate(
    text='Lens 2',
    textcoords="offset points",
    xy=(50, 0.8),
    xytext=(40, 30),
    fontsize=10,
    color='black'
)


plt.xlabel("z (mm)")
plt.ylabel("spot radius (mm)")
plt.title("Two-lens focii")
plt.grid(True)
plt.show()