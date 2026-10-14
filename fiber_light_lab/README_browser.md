# Fiber Light Lab

Interactive Python simulation for the MI-152 halogen illuminator and GH4001
plastic optical fiber. The default setup is your **10-inch cable** and
**2.5-inch tip-to-sample distance**.

## Open it

Double-click **Launch.cmd**. The local app opens at
[127.0.0.1:8057](http://127.0.0.1:8057). Keep the Python process running while
using it. This installation already has its own environment in `.venv`.

Or, in a terminal in this folder:

```powershell
.\.venv\Scripts\python.exe app.py
```

For a fresh installation (Python 3.11+):

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe app.py
```

The launcher automatically performs those steps if `.venv` is absent. It binds
only to localhost. `Stop.ps1` stops only the app process recorded by this launcher.

## Explore

- **Output field:** orbit the 3D fiber/output cone, zoom into the tip, inspect
  the sample irradiance and line profile, or switch to logarithmic intensity.
- **Inside the fiber:** inspect total internal reflection, launch acceptance,
  spectral attenuation and guided power along the cable.
- **Interference:** start with the halogen incoherent mean; switch to the
  coherent comparison to see speckle and signed constructive/destructive terms.
  Average independent patterns to see contrast fall. The phase slider changes
  local propagation phase, not the physical sample location.
- **Model & provenance:** read the assumptions and source status.
- **Export:** saves settings, irradiance and spectrum CSVs, the illustrative wave
  patch, documentation, and an offline interactive 3D HTML view as one ZIP.

All geometry and optical input controls are adjustable. Launch NA smaller than
0.50 underfills the fiber and narrows the beam. Input spot size and misalignment
change accepted power. Additional loss accounts for a measured loss you supply.

**The default 1 mW is an assumption, not a measurement.** The lamp's 150 W rating
is electrical. Spectral curves are approximate; absolute irradiance needs
calibration. The halogen and coherent comparison models are deliberately
distinguished. The model ends at the incident sample plane, before regolith
scattering, absorption or heating. See [MODEL.md](MODEL.md) for equations and limits.

## Files

- `app.py`: Python Dash interface and export.
- `physics.py`: reusable calculation functions, independent of the interface.
- `figures.py`: Plotly scenes and plots.
- `data/`: supplied references and transparent, editable spectral CSVs.
- `tests/test_physics.py`: conservation, analytic-limit, interference and export tests.
- `verification/`: checks and a default simulation export.

Run validation:

```powershell
.\.venv\Scripts\python.exe -m pytest tests -q
```

Reusable calculation example:

```python
from physics import Settings, transport, footprint
s = Settings(length_m=0.254, distance_mm=63.5, incident_mw=1.0)
print(transport(s)['output_mw'])
axis_mm, irradiance_mw_mm2, envelope_radius_mm = footprint(s, n=513)
```
