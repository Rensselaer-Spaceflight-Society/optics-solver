# Fiber optical simulation

> **New: native simulator.** [`native/`](native/README.md) contains **Fiber Light Lab**, a
> single-file desktop program (`FiberLightLab.exe`, no Python needed) that ray-traces the
> full goniometer bench: MI-152 → GH4001 → Thorlabs LB1092 lens → regolith-epoxy coupon →
> detector. It reproduces this Python model's fiber transport exactly and adds the lens, the
> mirror/instrument check, coupon reflectance models and the experiment matrix from the
> research plan. Download it from the repository's Releases, or build it with
> `cargo build --release` in `native/`. Equations and validation are in
> [native/MODEL_NATIVE.md](native/MODEL_NATIVE.md).

The Python versions below remain as the reference implementation.

## Desktop Python version

Run [fiber_simulation.py](fiber_simulation.py) for the **PySide6 desktop application
with OpenGL 3D views and Matplotlib charts**. It opens a normal light-mode window;
it does not start a browser or web server.

```powershell
cd C:\RSD_Optics_Research\fiber_light_lab
python fiber_simulation.py
```

Or double-click **Launch.cmd**, which now starts the desktop version using the
project's Python environment and installs missing desktop dependencies if needed.

For a new computer, Python 3.11+ is required:

```powershell
python -m pip install -r requirements-desktop.txt
python fiber_simulation.py
```

Keep `fiber_simulation.py`, `fast_graphics.py`, `physics.py`, `MODEL.md` and the `data/` folder together.
The application can be launched from any working directory. The ZIP distribution
contains these files and the Windows launcher; no Dash or Plotly installation is
required for the desktop application.

## Using the desktop app

- **Beam:** rotate the 3D geometry, inspect the sample irradiance and center
  profile, or select **Magnify fiber tip** for the last 8 mm of cable.
- **Fiber & spectrum:** view input/output spectra, attenuation, guided power
  along the cable and temporal coherence.
- **Interference:** choose the halogen mean, coherent comparison, or averaged
  independent patterns. The signed contribution separates constructive and
  destructive interference from the incoherent reference.
- The four main controls are cable length, sample distance, input power and
  launch NA. **Additional settings** reveals spot radius, offset, extra loss
  and wavelength band. Press Enter or leave a numeric field to apply edits.
- Drag the 3D view to rotate; scroll to zoom; middle-drag to pan. **Reset view**
  restores its camera. The divider between 3D and charts can be dragged.
- **Save data…** writes settings and CSVs to a ZIP archive. **Save 3D image…**
  saves the current camera view as PNG. The disk icon below the 2D charts saves
  those charts as PNG, PDF, SVG or another Matplotlib format.
- **Model notes** opens the equations, sources and limitations inside the app.

Defaults remain **10 inches of cable**, **2.5 inches to the sample** and an
**assumed 1 mW optical input**. Halogen illumination and the illustrative coherent
comparison remain separate. This version models incident illumination, not
regolith scattering, absorption or heating.

Calculations run in one background worker; rapid edits replace pending work.
Only the selected tab is updated, and unchanged tabs retain their plots and camera.
Bounded caches reuse calculations. The halogen mean skips complex-field allocation.
OpenGL renders the full 161 × 161 interference surface, with smooth vertex colors
and four-sample edge antialiasing where supported. The 2D image interpolation is
display-only: exported values and numerical resolution are unchanged. Rotation
updates only the 3D camera, without rerunning physics or redrawing Matplotlib.

If OpenGL or its dependencies are unavailable, the application falls back to
Matplotlib 3D. To request that renderer explicitly:

```powershell
python fiber_simulation.py --software
```

The status bar identifies the active renderer. Software mode uses a reduced
display mesh and rotates more slowly. On the tested Intel Iris Xe GPU, warm 3D
draws took about 4–8 ms, compared with 277–302 ms for the previous Matplotlib
rotation. These are rendering timings, not an end-to-end frame-rate guarantee;
display refresh, drivers, window size and hardware matter. See
`verification/performance_before.json`, `verification/performance_after.json`
and `verification/benchmark_desktop.py` for the measurements and repeatable check.

Desktop verification:

```powershell
python -m pytest tests/test_desktop.py -q
```

## Original browser version (workspace only)

The earlier browser application remains available as `app.py`. Use
**Launch browser version.cmd** to open it; see README_browser.md for its instructions. Its separate
dependencies are in `requirements.txt`. The desktop app is the default launcher.

