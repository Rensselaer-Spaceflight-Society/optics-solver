# Verification — 3 October 2026

- 13 pytest cases passed after the final physics/figure/export changes.
- Numerical results and measured calculation times are in `numerical_checks.json`.
- Default irradiance integral matched output power to 2.2e-16 relative error.
- Peak irradiance changed 0.0127% between 257² and 513² grids.
- Second radial moment agreed with an independent 262,144-ray Sobol calculation
  to 0.0201% at an 8 mm plane.
- Coherent illustrative contrast was 0.994; averaging 64 independent patterns
  reduced it to 0.121. The incoherent limit has exactly zero cross terms.
- Browser review covered output maps, the fiber cutaway, the coherent signed
  interference view, power changes, reset and download. Doubling optical input
  doubled output power and peak irradiance while preserving footprint size.
- Final browser error log was empty. The export button produced a ZIP in the
  user's Downloads folder; archive integrity and settings were checked.
- `output_field.jpg` and `interference.jpg` are screenshots from the running app.
- `default_simulation.zip` is a reproducible default export, including data,
  parameters, reference curves, model notes and an offline interactive 3D view.

These are checks of numerical consistency and software behavior, not validation
against measurements of the physical apparatus. Timing depends on machine load;
observed full figure-generation callbacks were about 0.7–1.1 s. Browser orbiting
is local and does not call the Python calculations. No frame-rate claim was tested.


## Native desktop optimization (4 October 2026)

The native PySide6 application now uses pyqtgraph/OpenGL for 3D, with Matplotlib
for 2D charts and a software 3D fallback. The original browser screenshots above
are not screenshots of the desktop application.

- `desktop_window.png` and `desktop_interference.png`: current desktop previews.
- `performance_before.json`: CPU baseline, median warm calculations and Agg
  rendering at a 1200 × 720 canvas; rotation included the surrounding 2D charts.
- `performance_after.json`: current local GPU timings and hardware details.
  GPU rendering uses an actual forced paint followed by `glFinish`; it isolates
  the 3D pane. These compare the application's old and new drawing paths, not
  identical primitive counts, resolutions or end-to-end displayed FPS.
- `benchmark_desktop.py`: repeatable native GPU smoke test and benchmark. Needs
  an interactive desktop and the desktop requirements. It opens a temporary
  test window, checks parameter changes and cached tab revisits, exercises all
  interference modes and zero-power / zero-distance cases, then closes it.
- The numerical and desktop suite passed 17 tests, covering energy conservation,
  reference transport, grid convergence, independent ray tracing, coherence,
  cache correctness under changed parameters, direct CSV values and software
  rendering. Four-sample framebuffer antialiasing was verified on Intel Iris Xe.

Run `python verification/benchmark_desktop.py` from the project directory.
Camera movement performs no model calculation. The full 161² interference grid
is rendered. First-time imports, new wave-basis construction, export compression,
monitor refresh and desktop composition are outside the warm draw measurements.
