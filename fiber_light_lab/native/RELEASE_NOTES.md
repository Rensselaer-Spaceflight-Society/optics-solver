**Fiber Light Lab** is a single-file simulator for the regolith-epoxy goniometer bench:
MI-152 halogen → GH4001 fiber → Thorlabs LB1092 lens → coupon → detector.

### Download

- **FiberLightLab.exe**: Windows 10/11 x64. Download and double-click; nothing to install.
- **FiberLightLab-windows-x64.zip**: the app, the `fiberlab-cli.exe` batch runner and the model notes.
- **FiberLightLab-linux-x64.tar.gz**: the same programs for Linux.

The executables are not code-signed, so Windows SmartScreen may show "Windows protected your
PC". Choose *More info → Run anyway*. If the window does not open on a machine without
DirectX 12 or Vulkan, start it with `FiberLightLab.exe --glow` to use OpenGL.

### What it does

- Reproduces the Python fiber model exactly (0.7843 mW default output; parity tests to 1e-12).
- Ray-traces the LB1092 (N-BK7 dispersion, exact spherical surfaces, Fresnel losses, aperture
  vignetting), the tilted coupon and a detector of finite size on the goniometer arm.
- Mirror / instrument-function check, diffuse reference, and a composite reflectance model
  (resin surface + body + Hapke exposed grains from Curtis et al. 2024) for LHS-1 and LMS-1
  at any regolith loading, as-cast or abraded.
- Near-specular vs away-from-specular metrics and the full experiment matrix from the
  research plan (0/30/45/60 wt% × as-cast/abraded × 30°/60°).
- 3D bench view, measured-scan overlay, CSV/JSON export, settings files, and a
  command-line runner.

Geometry, detector and coupon values are editable placeholders. Enter your measured bench
before drawing conclusions. See *Model notes* in the app for equations, sources,
verification and limitations.
