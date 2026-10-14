# Fiber Light Lab — model definition

## Scope

Python spectral/radiometric transport from an assumed incident lamp field to a
straight GH4001 fiber and a perpendicular plane in air. A separate local scalar
wave explorer demonstrates interference. The simulation ends before the light
interacts with regolith or composite material. It is not a calibrated prediction
of the user's apparatus, a thermal calculation, or a full vector Maxwell solver.

Defaults requested by the user: cable length **10 in = 0.254 m**, tip-to-sample
distance **2.5 in = 63.5 mm**. Default band-limited optical input **1 mW** is an
illustrative low-power scenario, not a measured or inferred lamp output. The
150 W electrical lamp rating cannot determine the optical power coupled into a
980 µm core without radiometry and coupling optics. Absolute results scale
linearly with incident power. Normalize to a power-meter reading before drawing
quantitative conclusions about the experiment.

## Provenance and parameter status

| Parameter | Value | Status |
|---|---:|---|
| Core diameter | 980 µm nominal (920–1040 µm) | Supplied GH-4001-P sheet, p. 1 |
| Cladding diameter | 1000 µm nominal | Supplied sheet, p. 1 |
| Jacket diameter | 2.20 mm nominal | Supplied sheet, p. 1 |
| Core index | 1.49 | Supplied sheet, p. 1; held constant versus wavelength |
| Numerical aperture | 0.50 | Supplied sheet, p. 1; held constant versus wavelength |
| Profile/material | Step index / PMMA / fluorinated polymer cladding | Supplied sheet |
| Loss at 650 nm | maximum 0.17 dB/m at 25 °C, 50% RH | Supplied sheet, p. 2; not a broadband constant |
| Minimum bend radius | 25 mm under specified quarter-bend loss test | Supplied sheet; not a general loss curve |
| Lamp | EKE, 150 W, 21 V | Dolan-Jenner manufacturer product page |
| Source spectrum | Manually estimated relative EKE curve anchors | User image, approximately 350–2100 nm |
| Attenuation spectrum | Manually estimated GH4001 curve anchors, 400–1000 nm | IFO typical ESKA graph, launch NA 0.1 |
| Input spot radius | 0.40 mm at 1/e² intensity | Adjustable assumption |
| Launch NA | 0.50 | Adjustable assumption; not inferred from the lamp |
| Optical input | 1 mW in selected spectral band | Adjustable assumption |

The graph source shows a typical loss curve, whereas the local datasheet gives
a maximum at one wavelength. They are not interchangeable. The sampled typical
curve uses 0.16 dB/m at 650 nm. Both CSV curves are coarse manual approximations,
not instrument measurements or high-precision digitization. Linear interpolation
is used for source radiance; log-linear interpolation for loss. The lamp image's
temperature labels are not used to invent a precise blackbody temperature. Its
dichroic-reflector curve is used directly instead.

Only 400–1000 nm or 400–700 nm contributes to calculated power. The source has
light beyond these bands; the program does not assume it is absent. The fiber
loss data provided here does not justify extrapolating into that infrared range.
Changing the band **holds entered optical input fixed in the newly selected
band**. This is not a physical filter-insertion calculation.

## Launch and transport equations

1. Normalize the chosen relative spectrum S(λ) over the selected band to the
   entered incident power. This is a user-defined power at the fiber-face plane,
   integrated over the entire assumed focused Gaussian, before spatial/angular
   acceptance and Fresnel loss. Relative radiance is assumed proportional to
   power spectrum under wavelength-independent coupling optics.
2. Launch intensity is a circular Gaussian exp(-2 r²/w²), displaced by d. The
   accepted fraction over the core disk is a noncentral chi-square CDF with two
   degrees of freedom, σ=w/2, threshold (a/σ)², noncentrality (d/σ)². Angular and
   spatial launch distributions are assumed independent.
3. The launch is uniform in projected solid angle (uniform radiance), or in
   transverse direction cosines sx,sy inside a disk of radius NA_launch. Its
   angular acceptance is min[1,(NA_fiber/NA_launch)²]. Rays beyond NA are rejected;
   leaky and skew-ray acceptance effects are omitted. Accepted output NA is
   min(NA_launch, NA_fiber); θmax = asin(NA_output) in air.
4. Apply unpolarized Fresnel power transmission at the air/PMMA entrance and
   reciprocal PMMA/air exit. Use the angle-dependent average of T², not the
   square of the average. The guided spectrum uses the single-face average T.
   The output spatial/angular distribution is weighted by the same T².
5. Bulk transport is 10^[-α(λ)L/10]. Extra user loss contributes a wavelength-
   independent multiplier 10^[-extra_dB/10]. Measured α is treated as loss per
   axial cable length. No extra ray-path multiplier is applied to that empirical
   loss. Application to launches unlike the source graph's NA=0.1 is approximate.
6. Spatial mixing is assumed to yield a uniform exit disk regardless of the
   accepted input spot. Angular distribution is conserved. These are modeling
   choices, including at zero/very short cable length, not measured mode mixing.

The nominal index relation gives n_clad = sqrt(1.49²−0.50²) ≈ 1.404, but there is
no vector eigenmode calculation. For context only, V=2πa NA/λ and M≈V²/2 give
approximately 2.8 million supported modes at 650 nm including two polarizations.
Supported capacity does not specify excited mode count or output coherence.

## Output irradiance

For axial distance z>0 and radial displacement r from a point on the output face,
uniform projected-solid-angle density maps to the plane kernel

    K(r,z) = z² / [π NA² (z²+r²)²],  r ≤ z tan(asin(NA)).

The r-dependent two-face Fresnel factor multiplies K before normalization. The
exit face is a uniform disk A(r), radius a=0.490 mm. The resulting irradiance is

    E(x,y,z) = P_out · (A * K)(x,y),

where A and K have unit integral. At z=0 use the uniform disk. A centered odd
257×257 grid contains the complete geometric envelope with 12% margin. Discrete
FFT convolution is padded by scipy.signal.fftconvolve, not circular convolution.
Small roundoff negatives are removed; the raster integral is normalized to
P_out. This prevents numerical loss from masquerading as absorption.

The finite grid smooths/quantizes boundaries; use a larger odd grid for exported
convergence studies (`footprint(settings,n=513)`). At large distances, the core
may occupy only a few cells but is already small compared with the footprint.
This is geometric radiometry: it does not claim diffraction accuracy at the
sharp aperture edge. The incoherent macroscopic model is appropriate for a large
core and broad thermal source under these assumptions.

Units: power mW; spectra mW/nm; geometric axes mm; irradiance mW/mm² internally.
**1 mW/mm² = 1000 W/m²**. Maps and UI report W/m². Export includes both units.
Geometric envelope diameter is 2[a+z tan(θmax)]; the reported 90%-power diameter
comes from an encircled-power integral and is smaller. The observation plane is
normal to the optical axis, with no sample tilt or lens.

## Local wave and interference explorer

The wave explorer is a separate, normalized educational model, not modulation
applied to the macroscopic irradiance. It uses 48 plane-wave components with
equal incoherent power. Their random phases and directions are seeded and
repeatable. Fields are summed as complex amplitudes before taking |E|²:

    I_coherent = |Σ a_j exp[i(kx_j x+ky_j y+φ_j)]|².
    I_incoherent = Σ |a_j|² = 1.
    Cross term = I_coherent − I_incoherent.

Negative cross terms indicate destructive interference relative to incoherent
addition; they are not negative physical irradiance. No independent per-frame
mean normalization is applied. Contrast is spatial standard deviation / mean
over the patch, not a measured speckle statistic of the real fiber.

Local angular spread is bounded by min(NA_output, a/sqrt(a²+z²)). This captures
the aperture subtended by the fiber at the on-axis sample in the far field,
instead of incorrectly resolving speckle everywhere at the full fiber NA.
Patch width = 12 λ/spread; 161 samples provide over six samples per shortest
intensity fringe. The circular aperture's exact angular correlations, full
eigenmode shapes, polarization structure and exact finite-aperture diffraction
are not solved. Geometrical scaling is approximate; patterns are illustrative.

Cable-length phase uses a sampled relative propagation constant
β_j−β_0 = k n [cos(θ_internal,j)−1]. It illustrates phase sensitivity, not a
fitted transmission matrix. The axial phase control adds (kz_j−k) Δz at fixed
reference geometry. It is not optical-cycle animation and does not change the
macroscopic observation plane. The finite ensemble option averages independent
random-phase intensity realizations (1, 4, 16 or 64). The halogen option instead
returns the exact ideal incoherent ensemble mean, with cross terms zero; its
basis wavelength has no effect on the intensity, only the illustrative axes.
These controls are not calibrated exposure time, source coherence, or bandwidth.

The temporal coherence chart independently integrates the actual modeled
output power spectrum:

    γ(δ) = ∫ S_out(λ) exp(i 2πδ/λ) dλ / ∫ S_out(λ) dλ.

δ is optical path difference. This is a visibility ceiling for equally powered,
equally polarized copies of the SAME spatial mode; it does not establish spatial
coherence across an extended halogen source. A narrowband filter alone cannot
make the whole lamp/fiber system coherent. The 0–15 µm chart uses 1 nm spectral
sampling; it is not extrapolated to macroscopic path differences, where coarse
spectral samples would alias.

## Rendering and performance

The default interface is now `fiber_simulation.py`, a native PySide6 window with
Matplotlib plots and an explicitly light palette. It uses the same `physics.py`
calculations as the browser version, with no change to numerical assumptions.
Numerical work runs in a single background worker, updates are debounced, and
stale results are dropped. Only the selected tab is drawn. The output plane is
displayed on a small polar mesh (18 radial by 49 angular samples); the wave
surface displays every fourth grid point. The full 257² irradiance and 161² wave
arrays are retained for calculation and CSV export. Matplotlib rotation does not
rerun the solver. Plot colors encode intensity, not perceived spectral color.

The following notes describe the retained browser interface (`app.py`):

Plotly/WebGL handles orbiting locally without rerunning Python. Representative
rays and a faint envelope explain geometry; they are not visible volumetric
scattering in clean air. The sample map is a cross-section. The fiber cutaway
shows only the final millimeters with true diameters; the full cable is included
in the transport calculation. The intensity surface's height encodes intensity,
not a physical sample surface. Colormaps encode intensity, not the actual color
the eye would see. The app serves bundled assets locally and does not need an
internet connection after installation. No GPU-specific Python stack is needed.

FFT footprints, spectra, wave bases and desktop calculation results use bounded
caches. Sliders update on release; geometry controls are debounced; camera
movement performs no physics calculations. The desktop OpenGL renderer uses the
full 161² local-field grid, indexed triangles and interpolated vertex colors.
Its detector display is a 65 × 193 polar mesh interpolated from the 257² footprint;
the 2D images use bilinear display interpolation. These are display operations,
not additional physical smoothing. CSV exports retain calculated grid values.
The vertical scale of the interference surface is relative intensity, with an
adjusted display aspect ratio; it is not a physical height. Native 3D rotation
does not redraw the Matplotlib charts. Unchanged tabs retain their views, while
changed charts reuse existing artists when their scale type stays the same.
The software fallback uses a reduced surface mesh. The original browser version
uses a 129² detector plane and 81² wave surface. No 3D voxel volume is recomputed.

## References

- User-provided Mitsubishi Chemical Corporation, *GH-4001-P Data Sheet*,
  14 November 2019, pp. 1–2. Local copy: `data/gh-4001-p-new.pdf`.
- User-provided *Quartz Halogen Lamps with Dichroic Reflectors* chart. Local
  image: `data/eke_source_reference.png`.
- Dolan-Jenner, [MI-152 illuminator specification](https://dolan-jenner.com/products/fiber-lite-mi-152-mi-157-illuminator-series-150w-halogen-light-source-mi152-mi157), accessed 3 October 2026.
- Industrial Fiber Optics, [Typical Transmission Loss Spectra of POF ESKA](https://i-fiberoptics.com/pdf/attenuation.pdf), accessed 3 October 2026. Local source PDF retained; curve is approximate and was obtained at launch NA 0.1.
- Rüdiger Paschotta, RP Photonics, [Case Study: Number of Modes of a Highly Multimode Fiber](https://www.rp-photonics.com/case_study_number_of_modes.html), for the large-V mode-capacity estimate.

## Measurements that would make this predictive

Measure output power in a known spectral band, the emitted spectrum, exit-face
image, and far-field angular distribution. Record actual lamp/filter/iris
settings, launch spot/NA/offset, end finish and bend routing. The easiest initial
calibration is multiplying incident power by measured_output / predicted_output
under matching spectral and geometry conditions. Replace CSV curves with measured
spectral data when available. Then add the sample's measured BRDF/absorption and
surface geometry to propagate the field into the regolith/composite experiment.
