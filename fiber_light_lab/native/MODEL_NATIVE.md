# Fiber Light Lab (native) — model notes

## Scope

This program follows light from the halogen illuminator to the goniometer detector:

    MI-152 (EKE lamp) -> GH4001 plastic fiber (10 in) -> Thorlabs LB1092 lens ->
    regolith-epoxy coupon at the rotation center -> detector on the moving arm

It answers instrument questions with real ray tracing: how big the spot is, how much
light reaches the coupon, how wide a mirror scan is, how much signal the detector gets,
and which angles the arm can reach. It also predicts angular scans of the coupons from a
physically based *hypothesis* model of the composite surface. The coupon model is not a
measurement. The experiment decides whether finishing or regolith content dominates.

The fiber transport is the original Python model (`physics.py`), reproduced exactly
(see Verification). Everything after the fiber exit face is new.

## Coordinates and geometry

Units are millimeters and degrees. The origin is the goniometer rotation center on the
nominal sample surface. +z is the sample normal, +x lies in the plane of incidence on the
specular side, and +y is the rotation axis. The source arm (fiber + lens) is fixed at
angle −θi. The detector moves on a circle of radius R in the x-z plane, so detector angle
+θi is the mirror direction. This is the geometry suggested in the research plan
(fiber + one lens fixed, detector moving around the sample).

Placeholders (edit them in the app): lens center 63.5 mm (2.5 in) from the rotation
center, detector radius 100 mm, coupon 25.4 × 25.4 × 5 mm, lens tube outer radius 8.9 mm.

## Parameters and provenance

| Item | Value | Source / status |
|---|---|---|
| Fiber | GH-4001-P, PMMA step index, core 980 µm (920–1040), n = 1.49, NA 0.50 | Mitsubishi data sheet (data/gh-4001-p-new.pdf) |
| Fiber loss | Typical ESKA curve, 400–1000 nm, launch NA 0.1; 0.17 dB/m max @ 650 nm | data/gh4001_loss.csv (manual digitization) |
| Lamp | Dolan-Jenner MI-152, EKE 150 W 21 V, dichroic reflector | data/eke_relative_spectrum.csv (manual digitization) |
| Lamp rating | 350,000 fc at the light-guide port | Dolan-Jenner product page; used only in "MI-152 rating" mode |
| Input power | 1 mW in the band (default) | Assumption, as in the Python app |
| Lens | Thorlabs LB1092: N-BK7 bi-convex, Ø12.7 mm, R1 = +14.6, R2 = −14.6 mm, CT 4.7 mm, ET 1.8 mm, f = 15.0 mm, BFL 13.4 mm (at 587.6 nm), uncoated | Thorlabs catalog listing |
| N-BK7 | Sellmeier B = 1.03961212, 0.231792344, 1.01046945; C = 0.00600069867, 0.0200179144, 103.560653 µm² | Schott data sheet |
| Detector presets | Thorlabs SM05PD1A / PDA36A2: Si 3.6 × 3.6 mm; S120C: Si Ø9.5 mm; BPW34: 7.5 mm² | Manufacturer specs; actual detector not yet verified |
| Si responsivity | Typical curve, peak ≈ 0.63 A/W near 950 nm | data/si_photodiode_responsivity_typical.csv (approximate; load your calibration) |
| CIE V(λ) | 1924 photopic | data/cie_photopic_vlambda.csv |
| LHS-1 grains | Hapke w 0.388, b 0.063, h 0.221, θ̄ 14.27°, φ 0.40; density 2.75 g/cm³ | Curtis et al. 2024 Apollo 16 (68810) smooth fit as highland analog; density estimated |
| LMS-1 grains | Hapke w 0.281, b 0.238, h 0.032, θ̄ 13.80°, φ 0.60; density 2.92 g/cm³ | Curtis et al. 2024 Apollo 11 (10084) smooth fit as mare analog; density estimated |
| Epoxy | n = 1.55, density 1.15 g/cm³ | Typical cured bisphenol-A epoxy; replace with your resin's data |
| Surface | RMS slope 1° as-cast, 12° abraded | Assumptions; measure Rdq if possible |

The simulants are immature (no space weathering), so they are probably brighter than the
Apollo soils used as analogs. Treat the grain parameters as starting values.

## 1. Fiber transport (unchanged from the Python app)

The EKE spectrum S(λ) is normalized over the band to the input power. Launch acceptance is
min[1, (0.50/NA_launch)²] times the fraction of a Gaussian spot on the core (noncentral
chi-square CDF). The two-face Fresnel factor T² is averaged over a uniformly filled NA disk.
Bulk loss is 10^(−αL/10) with α(λ) from the typical curve. The output is a uniform 980 µm
disk radiating uniformly in direction cosines inside NA 0.50 (uniform radiance).

In "MI-152 rating" mode, 350,000 fc (3.77 lm/mm²) over the 0.754 mm² core face gives
2.84 lm. Dividing by the luminous efficacy of the EKE spectrum in the band gives the
optical power on the core face (about 15 mW for 400–1000 nm). This is a manufacturer
claim at full intensity, not a calibration.

## 2. Ray source

Rays are drawn with a 5-D low-discrepancy sequence: wavelength from the output spectrum,
position uniform on the core, direction uniform in the NA disk. Each ray is weighted by
the angle-dependent two-face Fresnel factor, and the total equals the fiber output power
exactly. Because every ray has its own wavelength, lens dispersion is included.

## 3. LB1092 lens (exact ray trace)

Each ray is intersected with the two spherical surfaces and refracted with vector Snell's
law, using n_BK7(λ). Each surface transmits 1 − R, where R is the unpolarized Fresnel
reflectance at the actual angle (uncoated), or 0.5 % scaled by the Fresnel angle
dependence inside an AR band (-A, -B, -AB options). Rays outside the Ø12.7 mm aperture
(or a smaller mount aperture) are blocked. Rays leaving through the ground edge are
absorbed, and total internal reflection is tracked.

The paraxial focal lengths come from the thick-lens equations. With the catalog radii,
EFL = 14.95 mm and BFL = 13.31 mm at 587.6 nm, within the ±1 % focal-length tolerance of
the catalog 15.0 / 13.4 mm. "Collimate" puts the fiber at the paraxial front focus.
"Image fiber on sample" uses the paraxial conjugate of the rotation center.

## 4. Sample and detector

Rays are intersected with the (optionally tilted) sample plane and kept if they land on
the coupon. For every detector angle, the reflected power is integrated over the active
area:

    dP = P_ray · f_r(ω_i, ω_o) · cos θ_o · cos θ_det · dA / d²

The integral uses a stratified grid on the active area whose density follows the
narrowest BRDF lobe. The ideal mirror is handled by reflecting each ray and testing
whether it lands on the active area. Spot size, beam divergence, lens aberrations,
detector size and arm radius therefore all shape the scan, as on the real bench.

Outputs: optical power on the detector, and photocurrent from the responsivity curve
(per ray wavelength). A power meter set to λ_set reports photocurrent / R(λ_set), which is
what a broadband halogen measurement would show. The amplifier voltage is
current × transimpedance. The inferred BRDF is P / (P_on_coupon · Ω · cos θ_det).

The detector cannot reach angles within the "blocked" zone around −θi. Its half-width is
atan(r_source/L_source) + atan(r_detector/R).

Metrics follow the research plan. "Near" is the mean signal within ±5° of the specular
angle. "Away" is the mean signal more than 20° from it. Both values and their ratio are
reported, with blocked angles excluded. The mirror scan (R = 1) shown with every result is
the instrument function.

## 5. Coupon reflectance model

    f = (1 − f_exp) · [f_resin(σ) + f_body] + f_exp · f_Hapke  (+ f_back for neat epoxy)

Resin surface: Beckmann microfacet lobe, α = √2 · tan(RMS slope), with the Fresnel
reflectance of the resin (n = 1.55) and Smith shadowing (Walter et al. 2007).

Body (light that enters the resin and is returned by embedded grains): a
Saunderson-corrected Lambertian term,

    f_body = (1 − F(θi)) (1 − F(θo)) R_b / [π n² (1 − ρ_i R_b)],

where ρ_i ≈ 0.62 is the internal diffuse reflectance of the resin surface. R_b equals the
dry-powder reflectance (1 − γ)/(1 + γ), with γ = √(1 − w), times an "embedded darkening"
factor (default 0.5). Grains surrounded by resin scatter less than grains in air (Lekner &
Dorf 1988). Plagioclase, the main LHS-1 mineral, has nearly the same index as epoxy, so
this factor is very uncertain. For neat epoxy, R_b is the backing reflectance.

Exposed grains: abrasion cuts through grains. By Delesse's principle, the area fraction of
grains on a cut face equals their volume fraction:

    φ_v = (m/ρ_g) / [m/ρ_g + (1 − m)/ρ_resin],

so 30/45/60 wt% gives φ_v = 0.15/0.25/0.38 (LHS-1). The abraded exposed fraction is
f_exp = factor · φ_v. The as-cast face is a resin skin (f_exp = 0 by default). Exposed
grains follow the Hapke (2012) BRDF with porosity correction, shadow-hiding opposition
effect and macroscopic roughness (Hapke 1984). The double Henyey–Greenstein phase uses c
from Hapke's hockey-stick relation, because the published c values were not available.

Neat-epoxy back face: a transparent coupon also reflects from its back surface. This is
modeled as a second microfacet lobe with reflectance (1 − F)² · F_internal and a width
widened by the front-surface roughness on the way in and out.

Consequences built into the assumptions (things the experiment can falsify):
- As-cast composites are independent of regolith content (resin skin over an opaque
  body). Measured content dependence before abrasion would mean grains near the surface
  or a translucent body.
- Abrasion lowers the near/away ratio through two effects: a rougher resin surface
  (wider lobe) and exposed grains (diffuse, backscattering).

## What the model predicts for the current bench

These come from the default settings (collimated LB1092, lens 63.5 mm from the sample,
3.6 mm detector at 100 mm). They are predictions to check with the mirror test, not facts.

- The fiber's NA 0.50 cone overfills the Ø12.7 mm lens: about 38 % of the fiber output hits
  the mount, 5 % is reflected by the two bare surfaces, and about 55 % reaches the coupon.
  Underfilling the fiber (launch NA ≈ 0.25) would put about 92 % on the coupon.
- The beam is not well collimated. Spherical aberration of the fast bi-convex lens makes
  marginal rays cross the axis, so the beam narrows to about 6 mm near 60 mm and 90 % of
  the power is within ±6.5°. The source-size limit from the research plan (core/f ≈ 3.7°
  full) holds only for the paraxial core of the beam.
- A mirror scan is about 5.8° wide (FWHM), with wings out to about ±13° on a log scale.
  At R = 100 mm the width comes mostly from the beam footprint plus detector size divided
  by R. At R = 300–500 mm it falls to about 4.5–4.8°, set by the source-size divergence.
  Stopping the lens down to 6 mm removes the aberration wings (90 % spread 3.7°) but cuts
  the power about 4×.
- With the uncoated lens and 1 mW into the fiber, the mirror peak is about 40 µW
  (about 17 µA). An as-cast composite peaks at about 1.5 µW near specular and about 2 nW
  away from it (about 1 nA). An abraded one gives tens of nW in both windows. That is
  about four decades of photocurrent, so plan for switchable gain (for example the
  PDA36A2 gain settings) or a picoammeter, and confirm linearity across it.
- At 60° illumination the spot stretches to about 10 × 5 mm. That still fits a 25 mm
  coupon, but check against your real coupon size.

## Verification

All checks run with `cargo test` in `native/`.
- **Parity with the Python model:** every transport output agrees to 1e-12 relative for
  seven cases, including the published default output of 0.7843244 mW
  (tests/parity_python.rs). The fixture is written by
  `verification/make_native_reference.py`.
- **Spatial coupling:** agrees with scipy's ncx2 to 1e-12.
- **Bare-fiber footprint:** the exact radial integral matches Python's FFT grid (peak
  within 0.3 %), and the Monte Carlo trace matches the exact integral to 0.2 % in
  encircled power.
- **N-BK7:** within 2e-5 of Schott's nd, nF and nC.
- **LB1092:** EFL within the ±1 % catalog tolerance and BFL within 0.15 mm of the
  catalog value. The exact trace reproduces the paraxial BFL, and normal-incidence
  transmission is (1 − R0)².
- **Power budget:** closes to 1e-9.
- **PDF check:** with a small launch NA, the divergence equals core/f (3.7°).
- **Lambertian sample:** the scan follows P·R·cos θ·Ω/π within 5 % (finite spot).
- **Mirror:** reflected power is conserved, and a very smooth resin lobe converges to a
  mirror with R = F(θ).
- **BRDFs:** reciprocity holds, every finish/simulant/loading has albedo below 1, a
  smooth-resin lobe has albedo F(θ), and the Saunderson closed form is reproduced.

These are numerical self-consistency checks. They are not validation against your
apparatus.

## Limitations

- Fiber: straight, uniform output radiance, no modal or bend effects beyond the user loss.
  The fiber loss and lamp curves are approximate digitizations.
- Lens: no ghost reflections, no bulk absorption (under 0.3 % for 4.7 mm of N-BK7), and
  surface reflections averaged over polarization at each surface.
- Unmodeled light: rays that miss the lens or coupon are counted as lost; scattering from
  the mount, the stage and the coupon edges is not modeled.
- Detector: assumed to face the rotation center, with no window or angular response.
  Responsivity is a typical curve.
- Coupon: gray (no wavelength dependence of grain albedo) and in-plane scans only. The
  model parameters are literature analogs and assumptions.
- Composite body: the exposed-grain area and body reflectance are the least certain
  parameters. The model cannot predict coupon translucency, grain settling, dust, or
  polishing residue.
- Scope: the scan measures angular reflected brightness for this apparatus. It does not
  give solar albedo, emissivity or durability.

## Measurements that would make this predictive

1. Power out of the fiber (power meter at the tip, with the band noted). Use it as the input
   power or as a scale factor.
2. A mirror scan at 30° and 60° compared with the instrument function shown here. This
   checks alignment, spot size and resolution. Add a fiber offset or sample tilt in the
   app to reproduce any asymmetry.
3. A diffuse-reference scan to check the cos θ law and detector linearity.
4. The real detector model, its active area and its calibration curve (load the CSV).
5. Coupon size, resin index and density, and the surface slope (Rdq) before and after
   abrasion.
6. Load measured scans (angle, signal CSV) on the scan tab to compare shape with the
   model, then adjust exposure factor, slopes and darkening.

## References

- Mitsubishi Chemical, GH-4001-P data sheet (2019). Local copy: data/gh-4001-p-new.pdf.
- Industrial Fiber Optics, [Typical Transmission Loss Spectra of POF ESKA](https://i-fiberoptics.com/pdf/attenuation.pdf).
- Dolan-Jenner, [Fiber-Lite MI-150 series](https://dolan-jenner.com/products/fiber-lite-mi-152-mi-157-illuminator-series-150w-halogen-light-source-mi152-mi157).
- Thorlabs, [N-BK7 Bi-Convex Lenses, Uncoated](https://www.thorlabs.com/n-bk7-bi-convex-lenses-uncoated) (LB1092).
- Schott, N-BK7 optical glass data sheet (Sellmeier coefficients).
- Curtis, Warren, Shirley, Paige & Bowles (2024), [BRDF measurements of characterized Apollo regolith samples](https://onlinelibrary.wiley.com/doi/full/10.1111/maps.14266), Meteoritics & Planetary Science.
- Calisi, Caporali & Taurino (2026), [Resin-regolith composites for non-structural applications](https://flore.unifi.it/handle/2158/1470953), Materials 19, 2066.
- Hapke, B. (1984) Icarus 59, 41 (roughness); Hapke (2008) Icarus 195, 918 (porosity); Hapke (2012) *Theory of Reflectance and Emittance Spectroscopy*, 2nd ed.; Hapke (2012) Icarus 221, 1079 (hockey stick).
- Walter, Marschner, Li & Torrance (2007), Microfacet models for refraction through rough surfaces, EGSR.
- Saunderson (1942), J. Opt. Soc. Am. 32, 727; Lekner & Dorf (1988), Appl. Opt. 27, 1278.
- Newport, [Beam Focusing and Collimating](https://www.newport.com/n/focusing-and-collimating).
