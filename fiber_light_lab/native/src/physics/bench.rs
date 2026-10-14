//! Goniometer bench geometry and the Monte Carlo beam trace.
//!
//! World frame (mm): origin at the goniometer rotation center on the nominal
//! sample surface; +z is the nominal sample normal; +x lies in the plane of
//! incidence toward the specular side; +y is the rotation axis. The source arm
//! (fiber + optional LB1092) sits at angle -θi and stays fixed during a scan;
//! the detector moves on a circle of radius R in the x-z plane.
//!
//! Rays leave the GH4001 exit face uniformly over the 0.98 mm core and
//! uniformly in direction cosines inside the output NA (uniform radiance),
//! weighted by the two-face Fresnel factor — the same source model as the
//! Python app. Each ray carries a wavelength drawn from the modeled output
//! spectrum, so N-BK7 dispersion (chromatic aberration) is included.

use super::fiber::{CORE_RADIUS_MM, FiberSettings, Transport, face_transmission_averages, fresnel_transmission};
use super::lens::{LensLoss, LensSpec};
use super::math::{Frame, PI, RSequence, TAU, Vec3};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize, Default)]
pub enum CouponShape {
    #[default]
    Rectangle,
    Disk,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct Coupon {
    pub shape: CouponShape,
    /// Rectangle: extent along the plane of incidence (x). Disk: diameter.
    pub length_mm: f64,
    /// Rectangle: extent along the rotation axis (y).
    pub width_mm: f64,
    pub thickness_mm: f64,
}

impl Default for Coupon {
    fn default() -> Self {
        // Placeholder: 1 in × 1 in × 5 mm. Edit to the cast coupon size.
        Coupon { shape: CouponShape::Rectangle, length_mm: 25.4, width_mm: 25.4, thickness_mm: 5.0 }
    }
}

impl Coupon {
    pub fn contains(&self, u: f64, v: f64) -> bool {
        match self.shape {
            CouponShape::Rectangle => u.abs() <= 0.5 * self.length_mm && v.abs() <= 0.5 * self.width_mm,
            CouponShape::Disk => u.hypot(v) <= 0.5 * self.length_mm,
        }
    }
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct BenchSettings {
    /// Place the LB1092 between fiber and sample (false = bare fiber, as the
    /// original Python app).
    pub use_lens: bool,
    pub lens: LensSpec,
    /// Fiber exit face to lens front vertex (mm).
    pub fiber_to_lens_mm: f64,
    /// Lateral misalignment of the fiber tip relative to the lens axis (mm).
    pub fiber_offset_x_mm: f64,
    pub fiber_offset_y_mm: f64,
    /// With lens: lens center to rotation center. Without: fiber tip to center.
    pub source_distance_mm: f64,
    /// Illumination angle from the sample normal (deg).
    pub incidence_deg: f64,
    /// Sample misalignment: tilt about the rotation axis (in-plane, deg).
    pub sample_tilt_deg: f64,
    /// Sample misalignment: tilt about the x axis (out of plane, deg).
    pub sample_roll_deg: f64,
    pub coupon: Coupon,
    /// Outer radius of the source assembly (lens tube), for arm collisions.
    pub source_housing_radius_mm: f64,
}

impl Default for BenchSettings {
    fn default() -> Self {
        let lens = LensSpec::lb1092();
        BenchSettings {
            use_lens: true,
            fiber_to_lens_mm: lens.collimating_distance(587.6),
            lens,
            fiber_offset_x_mm: 0.0,
            fiber_offset_y_mm: 0.0,
            source_distance_mm: 63.5,
            incidence_deg: 30.0,
            sample_tilt_deg: 0.0,
            sample_roll_deg: 0.0,
            coupon: Coupon::default(),
            // SM05 lens tube: 0.535-40 thread, ~0.7 in outer diameter.
            source_housing_radius_mm: 8.9,
        }
    }
}

impl BenchSettings {
    pub fn validate(&self) -> Result<(), String> {
        if self.use_lens {
            self.lens.validate()?;
            if self.fiber_to_lens_mm <= 0.0 {
                return Err("Fiber-to-lens distance must be positive.".into());
            }
            if self.source_distance_mm <= 0.5 * self.lens.center_thickness_mm + 0.5 {
                return Err("Lens must sit in front of the sample.".into());
            }
        }
        if !(0.0..89.0).contains(&self.incidence_deg) {
            return Err("Incidence angle must be within 0-89°.".into());
        }
        if self.source_distance_mm < 0.0 || self.coupon.length_mm <= 0.0 || self.coupon.width_mm <= 0.0 {
            return Err("Distances and coupon size must be positive.".into());
        }
        Ok(())
    }

    /// Unit vector from the rotation center toward the source.
    pub fn source_direction(&self) -> Vec3 {
        let t = self.incidence_deg.to_radians();
        Vec3::new(-t.sin(), 0.0, t.cos())
    }

    /// Optical frame: w = optical axis toward the sample, v = rotation axis.
    pub fn optical_frame(&self) -> Frame {
        Frame::from_w_and_v(-self.source_direction(), Vec3::Y)
    }

    /// Distance along the axis from the fiber tip to the rotation center.
    pub fn tip_to_center_mm(&self) -> f64 {
        if self.use_lens {
            self.fiber_to_lens_mm + 0.5 * self.lens.center_thickness_mm + self.source_distance_mm
        } else {
            self.source_distance_mm
        }
    }

    /// World position of the nominal (unshifted) fiber tip.
    pub fn fiber_tip(&self) -> Vec3 {
        self.source_direction() * self.tip_to_center_mm()
    }

    /// World position of the lens front vertex on the axis.
    pub fn lens_front_vertex(&self) -> Vec3 {
        self.fiber_tip() - self.source_direction() * self.fiber_to_lens_mm
    }

    /// Sample surface normal including misalignment.
    pub fn sample_normal(&self) -> Vec3 {
        Vec3::Z.rotate_x(self.sample_roll_deg.to_radians()).rotate_y(self.sample_tilt_deg.to_radians())
    }

    /// Frame on the sample surface: u along the plane of incidence, v along y.
    pub fn sample_frame(&self) -> Frame {
        let n = self.sample_normal();
        let u = Vec3::Y.cross(n).normalized();
        let v = n.cross(u);
        Frame { u, v, w: n }
    }

    /// Direction of the mirror reflection of the central ray.
    pub fn specular_direction(&self) -> Vec3 {
        super::lens::reflect(-self.source_direction(), self.sample_normal())
    }

    /// Detector angle (deg, from the nominal normal) of the specular ray.
    pub fn specular_angle_deg(&self) -> f64 {
        let r = self.specular_direction();
        r.x.atan2(r.z).to_degrees()
    }
}

/// One ray that reached the sample plane.
#[derive(Clone, Copy, Debug)]
pub struct Hit {
    pub position: Vec3,
    pub direction: Vec3,
    /// Optical power carried (mW).
    pub power_mw: f64,
    pub wavelength_nm: f32,
    /// Sample-frame coordinates (mm).
    pub u: f32,
    pub v: f32,
    pub on_coupon: bool,
}

#[derive(Clone, Copy, Debug, Default, Serialize)]
pub struct PowerBudget {
    /// Fiber output power in the band (mW).
    pub fiber_output_mw: f64,
    pub blocked_by_mount_mw: f64,
    pub lens_reflection_mw: f64,
    pub lens_edge_mw: f64,
    pub lens_tir_mw: f64,
    pub past_lens_mw: f64,
    pub off_coupon_mw: f64,
    pub on_coupon_mw: f64,
}

#[derive(Clone, Debug)]
pub struct RayPath {
    pub points: Vec<Vec3>,
    pub wavelength_nm: f32,
    pub power_mw: f32,
    pub reached_sample: bool,
}

#[derive(Clone, Debug, Default)]
pub struct BeamResult {
    pub hits: Vec<Hit>,
    pub paths: Vec<RayPath>,
    pub budget: PowerBudget,
    pub rays_traced: usize,
}

/// Sample wavelengths from the output spectrum by inverse CDF (1 nm bins).
pub struct SpectrumSampler {
    wl: Vec<f64>,
    cdf: Vec<f64>,
}

impl SpectrumSampler {
    pub fn new(t: &Transport) -> Self {
        let mut cdf = vec![0.0; t.wavelength_nm.len()];
        for i in 1..cdf.len() {
            cdf[i] = cdf[i - 1]
                + 0.5 * (t.output_spectrum[i] + t.output_spectrum[i - 1]) * (t.wavelength_nm[i] - t.wavelength_nm[i - 1]);
        }
        let total = *cdf.last().unwrap();
        if total > 0.0 {
            cdf.iter_mut().for_each(|c| *c /= total);
        } else {
            cdf = (0..cdf.len()).map(|i| i as f64 / (cdf.len() - 1) as f64).collect();
        }
        SpectrumSampler { wl: t.wavelength_nm.clone(), cdf }
    }

    pub fn sample(&self, u: f64) -> f64 {
        let i = self.cdf.partition_point(|&c| c < u).clamp(1, self.cdf.len() - 1);
        let span = (self.cdf[i] - self.cdf[i - 1]).max(1e-300);
        let f = ((u - self.cdf[i - 1]) / span).clamp(0.0, 1.0);
        self.wl[i - 1] + f * (self.wl[i] - self.wl[i - 1])
    }
}

const PATHS_KEPT: usize = 900;

enum Fate {
    Mount,
    Edge,
    Tir,
    /// Reached the sample plane; power after the optics.
    Plane(Hit),
    /// Left the bench without meeting the sample plane.
    Escaped(f64),
}

/// Trace `n_rays` rays from the fiber exit face to the sample plane.
pub fn trace_beam(bench: &BenchSettings, fiber: &FiberSettings, t: &Transport, n_rays: usize, seed: u64) -> BeamResult {
    let n_rays = n_rays.max(1);
    let na = fiber.output_na();
    let (_, t2) = face_transmission_averages(na);
    let sampler = SpectrumSampler::new(t);
    let seq = RSequence::new(5, seed);
    let frame = bench.optical_frame();
    let tip = bench.fiber_tip();
    let sframe = bench.sample_frame();
    let normal = sframe.w;
    let base_power = t.output_mw / n_rays as f64;
    // Misalignment shifts the fiber face in the lens-local transverse plane.
    let tip_offset = Vec3::new(bench.fiber_offset_x_mm, bench.fiber_offset_y_mm, 0.0);
    let to_world = |p: Vec3| tip + frame.to_world(p);

    let results: Vec<(Fate, f64, Option<RayPath>)> = (0..n_rays)
        .into_par_iter()
        .map(|i| {
            let mut u = [0.0; 5];
            seq.point(i as u64, &mut u);
            let wl = sampler.sample(u[0]);
            let r = CORE_RADIUS_MM * u[1].sqrt();
            let phi = TAU * u[2];
            let s = na * u[3].sqrt();
            let psi = TAU * u[4];
            let ft = fresnel_transmission(s);
            let launched = base_power * ft * ft / t2;
            // Lens-local coordinates: fiber face at z = 0, optical axis +z.
            let o_local = Vec3::new(r * phi.cos(), r * phi.sin(), 0.0) + tip_offset;
            let d_local = Vec3::new(s * psi.cos(), s * psi.sin(), (1.0 - s * s).sqrt());
            let keep = i < PATHS_KEPT;
            let mut pts = Vec::new();
            if keep {
                pts.push(to_world(o_local));
            }
            let path = |pts: Vec<Vec3>, power: f64, reached: bool| {
                keep.then(|| RayPath { points: pts, wavelength_nm: wl as f32, power_mw: power as f32, reached_sample: reached })
            };
            let (origin, dir, power) = if bench.use_lens {
                match bench.lens.trace(o_local, d_local, wl, bench.fiber_to_lens_mm) {
                    Ok(tr) => {
                        if keep {
                            pts.push(to_world(tr.front_hit));
                            pts.push(to_world(tr.rear_hit));
                        }
                        (to_world(tr.rear_hit), frame.to_world(tr.direction), launched * tr.transmission)
                    }
                    Err(loss) => {
                        if keep {
                            // Draw the lost ray up to the lens front plane.
                            let tt = (bench.fiber_to_lens_mm - o_local.z) / d_local.z;
                            pts.push(to_world(o_local + d_local * tt));
                        }
                        let fate = match loss {
                            LensLoss::Aperture => Fate::Mount,
                            LensLoss::Edge => Fate::Edge,
                            LensLoss::Tir => Fate::Tir,
                        };
                        return (fate, launched, path(pts, launched, false));
                    }
                }
            } else {
                (to_world(o_local), frame.to_world(d_local), launched)
            };
            let denom = dir.dot(normal);
            let tt = if denom < -1e-12 { -origin.dot(normal) / denom } else { -1.0 };
            if tt <= 0.0 {
                if keep {
                    pts.push(origin + dir * 50.0);
                }
                return (Fate::Escaped(power), launched, path(pts, power, false));
            }
            let p = origin + dir * tt;
            let (su, sv) = (p.dot(sframe.u), p.dot(sframe.v));
            let on = bench.coupon.contains(su, sv);
            if keep {
                pts.push(p);
            }
            let hit = Hit {
                position: p,
                direction: dir,
                power_mw: power,
                wavelength_nm: wl as f32,
                u: su as f32,
                v: sv as f32,
                on_coupon: on,
            };
            (Fate::Plane(hit), launched, path(pts, power, on))
        })
        .collect();

    let mut out = BeamResult { rays_traced: n_rays, ..Default::default() };
    out.budget.fiber_output_mw = t.output_mw;
    // Renormalize so the launched rays carry exactly the fiber output power
    // (the sample mean of the Fresnel weight differs from t2 at ~1e-8).
    let launched_total: f64 = results.iter().map(|r| r.1).sum();
    let scale = if launched_total > 0.0 { t.output_mw / launched_total } else { 0.0 };
    let mut accounted = 0.0;
    for (mut fate, launched, path) in results {
        let launched = launched * scale;
        match &mut fate {
            Fate::Plane(h) => h.power_mw *= scale,
            Fate::Escaped(p) => *p *= scale,
            _ => {}
        }
        if let Some(p) = path {
            out.paths.push(p);
        }
        match fate {
            Fate::Mount => out.budget.blocked_by_mount_mw += launched,
            Fate::Edge => out.budget.lens_edge_mw += launched,
            Fate::Tir => out.budget.lens_tir_mw += launched,
            Fate::Plane(h) => {
                out.budget.past_lens_mw += h.power_mw;
                if h.on_coupon {
                    out.budget.on_coupon_mw += h.power_mw;
                } else {
                    out.budget.off_coupon_mw += h.power_mw;
                }
                accounted += h.power_mw - launched;
                out.hits.push(h);
            }
            Fate::Escaped(p) => {
                out.budget.past_lens_mw += p;
                out.budget.off_coupon_mw += p;
                accounted += p - launched;
            }
        }
    }
    // Launched minus transmitted power of the rays that went through the
    // lens is what its two surfaces reflected.
    out.budget.lens_reflection_mw = (-accounted).max(0.0);
    out
}

/// Beam diagnostics on planes normal to the optical axis.
#[derive(Clone, Debug, Default, Serialize)]
pub struct BeamProfile {
    /// Distance from the lens center (or bare fiber tip) along the axis (mm).
    pub distance_mm: Vec<f64>,
    /// Second-moment (D4σ-equivalent) diameter about the centroid (mm).
    pub d_second_moment_mm: Vec<f64>,
    /// Diameter enclosing 90 % of the power about the centroid (mm).
    pub d90_mm: Vec<f64>,
    /// Full cone angle enclosing 90 % of the power after the optics (deg).
    pub divergence90_full_deg: f64,
    /// Same, for the bare-geometry estimate core/f (deg).
    pub divergence_estimate_deg: f64,
    pub angle_hist_deg: Vec<f64>,
    pub angle_hist_power: Vec<f64>,
}

pub fn beam_profile(bench: &BenchSettings, beam: &BeamResult, wavelength_nm: f64, max_distance_mm: f64) -> BeamProfile {
    let frame = bench.optical_frame();
    let axis = frame.w;
    // Reference: lens center (or fiber tip) on the axis.
    let center = if bench.use_lens {
        bench.fiber_tip() + axis * (bench.fiber_to_lens_mm + 0.5 * bench.lens.center_thickness_mm)
    } else {
        bench.fiber_tip()
    };
    // Reconstruct rays just after the optics: back-propagate sample hits to
    // the reference plane (straight lines in air).
    let rays: Vec<(Vec3, Vec3, f64)> = beam
        .hits
        .iter()
        .step_by((beam.hits.len() / 60_000).max(1))
        .map(|h| (h.position, h.direction, h.power_mw))
        .collect();
    let mut prof = BeamProfile::default();
    let n = 61;
    for k in 0..n {
        let z = max_distance_mm * k as f64 / (n - 1) as f64;
        let plane = center + axis * z;
        let mut pts: Vec<(f64, f64, f64)> = Vec::with_capacity(rays.len());
        for &(p, d, w) in &rays {
            let dz = d.dot(axis);
            if dz <= 0.0 {
                continue;
            }
            let t = (plane - p).dot(axis) / dz;
            let q = p + d * t;
            let rel = q - plane;
            pts.push((rel.dot(frame.u), rel.dot(frame.v), w));
        }
        let wsum: f64 = pts.iter().map(|p| p.2).sum();
        if wsum <= 0.0 {
            continue;
        }
        let cx = pts.iter().map(|p| p.0 * p.2).sum::<f64>() / wsum;
        let cy = pts.iter().map(|p| p.1 * p.2).sum::<f64>() / wsum;
        let mut r: Vec<(f64, f64)> = pts.iter().map(|p| ((p.0 - cx).hypot(p.1 - cy), p.2)).collect();
        let m2 = r.iter().map(|(rr, w)| rr * rr * w).sum::<f64>() / wsum;
        r.sort_by(|a, b| a.0.total_cmp(&b.0));
        let mut acc = 0.0;
        let mut r90 = 0.0;
        for (rr, w) in &r {
            acc += w;
            if acc >= 0.9 * wsum {
                r90 = *rr;
                break;
            }
        }
        prof.distance_mm.push(z);
        prof.d_second_moment_mm.push(2.0 * (2.0 * m2).sqrt());
        prof.d90_mm.push(2.0 * r90);
    }
    // Angular distribution relative to the axis.
    let mut ang: Vec<(f64, f64)> = rays.iter().map(|&(_, d, w)| (d.dot(axis).clamp(-1.0, 1.0).acos().to_degrees(), w)).collect();
    let total: f64 = ang.iter().map(|a| a.1).sum();
    ang.sort_by(|a, b| a.0.total_cmp(&b.0));
    let mut acc = 0.0;
    for (a, w) in &ang {
        acc += w;
        if acc >= 0.9 * total {
            prof.divergence90_full_deg = 2.0 * a;
            break;
        }
    }
    let max_angle = if bench.use_lens { 15.0 } else { 35.0 };
    let bins = 150;
    prof.angle_hist_deg = (0..bins).map(|i| (i as f64 + 0.5) * max_angle / bins as f64).collect();
    prof.angle_hist_power = vec![0.0; bins];
    for (a, w) in &ang {
        let i = (a / max_angle * bins as f64) as usize;
        if i < bins {
            prof.angle_hist_power[i] += w;
        }
    }
    prof.divergence_estimate_deg = if bench.use_lens {
        (2.0 * CORE_RADIUS_MM / bench.lens.efl(wavelength_nm)).to_degrees()
    } else {
        2.0 * (PI / 6.0).to_degrees()
    };
    prof
}

/// Irradiance on the sample plane (sample frame u, v), W/m².
#[derive(Clone, Debug, Default, Serialize)]
pub struct SpotMap {
    pub n: usize,
    pub half_extent_mm: f64,
    /// Row-major [v][u] irradiance (W/m²).
    pub irradiance_w_m2: Vec<f64>,
    pub centroid_u_mm: f64,
    pub centroid_v_mm: f64,
    /// Widths enclosing 90 % of the power along u and v (mm).
    pub width90_u_mm: f64,
    pub width90_v_mm: f64,
    pub peak_w_m2: f64,
}

pub fn spot_map(beam: &BeamResult, n: usize, half_extent_mm: Option<f64>) -> SpotMap {
    let hits = &beam.hits;
    let total: f64 = hits.iter().map(|h| h.power_mw).sum();
    let mut map = SpotMap { n, ..Default::default() };
    if hits.is_empty() || total <= 0.0 {
        map.half_extent_mm = half_extent_mm.unwrap_or(10.0);
        map.irradiance_w_m2 = vec![0.0; n * n];
        return map;
    }
    map.centroid_u_mm = hits.iter().map(|h| h.u as f64 * h.power_mw).sum::<f64>() / total;
    map.centroid_v_mm = hits.iter().map(|h| h.v as f64 * h.power_mw).sum::<f64>() / total;
    let width90 = |f: &dyn Fn(&Hit) -> f64| {
        let mut v: Vec<(f64, f64)> = hits.iter().map(|h| (f(h), h.power_mw)).collect();
        v.sort_by(|a, b| a.0.total_cmp(&b.0));
        let (mut acc, mut lo, mut hi) = (0.0, v[0].0, v[v.len() - 1].0);
        let mut lo_set = false;
        for (x, w) in &v {
            acc += w;
            if !lo_set && acc >= 0.05 * total {
                lo = *x;
                lo_set = true;
            }
            if acc >= 0.95 * total {
                hi = *x;
                break;
            }
        }
        hi - lo
    };
    map.width90_u_mm = width90(&|h| h.u as f64);
    map.width90_v_mm = width90(&|h| h.v as f64);
    let extent = half_extent_mm.unwrap_or_else(|| {
        let mut r: Vec<(f64, f64)> =
            hits.iter().map(|h| ((h.u as f64).abs().max((h.v as f64).abs()), h.power_mw)).collect();
        r.sort_by(|a, b| a.0.total_cmp(&b.0));
        let mut acc = 0.0;
        let mut e = r[r.len() - 1].0;
        for (x, w) in &r {
            acc += w;
            if acc >= 0.998 * total {
                e = *x;
                break;
            }
        }
        (e * 1.08).max(1.0)
    });
    map.half_extent_mm = extent;
    let cell = 2.0 * extent / n as f64;
    let mut grid = vec![0.0; n * n];
    for h in hits {
        let iu = ((h.u as f64 + extent) / cell).floor();
        let iv = ((h.v as f64 + extent) / cell).floor();
        if iu >= 0.0 && iv >= 0.0 && (iu as usize) < n && (iv as usize) < n {
            grid[iv as usize * n + iu as usize] += h.power_mw;
        }
    }
    // mW per cell area (mm²) -> W/m² (×1000).
    let scale = 1000.0 / (cell * cell);
    grid.iter_mut().for_each(|g| *g *= scale);
    map.peak_w_m2 = grid.iter().cloned().fold(0.0, f64::max);
    map.irradiance_w_m2 = grid;
    map
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::physics::fiber::{BareFootprint, transport};

    #[test]
    fn bare_fiber_trace_matches_analytic_footprint() {
        // Normal incidence, no lens, large coupon: radial histogram of hits
        // must reproduce the exact convolution footprint.
        let fiber = FiberSettings::default();
        let t = transport(&fiber);
        let bench = BenchSettings {
            use_lens: false,
            incidence_deg: 0.0,
            source_distance_mm: 63.5,
            coupon: Coupon { length_mm: 500.0, width_mm: 500.0, ..Coupon::default() },
            ..BenchSettings::default()
        };
        let beam = trace_beam(&bench, &fiber, &t, 400_000, 1);
        let total: f64 = beam.hits.iter().map(|h| h.power_mw).sum();
        assert!((total / t.output_mw - 1.0).abs() < 1e-9);
        let model = BareFootprint::new(63.5, 0.5, t.output_mw);
        // Encircled power within several radii.
        for &rr in &[5.0, 15.0, 25.0, 33.0] {
            let inside: f64 =
                beam.hits.iter().filter(|h| (h.u as f64).hypot(h.v as f64) <= rr).map(|h| h.power_mw).sum();
            let gl = crate::physics::math::gauss_legendre(400);
            let analytic: f64 = gl
                .0
                .iter()
                .zip(&gl.1)
                .map(|(&x, &w)| {
                    let r = 0.5 * rr * (x + 1.0);
                    0.5 * rr * w * model.irradiance(r) * TAU * r
                })
                .sum();
            assert!((inside / analytic - 1.0).abs() < 2e-3, "r={rr}: {inside} vs {analytic}");
        }
    }

    #[test]
    fn lens_budget_closes_and_collimates() {
        let fiber = FiberSettings::default();
        let t = transport(&fiber);
        let bench = BenchSettings::default();
        let beam = trace_beam(&bench, &fiber, &t, 200_000, 2);
        let b = beam.budget;
        let sum = b.blocked_by_mount_mw + b.lens_reflection_mw + b.lens_edge_mw + b.lens_tir_mw + b.past_lens_mw;
        assert!((sum / b.fiber_output_mw - 1.0).abs() < 1e-9, "{b:?}");
        // The Ø12.7 mm lens at ~13.4 mm captures roughly sin²(23°)/0.25 of the NA 0.5 cone.
        let captured = (b.past_lens_mw + b.lens_reflection_mw) / b.fiber_output_mw;
        assert!(captured > 0.45 && captured < 0.75, "captured {captured}");
        // Bare-glass loss ≈ 2×4 % at the surfaces.
        let refl = b.lens_reflection_mw / (b.past_lens_mw + b.lens_reflection_mw);
        assert!(refl > 0.07 && refl < 0.12, "reflected {refl}");
    }

    #[test]
    fn small_launch_na_divergence_equals_core_over_f() {
        // Paraxial check of the PDF statement: 1 mm core, 15 mm lens -> ~3.8°.
        let fiber = FiberSettings { launch_na: 0.05, ..FiberSettings::default() };
        let t = transport(&fiber);
        let bench = BenchSettings::default();
        let beam = trace_beam(&bench, &fiber, &t, 100_000, 3);
        let prof = beam_profile(&bench, &beam, 587.6, 100.0);
        // Uniform core: the full angle is 2a/f; 90 % of a uniform disk sits
        // within sqrt(0.9) of the edge radius.
        let full = (2.0 * CORE_RADIUS_MM / bench.lens.efl(587.6)).to_degrees();
        assert!((prof.divergence90_full_deg / (full * 0.9f64.sqrt()) - 1.0).abs() < 0.06, "{} vs {full}", prof.divergence90_full_deg);
    }
}
