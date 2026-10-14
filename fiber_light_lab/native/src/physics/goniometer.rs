//! Detector on the moving arm, the angular scan and the comparison metrics.
//!
//! For each traced ray that lands on the coupon, the reflected power reaching
//! every detector position is integrated over the detector's active area
//! (next-event estimation): dP = P_ray · f_r(ω_i, ω_o) · cosθ_o · cosθ_det ·
//! dA / d². The ideal mirror is handled by reflecting each ray and testing
//! whether it lands inside the active area. Finite source size, lens
//! aberrations, spot size and detector size therefore all broaden the scan,
//! exactly as they would on the real bench.

use super::bench::{BeamResult, BenchSettings, Hit};
use super::data;
use super::lens::reflect;
use super::math::{Curve, PI, Vec3};
use super::sample::SampleModel;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize, Default)]
pub enum DetectorPreset {
    /// Thorlabs SM05PD1A / PDA36A2 class: Si, 3.6 × 3.6 mm (13 mm²).
    #[default]
    ThorlabsSi36,
    /// Thorlabs S120C power-meter head: Si, Ø9.5 mm aperture.
    ThorlabsS120C,
    /// Vishay BPW34: Si, 7.5 mm² radiant-sensitive area.
    Bpw34,
    Custom,
}

impl DetectorPreset {
    pub fn label(self) -> &'static str {
        match self {
            DetectorPreset::ThorlabsSi36 => "Thorlabs SM05PD1A / PDA36A2 — Si 3.6 × 3.6 mm",
            DetectorPreset::ThorlabsS120C => "Thorlabs S120C — Si Ø9.5 mm",
            DetectorPreset::Bpw34 => "Vishay BPW34 — Si 7.5 mm²",
            DetectorPreset::Custom => "Custom",
        }
    }
    pub const ALL: [DetectorPreset; 4] =
        [DetectorPreset::ThorlabsSi36, DetectorPreset::ThorlabsS120C, DetectorPreset::Bpw34, DetectorPreset::Custom];
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize, Default)]
pub enum DetectorShape {
    #[default]
    Square,
    Circle,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct DetectorSettings {
    pub preset: DetectorPreset,
    pub shape: DetectorShape,
    /// Side (square) or diameter (circle) of the active area (mm).
    pub size_mm: f64,
    /// Rotation center to detector active surface (mm).
    pub arm_radius_mm: f64,
    /// Outer radius of the detector housing, for arm collisions (mm).
    pub housing_radius_mm: f64,
    /// Responsivity curve (A/W). Default: typical Si photodiode shape.
    pub responsivity: Curve,
    pub responsivity_label: String,
    /// Wavelength a power meter is set to (it divides current by R(λset)).
    pub meter_wavelength_nm: f64,
    /// Transimpedance / load for a voltage estimate (V/A).
    pub transimpedance_v_per_a: f64,
}

impl Default for DetectorSettings {
    fn default() -> Self {
        let mut d = DetectorSettings {
            preset: DetectorPreset::ThorlabsSi36,
            shape: DetectorShape::Square,
            size_mm: 3.6,
            arm_radius_mm: 100.0,
            housing_radius_mm: 9.0,
            responsivity: data::si_responsivity().clone(),
            responsivity_label: "Typical Si photodiode (approximate)".into(),
            meter_wavelength_nm: 650.0,
            transimpedance_v_per_a: 1.0e6,
        };
        d.apply_preset(DetectorPreset::ThorlabsSi36);
        d
    }
}

impl DetectorSettings {
    pub fn apply_preset(&mut self, p: DetectorPreset) {
        self.preset = p;
        match p {
            DetectorPreset::ThorlabsSi36 => {
                self.shape = DetectorShape::Square;
                self.size_mm = 3.6;
                self.housing_radius_mm = 9.0;
            }
            DetectorPreset::ThorlabsS120C => {
                self.shape = DetectorShape::Circle;
                self.size_mm = 9.5;
                self.housing_radius_mm = 15.0;
            }
            DetectorPreset::Bpw34 => {
                self.shape = DetectorShape::Square;
                self.size_mm = 7.5f64.sqrt();
                self.housing_radius_mm = 3.0;
            }
            DetectorPreset::Custom => {}
        }
    }

    pub fn area_mm2(&self) -> f64 {
        match self.shape {
            DetectorShape::Square => self.size_mm * self.size_mm,
            DetectorShape::Circle => PI * 0.25 * self.size_mm * self.size_mm,
        }
    }

    /// Solid angle subtended at the rotation center (sr).
    pub fn solid_angle_sr(&self) -> f64 {
        self.area_mm2() / (self.arm_radius_mm * self.arm_radius_mm)
    }

    /// Full angular width seen from the rotation center (deg).
    pub fn angular_width_deg(&self) -> f64 {
        (2.0 * (0.5 * self.size_mm / self.arm_radius_mm).atan()).to_degrees()
    }

    pub fn validate(&self) -> Result<(), String> {
        if self.size_mm <= 0.0 || self.arm_radius_mm <= self.size_mm {
            return Err("Detector size must be positive and smaller than the arm radius.".into());
        }
        Ok(())
    }

    /// Stratified points on the active area in detector coordinates (mm),
    /// with the area element of each.
    fn area_points(&self, n: usize) -> (Vec<(f64, f64)>, f64) {
        let mut pts = Vec::with_capacity(n * n);
        let h = 0.5 * self.size_mm;
        for i in 0..n {
            for j in 0..n {
                let a = 2.0 * (i as f64 + 0.5) / n as f64 - 1.0;
                let b = 2.0 * (j as f64 + 0.5) / n as f64 - 1.0;
                match self.shape {
                    DetectorShape::Square => pts.push((a * h, b * h)),
                    DetectorShape::Circle => {
                        // Shirley–Chiu concentric map: equal-area square -> disk.
                        let (r, phi) = if a.abs() > b.abs() {
                            (a, PI / 4.0 * (b / a))
                        } else {
                            (b, PI / 2.0 - PI / 4.0 * (a / b))
                        };
                        pts.push((h * r * phi.cos(), h * r * phi.sin()));
                    }
                }
            }
        }
        (pts, self.area_mm2() / (n * n) as f64)
    }

    fn inside(&self, a: f64, b: f64) -> bool {
        let h = 0.5 * self.size_mm;
        match self.shape {
            DetectorShape::Square => a.abs() <= h && b.abs() <= h,
            DetectorShape::Circle => a.hypot(b) <= h,
        }
    }
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct ScanSettings {
    pub start_deg: f64,
    pub stop_deg: f64,
    pub step_deg: f64,
    /// Near-specular window: |θ − θspec| ≤ this (deg).
    pub near_halfwidth_deg: f64,
    /// Away-from-specular window: |θ − θspec| ≥ this (deg).
    pub off_min_deg: f64,
    /// Rays (of those on the coupon) used for the scan integration.
    pub scan_rays: usize,
}

impl Default for ScanSettings {
    fn default() -> Self {
        ScanSettings { start_deg: -85.0, stop_deg: 85.0, step_deg: 1.0, near_halfwidth_deg: 5.0, off_min_deg: 20.0, scan_rays: 4096 }
    }
}

impl ScanSettings {
    pub fn angles(&self) -> Vec<f64> {
        let step = self.step_deg.max(0.01);
        let n = ((self.stop_deg - self.start_deg) / step).floor() as usize + 1;
        (0..n).map(|i| self.start_deg + i as f64 * step).collect()
    }
    pub fn validate(&self) -> Result<(), String> {
        if self.stop_deg <= self.start_deg || self.step_deg <= 0.0 {
            return Err("Scan range must increase with a positive step.".into());
        }
        if self.start_deg < -89.9 || self.stop_deg > 89.9 {
            return Err("Detector angles must stay within ±89.9°.".into());
        }
        if self.angles().len() > 20_000 {
            return Err("Too many scan positions (reduce range or increase step).".into());
        }
        Ok(())
    }
}

#[derive(Clone, Debug, Default, Serialize)]
pub struct ScanMetrics {
    pub specular_angle_deg: f64,
    pub near_mean_mw: f64,
    pub off_mean_mw: f64,
    pub near_to_off_ratio: f64,
    pub peak_mw: f64,
    pub peak_angle_deg: f64,
    pub fwhm_deg: f64,
    pub near_points: usize,
    pub off_points: usize,
}

#[derive(Clone, Debug, Default, Serialize)]
pub struct ScanResult {
    pub angle_deg: Vec<f64>,
    /// Optical power on the active area (mW).
    pub power_mw: Vec<f64>,
    /// Photocurrent (mA) using the responsivity curve.
    pub current_ma: Vec<f64>,
    /// Reading of a power meter set to `meter_wavelength_nm` (mW).
    pub meter_mw: Vec<f64>,
    /// BRDF an experimenter would infer: P/(P_on_coupon·Ω·cosθd) (sr⁻¹).
    pub brdf_measured: Vec<f64>,
    /// Model BRDF at the nominal geometry (no instrument effects) (sr⁻¹).
    pub brdf_true: Vec<f64>,
    /// Detector or arm collides with the source arm / beam here.
    pub blocked: Vec<bool>,
    pub blocked_center_deg: f64,
    pub blocked_halfwidth_deg: f64,
    pub metrics: ScanMetrics,
    pub incident_on_coupon_mw: f64,
    pub detector_subsamples: usize,
}

/// Angular half-width (deg) of the zone where the detector cannot go.
pub fn blocked_zone(bench: &BenchSettings, det: &DetectorSettings) -> (f64, f64) {
    let center = -bench.incidence_deg;
    let src_d = if bench.use_lens {
        (bench.source_distance_mm - 0.5 * bench.lens.center_thickness_mm).max(1.0)
    } else {
        bench.source_distance_mm.max(1.0)
    };
    let half = (bench.source_housing_radius_mm / src_d).atan() + (det.housing_radius_mm / det.arm_radius_mm).atan();
    (center, half.to_degrees())
}

/// Pick a power-preserving, stratified subset of on-coupon hits.
fn scan_subset(beam: &BeamResult, n: usize) -> (Vec<Hit>, f64) {
    let on: Vec<&Hit> = beam.hits.iter().filter(|h| h.on_coupon).collect();
    let total: f64 = on.iter().map(|h| h.power_mw).sum();
    if on.is_empty() {
        return (Vec::new(), 0.0);
    }
    let step = (on.len() as f64 / n.max(1) as f64).max(1.0);
    let mut sub = Vec::with_capacity(n);
    let mut x = 0.0;
    while (x as usize) < on.len() {
        sub.push(*on[x as usize]);
        x += step;
    }
    let sub_total: f64 = sub.iter().map(|h| h.power_mw).sum();
    let k = if sub_total > 0.0 { total / sub_total } else { 0.0 };
    sub.iter_mut().for_each(|h| h.power_mw *= k);
    (sub, total)
}

pub fn compute_scan(
    bench: &BenchSettings,
    beam: &BeamResult,
    model: &SampleModel,
    det: &DetectorSettings,
    scan: &ScanSettings,
) -> ScanResult {
    let angles = scan.angles();
    // The mirror path is cheap per ray, so it uses many more rays.
    let n_hits = if model.is_delta() { scan.scan_rays.max(60_000) } else { scan.scan_rays };
    let (hits, on_coupon) = scan_subset(beam, n_hits);
    let sframe = bench.sample_frame();
    let normal = sframe.w;
    let r_det = det.arm_radius_mm;
    // Detector sampling density follows the narrowest BRDF lobe.
    let det_ang = det.size_mm / r_det;
    let n_sub = if model.is_delta() {
        1
    } else {
        ((det_ang / (0.35 * model.narrowest_lobe().max(1e-3))).ceil() as usize).clamp(2, 10)
    };
    let (area_pts, d_area) = det.area_points(n_sub);
    let resp: Vec<f64> = hits.iter().map(|h| det.responsivity.eval(h.wavelength_nm as f64)).collect();
    let local_in: Vec<Vec3> = hits.iter().map(|h| sframe.to_local(-h.direction)).collect();

    let per_angle: Vec<(f64, f64)> = angles
        .par_iter()
        .map(|&theta| {
            let t = theta.to_radians();
            let c = Vec3::new(t.sin(), 0.0, t.cos()) * r_det;
            let n_d = -(c / r_det); // faces the sample
            let tu = Vec3::new(t.cos(), 0.0, -t.sin());
            let tv = Vec3::Y;
            let (mut p_sum, mut i_sum) = (0.0, 0.0);
            if model.is_delta() {
                for (k, h) in hits.iter().enumerate() {
                    let r = reflect(h.direction, normal);
                    let denom = n_d.dot(r);
                    if denom >= 0.0 {
                        continue;
                    }
                    let tt = n_d.dot(c - h.position) / denom;
                    if tt <= 0.0 {
                        continue;
                    }
                    let q = h.position + r * tt - c;
                    if det.inside(q.dot(tu), q.dot(tv)) {
                        let w = h.power_mw * model.mirror_reflectance;
                        p_sum += w;
                        i_sum += w * resp[k];
                    }
                }
            } else {
                let world_pts: Vec<Vec3> = area_pts.iter().map(|&(a, b)| c + tu * a + tv * b).collect();
                for (k, h) in hits.iter().enumerate() {
                    let wi = local_in[k];
                    if wi.z <= 0.0 {
                        continue;
                    }
                    let mut acc = 0.0;
                    for q in &world_pts {
                        let d = *q - h.position;
                        let dist2 = d.norm2();
                        let wo_world = d / dist2.sqrt();
                        let cos_det = -wo_world.dot(n_d);
                        if cos_det <= 0.0 {
                            continue;
                        }
                        let wo = sframe.to_local(wo_world);
                        if wo.z <= 0.0 {
                            continue;
                        }
                        acc += model.brdf(wi, wo) * wo.z * cos_det / dist2;
                    }
                    let w = h.power_mw * acc * d_area;
                    p_sum += w;
                    i_sum += w * resp[k];
                }
            }
            (p_sum, i_sum)
        })
        .collect();

    let (block_c, block_hw) = blocked_zone(bench, det);
    let omega = det.solid_angle_sr();
    let r_meter = det.responsivity.eval(det.meter_wavelength_nm).max(1e-12);
    let wi_nom = sframe.to_local(bench.source_direction());
    let mut res = ScanResult {
        blocked_center_deg: block_c,
        blocked_halfwidth_deg: block_hw,
        incident_on_coupon_mw: on_coupon,
        detector_subsamples: n_sub * n_sub,
        ..Default::default()
    };
    for (&theta, &(p, i)) in angles.iter().zip(&per_angle) {
        let t = theta.to_radians();
        res.angle_deg.push(theta);
        res.power_mw.push(p);
        res.current_ma.push(i);
        res.meter_mw.push(i / r_meter);
        let denom = on_coupon * omega * t.cos();
        res.brdf_measured.push(if denom > 0.0 { p / denom } else { 0.0 });
        let wo_nom = sframe.to_local(Vec3::new(t.sin(), 0.0, t.cos()));
        res.brdf_true.push(if model.is_delta() { f64::NAN } else { model.brdf(wi_nom, wo_nom) });
        res.blocked.push((theta - block_c).abs() < block_hw);
    }
    res.metrics = metrics(&res, bench.specular_angle_deg(), scan);
    res
}

pub fn metrics(res: &ScanResult, spec: f64, scan: &ScanSettings) -> ScanMetrics {
    let mut m = ScanMetrics { specular_angle_deg: spec, ..Default::default() };
    let usable = |i: usize| !res.blocked[i];
    let (mut near, mut off) = (Vec::new(), Vec::new());
    for (i, &a) in res.angle_deg.iter().enumerate() {
        if !usable(i) {
            continue;
        }
        let d = (a - spec).abs();
        if d <= scan.near_halfwidth_deg {
            near.push(res.power_mw[i]);
        } else if d >= scan.off_min_deg {
            off.push(res.power_mw[i]);
        }
    }
    let mean = |v: &[f64]| if v.is_empty() { 0.0 } else { v.iter().sum::<f64>() / v.len() as f64 };
    m.near_mean_mw = mean(&near);
    m.off_mean_mw = mean(&off);
    m.near_points = near.len();
    m.off_points = off.len();
    m.near_to_off_ratio = if m.off_mean_mw > 0.0 { m.near_mean_mw / m.off_mean_mw } else { f64::INFINITY };
    // Peak and FWHM above the median off-specular baseline.
    let mut peak_i = None;
    for i in 0..res.angle_deg.len() {
        if usable(i) && peak_i.is_none_or(|j: usize| res.power_mw[i] > res.power_mw[j]) {
            peak_i = Some(i);
        }
    }
    if let Some(pi) = peak_i {
        m.peak_mw = res.power_mw[pi];
        m.peak_angle_deg = res.angle_deg[pi];
        let mut sorted = off.clone();
        sorted.sort_by(|a, b| a.total_cmp(b));
        let base = if sorted.is_empty() { 0.0 } else { sorted[sorted.len() / 2] };
        let half = base + 0.5 * (m.peak_mw - base);
        let cross = |dir: isize| -> Option<f64> {
            let mut i = pi as isize;
            loop {
                let j = i + dir;
                if j < 0 || j as usize >= res.angle_deg.len() || !usable(j as usize) {
                    return None;
                }
                let (pa, pb) = (res.power_mw[i as usize], res.power_mw[j as usize]);
                if pb <= half {
                    let f = (pa - half) / (pa - pb).max(1e-300);
                    let (aa, ab) = (res.angle_deg[i as usize], res.angle_deg[j as usize]);
                    return Some(aa + f * (ab - aa));
                }
                i = j;
            }
        };
        if let (Some(lo), Some(hi)) = (cross(-1), cross(1)) {
            m.fwhm_deg = (hi - lo).abs();
        }
    }
    m
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::physics::bench::{Coupon, trace_beam};
    use crate::physics::fiber::{FiberSettings, transport};
    use crate::physics::sample::{SampleKind, SampleSettings};

    fn setup(use_lens: bool, incidence: f64) -> (BenchSettings, BeamResult) {
        let fiber = FiberSettings::default();
        let t = transport(&fiber);
        let bench = BenchSettings {
            use_lens,
            incidence_deg: incidence,
            coupon: Coupon { length_mm: 200.0, width_mm: 200.0, ..Coupon::default() },
            ..BenchSettings::default()
        };
        let beam = trace_beam(&bench, &fiber, &t, 60_000, 11);
        (bench, beam)
    }

    #[test]
    fn lambertian_scan_follows_cosine_law() {
        let (bench, beam) = setup(true, 30.0);
        let model = SampleModel::new(&SampleSettings { kind: SampleKind::LambertianReference, ..Default::default() });
        let det = DetectorSettings::default();
        let scan = ScanSettings { start_deg: -60.0, stop_deg: 70.0, step_deg: 10.0, scan_rays: 8000, ..Default::default() };
        let r = compute_scan(&bench, &beam, &model, &det, &scan);
        for (i, &a) in r.angle_deg.iter().enumerate() {
            let expect = r.incident_on_coupon_mw * 0.99 / PI * a.to_radians().cos() * det.solid_angle_sr();
            // Finite spot (~15 mm) at 100 mm changes distances by a few %.
            assert!((r.power_mw[i] / expect - 1.0).abs() < 0.05, "θ={a}: {} vs {expect}", r.power_mw[i]);
            assert!((r.brdf_measured[i] * PI / 0.99 - 1.0).abs() < 0.05);
        }
    }

    #[test]
    fn mirror_scan_peaks_at_specular_and_conserves_power() {
        let (bench, beam) = setup(true, 30.0);
        let model = SampleModel::new(&SampleSettings { kind: SampleKind::Mirror, ..Default::default() });
        // A big detector catches the whole reflected beam.
        let mut det = DetectorSettings::default();
        det.apply_preset(DetectorPreset::Custom);
        det.size_mm = 60.0;
        let scan = ScanSettings { start_deg: 20.0, stop_deg: 40.0, step_deg: 0.5, scan_rays: 20_000, ..Default::default() };
        let r = compute_scan(&bench, &beam, &model, &det, &scan);
        // The big detector sees a flat plateau centered on the specular angle.
        let at30 = r.power_mw[r.angle_deg.iter().position(|&a| a == 30.0).unwrap()];
        assert!((at30 / (0.9 * r.incident_on_coupon_mw) - 1.0).abs() < 1e-9, "{at30}");
        assert!(r.metrics.peak_mw <= at30 * (1.0 + 1e-12));
    }

    #[test]
    fn mirror_instrument_function_is_degrees_wide() {
        // Collimated LB1092 beam (3.8° source-size divergence) and a 3.6 mm
        // detector at 100 mm (2.1°): the mirror scan must be a few degrees wide.
        let (bench, beam) = setup(true, 30.0);
        let model = SampleModel::new(&SampleSettings { kind: SampleKind::Mirror, ..Default::default() });
        let det = DetectorSettings::default();
        let scan = ScanSettings { start_deg: 10.0, stop_deg: 50.0, step_deg: 0.25, scan_rays: 40_000, ..Default::default() };
        let r = compute_scan(&bench, &beam, &model, &det, &scan);
        assert!(r.metrics.fwhm_deg > 3.0 && r.metrics.fwhm_deg < 12.0, "{}", r.metrics.fwhm_deg);
    }

    #[test]
    fn glossy_limit_approaches_mirror() {
        // A very smooth resin lobe should look like a mirror with R = F(θ).
        use crate::physics::sample::{CompositeSettings, dielectric_fresnel};
        let (bench, beam) = setup(true, 30.0);
        let mut c = CompositeSettings { regolith_mass_percent: 30.0, slope_as_cast_deg: 0.05, ..Default::default() };
        c.embedded_darkening = 0.0;
        let glossy = SampleModel::new(&SampleSettings { composite: c, ..Default::default() });
        let f = dielectric_fresnel(1.55, 30f64.to_radians().cos(), false);
        let mirror = SampleModel::new(&SampleSettings { kind: SampleKind::Mirror, mirror_reflectance: f, ..Default::default() });
        let det = DetectorSettings::default();
        let scan = ScanSettings { start_deg: 20.0, stop_deg: 40.0, step_deg: 1.0, scan_rays: 20_000, ..Default::default() };
        let a = compute_scan(&bench, &beam, &glossy, &det, &scan);
        let b = compute_scan(&bench, &beam, &mirror, &det, &scan);
        let sa: f64 = a.power_mw.iter().sum();
        let sb: f64 = b.power_mw.iter().sum();
        assert!((sa / sb - 1.0).abs() < 0.04, "{sa} vs {sb}");
    }

    #[test]
    fn blocked_zone_surrounds_the_source() {
        let (bench, _) = setup(true, 30.0);
        let det = DetectorSettings::default();
        let (c, hw) = blocked_zone(&bench, &det);
        assert_eq!(c, -30.0);
        assert!(hw > 5.0 && hw < 20.0, "{hw}");
    }
}
