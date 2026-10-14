//! End-to-end pipeline: fiber transport → beam trace → sample → scan.

use super::bench::{BeamProfile, BeamResult, BenchSettings, SpotMap, beam_profile, spot_map, trace_beam};
use super::fiber::{self, FiberSettings, PowerMode, Transport};
use super::goniometer::{DetectorSettings, ScanResult, ScanSettings, compute_scan};
use super::math::Vec3;
use super::sample::{BrdfParts, Finish, SampleKind, SampleModel, SampleSettings, Simulant};
use super::wave::WaveSettings;
use serde::{Deserialize, Serialize};
use std::sync::Arc;

pub const APP_VERSION: &str = env!("CARGO_PKG_VERSION");

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize, Default)]
pub enum Quality {
    Draft,
    #[default]
    Standard,
    High,
}

impl Quality {
    /// (rays traced from the fiber, rays used for scan integration)
    pub fn rays(self) -> (usize, usize) {
        match self {
            Quality::Draft => (60_000, 1_500),
            Quality::Standard => (250_000, 4_096),
            Quality::High => (1_000_000, 12_000),
        }
    }
    pub fn label(self) -> &'static str {
        match self {
            Quality::Draft => "Draft",
            Quality::Standard => "Standard",
            Quality::High => "High",
        }
    }
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct SimSettings {
    pub fiber: FiberSettings,
    pub bench: BenchSettings,
    pub sample: SampleSettings,
    pub detector: DetectorSettings,
    pub scan: ScanSettings,
    pub wave: WaveSettings,
    pub quality: Quality,
    pub seed: u64,
}

impl Default for SimSettings {
    fn default() -> Self {
        SimSettings {
            fiber: FiberSettings::default(),
            bench: BenchSettings::default(),
            sample: SampleSettings::default(),
            detector: DetectorSettings::default(),
            scan: ScanSettings::default(),
            wave: WaveSettings::default(),
            quality: Quality::Standard,
            seed: 7,
        }
    }
}

impl SimSettings {
    pub fn validate(&self) -> Result<(), String> {
        self.fiber.validate()?;
        self.bench.validate()?;
        self.sample.composite.validate()?;
        self.detector.validate()?;
        self.scan.validate()?;
        Ok(())
    }

    /// Key of everything the beam trace depends on (for caching).
    pub fn beam_key(&self) -> String {
        serde_json::to_string(&(&self.fiber, &self.bench, self.quality, self.seed)).unwrap_or_default()
    }

    pub fn from_json(text: &str) -> Result<Self, String> {
        let s: SimSettings = serde_json::from_str(text).map_err(|e| e.to_string())?;
        s.validate()?;
        Ok(s)
    }
}

/// In-plane BRDF of the model at the nominal incidence, per component.
#[derive(Clone, Debug, Default, Serialize)]
pub struct BrdfCurve {
    pub view_deg: Vec<f64>,
    pub parts: Vec<BrdfParts>,
    pub albedo: f64,
}

#[derive(Clone, Debug)]
pub struct SimResult {
    pub settings: SimSettings,
    pub transport: Transport,
    pub beam: Arc<BeamResult>,
    pub profile: BeamProfile,
    pub spot: SpotMap,
    pub scan: ScanResult,
    /// Mirror (R = 1) scan with the same bench: the instrument function.
    pub instrument: ScanResult,
    pub brdf: BrdfCurve,
    pub warnings: Vec<String>,
    pub mean_wavelength_nm: f64,
    pub compute_ms: f64,
}

pub fn brdf_curve(model: &SampleModel, incidence_deg: f64) -> BrdfCurve {
    let t = incidence_deg.to_radians();
    let wi = Vec3::new(-t.sin(), 0.0, t.cos());
    let view: Vec<f64> = (0..=356).map(|i| -89.0 + 0.5 * i as f64).collect();
    let parts = view
        .iter()
        .map(|&a| {
            let r = f64::to_radians(a);
            model.parts(wi, Vec3::new(r.sin(), 0.0, r.cos()))
        })
        .collect();
    BrdfCurve { view_deg: view, parts, albedo: model.albedo(t) }
}

/// Run everything. `cached_beam` lets the GUI skip the trace when only the
/// sample, detector or scan changed.
pub fn simulate(settings: &SimSettings, cached_beam: Option<Arc<BeamResult>>) -> Result<SimResult, String> {
    let start = std::time::Instant::now();
    settings.validate()?;
    let t = fiber::transport(&settings.fiber);
    let (n_rays, scan_rays) = settings.quality.rays();
    let beam = match cached_beam {
        Some(b) => b,
        None => Arc::new(trace_beam(&settings.bench, &settings.fiber, &t, n_rays, settings.seed)),
    };
    let mean_wl = t.weighted_mean(|w| w);
    let reach = (settings.bench.source_distance_mm * 1.6).max(40.0);
    let profile = beam_profile(&settings.bench, &beam, mean_wl, reach);
    let spot = spot_map(&beam, 96, None);
    let model = SampleModel::new(&settings.sample);
    let scan_settings = ScanSettings { scan_rays, ..settings.scan.clone() };
    let scan = compute_scan(&settings.bench, &beam, &model, &settings.detector, &scan_settings);
    let mirror = SampleModel::new(&SampleSettings { kind: SampleKind::Mirror, mirror_reflectance: 1.0, ..settings.sample.clone() });
    let instrument = compute_scan(&settings.bench, &beam, &mirror, &settings.detector, &scan_settings);
    let brdf = brdf_curve(&model, settings.bench.incidence_deg);
    let warnings = warnings(settings, &t, &beam, &spot, &scan);
    Ok(SimResult {
        settings: settings.clone(),
        transport: t,
        beam,
        profile,
        spot,
        scan,
        instrument,
        brdf,
        warnings,
        mean_wavelength_nm: mean_wl,
        compute_ms: start.elapsed().as_secs_f64() * 1000.0,
    })
}

fn warnings(s: &SimSettings, t: &Transport, beam: &BeamResult, spot: &SpotMap, scan: &ScanResult) -> Vec<String> {
    let mut w = Vec::new();
    let b = &beam.budget;
    if s.bench.use_lens && b.fiber_output_mw > 0.0 {
        let captured = (b.past_lens_mw + b.lens_reflection_mw) / b.fiber_output_mw;
        if captured < 0.8 {
            w.push(format!(
                "The Ø{:.1} mm lens captures {:.0} % of the fiber's NA {:.2} cone at {:.2} mm; {:.0} % hits the mount.",
                s.bench.lens.diameter_mm.min(s.bench.lens.mount_aperture_mm),
                100.0 * captured,
                t.output_na,
                s.bench.fiber_to_lens_mm,
                100.0 * b.blocked_by_mount_mw / b.fiber_output_mw
            ));
        }
    }
    if b.past_lens_mw > 0.0 && b.off_coupon_mw / b.past_lens_mw > 0.02 {
        w.push(format!(
            "{:.0} % of the beam misses the {:.1} × {:.1} mm coupon (spot {:.1} × {:.1} mm at θi = {:.0}°).",
            100.0 * b.off_coupon_mw / b.past_lens_mw,
            s.bench.coupon.length_mm,
            s.bench.coupon.width_mm,
            spot.width90_u_mm,
            spot.width90_v_mm,
            s.bench.incidence_deg
        ));
    }
    let det_w = s.detector.angular_width_deg();
    if det_w > 2.0 * s.scan.near_halfwidth_deg {
        w.push(format!(
            "The detector spans {det_w:.1}°, wider than the ±{:.1}° near-specular window.",
            s.scan.near_halfwidth_deg
        ));
    }
    let spec = s.bench.specular_angle_deg();
    if (spec - scan.blocked_center_deg).abs() < scan.blocked_halfwidth_deg {
        w.push("The specular direction falls inside the zone blocked by the source arm.".into());
    }
    if spot.centroid_u_mm.hypot(spot.centroid_v_mm) > 1.0 {
        w.push(format!(
            "Beam center is {:.1} mm from the rotation center (fiber offset or tilt).",
            spot.centroid_u_mm.hypot(spot.centroid_v_mm)
        ));
    }
    if s.fiber.power_mode == PowerMode::Assumed {
        w.push(format!(
            "Absolute powers scale with the assumed {:.3} mW input; normalize to a power-meter reading.",
            s.fiber.incident_mw
        ));
    }
    w
}

#[derive(Clone, Debug, Serialize)]
pub struct StudyRow {
    pub simulant: Simulant,
    pub regolith_mass_percent: f64,
    pub finish: Finish,
    pub incidence_deg: f64,
    pub near_mw: f64,
    pub off_mw: f64,
    pub ratio: f64,
    pub peak_angle_deg: f64,
    pub fwhm_deg: f64,
    pub exposed_fraction: f64,
    pub albedo: f64,
}

/// The experiment matrix from the research plan: regolith content ×
/// finish × illumination angle (scans of the same instrument model).
pub fn run_study(
    base: &SimSettings,
    simulants: &[Simulant],
    percents: &[f64],
    incidences: &[f64],
    progress: &(dyn Fn(f64) + Sync),
) -> Result<Vec<StudyRow>, String> {
    base.validate()?;
    let t = fiber::transport(&base.fiber);
    let (n_rays, scan_rays) = base.quality.rays();
    let total = (incidences.len() * simulants.len() * percents.len() * 2).max(1);
    let mut done = 0usize;
    let mut rows = Vec::new();
    for &inc in incidences {
        let bench = BenchSettings { incidence_deg: inc, ..base.bench.clone() };
        bench.validate()?;
        let beam = trace_beam(&bench, &base.fiber, &t, n_rays, base.seed);
        for &sim in simulants {
            for &pct in percents {
                for finish in [Finish::AsCast, Finish::Abraded] {
                    let mut sample = base.sample.clone();
                    sample.kind = SampleKind::Composite;
                    if sample.composite.simulant != sim {
                        sample.composite.apply_simulant(sim);
                    }
                    sample.composite.regolith_mass_percent = pct;
                    sample.composite.finish = finish;
                    let model = SampleModel::new(&sample);
                    let scan_settings = ScanSettings { scan_rays, ..base.scan.clone() };
                    let scan = compute_scan(&bench, &beam, &model, &base.detector, &scan_settings);
                    rows.push(StudyRow {
                        simulant: sim,
                        regolith_mass_percent: pct,
                        finish,
                        incidence_deg: inc,
                        near_mw: scan.metrics.near_mean_mw,
                        off_mw: scan.metrics.off_mean_mw,
                        ratio: scan.metrics.near_to_off_ratio,
                        peak_angle_deg: scan.metrics.peak_angle_deg,
                        fwhm_deg: scan.metrics.fwhm_deg,
                        exposed_fraction: sample.composite.exposed_fraction(),
                        albedo: model.albedo(inc.to_radians()),
                    });
                    done += 1;
                    progress(done as f64 / total as f64);
                }
            }
        }
    }
    Ok(rows)
}

/// Effect sizes in log10(near/off ratio): finishing vs. regolith content.
#[derive(Clone, Debug, Default, Serialize)]
pub struct EffectSummary {
    pub incidence_deg: f64,
    pub simulant: Option<Simulant>,
    /// Mean |log10 ratio(abraded) − log10 ratio(as-cast)| across contents.
    pub finish_effect: f64,
    /// Mean over finishes of the log10-ratio range across regolith contents.
    pub content_effect: f64,
}

pub fn summarize_effects(rows: &[StudyRow]) -> Vec<EffectSummary> {
    let mut keys: Vec<(f64, Simulant)> = Vec::new();
    for r in rows {
        if !keys.iter().any(|k| k.0 == r.incidence_deg && k.1 == r.simulant) {
            keys.push((r.incidence_deg, r.simulant));
        }
    }
    keys.iter()
        .map(|&(inc, sim)| {
            let sel: Vec<&StudyRow> = rows.iter().filter(|r| r.incidence_deg == inc && r.simulant == sim).collect();
            let lr = |r: &StudyRow| r.ratio.max(1e-12).log10();
            let mut finish = Vec::new();
            let mut pcts: Vec<f64> = sel.iter().map(|r| r.regolith_mass_percent).collect();
            pcts.sort_by(|a, b| a.total_cmp(b));
            pcts.dedup();
            for p in &pcts {
                let c = sel.iter().find(|r| r.regolith_mass_percent == *p && r.finish == Finish::AsCast);
                let a = sel.iter().find(|r| r.regolith_mass_percent == *p && r.finish == Finish::Abraded);
                if let (Some(c), Some(a)) = (c, a) {
                    finish.push((lr(a) - lr(c)).abs());
                }
            }
            let mut content = Vec::new();
            for f in [Finish::AsCast, Finish::Abraded] {
                let v: Vec<f64> = sel.iter().filter(|r| r.finish == f && r.regolith_mass_percent > 0.0).map(|r| lr(r)).collect();
                if v.len() >= 2 {
                    let (lo, hi) = v.iter().fold((f64::MAX, f64::MIN), |(a, b), &x| (a.min(x), b.max(x)));
                    content.push(hi - lo);
                }
            }
            let mean = |v: &[f64]| if v.is_empty() { 0.0 } else { v.iter().sum::<f64>() / v.len() as f64 };
            EffectSummary { incidence_deg: inc, simulant: Some(sim), finish_effect: mean(&finish), content_effect: mean(&content) }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_pipeline_runs_and_is_consistent() {
        let s = SimSettings { quality: Quality::Draft, ..Default::default() };
        let r = simulate(&s, None).unwrap();
        assert!(r.scan.metrics.near_to_off_ratio > 1.0);
        assert!(r.transport.output_mw > 0.78 && r.transport.output_mw < 0.79);
        assert!(r.beam.budget.on_coupon_mw > 0.0);
        // Reusing the beam gives identical scans.
        let again = simulate(&s, Some(r.beam.clone())).unwrap();
        assert_eq!(again.scan.power_mw, r.scan.power_mw);
        let json = serde_json::to_string(&s).unwrap();
        assert_eq!(SimSettings::from_json(&json).unwrap(), s);
    }

    #[test]
    fn study_matrix_has_expected_shape() {
        let s = SimSettings { quality: Quality::Draft, ..Default::default() };
        let rows = run_study(&s, &[Simulant::Lhs1], &[0.0, 30.0], &[30.0, 60.0], &|_| {}).unwrap();
        assert_eq!(rows.len(), 2 * 2 * 2);
        // Abrasion lowers the near/off ratio for the composite.
        let cast = rows.iter().find(|r| r.regolith_mass_percent == 30.0 && r.finish == Finish::AsCast && r.incidence_deg == 30.0).unwrap();
        let abr = rows.iter().find(|r| r.regolith_mass_percent == 30.0 && r.finish == Finish::Abraded && r.incidence_deg == 30.0).unwrap();
        assert!(abr.ratio < cast.ratio);
        let eff = summarize_effects(&rows);
        assert_eq!(eff.len(), 2);
    }
}
