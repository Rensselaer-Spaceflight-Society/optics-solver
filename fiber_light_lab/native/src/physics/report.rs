//! CSV / JSON export shared by the desktop app and the command-line runner.

use super::sim::{APP_VERSION, SimResult, StudyRow, summarize_effects};
use serde_json::json;
use std::fmt::Write as _;
use std::fs;
use std::path::Path;

pub const MODEL_NOTES: &str = include_str!("../../MODEL_NATIVE.md");

pub fn summary_json(r: &SimResult) -> serde_json::Value {
    let t = &r.transport;
    let b = &r.beam.budget;
    let s = &r.settings;
    let lens = &s.bench.lens;
    json!({
        "app": "Fiber Light Lab (native)",
        "version": APP_VERSION,
        "power_status": match s.fiber.power_mode {
            super::fiber::PowerMode::Assumed => "Assumed input power; not calibrated",
            super::fiber::PowerMode::Mi152Rating => "Estimated from MI-152 rating (350,000 fc); not calibrated",
        },
        "fiber": {
            "incident_mw": t.incident_mw,
            "output_mw": t.output_mw,
            "guided_mw": t.guided_mw,
            "output_na": t.output_na,
            "half_angle_deg": t.half_angle_deg,
            "mean_wavelength_nm": r.mean_wavelength_nm,
        },
        "lens": if s.bench.use_lens { json!({
            "part": lens.part,
            "efl_mm_at_587_6nm": lens.efl(587.6),
            "bfl_mm_at_587_6nm": lens.bfl(587.6),
            "efl_mm_at_mean_wavelength": lens.efl(r.mean_wavelength_nm),
            "fiber_to_front_vertex_mm": s.bench.fiber_to_lens_mm,
            "coating": lens.coating.label(),
        }) } else { json!(null) },
        "power_budget_mw": b,
        "beam": {
            "divergence90_full_deg": r.profile.divergence90_full_deg,
            "source_size_estimate_deg": r.profile.divergence_estimate_deg,
            "spot_width90_in_plane_mm": r.spot.width90_u_mm,
            "spot_width90_out_of_plane_mm": r.spot.width90_v_mm,
            "spot_centroid_mm": [r.spot.centroid_u_mm, r.spot.centroid_v_mm],
            "peak_irradiance_w_m2": r.spot.peak_w_m2,
        },
        "detector": {
            "label": s.detector.preset.label(),
            "area_mm2": s.detector.area_mm2(),
            "arm_radius_mm": s.detector.arm_radius_mm,
            "angular_width_deg": s.detector.angular_width_deg(),
            "solid_angle_sr": s.detector.solid_angle_sr(),
            "responsivity": s.detector.responsivity_label,
        },
        "scan_metrics": r.scan.metrics,
        "instrument_function": {
            "mirror_fwhm_deg": r.instrument.metrics.fwhm_deg,
            "mirror_peak_angle_deg": r.instrument.metrics.peak_angle_deg,
        },
        "sample": {
            "albedo_at_incidence": r.brdf.albedo,
            "regolith_volume_fraction": s.sample.composite.volume_fraction(),
            "exposed_grain_fraction": s.sample.composite.exposed_fraction(),
            "body_reflectance": s.sample.composite.body_reflectance(),
        },
        "blocked_zone_deg": [r.scan.blocked_center_deg - r.scan.blocked_halfwidth_deg,
                              r.scan.blocked_center_deg + r.scan.blocked_halfwidth_deg],
        "warnings": r.warnings,
        "limitations": "Model results for the stated assumptions; not measurements. See MODEL_NATIVE.md.",
    })
}

pub fn scan_csv(r: &SimResult) -> String {
    let mut out = String::from(
        "detector_angle_deg,power_W,photocurrent_A,meter_reading_W,voltage_V,brdf_measured_per_sr,brdf_model_per_sr,instrument_mirror_W,blocked\n",
    );
    let sc = &r.scan;
    let g = r.settings.detector.transimpedance_v_per_a;
    for i in 0..sc.angle_deg.len() {
        let _ = writeln!(
            out,
            "{},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{},{:.6e},{}",
            sc.angle_deg[i],
            sc.power_mw[i] * 1e-3,
            sc.current_ma[i] * 1e-3,
            sc.meter_mw[i] * 1e-3,
            sc.current_ma[i] * 1e-3 * g,
            sc.brdf_measured[i],
            if sc.brdf_true[i].is_finite() { format!("{:.6e}", sc.brdf_true[i]) } else { String::new() },
            r.instrument.power_mw.get(i).copied().unwrap_or(0.0) * 1e-3,
            sc.blocked[i] as u8
        );
    }
    out
}

pub fn spectrum_csv(r: &SimResult) -> String {
    let t = &r.transport;
    let mut out = String::from("wavelength_nm,incident_mW_per_nm,guided_mW_per_nm,output_mW_per_nm,loss_dB_per_m\n");
    for i in 0..t.wavelength_nm.len() {
        let _ = writeln!(
            out,
            "{},{:.6e},{:.6e},{:.6e},{:.6e}",
            t.wavelength_nm[i], t.incident_spectrum[i], t.guided_spectrum[i], t.output_spectrum[i], t.loss_db_per_m[i]
        );
    }
    out
}

pub fn beam_profile_csv(r: &SimResult) -> String {
    let p = &r.profile;
    let mut out = String::from("distance_from_lens_mm,second_moment_diameter_mm,d90_mm\n");
    for i in 0..p.distance_mm.len() {
        let _ = writeln!(out, "{:.3},{:.4},{:.4}", p.distance_mm[i], p.d_second_moment_mm[i], p.d90_mm[i]);
    }
    out
}

pub fn spot_csv(r: &SimResult) -> String {
    let m = &r.spot;
    let mut out = String::from("u_in_plane_mm,v_rotation_axis_mm,irradiance_W_per_m2\n");
    let cell = 2.0 * m.half_extent_mm / m.n as f64;
    for iv in 0..m.n {
        for iu in 0..m.n {
            let u = -m.half_extent_mm + (iu as f64 + 0.5) * cell;
            let v = -m.half_extent_mm + (iv as f64 + 0.5) * cell;
            let _ = writeln!(out, "{u:.4},{v:.4},{:.6e}", m.irradiance_w_m2[iv * m.n + iu]);
        }
    }
    out
}

pub fn brdf_csv(r: &SimResult) -> String {
    let mut out = String::from("view_angle_deg,resin_specular,body,exposed_grains,back_face,total_per_sr\n");
    for (a, p) in r.brdf.view_deg.iter().zip(&r.brdf.parts) {
        let _ = writeln!(
            out,
            "{a},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e}",
            p.resin_specular,
            p.body,
            p.grains,
            p.back_face,
            p.total()
        );
    }
    out
}

pub fn study_csv(rows: &[StudyRow]) -> String {
    let mut out = String::from(
        "simulant,regolith_mass_percent,finish,incidence_deg,near_specular_W,off_specular_W,near_to_off_ratio,peak_angle_deg,fwhm_deg,exposed_grain_fraction,albedo\n",
    );
    for r in rows {
        let _ = writeln!(
            out,
            "{:?},{},{:?},{},{:.6e},{:.6e},{:.6e},{:.2},{:.2},{:.4},{:.4}",
            r.simulant,
            r.regolith_mass_percent,
            r.finish,
            r.incidence_deg,
            r.near_mw * 1e-3,
            r.off_mw * 1e-3,
            r.ratio,
            r.peak_angle_deg,
            r.fwhm_deg,
            r.exposed_fraction,
            r.albedo
        );
    }
    out
}

/// Write a complete result folder. Returns the list of files written.
pub fn write_folder(dir: &Path, r: &SimResult, study: Option<&[StudyRow]>) -> std::io::Result<Vec<String>> {
    fs::create_dir_all(dir)?;
    let mut files = vec![
        ("settings.json".to_string(), serde_json::to_string_pretty(&r.settings).unwrap_or_default()),
        ("summary.json".to_string(), serde_json::to_string_pretty(&summary_json(r)).unwrap_or_default()),
        ("scan.csv".to_string(), scan_csv(r)),
        ("spectrum.csv".to_string(), spectrum_csv(r)),
        ("beam_profile.csv".to_string(), beam_profile_csv(r)),
        ("sample_irradiance.csv".to_string(), spot_csv(r)),
        ("sample_brdf_in_plane.csv".to_string(), brdf_csv(r)),
        ("MODEL_NATIVE.md".to_string(), MODEL_NOTES.to_string()),
    ];
    if let Some(rows) = study {
        files.push(("study_matrix.csv".into(), study_csv(rows)));
        files.push((
            "study_effects.json".into(),
            serde_json::to_string_pretty(&summarize_effects(rows)).unwrap_or_default(),
        ));
    }
    let mut names = Vec::new();
    for (name, body) in files {
        fs::write(dir.join(&name), body)?;
        names.push(name);
    }
    Ok(names)
}
