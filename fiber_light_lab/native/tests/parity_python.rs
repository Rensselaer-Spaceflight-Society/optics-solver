//! The native transport must reproduce the original Python model.
//! Fixture: tests/fixtures/python_reference.json from
//! `python verification/make_native_reference.py` (numpy/scipy).

use fiber_light_lab::physics::fiber::{
    self, Band, BareFootprint, FiberSettings, fresnel_transmission, spatial_coupling,
};
use serde_json::Value;

fn fixture() -> Value {
    let text = include_str!("fixtures/python_reference.json");
    serde_json::from_str(text).unwrap()
}

fn settings_from(v: &Value) -> FiberSettings {
    FiberSettings {
        length_m: v["length_m"].as_f64().unwrap(),
        incident_mw: v["incident_mw"].as_f64().unwrap(),
        launch_na: v["launch_na"].as_f64().unwrap(),
        spot_radius_mm: v["spot_radius_mm"].as_f64().unwrap(),
        offset_mm: v["offset_mm"].as_f64().unwrap(),
        extra_loss_db: v["extra_loss_db"].as_f64().unwrap(),
        band: match v["band"].as_str().unwrap() {
            "400–700 nm" => Band::Nm400To700,
            _ => Band::Nm400To1000,
        },
        ..FiberSettings::default()
    }
}

fn close(a: f64, b: f64, rel: f64, what: &str) {
    let scale = a.abs().max(b.abs()).max(1e-300);
    assert!((a - b).abs() / scale <= rel, "{what}: native {a} vs python {b}");
}

#[test]
fn transport_matches_python_for_every_case() {
    let fx = fixture();
    for (name, case) in fx["transport"].as_object().unwrap() {
        let s = settings_from(&case["settings"]);
        let t = fiber::transport(&s);
        for key in [
            "output_mw",
            "guided_mw",
            "coupling",
            "spatial_fraction",
            "angular_fraction",
            "face_pair_transmission",
            "half_angle_deg",
        ] {
            let native = match key {
                "output_mw" => t.output_mw,
                "guided_mw" => t.guided_mw,
                "coupling" => t.coupling,
                "spatial_fraction" => t.spatial_fraction,
                "angular_fraction" => t.angular_fraction,
                "face_pair_transmission" => t.face_pair_transmission,
                _ => t.half_angle_deg,
            };
            close(native, case[key].as_f64().unwrap(), 1e-12, &format!("{name}.{key}"));
        }
        for sample in case["samples"].as_array().unwrap() {
            let w = sample["wavelength_nm"].as_f64().unwrap();
            let i = t.wavelength_nm.iter().position(|&x| x == w).unwrap();
            close(t.incident_spectrum[i], sample["incident"].as_f64().unwrap(), 1e-12, "incident");
            close(t.guided_spectrum[i], sample["guided"].as_f64().unwrap(), 1e-12, "guided");
            close(t.output_spectrum[i], sample["output"].as_f64().unwrap(), 1e-12, "output");
            close(t.loss_db_per_m[i], sample["loss_db_per_m"].as_f64().unwrap(), 1e-12, "loss");
        }
    }
}

#[test]
fn default_output_power_is_the_published_value() {
    // verification/numerical_checks.json: output_mw = 0.7843244155491175
    let t = fiber::transport(&FiberSettings::default());
    close(t.output_mw, 0.7843244155491175, 1e-13, "default output");
}

#[test]
fn spatial_coupling_matches_scipy_ncx2() {
    for c in fixture()["spatial_coupling"].as_array().unwrap() {
        let (w, d) = (c["radius_mm"].as_f64().unwrap(), c["offset_mm"].as_f64().unwrap());
        let expect = c["fraction"].as_f64().unwrap();
        let got = spatial_coupling(w, d);
        assert!((got - expect).abs() < 1e-12, "w={w} d={d}: {got} vs {expect}");
    }
}

#[test]
fn fresnel_matches() {
    for c in fixture()["fresnel"].as_array().unwrap() {
        let s = c["sin_external"].as_f64().unwrap();
        close(fresnel_transmission(s), c["transmission"].as_f64().unwrap(), 1e-14, "fresnel");
    }
}

#[test]
fn analytic_footprint_matches_python_fft_grid() {
    // Python uses a 1025² FFT convolution grid; the native model integrates
    // the same convolution exactly in radius. Agreement is limited by the grid.
    for f in fixture()["footprint"].as_array().unwrap() {
        let z = f["distance_mm"].as_f64().unwrap();
        let na = f["launch_na"].as_f64().unwrap();
        let p = f["output_mw"].as_f64().unwrap();
        let model = BareFootprint::new(z, na.min(0.5), p);
        let peak = f["peak_mw_per_mm2"].as_f64().unwrap();
        close(model.irradiance(0.0), peak, 3e-3, &format!("peak z={z}"));
        close(model.envelope_radius_mm(), f["envelope_radius_mm"].as_f64().unwrap(), 1e-12, "envelope");
        // The pixel-sorted encircled radius in Python is biased by ~half a pixel.
        let px = 2.24 * f["envelope_radius_mm"].as_f64().unwrap() / (f["grid"].as_f64().unwrap() - 1.0);
        let r90 = f["encircled90_radius_mm"].as_f64().unwrap();
        assert!((model.encircled_radius(0.9) - r90).abs() < px.max(1e-3 * r90), "r90 z={z}");
        let r50 = f["encircled50_radius_mm"].as_f64().unwrap();
        assert!((model.encircled_radius(0.5) - r50).abs() < px.max(1e-3 * r50), "r50 z={z}");
        // Radial profile away from the cone edge.
        let radii = f["profile_mm"].as_array().unwrap();
        let values = f["profile_mw_per_mm2"].as_array().unwrap();
        let edge = model.envelope_radius_mm();
        for (r, v) in radii.iter().zip(values) {
            let (r, v) = (r.as_f64().unwrap(), v.as_f64().unwrap());
            // Skip the sloped rim [edge-2a, edge], which the grid resolves coarsely.
            if z > 0.0 && r < edge - 2.0 * fiber::CORE_RADIUS_MM - 2.0 * px {
                close(model.irradiance(r), v, 6e-3, &format!("profile z={z} r={r}"));
            }
        }
    }
}

#[test]
fn temporal_coherence_matches() {
    let fx = fixture();
    let t = fiber::transport(&FiberSettings::default());
    let opd: Vec<f64> = fx["coherence"]["opd_um"].as_array().unwrap().iter().map(|v| v.as_f64().unwrap()).collect();
    let gamma: Vec<f64> = fx["coherence"]["gamma"].as_array().unwrap().iter().map(|v| v.as_f64().unwrap()).collect();
    for (d, g) in opd.iter().zip(gamma) {
        let got = fiber_light_lab::physics::wave::coherence_at(&t, *d);
        assert!((got - g).abs() < 1e-10, "opd {d}: {got} vs {g}");
    }
}
