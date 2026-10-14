//! Temporal coherence of the modeled output spectrum and the illustrative
//! local interference explorer (port of `physics.py` `coherence_curve`,
//! `wave_basis`, `interference`). The explorer is educational: it is not the
//! speckle pattern of the real cable. Halogen illumination is incoherent.

use super::fiber::{CORE_INDEX, CORE_RADIUS_MM, Transport};
use super::math::{Rng, TAU, linspace, trapezoid};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

/// |γ(δ)| for optical path difference δ (µm) from the output spectrum.
pub fn coherence_at(t: &Transport, opd_um: f64) -> f64 {
    let power = trapezoid(&t.output_spectrum, &t.wavelength_nm);
    if power <= 0.0 {
        return 0.0;
    }
    let (re, im): (Vec<f64>, Vec<f64>) = t
        .wavelength_nm
        .iter()
        .zip(&t.output_spectrum)
        .map(|(&wl, &s)| {
            let ph = TAU * opd_um * 1000.0 / wl;
            (s * ph.cos(), s * ph.sin())
        })
        .unzip();
    let re = trapezoid(&re, &t.wavelength_nm);
    let im = trapezoid(&im, &t.wavelength_nm);
    (re * re + im * im).sqrt() / power
}

/// 301-point curve over 0..max_opd_um, as the Python chart.
pub fn coherence_curve(t: &Transport, max_opd_um: f64) -> (Vec<f64>, Vec<f64>) {
    let opd = linspace(0.0, max_opd_um, 301);
    let g = opd.iter().map(|&d| coherence_at(t, d)).collect();
    (opd, g)
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize, Default)]
pub enum WaveMode {
    /// Exact incoherent ensemble mean (what the halogen source delivers).
    #[default]
    Halogen,
    /// One coherent realization (illustrative speckle).
    Coherent,
    /// Mean of N independent random-phase realizations.
    Ensemble,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct WaveSettings {
    pub mode: WaveMode,
    pub wavelength_nm: f64,
    pub averages: usize,
    pub phase_shift_um: f64,
    pub seed: u64,
}

impl Default for WaveSettings {
    fn default() -> Self {
        WaveSettings { mode: WaveMode::Halogen, wavelength_nm: 650.0, averages: 16, phase_shift_um: 0.0, seed: 4 }
    }
}

#[derive(Clone, Debug, Default)]
pub struct WavePatch {
    /// Patch axis (µm), n points.
    pub axis_um: Vec<f64>,
    /// I / I_incoherent, row-major n×n.
    pub intensity: Vec<f32>,
    pub n: usize,
    pub contrast: f64,
}

const COMPONENTS: usize = 48;
const GRID: usize = 161;

/// Local angular spread: limited by NA near the tip and by the aperture
/// subtended by the core at the sample in the far field.
pub fn local_spread(distance_mm: f64, output_na: f64) -> f64 {
    output_na.min(CORE_RADIUS_MM / CORE_RADIUS_MM.hypot(distance_mm))
}

/// `spread_override`: local angular spread (sin) when optics sit between the
/// fiber and the sample (e.g. from the traced beam after the LB1092).
pub fn interference(
    w: &WaveSettings,
    cable_length_m: f64,
    distance_mm: f64,
    output_na: f64,
    spread_override: Option<f64>,
) -> WavePatch {
    let spread = spread_override.unwrap_or_else(|| local_spread(distance_mm, output_na)).clamp(1e-4, 0.99);
    let lam = w.wavelength_nm / 1000.0;
    let width = 12.0 * lam / spread;
    let axis = linspace(-width / 2.0, width / 2.0, GRID);
    if w.mode == WaveMode::Halogen {
        return WavePatch { axis_um: axis, intensity: vec![1.0; GRID * GRID], n: GRID, contrast: 0.0 };
    }
    let mut rng = Rng::new(w.seed);
    let k = TAU / lam;
    let mut sx = [0.0; COMPONENTS];
    let mut sy = [0.0; COMPONENTS];
    let mut radii = [0.0; COMPONENTS];
    for r in radii.iter_mut() {
        *r = spread * rng.uniform().sqrt();
    }
    for j in 0..COMPONENTS {
        let a = TAU * rng.uniform();
        sx[j] = radii[j] * a.cos();
        sy[j] = radii[j] * a.sin();
    }
    let base: Vec<f64> = (0..COMPONENTS).map(|_| TAU * rng.uniform()).collect();
    let beta: Vec<f64> = (0..COMPONENTS)
        .map(|_| {
            let s = output_na / CORE_INDEX * rng.uniform().sqrt();
            k * CORE_INDEX * ((1.0 - s * s).sqrt() - 1.0)
        })
        .collect();
    let kz: Vec<f64> = radii.iter().map(|r| k * ((1.0 - r * r).sqrt() - 1.0)).collect();
    let count = if w.mode == WaveMode::Ensemble { w.averages.max(1) } else { 1 };
    let mut rng2 = Rng::new(w.seed.wrapping_add(77));
    let realizations: Vec<Vec<f64>> = (0..count)
        .map(|r| {
            (0..COMPONENTS)
                .map(|j| {
                    let offset = if count == 1 && r == 0 { 0.0 } else { TAU * rng2.uniform() };
                    base[j] + (beta[j] * cable_length_m * 1e6).rem_euclid(TAU) + kz[j] * w.phase_shift_um + offset
                })
                .collect()
        })
        .collect();
    // Separable phasors exp(i k sx x) and exp(i k sy y).
    let colx: Vec<[(f32, f32); COMPONENTS]> = axis
        .iter()
        .map(|&x| {
            let mut row = [(0.0f32, 0.0f32); COMPONENTS];
            for j in 0..COMPONENTS {
                let p = k * sx[j] * x;
                row[j] = (p.cos() as f32, p.sin() as f32);
            }
            row
        })
        .collect();
    let rowy = colx
        .iter()
        .enumerate()
        .map(|(i, _)| {
            let y = axis[i];
            let mut row = [(0.0f32, 0.0f32); COMPONENTS];
            for j in 0..COMPONENTS {
                let p = k * sy[j] * y;
                row[j] = (p.cos() as f32, p.sin() as f32);
            }
            row
        })
        .collect::<Vec<_>>();
    let amp = 1.0 / (COMPONENTS as f32).sqrt();
    let phasors: Vec<Vec<(f32, f32)>> = realizations
        .iter()
        .map(|ph| ph.iter().map(|&p| ((p.cos() as f32) * amp, (p.sin() as f32) * amp)).collect())
        .collect();
    let intensity: Vec<f32> = (0..GRID * GRID)
        .into_par_iter()
        .map(|idx| {
            let (iy, ix) = (idx / GRID, idx % GRID);
            let mut acc = 0.0f32;
            for a in &phasors {
                let (mut re, mut im) = (0.0f32, 0.0f32);
                for j in 0..COMPONENTS {
                    let (xr, xi) = colx[ix][j];
                    let (yr, yi) = rowy[iy][j];
                    let (pr, pi) = (xr * yr - xi * yi, xr * yi + xi * yr);
                    let (ar, ai) = a[j];
                    re += ar * pr - ai * pi;
                    im += ar * pi + ai * pr;
                }
                acc += re * re + im * im;
            }
            acc / phasors.len() as f32
        })
        .collect();
    let mean = intensity.iter().map(|&v| v as f64).sum::<f64>() / intensity.len() as f64;
    let var = intensity.iter().map(|&v| (v as f64 - mean).powi(2)).sum::<f64>() / intensity.len() as f64;
    WavePatch { axis_um: axis, intensity, n: GRID, contrast: if mean > 0.0 { var.sqrt() / mean } else { 0.0 } }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::physics::fiber::{FiberSettings, transport};

    #[test]
    fn coherence_limits() {
        let t = transport(&FiberSettings::default());
        let (_, g) = coherence_curve(&t, 15.0);
        assert!((g[0] - 1.0).abs() < 1e-12);
        assert!(g.iter().all(|&v| v <= 1.0 + 1e-12));
        assert!(g[300] < 0.1);
    }

    #[test]
    fn interference_statistics() {
        let mut w = WaveSettings { mode: WaveMode::Coherent, ..Default::default() };
        let p = interference(&w, 0.254, 63.5, 0.5, None);
        assert!(p.intensity.iter().all(|&v| v >= 0.0));
        assert!(p.contrast > 0.7, "{}", p.contrast);
        w.mode = WaveMode::Ensemble;
        w.averages = 64;
        let avg = interference(&w, 0.254, 63.5, 0.5, None);
        assert!(avg.contrast < 0.25, "{}", avg.contrast);
        w.mode = WaveMode::Halogen;
        let h = interference(&w, 0.254, 63.5, 0.5, None);
        assert!(h.intensity.iter().all(|&v| v == 1.0));
        // Resolution: > 4 samples per shortest fringe (2·spread/λ).
        let spread = local_spread(63.5, 0.5);
        assert!(p.axis_um[1] - p.axis_um[0] < 0.65 / (4.0 * spread));
    }
}
