//! MI-152 halogen -> GH4001 plastic fiber spectral transport.
//!
//! This is a line-for-line port of `fiber_light_lab/physics.py` (`transport`,
//! `fresnel_transmission`, `spatial_coupling`, `footprint`), checked against
//! values written by `verification/make_native_reference.py`.
//! Units: power mW, spectra mW/nm, lengths mm (geometry) or m (cable).

use super::data;
use super::math::{PI, Curve, gauss_legendre, integrate_clustered, ln_factorial, trapezoid};
use serde::{Deserialize, Serialize};

/// GH-4001-P nominal core radius (980 µm diameter), Mitsubishi data sheet p.1.
pub const CORE_RADIUS_MM: f64 = 0.490;
/// PMMA core index (data sheet p.1), held constant versus wavelength.
pub const CORE_INDEX: f64 = 1.49;
/// Fiber numerical aperture (data sheet p.1).
pub const FIBER_NA: f64 = 0.50;
pub const CLADDING_DIAMETER_MM: f64 = 1.000;
pub const JACKET_DIAMETER_MM: f64 = 2.20;
/// Data sheet maximum loss at 650 nm, 25 °C 50 %RH (170 dB/km).
pub const SPEC_LOSS_650_DB_PER_M: f64 = 0.17;

/// Manufacturer claim: 350,000 foot-candles delivered to the light-guide port.
pub const MI152_RATED_FOOTCANDLES: f64 = 350_000.0;
pub const LUX_PER_FOOTCANDLE: f64 = 10.763_910_416_709_722;
pub const LUMENS_PER_WATT_555: f64 = 683.002;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize, Default)]
pub enum Band {
    #[default]
    #[serde(rename = "400–1000 nm")]
    Nm400To1000,
    #[serde(rename = "400–700 nm")]
    Nm400To700,
}

impl Band {
    pub fn upper_nm(self) -> f64 {
        match self {
            Band::Nm400To1000 => 1000.0,
            Band::Nm400To700 => 700.0,
        }
    }
    pub fn label(self) -> &'static str {
        match self {
            Band::Nm400To1000 => "400–1000 nm",
            Band::Nm400To700 => "400–700 nm",
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize, Default)]
pub enum PowerMode {
    /// Optical power in the band at the fiber face, entered by the user (mW).
    #[default]
    Assumed,
    /// Estimate from the MI-152 rating (350,000 fc at the port), uniform over
    /// the core face. An upper-bound style estimate, not a calibration.
    Mi152Rating,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct FiberSettings {
    pub length_m: f64,
    pub incident_mw: f64,
    pub launch_na: f64,
    pub spot_radius_mm: f64,
    pub offset_mm: f64,
    pub extra_loss_db: f64,
    pub band: Band,
    pub power_mode: PowerMode,
    /// Fraction of rated output when using the MI-152 rating (0-1).
    pub rating_fraction: f64,
}

impl Default for FiberSettings {
    fn default() -> Self {
        // Defaults of the existing Python app: 10 in cable, assumed 1 mW.
        FiberSettings {
            length_m: 0.254,
            incident_mw: 1.0,
            launch_na: 0.50,
            spot_radius_mm: 0.40,
            offset_mm: 0.0,
            extra_loss_db: 0.0,
            band: Band::Nm400To1000,
            power_mode: PowerMode::Assumed,
            rating_fraction: 1.0,
        }
    }
}

impl FiberSettings {
    pub fn validate(&self) -> Result<(), String> {
        let vals = [
            self.length_m,
            self.incident_mw,
            self.launch_na,
            self.spot_radius_mm,
            self.offset_mm,
            self.extra_loss_db,
            self.rating_fraction,
        ];
        if vals.iter().any(|v| !v.is_finite()) {
            return Err("Fiber settings must be finite.".into());
        }
        if self.length_m < 0.0 || self.incident_mw < 0.0 {
            return Err("Cable length and power must be nonnegative.".into());
        }
        if !(self.launch_na > 0.0 && self.launch_na < 1.0) {
            return Err("Launch NA must be between 0 and 1.".into());
        }
        if self.spot_radius_mm <= 0.0 || self.extra_loss_db < 0.0 || self.offset_mm < 0.0 {
            return Err("Spot radius must be positive; loss and offset nonnegative.".into());
        }
        if !(0.0..=1.0).contains(&self.rating_fraction) {
            return Err("Rating fraction must be within 0-1.".into());
        }
        Ok(())
    }

    pub fn output_na(&self) -> f64 {
        self.launch_na.min(FIBER_NA)
    }
}

/// Unpolarized power transmission of one lossless air/PMMA interface.
pub fn fresnel_transmission(sin_external: f64) -> f64 {
    let s = sin_external;
    let c1 = (1.0 - s * s).max(0.0).sqrt();
    let c2 = (1.0 - (s / CORE_INDEX).powi(2)).max(0.0).sqrt();
    let rs = ((c1 - CORE_INDEX * c2) / (c1 + CORE_INDEX * c2)).powi(2);
    let rp = ((CORE_INDEX * c1 - c2) / (CORE_INDEX * c1 + c2)).powi(2);
    1.0 - 0.5 * (rs + rp)
}

/// Fraction of a displaced circular Gaussian (1/e² radius `radius_mm`) that
/// falls on the core disk: the noncentral chi-square CDF with 2 degrees of
/// freedom, as `scipy.stats.ncx2.cdf((a/σ)², 2, (d/σ)²)` with σ = w/2.
pub fn spatial_coupling(radius_mm: f64, offset_mm: f64) -> f64 {
    let sigma = radius_mm / 2.0;
    let x = (CORE_RADIUS_MM / sigma).powi(2);
    let lambda = (offset_mm / sigma).powi(2);
    noncentral_chi2_cdf_2dof(x, lambda)
}

/// CDF of the noncentral chi-square distribution with two degrees of freedom,
/// summed as a Poisson mixture of central chi-square(2+2j) distributions.
pub fn noncentral_chi2_cdf_2dof(x: f64, lambda: f64) -> f64 {
    if x <= 0.0 {
        return 0.0;
    }
    let y = 0.5 * x;
    if lambda <= 0.0 {
        return 1.0 - (-y).exp();
    }
    let mu = 0.5 * lambda;
    let jmax = (mu + 15.0 * mu.sqrt() + 60.0).ceil() as u64;
    let ln_mu = mu.ln();
    let ln_y = y.ln();
    // Q(s, y) = e^-y Σ_{i<s} y^i/i!  (regularized upper gamma, integer s).
    let mut q = (-y).exp(); // s = 1
    let mut total = 0.0;
    for j in 0..=jmax {
        let s = j + 1;
        let weight = (-mu + j as f64 * ln_mu - ln_factorial(j)).exp();
        total += weight * (1.0 - q).max(0.0);
        // advance Q(s, y) -> Q(s+1, y)
        q += (-y + s as f64 * ln_y - ln_factorial(s)).exp();
        q = q.min(1.0);
    }
    total.clamp(0.0, 1.0)
}

/// Midpoint average of T and T² over a uniformly filled NA disk, exactly as
/// the Python model: u = (k+0.5)/2048, sin = NA·sqrt(u).
pub fn face_transmission_averages(output_na: f64) -> (f64, f64) {
    let n = 2048;
    let (mut t1, mut t2) = (0.0, 0.0);
    for k in 0..n {
        let u = (k as f64 + 0.5) / n as f64;
        let t = fresnel_transmission(output_na * u.sqrt());
        t1 += t;
        t2 += t * t;
    }
    (t1 / n as f64, t2 / n as f64)
}

#[derive(Clone, Debug, Default, PartialEq, Serialize)]
pub struct Transport {
    pub wavelength_nm: Vec<f64>,
    pub incident_spectrum: Vec<f64>,
    pub guided_spectrum: Vec<f64>,
    pub output_spectrum: Vec<f64>,
    pub loss_db_per_m: Vec<f64>,
    pub incident_mw: f64,
    pub output_mw: f64,
    pub guided_mw: f64,
    pub coupling: f64,
    pub spatial_fraction: f64,
    pub angular_fraction: f64,
    pub face_transmission: f64,
    pub face_pair_transmission: f64,
    pub half_angle_deg: f64,
    pub output_na: f64,
}

/// Optical power on the 0.98 mm core face implied by the MI-152 rating.
pub fn mi152_rated_core_power_mw(band: Band) -> f64 {
    let illuminance_lm_per_mm2 = MI152_RATED_FOOTCANDLES * LUX_PER_FOOTCANDLE * 1e-6;
    let lumens = illuminance_lm_per_mm2 * PI * CORE_RADIUS_MM * CORE_RADIUS_MM;
    lumens / band_luminous_efficacy(band) * 1000.0
}

/// Luminous efficacy (lm/W) of the EKE spectrum counted over the band.
pub fn band_luminous_efficacy(band: Band) -> f64 {
    let wl: Vec<f64> = (400..=band.upper_nm() as usize).map(|w| w as f64).collect();
    let s: Vec<f64> = wl.iter().map(|&w| data::eke_spectrum().eval(w)).collect();
    let vs: Vec<f64> = wl.iter().zip(&s).map(|(&w, &v)| v * data::cie_v_lambda().eval(w)).collect();
    LUMENS_PER_WATT_555 * trapezoid(&vs, &wl) / trapezoid(&s, &wl)
}

pub fn transport(s: &FiberSettings) -> Transport {
    transport_with(s, data::eke_spectrum(), data::gh4001_loss())
}

pub fn transport_with(s: &FiberSettings, source: &Curve, loss: &Curve) -> Transport {
    let upper = s.band.upper_nm() as usize;
    let wl: Vec<f64> = (400..=upper).map(|w| w as f64).collect();
    let relative: Vec<f64> = wl.iter().map(|&w| source.eval(w)).collect();
    let norm = trapezoid(&relative, &wl);
    let (incident_mw, spatial_fraction) = match s.power_mode {
        PowerMode::Assumed => (s.incident_mw, spatial_coupling(s.spot_radius_mm, s.offset_mm)),
        PowerMode::Mi152Rating => (mi152_rated_core_power_mw(s.band) * s.rating_fraction, 1.0),
    };
    let input: Vec<f64> = relative.iter().map(|r| r / norm * incident_mw).collect();
    let alpha: Vec<f64> = wl.iter().map(|&w| loss.eval_log(w)).collect();
    let angular_fraction = (FIBER_NA / s.launch_na).powi(2).min(1.0);
    let output_na = s.output_na();
    let (t1, t2) = face_transmission_averages(output_na);
    let acceptance = angular_fraction * spatial_fraction;
    let guided: Vec<f64> = input.iter().map(|v| v * acceptance * t1).collect();
    let output: Vec<f64> = input
        .iter()
        .zip(&alpha)
        .map(|(v, a)| v * acceptance * t2 * 10f64.powf(-(a * s.length_m + s.extra_loss_db) / 10.0))
        .collect();
    Transport {
        output_mw: trapezoid(&output, &wl),
        guided_mw: trapezoid(&guided, &wl),
        incident_mw,
        coupling: acceptance * t1,
        spatial_fraction,
        angular_fraction,
        face_transmission: t1,
        face_pair_transmission: t2,
        half_angle_deg: output_na.asin().to_degrees(),
        output_na,
        wavelength_nm: wl,
        incident_spectrum: input,
        guided_spectrum: guided,
        output_spectrum: output,
        loss_db_per_m: alpha,
    }
}

impl Transport {
    /// Guided power remaining after `length_m` of cable (mW).
    pub fn power_along_cable(&self, length_m: f64) -> f64 {
        let v: Vec<f64> = self
            .guided_spectrum
            .iter()
            .zip(&self.loss_db_per_m)
            .map(|(g, a)| g * 10f64.powf(-a * length_m / 10.0))
            .collect();
        trapezoid(&v, &self.wavelength_nm)
    }

    /// Output-spectrum-weighted mean of `f(lambda)`.
    pub fn weighted_mean(&self, f: impl Fn(f64) -> f64) -> f64 {
        let v: Vec<f64> = self.wavelength_nm.iter().zip(&self.output_spectrum).map(|(&w, &p)| p * f(w)).collect();
        let p = trapezoid(&self.output_spectrum, &self.wavelength_nm);
        if p > 0.0 { trapezoid(&v, &self.wavelength_nm) / p } else { 0.0 }
    }
}

/// Bare-fiber output irradiance on a plane normal to the fiber axis at
/// distance `z_mm`: a uniform exit disk convolved with the uniform-radiance
/// angular kernel K = z²/[π NA² (z²+ρ²)²] weighted by two Fresnel faces.
/// Radially symmetric, so the 2-D convolution reduces to one integral per
/// radius (exact geometry, no grid). Returns mW/mm².
pub struct BareFootprint {
    pub z_mm: f64,
    pub output_na: f64,
    pub output_mw: f64,
    t2: f64,
    gl: (Vec<f64>, Vec<f64>),
}

impl BareFootprint {
    pub fn new(z_mm: f64, output_na: f64, output_mw: f64) -> Self {
        BareFootprint {
            z_mm,
            output_na,
            output_mw,
            t2: face_transmission_averages(output_na).1,
            gl: gauss_legendre(96),
        }
    }

    pub fn envelope_radius_mm(&self) -> f64 {
        CORE_RADIUS_MM + self.z_mm * self.output_na.asin().tan()
    }

    fn kernel(&self, rho: f64) -> f64 {
        let z = self.z_mm;
        let rmax = z * self.output_na.asin().tan();
        if rho > rmax * (1.0 + 1e-12) {
            return 0.0;
        }
        let s = rho / (rho * rho + z * z).sqrt();
        let t = fresnel_transmission(s);
        z * z / (PI * self.output_na.powi(2) * (z * z + rho * rho).powi(2)) * t * t / self.t2
    }

    pub fn irradiance(&self, r: f64) -> f64 {
        let a = CORE_RADIUS_MM;
        let disk = self.output_mw / (PI * a * a);
        if self.z_mm <= 0.0 {
            return if r < a { disk } else { 0.0 };
        }
        let rmax = self.z_mm * self.output_na.asin().tan();
        let mut total = 0.0;
        // Fully inside: circle of radius rho around the point lies in the disk.
        let full_end = (a - r).min(rmax);
        if full_end > 0.0 {
            total += integrate_clustered(0.0, full_end, &self.gl, |rho| self.kernel(rho) * 2.0 * PI * rho);
        }
        let lo = (r - a).abs().max(if r < a { full_end } else { 0.0 });
        let hi = (r + a).min(rmax);
        if hi > lo {
            total += integrate_clustered(lo, hi, &self.gl, |rho| {
                if rho <= 0.0 || r <= 0.0 {
                    return 0.0;
                }
                let c = ((rho * rho + r * r - a * a) / (2.0 * rho * r)).clamp(-1.0, 1.0);
                self.kernel(rho) * 2.0 * c.acos() * rho
            });
        }
        disk * total
    }

    /// Radial samples (radius mm, irradiance mW/mm²) out to the envelope.
    pub fn profile(&self, n: usize) -> (Vec<f64>, Vec<f64>) {
        let rmax = self.envelope_radius_mm() * 1.02;
        let r: Vec<f64> = (0..n).map(|i| rmax * i as f64 / (n - 1) as f64).collect();
        let e = r.iter().map(|&v| self.irradiance(v)).collect();
        (r, e)
    }

    /// Radius containing `fraction` of the power (from the radial profile).
    pub fn encircled_radius(&self, fraction: f64) -> f64 {
        let (r, e) = self.profile(1200);
        let mut cum = vec![0.0; r.len()];
        for i in 1..r.len() {
            cum[i] = cum[i - 1] + 0.5 * (e[i] * r[i] + e[i - 1] * r[i - 1]) * 2.0 * PI * (r[i] - r[i - 1]);
        }
        let total = *cum.last().unwrap();
        if total <= 0.0 {
            return 0.0;
        }
        let target = fraction * total;
        let i = cum.partition_point(|&c| c < target).clamp(1, r.len() - 1);
        let f = (target - cum[i - 1]) / (cum[i] - cum[i - 1]).max(1e-300);
        r[i - 1] + f * (r[i] - r[i - 1])
    }
}

/// Supported mode capacity (context only): V = 2πa·NA/λ, M ≈ V²/2.
pub fn mode_count(wavelength_nm: f64) -> (f64, f64) {
    let v = 2.0 * PI * (CORE_RADIUS_MM * 1e6) * FIBER_NA / wavelength_nm;
    (v, v * v / 2.0)
}

/// Cladding index implied by the nominal core index and NA.
pub fn cladding_index() -> f64 {
    (CORE_INDEX * CORE_INDEX - FIBER_NA * FIBER_NA).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn analytic_launch_and_fresnel() {
        let expect = 1.0 - (-2.0 * CORE_RADIUS_MM.powi(2) / 0.16f64).exp();
        assert!((spatial_coupling(0.4, 0.0) - expect).abs() < 1e-14);
        let r0 = ((1.49 - 1.0) / (1.49 + 1.0)) * ((1.49 - 1.0) / (1.49 + 1.0));
        assert!((fresnel_transmission(0.0) - (1.0 - r0)).abs() < 1e-14);
        assert!(mode_count(650.0).1 > 2e6);
        assert!((cladding_index() - 1.404).abs() < 1e-3);
    }

    #[test]
    fn rating_estimate_is_tens_of_milliwatts() {
        let p = mi152_rated_core_power_mw(Band::Nm400To1000);
        assert!(p > 3.0 && p < 60.0, "{p}");
        let k = band_luminous_efficacy(Band::Nm400To700);
        assert!(k > 150.0 && k < 400.0, "{k}");
    }

    #[test]
    fn bare_footprint_conserves_power() {
        for &(z, na) in &[(0.5, 0.5), (8.0, 0.5), (63.5, 0.5), (63.5, 0.2)] {
            let f = BareFootprint::new(z, na, 2.0);
            let (r, e) = f.profile(1500);
            let p: f64 = (1..r.len())
                .map(|i| 0.5 * (e[i] * r[i] + e[i - 1] * r[i - 1]) * 2.0 * PI * (r[i] - r[i - 1]))
                .sum();
            assert!((p / 2.0 - 1.0).abs() < 2e-3, "z={z} na={na} p={p}");
        }
    }
}
