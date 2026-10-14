//! Sample reflectance models (BRDF, sr⁻¹).
//!
//! * Mirror — alignment / resolution check (ideal specular, constant R).
//! * Lambertian reference — diffuse calibration standard.
//! * Regolith–epoxy composite — a physically based *hypothesis* model with
//!   parameters tied to the experiment variables (simulant, regolith mass
//!   fraction, surface finish):
//!
//!   f = (1 − f_exp)·[f_resin_specular(σ) + f_body] + f_exp·f_Hapke + f_back
//!
//!   - Resin/air interface: Beckmann microfacet lobe (RMS slope σ), Fresnel
//!     for the resin index, Smith shadowing (Walter et al. 2007).
//!   - Body: light that enters the resin and is diffusely returned by the
//!     embedded grains, with the Saunderson / Wolff interface correction.
//!   - Exposed grains (area fraction f_exp, set by abrasion; Delesse: area
//!     fraction on a section = volume fraction): Hapke (2012) BRDF with
//!     parameters from Curtis et al. (2024) fits to Apollo soils.
//!   - Neat epoxy only: the internal reflection from the coupon back face.
//!
//! None of these parameters are measured for your coupons. The model shows
//! expected magnitudes and how the instrument would see them; the experiment
//! decides which effect dominates.

use super::math::{PI, TAU, Vec3, gauss_legendre};
use serde::{Deserialize, Serialize};

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize, Default)]
pub enum SampleKind {
    Mirror,
    LambertianReference,
    #[default]
    Composite,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize, Default)]
pub enum Simulant {
    /// Exolith / Space Resource Technologies LHS-1 lunar highlands simulant.
    #[default]
    Lhs1,
    /// Exolith / Space Resource Technologies LMS-1 lunar mare simulant.
    Lms1,
    Custom,
}

impl Simulant {
    pub fn label(self) -> &'static str {
        match self {
            Simulant::Lhs1 => "LHS-1 (lunar highlands)",
            Simulant::Lms1 => "LMS-1 (lunar mare)",
            Simulant::Custom => "Custom",
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize, Default)]
pub enum Finish {
    #[default]
    AsCast,
    Abraded,
}

impl Finish {
    pub fn label(self) -> &'static str {
        match self {
            Finish::AsCast => "As-cast",
            Finish::Abraded => "Controlled abrasion",
        }
    }
}

/// Hapke (2012) parameters for the exposed-grain component.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct HapkeParams {
    /// Single-scattering albedo.
    pub w: f64,
    /// Double Henyey–Greenstein lobe width b and backscatter fraction c.
    pub b: f64,
    pub c: f64,
    /// Shadow-hiding opposition effect amplitude and width.
    pub b_s0: f64,
    pub h_s: f64,
    /// Macroscopic roughness mean slope angle (deg).
    pub theta_bar_deg: f64,
    /// Filling factor φ (porosity correction K).
    pub filling_factor: f64,
}

impl Default for HapkeParams {
    fn default() -> Self {
        Self::apollo16_smooth()
    }
}

/// Hapke (2012) "hockey stick" relation c(b), clamped to [-1, 1].
pub fn hockey_stick_c(b: f64) -> f64 {
    (3.29 * (-17.4 * b * b).exp() - 0.908).clamp(-1.0, 1.0)
}

impl HapkeParams {
    /// Curtis et al. 2024 (MAPS 59, doi:10.1111/maps.14266) Apollo 16 68810
    /// "smooth" fit: highland analog for LHS-1. c from the hockey stick.
    pub fn apollo16_smooth() -> Self {
        HapkeParams { w: 0.388, b: 0.063, c: hockey_stick_c(0.063), b_s0: 1.0, h_s: 0.221, theta_bar_deg: 14.27, filling_factor: 0.40 }
    }
    /// Curtis et al. 2024 Apollo 11 10084 "smooth" fit: mare analog for LMS-1.
    pub fn apollo11_smooth() -> Self {
        HapkeParams { w: 0.281, b: 0.238, c: hockey_stick_c(0.238), b_s0: 1.0, h_s: 0.032, theta_bar_deg: 13.80, filling_factor: 0.60 }
    }
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct CompositeSettings {
    pub simulant: Simulant,
    /// Regolith content, percent by mass (0 = neat epoxy control).
    pub regolith_mass_percent: f64,
    pub finish: Finish,
    /// Cured resin refractive index (typical bisphenol-A epoxy ≈ 1.55).
    pub resin_index: f64,
    pub resin_density_g_cm3: f64,
    pub grain_density_g_cm3: f64,
    /// RMS surface slope of the resin skin (deg = atan(Rdq)).
    pub slope_as_cast_deg: f64,
    pub slope_abraded_deg: f64,
    /// Exposed-grain area fraction after abrasion = factor × volume fraction.
    pub exposure_factor: f64,
    /// Exposed-grain area fraction of the as-cast face (resin skin ⇒ ≈ 0).
    pub as_cast_exposed_fraction: f64,
    /// Embedded-grain body reflectance relative to the dry powder (wet-look
    /// darkening; Lekner & Dorf 1988).
    pub embedded_darkening: f64,
    /// Reflectance of whatever is behind a transparent neat-epoxy coupon.
    pub backing_reflectance: f64,
    /// Include the internal back-face reflection of neat-epoxy coupons.
    pub back_face_reflection: bool,
    pub back_face_slope_deg: f64,
    pub hapke: HapkeParams,
}

impl Default for CompositeSettings {
    fn default() -> Self {
        let mut s = CompositeSettings {
            simulant: Simulant::Lhs1,
            regolith_mass_percent: 30.0,
            finish: Finish::AsCast,
            resin_index: 1.55,
            resin_density_g_cm3: 1.15,
            grain_density_g_cm3: 2.75,
            slope_as_cast_deg: 1.0,
            slope_abraded_deg: 12.0,
            exposure_factor: 1.0,
            as_cast_exposed_fraction: 0.0,
            embedded_darkening: 0.5,
            backing_reflectance: 0.04,
            back_face_reflection: true,
            back_face_slope_deg: 1.0,
            hapke: HapkeParams::apollo16_smooth(),
        };
        s.apply_simulant(Simulant::Lhs1);
        s
    }
}

impl CompositeSettings {
    /// Load the preset grain properties for a simulant (custom keeps values).
    pub fn apply_simulant(&mut self, sim: Simulant) {
        self.simulant = sim;
        match sim {
            Simulant::Lhs1 => {
                // Anorthosite-dominated; grain density ≈ plagioclase. Estimate.
                self.grain_density_g_cm3 = 2.75;
                self.hapke = HapkeParams::apollo16_smooth();
            }
            Simulant::Lms1 => {
                // Basaltic (pyroxene, glass, ilmenite). Estimate.
                self.grain_density_g_cm3 = 2.92;
                self.hapke = HapkeParams::apollo11_smooth();
            }
            Simulant::Custom => {}
        }
    }

    pub fn mass_fraction(&self) -> f64 {
        (self.regolith_mass_percent / 100.0).clamp(0.0, 1.0)
    }

    /// Regolith volume fraction from mass fraction and densities.
    pub fn volume_fraction(&self) -> f64 {
        let m = self.mass_fraction();
        if m <= 0.0 {
            return 0.0;
        }
        let vg = m / self.grain_density_g_cm3;
        vg / (vg + (1.0 - m) / self.resin_density_g_cm3)
    }

    pub fn exposed_fraction(&self) -> f64 {
        if self.mass_fraction() <= 0.0 {
            return 0.0;
        }
        match self.finish {
            Finish::AsCast => self.as_cast_exposed_fraction,
            Finish::Abraded => self.exposure_factor * self.volume_fraction(),
        }
        .clamp(0.0, 1.0)
    }

    pub fn rms_slope_deg(&self) -> f64 {
        match self.finish {
            Finish::AsCast => self.slope_as_cast_deg,
            Finish::Abraded => self.slope_abraded_deg,
        }
    }

    /// Diffusive (bihemispherical) reflectance of the dry powder for the
    /// Hapke albedo: r0 = (1 − γ)/(1 + γ), γ = √(1 − w).
    pub fn dry_powder_reflectance(&self) -> f64 {
        let g = (1.0 - self.hapke.w).max(0.0).sqrt();
        (1.0 - g) / (1.0 + g)
    }

    /// Internal (under-resin) diffuse reflectance of the coupon body.
    pub fn body_reflectance(&self) -> f64 {
        if self.mass_fraction() <= 0.0 {
            self.backing_reflectance
        } else {
            (self.embedded_darkening * self.dry_powder_reflectance()).clamp(0.0, 0.999)
        }
    }

    pub fn validate(&self) -> Result<(), String> {
        if !(1.2..=2.0).contains(&self.resin_index) {
            return Err("Resin index should be within 1.2-2.0.".into());
        }
        if !(0.0..=100.0).contains(&self.regolith_mass_percent) {
            return Err("Regolith content must be 0-100 % by mass.".into());
        }
        if self.resin_density_g_cm3 <= 0.0 || self.grain_density_g_cm3 <= 0.0 {
            return Err("Densities must be positive.".into());
        }
        let h = &self.hapke;
        if !(0.0..1.0).contains(&h.w) || !(0.0..1.0).contains(&h.b) || !(-1.0..=1.0).contains(&h.c) {
            return Err("Hapke w, b in [0,1), c in [-1,1].".into());
        }
        if h.h_s <= 0.0 || !(0.0..0.75).contains(&h.filling_factor) || !(0.0..60.0).contains(&h.theta_bar_deg) {
            return Err("Hapke h > 0, filling factor < 0.75, θ̄ < 60°.".into());
        }
        if !(0.0..=45.0).contains(&self.slope_as_cast_deg) || !(0.0..=45.0).contains(&self.slope_abraded_deg) {
            return Err("RMS slopes must be within 0-45°.".into());
        }
        Ok(())
    }
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct SampleSettings {
    pub kind: SampleKind,
    /// Front-surface mirror reflectance (protected aluminum ≈ 0.90).
    pub mirror_reflectance: f64,
    /// Diffuse reference reflectance (PTFE/Spectralon-type ≈ 0.99).
    pub lambertian_reflectance: f64,
    pub composite: CompositeSettings,
}

impl Default for SampleSettings {
    fn default() -> Self {
        SampleSettings {
            kind: SampleKind::Composite,
            mirror_reflectance: 0.90,
            lambertian_reflectance: 0.99,
            composite: CompositeSettings::default(),
        }
    }
}

// ---------------------------------------------------------------- helpers

/// Unpolarized dielectric Fresnel reflectance, light incident from index 1
/// onto index n (or from n onto 1 when `from_inside`).
pub fn dielectric_fresnel(n: f64, cos_i: f64, from_inside: bool) -> f64 {
    let (n1, n2) = if from_inside { (n, 1.0) } else { (1.0, n) };
    super::lens::fresnel_reflectance(n1, n2, cos_i.abs())
}

/// Cosine-weighted hemispherical external reflectance of a flat interface.
pub fn hemispherical_fresnel(n: f64) -> f64 {
    let gl = gauss_legendre(64);
    // ∫0^1 F(μ) 2μ dμ
    gl.0.iter().zip(&gl.1).map(|(&x, &w)| {
        let mu = 0.5 * (x + 1.0);
        0.5 * w * dielectric_fresnel(n, mu, false) * 2.0 * mu
    }).sum()
}

fn beckmann_d(cos_h: f64, alpha: f64) -> f64 {
    if cos_h <= 0.0 {
        return 0.0;
    }
    let c2 = cos_h * cos_h;
    let tan2 = (1.0 - c2) / c2;
    (-tan2 / (alpha * alpha)).exp() / (PI * alpha * alpha * c2 * c2)
}

fn smith_g1(cos_v: f64, alpha: f64) -> f64 {
    if cos_v <= 0.0 {
        return 0.0;
    }
    let tan = (1.0 - cos_v * cos_v).max(0.0).sqrt() / cos_v;
    if tan <= 0.0 {
        return 1.0;
    }
    let a = 1.0 / (alpha * tan);
    if a >= 1.6 { 1.0 } else { (3.535 * a + 2.181 * a * a) / (1.0 + 2.276 * a + 2.577 * a * a) }
}

/// Beckmann α from the per-axis RMS slope angle (deg): α = √2·tan(slope).
pub fn beckmann_alpha(rms_slope_deg: f64) -> f64 {
    (std::f64::consts::SQRT_2 * rms_slope_deg.to_radians().tan()).max(1e-4)
}

/// Microfacet reflection lobe with a facet reflectance function F(cos_ih).
fn microfacet(wi: Vec3, wo: Vec3, alpha: f64, facet_r: impl Fn(f64) -> f64) -> f64 {
    let (ci, co) = (wi.z, wo.z);
    if ci <= 0.0 || co <= 0.0 {
        return 0.0;
    }
    let h = (wi + wo).normalized();
    let cih = wi.dot(h).max(0.0);
    facet_r(cih) * beckmann_d(h.z, alpha) * smith_g1(ci, alpha) * smith_g1(co, alpha) / (4.0 * ci * co)
}

// ---------------------------------------------------------------- Hapke

fn hapke_h(x: f64, w: f64) -> f64 {
    let g = (1.0 - w).max(0.0).sqrt();
    let r0 = (1.0 - g) / (1.0 + g);
    if x <= 0.0 {
        return 1.0;
    }
    1.0 / (1.0 - w * x * (r0 + 0.5 * (1.0 - 2.0 * r0 * x) * ((1.0 + x) / x).ln()))
}

fn hapke_porosity_k(phi: f64) -> f64 {
    if phi <= 1e-6 {
        return 1.0;
    }
    let a = 1.209 * phi.powf(2.0 / 3.0);
    -(1.0 - a).ln() / a
}

/// Hapke (1984) macroscopic roughness: returns (μ0e, μe, S).
pub fn hapke_roughness(mu0: f64, mu: f64, psi: f64, theta_bar: f64) -> (f64, f64, f64) {
    if theta_bar <= 1e-6 {
        return (mu0, mu, 1.0);
    }
    let i = mu0.clamp(-1.0, 1.0).acos();
    let e = mu.clamp(-1.0, 1.0).acos();
    let tb = theta_bar.tan();
    let cot_tb = 1.0 / tb;
    let xi = 1.0 / (1.0 + PI * tb * tb).sqrt();
    let cot = |x: f64| if x < 1e-9 { f64::INFINITY } else { 1.0 / x.tan() };
    let e1 = |x: f64| {
        let c = cot(x);
        if c.is_infinite() { 0.0 } else { (-(2.0 / PI) * cot_tb * c).exp() }
    };
    let e2 = |x: f64| {
        let c = cot(x);
        if c.is_infinite() { 0.0 } else { (-(1.0 / PI) * cot_tb * cot_tb * c * c).exp() }
    };
    let eta = |y: f64| xi * (y.cos() + y.sin() * tb * e2(y) / (2.0 - e1(y)));
    let psi = psi.clamp(0.0, PI);
    let f = if psi >= PI - 1e-9 { 0.0 } else { (-2.0 * (0.5 * psi).tan()).exp() };
    let s2 = (0.5 * psi).sin().powi(2);
    let (mu0e, mue, s) = if i <= e {
        let den = 2.0 - e1(e) - psi / PI * e1(i);
        let mu0e = xi * (i.cos() + i.sin() * tb * (psi.cos() * e2(e) + s2 * e2(i)) / den);
        let mue = xi * (e.cos() + e.sin() * tb * (e2(e) - s2 * e2(i)) / den);
        let s = mue / eta(e) * mu0 / eta(i) * xi / (1.0 - f + f * xi * mu0 / eta(i));
        (mu0e, mue, s)
    } else {
        let den = 2.0 - e1(i) - psi / PI * e1(e);
        let mu0e = xi * (i.cos() + i.sin() * tb * (e2(i) - s2 * e2(e)) / den);
        let mue = xi * (e.cos() + e.sin() * tb * (psi.cos() * e2(i) + s2 * e2(e)) / den);
        let s = mue / eta(e) * mu0 / eta(i) * xi / (1.0 - f + f * xi * mu / eta(e));
        (mu0e, mue, s)
    };
    (mu0e.max(1e-9), mue.max(1e-9), s.max(0.0))
}

/// Hapke BRDF (sr⁻¹) = r/μ0 with the porosity factor K (Hapke 2008/2012),
/// shadow-hiding opposition effect and macroscopic roughness. Directions
/// point away from the surface (toward source and detector), normal = +z.
pub fn hapke_brdf(p: &HapkeParams, wi: Vec3, wo: Vec3) -> f64 {
    let (mu0, mu) = (wi.z, wo.z);
    if mu0 <= 0.0 || mu <= 0.0 {
        return 0.0;
    }
    let cos_g = wi.dot(wo).clamp(-1.0, 1.0);
    let g = cos_g.acos();
    // Azimuth difference between the projections onto the surface.
    let (pi_x, pi_y) = (wi.x, wi.y);
    let (po_x, po_y) = (wo.x, wo.y);
    let ni = pi_x.hypot(pi_y);
    let no = po_x.hypot(po_y);
    let psi = if ni < 1e-12 || no < 1e-12 { 0.0 } else { ((pi_x * po_x + pi_y * po_y) / (ni * no)).clamp(-1.0, 1.0).acos() };
    let (mu0e, mue, s) = hapke_roughness(mu0, mu, psi, p.theta_bar_deg.to_radians());
    let k = hapke_porosity_k(p.filling_factor);
    let b = p.b;
    let hg = |cg: f64| (1.0 - b * b) / (1.0 - 2.0 * b * cg + b * b).powf(1.5);
    // g = 0 is backscatter: the (1+c)/2 lobe peaks toward the source.
    let phase = 0.5 * (1.0 + p.c) * hg(cos_g) + 0.5 * (1.0 - p.c) * hg(-cos_g);
    let bs = p.b_s0 / (1.0 + (0.5 * g).tan() / p.h_s);
    let m = hapke_h(mu0e / k, p.w) * hapke_h(mue / k, p.w) - 1.0;
    let r = k * p.w / (4.0 * PI) * mu0e / (mu0e + mue) * (phase * (1.0 + bs) + m) * s;
    (r / mu0).max(0.0)
}

// ---------------------------------------------------------------- model

/// Per-component BRDF values (sr⁻¹) for display.
#[derive(Clone, Copy, Debug, Default, Serialize)]
pub struct BrdfParts {
    pub resin_specular: f64,
    pub body: f64,
    pub grains: f64,
    pub back_face: f64,
}

impl BrdfParts {
    pub fn total(&self) -> f64 {
        self.resin_specular + self.body + self.grains + self.back_face
    }
}

/// A ready-to-evaluate sample model (derived quantities precomputed).
#[derive(Clone, Debug)]
pub struct SampleModel {
    pub kind: SampleKind,
    pub mirror_reflectance: f64,
    pub lambertian_reflectance: f64,
    pub n: f64,
    pub alpha: f64,
    pub exposed: f64,
    pub body_r: f64,
    pub rho_internal: f64,
    pub back_alpha: f64,
    pub back_face: bool,
    pub hapke: HapkeParams,
    pub volume_fraction: f64,
}

impl SampleModel {
    pub fn new(s: &SampleSettings) -> Self {
        let c = &s.composite;
        let n = c.resin_index;
        let rho_e = hemispherical_fresnel(n);
        let sigma_front = c.rms_slope_deg().to_radians().tan();
        let sigma_back = c.back_face_slope_deg.to_radians().tan();
        let sigma_eq = ((std::f64::consts::FRAC_1_SQRT_2 * (n - 1.0) * sigma_front).powi(2) + (n * sigma_back).powi(2)).sqrt();
        SampleModel {
            kind: s.kind,
            mirror_reflectance: s.mirror_reflectance,
            lambertian_reflectance: s.lambertian_reflectance,
            n,
            alpha: beckmann_alpha(c.rms_slope_deg()),
            exposed: c.exposed_fraction(),
            body_r: c.body_reflectance(),
            rho_internal: 1.0 - (1.0 - rho_e) / (n * n),
            back_alpha: (std::f64::consts::SQRT_2 * sigma_eq).max(1e-4),
            back_face: c.back_face_reflection && c.mass_fraction() <= 0.0,
            hapke: c.hapke,
            volume_fraction: c.volume_fraction(),
        }
    }

    /// True for the ideal mirror (handled by explicit ray reflection).
    pub fn is_delta(&self) -> bool {
        self.kind == SampleKind::Mirror
    }

    /// Angular width (rad) of the narrowest lobe, for detector sampling.
    pub fn narrowest_lobe(&self) -> f64 {
        match self.kind {
            SampleKind::Mirror => 0.0,
            SampleKind::LambertianReference => 1.0,
            SampleKind::Composite => {
                let mut a = self.alpha;
                if self.back_face {
                    a = a.min(self.back_alpha);
                }
                a
            }
        }
    }

    /// Components for directions in the sample frame (normal = +z).
    pub fn parts(&self, wi: Vec3, wo: Vec3) -> BrdfParts {
        if wi.z <= 0.0 || wo.z <= 0.0 {
            return BrdfParts::default();
        }
        match self.kind {
            SampleKind::Mirror => BrdfParts::default(),
            SampleKind::LambertianReference => BrdfParts { body: self.lambertian_reflectance / PI, ..Default::default() },
            SampleKind::Composite => {
                let n = self.n;
                let resin = 1.0 - self.exposed;
                let spec = microfacet(wi, wo, self.alpha, |c| dielectric_fresnel(n, c, false));
                let ti = 1.0 - dielectric_fresnel(n, wi.z, false);
                let to = 1.0 - dielectric_fresnel(n, wo.z, false);
                let body = ti * to * self.body_r / (PI * n * n * (1.0 - self.rho_internal * self.body_r));
                let grains = if self.exposed > 0.0 { hapke_brdf(&self.hapke, wi, wo) } else { 0.0 };
                let back = if self.back_face {
                    microfacet(wi, wo, self.back_alpha, |c| {
                        let t_ext = 1.0 - dielectric_fresnel(n, c, false);
                        let sin_t = (1.0 - c * c).max(0.0).sqrt() / n;
                        let cos_t = (1.0 - sin_t * sin_t).sqrt();
                        t_ext * t_ext * dielectric_fresnel(n, cos_t, true)
                    })
                } else {
                    0.0
                };
                BrdfParts {
                    resin_specular: resin * spec,
                    body: resin * body,
                    grains: self.exposed * grains,
                    back_face: resin * back,
                }
            }
        }
    }

    pub fn brdf(&self, wi: Vec3, wo: Vec3) -> f64 {
        self.parts(wi, wo).total()
    }

    /// Directional–hemispherical reflectance for incidence `theta_i` (rad).
    /// Lobes are integrated by stratified sampling of the half vector.
    pub fn albedo(&self, theta_i: f64) -> f64 {
        let wi = Vec3::new(-theta_i.sin(), 0.0, theta_i.cos());
        if self.kind == SampleKind::Mirror {
            return self.mirror_reflectance;
        }
        // Diffuse-like parts on a cosine-weighted grid.
        let (nt, np) = (96, 192);
        let mut diffuse = 0.0;
        for a in 0..nt {
            let u = (a as f64 + 0.5) / nt as f64;
            let ct = (1.0 - u).sqrt();
            let st = u.sqrt();
            for b in 0..np {
                let ph = TAU * (b as f64 + 0.5) / np as f64;
                let wo = Vec3::new(st * ph.cos(), st * ph.sin(), ct);
                let p = self.parts(wi, wo);
                // pdf = cos/π; integrand f·cos -> f·π
                diffuse += (p.body + p.grains) * PI / (nt * np) as f64;
            }
        }
        if self.kind != SampleKind::Composite {
            return diffuse;
        }
        let lobe = |alpha: f64, which: fn(&BrdfParts) -> f64| {
            let (nu, nv) = (256, 128);
            let mut acc = 0.0;
            for a in 0..nu {
                // Beckmann half-vector sampling: tan²θh = −α² ln(1−u)
                let u = (a as f64 + 0.5) / nu as f64;
                let tan2 = -alpha * alpha * (1.0 - u).ln();
                let ch = 1.0 / (1.0 + tan2).sqrt();
                let sh = (1.0 - ch * ch).sqrt();
                for b in 0..nv {
                    let ph = TAU * (b as f64 + 0.5) / nv as f64;
                    let h = Vec3::new(sh * ph.cos(), sh * ph.sin(), ch);
                    let wo = h * (2.0 * wi.dot(h)) - wi;
                    if wo.z <= 0.0 || wi.dot(h) <= 0.0 {
                        continue;
                    }
                    let pdf_h = beckmann_d(ch, alpha) * ch;
                    let pdf_o = pdf_h / (4.0 * wi.dot(h));
                    acc += which(&self.parts(wi, wo)) * wo.z / pdf_o / (nu * nv) as f64;
                }
            }
            acc
        };
        let mut total = diffuse + lobe(self.alpha, |p| p.resin_specular);
        if self.back_face {
            total += lobe(self.back_alpha, |p| p.back_face);
        }
        total
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn dir(theta_deg: f64, phi_deg: f64) -> Vec3 {
        let (t, p) = (theta_deg.to_radians(), phi_deg.to_radians());
        Vec3::new(t.sin() * p.cos(), t.sin() * p.sin(), t.cos())
    }

    fn composite(finish: Finish, percent: f64, sim: Simulant) -> SampleModel {
        let mut c = CompositeSettings::default();
        c.apply_simulant(sim);
        c.finish = finish;
        c.regolith_mass_percent = percent;
        SampleModel::new(&SampleSettings { composite: c, ..Default::default() })
    }

    #[test]
    fn volume_fraction_and_exposure() {
        let mut c = CompositeSettings::default();
        c.regolith_mass_percent = 60.0;
        let v = c.volume_fraction();
        let expect = (0.6 / 2.75) / (0.6 / 2.75 + 0.4 / 1.15);
        assert!((v - expect).abs() < 1e-12);
        c.finish = Finish::Abraded;
        assert!((c.exposed_fraction() - v).abs() < 1e-12);
        c.regolith_mass_percent = 0.0;
        assert_eq!(c.exposed_fraction(), 0.0);
    }

    #[test]
    fn reciprocity_of_resin_and_body_terms() {
        for m in [composite(Finish::AsCast, 30.0, Simulant::Lhs1), composite(Finish::Abraded, 0.0, Simulant::Lms1)] {
            for (a, b) in [((30.0, 180.0), (25.0, 0.0)), ((60.0, 180.0), (40.0, 10.0)), ((10.0, 90.0), (70.0, 0.0))] {
                let (wi, wo) = (dir(a.0, a.1), dir(b.0, b.1));
                let (p, q) = (m.parts(wi, wo), m.parts(wo, wi));
                assert!((p.resin_specular - q.resin_specular).abs() <= 1e-9 * p.resin_specular.max(1e-12));
                assert!((p.body - q.body).abs() <= 1e-12);
                assert!((p.back_face - q.back_face).abs() <= 1e-9 * p.back_face.max(1e-12));
            }
        }
    }

    #[test]
    fn smooth_hapke_is_reciprocal() {
        let p = HapkeParams { theta_bar_deg: 0.0, ..HapkeParams::apollo11_smooth() };
        let (wi, wo) = (dir(30.0, 180.0), dir(50.0, 20.0));
        let (a, b) = (hapke_brdf(&p, wi, wo), hapke_brdf(&p, wo, wi));
        assert!((a - b).abs() < 1e-12 * a.max(1.0));
    }

    #[test]
    fn energy_is_conserved() {
        let lambert = SampleModel::new(&SampleSettings { kind: SampleKind::LambertianReference, ..Default::default() });
        assert!((lambert.albedo(0.5) - 0.99).abs() < 2e-3);
        for sim in [Simulant::Lhs1, Simulant::Lms1] {
            for fin in [Finish::AsCast, Finish::Abraded] {
                for pct in [0.0, 30.0, 60.0] {
                    let m = composite(fin, pct, sim);
                    for t in [0.0f64, 30.0, 60.0, 80.0] {
                        let a = m.albedo(t.to_radians());
                        assert!(a > 0.0 && a < 1.0, "{sim:?} {fin:?} {pct}% θ={t}: albedo {a}");
                    }
                }
            }
        }
    }

    #[test]
    fn smooth_resin_lobe_albedo_is_fresnel() {
        // Narrow lobe of a near-flat interface reflects ≈ F(θ).
        let mut c = CompositeSettings { regolith_mass_percent: 30.0, slope_as_cast_deg: 0.5, ..Default::default() };
        c.embedded_darkening = 0.0;
        let m = SampleModel::new(&SampleSettings { composite: c, ..Default::default() });
        for t in [0.0f64, 30.0, 60.0] {
            let f = dielectric_fresnel(1.55, t.to_radians().cos(), false);
            let a = m.albedo(t.to_radians());
            assert!((a / f - 1.0).abs() < 0.03, "θ={t}: {a} vs F={f}");
        }
    }

    #[test]
    fn saunderson_body_matches_closed_form() {
        // Body albedo = (1−F(θi))·R(1−ρi)/(1−ρi R).
        let m = composite(Finish::AsCast, 45.0, Simulant::Lhs1);
        let t = 0.3f64;
        let wi = Vec3::new(-t.sin(), 0.0, t.cos());
        let (nt, np) = (200, 64);
        let mut a = 0.0;
        for i in 0..nt {
            let u = (i as f64 + 0.5) / nt as f64;
            let (ct, st) = ((1.0 - u).sqrt(), u.sqrt());
            for j in 0..np {
                let ph = TAU * (j as f64 + 0.5) / np as f64;
                a += m.parts(wi, Vec3::new(st * ph.cos(), st * ph.sin(), ct)).body * PI / (nt * np) as f64;
            }
        }
        let f = dielectric_fresnel(m.n, t.cos(), false);
        let expect = (1.0 - f) * m.body_r * (1.0 - m.rho_internal) / (1.0 - m.rho_internal * m.body_r);
        assert!((a / expect - 1.0).abs() < 2e-3, "{a} vs {expect}");
    }

    #[test]
    fn hapke_presets_are_backscattering_and_bounded() {
        for p in [HapkeParams::apollo16_smooth(), HapkeParams::apollo11_smooth()] {
            let wi = dir(30.0, 180.0);
            let back = hapke_brdf(&p, wi, dir(25.0, 180.0)); // near the source
            let fwd = hapke_brdf(&p, wi, dir(25.0, 0.0)); // specular side
            assert!(back > fwd, "{p:?}");
            assert!(back < 1.0 / PI * 3.0);
        }
    }

    #[test]
    fn abrasion_spreads_light_away_from_specular() {
        let cast = composite(Finish::AsCast, 45.0, Simulant::Lhs1);
        let abr = composite(Finish::Abraded, 45.0, Simulant::Lhs1);
        let wi = dir(30.0, 180.0);
        let spec = dir(30.0, 0.0);
        let off = dir(0.0, 0.0);
        assert!(cast.brdf(wi, spec) > 20.0 * abr.brdf(wi, spec));
        assert!(abr.brdf(wi, off) > cast.brdf(wi, off));
    }
}
