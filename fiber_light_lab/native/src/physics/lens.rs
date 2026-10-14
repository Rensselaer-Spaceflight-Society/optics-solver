//! Spherical singlet lenses, exact (non-paraxial) 3-D ray tracing.
//!
//! Default part: Thorlabs LB1092, N-BK7 bi-convex, Ø12.7 mm (Ø1/2"),
//! f = 15.0 mm at 587.6 nm, R1 = +14.6 mm, R2 = -14.6 mm, center thickness
//! 4.7 mm, edge thickness 1.8 mm, back focal length 13.4 mm, uncoated.
//! N-BK7 dispersion from the Schott Sellmeier formula.

use super::math::Vec3;
use serde::{Deserialize, Serialize};

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize, Default)]
pub enum Glass {
    #[default]
    NBK7,
}

impl Glass {
    /// Refractive index at `wavelength_nm` (Schott N-BK7 Sellmeier).
    pub fn index(self, wavelength_nm: f64) -> f64 {
        match self {
            Glass::NBK7 => {
                let l2 = (wavelength_nm / 1000.0).powi(2);
                let b = [1.039_612_12, 0.231_792_344, 1.010_469_45];
                let c = [0.006_000_698_67, 0.020_017_914_4, 103.560_653];
                (1.0 + (0..3).map(|i| b[i] * l2 / (l2 - c[i])).sum::<f64>()).sqrt()
            }
        }
    }
    pub fn label(self) -> &'static str {
        match self {
            Glass::NBK7 => "N-BK7 (Schott Sellmeier)",
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize, Default)]
pub enum Coating {
    /// Bare glass: exact Fresnel reflection at each surface.
    #[default]
    Uncoated,
    /// Thorlabs -A: Ravg < 0.5 % for 350-700 nm.
    ArA,
    /// Thorlabs -B: Ravg < 0.5 % for 650-1050 nm.
    ArB,
    /// Thorlabs -AB: Ravg < 0.5 % for 400-1100 nm.
    ArAB,
    /// Lossless surfaces (for comparison only).
    Ideal,
}

impl Coating {
    pub fn label(self) -> &'static str {
        match self {
            Coating::Uncoated => "Uncoated",
            Coating::ArA => "-A  AR 350–700 nm",
            Coating::ArB => "-B  AR 650–1050 nm",
            Coating::ArAB => "-AB AR 400–1100 nm",
            Coating::Ideal => "Ideal (no reflection)",
        }
    }
    fn band(self) -> Option<(f64, f64)> {
        match self {
            Coating::ArA => Some((350.0, 700.0)),
            Coating::ArB => Some((650.0, 1050.0)),
            Coating::ArAB => Some((400.0, 1100.0)),
            _ => None,
        }
    }
}

/// Unpolarized Fresnel power reflectance from index n1 into n2.
/// `cos_i` is the cosine of the incidence angle (positive). TIR returns 1.
pub fn fresnel_reflectance(n1: f64, n2: f64, cos_i: f64) -> f64 {
    let cos_i = cos_i.clamp(0.0, 1.0);
    let sin_t2 = (n1 / n2).powi(2) * (1.0 - cos_i * cos_i);
    if sin_t2 >= 1.0 {
        return 1.0;
    }
    let cos_t = (1.0 - sin_t2).sqrt();
    let rs = ((n1 * cos_i - n2 * cos_t) / (n1 * cos_i + n2 * cos_t)).powi(2);
    let rp = ((n1 * cos_t - n2 * cos_i) / (n1 * cos_t + n2 * cos_i)).powi(2);
    0.5 * (rs + rp)
}

/// Vector Snell refraction. `n` is the surface normal facing the incoming
/// ray (n·d < 0). Returns None on total internal reflection.
pub fn refract(d: Vec3, n: Vec3, n1: f64, n2: f64) -> Option<Vec3> {
    let eta = n1 / n2;
    let cos_i = -n.dot(d);
    let k = 1.0 - eta * eta * (1.0 - cos_i * cos_i);
    if k < 0.0 {
        return None;
    }
    Some((d * eta + n * (eta * cos_i - k.sqrt())).normalized())
}

pub fn reflect(d: Vec3, n: Vec3) -> Vec3 {
    d - n * (2.0 * d.dot(n))
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct LensSpec {
    pub part: String,
    /// Front radius (mm); positive when convex toward the fiber.
    pub r1_mm: f64,
    /// Rear radius (mm); negative when convex toward the sample.
    pub r2_mm: f64,
    pub center_thickness_mm: f64,
    pub diameter_mm: f64,
    /// Aperture of the mount / retaining ring (mm). Rays outside are blocked.
    pub mount_aperture_mm: f64,
    pub glass: Glass,
    pub coating: Coating,
    /// Catalog values for display and comparison only.
    pub catalog_efl_mm: f64,
    pub catalog_bfl_mm: f64,
    pub design_wavelength_nm: f64,
}

impl Default for LensSpec {
    fn default() -> Self {
        Self::lb1092()
    }
}

impl LensSpec {
    /// Thorlabs LB1092 (uncoated). Values from the Thorlabs catalog listing.
    pub fn lb1092() -> Self {
        LensSpec {
            part: "Thorlabs LB1092".into(),
            r1_mm: 14.6,
            r2_mm: -14.6,
            center_thickness_mm: 4.7,
            diameter_mm: 12.7,
            mount_aperture_mm: 12.7,
            glass: Glass::NBK7,
            coating: Coating::Uncoated,
            catalog_efl_mm: 15.0,
            catalog_bfl_mm: 13.4,
            design_wavelength_nm: 587.6,
        }
    }

    pub fn validate(&self) -> Result<(), String> {
        let a = 0.5 * self.diameter_mm;
        if !(self.r1_mm.abs() > a && self.r2_mm.abs() > a) {
            return Err("Lens radii must exceed the semi-diameter.".into());
        }
        if self.center_thickness_mm <= 0.0 || self.edge_thickness_mm() < 0.0 {
            return Err("Lens thickness is not physical for these radii.".into());
        }
        if self.mount_aperture_mm <= 0.0 {
            return Err("Mount aperture must be positive.".into());
        }
        Ok(())
    }

    pub fn index(&self, wavelength_nm: f64) -> f64 {
        self.glass.index(wavelength_nm)
    }

    /// Effective focal length from the thick-lens equation (mm).
    pub fn efl(&self, wavelength_nm: f64) -> f64 {
        let n = self.index(wavelength_nm);
        let (r1, r2, t) = (self.r1_mm, self.r2_mm, self.center_thickness_mm);
        1.0 / ((n - 1.0) * (1.0 / r1 - 1.0 / r2 + (n - 1.0) * t / (n * r1 * r2)))
    }

    /// Back focal length: rear vertex to rear focal point (mm).
    pub fn bfl(&self, wavelength_nm: f64) -> f64 {
        let n = self.index(wavelength_nm);
        let f = self.efl(wavelength_nm);
        f * (1.0 - (n - 1.0) * self.center_thickness_mm / (n * self.r1_mm))
    }

    /// Front focal length: front focal point to front vertex (mm).
    pub fn ffl(&self, wavelength_nm: f64) -> f64 {
        let n = self.index(wavelength_nm);
        let f = self.efl(wavelength_nm);
        f * (1.0 + (n - 1.0) * self.center_thickness_mm / (n * self.r2_mm))
    }

    /// Surface sag at radial height r for radius R (signed, along +z).
    pub fn sag(r: f64, radius: f64) -> f64 {
        radius - radius.signum() * (radius * radius - r * r).max(0.0).sqrt()
    }

    pub fn edge_thickness_mm(&self) -> f64 {
        let a = 0.5 * self.diameter_mm;
        self.center_thickness_mm - Self::sag(a, self.r1_mm) + Self::sag(a, self.r2_mm)
    }

    /// Fiber-to-front-vertex distance that collimates `wavelength_nm`.
    pub fn collimating_distance(&self, wavelength_nm: f64) -> f64 {
        self.ffl(wavelength_nm)
    }

    /// Paraxial fiber distance (to the front vertex) that images the fiber
    /// face onto a plane `image_from_center_mm` beyond the lens center.
    pub fn imaging_distance(&self, wavelength_nm: f64, image_from_center_mm: f64) -> Option<f64> {
        let f = self.efl(wavelength_nm);
        let half = 0.5 * self.center_thickness_mm;
        // Rear principal plane H' = rear vertex + (BFL - f).
        let s_i = image_from_center_mm - (half + self.bfl(wavelength_nm) - f);
        if s_i <= f {
            return None;
        }
        let s_o = f * s_i / (s_i - f);
        // Front principal plane H = front vertex + (f - FFL).
        let d = s_o - (f - self.ffl(wavelength_nm));
        (d > 0.0).then_some(d)
    }

    /// Surface power reflectance for the chosen coating.
    pub fn surface_reflectance(&self, wavelength_nm: f64, n1: f64, n2: f64, cos_i: f64) -> f64 {
        let bare = fresnel_reflectance(n1, n2, cos_i);
        match self.coating {
            Coating::Uncoated => bare,
            Coating::Ideal => 0.0,
            c => {
                let (lo, hi) = c.band().unwrap();
                if bare >= 1.0 {
                    return 1.0; // TIR is unaffected by a coating
                }
                if (lo..=hi).contains(&wavelength_nm) {
                    // Scale the angular dependence to Ravg = 0.5 % at normal incidence.
                    let r0 = fresnel_reflectance(n1, n2, 1.0);
                    (0.005 * bare / r0).min(bare)
                } else {
                    bare // approximate: outside the design band treat as bare glass
                }
            }
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LensLoss {
    /// Missed the front surface / mount aperture.
    Aperture,
    /// Left through the ground edge (between the polished surfaces).
    Edge,
    /// Total internal reflection at the rear surface.
    Tir,
}

#[derive(Clone, Copy, Debug)]
pub struct LensTrace {
    pub front_hit: Vec3,
    pub rear_hit: Vec3,
    pub direction: Vec3,
    /// Product of the two surface transmissions (Fresnel or coating).
    pub transmission: f64,
}

/// Intersect a ray with the spherical cap whose vertex is at z = `vertex_z`.
/// Returns the hit point and a normal facing -z (toward the fiber).
fn intersect_cap(o: Vec3, d: Vec3, vertex_z: f64, radius: f64) -> Option<(Vec3, Vec3)> {
    let c = Vec3::new(0.0, 0.0, vertex_z + radius);
    let oc = o - c;
    let b = oc.dot(d);
    let cc = oc.norm2() - radius * radius;
    let disc = b * b - cc;
    if disc < 0.0 {
        return None;
    }
    let sq = disc.sqrt();
    for t in [-b - sq, -b + sq] {
        if t > 1e-9 {
            let p = o + d * t;
            if (p.z - c.z) * radius < 0.0 {
                let n = (p - c) / radius;
                return Some((p, n));
            }
        }
    }
    None
}

impl LensSpec {
    /// Trace one ray (lens-local coordinates, optical axis +z, front vertex at
    /// `front_vertex_z`) through both surfaces.
    pub fn trace(
        &self,
        origin: Vec3,
        dir: Vec3,
        wavelength_nm: f64,
        front_vertex_z: f64,
    ) -> Result<LensTrace, LensLoss> {
        let n_glass = self.index(wavelength_nm);
        let half_ap = 0.5 * self.diameter_mm.min(self.mount_aperture_mm);
        let half_d = 0.5 * self.diameter_mm;
        let (p1, n1) = intersect_cap(origin, dir, front_vertex_z, self.r1_mm).ok_or(LensLoss::Aperture)?;
        if p1.x.hypot(p1.y) > half_ap {
            return Err(LensLoss::Aperture);
        }
        let cos1 = -n1.dot(dir);
        if cos1 <= 0.0 {
            return Err(LensLoss::Aperture);
        }
        let t1 = 1.0 - self.surface_reflectance(wavelength_nm, 1.0, n_glass, cos1);
        let d1 = refract(dir, n1, 1.0, n_glass).ok_or(LensLoss::Tir)?;
        let (p2, n2) =
            intersect_cap(p1, d1, front_vertex_z + self.center_thickness_mm, self.r2_mm).ok_or(LensLoss::Edge)?;
        let r2 = p2.x.hypot(p2.y);
        if r2 > half_d {
            return Err(LensLoss::Edge);
        }
        if r2 > half_ap {
            return Err(LensLoss::Aperture);
        }
        let cos2 = -n2.dot(d1);
        let r_rear = self.surface_reflectance(wavelength_nm, n_glass, 1.0, cos2);
        if r_rear >= 1.0 {
            return Err(LensLoss::Tir);
        }
        let d2 = refract(d1, n2, n_glass, 1.0).ok_or(LensLoss::Tir)?;
        Ok(LensTrace { front_hit: p1, rear_hit: p2, direction: d2, transmission: t1 * (1.0 - r_rear) })
    }

    /// Outline of the lens cross-section (z, r) for drawing, front vertex at 0.
    pub fn profile(&self, n: usize) -> Vec<(f64, f64)> {
        let a = 0.5 * self.diameter_mm;
        let mut pts = Vec::with_capacity(2 * n + 2);
        for i in 0..=n {
            let r = -a + 2.0 * a * i as f64 / n as f64;
            pts.push((Self::sag(r.abs(), self.r1_mm), r));
        }
        for i in 0..=n {
            let r = a - 2.0 * a * i as f64 / n as f64;
            pts.push((self.center_thickness_mm + Self::sag(r.abs(), self.r2_mm), r));
        }
        pts
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn nbk7_index_matches_schott() {
        // Schott catalog: nd = 1.51680 (587.56 nm), nF = 1.52238, nC = 1.51432.
        assert!((Glass::NBK7.index(587.56) - 1.51680).abs() < 2e-5);
        assert!((Glass::NBK7.index(486.13) - 1.52238).abs() < 2e-5);
        assert!((Glass::NBK7.index(656.27) - 1.51432).abs() < 2e-5);
    }

    #[test]
    fn lb1092_paraxial_data_match_catalog() {
        let l = LensSpec::lb1092();
        let f = l.efl(587.6);
        assert!((f - 15.0).abs() / 15.0 < 0.01, "EFL {f} not within the ±1% tolerance");
        assert!((l.bfl(587.6) - 13.4).abs() < 0.15, "BFL {}", l.bfl(587.6));
        assert!((l.edge_thickness_mm() - 1.8).abs() < 0.02, "ET {}", l.edge_thickness_mm());
        // Symmetric lens: FFL = BFL.
        assert!((l.ffl(587.6) - l.bfl(587.6)).abs() < 1e-12);
    }

    #[test]
    fn exact_trace_reproduces_paraxial_focus() {
        // A collimated near-axis ray from the sample side focuses at BFL.
        let l = LensSpec::lb1092();
        let h = 0.05;
        let tr = l.trace(Vec3::new(h, 0.0, -5.0), Vec3::Z, 587.6, 0.0).unwrap();
        let p = tr.rear_hit;
        let d = tr.direction;
        let z_cross = p.z - p.x * d.z / d.x;
        let bfl = z_cross - l.center_thickness_mm;
        assert!((bfl - l.bfl(587.6)).abs() < 2e-3, "{bfl} vs {}", l.bfl(587.6));
        // Normal-incidence uncoated transmission ≈ (1-R0)².
        let n = l.index(587.6);
        let r0 = ((n - 1.0) / (n + 1.0)).powi(2);
        assert!((tr.transmission - (1.0 - r0).powi(2)).abs() < 1e-5);
    }

    #[test]
    fn marginal_rays_show_spherical_aberration() {
        // Rays at larger heights cross the axis closer to the lens.
        let l = LensSpec::lb1092();
        let cross = |h: f64| {
            let tr = l.trace(Vec3::new(h, 0.0, -5.0), Vec3::Z, 587.6, 0.0).unwrap();
            tr.rear_hit.z - tr.rear_hit.x * tr.direction.z / tr.direction.x
        };
        assert!(cross(5.0) < cross(1.0) - 0.5);
    }

    #[test]
    fn rays_outside_the_aperture_are_blocked() {
        let l = LensSpec::lb1092();
        assert_eq!(l.trace(Vec3::new(7.0, 0.0, -1.0), Vec3::Z, 600.0, 0.0).unwrap_err(), LensLoss::Aperture);
    }

    #[test]
    fn imaging_distance_is_consistent() {
        let l = LensSpec::lb1092();
        let d = l.imaging_distance(587.6, 63.5).unwrap();
        assert!(d > l.ffl(587.6) && d < 25.0, "{d}");
        assert!(l.imaging_distance(587.6, 10.0).is_none());
    }
}
