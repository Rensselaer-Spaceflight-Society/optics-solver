//! Small numerical toolkit: 3-vectors, interpolation, quadrature and
//! deterministic random / low-discrepancy sequences. No external math crates.

use serde::{Deserialize, Serialize};
use std::ops::{Add, AddAssign, Div, Mul, Neg, Sub};

pub const PI: f64 = std::f64::consts::PI;
pub const TAU: f64 = std::f64::consts::TAU;

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Vec3 {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Vec3 {
    pub const ZERO: Vec3 = Vec3 { x: 0.0, y: 0.0, z: 0.0 };
    pub const X: Vec3 = Vec3 { x: 1.0, y: 0.0, z: 0.0 };
    pub const Y: Vec3 = Vec3 { x: 0.0, y: 1.0, z: 0.0 };
    pub const Z: Vec3 = Vec3 { x: 0.0, y: 0.0, z: 1.0 };

    pub const fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }
    pub fn dot(self, o: Vec3) -> f64 {
        self.x * o.x + self.y * o.y + self.z * o.z
    }
    pub fn cross(self, o: Vec3) -> Vec3 {
        Vec3::new(self.y * o.z - self.z * o.y, self.z * o.x - self.x * o.z, self.x * o.y - self.y * o.x)
    }
    pub fn norm(self) -> f64 {
        self.dot(self).sqrt()
    }
    pub fn norm2(self) -> f64 {
        self.dot(self)
    }
    pub fn normalized(self) -> Vec3 {
        let n = self.norm();
        if n > 0.0 { self / n } else { self }
    }
    /// Rotate about the y axis by `angle` (radians), right-handed.
    pub fn rotate_y(self, angle: f64) -> Vec3 {
        let (s, c) = angle.sin_cos();
        Vec3::new(c * self.x + s * self.z, self.y, -s * self.x + c * self.z)
    }
    /// Rotate about the x axis by `angle` (radians), right-handed.
    pub fn rotate_x(self, angle: f64) -> Vec3 {
        let (s, c) = angle.sin_cos();
        Vec3::new(self.x, c * self.y - s * self.z, s * self.y + c * self.z)
    }
    pub fn to_array(self) -> [f64; 3] {
        [self.x, self.y, self.z]
    }
}

impl Add for Vec3 {
    type Output = Vec3;
    fn add(self, o: Vec3) -> Vec3 {
        Vec3::new(self.x + o.x, self.y + o.y, self.z + o.z)
    }
}
impl AddAssign for Vec3 {
    fn add_assign(&mut self, o: Vec3) {
        *self = *self + o;
    }
}
impl Sub for Vec3 {
    type Output = Vec3;
    fn sub(self, o: Vec3) -> Vec3 {
        Vec3::new(self.x - o.x, self.y - o.y, self.z - o.z)
    }
}
impl Mul<f64> for Vec3 {
    type Output = Vec3;
    fn mul(self, s: f64) -> Vec3 {
        Vec3::new(self.x * s, self.y * s, self.z * s)
    }
}
impl Mul<Vec3> for f64 {
    type Output = Vec3;
    fn mul(self, v: Vec3) -> Vec3 {
        v * self
    }
}
impl Div<f64> for Vec3 {
    type Output = Vec3;
    fn div(self, s: f64) -> Vec3 {
        Vec3::new(self.x / s, self.y / s, self.z / s)
    }
}
impl Neg for Vec3 {
    type Output = Vec3;
    fn neg(self) -> Vec3 {
        Vec3::new(-self.x, -self.y, -self.z)
    }
}

/// Orthonormal basis whose third axis is `w`.
#[derive(Clone, Copy, Debug)]
pub struct Frame {
    pub u: Vec3,
    pub v: Vec3,
    pub w: Vec3,
}

impl Frame {
    /// `w` is the frame normal; `hint_v` fixes the rotation about it.
    pub fn from_w_and_v(w: Vec3, hint_v: Vec3) -> Self {
        let w = w.normalized();
        let u = hint_v.cross(w).normalized();
        let v = w.cross(u);
        Frame { u, v, w }
    }
    pub fn to_local(&self, a: Vec3) -> Vec3 {
        Vec3::new(a.dot(self.u), a.dot(self.v), a.dot(self.w))
    }
    pub fn to_world(&self, a: Vec3) -> Vec3 {
        self.u * a.x + self.v * a.y + self.w * a.z
    }
}

/// A tabulated curve with numpy.interp semantics (clamped at both ends).
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct Curve {
    pub x: Vec<f64>,
    pub y: Vec<f64>,
}

impl Curve {
    pub fn new(x: Vec<f64>, y: Vec<f64>) -> Self {
        assert_eq!(x.len(), y.len());
        assert!(!x.is_empty());
        assert!(x.windows(2).all(|w| w[1] > w[0]), "curve x must increase");
        Curve { x, y }
    }

    /// Parse a two-column CSV with one header line.
    pub fn from_csv(text: &str) -> Result<Self, String> {
        let mut x = Vec::new();
        let mut y = Vec::new();
        for (i, line) in text.lines().enumerate() {
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let mut parts = line.split([',', ';', '\t']).map(str::trim);
            let a = parts.next().unwrap_or("");
            let b = parts.next().unwrap_or("");
            match (a.parse::<f64>(), b.parse::<f64>()) {
                (Ok(a), Ok(b)) if a.is_finite() && b.is_finite() => {
                    x.push(a);
                    y.push(b);
                }
                _ if i == 0 => continue, // header
                _ => return Err(format!("line {}: expected two numbers, got '{line}'", i + 1)),
            }
        }
        if x.len() < 2 {
            return Err("need at least two rows of data".into());
        }
        let mut pairs: Vec<(f64, f64)> = x.into_iter().zip(y).collect();
        pairs.sort_by(|a, b| a.0.total_cmp(&b.0));
        pairs.dedup_by(|a, b| a.0 == b.0);
        let (x, y) = pairs.into_iter().unzip();
        Ok(Curve::new(x, y))
    }

    pub fn eval(&self, t: f64) -> f64 {
        interp(t, &self.x, &self.y)
    }

    /// Log-linear interpolation (used for attenuation in dB/m).
    pub fn eval_log(&self, t: f64) -> f64 {
        let n = self.x.len();
        if t <= self.x[0] {
            return self.y[0];
        }
        if t >= self.x[n - 1] {
            return self.y[n - 1];
        }
        let i = self.x.partition_point(|&v| v <= t) - 1;
        let f = (t - self.x[i]) / (self.x[i + 1] - self.x[i]);
        (self.y[i].ln() * (1.0 - f) + self.y[i + 1].ln() * f).exp()
    }
}

/// numpy.interp: linear interpolation, clamped outside the table.
pub fn interp(t: f64, xs: &[f64], ys: &[f64]) -> f64 {
    let n = xs.len();
    if t <= xs[0] {
        return ys[0];
    }
    if t >= xs[n - 1] {
        return ys[n - 1];
    }
    let i = xs.partition_point(|&v| v <= t) - 1;
    let f = (t - xs[i]) / (xs[i + 1] - xs[i]);
    ys[i] + f * (ys[i + 1] - ys[i])
}

/// numpy.trapezoid for samples `y` at abscissae `x`.
pub fn trapezoid(y: &[f64], x: &[f64]) -> f64 {
    y.windows(2).zip(x.windows(2)).map(|(yy, xx)| 0.5 * (yy[0] + yy[1]) * (xx[1] - xx[0])).sum()
}

/// Gauss-Legendre nodes and weights on [-1, 1].
pub fn gauss_legendre(n: usize) -> (Vec<f64>, Vec<f64>) {
    let mut nodes = vec![0.0; n];
    let mut weights = vec![0.0; n];
    for i in 0..n.div_ceil(2) {
        let mut z = (PI * (i as f64 + 0.75) / (n as f64 + 0.5)).cos();
        let mut dp = 0.0;
        for _ in 0..100 {
            let (mut p0, mut p1) = (1.0, z);
            for k in 2..=n {
                let p2 = ((2 * k - 1) as f64 * z * p1 - (k - 1) as f64 * p0) / k as f64;
                p0 = p1;
                p1 = p2;
            }
            let p = if n == 0 { 1.0 } else { p1 };
            let pm1 = p0;
            dp = n as f64 * (z * p - pm1) / (z * z - 1.0);
            let dz = p / dp;
            z -= dz;
            if dz.abs() < 1e-15 {
                break;
            }
        }
        nodes[i] = -z;
        nodes[n - 1 - i] = z;
        let w = 2.0 / ((1.0 - z * z) * dp * dp);
        weights[i] = w;
        weights[n - 1 - i] = w;
    }
    (nodes, weights)
}

/// Integrate `f` on [a, b] with a cosine-clustered Gauss-Legendre rule.
/// The substitution removes square-root endpoint singularities.
pub fn integrate_clustered(a: f64, b: f64, gl: &(Vec<f64>, Vec<f64>), f: impl Fn(f64) -> f64) -> f64 {
    if b <= a {
        return 0.0;
    }
    let half = 0.5 * (b - a);
    gl.0.iter()
        .zip(&gl.1)
        .map(|(&t, &w)| {
            // t in [-1,1] -> s = sin(pi t / 2) in [-1, 1]
            let s = (0.5 * PI * t).sin();
            let ds = 0.5 * PI * (0.5 * PI * t).cos();
            w * f(a + half * (s + 1.0)) * half * ds
        })
        .sum()
}

/// SplitMix64: tiny, fast, statistically sound for seeding and shuffles.
#[derive(Clone, Debug)]
pub struct Rng(u64);

impl Rng {
    pub fn new(seed: u64) -> Self {
        Rng(seed ^ 0x9E37_79B9_7F4A_7C15)
    }
    pub fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
    /// Uniform in [0, 1).
    pub fn uniform(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 * (1.0 / (1u64 << 53) as f64)
    }
}

/// Roberts' R_d low-discrepancy sequence with a Cranley-Patterson shift.
/// Stratifies all dimensions jointly; deterministic for a given seed.
#[derive(Clone, Debug)]
pub struct RSequence {
    alpha: Vec<f64>,
    shift: Vec<f64>,
}

impl RSequence {
    pub fn new(dims: usize, seed: u64) -> Self {
        // phi_d solves x^(d+1) = x + 1.
        let mut phi = 2.0f64;
        for _ in 0..64 {
            phi = (1.0 + phi).powf(1.0 / (dims as f64 + 1.0));
        }
        let alpha = (1..=dims).map(|k| (1.0 / phi.powi(k as i32)).fract()).collect();
        let mut rng = Rng::new(seed);
        let shift = (0..dims).map(|_| rng.uniform()).collect();
        RSequence { alpha, shift }
    }
    pub fn dims(&self) -> usize {
        self.alpha.len()
    }
    pub fn point(&self, index: u64, out: &mut [f64]) {
        let n = index as f64 + 1.0;
        for (k, o) in out.iter_mut().enumerate() {
            *o = (self.shift[k] + n * self.alpha[k]).fract();
        }
    }
}

/// ln(n!) for moderate n via lgamma-free summation cache-free approach.
pub fn ln_factorial(n: u64) -> f64 {
    if n < 2 {
        return 0.0;
    }
    if n < 64 {
        return (2..=n).map(|k| (k as f64).ln()).sum();
    }
    // Stirling series, accurate to ~1e-15 relative for n >= 64.
    let x = n as f64;
    x * x.ln() - x + 0.5 * (TAU * x).ln() + 1.0 / (12.0 * x) - 1.0 / (360.0 * x.powi(3)) + 1.0 / (1260.0 * x.powi(5))
}

/// Linear spacing like numpy.linspace (inclusive).
pub fn linspace(a: f64, b: f64, n: usize) -> Vec<f64> {
    if n == 1 {
        return vec![a];
    }
    (0..n).map(|i| a + (b - a) * i as f64 / (n - 1) as f64).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn legendre_integrates_polynomials_exactly() {
        let gl = gauss_legendre(12);
        let s: f64 = gl.0.iter().zip(&gl.1).map(|(x, w)| w * x.powi(10)).sum();
        assert!((s - 2.0 / 11.0).abs() < 1e-13);
        let total: f64 = gl.1.iter().sum();
        assert!((total - 2.0).abs() < 1e-13);
    }

    #[test]
    fn clustered_rule_handles_sqrt_singularity() {
        let gl = gauss_legendre(64);
        // ∫0^1 sqrt(x(1-x)) dx = pi/8
        let v = integrate_clustered(0.0, 1.0, &gl, |x| (x * (1.0 - x)).sqrt());
        assert!((v - PI / 8.0).abs() < 1e-7, "{v}");
    }

    #[test]
    fn interp_matches_numpy_semantics() {
        let xs = [1.0, 2.0, 4.0];
        let ys = [10.0, 20.0, 0.0];
        assert_eq!(interp(0.0, &xs, &ys), 10.0);
        assert_eq!(interp(5.0, &xs, &ys), 0.0);
        assert!((interp(3.0, &xs, &ys) - 10.0).abs() < 1e-12);
        assert!((interp(2.0, &xs, &ys) - 20.0).abs() < 1e-12);
    }

    #[test]
    fn r_sequence_is_equidistributed() {
        let seq = RSequence::new(5, 3);
        let mut p = [0.0; 5];
        let mut mean = [0.0; 5];
        let n = 20000;
        for i in 0..n {
            seq.point(i, &mut p);
            for k in 0..5 {
                mean[k] += p[k] / n as f64;
            }
        }
        for m in mean {
            assert!((m - 0.5).abs() < 2e-3);
        }
    }

    #[test]
    fn stirling_matches_sum() {
        let direct: f64 = (2..=100u64).map(|k| (k as f64).ln()).sum();
        assert!((ln_factorial(100) - direct).abs() < 1e-10);
    }
}
