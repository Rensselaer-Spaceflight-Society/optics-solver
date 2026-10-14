//! Interactive 3-D bench view (software projection, painter's algorithm).
//! Physical millimeters on every axis. Drag = orbit, right/middle drag = pan,
//! scroll = zoom, double-click = reset.

use super::style::{self, with_alpha};
use eframe::egui::{self, Align2, Color32, FontId, Pos2, Rect, Sense, Shape, Stroke, pos2};
use fiber_light_lab::physics::bench::CouponShape;
use fiber_light_lab::physics::goniometer::DetectorShape;
use fiber_light_lab::physics::lens::LensSpec;
use fiber_light_lab::physics::math::{Frame, TAU, Vec3};
use fiber_light_lab::physics::sample::{SampleKind, Simulant};
use fiber_light_lab::physics::sim::SimResult;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Camera {
    pub yaw_deg: f64,
    pub pitch_deg: f64,
    pub dist: f64,
    pub target: Vec3,
    pub fov_deg: f64,
}

impl Camera {
    pub fn home(scale: f64) -> Self {
        Camera {
            yaw_deg: -72.0,
            pitch_deg: 16.0,
            dist: scale * 3.4,
            target: Vec3::new(0.0, 0.0, scale * 0.36),
            fov_deg: 34.0,
        }
    }
    pub fn side(scale: f64) -> Self {
        Camera { yaw_deg: -90.0, pitch_deg: 0.0, ..Self::home(scale) }
    }
    pub fn top(scale: f64) -> Self {
        Camera { yaw_deg: -90.0, pitch_deg: 89.0, target: Vec3::ZERO, ..Self::home(scale) }
    }
    fn eye(&self) -> Vec3 {
        let (y, p) = (self.yaw_deg.to_radians(), self.pitch_deg.to_radians());
        self.target + Vec3::new(p.cos() * y.cos(), p.cos() * y.sin(), p.sin()) * self.dist
    }
    fn basis(&self) -> (Vec3, Vec3, Vec3, Vec3) {
        let eye = self.eye();
        let f = (self.target - eye).normalized();
        let mut r = f.cross(Vec3::Z);
        if r.norm() < 1e-9 {
            r = Vec3::X;
        }
        let r = r.normalized();
        let u = r.cross(f);
        (eye, f, r, u)
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct ViewOptions {
    pub show_rays: bool,
    pub ray_count: usize,
    pub show_blocked_rays: bool,
    pub show_spot: bool,
    pub show_lobe: bool,
    pub show_tube: bool,
    pub show_labels: bool,
    pub follow_specular: bool,
    pub detector_angle_deg: f64,
}

impl Default for ViewOptions {
    fn default() -> Self {
        ViewOptions {
            show_rays: true,
            ray_count: 260,
            show_blocked_rays: true,
            show_spot: true,
            show_lobe: true,
            show_tube: true,
            show_labels: true,
            follow_specular: true,
            detector_angle_deg: 30.0,
        }
    }
}

struct Tri {
    p: [Vec3; 3],
    color: Color32,
    lit: bool,
}

struct Seg {
    a: Vec3,
    b: Vec3,
    color: Color32,
    width: f32,
}

struct Label {
    p: Vec3,
    text: String,
    color: Color32,
}

#[derive(Default)]
pub struct Scene {
    tris: Vec<Tri>,
    segs: Vec<Seg>,
    labels: Vec<Label>,
    pub scale: f64,
}

impl Scene {
    fn tri(&mut self, a: Vec3, b: Vec3, c: Vec3, color: Color32, lit: bool) {
        self.tris.push(Tri { p: [a, b, c], color, lit });
    }
    fn quad(&mut self, a: Vec3, b: Vec3, c: Vec3, d: Vec3, color: Color32, lit: bool) {
        self.tri(a, b, c, color, lit);
        self.tri(a, c, d, color, lit);
    }
    fn seg(&mut self, a: Vec3, b: Vec3, color: Color32, width: f32) {
        self.segs.push(Seg { a, b, color, width });
    }
    fn polyline(&mut self, pts: &[Vec3], color: Color32, width: f32) {
        for w in pts.windows(2) {
            self.seg(w[0], w[1], color, width);
        }
    }
    fn dashed(&mut self, a: Vec3, b: Vec3, color: Color32, width: f32, dash: f64) {
        let len = (b - a).norm();
        let n = ((len / dash) as usize).max(1);
        for i in (0..n).step_by(2) {
            let t0 = i as f64 / n as f64;
            let t1 = ((i + 1) as f64 / n as f64).min(1.0);
            self.seg(a + (b - a) * t0, a + (b - a) * t1, color, width);
        }
    }
    fn label(&mut self, p: Vec3, text: impl Into<String>, color: Color32) {
        self.labels.push(Label { p, text: text.into(), color });
    }
    /// Closed cylinder along `axis` from `start` with `len` and `radius`.
    #[allow(clippy::too_many_arguments)]
    fn cylinder(&mut self, start: Vec3, axis: Vec3, len: f64, radius: f64, seg: usize, color: Color32, caps: bool) {
        let f = Frame::from_w_and_v(axis, if axis.dot(Vec3::Y).abs() > 0.9 { Vec3::X } else { Vec3::Y });
        let ring = |z: f64, k: usize| {
            let a = TAU * k as f64 / seg as f64;
            start + f.w * z + f.u * (radius * a.cos()) + f.v * (radius * a.sin())
        };
        for k in 0..seg {
            let (a, b) = (ring(0.0, k), ring(0.0, k + 1));
            let (c, d) = (ring(len, k + 1), ring(len, k));
            self.quad(a, b, c, d, color, true);
            if caps {
                self.tri(start, b, a, color, true);
                self.tri(start + f.w * len, d, c, color, true);
            }
        }
    }
}

fn sample_color(r: &SimResult) -> Color32 {
    let s = &r.settings.sample;
    match s.kind {
        SampleKind::Mirror => Color32::from_rgb(200, 206, 214),
        SampleKind::LambertianReference => Color32::from_rgb(244, 244, 240),
        SampleKind::Composite => {
            if s.composite.mass_fraction() <= 0.0 {
                Color32::from_rgba_unmultiplied(226, 196, 120, 150)
            } else {
                match s.composite.simulant {
                    Simulant::Lms1 => Color32::from_rgb(92, 90, 86),
                    _ => Color32::from_rgb(178, 172, 160),
                }
            }
        }
    }
}

/// Build the world-space scene for a result.
pub fn build_scene(r: &SimResult, opt: &ViewOptions) -> Scene {
    let s = &r.settings;
    let bench = &s.bench;
    let det = &s.detector;
    let big = det.arm_radius_mm.max(bench.tip_to_center_mm());
    let mut sc = Scene { scale: big, ..Default::default() };
    let frame = bench.optical_frame();
    let sf = bench.sample_frame();
    let label_c = Color32::from_gray(60);

    // Goniometer arc (detector path) with ticks.
    let rd = det.arm_radius_mm;
    let arc_pt = |deg: f64, r: f64| {
        let t = deg.to_radians();
        Vec3::new(r * t.sin(), 0.0, r * t.cos())
    };
    for k in (-90..90).step_by(2) {
        if k % 4 == 0 {
            sc.seg(arc_pt(k as f64, rd), arc_pt(k as f64 + 2.0, rd), Color32::from_gray(150), 1.0);
        }
    }
    for k in (-90..=90).step_by(10) {
        let len = if k % 30 == 0 { 0.06 } else { 0.03 } * rd;
        sc.seg(arc_pt(k as f64, rd), arc_pt(k as f64, rd + len), Color32::from_gray(130), 1.0);
        if opt.show_labels && k % 30 == 0 {
            sc.label(arc_pt(k as f64, rd * 1.12), format!("{k}°"), Color32::from_gray(110));
        }
    }
    // Blocked zone on the arc.
    let (bc, bh) = (r.scan.blocked_center_deg, r.scan.blocked_halfwidth_deg);
    let n = 24;
    let pts: Vec<Vec3> = (0..=n).map(|i| arc_pt(bc - bh + 2.0 * bh * i as f64 / n as f64, rd)).collect();
    sc.polyline(&pts, with_alpha(style::RED, 0.75), 3.0);

    // Rotation axis and sample normal.
    let ext = 0.5 * bench.coupon.width_mm.max(bench.coupon.length_mm) + 8.0;
    sc.dashed(Vec3::new(0.0, -ext, 0.0), Vec3::new(0.0, ext, 0.0), Color32::from_gray(120), 1.0, 2.0);
    sc.dashed(Vec3::ZERO, sf.w * (0.55 * rd), Color32::from_gray(120), 1.0, 2.5);
    if opt.show_labels {
        sc.label(sf.w * (0.58 * rd), "normal", label_c);
        sc.label(Vec3::new(0.0, ext + 2.0, 0.0), "rotation axis", Color32::from_gray(110));
    }

    // Coupon.
    let cpn = &bench.coupon;
    let col = sample_color(r);
    let side = Color32::from_rgba_unmultiplied(
        (col.r() as f32 * 0.8) as u8,
        (col.g() as f32 * 0.8) as u8,
        (col.b() as f32 * 0.8) as u8,
        col.a(),
    );
    let to_w = |u: f64, v: f64, w: f64| sf.u * u + sf.v * v + sf.w * w;
    match cpn.shape {
        CouponShape::Rectangle => {
            let (l, wd, t) = (0.5 * cpn.length_mm, 0.5 * cpn.width_mm, cpn.thickness_mm);
            let c = |a: f64, b: f64, z: f64| to_w(a, b, z);
            sc.quad(c(-l, -wd, 0.0), c(l, -wd, 0.0), c(l, wd, 0.0), c(-l, wd, 0.0), col, true);
            sc.quad(c(-l, -wd, -t), c(l, -wd, -t), c(l, -wd, 0.0), c(-l, -wd, 0.0), side, true);
            sc.quad(c(l, -wd, -t), c(l, wd, -t), c(l, wd, 0.0), c(l, -wd, 0.0), side, true);
            sc.quad(c(l, wd, -t), c(-l, wd, -t), c(-l, wd, 0.0), c(l, wd, 0.0), side, true);
            sc.quad(c(-l, wd, -t), c(-l, -wd, -t), c(-l, -wd, 0.0), c(-l, wd, 0.0), side, true);
        }
        CouponShape::Disk => {
            let rr = 0.5 * cpn.length_mm;
            let segs = 48;
            for k in 0..segs {
                let a0 = TAU * k as f64 / segs as f64;
                let a1 = TAU * (k + 1) as f64 / segs as f64;
                let p0 = to_w(rr * a0.cos(), rr * a0.sin(), 0.0);
                let p1 = to_w(rr * a1.cos(), rr * a1.sin(), 0.0);
                sc.tri(Vec3::ZERO, p0, p1, col, true);
                let d = sf.w * (-cpn.thickness_mm);
                sc.quad(p0 + d, p1 + d, p1, p0, side, true);
            }
        }
    }
    if opt.show_labels {
        let what = match s.sample.kind {
            SampleKind::Mirror => "Mirror".to_string(),
            SampleKind::LambertianReference => "Diffuse reference".to_string(),
            SampleKind::Composite => {
                let c = &s.sample.composite;
                if c.mass_fraction() <= 0.0 {
                    format!("Neat epoxy, {}", c.finish.label().to_lowercase())
                } else {
                    format!(
                        "{} {:.0} wt%, {}",
                        match c.simulant {
                            Simulant::Lhs1 => "LHS-1",
                            Simulant::Lms1 => "LMS-1",
                            Simulant::Custom => "Custom",
                        },
                        c.regolith_mass_percent,
                        c.finish.label().to_lowercase()
                    )
                }
            }
        };
        sc.label(to_w(0.0, -0.5 * cpn.width_mm - 4.0, -cpn.thickness_mm), what, label_c);
    }

    // Irradiance on the coupon surface.
    if opt.show_spot && r.spot.peak_w_m2 > 0.0 {
        let m = &r.spot;
        let nv = 48usize;
        let cell_src = 2.0 * m.half_extent_mm / m.n as f64;
        let cell = 2.0 * m.half_extent_mm / nv as f64;
        for iv in 0..nv {
            for iu in 0..nv {
                let u0 = -m.half_extent_mm + iu as f64 * cell;
                let v0 = -m.half_extent_mm + iv as f64 * cell;
                // Box-average the source grid into the display cell.
                let mut acc = 0.0;
                let mut cnt = 0.0;
                let su = ((u0 + m.half_extent_mm) / cell_src) as usize;
                let sv = ((v0 + m.half_extent_mm) / cell_src) as usize;
                let span = (cell / cell_src).ceil() as usize;
                for a in sv..(sv + span).min(m.n) {
                    for b in su..(su + span).min(m.n) {
                        acc += m.irradiance_w_m2[a * m.n + b];
                        cnt += 1.0;
                    }
                }
                let val = if cnt > 0.0 { acc / cnt / m.peak_w_m2 } else { 0.0 };
                if val < 0.02 || !cpn_contains(cpn, u0 + 0.5 * cell, v0 + 0.5 * cell) {
                    continue;
                }
                let c = with_alpha(style::viridis(val.sqrt()), 0.55 + 0.4 * val as f32);
                let z = 0.05;
                sc.quad(
                    to_w(u0, v0, z),
                    to_w(u0 + cell, v0, z),
                    to_w(u0 + cell, v0 + cell, z),
                    to_w(u0, v0 + cell, z),
                    c,
                    false,
                );
            }
        }
    }

    // Fiber, lens and tube.
    let tip = bench.fiber_tip() + frame.u * bench.fiber_offset_x_mm + frame.v * bench.fiber_offset_y_mm;
    let back = -frame.w;
    sc.cylinder(tip, back, 45.0, 1.1, 18, Color32::from_rgb(40, 40, 44), true);
    sc.cylinder(tip + frame.w * 0.01, back, 0.02, 0.49, 18, Color32::from_rgb(255, 210, 120), true);
    if opt.show_labels {
        sc.label(tip + back * 30.0 + frame.u * 4.0, "GH4001 fiber (Ø0.98 mm core)", label_c);
    }
    if bench.use_lens {
        let lens = &bench.lens;
        let fv = bench.lens_front_vertex() + frame.u * 0.0;
        add_lens(&mut sc, lens, fv, &frame);
        if opt.show_labels {
            sc.label(
                fv + frame.w * lens.center_thickness_mm + frame.u * (0.5 * lens.diameter_mm + 6.0),
                lens.part.clone(),
                label_c,
            );
        }
        if opt.show_tube {
            let tube_r = bench.source_housing_radius_mm;
            let z0 = -6.0;
            let z1 = bench.fiber_to_lens_mm + lens.center_thickness_mm + 3.0;
            let c = Color32::from_rgba_unmultiplied(90, 90, 100, 110);
            for z in [z0, 0.5 * (z0 + z1), z1] {
                let pts: Vec<Vec3> = (0..=40)
                    .map(|k| {
                        let a = TAU * k as f64 / 40.0;
                        tip - frame.u * bench.fiber_offset_x_mm - frame.v * bench.fiber_offset_y_mm
                            + frame.w * z
                            + frame.u * (tube_r * a.cos())
                            + frame.v * (tube_r * a.sin())
                    })
                    .collect();
                sc.polyline(&pts, c, 1.0);
            }
            for k in 0..8 {
                let a = TAU * k as f64 / 8.0;
                let off = frame.u * (tube_r * a.cos()) + frame.v * (tube_r * a.sin());
                let base = tip - frame.u * bench.fiber_offset_x_mm - frame.v * bench.fiber_offset_y_mm;
                sc.seg(base + frame.w * z0 + off, base + frame.w * z1 + off, c, 1.0);
            }
        }
    }

    // Rays.
    if opt.show_rays {
        let mut drawn = 0;
        for p in r.beam.paths.iter() {
            if drawn >= opt.ray_count {
                break;
            }
            if !p.reached_sample && !opt.show_blocked_rays {
                continue;
            }
            drawn += 1;
            let c = if p.reached_sample {
                with_alpha(style::wavelength_color(p.wavelength_nm as f64), 0.55)
            } else {
                Color32::from_rgba_unmultiplied(120, 120, 120, 70)
            };
            sc.polyline(&p.points, c, 1.0);
        }
    }

    // Reflected-light lobe from the simulated scan (in the plane of incidence).
    let pmax = r.scan.power_mw.iter().cloned().fold(0.0, f64::max);
    if opt.show_lobe && pmax > 0.0 {
        let lobe_r = 0.82 * rd;
        let pts: Vec<Vec3> = r
            .scan
            .angle_deg
            .iter()
            .zip(&r.scan.power_mw)
            .map(|(&a, &p)| arc_pt(a, lobe_r * (p / pmax).max(0.0).sqrt()))
            .collect();
        let fill = with_alpha(style::ORANGE, 0.20);
        for w in pts.windows(2) {
            sc.tri(Vec3::new(0.0, 0.0, 0.2), w[0], w[1], fill, false);
        }
        sc.polyline(&pts, with_alpha(style::ORANGE, 0.95), 1.8);
        if opt.show_labels {
            let ia = r.scan.power_mw.iter().position(|&p| p == pmax).unwrap_or(0);
            let peak = arc_pt(r.scan.angle_deg[ia], lobe_r * 1.02);
            sc.label(peak + Vec3::new(0.0, 0.0, 4.0), "measured lobe (√ scale)", style::ORANGE);
        }
    }
    // Specular direction.
    let spec = bench.specular_direction();
    sc.dashed(Vec3::ZERO, spec * (rd * 0.98), with_alpha(style::ACCENT, 0.8), 1.2, 3.0);

    // Detector at the preview angle.
    let th = if opt.follow_specular { bench.specular_angle_deg() } else { opt.detector_angle_deg };
    let t = th.to_radians();
    let c = arc_pt(th, rd);
    let n_d = -(c / rd);
    let tu = Vec3::new(t.cos(), 0.0, -t.sin());
    let tv = Vec3::Y;
    sc.seg(Vec3::ZERO, c, Color32::from_gray(150), 1.0);
    let hr = det.housing_radius_mm;
    let body = Color32::from_rgb(60, 66, 78);
    sc.cylinder(c - n_d * 0.5, -n_d, 12.0, hr, 20, body, true);
    let h = 0.5 * det.size_mm;
    let act = Color32::from_rgb(40, 90, 200);
    let face = c + n_d * 0.05;
    match det.shape {
        DetectorShape::Square => sc.quad(
            face + tu * -h + tv * -h,
            face + tu * h + tv * -h,
            face + tu * h + tv * h,
            face + tu * -h + tv * h,
            act,
            false,
        ),
        DetectorShape::Circle => {
            for k in 0..24 {
                let a0 = TAU * k as f64 / 24.0;
                let a1 = TAU * (k + 1) as f64 / 24.0;
                sc.tri(
                    face,
                    face + tu * (h * a0.cos()) + tv * (h * a0.sin()),
                    face + tu * (h * a1.cos()) + tv * (h * a1.sin()),
                    act,
                    false,
                );
            }
        }
    }
    if opt.show_labels {
        let i = r.scan.angle_deg.iter().position(|&a| (a - th).abs() < 0.5 * r.settings.scan.step_deg + 1e-9);
        let txt = match i {
            Some(i) => format!("Detector {th:.0}°: {}", style::si(r.scan.power_mw[i] * 1e-3, "W")),
            None => format!("Detector {th:.0}°"),
        };
        sc.label(c * 1.32 + Vec3::new(0.0, 0.0, hr), txt, style::ACCENT);
    }

    // Axis triad.
    let o = Vec3::new(0.6 * rd, -0.5 * rd, 0.0);
    let l = 0.08 * rd;
    sc.seg(o, o + Vec3::X * l, style::RED, 2.0);
    sc.seg(o, o + Vec3::Y * l, style::GREEN, 2.0);
    sc.seg(o, o + Vec3::Z * l, style::ACCENT, 2.0);
    if opt.show_labels {
        sc.label(o + Vec3::X * (l * 1.3), "x", style::RED);
        sc.label(o + Vec3::Y * (l * 1.3), "y", style::GREEN);
        sc.label(o + Vec3::Z * (l * 1.3), "z", style::ACCENT);
    }
    sc
}

fn cpn_contains(c: &fiber_light_lab::physics::bench::Coupon, u: f64, v: f64) -> bool {
    c.contains(u, v)
}

fn add_lens(sc: &mut Scene, lens: &LensSpec, front_vertex: Vec3, frame: &Frame) {
    let a = 0.5 * lens.diameter_mm;
    let rings = 10;
    let segs = 36;
    let glass = Color32::from_rgba_unmultiplied(150, 200, 235, 120);
    let edge = Color32::from_rgba_unmultiplied(120, 170, 210, 150);
    let p =
        |r: f64, ang: f64, z: f64| front_vertex + frame.w * z + frame.u * (r * ang.cos()) + frame.v * (r * ang.sin());
    for k in 0..segs {
        let a0 = TAU * k as f64 / segs as f64;
        let a1 = TAU * (k + 1) as f64 / segs as f64;
        for i in 0..rings {
            let r0 = a * i as f64 / rings as f64;
            let r1 = a * (i + 1) as f64 / rings as f64;
            let z0 = LensSpec::sag(r0, lens.r1_mm);
            let z1 = LensSpec::sag(r1, lens.r1_mm);
            sc.quad(p(r0, a0, z0), p(r1, a0, z1), p(r1, a1, z1), p(r0, a1, z0), glass, true);
            let ct = lens.center_thickness_mm;
            let w0 = ct + LensSpec::sag(r0, lens.r2_mm);
            let w1 = ct + LensSpec::sag(r1, lens.r2_mm);
            sc.quad(p(r0, a0, w0), p(r0, a1, w0), p(r1, a1, w1), p(r1, a0, w1), glass, true);
        }
        let zf = LensSpec::sag(a, lens.r1_mm);
        let zr = lens.center_thickness_mm + LensSpec::sag(a, lens.r2_mm);
        sc.quad(p(a, a0, zf), p(a, a1, zf), p(a, a1, zr), p(a, a0, zr), edge, true);
    }
}

/// Project, sort and paint the scene. Returns the response for interaction.
pub fn show(ui: &mut egui::Ui, scene: &Scene, cam: &mut Camera, home: Camera) -> egui::Response {
    let (rect, response) = ui.allocate_exact_size(ui.available_size(), Sense::click_and_drag());
    let painter = ui.painter_at(rect);
    painter.rect_filled(rect, 0.0, Color32::from_rgb(250, 250, 252));

    // Interaction.
    if response.dragged_by(egui::PointerButton::Primary) {
        let d = response.drag_delta();
        cam.yaw_deg -= d.x as f64 * 0.35;
        cam.pitch_deg = (cam.pitch_deg + d.y as f64 * 0.35).clamp(-89.0, 89.0);
    }
    let (_, _, right, up) = cam.basis();
    if response.dragged_by(egui::PointerButton::Secondary) || response.dragged_by(egui::PointerButton::Middle) {
        let d = response.drag_delta();
        let k = 2.0 * cam.dist * (cam.fov_deg.to_radians() * 0.5).tan() / rect.height().max(1.0) as f64;
        cam.target = cam.target - right * (d.x as f64 * k) + up * (d.y as f64 * k);
    }
    if response.hovered() {
        let (scroll, zoom) = ui.input(|i| (i.smooth_scroll_delta.y, i.zoom_delta()));
        if scroll != 0.0 {
            cam.dist *= (-scroll as f64 * 0.0015).exp();
        }
        if zoom != 1.0 {
            cam.dist /= zoom as f64;
        }
        cam.dist = cam.dist.clamp(scene.scale * 0.05, scene.scale * 20.0);
    }
    if response.double_clicked() {
        *cam = home;
    }

    let (eye, fwd, right, up) = cam.basis();
    let focal = (rect.height() as f64 * 0.5) / (cam.fov_deg.to_radians() * 0.5).tan();
    let center = rect.center();
    let near = cam.dist * 0.02;
    let project = |p: Vec3| -> Option<(Pos2, f64)> {
        let q = p - eye;
        let z = q.dot(fwd);
        if z < near {
            return None;
        }
        let x = q.dot(right) / z * focal;
        let y = q.dot(up) / z * focal;
        Some((pos2(center.x + x as f32, center.y - y as f32), z))
    };

    enum Prim {
        Tri([Pos2; 3], Color32),
        Seg(Pos2, Pos2, Color32, f32),
    }
    let mut prims: Vec<(f64, Prim)> = Vec::with_capacity(scene.tris.len() + scene.segs.len());
    for t in &scene.tris {
        let (Some(a), Some(b), Some(c)) = (project(t.p[0]), project(t.p[1]), project(t.p[2])) else { continue };
        let mut color = t.color;
        if t.lit {
            let n = (t.p[1] - t.p[0]).cross(t.p[2] - t.p[0]).normalized();
            let centroid = (t.p[0] + t.p[1] + t.p[2]) / 3.0;
            let v = (eye - centroid).normalized();
            let shade = (0.5 + 0.5 * n.dot(v).abs()) as f32;
            color = Color32::from_rgba_unmultiplied(
                (color.r() as f32 * shade) as u8,
                (color.g() as f32 * shade) as u8,
                (color.b() as f32 * shade) as u8,
                color.a(),
            );
        }
        prims.push(((a.1 + b.1 + c.1) / 3.0, Prim::Tri([a.0, b.0, c.0], color)));
    }
    for s in &scene.segs {
        let (Some(a), Some(b)) = (project(s.a), project(s.b)) else { continue };
        // Bias lines slightly toward the viewer so they sit on surfaces.
        prims.push(((a.1 + b.1) / 2.0 - 0.002 * cam.dist, Prim::Seg(a.0, b.0, s.color, s.width)));
    }
    prims.sort_by(|a, b| b.0.total_cmp(&a.0));

    let mut mesh = egui::Mesh::default();
    let flush = |mesh: &mut egui::Mesh, painter: &egui::Painter| {
        if !mesh.is_empty() {
            painter.add(Shape::mesh(std::mem::take(mesh)));
        }
    };
    for (_, p) in &prims {
        match p {
            Prim::Tri(pts, c) => {
                let i = mesh.vertices.len() as u32;
                for q in pts {
                    mesh.colored_vertex(*q, *c);
                }
                mesh.add_triangle(i, i + 1, i + 2);
            }
            Prim::Seg(a, b, c, w) => {
                flush(&mut mesh, &painter);
                painter.line_segment([*a, *b], Stroke::new(*w, *c));
            }
        }
    }
    flush(&mut mesh, &painter);
    for l in &scene.labels {
        if let Some((p, _)) = project(l.p)
            && rect.contains(p)
        {
            painter.text(p, Align2::CENTER_CENTER, &l.text, FontId::proportional(12.0), l.color);
        }
    }
    painter.text(
        rect.left_bottom() + egui::vec2(8.0, -8.0),
        Align2::LEFT_BOTTOM,
        "Drag: orbit · Right-drag: pan · Scroll: zoom · Double-click: reset · units: mm",
        FontId::proportional(11.0),
        Color32::from_gray(120),
    );
    response
}

/// Crop a full-window screenshot to `rect` (points → pixels).
pub fn crop(image: &egui::ColorImage, rect: Rect, pixels_per_point: f32) -> egui::ColorImage {
    let x0 = (rect.min.x * pixels_per_point).max(0.0) as usize;
    let y0 = (rect.min.y * pixels_per_point).max(0.0) as usize;
    let x1 = ((rect.max.x * pixels_per_point) as usize).min(image.size[0]);
    let y1 = ((rect.max.y * pixels_per_point) as usize).min(image.size[1]);
    let (w, h) = (x1.saturating_sub(x0), y1.saturating_sub(y0));
    let mut pixels = Vec::with_capacity(w * h);
    for y in y0..y1 {
        pixels.extend_from_slice(&image.pixels[y * image.size[0] + x0..y * image.size[0] + x1]);
    }
    egui::ColorImage::new([w, h], pixels)
}
