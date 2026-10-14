//! Result tabs: beam & spot, goniometer scan, sample BRDF, experiment
//! matrix, fiber & spectrum, interference.

use super::style::{self, ACCENT, GRAY, GREEN, ORANGE, PURPLE, RED, note, si, si_scale, with_alpha};
use eframe::egui::{self, Color32, RichText, TextureHandle, TextureOptions};
use egui_plot::{
    Bar, BarChart, HoverPosition, Legend, Line, LineStyle, MarkerShape, Plot, PlotImage, PlotPoint, PlotPoints, Points,
    Span, Text, VLine,
};
use fiber_light_lab::physics::bench::{Coupon, CouponShape};
use fiber_light_lab::physics::fiber::SPEC_LOSS_650_DB_PER_M;
use fiber_light_lab::physics::goniometer::ScanResult;
use fiber_light_lab::physics::sample as crate_parts;
use fiber_light_lab::physics::sample::{Finish, SampleKind, Simulant};
use fiber_light_lab::physics::sim::{EffectSummary, SimResult, StudyRow};
use fiber_light_lab::physics::wave::{self, WaveMode, WavePatch};
use fiber_light_lab::physics::{data, math::Curve};

fn pts(x: &[f64], y: &[f64]) -> PlotPoints<'static> {
    PlotPoints::from(x.iter().zip(y).map(|(&a, &b)| [a, b]).collect::<Vec<_>>())
}

fn log10_floor(v: f64, floor: f64) -> f64 {
    v.max(floor).log10()
}

/// Tick label for a log10 axis value.
fn log_tick(exp10: f64) -> String {
    let v = 10f64.powf(exp10);
    if (0.01..10_000.0).contains(&v) {
        let t = format!("{v:.3}");
        t.trim_end_matches('0').trim_end_matches('.').to_string()
    } else {
        format!("{v:.0e}")
    }
}

/// Strong grid marks at integer category positions (bar charts).
fn category_marks(n: usize) -> impl Fn(egui_plot::GridInput) -> Vec<egui_plot::GridMark> {
    move |_| (0..n).map(|i| egui_plot::GridMark { value: i as f64, step_size: 1.0 }).collect()
}

fn heat_texture(
    ctx: &egui::Context,
    name: &str,
    n_x: usize,
    n_y: usize,
    values: &[f64],
    map: impl Fn(f64) -> Color32,
) -> TextureHandle {
    // values row-major [y][x] with y increasing upward; images store the top row first.
    let mut px = Vec::with_capacity(n_x * n_y);
    for row in (0..n_y).rev() {
        for col in 0..n_x {
            px.push(map(values[row * n_x + col]));
        }
    }
    ctx.load_texture(name, egui::ColorImage::new([n_x, n_y], px), TextureOptions::LINEAR)
}

fn coupon_outline(c: &Coupon) -> Vec<[f64; 2]> {
    match c.shape {
        CouponShape::Rectangle => {
            let (l, w) = (0.5 * c.length_mm, 0.5 * c.width_mm);
            vec![[-l, -w], [l, -w], [l, w], [-l, w], [-l, -w]]
        }
        CouponShape::Disk => (0..=72)
            .map(|k| {
                let a = std::f64::consts::TAU * k as f64 / 72.0;
                [0.5 * c.length_mm * a.cos(), 0.5 * c.length_mm * a.sin()]
            })
            .collect(),
    }
}

// ------------------------------------------------------------ beam & spot

pub struct TextureCache {
    pub spot: Option<(u64, TextureHandle)>,
    pub wave: Option<(String, TextureHandle)>,
}

pub fn beam_tab(ui: &mut egui::Ui, r: &SimResult, generation: u64, tex: &mut TextureCache) {
    let ctx = ui.ctx().clone();
    if tex.spot.as_ref().is_none_or(|t| t.0 != generation) {
        let m = &r.spot;
        let peak = m.peak_w_m2.max(1e-300);
        let t = heat_texture(&ctx, "spot", m.n, m.n, &m.irradiance_w_m2, |v| style::viridis((v / peak).sqrt()));
        tex.spot = Some((generation, t));
    }
    let h = ui.available_height();
    ui.columns(2, |cols| {
        let ui = &mut cols[0];
        style::heading(ui, "Irradiance on the sample plane");
        note(ui, "Sample coordinates: u along the plane of incidence, v along the rotation axis. Color scale: square root of irradiance.");
        let m = &r.spot;
        let e = m.half_extent_mm;
        let tex_id = tex.spot.as_ref().unwrap().1.id();
        let outline = coupon_outline(&r.settings.bench.coupon);
        Plot::new("spot")
            .data_aspect(1.0)
            .height(h * 0.62)
            .x_axis_label("u (mm)")
            .y_axis_label("v (mm)")
            .legend(Legend::default())
            .show(ui, |p| {
                p.image(PlotImage::new("irradiance", tex_id, PlotPoint::new(0.0, 0.0), egui::vec2(2.0 * e as f32, 2.0 * e as f32)));
                p.line(Line::new("coupon edge", PlotPoints::from(outline)).color(Color32::from_rgb(255, 120, 60)).width(1.8));
                p.points(
                    Points::new("beam centroid", PlotPoints::from(vec![[m.centroid_u_mm, m.centroid_v_mm]]))
                        .shape(MarkerShape::Cross)
                        .radius(6.0)
                        .color(RED),
                );
            });
        egui::Grid::new("spotinfo").num_columns(2).show(ui, |ui| {
            ui.label("Peak irradiance");
            ui.label(si(m.peak_w_m2, "W/m²"));
            ui.end_row();
            ui.label("90 % width (u × v)");
            ui.label(format!("{:.2} × {:.2} mm", m.width90_u_mm, m.width90_v_mm));
            ui.end_row();
            ui.label("Centroid");
            ui.label(format!("({:.2}, {:.2}) mm", m.centroid_u_mm, m.centroid_v_mm));
            ui.end_row();
        });
        power_budget(ui, r);

        let ui = &mut cols[1];
        let p = &r.profile;
        style::heading(ui, if r.settings.bench.use_lens { "Beam diameter after the lens" } else { "Beam diameter from the fiber tip" });
        let sample_d = r.settings.bench.source_distance_mm;
        Plot::new("profile")
            .height(h * 0.40)
            .legend(Legend::default())
            .x_axis_label(if r.settings.bench.use_lens { "distance from lens center (mm)" } else { "distance from fiber tip (mm)" })
            .y_axis_label("diameter (mm)")
            .include_y(0.0)
            .show(ui, |pl| {
                pl.line(Line::new("90 % encircled", pts(&p.distance_mm, &p.d90_mm)).color(ACCENT).width(2.0));
                pl.line(Line::new("second-moment (D4σ)", pts(&p.distance_mm, &p.d_second_moment_mm)).color(ORANGE).width(1.5));
                pl.vline(VLine::new("rotation center", sample_d).color(GRAY).style(LineStyle::dashed_loose()));
            });
        style::heading(ui, "Angular distribution after the optics");
        let (k, unit) = si_scale(p.angle_hist_power.iter().cloned().fold(0.0, f64::max) * 1e-3, "W");
        let hist: Vec<f64> = p.angle_hist_power.iter().map(|v| v * 1e-3 * k).collect();
        Plot::new("angles")
            .height(h * 0.30)
            .legend(Legend::default())
            .x_axis_label("angle from optical axis (deg)")
            .y_axis_label(format!("{unit} per 0.1°"))
            .show(ui, |pl| {
                pl.line(Line::new("power", pts(&p.angle_hist_deg, &hist)).color(PURPLE).fill(0.0).width(1.5));
                pl.vline(VLine::new("source-size half-angle a/f", 0.5 * p.divergence_estimate_deg).color(GREEN).style(LineStyle::dashed_loose()));
                pl.vline(VLine::new("90 % half-angle", 0.5 * p.divergence90_full_deg).color(RED).style(LineStyle::dashed_loose()));
            });
        note(
            ui,
            &format!(
                "90 % of the power within ±{:.2}° (full {:.2}°). Source-size limit core/f = {:.2}° full. The difference is lens aberration and dispersion.",
                0.5 * p.divergence90_full_deg,
                p.divergence90_full_deg,
                p.divergence_estimate_deg
            ),
        );
    });
}

fn power_budget(ui: &mut egui::Ui, r: &SimResult) {
    let b = &r.beam.budget;
    style::heading(ui, "Power budget");
    let total = b.fiber_output_mw.max(1e-300);
    let items = [
        ("On coupon", b.on_coupon_mw, GREEN),
        ("Misses coupon", b.off_coupon_mw, ORANGE),
        ("Lens surface reflections", b.lens_reflection_mw, PURPLE),
        ("Lens edge / TIR", b.lens_edge_mw + b.lens_tir_mw, GRAY),
        ("Blocked by mount (outside Ø)", b.blocked_by_mount_mw, RED),
    ];
    let (rect, _) = ui.allocate_exact_size(egui::vec2(ui.available_width(), 18.0), egui::Sense::hover());
    let mut x = rect.left();
    for (_, v, c) in items {
        let w = (v / total) as f32 * rect.width();
        ui.painter().rect_filled(
            egui::Rect::from_min_size(egui::pos2(x, rect.top()), egui::vec2(w, rect.height())),
            0.0,
            c,
        );
        x += w;
    }
    egui::Grid::new("budget").num_columns(3).show(ui, |ui| {
        ui.label(RichText::new("Fiber output").strong());
        ui.label(si(b.fiber_output_mw * 1e-3, "W"));
        ui.label("100 %");
        ui.end_row();
        for (name, v, c) in items {
            ui.label(RichText::new(format!("■ {name}")).color(c));
            ui.label(si(v * 1e-3, "W"));
            ui.label(format!("{:.1} %", 100.0 * v / total));
            ui.end_row();
        }
    });
}

// ------------------------------------------------------------ scan

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ScanUnits {
    Power,
    Current,
    Voltage,
    Meter,
    Brdf,
}

impl ScanUnits {
    fn label(self) -> &'static str {
        match self {
            ScanUnits::Power => "Optical power on detector",
            ScanUnits::Current => "Photocurrent",
            ScanUnits::Voltage => "Amplifier voltage",
            ScanUnits::Meter => "Power-meter reading",
            ScanUnits::Brdf => "Inferred BRDF (1/sr)",
        }
    }
}

pub struct ScanView {
    pub units: ScanUnits,
    pub log: bool,
    pub show_instrument: bool,
    pub show_model: bool,
    pub show_windows: bool,
    pub compare_finish: bool,
    pub measured: Option<(String, Curve)>,
    pub scale_measured: bool,
}

impl Default for ScanView {
    fn default() -> Self {
        ScanView {
            units: ScanUnits::Power,
            log: true,
            show_instrument: true,
            show_model: true,
            show_windows: true,
            compare_finish: true,
            measured: None,
            scale_measured: true,
        }
    }
}

fn scan_values(s: &ScanResult, units: ScanUnits, gain: f64) -> (Vec<f64>, &'static str) {
    match units {
        ScanUnits::Power => (s.power_mw.iter().map(|v| v * 1e-3).collect(), "W"),
        ScanUnits::Current => (s.current_ma.iter().map(|v| v * 1e-3).collect(), "A"),
        ScanUnits::Voltage => (s.current_ma.iter().map(|v| v * 1e-3 * gain).collect(), "V"),
        ScanUnits::Meter => (s.meter_mw.iter().map(|v| v * 1e-3).collect(), "W"),
        ScanUnits::Brdf => (s.brdf_measured.clone(), "1/sr"),
    }
}

pub fn scan_tab(ui: &mut egui::Ui, r: &SimResult, view: &mut ScanView, compare: Option<&(String, ScanResult)>) -> bool {
    let mut want_compare_changed = false;
    ui.horizontal_wrapped(|ui| {
        egui::ComboBox::from_id_salt("units").selected_text(view.units.label()).show_ui(ui, |ui| {
            for u in [ScanUnits::Power, ScanUnits::Current, ScanUnits::Voltage, ScanUnits::Meter, ScanUnits::Brdf] {
                ui.selectable_value(&mut view.units, u, u.label());
            }
        });
        ui.checkbox(&mut view.log, "Log scale");
        ui.checkbox(&mut view.show_windows, "Windows");
        ui.checkbox(&mut view.show_instrument, "Mirror instrument function");
        if view.units == ScanUnits::Brdf {
            ui.checkbox(&mut view.show_model, "Model BRDF (no instrument)");
        }
        if r.settings.sample.kind == SampleKind::Composite {
            want_compare_changed |= ui.checkbox(&mut view.compare_finish, "Compare other finish").changed();
        }
        if ui.button("Load measured scan…").on_hover_text("CSV: detector_angle_deg, signal").clicked()
            && let Some(path) = rfd::FileDialog::new().add_filter("CSV", &["csv", "txt"]).pick_file()
        {
            match std::fs::read_to_string(&path).map_err(|e| e.to_string()).and_then(|t| Curve::from_csv(&t)) {
                Ok(c) => {
                    view.measured =
                        Some((path.file_name().map(|n| n.to_string_lossy().into_owned()).unwrap_or_default(), c))
                }
                Err(e) => {
                    rfd::MessageDialog::new().set_title("Could not read scan").set_description(e).show();
                }
            }
        }
        if view.measured.is_some() {
            ui.checkbox(&mut view.scale_measured, "Scale measured to model peak");
            if ui.button("Clear").clicked() {
                view.measured = None;
            }
        }
    });
    let s = &r.scan;
    let gain = r.settings.detector.transimpedance_v_per_a;
    let (vals, unit) = scan_values(s, view.units, gain);
    let vmax = vals.iter().cloned().fold(0.0, f64::max);
    let (k, unit_lbl) = if view.units == ScanUnits::Brdf { (1.0, unit.to_string()) } else { si_scale(vmax, unit) };
    let floor = (vmax * k * 1e-6).max(1e-30);
    let tr = |v: f64| if view.log { log10_floor(v * k, floor) } else { v * k };
    let m = &s.metrics;
    let sp = m.specular_angle_deg;
    let usable: Vec<[f64; 2]> =
        s.angle_deg.iter().zip(&vals).zip(&s.blocked).filter(|(_, b)| !**b).map(|((&a, &v), _)| [a, tr(v)]).collect();
    let blocked_pts: Vec<[f64; 2]> =
        s.angle_deg.iter().zip(&vals).zip(&s.blocked).filter(|(_, b)| **b).map(|((&a, &v), _)| [a, tr(v)]).collect();
    let unit_axis = unit_lbl.clone();
    let log = view.log;
    let h = ui.available_height() - 150.0;
    Plot::new("scan")
        .height(h.max(200.0))
        .legend(Legend::default())
        .x_axis_label("detector angle from sample normal (deg)  —  source arm at −θi")
        .y_axis_label(format!("{} ({unit_lbl})", view.units.label()))
        .y_axis_formatter(move |mark, _| if log { log_tick(mark.value) } else { format!("{:.3}", mark.value) })
        .label_formatter(move |pos| {
            let (name, p) = match pos {
                HoverPosition::NearDataPoint { plot_name, position, .. } => (*plot_name, *position),
                HoverPosition::Elsewhere { position } => ("", *position),
            };
            let v = if log { 10f64.powf(p.y) } else { p.y };
            Some(format!("{name}\nθ = {:.2}°\n{:.4} {unit_axis}", p.x, v))
        })
        .show(ui, |pl| {
            if view.show_windows {
                pl.span(Span::new("near", (sp - s_near(r))..=(sp + s_near(r))).fill(with_alpha(GREEN, 0.10)));
                let off = r.settings.scan.off_min_deg;
                pl.span(Span::new("away", (sp + off)..=90.0).fill(with_alpha(ACCENT, 0.06)));
                pl.span(Span::new("away", -90.0..=(sp - off)).fill(with_alpha(ACCENT, 0.06)));
            }
            pl.span(
                Span::new(
                    "blocked",
                    (s.blocked_center_deg - s.blocked_halfwidth_deg)..=(s.blocked_center_deg + s.blocked_halfwidth_deg),
                )
                .fill(with_alpha(RED, 0.10)),
            );
            pl.vline(VLine::new("specular angle", sp).color(GRAY).style(LineStyle::dashed_loose()));
            if view.show_instrument && r.settings.sample.kind != SampleKind::Mirror {
                let (iv, _) = scan_values(&r.instrument, view.units, gain);
                let imax = iv.iter().cloned().fold(0.0, f64::max);
                if imax > 0.0 && vmax > 0.0 {
                    let scale = vmax / imax;
                    let y: Vec<[f64; 2]> = r
                        .instrument
                        .angle_deg
                        .iter()
                        .zip(&iv)
                        .filter(|(_, v)| !log || **v * scale * k > floor)
                        .map(|(&a, &v)| [a, tr(v * scale)])
                        .collect();
                    pl.line(
                        Line::new("mirror instrument function (scaled)", PlotPoints::from(y))
                            .color(GRAY)
                            .style(LineStyle::dashed_dense())
                            .width(1.2),
                    );
                }
            }
            if let Some((name, other)) = compare {
                let (ov, _) = scan_values(other, view.units, gain);
                let y: Vec<[f64; 2]> = other
                    .angle_deg
                    .iter()
                    .zip(&ov)
                    .zip(&other.blocked)
                    .filter(|(_, b)| !**b)
                    .map(|((&a, &v), _)| [a, tr(v)])
                    .collect();
                pl.line(Line::new(name.clone(), PlotPoints::from(y)).color(PURPLE).width(1.8));
            }
            if view.units == ScanUnits::Brdf && view.show_model && r.settings.sample.kind != SampleKind::Mirror {
                let y: Vec<[f64; 2]> = s
                    .angle_deg
                    .iter()
                    .zip(&s.brdf_true)
                    .filter(|(_, v)| v.is_finite())
                    .map(|(&a, &v)| [a, tr(v)])
                    .collect();
                pl.line(
                    Line::new("model BRDF (ideal instrument)", PlotPoints::from(y))
                        .color(GREEN)
                        .style(LineStyle::dashed_loose())
                        .width(1.5),
                );
            }
            pl.line(Line::new("simulated scan", PlotPoints::from(usable)).color(ORANGE).width(2.2));
            if !blocked_pts.is_empty() {
                pl.points(
                    Points::new("unreachable (blocked)", PlotPoints::from(blocked_pts))
                        .color(with_alpha(RED, 0.5))
                        .radius(1.5),
                );
            }
            if let Some((name, c)) = &view.measured {
                let mmax = c.y.iter().cloned().fold(0.0, f64::max);
                let scale = if view.scale_measured && mmax > 0.0 { vmax / mmax } else { 1.0 };
                let y: Vec<[f64; 2]> = c.x.iter().zip(&c.y).map(|(&a, &v)| [a, tr(v * scale)]).collect();
                pl.points(
                    Points::new(format!("measured: {name}"), PlotPoints::from(y))
                        .color(Color32::BLACK)
                        .radius(2.5)
                        .shape(MarkerShape::Circle),
                );
            }
        });
    metrics_grid(ui, r, compare);
    want_compare_changed
}

fn s_near(r: &SimResult) -> f64 {
    r.settings.scan.near_halfwidth_deg
}

fn metrics_grid(ui: &mut egui::Ui, r: &SimResult, compare: Option<&(String, ScanResult)>) {
    let m = &r.scan.metrics;
    let peak_current_a = r.scan.current_at(r.scan.metrics.peak_angle_deg) * 1e-3;
    ui.horizontal_wrapped(|ui| {
        let chip = |ui: &mut egui::Ui, k: &str, v: String| {
            egui::Frame::group(ui.style()).show(ui, |ui| {
                ui.vertical(|ui| {
                    ui.label(RichText::new(k).small().color(Color32::from_gray(100)));
                    ui.label(RichText::new(v).strong());
                });
            });
        };
        chip(ui, "Near specular (mean)", si(m.near_mean_mw * 1e-3, "W"));
        chip(ui, "Away from specular (mean)", si(m.off_mean_mw * 1e-3, "W"));
        chip(ui, "Near / away ratio", format!("{:.3}", m.near_to_off_ratio));
        chip(ui, "Peak", format!("{} at {:.1}°", si(m.peak_mw * 1e-3, "W"), m.peak_angle_deg));
        chip(ui, "Peak FWHM", if m.fwhm_deg > 0.0 { format!("{:.2}°", m.fwhm_deg) } else { "—".into() });
        chip(ui, "Mirror FWHM (instrument)", format!("{:.2}°", r.instrument.metrics.fwhm_deg));
        chip(ui, "Peak photocurrent", si(peak_current_a, "A"));
        if let Some((name, other)) = compare {
            chip(ui, &format!("Ratio, {name}"), format!("{:.3}", other.metrics.near_to_off_ratio));
        }
    });
    note(
        ui,
        &format!(
            "Windows: near = |θ − {:.1}°| ≤ {:.1}°, away = |θ − {:.1}°| ≥ {:.1}° ({} and {} positions; blocked angles excluded). Light incident on the coupon: {}.",
            m.specular_angle_deg,
            r.settings.scan.near_halfwidth_deg,
            m.specular_angle_deg,
            r.settings.scan.off_min_deg,
            m.near_points,
            m.off_points,
            si(r.scan.incident_on_coupon_mw * 1e-3, "W")
        ),
    );
}

// ------------------------------------------------------------ sample BRDF

pub fn brdf_tab(ui: &mut egui::Ui, r: &SimResult) {
    let b = &r.brdf;
    let s = &r.settings.sample;
    style::heading(
        ui,
        &format!(
            "Model BRDF in the plane of incidence, θi = {:.1}° (ideal instrument)",
            r.settings.bench.incidence_deg
        ),
    );
    if s.kind == SampleKind::Mirror {
        note(
            ui,
            &format!(
                "An ideal mirror is a delta function: R = {:.2} into the specular direction only. See the scan tab for what the detector records.",
                s.mirror_reflectance
            ),
        );
        return;
    }
    let floor = 1e-5;
    let comp = |f: fn(&fiber_light_lab::physics::sample::BrdfParts) -> f64| -> Vec<[f64; 2]> {
        b.view_deg.iter().zip(&b.parts).map(|(&a, p)| [a, log10_floor(f(p), floor)]).collect()
    };
    let h = ui.available_height() - 160.0;
    Plot::new("brdf")
        .height(h.max(220.0))
        .legend(Legend::default())
        .x_axis_label("view angle from normal (deg; specular side positive)")
        .y_axis_label("BRDF (1/sr, log scale)")
        .y_axis_formatter(|m, _| log_tick(m.value))
        .show(ui, |pl| {
            pl.line(Line::new("total", PlotPoints::from(comp(|p| p.total()))).color(Color32::BLACK).width(2.2));
            if s.kind == SampleKind::Composite {
                let any = |f: fn(&crate_parts::BrdfParts) -> f64| b.parts.iter().any(|p| f(p) > floor);
                if any(|p| p.resin_specular) {
                    pl.line(
                        Line::new("resin surface (specular lobe)", PlotPoints::from(comp(|p| p.resin_specular)))
                            .color(ACCENT)
                            .width(1.4),
                    );
                }
                if any(|p| p.body) {
                    pl.line(Line::new("body under resin", PlotPoints::from(comp(|p| p.body))).color(GREEN).width(1.4));
                }
                if any(|p| p.grains) {
                    pl.line(
                        Line::new("exposed grains (Hapke)", PlotPoints::from(comp(|p| p.grains)))
                            .color(ORANGE)
                            .width(1.4),
                    );
                }
                if any(|p| p.back_face) {
                    pl.line(
                        Line::new("neat-epoxy back face", PlotPoints::from(comp(|p| p.back_face)))
                            .color(PURPLE)
                            .width(1.4),
                    );
                }
            }
            pl.vline(
                VLine::new("specular", r.settings.bench.specular_angle_deg())
                    .color(GRAY)
                    .style(LineStyle::dashed_loose()),
            );
        });
    let c = &s.composite;
    egui::Grid::new("brdfinfo").num_columns(4).spacing([18.0, 4.0]).show(ui, |ui| {
        ui.label("Directional-hemispherical reflectance");
        ui.label(RichText::new(format!("{:.4}", b.albedo)).strong());
        if s.kind == SampleKind::Composite {
            ui.label("Regolith volume fraction");
            ui.label(format!("{:.3}", c.volume_fraction()));
            ui.end_row();
            ui.label("Exposed-grain area fraction");
            ui.label(format!("{:.3}", c.exposed_fraction()));
            ui.label("Body reflectance under resin");
            ui.label(format!("{:.3}", c.body_reflectance()));
            ui.end_row();
            ui.label("Resin RMS slope");
            ui.label(format!("{:.2}°", c.rms_slope_deg()));
            ui.label("Dry-powder reflectance (Hapke w)");
            ui.label(format!("{:.3} (w = {:.3})", c.dry_powder_reflectance(), c.hapke.w));
        }
        ui.end_row();
    });
    note(
        ui,
        "Hypothesis model: parameters are literature analogs and assumptions, not measurements of your coupons. The experiment decides which effect dominates.",
    );
}

// ------------------------------------------------------------ experiment

pub struct StudyUi {
    pub lhs: bool,
    pub lms: bool,
    pub contents: String,
    pub incidences: String,
    pub rows: Vec<StudyRow>,
    pub effects: Vec<EffectSummary>,
    pub error: Option<String>,
}

impl Default for StudyUi {
    fn default() -> Self {
        StudyUi {
            lhs: true,
            lms: true,
            contents: "0, 30, 45, 60".into(),
            incidences: "30, 60".into(),
            rows: vec![],
            effects: vec![],
            error: None,
        }
    }
}

pub fn parse_list(s: &str) -> Result<Vec<f64>, String> {
    s.split([',', ' ', ';'])
        .filter(|t| !t.trim().is_empty())
        .map(|t| t.trim().parse::<f64>().map_err(|_| format!("'{t}' is not a number")))
        .collect()
}

/// Returns true when the user asks to run the study.
pub fn study_tab(ui: &mut egui::Ui, st: &mut StudyUi, running: Option<f64>) -> bool {
    let mut run = false;
    style::heading(ui, "Experiment matrix: does finishing change reflection more than regolith content?");
    note(
        ui,
        "Simulates every coupon type of the plan (as-cast vs abraded twins, 0 % = neat-epoxy control) with the current bench, detector and model parameters.",
    );
    ui.horizontal_wrapped(|ui| {
        ui.checkbox(&mut st.lhs, "LHS-1");
        ui.checkbox(&mut st.lms, "LMS-1");
        ui.label("Regolith wt%:");
        ui.add(egui::TextEdit::singleline(&mut st.contents).desired_width(110.0));
        ui.label("Illumination °:");
        ui.add(egui::TextEdit::singleline(&mut st.incidences).desired_width(70.0));
        match running {
            Some(p) => {
                ui.add(egui::ProgressBar::new(p as f32).desired_width(160.0).show_percentage());
            }
            None => {
                if ui.button(RichText::new("Run experiment matrix").strong()).clicked() {
                    run = true;
                }
            }
        }
    });
    if let Some(e) = &st.error {
        ui.colored_label(RED, e);
    }
    if st.rows.is_empty() {
        return run;
    }
    for e in &st.effects {
        let fin = 10f64.powf(e.finish_effect);
        let con = 10f64.powf(e.content_effect);
        let verdict =
            if e.finish_effect > e.content_effect { "finishing dominates" } else { "regolith content dominates" };
        ui.label(
            RichText::new(format!(
                "{} at {:.0}°: abrasion changes the near/away ratio by ×{:.1} on average; regolith content (30–60 %) by ×{:.1}: {verdict} (under these model assumptions).",
                e.simulant.map(|s| s.label()).unwrap_or(""),
                e.incidence_deg,
                fin,
                con
            ))
            .color(Color32::from_gray(40)),
        );
    }
    let mut incs: Vec<f64> = st.rows.iter().map(|r| r.incidence_deg).collect();
    incs.sort_by(|a, b| a.total_cmp(b));
    incs.dedup();
    let sims: Vec<Simulant> = {
        let mut v: Vec<Simulant> = Vec::new();
        for r in &st.rows {
            if !v.contains(&r.simulant) {
                v.push(r.simulant);
            }
        }
        v
    };
    let h = (ui.available_height() * 0.5).max(180.0);
    let ncols = (incs.len() * sims.len()).max(1);
    ui.columns(ncols, |cols| {
        let mut k = 0;
        for sim in &sims {
            for inc in &incs {
                let ui = &mut cols[k];
                k += 1;
                ui.label(RichText::new(format!("{} · θi = {inc:.0}°", sim.label())).strong());
                let sel: Vec<&StudyRow> =
                    st.rows.iter().filter(|r| r.incidence_deg == *inc && r.simulant == *sim).collect();
                let mut cast = Vec::new();
                let mut abr = Vec::new();
                for (i, r) in sel.iter().filter(|r| r.finish == Finish::AsCast).enumerate() {
                    cast.push(
                        Bar::new(i as f64 - 0.18, r.ratio.max(1e-9).log10())
                            .width(0.34)
                            .name(format!("{:.0} wt%", r.regolith_mass_percent)),
                    );
                }
                for (i, r) in sel.iter().filter(|r| r.finish == Finish::Abraded).enumerate() {
                    abr.push(
                        Bar::new(i as f64 + 0.18, r.ratio.max(1e-9).log10())
                            .width(0.34)
                            .name(format!("{:.0} wt%", r.regolith_mass_percent)),
                    );
                }
                let labels: Vec<String> = sel
                    .iter()
                    .filter(|r| r.finish == Finish::AsCast)
                    .map(|r| format!("{:.0}%", r.regolith_mass_percent))
                    .collect();
                Plot::new(format!("study{k}"))
                    .height(h)
                    .legend(Legend::default())
                    .y_axis_label("log10(near / away)")
                    .x_grid_spacer(category_marks(labels.len()))
                    .x_axis_formatter(move |m, _| {
                        let i = m.value.round();
                        if (m.value - i).abs() < 1e-6 && i >= 0.0 {
                            labels.get(i as usize).cloned().unwrap_or_default()
                        } else {
                            String::new()
                        }
                    })
                    .allow_drag(false)
                    .show(ui, |pl| {
                        pl.bar_chart(BarChart::new("as-cast", cast).color(ACCENT));
                        pl.bar_chart(BarChart::new("abraded", abr).color(ORANGE));
                    });
            }
        }
    });
    egui::ScrollArea::vertical().id_salt("studytable").show(ui, |ui| {
        egui::Grid::new("studyrows").striped(true).num_columns(9).show(ui, |ui| {
            for h in ["Simulant", "wt%", "Finish", "θi", "Near", "Away", "Ratio", "Exposed", "Albedo"] {
                ui.label(RichText::new(h).strong());
            }
            ui.end_row();
            for r in &st.rows {
                ui.label(match r.simulant {
                    Simulant::Lhs1 => "LHS-1",
                    Simulant::Lms1 => "LMS-1",
                    Simulant::Custom => "Custom",
                });
                ui.label(format!("{:.0}", r.regolith_mass_percent));
                ui.label(r.finish.label());
                ui.label(format!("{:.0}°", r.incidence_deg));
                ui.label(si(r.near_mw * 1e-3, "W"));
                ui.label(si(r.off_mw * 1e-3, "W"));
                ui.label(format!("{:.2}", r.ratio));
                ui.label(format!("{:.3}", r.exposed_fraction));
                ui.label(format!("{:.3}", r.albedo));
                ui.end_row();
            }
        });
    });
    run
}

// ------------------------------------------------------------ fiber & spectrum

pub fn fiber_tab(ui: &mut egui::Ui, r: &SimResult, coherence: &(Vec<f64>, Vec<f64>)) {
    let t = &r.transport;
    let h = (ui.available_height() - 60.0) / 2.0;
    ui.columns(2, |cols| {
        let ui = &mut cols[0];
        style::heading(ui, "Spectra (EKE lamp to GH4001 output)");
        let scale = |v: &Vec<f64>| v.iter().map(|x| x * 1000.0).collect::<Vec<_>>();
        Plot::new("spec")
            .height(h)
            .legend(Legend::default())
            .x_axis_label("wavelength (nm)")
            .y_axis_label("µW/nm")
            .show(ui, |pl| {
                pl.line(
                    Line::new("incident", pts(&t.wavelength_nm, &scale(&t.incident_spectrum))).color(GRAY).width(1.5),
                );
                pl.line(
                    Line::new("inside fiber", pts(&t.wavelength_nm, &scale(&t.guided_spectrum)))
                        .color(ACCENT)
                        .width(1.5),
                );
                pl.line(
                    Line::new("fiber output", pts(&t.wavelength_nm, &scale(&t.output_spectrum)))
                        .color(ORANGE)
                        .width(2.0),
                );
            });
        style::heading(ui, "GH4001 attenuation (typical, launch NA 0.1)");
        let loss = data::gh4001_loss();
        let ll: Vec<f64> = loss.y.iter().map(|v| v.log10()).collect();
        Plot::new("loss")
            .height(h)
            .legend(Legend::default())
            .x_axis_label("wavelength (nm)")
            .y_axis_label("dB/m (log)")
            .y_axis_formatter(|m, _| log_tick(m.value))
            .show(ui, |pl| {
                pl.line(Line::new("approximate typical curve", pts(&loss.x, &ll)).color(ACCENT).width(1.5));
                pl.points(
                    Points::new(
                        "650 nm datasheet maximum",
                        PlotPoints::from(vec![[650.0, SPEC_LOSS_650_DB_PER_M.log10()]]),
                    )
                    .color(ORANGE)
                    .radius(4.0),
                );
            });
        let ui = &mut cols[1];
        style::heading(ui, "Guided power along the cable");
        let lmax = (r.settings.fiber.length_m * 1.1).max(0.3);
        let xs: Vec<f64> = (0..100).map(|i| lmax * i as f64 / 99.0).collect();
        let ys: Vec<f64> = xs.iter().map(|&l| t.power_along_cable(l)).collect();
        Plot::new("along")
            .height(h)
            .legend(Legend::default())
            .x_axis_label("distance along cable (m)")
            .y_axis_label("mW")
            .include_y(0.0)
            .show(ui, |pl| {
                pl.line(Line::new("guided power", pts(&xs, &ys)).color(ACCENT).width(1.8));
                pl.vline(
                    VLine::new("cable end", r.settings.fiber.length_m).color(ORANGE).style(LineStyle::dashed_loose()),
                );
            });
        style::heading(ui, "Temporal coherence of the output");
        Plot::new("coh")
            .height(h)
            .legend(Legend::default())
            .x_axis_label("optical path difference (µm)")
            .y_axis_label("|γ|")
            .include_y(0.0)
            .include_y(1.05)
            .show(ui, |pl| {
                pl.line(Line::new("|γ(δ)|", pts(&coherence.0, &coherence.1)).color(PURPLE).width(1.8));
            });
    });
    note(
        ui,
        &format!(
            "Input {} · inside fiber {} · output {} · coupling {:.3} · output NA {:.2} (half-angle {:.1}°) · mean λ {:.0} nm. Approximate spectra; straight fiber in air.",
            si(t.incident_mw * 1e-3, "W"),
            si(t.guided_mw * 1e-3, "W"),
            si(t.output_mw * 1e-3, "W"),
            t.coupling,
            t.output_na,
            t.half_angle_deg,
            r.mean_wavelength_nm
        ),
    );
}

// ------------------------------------------------------------ interference

pub fn wave_tab(
    ui: &mut egui::Ui,
    r: &SimResult,
    w: &mut wave::WaveSettings,
    cross: &mut bool,
    patch: &mut Option<(String, WavePatch)>,
    tex: &mut TextureCache,
) {
    let mut changed = false;
    ui.horizontal_wrapped(|ui| {
        egui::ComboBox::from_id_salt("wmode")
            .selected_text(match w.mode {
                WaveMode::Halogen => "Halogen (incoherent mean)",
                WaveMode::Coherent => "Coherent comparison",
                WaveMode::Ensemble => "Average independent patterns",
            })
            .show_ui(ui, |ui| {
                changed |= ui.selectable_value(&mut w.mode, WaveMode::Halogen, "Halogen (incoherent mean)").changed();
                changed |= ui.selectable_value(&mut w.mode, WaveMode::Coherent, "Coherent comparison").changed();
                changed |=
                    ui.selectable_value(&mut w.mode, WaveMode::Ensemble, "Average independent patterns").changed();
            });
        ui.add_enabled_ui(w.mode != WaveMode::Halogen, |ui| {
            changed |= ui.add(egui::DragValue::new(&mut w.wavelength_nm).range(400.0..=1000.0).suffix(" nm")).changed();
        });
        ui.add_enabled_ui(w.mode == WaveMode::Ensemble, |ui| {
            egui::ComboBox::from_id_salt("avg").selected_text(format!("{} patterns", w.averages)).show_ui(ui, |ui| {
                for n in [1, 4, 16, 64] {
                    changed |= ui.selectable_value(&mut w.averages, n, format!("{n}")).changed();
                }
            });
        });
        ui.add_enabled_ui(w.mode != WaveMode::Halogen, |ui| {
            changed |=
                ui.add(egui::Slider::new(&mut w.phase_shift_um, 0.0..=1000.0).text("axial phase shift (µm)")).changed();
            if ui.button("New pattern").clicked() {
                w.seed += 1;
                changed = true;
            }
        });
        ui.checkbox(cross, "Show interference contribution (I − I_inc)");
    });
    let bench = &r.settings.bench;
    let distance = if bench.use_lens { bench.tip_to_center_mm() } else { bench.source_distance_mm };
    // With the lens, the local angular spread is that of the traced beam.
    let spread = bench.use_lens.then(|| (0.5 * r.profile.divergence90_full_deg).to_radians().sin());
    let key = format!("{:?}|{}|{}|{}|{:?}", w, r.settings.fiber.length_m, distance, r.transport.output_na, spread);
    if changed || patch.as_ref().is_none_or(|p| p.0 != key) {
        let p = wave::interference(w, r.settings.fiber.length_m, distance, r.transport.output_na, spread);
        *patch = Some((key.clone(), p));
    }
    let (_, p) = patch.as_ref().unwrap();
    let tkey = format!("{key}|{cross}");
    if tex.wave.as_ref().is_none_or(|t| t.0 != tkey) {
        let vals: Vec<f64> = p.intensity.iter().map(|&v| v as f64).collect();
        let max = vals.iter().cloned().fold(1.0, f64::max);
        let min = vals.iter().cloned().fold(f64::MAX, f64::min);
        let t = if *cross {
            let span = (max - 1.0).max(1.0 - min).max(1e-9);
            heat_texture(ui.ctx(), "wave", p.n, p.n, &vals, |v| style::diverging((v - 1.0) / span))
        } else {
            let top = max.max(2.0);
            heat_texture(ui.ctx(), "wave", p.n, p.n, &vals, |v| style::viridis(v / top))
        };
        tex.wave = Some((tkey, t));
    }
    let half = *p.axis_um.last().unwrap_or(&1.0);
    let tid = tex.wave.as_ref().unwrap().1.id();
    Plot::new("wave")
        .data_aspect(1.0)
        .height(ui.available_height() - 50.0)
        .x_axis_label("x (µm)")
        .y_axis_label("y (µm)")
        .show(ui, |pl| {
            pl.image(PlotImage::new(
                "pattern",
                tid,
                PlotPoint::new(0.0, 0.0),
                egui::vec2(2.0 * half as f32, 2.0 * half as f32),
            ));
            pl.text(Text::new("", PlotPoint::new(0.0, half * 1.06), format!("contrast σ/mean = {:.3}", p.contrast)));
        });
    note(
        ui,
        "Illustrative local field near the sample, not the exact speckle of this cable. The halogen lamp is incoherent: its mean has no fringes.",
    );
}
