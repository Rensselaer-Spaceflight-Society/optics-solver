//! Desktop application shell: layout, background worker, export.

mod controls;
pub mod style;
mod tabs;
pub mod view3d;

use eframe::egui::{self, Color32, RichText};
use fiber_light_lab::physics::bench::BeamResult;
use fiber_light_lab::physics::goniometer::{ScanResult, ScanSettings, compute_scan};
use fiber_light_lab::physics::report;
use fiber_light_lab::physics::sample::{Finish, SampleKind, SampleModel, Simulant};
use fiber_light_lab::physics::sim::{self, SimResult, SimSettings};
use fiber_light_lab::physics::wave::{self, WavePatch, WaveSettings};
use std::path::PathBuf;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::mpsc::{Receiver, Sender, channel};
use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant};
use style::{note, si};
use tabs::{ScanView, StudyUi, TextureCache};
use view3d::{Camera, ViewOptions};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Tab {
    Bench,
    Beam,
    Scan,
    Sample,
    Study,
    Fiber,
    Wave,
    Notes,
}

impl Tab {
    const ALL: [Tab; 8] =
        [Tab::Bench, Tab::Beam, Tab::Scan, Tab::Sample, Tab::Study, Tab::Fiber, Tab::Wave, Tab::Notes];
    fn label(self) -> &'static str {
        match self {
            Tab::Bench => "3D bench",
            Tab::Beam => "Beam & spot",
            Tab::Scan => "Goniometer scan",
            Tab::Sample => "Sample BRDF",
            Tab::Study => "Experiment matrix",
            Tab::Fiber => "Fiber & spectrum",
            Tab::Wave => "Interference",
            Tab::Notes => "Model notes",
        }
    }
}

/// (x, y) samples of a 1-D curve.
type Curve2 = (Vec<f64>, Vec<f64>);

struct Job {
    generation: u64,
    settings: SimSettings,
    beam: Option<Arc<BeamResult>>,
    compare: bool,
}

struct Done {
    generation: u64,
    result: Result<SimResult, String>,
    compare: Option<(String, ScanResult)>,
}

fn other_finish_scan(r: &SimResult) -> Option<(String, ScanResult)> {
    if r.settings.sample.kind != SampleKind::Composite {
        return None;
    }
    let mut sample = r.settings.sample.clone();
    sample.composite.finish = match sample.composite.finish {
        Finish::AsCast => Finish::Abraded,
        Finish::Abraded => Finish::AsCast,
    };
    let model = SampleModel::new(&sample);
    let scan = ScanSettings { scan_rays: r.settings.quality.rays().1, ..r.settings.scan.clone() };
    let s = compute_scan(&r.settings.bench, &r.beam, &model, &r.settings.detector, &scan);
    Some((format!("{} (same coupon)", sample.composite.finish.label().to_lowercase()), s))
}

fn spawn_worker(ctx: egui::Context) -> (Sender<Job>, Receiver<Done>) {
    let (job_tx, job_rx) = channel::<Job>();
    let (done_tx, done_rx) = channel::<Done>();
    std::thread::Builder::new()
        .name("fiberlab-worker".into())
        .spawn(move || {
            while let Ok(mut job) = job_rx.recv() {
                // Only the newest request matters.
                while let Ok(newer) = job_rx.try_recv() {
                    job = newer;
                }
                let result = sim::simulate(&job.settings, job.beam.take());
                let compare = match (&result, job.compare) {
                    (Ok(r), true) => other_finish_scan(r),
                    _ => None,
                };
                if done_tx.send(Done { generation: job.generation, result, compare }).is_err() {
                    break;
                }
                ctx.request_repaint();
            }
        })
        .expect("worker thread");
    (job_tx, done_rx)
}

type StudyOutcome = Arc<Mutex<Option<Result<Vec<sim::StudyRow>, String>>>>;

struct StudyRun {
    progress: Arc<AtomicU64>,
    result: StudyOutcome,
}

/// Command-line options for the desktop app.
#[derive(Default)]
pub struct Startup {
    pub settings: Option<SimSettings>,
    pub tab: Option<String>,
    /// Save a full-window PNG once results are shown, then quit.
    pub snapshot: Option<PathBuf>,
    /// Start the experiment matrix immediately.
    pub run_study: bool,
}

pub struct FiberLabApp {
    snapshot: Option<(PathBuf, u32)>,
    settings: SimSettings,
    result: Option<SimResult>,
    error: Option<String>,
    jobs: Sender<Job>,
    done: Receiver<Done>,
    generation: u64,
    in_flight: bool,
    sent: Option<SimSettings>,
    dirty_at: Option<Instant>,
    tab: Tab,
    camera: Option<Camera>,
    view: ViewOptions,
    scene: Option<(u64, ViewOptions, view3d::Scene)>,
    view_rect: egui::Rect,
    scan_view: ScanView,
    compare: Option<(String, ScanResult)>,
    study: StudyUi,
    study_run: Option<StudyRun>,
    wave_cross: bool,
    wave_patch: Option<(String, WavePatch)>,
    coherence: Option<(u64, Curve2)>,
    textures: TextureCache,
    pending_screenshot: Option<PathBuf>,
    status: String,
    renderer: String,
    notes: Vec<NoteLine>,
}

impl FiberLabApp {
    pub fn new(cc: &eframe::CreationContext<'_>, renderer: &str, startup: Startup) -> Self {
        style::set_light_theme(&cc.egui_ctx);
        let (jobs, done) = spawn_worker(cc.egui_ctx.clone());
        let tab = startup
            .tab
            .as_deref()
            .and_then(|t| {
                Tab::ALL.into_iter().find(|x| {
                    x.label()
                        .to_lowercase()
                        .replace([' ', '&'], "")
                        .starts_with(&t.to_lowercase().replace([' ', '&'], ""))
                })
            })
            .unwrap_or(Tab::Bench);
        let mut app = FiberLabApp {
            snapshot: startup.snapshot.map(|p| (p, 0)),
            settings: startup.settings.unwrap_or_default(),
            result: None,
            error: None,
            jobs,
            done,
            generation: 0,
            in_flight: false,
            sent: None,
            dirty_at: Some(Instant::now() - Duration::from_secs(1)),
            tab,
            camera: None,
            view: ViewOptions::default(),
            scene: None,
            view_rect: egui::Rect::NOTHING,
            scan_view: ScanView::default(),
            compare: None,
            study: StudyUi::default(),
            study_run: None,
            wave_cross: false,
            wave_patch: None,
            coherence: None,
            textures: TextureCache { spot: None, wave: None },
            pending_screenshot: None,
            status: "Starting…".into(),
            renderer: renderer.to_string(),
            notes: parse_notes(report::MODEL_NOTES),
        };
        app.submit();
        if startup.run_study {
            app.start_study();
        }
        app
    }

    fn mark_dirty(&mut self) {
        self.dirty_at = Some(Instant::now());
    }

    fn submit(&mut self) {
        self.generation += 1;
        let beam =
            self.result.as_ref().filter(|r| r.settings.beam_key() == self.settings.beam_key()).map(|r| r.beam.clone());
        let _ = self.jobs.send(Job {
            generation: self.generation,
            settings: self.settings.clone(),
            beam,
            compare: self.scan_view.compare_finish,
        });
        self.sent = Some(self.settings.clone());
        self.in_flight = true;
        self.dirty_at = None;
        self.status = "Calculating…".into();
    }

    fn poll(&mut self, ctx: &egui::Context) {
        while let Ok(done) = self.done.try_recv() {
            if done.generation != self.generation {
                continue;
            }
            self.in_flight = false;
            match done.result {
                Ok(r) => {
                    self.status = format!(
                        "Ready · {:.0} ms · {}k rays traced · {} detector sub-points",
                        r.compute_ms,
                        r.beam.rays_traced / 1000,
                        r.scan.detector_subsamples
                    );
                    if self.camera.is_none() {
                        self.camera = Some(Camera::home(home_scale(&r)));
                    }
                    self.result = Some(r);
                    self.compare = done.compare;
                    self.error = None;
                }
                Err(e) => {
                    self.error = Some(e);
                    self.status = "Settings need attention".into();
                }
            }
        }
        if self.sent.as_ref() != Some(&self.settings) && self.dirty_at.is_none() {
            self.mark_dirty();
        }
        if let Some(t) = self.dirty_at {
            if !self.in_flight && t.elapsed() > Duration::from_millis(120) {
                self.submit();
            } else {
                ctx.request_repaint_after(Duration::from_millis(40));
            }
        }
        if self.in_flight {
            ctx.request_repaint_after(Duration::from_millis(60));
        }
        // Experiment matrix in its own thread.
        if let Some(run) = &self.study_run {
            let done = run.result.lock().unwrap().take();
            if let Some(res) = done {
                match res {
                    Ok(rows) => {
                        self.study.effects = sim::summarize_effects(&rows);
                        self.study.rows = rows;
                        self.study.error = None;
                    }
                    Err(e) => self.study.error = Some(e),
                }
                self.study_run = None;
            } else {
                ctx.request_repaint_after(Duration::from_millis(100));
            }
        }
    }

    fn start_study(&mut self) {
        let contents = match tabs::parse_list(&self.study.contents) {
            Ok(v) if !v.is_empty() => v,
            Ok(_) => return self.study.error = Some("Enter at least one regolith content.".into()),
            Err(e) => return self.study.error = Some(e),
        };
        let incidences = match tabs::parse_list(&self.study.incidences) {
            Ok(v) if !v.is_empty() => v,
            Ok(_) => return self.study.error = Some("Enter at least one illumination angle.".into()),
            Err(e) => return self.study.error = Some(e),
        };
        let mut sims = Vec::new();
        if self.study.lhs {
            sims.push(Simulant::Lhs1);
        }
        if self.study.lms {
            sims.push(Simulant::Lms1);
        }
        if sims.is_empty() {
            return self.study.error = Some("Pick at least one simulant.".into());
        }
        let progress = Arc::new(AtomicU64::new(0f64.to_bits()));
        let result = Arc::new(Mutex::new(None));
        let (p2, r2, settings) = (progress.clone(), result.clone(), self.settings.clone());
        std::thread::spawn(move || {
            let res =
                sim::run_study(&settings, &sims, &contents, &incidences, &|f| p2.store(f.to_bits(), Ordering::Relaxed));
            *r2.lock().unwrap() = Some(res);
        });
        self.study_run = Some(StudyRun { progress, result });
    }

    fn export(&mut self) {
        let Some(r) = &self.result else { return };
        if let Some(dir) = rfd::FileDialog::new().set_title("Choose a folder for the results").pick_folder() {
            let stamp =
                std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).map(|d| d.as_secs()).unwrap_or(0);
            let target = dir.join(format!("fiberlab_results_{stamp}"));
            let rows = (!self.study.rows.is_empty()).then_some(self.study.rows.as_slice());
            match report::write_folder(&target, r, rows) {
                Ok(files) => self.status = format!("Saved {} files to {}", files.len(), target.display()),
                Err(e) => {
                    rfd::MessageDialog::new().set_title("Could not save").set_description(e.to_string()).show();
                }
            }
        }
    }

    fn save_settings(&mut self) {
        if let Some(path) =
            rfd::FileDialog::new().add_filter("JSON", &["json"]).set_file_name("fiberlab_settings.json").save_file()
        {
            let text = serde_json::to_string_pretty(&self.settings).unwrap_or_default();
            match std::fs::write(&path, text) {
                Ok(()) => self.status = format!("Settings saved to {}", path.display()),
                Err(e) => {
                    rfd::MessageDialog::new().set_title("Could not save").set_description(e.to_string()).show();
                }
            }
        }
    }

    fn load_settings(&mut self) {
        if let Some(path) = rfd::FileDialog::new().add_filter("JSON", &["json"]).pick_file() {
            match std::fs::read_to_string(&path).map_err(|e| e.to_string()).and_then(|t| SimSettings::from_json(&t)) {
                Ok(s) => {
                    self.settings = s;
                    self.mark_dirty();
                }
                Err(e) => {
                    rfd::MessageDialog::new().set_title("Could not load settings").set_description(e).show();
                }
            }
        }
    }

    fn handle_snapshot(&mut self, ui: &egui::Ui) {
        let ready = self.result.is_some() && !self.in_flight && self.dirty_at.is_none() && self.study_run.is_none();
        let Some((path, frames)) = self.snapshot.as_mut() else { return };
        if !ready {
            ui.ctx().request_repaint();
            return;
        }
        *frames += 1;
        if *frames == 8 {
            ui.ctx().send_viewport_cmd(egui::ViewportCommand::Screenshot(egui::UserData::default()));
        }
        if *frames < 8 {
            ui.ctx().request_repaint();
            return;
        }
        let image = ui.input(|i| {
            i.raw.events.iter().find_map(|e| match e {
                egui::Event::Screenshot { image, .. } => Some(image.clone()),
                _ => None,
            })
        });
        if let Some(image) = image {
            let bytes: Vec<u8> = image.pixels.iter().flat_map(|c| c.to_array()).collect();
            let _ = image::save_buffer(
                &*path,
                &bytes,
                image.size[0] as u32,
                image.size[1] as u32,
                image::ExtendedColorType::Rgba8,
            );
            ui.ctx().send_viewport_cmd(egui::ViewportCommand::Close);
        } else {
            ui.ctx().request_repaint();
        }
    }

    fn handle_screenshot(&mut self, ui: &egui::Ui) {
        let Some(path) = self.pending_screenshot.clone() else { return };
        let image = ui.input(|i| {
            i.raw.events.iter().find_map(|e| match e {
                egui::Event::Screenshot { image, .. } => Some(image.clone()),
                _ => None,
            })
        });
        if let Some(image) = image {
            self.pending_screenshot = None;
            let crop = view3d::crop(&image, self.view_rect, ui.ctx().pixels_per_point());
            let bytes: Vec<u8> = crop.pixels.iter().flat_map(|c| c.to_array()).collect();
            let res = image::save_buffer(
                &path,
                &bytes,
                crop.size[0] as u32,
                crop.size[1] as u32,
                image::ExtendedColorType::Rgba8,
            );
            self.status = match res {
                Ok(()) => format!("Saved 3D view to {}", path.display()),
                Err(e) => format!("Could not save image: {e}"),
            };
        }
    }

    fn top_bar(&mut self, ui: &mut egui::Ui) {
        egui::MenuBar::new().ui(ui, |ui| {
            ui.menu_button("File", |ui| {
                if ui.button("Export results…").clicked() {
                    ui.close();
                    self.export();
                }
                if ui.button("Save settings…").clicked() {
                    ui.close();
                    self.save_settings();
                }
                if ui.button("Load settings…").clicked() {
                    ui.close();
                    self.load_settings();
                }
                ui.separator();
                if ui.button("Reset all settings").clicked() {
                    ui.close();
                    self.settings = SimSettings::default();
                    self.mark_dirty();
                }
            });
            ui.menu_button("View", |ui| {
                for t in Tab::ALL {
                    if ui.selectable_label(self.tab == t, t.label()).clicked() {
                        self.tab = t;
                        ui.close();
                    }
                }
            });
            ui.separator();
            ui.label(RichText::new("Fiber Light Lab").strong().color(style::ACCENT));
            note(ui, "MI-152 halogen -> GH4001 fiber -> LB1092 lens -> regolith-epoxy coupon -> detector");
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                if ui.button("Export results…").clicked() {
                    self.export();
                }
                let mut pick = self.settings.quality;
                egui::ComboBox::from_id_salt("qtop").selected_text(format!("Quality: {}", pick.label())).show_ui(
                    ui,
                    |ui| {
                        for q in [sim::Quality::Draft, sim::Quality::Standard, sim::Quality::High] {
                            ui.selectable_value(&mut pick, q, q.label());
                        }
                    },
                );
                if pick != self.settings.quality {
                    self.settings.quality = pick;
                    self.mark_dirty();
                }
            });
        });
    }

    fn bench_tab(&mut self, ui: &mut egui::Ui) {
        let Some(r) = &self.result else { return };
        let scale = home_scale(r);
        ui.horizontal_wrapped(|ui| {
            let v = &mut self.view;
            ui.checkbox(&mut v.show_rays, "Rays");
            ui.add_enabled(v.show_rays, egui::Slider::new(&mut v.ray_count, 20..=900).text("rays"));
            ui.checkbox(&mut v.show_blocked_rays, "Lost rays");
            ui.checkbox(&mut v.show_spot, "Irradiance");
            ui.checkbox(&mut v.show_lobe, "Scan lobe");
            ui.checkbox(&mut v.show_tube, "Lens tube");
            ui.checkbox(&mut v.show_labels, "Labels");
            ui.separator();
            ui.checkbox(&mut v.follow_specular, "Detector at specular");
            ui.add_enabled(
                !v.follow_specular,
                egui::Slider::new(&mut v.detector_angle_deg, -89.0..=89.0).text("detector °"),
            );
            ui.separator();
            if ui.button("Iso").clicked() {
                self.camera = Some(Camera::home(scale));
            }
            if ui.button("Side").on_hover_text("Look along the rotation axis (plane of incidence)").clicked() {
                self.camera = Some(Camera::side(scale));
            }
            if ui.button("Top").clicked() {
                self.camera = Some(Camera::top(scale));
            }
            if ui.button("Save image…").clicked()
                && let Some(path) =
                    rfd::FileDialog::new().add_filter("PNG", &["png"]).set_file_name("fiberlab_3d.png").save_file()
            {
                self.pending_screenshot = Some(path);
                ui.ctx().send_viewport_cmd(egui::ViewportCommand::Screenshot(egui::UserData::default()));
            }
        });
        let generation = self.generation;
        let rebuild = self.scene.as_ref().is_none_or(|(g, v, _)| *g != generation || *v != self.view);
        if rebuild {
            self.scene = Some((generation, self.view.clone(), view3d::build_scene(r, &self.view)));
        }
        egui::Panel::right("bench_info").resizable(true).default_size(270.0).show(ui, |ui| {
            egui::ScrollArea::vertical().show(ui, |ui| key_results(ui, r));
        });
        let cam = self.camera.get_or_insert_with(|| Camera::home(scale));
        let home = Camera::home(scale);
        let scene = &self.scene.as_ref().unwrap().2;
        let resp = view3d::show(ui, scene, cam, home);
        self.view_rect = resp.rect;
    }
}

fn home_scale(r: &SimResult) -> f64 {
    r.settings.detector.arm_radius_mm.max(r.settings.bench.tip_to_center_mm())
}

fn key_results(ui: &mut egui::Ui, r: &SimResult) {
    style::heading(ui, "Key results");
    let t = &r.transport;
    let b = &r.beam.budget;
    let m = &r.scan.metrics;
    let resp = r.settings.detector.responsivity.eval(r.mean_wavelength_nm);
    let rows: Vec<(&str, String)> = vec![
        ("Fiber output", si(t.output_mw * 1e-3, "W")),
        (
            "On coupon",
            format!("{} ({:.0} %)", si(b.on_coupon_mw * 1e-3, "W"), 100.0 * b.on_coupon_mw / t.output_mw.max(1e-300)),
        ),
        ("Spot 90 % (u × v)", format!("{:.1} × {:.1} mm", r.spot.width90_u_mm, r.spot.width90_v_mm)),
        ("Beam spread (90 %)", format!("{:.1}° full", r.profile.divergence90_full_deg)),
        ("Source-size limit", format!("{:.2}° full", r.profile.divergence_estimate_deg)),
        ("Mirror FWHM", format!("{:.2}°", r.instrument.metrics.fwhm_deg)),
        ("Detector spans", format!("{:.2}°", r.settings.detector.angular_width_deg())),
        ("Specular angle", format!("{:.2}°", m.specular_angle_deg)),
        ("Peak signal", format!("{} at {:.1}°", si(m.peak_mw * 1e-3, "W"), m.peak_angle_deg)),
        ("Peak photocurrent", si(m.peak_mw * 1e-3 * resp, "A")),
        ("Near specular", si(m.near_mean_mw * 1e-3, "W")),
        ("Away from specular", si(m.off_mean_mw * 1e-3, "W")),
        ("Near / away", format!("{:.2}", m.near_to_off_ratio)),
        ("Sample albedo", format!("{:.3}", r.brdf.albedo)),
    ];
    egui::Grid::new("keys").num_columns(2).striped(true).show(ui, |ui| {
        for (k, v) in rows {
            ui.label(RichText::new(k).color(Color32::from_gray(90)));
            ui.label(RichText::new(v).strong());
            ui.end_row();
        }
    });
    ui.add_space(6.0);
    note(
        ui,
        "Orange lobe: simulated detector signal vs angle (√ scale). Blue dashes: specular ray. Red arc: angles the detector cannot reach.",
    );
}

impl eframe::App for FiberLabApp {
    fn ui(&mut self, ui: &mut egui::Ui, _frame: &mut eframe::Frame) {
        let ctx = ui.ctx().clone();
        self.poll(&ctx);
        self.handle_screenshot(ui);
        self.handle_snapshot(ui);

        egui::Panel::top("menu").show(ui, |ui| self.top_bar(ui));
        egui::Panel::bottom("status").show(ui, |ui| {
            ui.horizontal(|ui| {
                if self.in_flight {
                    ui.spinner();
                }
                ui.label(&self.status);
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    note(
                        ui,
                        &format!(
                            "v{} · {} renderer · optical power is {}",
                            sim::APP_VERSION,
                            self.renderer,
                            match self.settings.fiber.power_mode {
                                fiber_light_lab::physics::fiber::PowerMode::Assumed => "assumed",
                                fiber_light_lab::physics::fiber::PowerMode::Mi152Rating => "a rating estimate",
                            }
                        ),
                    );
                });
            });
        });
        egui::Panel::left("controls").resizable(true).default_size(330.0).show(ui, |ui| {
            egui::ScrollArea::vertical().show(ui, |ui| {
                let mean = self.result.as_ref().map(|r| r.mean_wavelength_nm).unwrap_or(650.0);
                if controls::panel(ui, &mut self.settings, mean) {
                    self.mark_dirty();
                }
                ui.add_space(8.0);
                ui.horizontal(|ui| {
                    if ui.button("Reset all").clicked() {
                        self.settings = SimSettings::default();
                        self.mark_dirty();
                    }
                    if ui.button("Save settings…").clicked() {
                        self.save_settings();
                    }
                    if ui.button("Load…").clicked() {
                        self.load_settings();
                    }
                });
            });
        });
        egui::CentralPanel::default().show(ui, |ui| {
            ui.horizontal_wrapped(|ui| {
                for t in Tab::ALL {
                    if ui.selectable_label(self.tab == t, RichText::new(t.label()).size(14.0)).clicked() {
                        self.tab = t;
                    }
                }
            });
            ui.separator();
            if let Some(e) = &self.error {
                ui.colored_label(style::RED, format!("⚠ {e}"));
            }
            if let Some(r) = &self.result {
                for w in &r.warnings {
                    egui::Frame::new().fill(style::WARN_BG).inner_margin(4.0).corner_radius(3.0).show(ui, |ui| {
                        ui.label(RichText::new(format!("⚠ {w}")).color(style::WARN_FG).small());
                    });
                }
            }
            if self.result.is_none() {
                ui.centered_and_justified(|ui| {
                    ui.spinner();
                });
                return;
            }
            match self.tab {
                Tab::Bench => self.bench_tab(ui),
                Tab::Beam => {
                    let r = self.result.as_ref().unwrap();
                    tabs::beam_tab(ui, r, self.generation, &mut self.textures)
                }
                Tab::Scan => {
                    let r = self.result.as_ref().unwrap();
                    if tabs::scan_tab(ui, r, &mut self.scan_view, self.compare.as_ref()) {
                        self.mark_dirty();
                        self.sent = None;
                    }
                }
                Tab::Sample => tabs::brdf_tab(ui, self.result.as_ref().unwrap()),
                Tab::Study => {
                    let running = self.study_run.as_ref().map(|r| f64::from_bits(r.progress.load(Ordering::Relaxed)));
                    if tabs::study_tab(ui, &mut self.study, running) {
                        self.start_study();
                    }
                }
                Tab::Fiber => {
                    let r = self.result.as_ref().unwrap();
                    if self.coherence.as_ref().is_none_or(|c| c.0 != self.generation) {
                        self.coherence = Some((self.generation, wave::coherence_curve(&r.transport, 15.0)));
                    }
                    tabs::fiber_tab(ui, r, &self.coherence.as_ref().unwrap().1)
                }
                Tab::Wave => {
                    let r = self.result.as_ref().unwrap();
                    let mut w: WaveSettings = self.settings.wave.clone();
                    tabs::wave_tab(ui, r, &mut w, &mut self.wave_cross, &mut self.wave_patch, &mut self.textures);
                    if w != self.settings.wave {
                        // Wave settings do not affect the main simulation.
                        self.settings.wave = w.clone();
                        if let Some(s) = self.sent.as_mut() {
                            s.wave = w;
                        }
                    }
                }
                Tab::Notes => notes_view(ui, &self.notes),
            }
        });
    }
}

// ------------------------------------------------------------ model notes

enum NoteLine {
    H1(String),
    H2(String),
    H3(String),
    Bullet(String),
    Mono(String),
    Para(String),
    Gap,
}

fn parse_notes(md: &str) -> Vec<NoteLine> {
    let mut out = Vec::new();
    let mut para = String::new();
    let mut in_code = false;
    let flush = |para: &mut String, out: &mut Vec<NoteLine>| {
        if !para.trim().is_empty() {
            out.push(NoteLine::Para(std::mem::take(para).trim().to_string()));
        }
        para.clear();
    };
    for line in md.lines() {
        if line.trim_start().starts_with("```") {
            flush(&mut para, &mut out);
            in_code = !in_code;
            continue;
        }
        if in_code || line.starts_with("    ") || line.starts_with('|') {
            flush(&mut para, &mut out);
            out.push(NoteLine::Mono(line.to_string()));
            continue;
        }
        let t = line.trim();
        if t.is_empty() {
            flush(&mut para, &mut out);
            out.push(NoteLine::Gap);
        } else if let Some(h) = t.strip_prefix("### ") {
            flush(&mut para, &mut out);
            out.push(NoteLine::H3(h.into()));
        } else if let Some(h) = t.strip_prefix("## ") {
            flush(&mut para, &mut out);
            out.push(NoteLine::H2(h.into()));
        } else if let Some(h) = t.strip_prefix("# ") {
            flush(&mut para, &mut out);
            out.push(NoteLine::H1(h.into()));
        } else if let Some(b) = t.strip_prefix("- ").or_else(|| t.strip_prefix("* ")) {
            flush(&mut para, &mut out);
            out.push(NoteLine::Bullet(b.into()));
        } else {
            if !para.is_empty() {
                para.push(' ');
            }
            para.push_str(t);
        }
    }
    flush(&mut para, &mut out);
    out
}

fn strip_md(s: &str) -> String {
    s.replace("**", "").replace('`', "")
}

fn rich_with_links(ui: &mut egui::Ui, text: &str) {
    // Render [label](url) as hyperlinks, the rest as wrapped text.
    let text = strip_md(text);
    let mut rest = text.as_str();
    ui.horizontal_wrapped(|ui| {
        ui.spacing_mut().item_spacing.x = 0.0;
        while let Some(start) = rest.find('[') {
            let Some(mid) = rest[start..].find("](") else { break };
            let Some(end) = rest[start + mid..].find(')') else { break };
            ui.label(&rest[..start]);
            let label = &rest[start + 1..start + mid];
            let url = &rest[start + mid + 2..start + mid + end];
            ui.hyperlink_to(label, url);
            rest = &rest[start + mid + end + 1..];
        }
        ui.label(rest);
    });
}

fn notes_view(ui: &mut egui::Ui, notes: &[NoteLine]) {
    egui::ScrollArea::vertical().auto_shrink(false).show(ui, |ui| {
        ui.set_max_width(980.0);
        for n in notes {
            match n {
                NoteLine::H1(t) => {
                    ui.label(RichText::new(strip_md(t)).size(22.0).strong());
                }
                NoteLine::H2(t) => {
                    ui.add_space(6.0);
                    ui.label(RichText::new(strip_md(t)).size(17.0).strong().color(style::ACCENT));
                }
                NoteLine::H3(t) => {
                    ui.label(RichText::new(strip_md(t)).size(14.5).strong());
                }
                NoteLine::Bullet(t) => {
                    ui.horizontal_wrapped(|ui| {
                        ui.label("  •");
                        rich_with_links(ui, t);
                    });
                }
                NoteLine::Mono(t) => {
                    ui.label(RichText::new(t).monospace().size(12.0));
                }
                NoteLine::Para(t) => rich_with_links(ui, t),
                NoteLine::Gap => ui.add_space(4.0),
            }
        }
    });
}
