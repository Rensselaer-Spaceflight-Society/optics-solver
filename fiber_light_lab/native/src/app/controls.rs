//! Left-hand parameter panel. Every edit returns `true` so the app can
//! schedule a recalculation.

use super::style::{self, note};
use eframe::egui::{self, DragValue, RichText};
use fiber_light_lab::physics::bench::CouponShape;
use fiber_light_lab::physics::fiber::{self, Band, PowerMode};
use fiber_light_lab::physics::goniometer::{DetectorPreset, DetectorShape};
use fiber_light_lab::physics::lens::{Coating, LensSpec};
use fiber_light_lab::physics::math::Curve;
use fiber_light_lab::physics::sample::{CompositeSettings, Finish, HapkeParams, SampleKind, Simulant, hockey_stick_c};
use fiber_light_lab::physics::sim::{Quality, SimSettings};

const IN: f64 = 25.4;

fn row(ui: &mut egui::Ui, label: &str, tip: &str, add: impl FnOnce(&mut egui::Ui) -> bool) -> bool {
    let r = ui.label(label);
    if !tip.is_empty() {
        r.on_hover_text(tip);
    }
    let changed = add(ui);
    ui.end_row();
    changed
}

fn drag(ui: &mut egui::Ui, v: &mut f64, lo: f64, hi: f64, speed: f64, suffix: &str, decimals: usize) -> bool {
    ui.add(DragValue::new(v).range(lo..=hi).speed(speed).suffix(suffix).max_decimals(decimals)).changed()
}

fn grid(ui: &mut egui::Ui, id: &str, f: impl FnOnce(&mut egui::Ui) -> bool) -> bool {
    egui::Grid::new(id).num_columns(2).spacing([10.0, 6.0]).striped(false).show(ui, f).inner
}

pub fn panel(ui: &mut egui::Ui, s: &mut SimSettings, mean_wl: f64) -> bool {
    let mut changed = false;
    changed |= source_section(ui, s);
    changed |= lens_section(ui, s, mean_wl);
    changed |= geometry_section(ui, s);
    changed |= detector_section(ui, s);
    changed |= sample_section(ui, s);
    changed |= scan_section(ui, s);
    changed
}

fn source_section(ui: &mut egui::Ui, s: &mut SimSettings) -> bool {
    let mut c = false;
    egui::CollapsingHeader::new(RichText::new("Source & fiber").strong()).default_open(true).show(ui, |ui| {
        note(ui, "Dolan-Jenner MI-152 (EKE 150 W, 21 V) -> Mitsubishi GH4001 POF, 980 µm core, NA 0.50, n = 1.49");
        let f = &mut s.fiber;
        ui.horizontal(|ui| {
            c |= ui.radio_value(&mut f.power_mode, PowerMode::Assumed, "Assumed power").changed();
            c |= ui.radio_value(&mut f.power_mode, PowerMode::Mi152Rating, "MI-152 rating").changed();
        });
        c |= grid(ui, "src", |ui| {
            let mut k = false;
            if f.power_mode == PowerMode::Assumed {
                k |= row(ui, "Input power", "Band-limited optical power at the fiber face (assumed).", |ui| {
                    drag(ui, &mut f.incident_mw, 0.0, 10_000.0, 0.01, " mW", 4)
                });
            } else {
                let mut pct = f.rating_fraction * 100.0;
                k |= row(ui, "Lamp setting", "Fraction of rated output (linear scaling only).", |ui| {
                    let ch = drag(ui, &mut pct, 0.0, 100.0, 0.5, " %", 1);
                    f.rating_fraction = pct / 100.0;
                    ch
                });
            }
            k |= row(ui, "Band", "Wavelength range included in the calculation.", |ui| {
                let mut ch = false;
                egui::ComboBox::from_id_salt("band").selected_text(f.band.label()).show_ui(ui, |ui| {
                    ch |= ui.selectable_value(&mut f.band, Band::Nm400To1000, Band::Nm400To1000.label()).changed();
                    ch |= ui.selectable_value(&mut f.band, Band::Nm400To700, Band::Nm400To700.label()).changed();
                });
                ch
            });
            let mut inches = f.length_m * 1000.0 / IN;
            k |= row(ui, "Cable length", "Fiber length (the original setup is 10 in).", |ui| {
                let ch = drag(ui, &mut inches, 0.0, 400.0, 0.05, " in", 2);
                f.length_m = inches * IN / 1000.0;
                ch
            });
            k |= row(ui, "Launch NA", "Cone of light entering the fiber. Above 0.50 overfills it.", |ui| {
                drag(ui, &mut f.launch_na, 0.02, 0.95, 0.005, "", 3)
            });
            k
        });
        if f.power_mode == PowerMode::Mi152Rating {
            note(
                ui,
                &format!(
                    "350,000 fc at the port gives about {:.1} mW on the core face in {} (uniform overfill). Manufacturer claim, not a calibration.",
                    fiber::mi152_rated_core_power_mw(f.band),
                    f.band.label()
                ),
            );
        }
        egui::CollapsingHeader::new("Launch details").id_salt("launch").show(ui, |ui| {
            c |= grid(ui, "launch_g", |ui| {
                let mut k = false;
                k |= row(ui, "Spot radius", "Gaussian 1/e² radius at the input face.", |ui| {
                    drag(ui, &mut f.spot_radius_mm, 0.02, 3.0, 0.005, " mm", 3)
                });
                k |= row(ui, "Spot offset", "Lateral offset of the launch spot.", |ui| {
                    drag(ui, &mut f.offset_mm, 0.0, 2.0, 0.005, " mm", 3)
                });
                k |= row(ui, "Extra loss", "Measured connector / bend loss to add.", |ui| {
                    drag(ui, &mut f.extra_loss_db, 0.0, 20.0, 0.05, " dB", 2)
                });
                k
            });
        });
    });
    c
}

fn lens_section(ui: &mut egui::Ui, s: &mut SimSettings, mean_wl: f64) -> bool {
    let mut c = false;
    egui::CollapsingHeader::new(RichText::new("Lens (Thorlabs LB1092)").strong()).default_open(true).show(ui, |ui| {
        let b = &mut s.bench;
        c |= ui.checkbox(&mut b.use_lens, "Lens between fiber and sample").changed();
        if !b.use_lens {
            note(ui, "Bare fiber: as in the original Python app.");
            return;
        }
        note(ui, "N-BK7 bi-convex, Ø12.7 mm, f = 15.0 mm @ 587.6 nm, R = ±14.6 mm, CT 4.7 mm");
        c |= grid(ui, "lens", |ui| {
            let mut k = false;
            k |= row(ui, "Fiber to lens", "Fiber exit face to the lens front vertex.", |ui| {
                drag(ui, &mut b.fiber_to_lens_mm, 0.5, 200.0, 0.01, " mm", 3)
            });
            k |= row(ui, "Coating", "", |ui| {
                let mut ch = false;
                egui::ComboBox::from_id_salt("coat").selected_text(b.lens.coating.label()).show_ui(ui, |ui| {
                    for co in [Coating::Uncoated, Coating::ArA, Coating::ArB, Coating::ArAB, Coating::Ideal] {
                        ch |= ui.selectable_value(&mut b.lens.coating, co, co.label()).changed();
                    }
                });
                ch
            });
            k |= row(ui, "Mount aperture", "Clear opening of the lens mount / retaining ring.", |ui| {
                drag(ui, &mut b.lens.mount_aperture_mm, 1.0, 50.0, 0.05, " mm", 2)
            });
            k |= row(ui, "Fiber offset x", "Tip misalignment across the axis (in the plane of incidence).", |ui| {
                drag(ui, &mut b.fiber_offset_x_mm, -5.0, 5.0, 0.01, " mm", 3)
            });
            k |= row(ui, "Fiber offset y", "Tip misalignment along the rotation axis.", |ui| {
                drag(ui, &mut b.fiber_offset_y_mm, -5.0, 5.0, 0.01, " mm", 3)
            });
            k
        });
        ui.horizontal(|ui| {
            if ui.button("Collimate").on_hover_text("Fiber at the paraxial front focus (587.6 nm).").clicked() {
                b.fiber_to_lens_mm = b.lens.collimating_distance(587.6);
                c = true;
            }
            if ui.button("Image fiber on sample").on_hover_text("Paraxial conjugate that images the core onto the rotation center.").clicked()
                && let Some(d) = b.lens.imaging_distance(587.6, b.source_distance_mm) {
                    b.fiber_to_lens_mm = d;
                    c = true;
                }
        });
        let l = &b.lens;
        note(
            ui,
            &format!(
                "EFL {:.3} mm, BFL {:.3} mm @ 587.6 nm · EFL {:.3} mm @ {:.0} nm (spectrum mean) · catalog f {:.1}, BFL {:.1}",
                l.efl(587.6),
                l.bfl(587.6),
                l.efl(mean_wl),
                mean_wl,
                l.catalog_efl_mm,
                l.catalog_bfl_mm
            ),
        );
        egui::CollapsingHeader::new("Lens prescription").id_salt("presc").show(ui, |ui| {
            c |= grid(ui, "presc_g", |ui| {
                let l = &mut b.lens;
                let mut k = false;
                k |= row(ui, "R1", "Front radius (+ convex toward fiber).", |ui| drag(ui, &mut l.r1_mm, -500.0, 500.0, 0.01, " mm", 3));
                k |= row(ui, "R2", "Rear radius (− convex toward sample).", |ui| drag(ui, &mut l.r2_mm, -500.0, 500.0, 0.01, " mm", 3));
                k |= row(ui, "Center thickness", "", |ui| drag(ui, &mut l.center_thickness_mm, 0.2, 30.0, 0.01, " mm", 3));
                k |= row(ui, "Diameter", "", |ui| drag(ui, &mut l.diameter_mm, 1.0, 50.8, 0.01, " mm", 2));
                k
            });
            if ui.button("Reset to LB1092 catalog values").clicked() {
                let coating = b.lens.coating;
                b.lens = LensSpec { coating, ..LensSpec::lb1092() };
                c = true;
            }
            note(ui, &format!("Edge thickness {:.2} mm (catalog 1.8 mm)", b.lens.edge_thickness_mm()));
        });
    });
    c
}

fn geometry_section(ui: &mut egui::Ui, s: &mut SimSettings) -> bool {
    let mut c = false;
    egui::CollapsingHeader::new(RichText::new("Goniometer geometry").strong()).default_open(true).show(ui, |ui| {
        let b = &mut s.bench;
        ui.horizontal(|ui| {
            ui.label("Illumination");
            for a in [30.0, 60.0] {
                if ui.selectable_label(b.incidence_deg == a, format!("{a:.0}°")).clicked() {
                    b.incidence_deg = a;
                    c = true;
                }
            }
            c |= drag(ui, &mut b.incidence_deg, 0.0, 85.0, 0.25, "°", 2);
        });
        c |= grid(ui, "geom", |ui| {
            let mut k = false;
            let what = if b.use_lens { "Lens to sample" } else { "Fiber tip to sample" };
            k |= row(ui, what, "Source arm length to the rotation center (lens center when a lens is used).", |ui| {
                drag(ui, &mut b.source_distance_mm, 3.0, 1000.0, 0.1, " mm", 2)
            });
            k |= row(ui, "Detector radius", "Rotation center to the detector's active surface.", |ui| {
                drag(ui, &mut s.detector.arm_radius_mm, 10.0, 2000.0, 0.5, " mm", 1)
            });
            k |= row(ui, "Source housing r", "Outer radius of the lens tube (sets the blocked zone).", |ui| {
                drag(ui, &mut b.source_housing_radius_mm, 0.5, 60.0, 0.1, " mm", 1)
            });
            k |= row(ui, "Sample tilt", "Misalignment about the rotation axis.", |ui| {
                drag(ui, &mut b.sample_tilt_deg, -20.0, 20.0, 0.05, "°", 2)
            });
            k |= row(ui, "Sample roll", "Misalignment out of the plane of incidence.", |ui| {
                drag(ui, &mut b.sample_roll_deg, -20.0, 20.0, 0.05, "°", 2)
            });
            k
        });
        note(
            ui,
            &format!(
                "{} = {:.2} in",
                if b.use_lens { "Lens to sample" } else { "Tip to sample" },
                b.source_distance_mm / IN
            ),
        );
        egui::CollapsingHeader::new("Coupon").id_salt("coupon").show(ui, |ui| {
            let cp = &mut b.coupon;
            ui.horizontal(|ui| {
                c |= ui.radio_value(&mut cp.shape, CouponShape::Rectangle, "Rectangle").changed();
                c |= ui.radio_value(&mut cp.shape, CouponShape::Disk, "Disk").changed();
            });
            c |= grid(ui, "cpn", |ui| {
                let mut k = false;
                let l1 = if cp.shape == CouponShape::Disk { "Diameter" } else { "Length (in plane)" };
                k |= row(ui, l1, "", |ui| drag(ui, &mut cp.length_mm, 1.0, 300.0, 0.1, " mm", 2));
                if cp.shape == CouponShape::Rectangle {
                    k |= row(ui, "Width (along axis)", "", |ui| drag(ui, &mut cp.width_mm, 1.0, 300.0, 0.1, " mm", 2));
                }
                k |= row(ui, "Thickness", "", |ui| drag(ui, &mut cp.thickness_mm, 0.1, 50.0, 0.1, " mm", 2));
                k
            });
            note(ui, "Placeholder 25.4 × 25.4 × 5 mm. Enter the cast coupon size.");
        });
    });
    c
}

fn detector_section(ui: &mut egui::Ui, s: &mut SimSettings) -> bool {
    let mut c = false;
    egui::CollapsingHeader::new(RichText::new("Detector").strong()).default_open(false).show(ui, |ui| {
        let d = &mut s.detector;
        let mut preset = d.preset;
        egui::ComboBox::from_id_salt("detp").width(230.0).selected_text(preset.label()).show_ui(ui, |ui| {
            for p in DetectorPreset::ALL {
                ui.selectable_value(&mut preset, p, p.label());
            }
        });
        if preset != d.preset {
            d.apply_preset(preset);
            c = true;
        }
        note(ui, "Model not yet verified — pick the real detector and confirm its active area.");
        c |= grid(ui, "det", |ui| {
            let mut k = false;
            k |= row(ui, "Shape", "", |ui| {
                let mut ch = ui.radio_value(&mut d.shape, DetectorShape::Square, "Square").changed();
                ch |= ui.radio_value(&mut d.shape, DetectorShape::Circle, "Circle").changed();
                if ch {
                    d.preset = DetectorPreset::Custom;
                }
                ch
            });
            let lbl = if d.shape == DetectorShape::Square { "Active side" } else { "Active diameter" };
            k |= row(ui, lbl, "", |ui| {
                let ch = drag(ui, &mut d.size_mm, 0.1, 50.0, 0.01, " mm", 3);
                if ch {
                    d.preset = DetectorPreset::Custom;
                }
                ch
            });
            k |= row(ui, "Housing radius", "Outer radius of the detector body (blocked zone).", |ui| {
                drag(ui, &mut d.housing_radius_mm, 0.5, 60.0, 0.1, " mm", 1)
            });
            k |= row(ui, "Meter λ setting", "Wavelength a power meter is set to.", |ui| {
                drag(ui, &mut d.meter_wavelength_nm, 350.0, 1100.0, 1.0, " nm", 0)
            });
            k |= row(ui, "Transimpedance", "Amplifier gain or load for the voltage estimate.", |ui| {
                ui.add(DragValue::new(&mut d.transimpedance_v_per_a).range(1.0..=1e10).speed(1000.0).suffix(" V/A"))
                    .changed()
            });
            k
        });
        note(
            ui,
            &format!(
                "Area {:.2} mm² · spans {:.2}° · Ω = {:.2e} sr · {}",
                d.area_mm2(),
                d.angular_width_deg(),
                d.solid_angle_sr(),
                d.responsivity_label
            ),
        );
        ui.horizontal(|ui| {
            if ui.button("Load responsivity CSV…").on_hover_text("Two columns: wavelength_nm, A/W").clicked()
                && let Some(path) = rfd::FileDialog::new().add_filter("CSV", &["csv", "txt"]).pick_file()
            {
                match std::fs::read_to_string(&path).map_err(|e| e.to_string()).and_then(|t| Curve::from_csv(&t)) {
                    Ok(curve) => {
                        d.responsivity = curve;
                        d.responsivity_label =
                            path.file_name().map(|n| n.to_string_lossy().into_owned()).unwrap_or_default();
                        c = true;
                    }
                    Err(e) => {
                        rfd::MessageDialog::new().set_title("Could not read responsivity").set_description(e).show();
                    }
                }
            }
            if ui.button("Typical Si").clicked() {
                d.responsivity = fiber_light_lab::physics::data::si_responsivity().clone();
                d.responsivity_label = "Typical Si photodiode (approximate)".into();
                c = true;
            }
        });
    });
    c
}

fn sample_section(ui: &mut egui::Ui, s: &mut SimSettings) -> bool {
    let mut c = false;
    egui::CollapsingHeader::new(RichText::new("Sample").strong()).default_open(true).show(ui, |ui| {
        let sm = &mut s.sample;
        ui.horizontal(|ui| {
            c |= ui
                .selectable_value(&mut sm.kind, SampleKind::Mirror, "Mirror")
                .on_hover_text("Alignment and resolution check")
                .changed();
            c |= ui.selectable_value(&mut sm.kind, SampleKind::LambertianReference, "Diffuse ref.").changed();
            c |= ui.selectable_value(&mut sm.kind, SampleKind::Composite, "Composite").changed();
        });
        match sm.kind {
            SampleKind::Mirror => {
                c |= grid(ui, "mir", |ui| {
                    row(ui, "Reflectance", "Protected aluminum ≈ 0.90", |ui| {
                        drag(ui, &mut sm.mirror_reflectance, 0.0, 1.0, 0.005, "", 3)
                    })
                });
            }
            SampleKind::LambertianReference => {
                c |= grid(ui, "lam", |ui| {
                    row(ui, "Reflectance", "PTFE / Spectralon ≈ 0.99", |ui| {
                        drag(ui, &mut sm.lambertian_reflectance, 0.0, 1.0, 0.005, "", 3)
                    })
                });
            }
            SampleKind::Composite => c |= composite_controls(ui, &mut sm.composite),
        }
    });
    c
}

fn composite_controls(ui: &mut egui::Ui, cp: &mut CompositeSettings) -> bool {
    let mut c = false;
    let mut sim = cp.simulant;
    egui::ComboBox::from_id_salt("sim").selected_text(sim.label()).show_ui(ui, |ui| {
        for s in [Simulant::Lhs1, Simulant::Lms1, Simulant::Custom] {
            ui.selectable_value(&mut sim, s, s.label());
        }
    });
    if sim != cp.simulant {
        cp.apply_simulant(sim);
        c = true;
    }
    ui.horizontal(|ui| {
        ui.label("Regolith");
        for p in [0.0, 30.0, 45.0, 60.0] {
            if ui.selectable_label(cp.regolith_mass_percent == p, format!("{p:.0}%")).clicked() {
                cp.regolith_mass_percent = p;
                c = true;
            }
        }
    });
    c |= ui.add(egui::Slider::new(&mut cp.regolith_mass_percent, 0.0..=70.0).suffix(" wt%").step_by(0.5)).changed();
    ui.horizontal(|ui| {
        c |= ui.radio_value(&mut cp.finish, Finish::AsCast, "As-cast").changed();
        c |= ui.radio_value(&mut cp.finish, Finish::Abraded, "Abraded").changed();
    });
    note(
        ui,
        &format!(
            "Volume fraction {:.3} · exposed grains {:.3} · body R {:.3} · RMS slope {:.1}°",
            cp.volume_fraction(),
            cp.exposed_fraction(),
            cp.body_reflectance(),
            cp.rms_slope_deg()
        ),
    );
    egui::CollapsingHeader::new("Composite model parameters").id_salt("compadv").show(ui, |ui| {
        note(ui, "Assumptions, not measurements. Measure roughness (Rdq) and resin index if you can.");
        c |= grid(ui, "comp", |ui| {
            let mut k = false;
            k |= row(ui, "Resin index", "Cured epoxy, typical 1.54–1.57.", |ui| drag(ui, &mut cp.resin_index, 1.3, 1.8, 0.001, "", 4));
            k |= row(ui, "Resin density", "", |ui| drag(ui, &mut cp.resin_density_g_cm3, 0.5, 3.0, 0.005, " g/cm³", 3));
            k |= row(ui, "Grain density", "Estimate from mineralogy.", |ui| drag(ui, &mut cp.grain_density_g_cm3, 1.0, 6.0, 0.005, " g/cm³", 3));
            k |= row(ui, "As-cast RMS slope", "atan(Rdq) of the cast resin skin.", |ui| drag(ui, &mut cp.slope_as_cast_deg, 0.0, 45.0, 0.05, "°", 2));
            k |= row(ui, "Abraded RMS slope", "atan(Rdq) after controlled abrasion.", |ui| drag(ui, &mut cp.slope_abraded_deg, 0.0, 45.0, 0.1, "°", 2));
            k |= row(ui, "Exposure factor", "Exposed-grain area = factor × volume fraction (Delesse).", |ui| drag(ui, &mut cp.exposure_factor, 0.0, 3.0, 0.01, "", 3));
            k |= row(ui, "As-cast exposed", "Grain area fraction on the as-cast face.", |ui| drag(ui, &mut cp.as_cast_exposed_fraction, 0.0, 1.0, 0.005, "", 3));
            k |= row(ui, "Embedded darkening", "Body reflectance / dry-powder reflectance.", |ui| drag(ui, &mut cp.embedded_darkening, 0.0, 1.0, 0.01, "", 3));
            k |= row(ui, "Backing reflectance", "Behind transparent neat epoxy.", |ui| drag(ui, &mut cp.backing_reflectance, 0.0, 1.0, 0.005, "", 3));
            k |= row(ui, "Back-face reflection", "Neat epoxy only.", |ui| ui.checkbox(&mut cp.back_face_reflection, "").changed());
            k |= row(ui, "Back-face slope", "", |ui| drag(ui, &mut cp.back_face_slope_deg, 0.0, 45.0, 0.05, "°", 2));
            k
        });
        style::heading(ui, "Exposed grains — Hapke");
        note(ui, "Presets: Curtis et al. 2024 fits (Apollo 16 for LHS-1, Apollo 11 for LMS-1); c from Hapke's hockey-stick relation.");
        c |= grid(ui, "hapke", |ui| {
            let h: &mut HapkeParams = &mut cp.hapke;
            let mut k = false;
            k |= row(ui, "w (albedo)", "", |ui| drag(ui, &mut h.w, 0.0, 0.999, 0.002, "", 3));
            k |= row(ui, "b", "", |ui| drag(ui, &mut h.b, 0.0, 0.95, 0.002, "", 3));
            k |= row(ui, "c", "Backscatter fraction (+ = backscattering).", |ui| drag(ui, &mut h.c, -1.0, 1.0, 0.005, "", 3));
            k |= row(ui, "B_S0", "", |ui| drag(ui, &mut h.b_s0, 0.0, 5.0, 0.01, "", 2));
            k |= row(ui, "h_S", "", |ui| drag(ui, &mut h.h_s, 0.001, 1.0, 0.001, "", 3));
            k |= row(ui, "theta-bar", "Macroscopic roughness.", |ui| drag(ui, &mut h.theta_bar_deg, 0.0, 45.0, 0.1, "°", 2));
            k |= row(ui, "φ (filling)", "", |ui| drag(ui, &mut h.filling_factor, 0.0, 0.74, 0.005, "", 3));
            k
        });
        ui.horizontal(|ui| {
            if ui.button("c from hockey stick").clicked() {
                cp.hapke.c = hockey_stick_c(cp.hapke.b);
                c = true;
            }
            if ui.button("Reset composite").clicked() {
                let keep = (cp.simulant, cp.regolith_mass_percent, cp.finish);
                *cp = CompositeSettings::default();
                cp.apply_simulant(keep.0);
                cp.regolith_mass_percent = keep.1;
                cp.finish = keep.2;
                c = true;
            }
        });
    });
    c
}

fn scan_section(ui: &mut egui::Ui, s: &mut SimSettings) -> bool {
    let mut c = false;
    egui::CollapsingHeader::new(RichText::new("Scan & accuracy").strong()).default_open(false).show(ui, |ui| {
        let sc = &mut s.scan;
        c |= grid(ui, "scan", |ui| {
            let mut k = false;
            k |= row(ui, "From", "", |ui| drag(ui, &mut sc.start_deg, -89.0, 89.0, 0.5, "°", 2));
            k |= row(ui, "To", "", |ui| drag(ui, &mut sc.stop_deg, -89.0, 89.0, 0.5, "°", 2));
            k |= row(ui, "Step", "", |ui| drag(ui, &mut sc.step_deg, 0.05, 10.0, 0.05, "°", 2));
            k |= row(ui, "Near window ±", "Around the specular angle.", |ui| {
                drag(ui, &mut sc.near_halfwidth_deg, 0.1, 30.0, 0.1, "°", 2)
            });
            k |= row(ui, "Away beyond", "Distance from specular for the 'away' window.", |ui| {
                drag(ui, &mut sc.off_min_deg, 1.0, 80.0, 0.5, "°", 2)
            });
            k |= row(ui, "Quality", "Rays traced / used in the scan.", |ui| {
                let mut ch = false;
                egui::ComboBox::from_id_salt("q").selected_text(s.quality.label()).show_ui(ui, |ui| {
                    for q in [Quality::Draft, Quality::Standard, Quality::High] {
                        let (a, b) = q.rays();
                        ch |= ui
                            .selectable_value(
                                &mut s.quality,
                                q,
                                format!("{} ({}k / {}k rays)", q.label(), a / 1000, b as f64 / 1000.0),
                            )
                            .changed();
                    }
                });
                ch
            });
            k |= row(ui, "Random seed", "", |ui| ui.add(DragValue::new(&mut s.seed).range(0..=1_000_000)).changed());
            k
        });
    });
    c
}
