//! Fiber Light Lab desktop application.
//!
//! Renderer: wgpu (DirectX 12 / Vulkan / Metal) with an automatic OpenGL
//! fallback. Force one with `--glow` / `--wgpu` or FIBERLAB_RENDERER.
#![cfg_attr(all(windows, not(debug_assertions)), windows_subsystem = "windows")]

mod app;

use eframe::egui;
use std::sync::Arc;

fn icon() -> egui::IconData {
    // 64×64: dark disk, a lens and a warm beam.
    let n = 64usize;
    let mut rgba = vec![0u8; n * n * 4];
    for y in 0..n {
        for x in 0..n {
            let (fx, fy) = (x as f32 + 0.5 - 32.0, y as f32 + 0.5 - 32.0);
            let r = (fx * fx + fy * fy).sqrt();
            let i = (y * n + x) * 4;
            let mut c = [0u8, 0, 0, 0];
            if r < 30.0 {
                c = [36, 52, 72, 255];
                // Beam: narrow on the left, spreading after the lens.
                let half = if fx < -2.0 { 1.5 } else { 1.5 + (fx + 2.0) * 0.35 };
                if fy.abs() < half && fx > -26.0 {
                    let t = 1.0 - fy.abs() / half;
                    c = [255, (180.0 + 60.0 * t) as u8, (60.0 + 80.0 * t) as u8, 255];
                }
                // Lens.
                if (fx + 2.0).abs() < 2.5 - (fy / 14.0).powi(2) * 2.0 && fy.abs() < 14.0 {
                    c = [150, 205, 240, 255];
                }
            } else if r < 31.5 {
                c = [36, 52, 72, ((31.5 - r) / 1.5 * 255.0) as u8];
            }
            rgba[i..i + 4].copy_from_slice(&c);
        }
    }
    egui::IconData { rgba, width: n as u32, height: n as u32 }
}

fn options(renderer: eframe::Renderer) -> eframe::NativeOptions {
    eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_title("Fiber Light Lab — regolith goniometer simulator")
            .with_inner_size([1480.0, 920.0])
            .with_min_inner_size([1000.0, 640.0])
            .with_icon(Arc::new(icon())),
        renderer,
        centered: true,
        ..Default::default()
    }
}

fn startup(args: &[String]) -> Result<app::Startup, String> {
    let mut s = app::Startup::default();
    let mut it = args.iter().skip(1);
    while let Some(a) = it.next() {
        match a.as_str() {
            "--settings" => {
                let p = it.next().ok_or("--settings needs a file")?;
                let text = std::fs::read_to_string(p).map_err(|e| format!("{p}: {e}"))?;
                s.settings = Some(fiber_light_lab::physics::sim::SimSettings::from_json(&text)?);
            }
            "--tab" => s.tab = it.next().cloned(),
            "--snapshot" => s.snapshot = it.next().map(std::path::PathBuf::from),
            "--run-study" => s.run_study = true,
            "--glow" | "--wgpu" => {}
            other if other.ends_with(".json") => {
                // Opening a settings file with the app (e.g. drag onto the .exe).
                let text = std::fs::read_to_string(other).map_err(|e| format!("{other}: {e}"))?;
                s.settings = Some(fiber_light_lab::physics::sim::SimSettings::from_json(&text)?);
            }
            other => return Err(format!("Unknown argument '{other}'.")),
        }
    }
    Ok(s)
}

fn run(renderer: eframe::Renderer, args: &[String]) -> eframe::Result {
    let name = match renderer {
        eframe::Renderer::Glow => "OpenGL",
        _ => "wgpu",
    };
    let start = match startup(args) {
        Ok(s) => s,
        Err(e) => {
            rfd::MessageDialog::new().set_title("Fiber Light Lab").set_description(e).show();
            app::Startup::default()
        }
    };
    eframe::run_native(
        "Fiber Light Lab",
        options(renderer),
        Box::new(move |cc| Ok(Box::new(app::FiberLabApp::new(cc, name, start)))),
    )
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let env = std::env::var("FIBERLAB_RENDERER").unwrap_or_default().to_lowercase();
    let prefer_glow = args.iter().any(|a| a == "--glow") || env == "glow" || env == "opengl";
    let order = if prefer_glow {
        [eframe::Renderer::Glow, eframe::Renderer::Wgpu]
    } else {
        [eframe::Renderer::Wgpu, eframe::Renderer::Glow]
    };
    let mut last_error = String::new();
    for r in order {
        match run(r, &args) {
            Ok(()) => return,
            Err(e) => last_error = e.to_string(),
        }
    }
    rfd::MessageDialog::new()
        .set_title("Fiber Light Lab could not start")
        .set_description(format!(
            "No graphics backend could be initialized.\n\n{last_error}\n\nThe command-line runner (fiberlab-cli) works without graphics."
        ))
        .set_level(rfd::MessageLevel::Error)
        .show();
}
