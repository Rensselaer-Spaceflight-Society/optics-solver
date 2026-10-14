//! Colors, colormaps and number formatting for the desktop app.

use eframe::egui::{self, Color32};

pub const ACCENT: Color32 = Color32::from_rgb(0x32, 0x6b, 0x91);
pub const ORANGE: Color32 = Color32::from_rgb(0xc2, 0x6a, 0x1b);
pub const GREEN: Color32 = Color32::from_rgb(0x2f, 0x8a, 0x4e);
pub const RED: Color32 = Color32::from_rgb(0xb8, 0x3b, 0x3b);
pub const PURPLE: Color32 = Color32::from_rgb(0x7a, 0x4e, 0xa8);
pub const GRAY: Color32 = Color32::from_rgb(0x88, 0x88, 0x88);
pub const WARN_BG: Color32 = Color32::from_rgb(0xfd, 0xf3, 0xd6);
pub const WARN_FG: Color32 = Color32::from_rgb(0x6b, 0x4a, 0x00);

/// Matplotlib viridis, 9 control points.
const VIRIDIS: [[f32; 3]; 9] = [
    [0.267, 0.005, 0.329],
    [0.283, 0.141, 0.458],
    [0.254, 0.265, 0.530],
    [0.207, 0.372, 0.553],
    [0.164, 0.471, 0.558],
    [0.128, 0.567, 0.551],
    [0.135, 0.659, 0.518],
    [0.478, 0.821, 0.319],
    [0.993, 0.906, 0.144],
];

pub fn viridis(t: f64) -> Color32 {
    let t = t.clamp(0.0, 1.0) as f32 * (VIRIDIS.len() - 1) as f32;
    let i = (t.floor() as usize).min(VIRIDIS.len() - 2);
    let f = t - i as f32;
    let a = VIRIDIS[i];
    let b = VIRIDIS[i + 1];
    let c = |k: usize| ((a[k] + f * (b[k] - a[k])) * 255.0).round() as u8;
    Color32::from_rgb(c(0), c(1), c(2))
}

/// Diverging blue-white-red for signed values in [-1, 1].
pub fn diverging(t: f64) -> Color32 {
    let t = t.clamp(-1.0, 1.0) as f32;
    let (r, g, b) = if t < 0.0 {
        let s = -t;
        (1.0 - 0.8 * s, 1.0 - 0.6 * s, 1.0 - 0.25 * s)
    } else {
        (1.0 - 0.25 * t, 1.0 - 0.75 * t, 1.0 - 0.8 * t)
    };
    Color32::from_rgb((r * 255.0) as u8, (g * 255.0) as u8, (b * 255.0) as u8)
}

/// Approximate perceived color of monochromatic light (display only).
pub fn wavelength_color(nm: f64) -> Color32 {
    let (r, g, b) = match nm {
        w if w < 440.0 => (-(w - 440.0) / 60.0, 0.0, 1.0),
        w if w < 490.0 => (0.0, (w - 440.0) / 50.0, 1.0),
        w if w < 510.0 => (0.0, 1.0, -(w - 510.0) / 20.0),
        w if w < 580.0 => ((w - 510.0) / 70.0, 1.0, 0.0),
        w if w < 645.0 => (1.0, -(w - 645.0) / 65.0, 0.0),
        _ => (1.0, 0.0, 0.0),
    };
    // Fade toward the ends of vision; render near-IR as deep red-brown.
    let fade = match nm {
        w if w < 420.0 => 0.35 + 0.65 * (w - 380.0).max(0.0) / 40.0,
        w if w > 700.0 => (0.35 + 0.65 * (780.0 - w).max(0.0) / 80.0).max(0.45),
        _ => 1.0,
    };
    let c = |v: f64| ((v * fade).clamp(0.0, 1.0).powf(0.8) * 255.0) as u8;
    if nm > 760.0 {
        return Color32::from_rgb(120, 30, 30);
    }
    Color32::from_rgb(c(r), c(g), c(b))
}

pub fn with_alpha(c: Color32, a: f32) -> Color32 {
    Color32::from_rgba_unmultiplied(c.r(), c.g(), c.b(), (a.clamp(0.0, 1.0) * 255.0) as u8)
}

/// Format a value with an SI prefix: si(1.46e-6, "W") -> "1.46 µW".
pub fn si(value: f64, unit: &str) -> String {
    if !value.is_finite() {
        return format!("— {unit}");
    }
    if value == 0.0 {
        return format!("0 {unit}");
    }
    let prefixes = [(1e-15, "f"), (1e-12, "p"), (1e-9, "n"), (1e-6, "µ"), (1e-3, "m"), (1.0, ""), (1e3, "k")];
    let a = value.abs();
    let mut chosen = prefixes[0];
    for p in prefixes {
        if a >= p.0 * 0.9995 {
            chosen = p;
        }
    }
    let v = value / chosen.0;
    let digits = if v.abs() >= 100.0 {
        0
    } else if v.abs() >= 10.0 {
        1
    } else {
        2
    };
    format!("{v:.digits$} {}{unit}", chosen.1)
}

/// Pick a display scale for a series: returns (multiplier, unit label).
pub fn si_scale(max_abs: f64, unit: &str) -> (f64, String) {
    let prefixes = [(1e-15, "f"), (1e-12, "p"), (1e-9, "n"), (1e-6, "µ"), (1e-3, "m"), (1.0, "")];
    let mut chosen = (1.0, "");
    for p in prefixes {
        if max_abs >= p.0 {
            chosen = p;
        }
    }
    (1.0 / chosen.0, format!("{}{unit}", chosen.1))
}

pub fn heading(ui: &mut egui::Ui, text: &str) {
    ui.add_space(4.0);
    ui.label(egui::RichText::new(text).strong().size(13.5));
}

pub fn note(ui: &mut egui::Ui, text: &str) {
    ui.label(egui::RichText::new(text).small().color(Color32::from_gray(105)));
}

pub fn set_light_theme(ctx: &egui::Context) {
    ctx.set_visuals(egui::Visuals::light());
    let mut style = (*ctx.global_style()).clone();
    style.spacing.item_spacing = egui::vec2(8.0, 5.0);
    style.spacing.slider_width = 150.0;
    ctx.set_global_style(style);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn si_prefixes() {
        assert_eq!(si(1.46e-6, "W"), "1.46 µW");
        assert_eq!(si(0.000994, "W"), "994 µW");
        assert_eq!(si(25.0, "W/m²"), "25.0 W/m²");
        let (k, u) = si_scale(3e-8, "W");
        assert_eq!(u, "nW");
        assert!((k - 1e9).abs() < 1.0);
    }
}
