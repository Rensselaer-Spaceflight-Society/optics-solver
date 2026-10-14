//! Headless runner: `fiberlab-cli --config settings.json --out results/`.
//! Uses exactly the same physics as the desktop app.

use fiber_light_lab::physics::report;
use fiber_light_lab::physics::sample::Simulant;
use fiber_light_lab::physics::sim::{self, Quality, SimSettings};
use std::path::PathBuf;
use std::process::ExitCode;

const HELP: &str = "fiberlab-cli — Fiber Light Lab headless simulation

USAGE:
  fiberlab-cli [--config FILE] [--out DIR] [--quality draft|standard|high]
               [--incidence DEG] [--study] [--print-default-config]

  --config FILE           settings JSON (as saved by the desktop app)
  --out DIR               write CSV/JSON results into DIR
  --quality Q             override ray counts
  --incidence DEG         override the illumination angle
  --study                 also run the experiment matrix (0/30/45/60 wt%,
                          as-cast/abraded, 30°/60°, LHS-1 and LMS-1)
  --print-default-config  print the default settings JSON and exit
";

fn main() -> ExitCode {
    match run() {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("error: {e}");
            ExitCode::FAILURE
        }
    }
}

fn run() -> Result<(), String> {
    let mut args = std::env::args().skip(1);
    let mut settings = SimSettings::default();
    let mut out: Option<PathBuf> = None;
    let mut study = false;
    while let Some(a) = args.next() {
        match a.as_str() {
            "-h" | "--help" => {
                print!("{HELP}");
                return Ok(());
            }
            "--print-default-config" => {
                println!("{}", serde_json::to_string_pretty(&SimSettings::default()).unwrap());
                return Ok(());
            }
            "--config" => {
                let p = args.next().ok_or("--config needs a file")?;
                let text = std::fs::read_to_string(&p).map_err(|e| format!("{p}: {e}"))?;
                settings = SimSettings::from_json(&text)?;
            }
            "--out" => out = Some(PathBuf::from(args.next().ok_or("--out needs a directory")?)),
            "--quality" => {
                settings.quality = match args.next().as_deref() {
                    Some("draft") => Quality::Draft,
                    Some("standard") => Quality::Standard,
                    Some("high") => Quality::High,
                    _ => return Err("--quality must be draft, standard or high".into()),
                }
            }
            "--incidence" => {
                settings.bench.incidence_deg =
                    args.next().and_then(|v| v.parse().ok()).ok_or("--incidence needs a number")?
            }
            "--study" => study = true,
            other => return Err(format!("unknown argument '{other}' (see --help)")),
        }
    }
    let r = sim::simulate(&settings, None)?;
    let summary = report::summary_json(&r);
    println!("{}", serde_json::to_string_pretty(&summary).unwrap());
    let rows = if study {
        let rows = sim::run_study(
            &settings,
            &[Simulant::Lhs1, Simulant::Lms1],
            &[0.0, 30.0, 45.0, 60.0],
            &[30.0, 60.0],
            &|_| {},
        )?;
        eprintln!("{}", report::study_csv(&rows));
        Some(rows)
    } else {
        None
    };
    if let Some(dir) = out {
        let files = report::write_folder(&dir, &r, rows.as_deref()).map_err(|e| e.to_string())?;
        eprintln!("wrote {} files to {}", files.len(), dir.display());
    }
    Ok(())
}
