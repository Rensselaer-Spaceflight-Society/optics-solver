//! Reference curves compiled into the executable from `fiber_light_lab/data/`.
//! The CSV files stay the single editable source shared with the Python app.

use super::math::Curve;
use std::sync::OnceLock;

pub const EKE_SPECTRUM_CSV: &str = include_str!("../../../data/eke_relative_spectrum.csv");
pub const GH4001_LOSS_CSV: &str = include_str!("../../../data/gh4001_loss.csv");
pub const CIE_V_LAMBDA_CSV: &str = include_str!("../../../data/cie_photopic_vlambda.csv");
pub const SI_RESPONSIVITY_CSV: &str = include_str!("../../../data/si_photodiode_responsivity_typical.csv");

fn cached(cell: &'static OnceLock<Curve>, csv: &str) -> &'static Curve {
    cell.get_or_init(|| Curve::from_csv(csv).expect("embedded reference CSV is valid"))
}

/// Relative spectral radiance of the EKE lamp with dichroic reflector
/// (manual approximation of the supplied chart, 350-2100 nm).
pub fn eke_spectrum() -> &'static Curve {
    static C: OnceLock<Curve> = OnceLock::new();
    cached(&C, EKE_SPECTRUM_CSV)
}

/// Typical GH4001 transmission loss in dB/m (launch NA 0.1, ESKA chart).
pub fn gh4001_loss() -> &'static Curve {
    static C: OnceLock<Curve> = OnceLock::new();
    cached(&C, GH4001_LOSS_CSV)
}

/// CIE 1924 photopic luminous efficiency V(lambda).
pub fn cie_v_lambda() -> &'static Curve {
    static C: OnceLock<Curve> = OnceLock::new();
    cached(&C, CIE_V_LAMBDA_CSV)
}

/// Typical silicon PIN photodiode responsivity (A/W). Approximate shape only.
pub fn si_responsivity() -> &'static Curve {
    static C: OnceLock<Curve> = OnceLock::new();
    cached(&C, SI_RESPONSIVITY_CSV)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn embedded_curves_parse() {
        assert_eq!(eke_spectrum().x.len(), 25);
        assert_eq!(gh4001_loss().x.len(), 41);
        assert!((cie_v_lambda().eval(555.0) - 0.99975).abs() < 0.01);
        assert!(si_responsivity().eval(650.0) > 0.3);
    }
}
