//! Fiber Light Lab: native simulator for the regolith-epoxy goniometer bench.
//!
//! Light path: MI-152 EKE halogen -> GH4001 POF -> Thorlabs LB1092 lens ->
//! sample at the goniometer rotation center -> detector on a moving arm.
//! `physics` holds every calculation; the GUI lives in `src/app`.

pub mod physics;
