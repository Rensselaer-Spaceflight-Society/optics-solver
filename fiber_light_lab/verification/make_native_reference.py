"""Write reference values from physics.py for the native (Rust) parity tests.

Run from fiber_light_lab/:  python verification/make_native_reference.py
Output: native/tests/fixtures/python_reference.json
"""
from dataclasses import asdict
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from physics import (Settings, transport, footprint, encircled_radius,  # noqa: E402
                     spatial_coupling, fresnel_transmission, coherence_curve)

CASES = {
    'default': Settings(),
    'band_400_700': Settings(band='400–700 nm'),
    'underfilled_na_0_3': Settings(launch_na=0.3),
    'overfilled_na_0_8': Settings(launch_na=0.8),
    'spot_offset': Settings(spot_radius_mm=0.6, offset_mm=0.3),
    'long_lossy': Settings(length_m=2.0, extra_loss_db=3.0, incident_mw=5.0),
    'short_small_spot': Settings(length_m=0.0, spot_radius_mm=0.05),
}


def transport_record(s):
    t = transport(s)
    wl = t['wavelength_nm']
    pick = [i for i, w in enumerate(wl) if w in (400, 450, 550, 650, 700, 850, 1000)]
    return {
        'settings': asdict(s),
        'output_mw': t['output_mw'],
        'guided_mw': t['guided_mw'],
        'coupling': t['coupling'],
        'spatial_fraction': t['spatial_fraction'],
        'angular_fraction': t['angular_fraction'],
        'face_pair_transmission': t['face_pair_transmission'],
        'half_angle_deg': t['half_angle_deg'],
        'samples': [{'wavelength_nm': float(wl[i]),
                     'incident': float(t['incident_spectrum'][i]),
                     'guided': float(t['guided_spectrum'][i]),
                     'output': float(t['output_spectrum'][i]),
                     'loss_db_per_m': float(t['loss_db_per_m'][i])} for i in pick],
    }


def footprint_record(distance, na, n=1025):
    s = Settings(distance_mm=distance, launch_na=na)
    axis, e, rmax = footprint(s, n)
    c = len(axis)//2
    return {'distance_mm': distance, 'launch_na': na, 'grid': n,
            'output_mw': transport(s)['output_mw'],
            'peak_mw_per_mm2': float(e[c, c]),
            'encircled90_radius_mm': encircled_radius(axis, e, 0.9),
            'encircled50_radius_mm': encircled_radius(axis, e, 0.5),
            'envelope_radius_mm': float(rmax),
            'profile_mm': [float(v) for v in axis[c::max(1, (n//2)//16)]],
            'profile_mw_per_mm2': [float(v) for v in e[c, c::max(1, (n//2)//16)]]}


def main():
    opd, gamma = coherence_curve(Settings())
    record = {
        'generator': 'verification/make_native_reference.py',
        'numpy': np.__version__,
        'transport': {name: transport_record(s) for name, s in CASES.items()},
        'spatial_coupling': [{'radius_mm': w, 'offset_mm': d, 'fraction': spatial_coupling(w, d)}
                             for w, d in [(0.4, 0), (0.4, 0.3), (0.05, 0), (0.6, 0.3), (2.0, 0.0),
                                          (2.0, 1.5), (0.1, 0.45), (0.05, 1.5), (1.0, 0.7)]],
        'fresnel': [{'sin_external': s, 'transmission': float(fresnel_transmission(s))}
                    for s in (0.0, 0.1, 0.3, 0.5, 0.7, 0.9)],
        'footprint': [footprint_record(z, na) for z, na in
                      [(0.0, 0.5), (5.0, 0.5), (8.0, 0.5), (63.5, 0.5), (150.0, 0.5), (63.5, 0.2)]],
        'coherence': {'opd_um': [float(opd[i]) for i in (0, 20, 40, 100, 300)],
                      'gamma': [float(gamma[i]) for i in (0, 20, 40, 100, 300)]},
    }
    out = ROOT/'native'/'tests'/'fixtures'/'python_reference.json'
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(record, indent=2), encoding='utf-8')
    print('Wrote', out)


if __name__ == '__main__':
    main()
