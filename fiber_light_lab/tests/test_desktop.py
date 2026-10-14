"""Checks the desktop rendering and export without starting a browser/server."""
import io
import json
import zipfile
import numpy as np
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
from fiber_simulation import calculate, draw_output, draw_transport, draw_interference, save_data
from physics import Settings, transport
from dataclasses import replace
import physics


OPTIONS = dict(wavelength_nm=650, mode='coherent', averages=16, phase_shift=0, seed=4)


def test_cached_results_follow_changes_and_keep_physical_units():
    s = Settings()
    original = calculate(s, 2, OPTIONS)
    doubled = calculate(replace(s, incident_mw=2), 2, OPTIONS)
    np.testing.assert_allclose(doubled['footprint'][1], 2*original['footprint'][1])
    np.testing.assert_array_equal(doubled['wave'][1], original['wave'][1])
    for change in [dict(length_m=.5), dict(distance_mm=2), dict(launch_na=.2)]:
        changed = replace(s, **change)
        result = calculate(changed, 2, OPTIONS)
        for cached, direct in zip(result['wave'], physics.interference(changed, **OPTIONS)):
            np.testing.assert_array_equal(cached, direct)
    for change in [dict(phase_shift=900), dict(wavelength_nm=500), dict(seed=5), dict(mode='ensemble', averages=4)]:
        options = {**OPTIONS, **change}
        result = calculate(s, 2, options)
        np.testing.assert_array_equal(result['wave'][1], physics.interference(s, **options)[1])
        assert not np.allclose(result['wave'][1], original['wave'][1])
    assert calculate(replace(s, incident_mw=0), 0, OPTIONS)['diameter90'] == 0


def test_incoherent_mean_skips_complex_basis(monkeypatch):
    expected_axis = physics.wave_basis(650., 63.5, .5, 4)[0].copy()
    def forbidden(*args, **kwargs):
        raise AssertionError('Halogen mean must not allocate a complex wave basis')
    monkeypatch.setattr(physics, 'wave_basis', forbidden)
    axis, intensity, cross = physics.interference(Settings(), mode='halogen')
    np.testing.assert_array_equal(axis, expected_axis)
    assert np.all(intensity == 1) and np.all(cross == 0)


def test_desktop_preserves_physics_and_exports(tmp_path):
    s = Settings()
    result = calculate(s, 2, OPTIONS)
    assert result['transport']['output_mw'] == transport(s)['output_mw']
    path = tmp_path/'desktop.zip'
    save_data(path, s, OPTIONS)
    with zipfile.ZipFile(path) as z:
        assert z.testzip() is None
        config = json.loads(z.read('settings.json'))
        assert config['settings']['distance_mm'] == 63.5
        assert config['wave'] == OPTIONS
        values = np.loadtxt(io.BytesIO(z.read('irradiance.csv')), delimiter=',', skiprows=1)
        assert np.isclose(values[:,2].max(), result['footprint'][1].max()*1000)


def test_all_desktop_plots_render_and_zero_power_is_finite():
    for s in [Settings(), Settings(incident_mw=0, distance_mm=0)]:
        f = Figure(figsize=(10, 6), layout='constrained')
        canvas = FigureCanvasAgg(f)
        result = calculate(s, 0, OPTIONS)
        for tip, log in [(False, False), (True, True)]:
            draw_output(f, s, result, tip, log)
            canvas.draw()
        draw_transport(f, s, calculate(s, 1, OPTIONS))
        canvas.draw()
        for mode, cross in [('coherent', False), ('coherent', True), ('halogen', True)]:
            result = calculate(s, 2, {**OPTIONS, 'mode':mode})
            draw_interference(f, result, cross)
            canvas.draw()
            assert np.isfinite(result['wave'][1]).all()
