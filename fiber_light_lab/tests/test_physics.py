import io
import json
import zipfile
from dataclasses import replace
import numpy as np
import pytest
from scipy.stats import qmc
from physics import (Settings, CORE_RADIUS_MM, spatial_coupling, fresnel_transmission,
                     transport, footprint, interference, coherence_curve, mode_count)


def test_defaults_and_analytic_launch():
    s = Settings()
    assert s.length_m == pytest.approx(10*.0254)
    assert s.distance_mm == pytest.approx(2.5*25.4)
    assert spatial_coupling(.4, 0) == pytest.approx(1-np.exp(-2*CORE_RADIUS_MM**2/.4**2))
    assert transport(s)['half_angle_deg'] == pytest.approx(30)
    assert fresnel_transmission(0) == pytest.approx(1-((1.49-1)/(1.49+1))**2)
    assert mode_count()[1] > 2e6


@pytest.mark.parametrize('distance,na', [(0,.5), (.1,.5), (5,.1), (63.5,.5), (300,.5)])
def test_plane_conserves_spectral_power(distance, na):
    s = Settings(distance_mm=distance, launch_na=na)
    x, e, radius = footprint(s)
    dx = x[1]-x[0]
    assert e.min() >= 0
    assert e.sum()*dx*dx == pytest.approx(transport(s)['output_mw'], rel=1e-11)
    assert radius == pytest.approx(.49+distance*np.tan(np.arcsin(na)))
    assert np.allclose(e, e[::-1, :], atol=1e-12)
    assert np.allclose(e, e[:, ::-1], atol=1e-12)


def test_transport_limits():
    s = Settings(length_m=0, spot_radius_mm=.05)
    t = transport(s)
    assert 0 < t['output_mw'] < t['guided_mw'] < s.incident_mw
    assert np.allclose(transport(replace(s, extra_loss_db=10))['output_spectrum'], t['output_spectrum']/10)
    assert np.allclose(transport(replace(s, incident_mw=2))['output_spectrum'], t['output_spectrum']*2)
    assert transport(replace(s, length_m=10))['output_mw'] < t['output_mw']
    assert transport(replace(s, offset_mm=1))['output_mw'] < t['output_mw']*.01
    assert transport(replace(s, launch_na=.8))['angular_fraction'] == pytest.approx((.5/.8)**2)
    at650 = np.where(t['wavelength_nm']==650)[0][0]
    ratio = transport(replace(s, length_m=10))['output_spectrum'][at650]/t['output_spectrum'][at650]
    assert ratio == pytest.approx(10**(-.16*10/10))


def test_far_field_inverse_square_and_grid_convergence():
    s = Settings(distance_mm=100)
    _, e1, _ = footprint(s)
    _, e2, _ = footprint(replace(s, distance_mm=200))
    assert e1.max()/e2.max() == pytest.approx(4, rel=.03)
    _, fine, _ = footprint(Settings(), 513)
    _, coarse, _ = footprint(Settings(), 257)
    assert coarse.max() == pytest.approx(fine.max(), rel=.015)


def test_independent_ray_trace_second_moment():
    # Independent quasi-Monte-Carlo ray transport checks the convolution Jacobian.
    s = Settings(distance_mm=8)
    u = qmc.Sobol(4, scramble=True, seed=7).random_base2(18)
    r = CORE_RADIUS_MM*np.sqrt(u[:,0])
    phi = 2*np.pi*u[:,1]
    sinth = .5*np.sqrt(u[:,2])
    theta = 2*np.pi*u[:,3]
    disp = s.distance_mm*sinth/np.sqrt(1-sinth**2)
    rx = r*np.cos(phi)+disp*np.cos(theta)
    ry = r*np.sin(phi)+disp*np.sin(theta)
    weight = fresnel_transmission(sinth)**2
    ray_moment = np.average(rx*rx+ry*ry, weights=weight)
    ax, e, _ = footprint(s, 513)
    x,y = np.meshgrid(ax,ax)
    raster_moment = np.sum(e*(x*x+y*y))/e.sum()
    assert raster_moment == pytest.approx(ray_moment, rel=.008)


def test_interference_and_ensemble_limits():
    s = Settings()
    _, halogen, cross = interference(s, mode='halogen')
    assert np.all(halogen==1) and np.all(cross==0)
    axis, coherent, cross = interference(s)
    assert coherent.min() >= 0
    assert cross.min() < 0 < cross.max()
    assert coherent.std()/coherent.mean() > .7
    _, averaged, _ = interference(s, mode='ensemble', averages=64)
    assert averaged.std()/averaged.mean() < .22
    _, repeat, _ = interference(s)
    np.testing.assert_array_equal(coherent, repeat)
    _, changed, _ = interference(s, phase_shift=1000)
    assert not np.allclose(coherent, changed)
    # Maximum intensity spatial frequency = 2*spread/lambda.
    spread = min(.5, .49/np.hypot(.49,s.distance_mm))
    assert axis[1]-axis[0] < .65/(4*spread)


def test_temporal_coherence():
    _, gamma = coherence_curve(Settings())
    assert gamma[0] == pytest.approx(1)
    assert gamma.max() <= 1+1e-12
    assert gamma[-1] < .1


def test_zero_power_and_invalid_settings():
    s = Settings(incident_mw=0)
    assert transport(s)['output_mw'] == 0
    assert not np.any(footprint(s)[1])
    assert not np.any(coherence_curve(s)[1])
    for values in [dict(launch_na=0), dict(length_m=-1), dict(incident_mw=float('nan'))]:
        with pytest.raises(ValueError):
            Settings(**values)


def test_app_callback_export_and_local_routes():
    from app import app, update, export_bundle
    outputs = update(10, 2.5, 1, .5, .4, 0, 0, '400–1000 nm', [], 'beam')
    assert len(outputs)==15
    assert outputs[-1]['length_m'] == .254
    client = app.server.test_client()
    for path in ['/', '/_dash-layout', '/_dash-dependencies', '/assets/style.css']:
        assert client.get(path).status_code == 200
    archive = zipfile.ZipFile(io.BytesIO(export_bundle(Settings())))
    settings = json.loads(archive.read('settings.json'))
    assert settings['output_mw'] > 0
    csv = np.loadtxt(io.BytesIO(archive.read('irradiance.csv')), delimiter=',', skiprows=1)
    assert np.allclose(csv[:,3], csv[:,2]*1000)
    html = archive.read('output_3d.html').decode()
    assert 'plotly.js' in html and 'Plotly.newPlot' in html
