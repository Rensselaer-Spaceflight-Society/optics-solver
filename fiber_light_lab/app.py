"""Run `python app.py`, then open http://127.0.0.1:8057."""
import argparse
from dataclasses import asdict
import io
import json
from pathlib import Path
import time
import zipfile
import numpy as np
from dash import Dash, Input, Output, State, dcc, html, no_update
import plotly.io as pio
from physics import Settings, transport, footprint, encircled_radius, interference, mode_count
import figures as f

ROOT = Path(__file__).parent
app = Dash(__name__, title='Fiber Light Lab', update_title=None)
server = app.server


def graph(id, **kwargs):
    return dcc.Graph(id=id, config=f.CONFIG, **kwargs)


def card(title, subtitle, children, className='panel'):
    return html.Section([html.Div([html.H3(title), html.P(subtitle)], className='panel-heading'), children], className=className)


def number(label, id, value, unit, min=0, max=10000, step=.1):
    return html.Div([html.Label(label, htmlFor=id), html.Div([
        dcc.Input(id=id, type='number', value=value, min=min, max=max, step=step, debounce=True),
        html.Span(unit)], className='number-row')], className='control')


def slider(label, id, value, min, max, step, marks):
    return html.Div([html.Label(label), dcc.Slider(id=id, min=min, max=max, step=step, value=value,
        marks=marks, tooltip={'placement': 'bottom', 'always_visible': False}, updatemode='mouseup')], className='control slider')


def metric(label, id):
    return html.Div([html.Span(label), html.Strong(id=id)], className='metric')


app.layout = html.Div([
    dcc.Store(id='settings'), dcc.Download(id='download'),
    html.Aside([
        html.Div([html.Div('◌', className='logo'), html.Div(['FIBER', html.Br(), 'LIGHT LAB'], className='wordmark')], className='brand'),
        html.P('LUNAR MATERIALS / OPTICAL BENCH', className='eyebrow sidebar-caption'),
        html.H2('Your setup'),
        html.P('MI-152 · EKE halogen\nGH4001 · PMMA step index', className='equipment'),
        number('Cable length', 'length', 10, 'in', max=400, step=.25),
        number('Tip → sample', 'distance', 2.5, 'in', max=12, step=.05),
        number('Incident optical power', 'power', 1, 'mW', max=1000, step=.1),
        html.P('Assumed power in the selected band, before coupling. The 150 W lamp rating is electrical.', className='hint'),
        html.Label('Wavelength band'),
        dcc.Dropdown(id='band', options=['400–1000 nm', '400–700 nm'], value='400–1000 nm', clearable=False),
        html.Div(className='rule'),
        slider('Launch NA', 'na', .5, .05, .9, .01, {.05:'.05', .5:'.50', .9:'.90'}),
        slider('Focused spot · 1/e² radius (mm)', 'spot', .4, .05, 2, .01, {.05:'.05', 1:'1', 2:'2'}),
        slider('Lateral misalignment (mm)', 'offset', 0, 0, 1.5, .01, {0:'0', .75:'.75', 1.5:'1.5'}),
        slider('Additional loss (dB)', 'loss', 0, 0, 10, .1, {0:'0', 5:'5', 10:'10'}),
        html.Button('Restore your setup', id='reset', className='secondary'),
        html.Button('Export this simulation ↓', id='export', className='primary'),
        html.P('CSV data · settings · offline 3D view', className='hint'),
        html.Div(['980 µm core', html.Br(), 'NA 0.50  /  n = 1.49'], className='spec-stamp'),
    ], className='sidebar'),
    html.Main([
        html.Header([html.Div([html.P('SOURCE → FIBER → FREE SPACE', className='eyebrow'),
                              html.H1('Light, from launch to landing.'),
                              html.P('A spectral transport model for your regolith illumination bench.', className='subtitle')]),
                     html.Div([html.Span(className='status-dot'), 'LOCAL SIMULATION'], className='status')]),
        html.Div([metric('OUTPUT POWER · mW', 'pout'), metric('90% POWER DIAMETER · mm', 'diameter'),
                  metric('PEAK IRRADIANCE · W/m²', 'peak'), metric('OUTPUT HALF-ANGLE', 'angle')], className='metrics'),
        html.Div(id='summary', className='summary'),
        dcc.Tabs(id='tabs', value='beam', className='tabs', children=[
            dcc.Tab(label='01  Output field', value='beam', className='tab', selected_className='tab selected', children=[
                html.Div([card('The output volume', 'Drag to orbit · scroll to zoom · all geometry in millimeters',
                    html.Div([html.Div([dcc.RadioItems(id='scene-view', options=[{'label':'Full beam', 'value':'beam'}, {'label':'Tip detail', 'value':'tip'}],
                         value='beam', inline=True), dcc.Checklist(id='log', options=[{'label':'Log intensity', 'value':'log'}], value=[])], className='toolbar'),
                         dcc.Loading(graph('scene'), type='dot', color=f.AMBER),
                         html.P('Rays and the translucent envelope indicate propagation; clear air itself does not glow. The sample plane shows incident irradiance.', className='plot-note')])),
                    card('At the sample', 'A perpendicular observation plane in air', graph('footprint'))], className='beam-grid'),
                html.Div([card('Across the footprint', 'Central line profile · absolute values scale with assumed optical power', graph('profile')),
                          card('The spectrum in transit', 'EKE curve approximation × coupling × GH4001 transmission', graph('spectrum'))], className='two-grid'),
            ]),
            dcc.Tab(label='02  Inside the fiber', value='fiber', className='tab', selected_className='tab selected', children=[
                html.Div([card('Total internal reflection', 'Physical tip segment · representative meridional rays, not individual modes', graph('fiber-scene')),
                          card('Power along the cable', 'Bulk spectral loss; additional loss and the exit face are applied afterward', graph('length-power'))], className='two-grid'),
                html.Div([card('GH4001 attenuation', 'Approximate points from the manufacturer’s typical curve (launch NA 0.1)', graph('loss-plot', figure=f.loss_plot())),
                          card('Launch & acceptance', 'The light source optics control which rays enter the fiber.', html.Div(id='coupling-info', className='explanation'))], className='two-grid'),
            ]),
            dcc.Tab(label='03  Interference', value='wave', className='tab', selected_className='tab selected', children=[
                html.Div([html.H3('Look closer. Coherence changes the picture.'),
                    html.P('A local, on-axis wave patch illustrates interference. It is a 48-component scalar surrogate, not a measured GH4001 speckle pattern or a full fiber mode solution.')], className='intro'),
                html.Div([
                    html.Div([html.Label('Illumination model'), dcc.Dropdown(id='wave-mode', options=[
                        {'label':'Halogen · incoherent mean', 'value':'halogen'},
                        {'label':'Coherent monochromatic · comparison', 'value':'coherent'},
                        {'label':'Independent patterns · averaging demo', 'value':'ensemble'}], value='halogen', clearable=False)], className='wave-select'),
                    html.Div([html.Label('Display'), dcc.RadioItems(id='wave-view', options=[{'label':'Intensity', 'value':'intensity'},
                        {'label':'Constructive / destructive', 'value':'cross'}], value='intensity')]),
                    html.Div([html.Label('Wavelength · nm'), dcc.Slider(id='wavelength', min=400, max=1000, value=650, step=5,
                        marks={400:'400', 650:'650', 1000:'1000'}, updatemode='mouseup')], className='wave-slider'),
                ], className='wave-controls'),
                html.Div([html.Div([html.Label('Axial phase exploration · µm'), dcc.Slider(id='phase', min=0, max=1000, value=0, step=1,
                         marks={0:'0', 500:'500', 1000:'1000'}, updatemode='mouseup')]),
                    html.Div([html.Label('Independent realizations'), dcc.RadioItems(id='averages', options=[1, 4, 16, 64], value=16, inline=True)]),
                    html.Button('New phase realization', id='seed', className='secondary')], className='wave-controls'),
                html.Div(id='wave-summary', className='summary'),
                html.Div([card('Resolved wave patch', 'Coordinates adapt to wavelength and the aperture seen from the sample.', dcc.Loading(graph('wave-map'), color=f.AMBER)),
                          card('Intensity landscape', 'Height is intensity, not physical surface height.', graph('wave-surface'))], className='two-grid'),
                html.Div([card('Temporal coherence of the output spectrum', 'Visibility ceiling for two equally polarized copies of the same spatial field.', graph('coherence')),
                          card('Reading the interference', 'Blue removes intensity; red adds it relative to incoherent addition.', html.Div([
                              html.P('For coherent light: I = |Σ Eⱼ|². The cross terms create bright and dark structure.'),
                              html.P('For the ideal incoherent limit: ⟨I⟩ = Σ ⟨|Eⱼ|²⟩. Cross terms vanish. The default halogen view is this ensemble mean.'),
                              html.P('A narrow filter alone does not make an extended halogen source spatially coherent. The averaging control is an educational ensemble count, not a calibrated exposure time.'),
                              html.P('The phase slider propagates the local plane-wave basis relative to its reference plane; it holds the beam envelope and geometry fixed.')], className='explanation'))], className='two-grid'),
            ]),
            dcc.Tab(label='04  Model & provenance', value='model', className='tab', selected_className='tab selected', children=[
                html.Div([
                    html.H2('An inspectable model, with explicit assumptions.'),
                    html.H3('What comes from the specifications'),
                    html.P('GH-4001-P: 980 µm nominal core, 1000 µm cladding diameter, 2.20 mm jacket, PMMA core index 1.49, step index and NA 0.50. Maximum specified loss is 0.17 dB/m at 650 nm under the stated room-temperature test conditions.'),
                    html.P('The MI-152 uses a 150 W / 21 V EKE lamp. The uploaded spectrum is relative radiance, not calibrated optical power. It includes the dichroic reflector response; a pure blackbody curve would miss that shape.'),
                    html.H3('What is assumed'),
                    html.P('A straight, clean-ended fiber; air outside; no downstream lens; centered normal observation plane; unpolarized light; Gaussian input spot; uniform radiance within the launch cone; a uniform, spatially mixed exit face; conserved angular distribution. Wavelength-independent NA and refractive index. Offset changes accepted power, not the assumed mixed exit profile.'),
                    html.P('The default 1 mW is an illustrative optical input in the selected band, not a conversion of 150 W electrical power. Changing the band holds that entered band power fixed and renormalizes the source; it does not simulate insertion of a filter.'),
                    html.H3('What is approximate'),
                    html.P('Both spectral curves use manually estimated anchor points. GH4001 loss was measured at launch NA 0.1; applying it to other launch angles is an approximation. Extra loss is a user-defined lumped term, not a bend-radius model. No spectral extrapolation is made above 1000 nm or below 400 nm; power outside the chosen band is excluded.'),
                    html.P('The irradiance model is geometric radiometry, appropriate for this large-core, low-coherence source. It does not resolve aperture diffraction or individual fiber modes. The wave explorer illustrates local interference with a sampled scalar angular spectrum; it cannot predict the exact speckle of your cable. A full-field high-NA Maxwell simulation across a 980 µm core and inches of fiber is intentionally not attempted.'),
                    html.H3('About the Desmos disk'),
                    html.P('The supplied square-root surface describes a geometric cap, not an optical irradiance law. Intensity requires optical power, an angular distribution and propagation geometry. Interference requires adding complex fields before squaring. This simulation computes those separately.'),
                    html.H3('For the regolith experiment'),
                    html.P('This version ends at the incident sample plane. It does not yet model regolith absorption, multiple scattering, BRDF, surface roughness, composite inclusions or heating. Exported irradiance and spectra provide boundary conditions for that next stage.'),
                    html.H3('Sources and documentation'),
                    html.Ul([html.Li(html.A('Dolan-Jenner MI-152 product specification', href='https://dolan-jenner.com/products/fiber-lite-mi-152-mi-157-illuminator-series-150w-halogen-light-source-mi152-mi157', target='_blank')),
                             html.Li(html.A('Industrial Fiber Optics: typical ESKA loss spectra', href='https://i-fiberoptics.com/pdf/attenuation.pdf', target='_blank')),
                             html.Li('User-supplied Mitsubishi GH-4001-P datasheet, 14 November 2019, and EKE spectrum image (copied into data/).'),
                             html.Li('MODEL.md contains equations, units, approximation boundaries and the full source list.')]),
                ], className='model-text'),
            ]),
        ]),
        html.Footer([html.Span('FIBER LIGHT LAB / RESEARCH PROTOTYPE'), html.Span(id='timing')]),
    ], className='main'),
], className='shell')


@app.callback([Output(k, 'value') for k in ['length', 'distance', 'power', 'na', 'spot', 'offset', 'loss', 'band']],
              Input('reset', 'n_clicks'), prevent_initial_call=True)
def reset(_):
    return 10, 2.5, 1, .5, .4, 0, 0, '400–1000 nm'


@app.callback(
    [Output(k, 'figure') for k in ['scene', 'footprint', 'profile', 'spectrum', 'fiber-scene', 'length-power', 'coherence']]+
    [Output(k, 'children') for k in ['pout', 'diameter', 'peak', 'angle', 'summary', 'coupling-info', 'timing']]+
    [Output('settings', 'data')],
    [Input(k, 'value') for k in ['length', 'distance', 'power', 'na', 'spot', 'offset', 'loss', 'band', 'log', 'scene-view']])
def update(length, distance, power, na, spot, offset, loss, band, log, scene_view):
    if any(v is None for v in [length, distance, power, na, spot, offset, loss]):
        return [no_update]*15
    start = time.perf_counter()
    try:
        s = Settings(length*.0254, distance*25.4, power, na, spot, offset, loss, band)
    except ValueError:
        return [no_update]*15
    t = transport(s)
    axis, irr, bound = footprint(s)
    r90 = encircled_radius(axis, irr)
    stats = [f"{t['output_mw']:.4g}", f'{2*r90:.2f}', f'{irr.max()*1000:.4g}', f"{t['half_angle_deg']:.1f}°"]
    summary = f'{s.length_m*1000:g} mm cable · {s.distance_mm:g} mm working distance · {2*bound:.1f} mm geometric envelope · {band} · assumed input power'
    _, modes = mode_count()
    info = [html.P(f"Spatial overlap: {t['spatial_fraction']:.1%}. Angular acceptance: {t['angular_fraction']:.1%}."),
            html.P(f"Guided power after the entrance: {t['guided_mw']:.4g} mW. Output after the exit: {t['output_mw']:.4g} mW."),
            html.P(f'High-V mode-count estimate at 650 nm: {modes/1e6:.2f} million, counting both polarizations. This is capacity, not the number actually excited.'),
            html.P('The cutaway shows only the final 8 mm. Reflections use Snell’s law and the nominal core radius; representative ray count is for clarity.'),
            html.P('The exit face is assumed spatially mixed. Actual launch optics, cable bends and polish can change the near-field and far-field profiles.')]
    figures = [f.scene(s, bool(log), scene_view=='tip'), f.irradiance_map(s, bool(log)), f.radial_profile(s),
               f.spectrum(s), f.scene(s, False, True), f.length_power(s), f.coherence_plot(s)]
    return figures+stats+[summary, info, f'Calculated in {(time.perf_counter()-start)*1000:.0f} ms', asdict(s)]


@app.callback(Output('wave-map', 'figure'), Output('wave-surface', 'figure'), Output('wave-summary', 'children'),
              Input('settings', 'data'), Input('wave-mode', 'value'), Input('wavelength', 'value'),
              Input('averages', 'value'), Input('phase', 'value'), Input('seed', 'n_clicks'), Input('wave-view', 'value'))
def update_wave(settings, mode, wavelength, averages, phase, seed, view):
    s = Settings(**settings) if settings else Settings()
    axis, intensity, cross = interference(s, wavelength, mode, averages, phase, (seed or 0)+4)
    figs = f.wave_figures(axis, intensity, cross, view)
    contrast = float(intensity.std()/intensity.mean())
    text = f'{axis[-1]-axis[0]:.1f} µm patch · {wavelength} nm basis · contrast σ/mean = {contrast:.3f} · '
    text += 'ideal incoherent mean; zero cross terms' if mode=='halogen' else 'illustrative local field, normalized to incoherent intensity'
    return *figs, text


def csv_bytes(columns, header):
    out = io.StringIO()
    np.savetxt(out, columns, delimiter=',', header=header, comments='')
    return out.getvalue()


def export_bundle(s, mode='halogen', wavelength=650, averages=16, phase=0, seed=4):
    t = transport(s)
    axis, irr, _ = footprint(s)
    x, y = np.meshgrid(axis, axis)
    wa, wi, wc = interference(s, wavelength, mode, averages, phase, seed)
    wx, wy = np.meshgrid(wa, wa)
    buffer = io.BytesIO()
    metadata = dict(settings=asdict(s), output_mw=t['output_mw'],
                    wave=dict(mode=mode, wavelength_nm=wavelength, averages=averages, phase_um=phase, seed=seed),
                    assumptions='See MODEL.md. Input optical power is assumed. Spectrum and loss curves are approximate. Wave patch is an illustrative surrogate.')
    with zipfile.ZipFile(buffer, 'w', zipfile.ZIP_DEFLATED) as z:
        z.writestr('settings.json', json.dumps(metadata, indent=2))
        z.writestr('irradiance.csv', csv_bytes(np.c_[x.ravel(), y.ravel(), irr.ravel(), irr.ravel()*1000], 'x_mm,y_mm,irradiance_mW_per_mm2,irradiance_W_per_m2'))
        z.writestr('spectrum.csv', csv_bytes(np.c_[t['wavelength_nm'], t['incident_spectrum'], t['guided_spectrum'], t['output_spectrum'], t['loss_db_per_m']], 'wavelength_nm,input_mW_per_nm,guided_mW_per_nm,output_mW_per_nm,loss_dB_per_m'))
        z.writestr('wave_patch_illustrative.csv', csv_bytes(np.c_[wx.ravel(), wy.ravel(), wi.ravel(), wc.ravel()], 'x_um,y_um,I_over_incoherent,cross_term_over_incoherent'))
        z.writestr('output_3d.html', pio.to_html(f.scene(s), include_plotlyjs=True, full_html=True, config=f.CONFIG))
        z.writestr('MODEL.md', (ROOT/'MODEL.md').read_text(encoding='utf-8'))
        for name in ['eke_relative_spectrum.csv', 'gh4001_loss.csv']:
            z.write(ROOT/'data'/name, 'reference_curves/'+name)
    return buffer.getvalue()


@app.callback(Output('download', 'data'), Input('export', 'n_clicks'), State('settings', 'data'),
              State('wave-mode', 'value'), State('wavelength', 'value'), State('averages', 'value'),
              State('phase', 'value'), State('seed', 'n_clicks'), prevent_initial_call=True)
def download(_, settings, mode, wavelength, averages, phase, seed):
    s = Settings(**settings) if settings else Settings()
    return dcc.send_bytes(export_bundle(s, mode, wavelength, averages, phase, (seed or 0)+4), 'fiber_light_simulation.zip')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--port', type=int, default=8057)
    args = parser.parse_args()
    app.run(host='127.0.0.1', port=args.port, debug=False)
