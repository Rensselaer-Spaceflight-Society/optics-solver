"""Small Plotly scenes; all physical calculations live in physics.py."""
import numpy as np
import plotly.graph_objects as go
from physics import (CORE_RADIUS_MM, CORE_INDEX, footprint, transport, reference_data,
                     normalized_footprint, coherence_curve, fresnel_transmission)

BG = '#111a21'
TEXT = '#d8e2e5'
MUTED = '#8eaaaF'
AMBER = '#e4b66f'
TEAL = '#71c5c4'
HEAT = [[0, '#111a21'], [.12, '#263b46'], [.38, '#4c7275'], [.65, '#b99862'], [1, '#ffe0a3']]
CONFIG = dict(displaylogo=False, scrollZoom=True,
              modeBarButtonsToRemove=['lasso2d', 'select2d'],
              toImageButtonOptions=dict(format='png', scale=2))


def style(fig, height=320, **kwargs):
    fig.update_layout(template='plotly_dark', paper_bgcolor=BG, plot_bgcolor=BG,
                      font=dict(family='Segoe UI, Arial, sans-serif', color=TEXT, size=12),
                      height=height, margin=dict(l=58, r=25, t=25, b=48),
                      legend=dict(orientation='h', y=1.14, x=0),
                      uirevision='preserve-view', **kwargs)
    fig.update_xaxes(gridcolor='#25323b', zeroline=False)
    fig.update_yaxes(gridcolor='#25323b', zeroline=False)
    return fig


def cylinder(fig, start, stop, radius, color, opacity=1, cutaway=False):
    theta = np.linspace(.4 if cutaway else 0, 5.0 if cutaway else 2*np.pi, 45)
    axial, theta = np.meshgrid([start, stop], theta)
    fig.add_trace(go.Surface(x=axial, y=radius*np.cos(theta), z=radius*np.sin(theta),
                            surfacecolor=np.zeros_like(axial), colorscale=[[0, color], [1, color]],
                            showscale=False, opacity=opacity, hoverinfo='skip',
                            lighting=dict(ambient=.6, diffuse=.6)))


def rays_inside(fig, segment=8, launch_na=.5, count=10):
    # Meridional representative rays; triangle-wave folding is exact for these rays.
    for i in range(count):
        phi = 2*np.pi*i/count
        sint = launch_na*np.sqrt((i+.5)/count)/CORE_INDEX
        slope = np.tan(np.arcsin(sint))
        z = np.linspace(-segment, 0, 250)
        unfolded = slope*z + .12*np.sin(i)
        a = CORE_RADIUS_MM
        r = a-np.abs(np.mod(unfolded+a, 4*a)-2*a)
        fig.add_trace(go.Scatter3d(x=z, y=r*np.cos(phi), z=r*np.sin(phi), mode='lines',
                                  line=dict(color=AMBER, width=2), opacity=.75,
                                  hoverinfo='skip', showlegend=False))


def scene(s, logarithmic=False, detail=False):
    fig = go.Figure()
    length = 8 if detail else min(12, max(4, s.distance_mm*.16))
    cylinder(fig, -length, -.45, 1.10, '#28363e', cutaway=True)
    cylinder(fig, -length, 0, .50, '#729394', .13)
    rays_inside(fig, length, s.output_na, 8)
    axis, irr, bound = footprint(s, 129)
    xx, yy = np.meshgrid(axis, axis)
    peak = max(float(irr.max()), 1e-30)
    display = np.log10(np.maximum(irr/peak, 1e-4)) if logarithmic else irr/peak
    fig.add_trace(go.Surface(x=np.full_like(xx, s.distance_mm), y=xx, z=yy,
                            surfacecolor=display, colorscale=HEAT, cmin=-4 if logarithmic else 0,
                            cmax=0 if logarithmic else 1, showscale=False,
                            customdata=irr*1000,
                            hovertemplate='x %{y:.2f} mm<br>y %{z:.2f} mm<br>%{customdata:.3f} W/m²<extra></extra>',
                            lighting=dict(ambient=1, diffuse=0)))
    theta = np.linspace(0, 2*np.pi, 65)
    zz, tt = np.meshgrid([0, s.distance_mm], theta)
    rr = CORE_RADIUS_MM + zz*np.tan(np.arcsin(s.output_na))
    fig.add_trace(go.Surface(x=zz, y=rr*np.cos(tt), z=rr*np.sin(tt), opacity=.055,
                            colorscale=[[0, AMBER], [1, AMBER]], showscale=False, hoverinfo='skip'))
    for i in range(22):
        phi = i*np.pi*(3-np.sqrt(5))
        sine = s.output_na*np.sqrt((i+.5)/22)
        slope = np.tan(np.arcsin(sine))
        r0 = CORE_RADIUS_MM*.6*np.sin(i*2.3)
        end = r0+s.distance_mm*slope
        fig.add_trace(go.Scatter3d(x=[0, s.distance_mm], y=np.cos(phi)*np.array([r0, end]),
                                  z=np.sin(phi)*np.array([r0, end]), mode='lines',
                                  line=dict(color=AMBER, width=1.5), opacity=.32,
                                  showlegend=False, hoverinfo='skip'))
    fig.add_trace(go.Scatter3d(x=[0, s.distance_mm], y=[0, 0], z=[0, 0],
                              mode='markers+text', marker=dict(size=3, color=TEAL),
                              text=['Fiber tip', 'Observation plane'], textposition='top center',
                              textfont=dict(size=12), showlegend=False, hoverinfo='skip'))
    ax = dict(showbackground=False, gridcolor='#2a3942', zeroline=False,
              color=MUTED, showspikes=False)
    scene_settings = dict(xaxis=dict(**ax, title='z · mm'),
                          yaxis=dict(**ax, title='x · mm'), zaxis=dict(**ax, title='y · mm'),
                          aspectmode='data', camera=dict(eye=dict(x=-1.45, y=1.65, z=1.05)))
    if detail:
        scene_settings.update(xaxis=dict(**ax, title='z · mm', range=[-8, 5]),
                              yaxis=dict(**ax, title='x · mm', range=[-3, 3]),
                              zaxis=dict(**ax, title='y · mm', range=[-3, 3]))
    return style(fig, 500, scene=scene_settings)


def irradiance_map(s, log=False):
    axis, irr, bound = footprint(s)
    irradiance_w = irr*1000  # 1 mW/mm² = 1000 W/m².
    peak = max(float(irr.max()), 1e-30)
    value = np.log10(np.maximum(irr/peak, 1e-4)) if log else irradiance_w
    fig = go.Figure(go.Heatmap(x=axis, y=axis, z=value, colorscale=HEAT,
                              zmin=-4 if log else 0, zmax=0 if log else None,
                              customdata=irradiance_w,
                              colorbar=dict(title='log₁₀(I/peak)' if log else 'W/m²', thickness=11),
                              hovertemplate='x %{x:.2f} mm<br>y %{y:.2f} mm<br>%{customdata:.4f} W/m²<extra></extra>'))
    fig.update_xaxes(title='x · mm', constrain='domain')
    fig.update_yaxes(title='y · mm', scaleanchor='x', scaleratio=1)
    return style(fig, 430)


def radial_profile(s):
    axis, irr, _ = footprint(s)
    fig = go.Figure(go.Scatter(x=axis, y=irr[len(axis)//2]*1000,
                               line=dict(color=AMBER, width=2), fill='tozeroy',
                               fillcolor='rgba(228,182,111,.09)', name='Central section'))
    fig.update_xaxes(title='x through beam center · mm')
    fig.update_yaxes(title='Irradiance · W/m²', rangemode='tozero')
    return style(fig, 260)


def spectrum(s):
    t = transport(s)
    fig = go.Figure()
    for key, name, color in [('incident_spectrum', 'Incident', '#748e9a'),
                             ('guided_spectrum', 'Inside fiber', TEAL),
                             ('output_spectrum', 'Output', AMBER)]:
        fig.add_trace(go.Scatter(x=t['wavelength_nm'], y=t[key]*1000,
                                line=dict(color=color, width=2), name=name))
    fig.update_xaxes(title='Wavelength · nm')
    fig.update_yaxes(title='Spectral power · µW/nm', rangemode='tozero')
    return style(fig, 300)


def length_power(s):
    t = transport(s)
    length = np.linspace(0, max(s.length_m*1.1, .3), 120)
    powers = np.trapezoid(t['guided_spectrum'][:, None]*10**(-t['loss_db_per_m'][:, None]*length/10),
                          t['wavelength_nm'], axis=0)
    fig = go.Figure(go.Scatter(x=length, y=powers, line=dict(color=TEAL, width=2), name='Inside fiber'))
    fig.add_vline(x=s.length_m, line_dash='dot', line_color=AMBER)
    fig.update_xaxes(title='Distance along cable · m')
    fig.update_yaxes(title='Guided power · mW', rangemode='tozero')
    return style(fig, 310)


def loss_plot():
    _, loss = reference_data()
    fig = go.Figure(go.Scatter(x=loss[:, 0], y=loss[:, 1], line=dict(color=TEAL, width=2), name='Typical curve'))
    fig.add_trace(go.Scatter(x=[650], y=[.17], mode='markers', name='650 nm specification maximum',
                            marker=dict(color=AMBER, size=8, symbol='diamond')))
    fig.update_xaxes(title='Wavelength · nm')
    fig.update_yaxes(title='Attenuation · dB/m', type='log', tickvals=[.1, 1, 10, 100], ticktext=['0.1', '1', '10', '100'])
    return style(fig, 310)


def wave_figures(axis, intensity, cross, view='intensity'):
    data = cross if view == 'cross' else intensity
    color = 'RdBu_r' if view == 'cross' else HEAT
    limit = max(1., abs(float(data.min())), abs(float(data.max())))
    lo, hi = (-limit, limit) if view == 'cross' else (0, max(float(data.max()), 1))
    title = 'I − Iᵢₙ꜀' if view == 'cross' else 'I / Iᵢₙ꜀'
    h = go.Figure(go.Heatmap(x=axis, y=axis, z=data, colorscale=color, zmin=lo, zmax=hi,
                            colorbar=dict(title=title, thickness=12)))
    h.update_xaxes(title='x · µm', constrain='domain')
    h.update_yaxes(title='y · µm', scaleanchor='x')
    q = slice(None, None, 2)
    surface = go.Figure(go.Surface(x=axis[q], y=axis[q], z=data[q, q], colorscale=color,
                                  cmin=lo, cmax=hi, showscale=False,
                                  lighting=dict(ambient=.8, diffuse=.5, roughness=.8)))
    ax = dict(gridcolor='#273640', showbackground=False)
    scene = dict(xaxis=dict(**ax, title='x · µm'), yaxis=dict(**ax, title='y · µm'),
                 zaxis=dict(**ax, title=title), aspectratio=dict(x=1, y=1, z=.45),
                 camera=dict(eye=dict(x=1.3, y=1.5, z=1.1)))
    return style(h, 410), style(surface, 410, scene=scene)


def coherence_plot(s):
    opd, gamma = coherence_curve(s)
    fig = go.Figure(go.Scatter(x=opd, y=gamma, line=dict(color=TEAL, width=2)))
    fig.update_xaxes(title='Optical path difference · µm')
    fig.update_yaxes(title='Temporal coherence |γ|', range=[0, 1.05])
    return style(fig, 270)
