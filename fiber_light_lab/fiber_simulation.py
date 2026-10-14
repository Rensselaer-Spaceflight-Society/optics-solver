"""Native desktop simulation. Run: python fiber_simulation.py

Dependencies are listed in requirements-desktop.txt.
Keep physics.py and data/ alongside this file. No web server is used.
"""
import os
os.environ.setdefault('QT_API', 'pyside6')
os.environ.setdefault('OPENBLAS_NUM_THREADS', '2')

from dataclasses import asdict
from functools import lru_cache
import io
import json
from pathlib import Path
import sys
import time
import zipfile
import numpy as np
from PySide6 import QtCore, QtGui, QtWidgets as Qt
from matplotlib import colormaps, rcParams
from matplotlib.colors import Normalize, LogNorm, TwoSlopeNorm
from matplotlib.figure import Figure
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from physics import (Settings, CORE_RADIUS_MM, CORE_INDEX, transport, footprint,
                     encircled_radius, interference, coherence_curve, reference_data)

ROOT = Path(__file__).resolve().parent
BLUE, ORANGE = '#326b91', '#b56a21'
rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9, 'axes.titlesize': 11,
                 'axes.labelsize': 9, 'axes.spines.top': False, 'axes.spines.right': False,
                 'figure.facecolor': 'white', 'axes.facecolor': 'white',
                 'text.color': '#222222', 'axes.labelcolor': '#333333',
                 'xtick.color': '#555555', 'ytick.color': '#555555',
                 'grid.color': '#dedede', 'grid.linewidth': .6,
                 'savefig.facecolor': 'white'})


def set_light_theme(app):
    """Explicit palette keeps the app light even when Windows uses dark mode."""
    app.setStyle('Fusion')
    palette = QtGui.QPalette()
    for role, value in {
        'Window': '#f4f4f4', 'WindowText': '#202020', 'Base': '#ffffff',
        'AlternateBase': '#f2f2f2', 'Text': '#202020', 'Button': '#eeeeee',
        'ButtonText': '#202020', 'BrightText': '#ffffff', 'Highlight': '#326b91',
        'HighlightedText': '#ffffff', 'ToolTipBase': '#fffbe9',
        'ToolTipText': '#202020', 'Light': '#ffffff', 'Mid': '#b8b8b8',
        'Dark': '#808080', 'Link': '#326b91', 'PlaceholderText': '#737373'
    }.items():
        palette.setColor(getattr(QtGui.QPalette.ColorRole, role), QtGui.QColor(value))
    palette.setColor(QtGui.QPalette.ColorGroup.Disabled, QtGui.QPalette.ColorRole.Text, QtGui.QColor('#888888'))
    palette.setColor(QtGui.QPalette.ColorGroup.Disabled, QtGui.QPalette.ColorRole.ButtonText, QtGui.QColor('#888888'))
    app.setPalette(palette)
    app.setFont(QtGui.QFont('Segoe UI', 10))


class PlotPage(Qt.QWidget):
    def __init__(self, parent=None, accelerated=False):
        super().__init__(parent)
        self.figure = Figure(figsize=(9, 6), layout='constrained')
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.toolbar = NavigationToolbar2QT(self.canvas, self)
        self.extra = Qt.QHBoxLayout()
        layout = Qt.QVBoxLayout(self)
        layout.setContentsMargins(4, 4, 4, 4)
        layout.addLayout(self.extra)
        self.scene = None
        self.chart_key = None
        if accelerated:
            from fast_graphics import Scene
            self.scene = Scene()
            self.splitter = Qt.QSplitter(QtCore.Qt.Orientation.Horizontal)
            scene_panel = Qt.QWidget()
            scene_layout = Qt.QVBoxLayout(scene_panel)
            scene_layout.setContentsMargins(0, 0, 0, 0)
            self.scene_title = Qt.QLabel('Output geometry')
            self.scene_title.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
            scene_layout.addWidget(self.scene_title)
            scene_layout.addWidget(self.scene, 1)
            buttons = Qt.QHBoxLayout()
            home = Qt.QPushButton('Reset view')
            home.clicked.connect(self.scene.home)
            save = Qt.QPushButton('Save 3D image…')
            save.clicked.connect(self.save_scene)
            buttons.addWidget(home); buttons.addWidget(save); buttons.addStretch()
            scene_layout.addLayout(buttons)
            plot_panel = Qt.QWidget()
            plot_layout = Qt.QVBoxLayout(plot_panel)
            plot_layout.setContentsMargins(0, 0, 0, 0)
            plot_layout.addWidget(self.canvas, 1)
            plot_layout.addWidget(self.toolbar)
            self.splitter.addWidget(scene_panel)
            self.splitter.addWidget(plot_panel)
            self.splitter.setSizes([550, 450])
            self.splitter.setChildrenCollapsible(False)
            layout.addWidget(self.splitter, 1)
        else:
            layout.addWidget(self.canvas, 1)
            layout.addWidget(self.toolbar)

    def save_scene(self):
        path, _ = Qt.QFileDialog.getSaveFileName(self, 'Save 3D view', 'fiber_3d.png', 'PNG image (*.png)')
        if path and not self.scene.grab().save(path):
            Qt.QMessageBox.warning(self, 'Could not save', 'The image could not be written.')

    def update_output_chart(self, s, result, logarithmic):
        key = (s, logarithmic)
        if self.chart_key == key:
            return False
        if self.chart_key is None or self.chart_key[1] != logarithmic:
            draw_output(self.figure, s, result, False, logarithmic, charts_only=True)
        else:
            axis, e, _ = result['footprint']
            data = e*1000
            maximum = max(float(data.max()), 1e-12)
            heat, _, profile = self.figure.axes
            image = heat.images[0]
            image.set_data(np.maximum(data, maximum*1e-4) if logarithmic else data)
            image.set_extent([axis[0], axis[-1], axis[0], axis[-1]])
            image.set_clim(maximum*1e-4 if logarithmic else 0, maximum)
            profile.lines[0].set_data(axis, data[len(axis)//2])
            profile.set(xlim=(axis[0], axis[-1]), ylim=(0, maximum*1.08))
        self.chart_key = key
        return True

    def update_wave_chart(self, result, cross):
        if self.chart_key != cross:
            draw_interference(self.figure, result, cross, charts_only=True)
        else:
            axis, intensity, terms = result['wave']
            data = terms if cross else intensity
            image = self.figure.axes[0].images[0]
            image.set_data(data)
            image.set_extent([axis[0], axis[-1], axis[0], axis[-1]])
            image.set_clim(-max(1., abs(float(data.min()))) if cross else 0, max(1., float(data.max())))
        self.chart_key = cross


def configure_3d(ax):
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.set_pane_color((.98, .98, .98, 1))
    ax.tick_params(labelsize=8, pad=1)
    ax.set_xlabel('Distance from tip (mm)', labelpad=6)
    ax.set_ylabel('x (mm)', labelpad=2)
    ax.set_zlabel('y (mm)', labelpad=2)
    ax.view_init(elev=22, azim=-62)


def draw_beam(ax, settings, tip=False, logarithmic=False):
    """Low polygon scene; dragging rotates it without invoking the solver."""
    z = settings.distance_mm
    angle = np.arcsin(settings.output_na)
    a = CORE_RADIUS_MM
    bound = a+z*np.tan(angle)
    segment = 8 if tip else min(12, max(3, z*.16))
    th = np.linspace(.6, 5.4, 25)
    axial, theta = np.meshgrid([-segment, -.25], th)
    ax.plot_surface(axial, 1.1*np.cos(theta), 1.1*np.sin(theta),
                    color='#999999', alpha=.7, linewidth=0, shade=True)
    # A few exactly reflected meridional rays in the straight core.
    rays = []
    for i in range(8):
        phi = 2*np.pi*i/8
        slope = np.tan(np.arcsin(settings.output_na*np.sqrt((i+.5)/8)/CORE_INDEX))
        zi = np.linspace(-segment, 0, 100)
        unfolded = slope*zi+.12*np.sin(i)
        ri = a-np.abs(np.mod(unfolded+a, 4*a)-2*a)
        rays.append(np.c_[zi, ri*np.cos(phi), ri*np.sin(phi)])
    ax.add_collection3d(Line3DCollection(rays, colors=ORANGE, linewidths=.8, alpha=.8))
    th = np.linspace(0, 2*np.pi, 64)
    ax.plot(np.zeros_like(th), a*np.cos(th), a*np.sin(th), color=BLUE, linewidth=1.5)
    length = min(z, 5) if tip else z
    segments = []
    for i in range(18):
        phi = i*np.pi*(3-np.sqrt(5))
        slope = np.tan(np.arcsin(settings.output_na*np.sqrt((i+.5)/18)))
        r0 = .6*a*np.sin(i*2.3)
        r1 = r0+length*slope
        segments.append([[0, r0*np.cos(phi), r0*np.sin(phi)],
                         [length, r1*np.cos(phi), r1*np.sin(phi)]])
    ax.add_collection3d(Line3DCollection(segments, colors=ORANGE, linewidths=.7, alpha=.55))
    if not tip:
        # Polar grid makes the boundary circular without thousands of polygons.
        axis, e, _ = footprint(settings, 129)
        radius, theta = np.meshgrid(np.linspace(0, bound, 18), np.linspace(0, 2*np.pi, 49))
        radial = np.interp(radius, axis[len(axis)//2:], e[len(axis)//2, len(axis)//2:])*1000
        maximum = max(float(e.max()*1000), 1e-12)
        norm = LogNorm(maximum*1e-4, maximum) if logarithmic else Normalize(0, maximum)
        ax.plot_surface(np.full_like(radius, z), radius*np.cos(theta), radius*np.sin(theta),
                        facecolors=colormaps['viridis'](norm(np.maximum(radial, maximum*1e-4))),
                        rstride=1, cstride=1, shade=False, antialiased=False, linewidth=0)
        ax.plot(np.full_like(th, z), bound*np.cos(th), bound*np.sin(th), color='#707070', linewidth=.6)
        ax.text(z, 0, bound*1.1, 'Sample plane', fontsize=9)
        ax.set_xlim(-segment, max(z, 1))
        ax.set_ylim(-bound*1.15, bound*1.15)
        ax.set_zlim(-bound*1.15, bound*1.15)
        ax.set_box_aspect((z+segment, 2.3*bound, 2.3*bound))
    else:
        ax.set_xlim(-8, 5); ax.set_ylim(-3, 3); ax.set_zlim(-3, 3)
        ax.set_box_aspect((13, 6, 6))
    ax.text(0, 0, .9 if tip else bound*.14, 'Tip', fontsize=9)
    ax.set_title('Fiber tip (last 8 mm)' if tip else 'Output geometry')
    configure_3d(ax)


def draw_output(figure, s, result, tip, logarithmic, camera=None, charts_only=False):
    figure.clear()
    figure.set_layout_engine('constrained')
    grid = figure.add_gridspec(2, 1 if charts_only else 2, height_ratios=[3, 1.25])
    ax = None
    if not charts_only:
        ax = figure.add_subplot(grid[:, 0], projection='3d')
        draw_beam(ax, s, tip, logarithmic)
        if camera is not None:
            ax.view_init(*camera)
    x, e, _ = result['footprint']
    mapax = figure.add_subplot(grid[0, -1])
    data = e*1000
    maximum = max(float(data.max()), 1e-12)
    norm = LogNorm(maximum*1e-4, maximum) if logarithmic else Normalize(0, maximum)
    image = mapax.imshow(np.maximum(data, maximum*1e-4) if logarithmic else data,
                        origin='lower', extent=[x[0], x[-1], x[0], x[-1]],
                        cmap='viridis', norm=norm, interpolation='bilinear')
    mapax.set(title='Irradiance at sample', xlabel='x (mm)', ylabel='y (mm)')
    figure.colorbar(image, ax=mapax, shrink=.87, label='W/m²', pad=.03)
    profile = figure.add_subplot(grid[1, -1])
    profile.plot(x, data[len(x)//2], color=BLUE)
    profile.set(xlabel='x through center (mm)', ylabel='W/m²', ylim=(0, maximum*1.08))
    profile.grid(True, alpha=.6)
    return ax


def draw_transport(figure, s, result):
    figure.clear()
    figure.set_layout_engine('constrained')
    axes = figure.subplots(2, 2)
    t = result['transport']
    ax = axes[0, 0]
    for key, label, color in [('incident_spectrum', 'Incident', '#888888'),
                               ('guided_spectrum', 'Inside fiber', BLUE), ('output_spectrum', 'Output', ORANGE)]:
        ax.plot(t['wavelength_nm'], t[key]*1000, label=label, color=color)
    ax.set(title='Source and output spectrum', xlabel='Wavelength (nm)', ylabel='µW/nm')
    ax.legend(frameon=False, fontsize=8)
    ax = axes[0, 1]
    lengths = np.linspace(0, max(.3, s.length_m*1.1), 100)
    power = np.trapezoid(t['guided_spectrum'][:, None]*10**(-t['loss_db_per_m'][:, None]*lengths/10), t['wavelength_nm'], axis=0)
    ax.plot(lengths, power, color=BLUE)
    ax.axvline(s.length_m, color=ORANGE, linestyle=':', label='Cable end')
    ax.set(title='Power inside the fiber', xlabel='Distance along cable (m)', ylabel='mW', ylim=(0, max(float(power.max())*1.1, 1e-9)))
    ax.legend(frameon=False, fontsize=8)
    ax = axes[1, 0]
    _, loss = reference_data()
    ax.semilogy(loss[:, 0], loss[:, 1], color=BLUE, label='Approximate typical curve')
    ax.plot(650, .17, 'o', color=ORANGE, label='650 nm specification maximum')
    ax.set(title='GH4001 attenuation', xlabel='Wavelength (nm)', ylabel='dB/m')
    ax.legend(frameon=False, fontsize=8)
    ax = axes[1, 1]
    opd, gamma = result['coherence']
    ax.plot(opd, gamma, color=BLUE)
    ax.set(title='Temporal coherence of output', xlabel='Optical path difference (µm)', ylabel='|γ|', ylim=(0, 1.05))
    for ax in axes.flat:
        ax.grid(True, alpha=.6)


def draw_interference(figure, result, cross, camera=None, charts_only=False):
    figure.clear()
    figure.set_layout_engine('constrained')
    axis, intensity, terms = result['wave']
    data = terms if cross else intensity
    label = '(I − Iinc) / Iinc' if cross else 'I / Iinc'
    norm = TwoSlopeNorm(0, -max(1., abs(float(data.min()))), max(1., float(data.max()))) if cross else Normalize(0, max(1., float(data.max())))
    cmap = 'RdBu_r' if cross else 'viridis'
    heat = figure.add_subplot(111 if charts_only else 121)
    image = heat.imshow(data, origin='lower', extent=[axis[0], axis[-1], axis[0], axis[-1]], cmap=cmap, norm=norm, interpolation='bilinear')
    heat.set(title='Interference contribution' if cross else 'Local intensity', xlabel='x (µm)', ylabel='y (µm)')
    figure.colorbar(image, ax=heat, shrink=.65, label=label, pad=.02)
    if charts_only:
        return None
    ax = figure.add_subplot(122, projection='3d')
    # ~1,600 quads, rather than the 25,921-cell computational grid.
    x, y = np.meshgrid(axis[::4], axis[::4])
    ax.plot_surface(x, y, data[::4, ::4], cmap=cmap, norm=norm, rstride=1, cstride=1,
                    linewidth=0, antialiased=False)
    configure_3d(ax)
    ax.set(xlabel='x (µm)', ylabel='y (µm)', zlabel=label, title='Intensity surface')
    ax.set_box_aspect((1, 1, .6))
    if camera is not None:
        ax.view_init(*camera)
    return ax


@lru_cache(maxsize=12)
def base_calculation(s):
    t = transport(s)
    axis, e, radius = footprint(s)
    return {'transport': t, 'footprint': (axis, e, radius),
            'diameter90': 2*encircled_radius(axis, e)}


@lru_cache(maxsize=8)
def cached_coherence(s):
    return coherence_curve(s)


@lru_cache(maxsize=8)
def cached_wave(length, distance, na, options):
    return interference(Settings(length_m=length, distance_mm=distance, launch_na=na), **dict(options))


def calculate(s, tab, wave_options):
    result = dict(base_calculation(s))
    if tab == 1:
        result['coherence'] = cached_coherence(s)
    if tab == 2:
        result['wave'] = cached_wave(s.length_m, s.distance_mm, s.output_na, tuple(sorted(wave_options.items())))
    return result


class WorkerSignals(QtCore.QObject):
    ready = QtCore.Signal(int, object, float)
    failed = QtCore.Signal(int, str)


class Calculation(QtCore.QRunnable):
    def __init__(self, generation, settings, tab, options):
        super().__init__()
        self.signals = WorkerSignals()
        self.generation, self.settings, self.tab, self.options = generation, settings, tab, options

    @QtCore.Slot()
    def run(self):
        start = time.perf_counter()
        try:
            result = calculate(self.settings, self.tab, self.options)
            self.signals.ready.emit(self.generation, result, time.perf_counter()-start)
        except Exception as exc:
            self.signals.failed.emit(self.generation, str(exc))


def save_data(path, settings, wave_options):
    """Desktop export is independent of Dash and Plotly."""
    result = calculate(settings, 2, wave_options)
    t = result['transport']
    ax, e, _ = result['footprint']; x, y = np.meshgrid(ax, ax)
    wa, wi, wc = result['wave']; wx, wy = np.meshgrid(wa, wa)

    def csv(data, header):
        buffer = io.StringIO()
        np.savetxt(buffer, data, delimiter=',', header=header, comments='')
        return buffer.getvalue()

    with zipfile.ZipFile(path, 'w', zipfile.ZIP_DEFLATED) as z:
        z.writestr('settings.json', json.dumps({'settings': asdict(settings), 'wave': wave_options,
                        'output_mw': t['output_mw'], 'power_status': 'Assumed; not calibrated'}, indent=2))
        z.writestr('irradiance.csv', csv(np.c_[x.ravel(), y.ravel(), e.ravel()*1000], 'x_mm,y_mm,irradiance_W_per_m2'))
        z.writestr('spectrum.csv', csv(np.c_[t['wavelength_nm'], t['incident_spectrum'], t['guided_spectrum'], t['output_spectrum']],
                                      'wavelength_nm,incident_mW_per_nm,guided_mW_per_nm,output_mW_per_nm'))
        z.writestr('illustrative_wave.csv', csv(np.c_[wx.ravel(), wy.ravel(), wi.ravel(), wc.ravel()],
                                              'x_um,y_um,I_over_incoherent,cross_term_over_incoherent'))
        for name in ['eke_relative_spectrum.csv', 'gh4001_loss.csv']:
            z.write(ROOT/'data'/name, 'reference_curves/'+name)
        z.write(ROOT/'MODEL.md', 'MODEL.md')


class SimulationWindow(Qt.QMainWindow):
    def __init__(self, accelerated=False):
        super().__init__()
        self.setWindowTitle('Fiber optical simulation')
        self.resize(1280, 820)
        self.setMinimumSize(960, 640)
        self._generation, self._busy, self._pending = 0, False, False
        self._closing = False
        self.accelerated = accelerated
        self._rendered = {}
        self.seed = 4
        self.result = None
        self.jobs = QtCore.QThreadPool(self)
        self.jobs.setMaxThreadCount(1)
        self.debounce = QtCore.QTimer(self)
        self.debounce.setSingleShot(True)
        self.debounce.setInterval(70)
        self.debounce.timeout.connect(self.start_calculation)
        self.cameras = {}
        self.axes3d = {}
        root = Qt.QWidget()
        self.setCentralWidget(root)
        layout = Qt.QHBoxLayout(root)
        sidebar = Qt.QWidget()
        sidebar.setFixedWidth(245)
        controls = Qt.QVBoxLayout(sidebar)
        controls.setContentsMargins(8, 8, 12, 8)
        title = Qt.QLabel('MI-152 / GH4001')
        title.setFont(QtGui.QFont('Segoe UI', 12, QtGui.QFont.Weight.DemiBold))
        controls.addWidget(title)
        controls.addWidget(Qt.QLabel('980 µm core · NA 0.50'))
        controls.addSpacing(16)
        form = Qt.QFormLayout()
        form.setVerticalSpacing(14)
        self.length = self.spin(10, 0, 400, .25, 2, 'Cable length (in)')
        self.distance = self.spin(2.5, 0, 12, .05, 2, 'Tip to sample (in)')
        self.power = self.spin(1, 0, 1000, .1, 3, 'Input power (mW)')
        self.na = self.spin(.5, .05, .9, .01, 2, 'Launch NA')
        for label, widget in [('Cable (in)', self.length), ('Distance (in)', self.distance),
                              ('Input (mW)', self.power), ('Launch NA', self.na)]:
            form.addRow(label, widget)
        controls.addLayout(form)
        note = Qt.QLabel('Input power is assumed, in the selected wavelength band. The lamp’s 150 W rating is electrical.')
        note.setWordWrap(True)
        note.setStyleSheet('color: #606060; font-size: 11px;')
        controls.addSpacing(10)
        controls.addWidget(note)
        controls.addSpacing(12)
        advanced = Qt.QToolButton()
        advanced.setText('Additional settings')
        advanced.setCheckable(True)
        advanced.setToolButtonStyle(QtCore.Qt.ToolButtonStyle.ToolButtonTextBesideIcon)
        advanced.setArrowType(QtCore.Qt.ArrowType.RightArrow)
        controls.addWidget(advanced)
        extra = Qt.QWidget()
        extraform = Qt.QFormLayout(extra)
        extraform.setContentsMargins(0, 10, 0, 0)
        self.spot = self.spin(.4, .05, 2, .05, 2, 'Spot radius (mm)')
        self.offset = self.spin(0, 0, 1.5, .05, 2, 'Offset (mm)')
        self.loss = self.spin(0, 0, 10, .1, 1, 'Extra loss (dB)')
        self.band = Qt.QComboBox()
        self.band.addItems(['400–1000 nm', '400–700 nm'])
        self.band.setAccessibleName('Wavelength band')
        self.band.currentIndexChanged.connect(self.schedule)
        for label, widget in [('Spot radius (mm)', self.spot), ('Offset (mm)', self.offset),
                              ('Extra loss (dB)', self.loss), ('Band', self.band)]:
            extraform.addRow(label, widget)
        self.spot.setToolTip('Gaussian 1/e² intensity radius at the input face.')
        extra.hide()
        advanced.toggled.connect(extra.setVisible)
        advanced.toggled.connect(lambda checked: advanced.setArrowType(QtCore.Qt.ArrowType.DownArrow if checked else QtCore.Qt.ArrowType.RightArrow))
        controls.addWidget(extra)
        controls.addStretch()
        self.readout = Qt.QLabel('Calculating…')
        self.readout.setWordWrap(True)
        self.readout.setTextInteractionFlags(QtCore.Qt.TextInteractionFlag.TextSelectableByMouse)
        controls.addWidget(self.readout)
        controls.addSpacing(12)
        for text, method in [('Reset settings', self.reset), ('Save data…', self.export), ('Model notes', self.notes)]:
            button = Qt.QPushButton(text)
            button.clicked.connect(method)
            controls.addWidget(button)
        layout.addWidget(sidebar)
        self.tabs = Qt.QTabWidget()
        layout.addWidget(self.tabs, 1)
        self.pages = [PlotPage(accelerated=accelerated and i != 1) for i in range(3)]
        for name, page in zip(['Beam', 'Fiber & spectrum', 'Interference'], self.pages):
            self.tabs.addTab(page, name)
        self.tip = Qt.QCheckBox('Magnify fiber tip')
        self.log = Qt.QCheckBox('Log intensity')
        for widget in [self.tip, self.log]:
            widget.toggled.connect(self.schedule)
            self.pages[0].extra.addWidget(widget)
        self.pages[0].extra.addStretch()
        self.pages[0].extra.addWidget(Qt.QLabel('Drag 3D view to rotate'))
        caption = Qt.QLabel('Approximate spectra; straight fiber in air. See Model notes for assumptions.')
        self.pages[1].extra.addWidget(caption)
        self.mode = Qt.QComboBox()
        self.mode.addItem('Halogen (incoherent mean)', 'halogen')
        self.mode.addItem('Coherent comparison', 'coherent')
        self.mode.addItem('Average independent patterns', 'ensemble')
        self.mode.setAccessibleName('Interference model')
        self.mode.currentIndexChanged.connect(self.wave_controls_changed)
        self.wavelength = self.spin(650, 400, 1000, 5, 0, 'Wavelength (nm)')
        self.wavelength.setSuffix(' nm')
        self.wavelength.setMaximumWidth(105)
        self.averages = Qt.QComboBox()
        self.averages.addItems(['1', '4', '16', '64'])
        self.averages.setCurrentText('16')
        self.averages.currentIndexChanged.connect(self.schedule)
        self.averages.setAccessibleName('Independent patterns')
        self.cross = Qt.QCheckBox('Show interference contribution')
        self.cross.toggled.connect(self.schedule)
        wavebar = self.pages[2].extra
        wavebar.addWidget(self.mode)
        wavebar.addWidget(self.wavelength)
        wavebar.addWidget(Qt.QLabel('Patterns:'))
        wavebar.addWidget(self.averages)
        wavebar.addStretch()
        wavebar.addWidget(self.cross)
        wavebottom = Qt.QHBoxLayout()
        self.phase = Qt.QSlider(QtCore.Qt.Orientation.Horizontal)
        self.phase.setRange(0, 1000)
        self.phase.setAccessibleName('Axial phase shift (micrometers)')
        self.phase.setTracking(False)
        self.phase.valueChanged.connect(self.schedule)
        self.phase.valueChanged.connect(lambda value: self.phase_label.setText(f'Phase shift: {value} µm'))
        self.phase_label = Qt.QLabel('Phase shift: 0 µm')
        wavebottom.addWidget(self.phase_label)
        wavebottom.addWidget(self.phase, 1)
        self.new_phase = Qt.QPushButton('New pattern')
        self.new_phase.clicked.connect(self.reseed)
        wavebottom.addWidget(self.new_phase)
        self.pages[2].layout().insertLayout(1, wavebottom)
        self.wave_note = Qt.QLabel('Illustrative local field; not the exact speckle pattern of this cable.')
        self.wave_note.setWordWrap(True)
        self.pages[2].layout().insertWidget(2, self.wave_note)
        self.tabs.currentChanged.connect(self.schedule)
        self.wave_controls_changed()
        self.schedule()

    def spin(self, value, lo, hi, step, decimals, accessible_name):
        widget = Qt.QDoubleSpinBox()
        widget.setRange(lo, hi); widget.setDecimals(decimals)
        widget.setSingleStep(step); widget.setValue(value)
        widget.setKeyboardTracking(False)
        widget.setAccessibleName(accessible_name)
        widget.valueChanged.connect(self.schedule)
        return widget

    def settings(self):
        return Settings(self.length.value()*.0254, self.distance.value()*25.4,
                        self.power.value(), self.na.value(), self.spot.value(),
                        self.offset.value(), self.loss.value(), self.band.currentText())

    def wave_options(self):
        return dict(wavelength_nm=self.wavelength.value(), mode=self.mode.currentData(),
                    averages=int(self.averages.currentText()), phase_shift=self.phase.value(), seed=self.seed)

    def wave_controls_changed(self, *_):
        coherent = self.mode.currentData() != 'halogen'
        self.averages.setEnabled(self.mode.currentData() == 'ensemble')
        self.phase.setEnabled(coherent)
        self.new_phase.setEnabled(coherent)
        self.wavelength.setEnabled(coherent)
        self.schedule()

    def reseed(self):
        self.seed += 1
        self.schedule()

    def reset(self):
        for widget, value in [(self.length, 10), (self.distance, 2.5), (self.power, 1),
                               (self.na, .5), (self.spot, .4), (self.offset, 0), (self.loss, 0),
                               (self.wavelength, 650), (self.phase, 0)]:
            widget.setValue(value)
        self.band.setCurrentIndex(0); self.mode.setCurrentIndex(0)
        self.averages.setCurrentText('16'); self.cross.setChecked(False)
        self.tip.setChecked(False); self.log.setChecked(False)
        self.cameras.clear(); self.axes3d.clear(); self.seed = 4
        self._rendered.clear()
        for page in self.pages:
            if page.scene is not None:
                page.scene._geometry = None
        self.schedule()

    def schedule(self, *_):
        self._generation += 1
        self._pending = True
        self.debounce.start()

    def start_calculation(self):
        if self._busy or self._closing:
            return
        # Capture cameras before clearing any axes, preserving user rotation.
        for key, ax in self.axes3d.items():
            self.cameras[key] = (ax.elev, ax.azim, ax.roll)
        self._pending, self._busy = False, True
        s, tab = self.settings(), self.tabs.currentIndex()
        key = (s, self.tip.isChecked(), self.log.isChecked()) if tab == 0 else ((s,) if tab == 1 else (s, tuple(sorted(self.wave_options().items())), self.cross.isChecked()))
        if self._rendered.get(tab) == key:
            self._busy = False
            self.statusBar().showMessage('Ready · view retained   |   Optical power is assumed')
            return
        self._request_key = key
        self._request = (s, tab, self.tip.isChecked(), self.log.isChecked(), self.cross.isChecked())
        self.statusBar().showMessage('Calculating…')
        self.worker = Calculation(self._generation, s, tab, self.wave_options())
        self.worker.signals.ready.connect(self.finished)
        self.worker.signals.failed.connect(self.failed)
        self.jobs.start(self.worker)

    @QtCore.Slot(int, object, float)
    def finished(self, generation, result, calculation_time):
        self._busy = False
        if self._closing:
            return
        if generation != self._generation:
            self.start_calculation()
            return
        s, tab, tip, log, cross = self._request
        self.result = result
        t, (_, e, _) = result['transport'], result['footprint']
        self.readout.setText(f"Output: {t['output_mw']:.4g} mW\nPeak: {e.max()*1000:.4g} W/m²\n"
                             f"90% diameter: {result['diameter90']:.2f} mm\nHalf-angle: {t['half_angle_deg']:.1f}°\n\n"
                             f'{s.length_m*1000:g} mm cable\n{s.distance_mm:g} mm to sample')
        page = self.pages[tab]
        start = time.perf_counter()
        draw_chart = True
        if tab == 0:
            if page.scene is not None:
                page.scene.beam(s, result, tip, log)
                page.scene_title.setText('Fiber tip (last 8 mm)' if tip else 'Output geometry · mm')
                draw_chart = page.update_output_chart(s, result, log)
            else:
                self.axes3d[0] = draw_output(page.figure, s, result, tip, log, self.cameras.get(0))
        elif tab == 1:
            draw_transport(page.figure, s, result)
        else:
            if page.scene is not None:
                page.scene.wave(result, cross)
                page.scene_title.setText('Interference contribution · (I − Iinc) / Iinc' if cross else 'Intensity surface · I / Iinc')
                page.update_wave_chart(result, cross)
            else:
                self.axes3d[2] = draw_interference(page.figure, result, cross, self.cameras.get(2))
            axis, intensity, _ = result['wave']
            contrast = intensity.std()/intensity.mean()
            self.wave_note.setText(f'Illustrative local field; not the exact pattern of this cable.  '
                                  f'Patch: {np.ptp(axis):.0f} µm · contrast: {contrast:.3f}\n'
                                  'Surface height is relative intensity. The phase control holds the beam envelope fixed.')
        # Clear old toolbar view limits when new physical geometry is drawn.
        if draw_chart:
            page.toolbar.update()
            page.canvas.draw()
        # Keep the completed layout during mouse rotation. Re-solving label
        # placement on every mouse move adds unnecessary work to mplot3d.
        page.figure.set_layout_engine(None)
        elapsed = time.perf_counter()-start
        self._rendered[tab] = self._request_key
        backend = 'OpenGL 3D' if self.accelerated else 'Software 3D'
        self.statusBar().showMessage(f'{backend} · calculation {calculation_time*1000:.0f} ms · plot update {elapsed*1000:.0f} ms   |   Optical power is assumed')
        if self._pending:
            self.debounce.start()

    @QtCore.Slot(int, str)
    def failed(self, generation, message):
        self._busy = False
        if generation == self._generation:
            self.statusBar().showMessage('Calculation failed: '+message)
        if self._pending:
            self.debounce.start()

    def export(self):
        name, _ = Qt.QFileDialog.getSaveFileName(self, 'Save simulation data', str(ROOT/'simulation_data.zip'), 'ZIP archive (*.zip)')
        if not name:
            return
        try:
            save_data(name, self.settings(), self.wave_options())
            self.statusBar().showMessage('Saved '+name)
        except Exception as exc:
            Qt.QMessageBox.warning(self, 'Could not save', str(exc))

    def notes(self):
        dialog = Qt.QDialog(self)
        dialog.setWindowTitle('Model notes')
        dialog.resize(820, 680)
        layout = Qt.QVBoxLayout(dialog)
        text = Qt.QTextBrowser()
        text.setOpenExternalLinks(True)
        text.setMarkdown((ROOT/'MODEL.md').read_text(encoding='utf-8'))
        layout.addWidget(text)
        close = Qt.QDialogButtonBox(Qt.QDialogButtonBox.StandardButton.Close)
        close.rejected.connect(dialog.reject)
        layout.addWidget(close)
        dialog.exec()

    def closeEvent(self, event):
        self._closing = True
        self.debounce.stop()
        self.jobs.waitForDone()
        event.accept()


def main():
    accelerated = False
    if '--software' not in sys.argv and os.environ.get('FIBER_RENDERER') != 'software':
        try:
            from fast_graphics import configure_surface, opengl_available
            configure_surface()
            accelerated = True
        except ImportError:
            pass
    app = Qt.QApplication(sys.argv)
    app.setApplicationName('Fiber optical simulation')
    set_light_theme(app)
    if accelerated:
        accelerated = opengl_available()
    window = SimulationWindow(accelerated=accelerated)
    window.show()
    return app.exec()


if __name__ == '__main__':
    sys.exit(main())
