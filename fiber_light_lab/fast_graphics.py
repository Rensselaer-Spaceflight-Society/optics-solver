"""GPU display geometry only. Numerical results remain in physics.py.

No solver runs during camera movement. Indexed meshes and float32 buffers keep
the full interference grid inexpensive to draw. OpenGL is optional at startup.
"""
from functools import lru_cache
import numpy as np
from PySide6 import QtCore, QtGui, QtWidgets
from matplotlib import colormaps
from matplotlib.colors import Normalize, LogNorm, TwoSlopeNorm
import pyqtgraph.opengl as gl
from physics import CORE_RADIUS_MM, CORE_INDEX


def configure_surface():
    fmt = QtGui.QSurfaceFormat()
    fmt.setDepthBufferSize(24)
    fmt.setSamples(4)
    fmt.setSwapInterval(1)
    QtGui.QSurfaceFormat.setDefaultFormat(fmt)


def opengl_available():
    """Probe before creating visible widgets, including remote desktop sessions."""
    surface = QtGui.QOffscreenSurface()
    surface.setFormat(QtGui.QSurfaceFormat.defaultFormat())
    surface.create()
    context = QtGui.QOpenGLContext()
    context.setFormat(surface.format())
    if not context.create() or not context.makeCurrent(surface):
        return False
    version = context.format().version()
    context.doneCurrent()
    return version >= (2, 1)


@lru_cache(maxsize=6)
def grid_faces(rows, columns):
    cells = np.arange(rows*columns, dtype=np.uint32).reshape(rows, columns)[:-1, :-1].ravel()
    return np.vstack((np.c_[cells, cells+1, cells+columns],
                      np.c_[cells+1, cells+columns+1, cells+columns])).astype(np.uint32)


class Scene(gl.GLViewWidget):
    def __init__(self, parent=None):
        super().__init__(parent=parent)
        self.setBackgroundColor('white')
        self.setMinimumSize(260, 260)
        self.setSizePolicy(QtWidgets.QSizePolicy.Policy.Expanding, QtWidgets.QSizePolicy.Policy.Expanding)
        self._geometry = None
        self._labels = []
        self._home = dict(pos=QtGui.QVector3D(0, 0, 0), distance=5, elevation=22, azimuth=-55)
        self.setToolTip('Drag to rotate. Scroll to zoom. Middle-drag to pan.')

    def home(self):
        self.setCameraPosition(**self._home)

    def clear(self):
        for _, label in self._labels:
            label.hide()
            label.deleteLater()
        self._labels.clear()
        super().clear()

    def paintGL(self):
        super().paintGL()
        if not self._labels:
            return
        # Ordinary Qt text overlays stay legible across OpenGL drivers and avoid
        # switching between QPainter and GL inside each frame.
        viewport = self.getViewport()
        projection = self.projectionMatrix(viewport, viewport)*self.viewMatrix()
        for pos, label in self._labels:
            point = projection.map(QtGui.QVector3D(*pos))
            if -1 <= point.z() <= 1:
                label.move(round((point.x()+1)*self.width()/2),
                           round((1-point.y())*self.height()/2)-label.height())
                label.show()
            else:
                label.hide()

    def frame(self, key, center, distance):
        self._home = dict(pos=QtGui.QVector3D(*map(float, center)), distance=float(distance),
                          elevation=22, azimuth=-55)
        if self._geometry != key:
            self._geometry = key
            self.home()

    def line(self, points, color=(.45, .45, .45, 1), mode='line_strip', width=1):
        item = gl.GLLinePlotItem(pos=np.asarray(points, dtype=np.float32), color=color,
                                width=width, mode=mode, antialias=True)
        item.setGLOptions('opaque')
        self.addItem(item)

    def label(self, pos, text):
        label = QtWidgets.QLabel(text, self)
        label.setFont(QtGui.QFont('Segoe UI', 9))
        label.setStyleSheet('color: #555555; background: transparent;')
        label.setAttribute(QtCore.Qt.WidgetAttribute.WA_TransparentForMouseEvents)
        label.adjustSize()
        label.hide()
        self._labels.append((tuple(map(float, pos)), label))

    def mesh(self, vertices, faces, colors):
        item = gl.GLMeshItem(vertexes=np.asarray(vertices, dtype=np.float32), faces=faces,
                            vertexColors=np.asarray(colors, dtype=np.float32),
                            smooth=True, computeNormals=False, shader=None)
        self.addItem(item)

    def beam(self, s, result, tip=False, logarithmic=False):
        self.clear()
        a, z = CORE_RADIUS_MM, s.distance_mm
        axis, irradiance, bound = result['footprint']
        segment = 8 if tip else min(12, max(3, z*.16))
        length = min(z, 5) if tip else z
        radius = 3 if tip else max(bound, 1.2)
        # The scene uses physical millimeters on all three axes.
        self.frame(('beam', tip, z, s.output_na), ((length-segment)/2, 0, -radius*.18),
                   max(length+segment, radius*2.3)*1.5)
        theta = np.linspace(.6, 5.4, 65)
        axial, angle = np.meshgrid([-segment, -.25], theta)
        vertices = np.c_[axial.ravel(), 1.1*np.cos(angle).ravel(), 1.1*np.sin(angle).ravel()]
        shade = .65+.16*np.cos(angle.ravel())
        self.mesh(vertices, grid_faces(65, 2), np.c_[shade, shade, shade, np.ones_like(shade)])
        orange = (.71, .42, .13, 1)
        internal = []
        for i in range(8):
            phi = 2*np.pi*i/8
            slope = np.tan(np.arcsin(s.output_na*np.sqrt((i+.5)/8)/CORE_INDEX))
            zi = np.linspace(-segment, 0, 100)
            ri = a-np.abs(np.mod(slope*zi+.12*np.sin(i)+a, 4*a)-2*a)
            p = np.c_[zi, ri*np.cos(phi), ri*np.sin(phi)]
            internal.extend(np.stack((p[:-1], p[1:]), axis=1).reshape(-1, 3))
        self.line(internal, orange, 'lines')
        rays = []
        for i in range(18):
            phi = i*np.pi*(3-np.sqrt(5))
            r0 = .6*a*np.sin(i*2.3)
            r1 = r0+length*np.tan(np.arcsin(s.output_na*np.sqrt((i+.5)/18)))
            rays.extend([[0, r0*np.cos(phi), r0*np.sin(phi)], [length, r1*np.cos(phi), r1*np.sin(phi)]])
        self.line(rays, (.72, .52, .28, 1), 'lines')
        theta = np.linspace(0, 2*np.pi, 193)
        self.line(np.c_[theta*0, a*np.cos(theta), a*np.sin(theta)], (.2, .42, .57, 1), width=2)
        if not tip:
            # Dense polar mesh, interpolated vertex colors; no polygon bands.
            rr, tt = np.meshgrid(np.linspace(0, bound, 65), theta)
            data = np.interp(rr, axis[len(axis)//2:], irradiance[len(axis)//2, len(axis)//2:])*1000
            maximum = max(float(irradiance.max()*1000), 1e-12)
            norm = LogNorm(maximum*1e-4, maximum) if logarithmic else Normalize(0, maximum)
            colors = colormaps['viridis'](norm(np.maximum(data, maximum*1e-4))).reshape(-1, 4)
            self.mesh(np.c_[np.full(rr.size, z), (rr*np.cos(tt)).ravel(), (rr*np.sin(tt)).ravel()],
                      grid_faces(*rr.shape), colors)
            self.line(np.c_[np.full(theta.size, z), bound*np.cos(theta), bound*np.sin(theta)])
            self.label((z, 0, radius*1.13), 'Sample plane')
        self.label((0, 0, max(.9, radius*.14)), 'Tip')
        self.line([[-segment, -radius*1.15, -radius], [length, -radius*1.15, -radius]])
        self.label((0, -radius*1.15, -radius), '0')
        self.label((length, -radius*1.15, -radius), f'{length:g} mm')
        self.line([[length, -radius, -radius*1.2], [length, radius, -radius*1.2]])
        self.label((length, 0, -radius*1.35), f'Width: {2*radius:.1f} mm')
        self.update()

    def wave(self, result, cross=False):
        self.clear()
        axis, intensity, terms = result['wave']
        data = terms if cross else intensity
        norm = TwoSlopeNorm(0, -max(1., abs(float(data.min()))), max(1., float(data.max()))) if cross else Normalize(0, max(1., float(data.max())))
        colors = colormaps['RdBu_r' if cross else 'viridis'](norm(data)).reshape(-1, 4)
        # Display aspect ratio only: x/y normalized to patch width; z is intensity.
        # No filtering or decimation of the 161 x 161 calculated field.
        coords = axis/(np.ptp(axis)/2)
        x, y = np.meshgrid(coords, coords)
        scale = .95/max(1., float(np.ptp(data)), float(np.abs(data).max()))
        self.mesh(np.c_[x.ravel(), y.ravel(), data.ravel()*scale], grid_faces(*data.shape), colors)
        self.frame(('wave',), (0, 0, .35), 3.7)
        self.line([[-1, -1, 0], [1, -1, 0], [1, 1, 0]])
        self.label((0, -1.2, 0), f'x (µm), ±{axis[-1]:.0f}')
        self.label((1.12, 0, 0), f'y (µm), ±{axis[-1]:.0f}')
        lo, hi = min(0., float(data.min())), max(1., float(data.max()))
        self.line([[-1.13, 1, lo*scale], [-1.13, 1, hi*scale]])
        self.label((-1.13, 1, lo*scale), f'{lo:.2g}')
        self.label((-1.13, 1, hi*scale), f'{hi:.2g}')
        self.update()
