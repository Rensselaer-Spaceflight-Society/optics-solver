"""Run locally with a display: python verification/benchmark_desktop.py.

Uses the application's own widgets; saves screenshots and timings. GPU frame
timing includes glFinish, excludes input delivery / compositor / monitor refresh.
"""
import sys
import json
import statistics
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from fast_graphics import configure_surface, opengl_available
from fiber_simulation import Qt, SimulationWindow, set_light_theme, calculate
from physics import Settings
from OpenGL import GL


def settle(window):
    deadline = time.monotonic()+30
    while window._busy or window._pending or window.debounce.isActive():
        app.processEvents()
        time.sleep(.005)
        if time.monotonic() > deadline:
            raise TimeoutError(window.statusBar().currentMessage())
    for _ in range(5):
        app.processEvents()
    assert 'failed' not in window.statusBar().currentMessage().lower()


def median_ms(operation, repeats=20):
    operation()
    times = []
    for _ in range(repeats):
        start = time.perf_counter()
        operation()
        times.append((time.perf_counter()-start)*1000)
    return round(statistics.median(times), 3)


def orbit_ms(scene):
    def render():
        scene.opts['azimuth'] += 1
        scene.makeCurrent()
        # Force an actual draw: QWidget.repaint may reuse its backing texture
        # or skip an occluded window, producing misleadingly small timings.
        ratio = scene.devicePixelRatioF()
        GL.glViewport(0, 0, round(scene.width()*ratio), round(scene.height()*ratio))
        scene.paintGL()
        GL.glFinish()
        scene.doneCurrent()
    result = median_ms(render)
    return result


configure_surface()
app = Qt.QApplication([])
set_light_theme(app)
assert opengl_available(), 'No supported OpenGL context'
window = SimulationWindow(accelerated=True)
window.show()
settle(window)
scene = window.pages[0].scene
assert scene.isValid()
scene.makeCurrent()
report = {name: GL.glGetString(enum).decode() for name, enum in [
    ('vendor', GL.GL_VENDOR), ('renderer', GL.GL_RENDERER), ('version', GL.GL_VERSION)]}
report['multisamples'] = scene.format().samples()
report['framebuffer_samples'] = int(GL.glGetIntegerv(GL.GL_SAMPLES))
scene.doneCurrent()
report['beam_view_pixels'] = [scene.width(), scene.height()]
report['beam_rotation_draw_ms'] = orbit_ms(scene)
scene.home()
app.processEvents()
scene.repaint()
assert scene.grab().save(str(ROOT/'verification/desktop_beam.png'))
assert window.grab().save(str(ROOT/'verification/desktop_window.png'))
window.tabs.setCurrentIndex(2)
window.mode.setCurrentIndex(1)
settle(window)
scene = window.pages[2].scene
report['wave_rotation_draw_ms'] = orbit_ms(scene)
scene.home()
app.processEvents()
scene.repaint()
assert scene.grab().save(str(ROOT/'verification/desktop_wave.png'))
assert window.grab().save(str(ROOT/'verification/desktop_interference.png'))
report['wave_chart_update_ms'] = median_ms(lambda: (window.pages[2].update_wave_chart(window.result, False), window.pages[2].canvas.draw()), repeats=6)
# Same-state revisits must preserve existing charts and camera.
old_axes = window.pages[0].figure.axes[0]
window.tabs.setCurrentIndex(0)
settle(window)
assert window.pages[0].figure.axes[0] is old_axes
# Exercise geometry extremes and every interference representation.
window.distance.setValue(0)
window.power.setValue(0)
settle(window)
assert not window.pages[0].figure.axes[0].images[0].get_array().any()
window.log.setChecked(True)
settle(window)
window.tip.setChecked(True)
settle(window)
window.tabs.setCurrentIndex(1)
settle(window)
window.tabs.setCurrentIndex(2)
window.cross.setChecked(True)
settle(window)
window.mode.setCurrentIndex(0)
settle(window)
window.mode.setCurrentIndex(2)
window.averages.setCurrentText('64')
settle(window)
options = dict(wavelength_nm=650, mode='coherent', averages=16, phase_shift=0, seed=4)
for name, tab, mode in [('beam', 0, 'coherent'), ('coherent', 2, 'coherent'), ('halogen', 2, 'halogen')]:
    report[name+'_calculation_ms'] = median_ms(lambda: calculate(Settings(), tab, {**options, 'mode': mode}))
report['timing_scope'] = 'Warm-cache medians. Forced GL paint plus glFinish; excludes compositor, monitor refresh and input delivery. Rotation performs no physics or 2D redraw. CPU baseline used a 1200x720 full Matplotlib canvas; the new renderer isolates the 3D pane.'
(ROOT/'verification/performance_after.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
print(json.dumps(report, indent=2))
window.close()
