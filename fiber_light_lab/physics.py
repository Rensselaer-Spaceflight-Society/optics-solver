"""Radiometric transport and a deliberately separate local coherence surrogate.

Length units are mm in the output geometry, m in attenuation, um in the wave
patch, nm in spectra. Optical powers are mW; irradiance is mW/mm².
See MODEL.md for equations, provenance, and limitations.
"""
from dataclasses import dataclass, asdict
from functools import lru_cache
from pathlib import Path
import numpy as np
from scipy.signal import fftconvolve
from scipy.stats import ncx2

DATA = Path(__file__).parent / 'data'
CORE_RADIUS_MM = 0.490
CORE_INDEX = 1.49
FIBER_NA = 0.50


@dataclass(frozen=True)
class Settings:
    length_m: float = 0.254
    distance_mm: float = 63.5
    incident_mw: float = 1.0
    launch_na: float = 0.50
    spot_radius_mm: float = 0.40
    offset_mm: float = 0.0
    extra_loss_db: float = 0.0
    band: str = '400–1000 nm'

    def __post_init__(self):
        vals = [v for v in asdict(self).values() if isinstance(v, (float, int))]
        if not all(np.isfinite(v) for v in vals):
            raise ValueError('Settings must be finite.')
        if self.length_m < 0 or self.distance_mm < 0 or self.incident_mw < 0:
            raise ValueError('Length, distance and power must be nonnegative.')
        if not 0 < self.launch_na < 1 or self.spot_radius_mm <= 0 or self.extra_loss_db < 0:
            raise ValueError('Invalid launch NA, spot radius or loss.')
        if self.band not in ('400–1000 nm', '400–700 nm'):
            raise ValueError('Unsupported spectral band.')

    @property
    def output_na(self):
        return min(self.launch_na, FIBER_NA)


def fresnel_transmission(sin_external):
    """Unpolarized power transmission of one lossless air/PMMA interface."""
    s = np.asarray(sin_external)
    c1, c2 = np.sqrt(1-s*s), np.sqrt(1-(s/CORE_INDEX)**2)
    rs = ((c1-CORE_INDEX*c2)/(c1+CORE_INDEX*c2))**2
    rp = ((CORE_INDEX*c1-c2)/(CORE_INDEX*c1+c2))**2
    return 1 - (rs+rp)/2


def spatial_coupling(radius_mm, offset_mm):
    # A displaced circular Gaussian: sigma=w/2, w is the 1/e² radius.
    sigma = radius_mm / 2
    return float(ncx2.cdf((CORE_RADIUS_MM/sigma)**2, 2, (offset_mm/sigma)**2))


@lru_cache(maxsize=1)
def reference_data():
    return (np.loadtxt(DATA/'eke_relative_spectrum.csv', delimiter=',', skiprows=1),
            np.loadtxt(DATA/'gh4001_loss.csv', delimiter=',', skiprows=1))


@lru_cache(maxsize=64)
def transport(s: Settings):
    source, loss = reference_data()
    upper = 700 if s.band == '400–700 nm' else 1000
    wl = np.arange(400., upper+0.1, 1.)
    relative = np.interp(wl, source[:, 0], source[:, 1])
    input_density = relative / np.trapezoid(relative, wl) * s.incident_mw
    alpha = np.exp(np.interp(wl, loss[:, 0], np.log(loss[:, 1])))
    angular_fraction = min(1., (FIBER_NA/s.launch_na)**2)
    spatial_fraction = spatial_coupling(s.spot_radius_mm, s.offset_mm)
    u = (np.arange(2048)+0.5)/2048
    t = fresnel_transmission(s.output_na*np.sqrt(u))
    t1, t2 = float(t.mean()), float((t*t).mean())
    acceptance = angular_fraction*spatial_fraction
    guided = input_density*acceptance*t1
    bulk_transmission = 10**(-(alpha*s.length_m+s.extra_loss_db)/10)
    output = input_density*acceptance*t2*bulk_transmission
    pout = float(np.trapezoid(output, wl))
    return dict(wavelength_nm=wl, incident_spectrum=input_density,
                guided_spectrum=guided, output_spectrum=output, loss_db_per_m=alpha,
                output_mw=pout, guided_mw=float(np.trapezoid(guided, wl)),
                coupling=acceptance*t1, spatial_fraction=spatial_fraction,
                angular_fraction=angular_fraction, face_pair_transmission=t2,
                half_angle_deg=float(np.degrees(np.arcsin(s.output_na))))


@lru_cache(maxsize=48)
def normalized_footprint(distance_mm, output_na, n=257):
    """Convolve a uniform exit disk with a radiance-preserving angular kernel.

    Uniform projected solid-angle density: p(sx,sy)=constant in an NA disk.
    Its plane density is z²/[pi NA² (z²+r²)²]. Two Fresnel transmissions
    weight this kernel. Raster probabilities are normalized before convolution.
    """
    if n < 33 or n % 2 == 0:
        raise ValueError('Use an odd grid size >=33.')
    if distance_mm < 0 or not 0 < output_na < 1:
        raise ValueError('Invalid propagation geometry.')
    rmax = CORE_RADIUS_MM + distance_mm*np.tan(np.arcsin(output_na))
    axis = np.linspace(-1.12*rmax, 1.12*rmax, n)
    dx = axis[1]-axis[0]
    x, y = np.meshgrid(axis, axis)
    rr = x*x+y*y
    disk = (rr <= CORE_RADIUS_MM**2).astype(float)
    disk /= disk.sum()
    if distance_mm == 0:
        density = disk/(dx*dx)
    else:
        sin_angle = np.sqrt(rr/(rr+distance_mm**2))
        kernel = distance_mm**2/(rr+distance_mm**2)**2
        kernel *= fresnel_transmission(sin_angle)**2
        kernel[sin_angle > output_na] = 0
        kernel /= kernel.sum()
        density = np.maximum(fftconvolve(disk, kernel, mode='same'), 0)/(dx*dx)
        density /= density.sum()*dx*dx
    return axis, density, rmax


def footprint(s, n=257):
    axis, density, rmax = normalized_footprint(s.distance_mm, s.output_na, n)
    return axis, density*transport(s)['output_mw'], rmax


def encircled_radius(axis, irradiance, fraction=0.9):
    x, y = np.meshgrid(axis, axis)
    radius = np.hypot(x, y).ravel()
    order = np.argsort(radius)
    powers = irradiance.ravel()[order]
    if powers.sum() <= 0:
        return 0.
    index = np.searchsorted(np.cumsum(powers), fraction*powers.sum())
    return float(radius[order[min(index, len(order)-1)]])


def mode_count(wavelength_nm=650.):
    v = 2*np.pi*(CORE_RADIUS_MM*1e6)*FIBER_NA/wavelength_nm
    return v, v*v/2


@lru_cache(maxsize=8)
def wave_basis(wavelength_nm, distance_mm, output_na, seed=4, n=161):
    """48 angular components in a local on-axis observation patch.

    This is NOT an eigenmode solver. The aperture limits local angular spread
    in the far field; the NA limits it near the tip. All axes are resolved.
    """
    rng = np.random.default_rng(seed)
    spread = min(output_na, CORE_RADIUS_MM/np.hypot(CORE_RADIUS_MM, distance_mm))
    lam = wavelength_nm/1000
    width = 12*lam/spread
    axis = np.linspace(-width/2, width/2, n)
    x, y = np.meshgrid(axis, axis)
    radii = spread*np.sqrt(rng.uniform(0, 1, 48))
    angle = rng.uniform(0, 2*np.pi, 48)
    sx, sy = radii*np.cos(angle), radii*np.sin(angle)
    k = 2*np.pi/lam
    phase = k*(sx[:, None, None]*x+sy[:, None, None]*y)
    basis = np.exp(1j*phase).reshape(48, -1).astype(np.complex64)
    internal_sin = output_na/CORE_INDEX*np.sqrt(rng.uniform(0, 1, 48))
    beta_relative = k*CORE_INDEX*(np.sqrt(1-internal_sin**2)-1)
    kz_relative = k*(np.sqrt(1-radii*radii)-1)
    return axis, basis, rng.uniform(0, 2*np.pi, 48), beta_relative, kz_relative


def interference(s, wavelength_nm=650., mode='coherent', averages=1, phase_shift=0., seed=4):
    if mode == 'halogen':
        # The exact incoherent expectation needs no complex wave basis.
        spread = min(s.output_na, CORE_RADIUS_MM/np.hypot(CORE_RADIUS_MM, s.distance_mm))
        width = 12*(float(wavelength_nm)/1000)/spread
        axis = np.linspace(-width/2, width/2, 161)
        intensity = np.ones((161, 161))
        return axis, intensity, intensity-1
    axis, basis, phase, beta, kz = wave_basis(float(wavelength_nm), s.distance_mm, s.output_na, int(seed))
    n = len(axis)
    count = int(averages) if mode == 'ensemble' else 1
    rng = np.random.default_rng(int(seed)+77)
    phases = phase + np.remainder(beta*s.length_m*1e6, 2*np.pi) + kz*phase_shift
    offsets = np.zeros((count, 48)) if count == 1 else rng.uniform(0, 2*np.pi, (count, 48))
    amplitudes = np.exp(1j*(phases[None, :]+offsets)).astype(np.complex64)/np.sqrt(48)
    fields = amplitudes@basis
    intensity = np.mean(np.abs(fields)**2, axis=0).reshape(n, n)
    # Incoherent reference is analytically 1; do not normalize individual frames.
    return axis, intensity, intensity-1


def coherence_curve(s, max_opd_um=15.):
    t = transport(s)
    wl, spectrum = t['wavelength_nm'], t['output_spectrum']
    opd = np.linspace(0, max_opd_um, 301)
    power = np.trapezoid(spectrum, wl)
    if power == 0:
        return opd, np.zeros_like(opd)
    gamma = np.trapezoid(spectrum[:, None]*np.exp(2j*np.pi*opd[None, :]*1000/wl[:, None]), wl, axis=0)/power
    return opd, np.abs(gamma)
