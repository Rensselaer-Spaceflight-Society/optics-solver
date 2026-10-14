"""Create a reproducible numerical check summary and default export bundle."""
from dataclasses import replace
import json
from pathlib import Path
import platform
import time
import numpy as np
from scipy.stats import qmc
from physics import Settings, transport, footprint, encircled_radius, fresnel_transmission, interference
from app import export_bundle

out = Path(__file__).parent/'verification'
out.mkdir(exist_ok=True)
s = Settings()
t = transport(s)
start = time.perf_counter()
x,e,_ = footprint(s)
first_ms = (time.perf_counter()-start)*1000
start = time.perf_counter()
footprint(s)
cached_ms = (time.perf_counter()-start)*1000
_,fine,_ = footprint(s,513)
_,coherent,_ = interference(s)
_,averaged,_ = interference(s,mode='ensemble',averages=64)
dx=x[1]-x[0]
u=qmc.Sobol(4,scramble=True,seed=7).random_base2(18)
r=.49*np.sqrt(u[:,0]); p=2*np.pi*u[:,1]
sinth=.5*np.sqrt(u[:,2]); theta=2*np.pi*u[:,3]
disp=8*sinth/np.sqrt(1-sinth*sinth)
raymoment=np.average((r*np.cos(p)+disp*np.cos(theta))**2+(r*np.sin(p)+disp*np.sin(theta))**2,weights=fresnel_transmission(sinth)**2)
ax,ee,_=footprint(replace(s,distance_mm=8),513)
xx,yy=np.meshgrid(ax,ax)
rastermoment=np.sum(ee*(xx*xx+yy*yy))/ee.sum()
report=dict(python=platform.python_version(),output_mw=t['output_mw'],peak_W_per_m2=float(e.max()*1000),
            ninety_percent_diameter_mm=2*encircled_radius(x,e),
            relative_power_conservation_error=float(abs(e.sum()*dx*dx/t['output_mw']-1)),
            relative_peak_grid_change_257_to_513=float(abs(e.max()/fine.max()-1)),
            independent_ray_second_moment_relative_error=float(abs(rastermoment/raymoment-1)),
            coherent_contrast=float(coherent.std()/coherent.mean()),
            averaged_64_contrast=float(averaged.std()/averaged.mean()),
            footprint_first_ms=first_ms,footprint_cached_ms=cached_ms,
            limitations='Numerical self-consistency checks; not experimental validation. See MODEL.md.')
(out/'numerical_checks.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
(out/'default_simulation.zip').write_bytes(export_bundle(s))
print(json.dumps(report,indent=2))
