from models.space import CoordinateVector, UnOrderedRotation
from itertools import product
import math

# Helper bounds from vertices
def bounds_from_points(points: list[CoordinateVector]):
    coords = [p.coords for p in points]
    mins = tuple(min(axis) for axis in zip(*coords))
    maxs = tuple(max(axis) for axis in zip(*coords))
    return mins, maxs

def rotated_box_bounds(center: CoordinateVector, rotation: UnOrderedRotation, *lengths):
    n = center.dim
    half = [l / 2 for l in lengths]
    # +- combinations of the 1/2 of the lengths to get the corners around a center pt
    verts = [
        center + rotation.transform(CoordinateVector(n, *(s * h for s, h in zip(signs, half))))
        for signs in product((-1, 1), repeat=n)
    ]
    return bounds_from_points(verts)

# rotated_box_bounds if it wasn't rotated
def axis_aligned_box_bounds(center: CoordinateVector, *lengths):
    half = [l / 2 for l in lengths]
    mins = tuple(c - h for c, h in zip(center.coords, half))
    maxs = tuple(c + h for c, h in zip(center.coords, half))
    return mins, maxs


def sphere_bounds(center: CoordinateVector, radius):
    mins = tuple(c - radius for c in center.coords)
    maxs = tuple(c + radius for c in center.coords)
    return mins, maxs

# okay this one took a while -> but here is a bounding box for a lens able to have both 
# concave and convex

class Hull:
    def __init__(self, bounds_map, *params):
        self.mins, self.maxs = bounds_map(*params)
        self.dim = len(self.mins)

    def overlap(self, target):
        if isinstance(target, CoordinateVector):
            return all(lo <= c <= hi for lo, c, hi in zip(self.mins, target.coords, self.maxs))
        if isinstance(target, Hull):
            return all(
                a_min <= b_max and b_min <= a_max
                for a_min, a_max, b_min, b_max in zip(self.mins, self.maxs, target.mins, target.maxs)
            )
        raise TypeError(f"Cannot check overlap with {type(target).__name__}")

# Additional class for curved hulls -> specifically the lens

import math

def _scalar_bezier(vals, t):
    """de Casteljau on plain numbers."""
    vals = list(vals)
    while len(vals) > 1:
        vals = [(1 - t) * a + t * b for a, b in zip(vals, vals[1:])]
    return vals[0]


def _scalar_range(vals):
    """Exact (min, max) of a 1-D Bezier with these control values (degree <= 3)."""
    d = len(vals) - 1
    candidates = [vals[0], vals[-1]]
    if d in (2, 3):
        a = vals[1] - vals[0]
        b = vals[2] - vals[1]
        if d == 2:
            A, B, C = 0.0, b - a, a          # derivative ∝ C + B t
        else:
            c = vals[3] - vals[2]
            A, B, C = a - 2 * b + c, 2 * (b - a), a   # derivative ∝ C + B t + A t²
        roots = []
        if abs(A) < 1e-12:
            if abs(B) > 1e-12:
                roots.append(-C / B)
        else:
            disc = B * B - 4 * A * C
            if disc >= 0:
                s = math.sqrt(disc)
                roots += [(-B + s) / (2 * A), (-B - s) / (2 * A)]
        candidates += [_scalar_bezier(vals, t) for t in roots if 0 < t < 1]
    elif d > 3:
        candidates = list(vals)  # control-point bound: conservative, not tight
    return min(candidates), max(candidates)

# lens hull class (a lot of work for one elem imo)
# need to switch between spherical and xyz essentially -> xyz is default basis

def _arc(xv, R, radius, m):
    """Exact cap profile as (x, rho) samples, vertex -> rim, stepping evenly in angle."""
    if math.isinf(R):
        return [(xv, radius * i / (m - 1)) for i in range(m)]
    tmax = math.asin(radius / abs(R))
    pts = []
    for i in range(m):
        t = tmax * i / (m - 1)
        pts.append((xv + R * (1 - math.cos(t)), abs(R) * math.sin(t)))
    return pts




# LENS CLASS + HELPERS
# creates a hull (similar to above) of a lens
# this shape is essentially a cylinder with curved endcaps
# might be too much -> consider redoing?
# one issue with all these classes is they're good for giving rendering instructions
# but not actually useful for simulation -> might be worth threadpooling

# like the offset created by having a curved endcap on a cylinder
def _sag(R, h):
    if math.isinf(R):
        return 0.0
    if abs(R) < h:
        raise ValueError(f"Radius of curvature {abs(R)} must be >= cylinder radius {h}")
    return R - math.copysign(math.sqrt(R * R - h * h), R)

def _support_local(u, w, radius, xv1, xe1, xv2, xe2, R1, R2):
    cands = [(xv1, 0.0), (xe1, radius), (xv2, 0.0), (xe2, radius)]

    norm = math.hypot(u, w)
    if norm > 0:
        for xv, R, convex in ((xv1, R1, R1 > 0), (xv2, R2, R2 < 0)):
            if math.isinf(R) or not convex:
                continue
            cx = xv + R
            px = cx + abs(R) * u / norm
            rho = abs(R) * w / norm
            on_cap = (px <= cx) if R > 0 else (px >= cx)
            if on_cap and rho <= radius:
                cands.append((px, rho))

    return max(u * x + w * rho for x, rho in cands)


def capped_cylinder_bounds(center, rotation, thickness, radius, R1, R2, thickness_at="center"):
    n = center.dim
    s1, s2 = _sag(R1, radius), _sag(R2, radius)

    if thickness_at == "center":
        xv1, xv2 = -thickness / 2, thickness / 2
        xe1, xe2 = xv1 + s1, xv2 + s2
    else:
        xe1, xe2 = -thickness / 2, thickness / 2
        xv1, xv2 = xe1 - s1, xe2 - s2

    if xv1 > xv2 or xe1 > xe2:
        raise ValueError("Caps cross")

    # direction of each local axis (columns of the rotation mtx) found via transform
    cols = [
        rotation.transform(CoordinateVector(n, *(1.0 if j == k else 0.0 for j in range(n)))).coords
        for k in range(n)
    ]

    mins, maxs = [], []
    for i in range(n):
        u = cols[0][i] # axis i, axial part
        w = math.sqrt(sum(cols[k][i] ** 2 for k in range(1, n)))  # side part
        hi = _support_local(u, w, radius, xv1, xe1, xv2, xe2, R1, R2)
        lo = -_support_local(-u, w, radius, xv1, xe1, xv2, xe2, R1, R2)
        mins.append(center.coords[i] + lo)
        maxs.append(center.coords[i] + hi)
    return tuple(mins), tuple(maxs)

class LensHull:
    def __init__(self, center, rotation, thickness, radius, R1, R2, thickness_at="center"):
        self.center, self.rotation = center, rotation
        self.n, self.radius, self.R1, self.R2 = center.dim, radius, R1, R2

        s1, s2 = _sag(R1, radius), _sag(R2, radius) # also validates |R| >= radius
        if thickness_at == "center":
            self.xv1, self.xv2 = -thickness / 2, thickness / 2
        else:
            self.xv1, self.xv2 = -thickness / 2 - s1, thickness / 2 - s2
        self.xe1, self.xe2 = self.xv1 + s1, self.xv2 + s2
        if self.xv1 > self.xv2 or self.xe1 > self.xe2:
            raise ValueError("uh oh caps overlap")

        # stage direction of each local axis
        self.axes = [
            rotation.transform(CoordinateVector(self.n, *(1.0 if j == k else 0.0 for j in range(self.n)))).coords
            for k in range(self.n)
        ]

    def front(self, rho):
        return self.xv1 + _sag(self.R1, rho)
    
    def back(self, rho):
        return self.xv2 + _sag(self.R2, rho)

    def to_local(self, point):
        off = [p - c for p, c in zip(point.coords, self.center.coords)]
        return tuple(sum(o * a for o, a in zip(off, axis)) for axis in self.axes)

    def contains(self, point, eps=1e-12):
        local = self.to_local(point)
        rho = math.hypot(*local[1:])
        if rho > self.radius + eps:
            return False
        rho = min(rho, self.radius)
        return self.front(rho) - eps <= local[0] <= self.back(rho) + eps

    def outline(self, m=32):
        top = _arc(self.xv1, self.R1, self.radius, m) + _arc(self.xv2, self.R2, self.radius, m)[::-1]
        pts = top + [(x, -rho) for x, rho in reversed(top)][1:-1]  # mirror for the other side
        out = []
        for x, s in pts:
            out.append(CoordinateVector(self.n, *(
                self.center.coords[i] + x * self.axes[0][i] + s * self.axes[1][i]
                for i in range(self.n)
            )))
        return out

# similar to what the other "bounds" functions do
    def bounds(self):
        mins, maxs = [], []
        for i in range(self.n):
            u = self.axes[0][i]
            w = math.sqrt(sum(self.axes[k][i] ** 2 for k in range(1, self.n)))
            args = (self.radius, self.xv1, self.xe1, self.xv2, self.xe2, self.R1, self.R2)
            mins.append(self.center.coords[i] - _support_local(-u, w, *args))
            maxs.append(self.center.coords[i] + _support_local(u, w, *args))
        return tuple(mins), tuple(maxs)