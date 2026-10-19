# STAGE CLASS
# the stage is a container for all the optics: arms (straight optical paths) and the
# elements in them. a Port is an element that sits in several arms and passes beams between them.
# run() traces every beam and returns what the diagnostics saw.
# possibly add a child class for transients
# - Possible "OptimizationStage" child class that lets you optimize a setup for a certain goal
# ^ This is where something like element overlap might matter
# Another Stage note -> if you use a transient stage, you'd could use a minkowski metric with rotations?

# RUN SCHEME
# every beam is traced on its own, from where it's born to where it dies.
# a source arm births a beam -> Arm.run() carries it down the arm -> a Port ends it and
# births one child beam per link -> repeat. Nothing is merged, so order never matters.
# A diagnostic hit by several beams gets several records; combining them (sum the
# powers, add the fields, ...) is post-processing, not part of tracing.

import copy
from collections import defaultdict, deque

import numpy as np

from models.diagnostic import Diagnostic
from models.elements import Element

EPS = 1e-12


def _coords(vector):
    return np.asarray(getattr(vector, "coords", vector), dtype=float)


def _unit(vector):
    vector = _coords(vector)
    norm = np.linalg.norm(vector)
    if vector.ndim != 1 or norm < EPS:
        raise ValueError("optical_axis must be a non-zero vector")
    return vector / norm


def _make_transverse_basis(axis):
    axis = _unit(axis)
    if axis.size != 3:
        raise ValueError("optical_axis must have exactly three coordinates")
    # Choose the global axis least parallel to the optical axis for stability.
    reference = np.eye(3)[np.argmin(np.abs(axis))]
    e1 = np.cross(axis, reference)
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(axis, e1)
    return np.column_stack((e1, e2))


def _rotation_between(a, b, fallback_axis):
    a, b = _unit(a), _unit(b)
    cross = np.cross(a, b)
    s = np.linalg.norm(cross)
    c = float(np.clip(np.dot(a, b), -1.0, 1.0))
    if s < EPS:
        if c > 0:
            return np.eye(3)
        # Opposite directions: choose a stable transverse axis as the rotation axis.
        axis = np.asarray(fallback_axis, dtype=float)
        axis = axis - np.dot(axis, a) * a
        axis /= np.linalg.norm(axis)
        return 2.0 * np.outer(axis, axis) - np.eye(3)
    k = cross / s
    K = np.array([[0.0, -k[2], k[1]],
                  [k[2], 0.0, -k[0]],
                  [-k[1], k[0], 0.0]])
    return np.eye(3) + s * K + (1.0 - c) * (K @ K)


# STAGE ELEMENTS
class Port(Element):
    def __init__(self, position, rotation, links, bases=None):
        super().__init__(position, rotation)
        self.links = links
        # Optional 2x2 transverse coordinate transforms keyed by (in_arm, out_arm).
        self.bases = {} if bases is None else bases

    def matrix(self, src):
        # override (and/or transmission) for a port with real optics in it, e.g. a plate's thickness
        return np.eye(2), src.n

    def arm_ids(self):
        return {arm for pair in self.links for arm in pair}

    def _frame_transform(self, in_arm, out_arm):
        # Rotate the incoming axis onto the outgoing axis, then express its transverse
        # basis vectors in the outgoing arm's local coordinates.
        Q = _rotation_between(in_arm.axis, out_arm.axis, in_arm.basis[:, 0])
        return out_arm.basis.T @ Q @ in_arm.basis

    def apply(self, src, in_arm, arms=None):
        outs = {out: p for (i, out), p in self.links.items() if i == in_arm}
        if not outs:
            raise KeyError(f"port has no links for beams arriving from arm {in_arm!r}")
        if sum(outs.values()) > 1.0 + 1e-9:
            raise ValueError(f"links from arm {in_arm!r} add up to {sum(outs.values()):.4f} > 1 (port would create power)")
        src = self.transfer(src)
        children = []
        for out, fraction in outs.items():
            basis = self.bases.get((in_arm, out))
            if basis is None and arms is not None:
                basis = self._frame_transform(arms[in_arm], arms[out])
            children.append((src.fork(fraction, basis), out))
        return children


class Arm:
    def __init__(self, optical_axis, source=None, start=0.0, transverse_basis=None):
        self.axis = _unit(optical_axis)
        if self.axis.size != 3:
            raise ValueError("optical_axis must have exactly three coordinates")
        self.basis = (_make_transverse_basis(self.axis) if transverse_basis is None
                      else np.asarray(transverse_basis, dtype=float))
        if self.basis.shape != (3, 2):
            raise ValueError("transverse_basis must have shape (3, 2)")
        if not np.allclose(self.basis.T @ self.basis, np.eye(2), atol=1e-8):
            raise ValueError("transverse_basis columns must be orthonormal")
        if not np.allclose(self.basis.T @ self.axis, np.zeros(2), atol=1e-8):
            raise ValueError("transverse_basis must be perpendicular to optical_axis")
        self.source = source  # template for a source arm, never mutated (run copies it)
        self.start = start
        self.elements = []

    def add_elem(self, elem):
        self.elements.append(elem)

    def coord(self, elem):
        return float(np.dot(_coords(elem.position), self.axis))

    def run(self, src, entry=None):
        # sorted every run so positions can change freely -> stable, so ties keep the order added
        order = sorted(((self.coord(e), e) for e in self.elements), key=lambda ce: ce[0])
        todo = [(c, e) for c, e in order if e is not entry]
        z = self.start if entry is None else self.coord(entry)
        hits = []
        for i, (c, e) in enumerate(todo):
            if c < z - EPS and not isinstance(e, Diagnostic):
                raise ValueError(f"{type(e).__name__} at z={c} overlaps the previous element " f"or sits behind where the beam enters (z={z})")
            if c > z:
                src.apply(src.translation(c - z, src.n))  # free space up to the element
            if isinstance(e, Port):
                if todo[i + 1:]:
                    raise ValueError(f"{len(todo) - i - 1} element(s) sit beyond the exit port at z={c}")
                return src, e, hits
            src = e.transfer(src)
            if isinstance(e, Diagnostic):
                hits.append((e, e.data[-1]))
            z = c + getattr(e, "thickness", 0.0)  # lens matrix already includes its thickness
        return src, None, hits


class Stage:
    def __init__(self, dispersion=None):
        self.dispersion = dispersion  # placeholder (metric / dispersion hook later)
        self.arms = {}

    def add_arm(self, arm_id, optical_axis, source=None, start=0.0, transverse_basis=None):
        self.arms[arm_id] = Arm(optical_axis, source, start, transverse_basis)
        return self.arms[arm_id]

    def add_elem(self, arm_id, elem):
        self.arms[arm_id].add_elem(elem)

    def add_port(self, port):
        missing = port.arm_ids() - self.arms.keys()
        if missing:
            raise KeyError(f"port links to arms that don't exist: {list(missing)}")
        for arm_id in port.arm_ids():
            self.arms[arm_id].add_elem(port)

    def run(self):
        for arm in self.arms.values():
            for e in arm.elements:
                if isinstance(e, Diagnostic):
                    e.data.clear()

        readings = defaultdict(list)
        # a beam is (arm_id, source, entry_port, path)
        work = deque((a_id, copy.deepcopy(a.source), None, (a_id,))
                     for a_id, a in self.arms.items() if a.source is not None)
        while work:
            arm_id, src, entry, path = work.popleft()
            if len(path) > len(self.arms):  # a beam can't visit more arms than exist unless the arms loop
                raise ValueError(f"beam went through {len(path)} arms but there are only {len(self.arms)}: loop? {path}")
            src, port, hits = self.arms[arm_id].run(src, entry)
            for diag, record in hits:
                readings[diag].append({**record, "path": path})
            if port is not None:
                for child, out in port.apply(src, arm_id, self.arms):
                    work.append((out, child, port, path + (out,)))
        return dict(readings)
