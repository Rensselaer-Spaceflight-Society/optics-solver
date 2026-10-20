from itertools import combinations
import numpy as np

# Note: I've kinda given up on generalizing N-dimensional coordinates, but with the base classes
# it's possible to write your own but it doesn't end up being physically useful

# COORDINATE CLASSES
# generic n-dimensional coordinate vector -> useful if we need more generalized coordinates
class CoordinateVector:
    def __init__(self, dim, *args):
        self.dim = dim
        self.coords = args

    def __call__(self):
        return self

    def __add__(self, other):
        if self.dim != other.dim:
            raise ValueError(f"Cannot add vectors of different dimensions ({self.dim} and {other.dim})")
        return CoordinateVector(self.dim, *(a + b for a, b in zip(self.coords, other.coords)))

    def __sub__(self, other):
            if self.dim != other.dim:
                raise ValueError(f"Cannot add vectors of different dimensions ({self.dim} and {other.dim})")
            return CoordinateVector(self.dim, *(a - b for a, b in zip(self.coords, other.coords)))

    def __eq__(self, other):
        return self.dim == other.dim and self.coords == other.coords

    def __hash__(self):
        return hash((self.dim, self.coords))

    def is_less_than_on_axis(self, other: CoordinateVector, axis: int) -> bool:
        return self.coords[axis] < other.coords[axis]

    def proj(self, AXIS: CoordinateVector):
        return np.dot(self.coords, AXIS.coords)



# BELOW IS THE MOST OVERDONE ROTATION CLASS OF ALL TIME
# - I'm thinking of adding a "Fast3Rotation" or something if this proves too slow
# - it's useful because you don't have to worry about euler angles/deal with quaternions

# UnOrderedRotation is just a rotation in n dimensions furnished by a position vector ->
# maybe look into different representations
# e.g. a higher/lower dimensional representation might be less computation-heavy if 
# it's furnished by a higher rank tensor?
# ok this is kind of group theory semantics:
# by default assume the N-dimensional UnOrderedRotation only acts upon rank 1 tensors of dimension N
# on init it creates a matrix representation of that transformation -> this may end up computationally heavy - consider a Fast3Rotation class alternative?
class UnOrderedRotation(CoordinateVector):
    def __init__(self, dim, *args):
        super().__init__(dim, *args)

        # check for the right number of parameters
        self.param_dim = self._param_count(len(args))

        # we want the generators for so(N) via the lie algebra
        # antisymmetric
        self._generators = []
        for i, j in combinations(range(self.dim), 2):
            generator = np.zeros((self.dim, self.dim))
            # antisymmetric generator "base" for N dimensions
            generator[i, j] = -1.0
            generator[j, i] = 1.0
            self._generators.append(generator)

        # get the final matrix!!!!
        self._matrix = self._mtx_rep()

    # actually create the matrix representation
    # you can think of this as e^J where J is a weighted combo of the "base" generators
    # base generators will generate infinitesimal rotations along one axis when exponentiated
    # so if you have the exponenetial of a weighted combo (i'm hesitant to say linear combo lest I get stabbed walking home at night)
    # it creates the rotation matrix we all know and (maybe don't) love
    def _mtx_rep(self):
        gens = self._generators
        params = np.asarray(self.coords, dtype=float)
        # M is essentially a linear combo of the generators
        # where the coefficients are like a weighted sum based on the taylor expansion
        M = np.zeros((self.dim, self.dim))
        for i in range(self.param_dim):
            M = M + params[i] * gens[i]
        return self._taylor_exp(M)

    # series expansion of e^M where M is a matrix -> don't worry about why it does that
    # you can exponentiate a matrix if you taylor expand it and physics says we can do that
    # by default I made the order 20 -> just wanted to slow it down -> exactness doesn't rly matter
    # this is also a neat algorithm for fast-converging taylor series (using lin alg and the s param) shout out chijan lim
    # see: https://math.stackexchange.com/questions/4030138/taylor-series-approximation-squaring-and-scaling
    def _taylor_exp(self, M, order=20):
        norm = np.linalg.norm(M,1)
        # how many times to halve M so its norm is close to 1/2
        scalsquar_count = (
            max(0, int(np.ceil(np.log2(norm / 0.5))))
            if norm > 0.5 else 0
        )
        M = M / (2 ** scalsquar_count)
        # np.eye creates a matrix with 1s along a specified diagonal
        R, term = np.eye(len(M)), np.eye(len(M))
        # get the terms
        # TIL @ is the matrix multiplication operator in python
        for i in range(1, order + 1):
            term = term @ M / i
            R = R + term

        for i in range(scalsquar_count):
            R = R @ R
        return R

    def transform(self, vec: CoordinateVector):
        coords = vec.coords
        v = np.asarray(coords, dtype=float)
        return CoordinateVector(self.dim, *(self._matrix @ v))

    def _param_count(self, m):
        if m < 1:
            raise ValueError(f"{m} parameters not allowed, need at least 1.")

        n = 2
        while n * (n - 1) // 2 < m:
            n += 1

        if n * (n - 1) // 2 != m:
            lo = (n - 1) * (n - 2) // 2
            hi = n * (n - 1) // 2

            raise ValueError(f"Invalid parameter count {m}. Valid neighboring counts are {lo} (SO({n-1})) or {hi} (SO({n})).")
        return n