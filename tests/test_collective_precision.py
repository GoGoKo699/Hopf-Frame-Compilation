"""Bounded checks of collective precision correction, with literal phases.

History is checked exactly on every physical input. Other checks use complete
small matrices, including occupied dirty columns, and test identities rather
than fitting a synthesis. Native width and source-call bounds are analytic;
their explicit sufficient reservations are checked in integer arithmetic.
"""
from dataclasses import dataclass
from fractions import Fraction
import unittest

import numpy as np

from tests.test_operator_source_compiler import X, Y, Z


ATOL = 3e-12


@dataclass(frozen=True)
class _Quadratic:
    """Exact a + b sqrt(2) amplitudes of the geometric history network."""

    a: Fraction = Fraction(0)
    b: Fraction = Fraction(0)

    def __add__(self, other):
        return _Quadratic(self.a + other.a, self.b + other.b)

    def __neg__(self):
        return _Quadratic(-self.a, -self.b)

    def __mul__(self, other):
        return _Quadratic(self.a * other.a + 2 * self.b * other.b,
                          self.a * other.b + self.b * other.a)


_ZERO = _Quadratic()
_ONE = _Quadratic(Fraction(1))
_HALFROOT = _Quadratic(Fraction(0), Fraction(1, 2))


def _history(n, basis):
    size = 1 << n
    vector = {basis ^ (2 * size): _ONE}
    for level in range(n, 0, -1):
        # Route the heap delimiter to a, path to high x, and flip a.
        def route(z):
            if z & (2 * size):
                return z
            return ((((z >> level) & 1) << n)
                    | ((z & ((1 << level) - 1)) << (n - level))
                    | (z >> (level + 1))) ^ size

        permutation = [route(z) for z in range(4 * size)]
        if sorted(permutation) != list(range(4 * size)):
            raise AssertionError('history routing is not a permutation')
        inverse = {v: k for k, v in enumerate(permutation)}
        routed = {permutation[z]: value for z, value in vector.items()}
        output = dict(routed)
        for address in range(1 << level):
            low = address << (n - level)
            high = low + 2 * size
            u, v = routed.get(low, _ZERO), routed.get(high, _ZERO)
            # Literal Ry(-pi/4), not a half-angle convention.
            output[low] = _HALFROOT * (u + v)
            output[high] = _HALFROOT * ((-u) + v)
        vector = {inverse[z]: value for z, value in output.items() if value != _ZERO}
    return vector


def _amplitude(exponent):
    if exponent % 2:
        return _Quadratic(Fraction(0), Fraction(1, 2 ** ((exponent + 1) // 2)))
    return _Quadratic(Fraction(1, 2 ** (exponent // 2)))


def _exp(h):
    eigenvalues, eigenvectors = np.linalg.eigh(h)
    return (eigenvectors * np.exp(1j * eigenvalues)) @ eigenvectors.conj().T


def _norm(matrix):
    return np.linalg.norm(matrix, 2)


def _coin(alpha):
    t = np.tan(alpha)
    g = 1 / (1 + 1j * t)
    h = t * g
    return np.array([[1j * h, -g.conjugate()], [g, -1j * h.conjugate()]])


def _layer(n, depth, alpha):
    out = np.eye(1 << n, dtype=complex)
    for p in range(1 << depth):
        ids = [p << (n - depth), (2 * p + 1) << (n - depth - 1)]
        out[np.ix_(ids, ids)] = _coin(alpha[(1 << depth) + p])
    return out


def _frame(n, alpha):
    out = np.eye(1 << n, dtype=complex)
    for depth in range(n):
        out = _layer(n, depth, alpha) @ out
    return out


def _path(n, scale):
    base, difference = np.zeros(1 << n), np.zeros(1 << n)
    for v in range(1, 1 << n):
        base[v] = 0 if v % 5 == 0 else .61 * np.sin(1.17 * v)
        difference[v] = scale * 2. ** (v.bit_length() - 1 - n) * np.cos(.71 * v)
    return base, difference


def _velocity(n, base, direction):
    prefix = np.eye(1 << n, dtype=complex)
    velocity = np.zeros_like(prefix)
    for depth in range(n):
        local = np.zeros_like(prefix)
        for p in range(1 << depth):
            v = (1 << depth) + p
            ids = [p << (n - depth), (2 * p + 1) << (n - depth - 1)]
            local[np.ix_(ids, ids)] = -direction[v] * (
                np.cos(2 * base[v]) * X - np.sin(2 * base[v]) * Y + Z)
        velocity += prefix.conj().T @ local @ prefix
        prefix = _layer(n, depth, base) @ prefix
    return velocity


def _filter(unitary, embedding, sign=1):
    dimension = len(unitary)
    projector = embedding @ embedding.conj().T
    phase = _exp(sign * (np.pi / 3) * (projector - np.eye(dimension) / 4))
    return (phase.conj().T @ phase.conj().T @ unitary @ phase
            @ unitary.conj().T @ phase @ unitary)


class CollectivePrecisionTests(unittest.TestCase):
    def test_exact_history_extension_on_all_occupied_flag_inputs(self):
        for n in range(1, 5):
            size = 1 << n
            columns = [_history(n, j) for j in range(4 * size)]
            for i, u in enumerate(columns):
                for j in range(i + 1):
                    v = columns[j]
                    inner = _ZERO
                    for k in u.keys() & v.keys():
                        inner = inner + u[k] * v[k]
                    self.assertEqual(inner, _ONE if i == j else _ZERO)
            for x in range(size):
                birth = 0 if x == 0 else n - (x & -x).bit_length()
                expected = {2 * size + x: _amplitude(n - birth)}
                for level in range(birth + 1, n + 1):
                    expected[(1 << level) + (x >> (n - level))] = _amplitude(n - level + 1)
                self.assertEqual(columns[x], expected)
            self.assertLess((2 * size - 2) + size, 4 * size)

    def test_collective_tangent_compression_has_the_correct_relative_orientation(self):
        for n in (2, 3, 4):
            size = 1 << n
            errors = []
            for scale in (.02, .01):
                base, difference = _path(n, scale)
                s = np.eye(4 * size, dtype=complex)
                middle = s.copy()
                embedding = np.zeros((4 * size, size), dtype=complex)
                weights = np.zeros(size)
                relative_product = np.eye(size, dtype=complex)
                kappa = 0
                for level in range(1, n + 1):
                    width, weight = 1 << level, 2. ** (level - n - 1)
                    ids = np.arange(width, 2 * width)
                    prefix = _layer(level, level - 1, base).conj().T @ _frame(level, base)
                    s[np.ix_(ids, ids)] = prefix
                    h = np.zeros((width, width), dtype=complex)
                    for p in range(width // 2):
                        v = width // 2 + p
                        r = _coin(base[v]).conj().T @ _coin(base[v] + difference[v])
                        values, vectors = np.linalg.eigh((r - r.conj().T) / (2j))
                        local = (vectors * np.arcsin(values)) @ vectors.conj().T
                        h[2 * p:2 * p + 2, 2 * p:2 * p + 2] = local
                        kappa = max(kappa, _norm(local) / weight)
                    middle[np.ix_(ids, ids)] = _exp(h / weight)
                    injection = np.zeros((width, size))
                    for m in range(width):
                        x = m << (n - level)
                        injection[m, x] = 1
                        embedding[width + m, x] = np.sqrt(weight)
                        weights[x] += weight
                    generator = injection.T @ prefix.conj().T @ h @ prefix @ injection
                    relative_product = _exp(generator) @ relative_product
                embedding[2 * size:3 * size] = np.diag(np.sqrt(1 - weights))
                relative = _frame(n, base).conj().T @ _frame(n, base + difference)
                np.testing.assert_allclose(relative_product, relative, atol=ATOL, rtol=0)
                np.testing.assert_allclose(embedding.conj().T @ embedding, np.eye(size), atol=ATOL, rtol=0)
                compression = embedding.conj().T @ s.conj().T @ middle @ s @ embedding
                errors.append(_norm(compression - relative))
                self.assertLessEqual(errors[-1], kappa ** 2 + ATOL)
            self.assertGreater(errors[0] / errors[1], 3.8)
            self.assertLess(errors[0] / errors[1], 4.2)

    def test_generic_two_flag_purifier_controls_dirty_columns_and_literal_phase(self):
        embedding = np.eye(16, 4, dtype=complex)
        target = np.kron(_exp(-.19 * Z) @ _exp(.23 * X), np.eye(2))
        for delta in (.02, .01):
            flag_x = np.kron(X, np.eye(2))
            c1 = _exp(delta * np.kron(flag_x, np.kron(X, Z))) @ np.kron(np.eye(4), np.kron(_exp(.23 * X), np.eye(2)))
            c2 = _exp(.7 * delta * np.kron(flag_x, np.kron(Z, X))) @ np.kron(np.eye(4), np.kron(_exp(-.19 * Z), np.eye(2)))
            actual = c2 @ c1
            accepted = embedding.conj().T @ actual @ embedding
            left, singulars, right = np.linalg.svd(accepted)
            polar = left @ right
            positive = (right.conj().T * singulars) @ right
            filtered = _filter(actual, embedding)
            f = positive @ (np.eye(4) + np.exp(-1j * np.pi / 3) * (np.eye(4) - positive @ positive))
            np.testing.assert_allclose(embedding.conj().T @ filtered @ embedding, polar @ f, atol=ATOL, rtol=0)
            expected = max((1 - singulars) * np.sqrt(singulars + 2))
            self.assertAlmostEqual(_norm(filtered @ embedding - embedding @ polar), expected, delta=ATOL)
            self.assertLessEqual(_norm(filtered @ embedding - embedding @ target), (2 + np.sqrt(3)) * _norm(accepted - target) + ATOL)
            np.testing.assert_allclose(filtered.conj().T @ filtered, np.eye(16), atol=ATOL, rtol=0)
        # Three Pauli phases realize D literally; no uncharged scalar.
        projector = embedding @ embedding.conj().T
        za, zc = np.kron(Z, np.eye(8)), np.kron(np.eye(2), np.kron(Z, np.eye(4)))
        np.testing.assert_allclose(_exp(np.pi / 12 * (za + zc + za @ zc)),
                                   _exp(np.pi / 3 * (projector - np.eye(16) / 4)), atol=ATOL, rtol=0)

    def test_midpoint_velocity_is_cubic_in_fixed_input_coordinates(self):
        for n in (2, 3, 4):
            errors = []
            for scale in (.02, .01):
                base, difference = _path(n, scale)
                target = _frame(n, base).conj().T @ _frame(n, base + difference)
                midpoint = _exp(_velocity(n, base + difference / 2, difference))
                errors.append(_norm(target - midpoint))
                self.assertLessEqual(errors[-1], 16 * scale ** 3)
            self.assertGreater(errors[0] / errors[1], 7.8)
            self.assertLess(errors[0] / errors[1], 8.2)

    def test_opposite_filters_cancel_quadratic_phase_including_dirty_leakage(self):
        rng = np.random.default_rng(734)
        raw = rng.normal(size=(16, 16)) + 1j * rng.normal(size=(16, 16))
        generator = (raw + raw.conj().T) / 2
        generator /= _norm(generator)
        embedding = np.eye(16, 4, dtype=complex)
        errors = []
        for scale in (.02, .01):
            actual = _exp(scale * generator)
            left, _, right = np.linalg.svd(embedding.conj().T @ actual @ embedding)
            polar = left @ right
            plus, minus = _filter(actual, embedding), _filter(actual, embedding, -1)
            output = minus @ plus @ embedding
            errors.append(_norm(output - embedding @ polar @ polar))
            self.assertLessEqual(errors[-1], 11 * scale ** 3)
            self.assertGreater(_norm((plus @ plus - minus @ plus) @ embedding), 10 * errors[-1])
            np.testing.assert_allclose(plus.conj().T @ minus.conj().T @ output, embedding, atol=ATOL, rtol=0)
        self.assertGreater(errors[0] / errors[1], 7.8)
        self.assertLess(errors[0] / errors[1], 8.2)

    def test_fixed_coarse_basis_misses_the_real_second_order_edge(self):
        base = np.zeros(4)
        prefix = np.eye(4, dtype=complex)
        for depth in range(2):
            for p in range(1 << depth):
                ids = [p << (2 - depth), (2 * p + 1) << (1 - depth)]
                for pauli in (X, Y, Z):
                    local = np.zeros((4, 4), complex)
                    local[np.ix_(ids, ids)] = pauli
                    self.assertEqual((prefix.conj().T @ local @ prefix)[1, 0], 0)
            prefix = _layer(2, depth, base) @ prefix
        for epsilon in (.01, .005):
            relative = _frame(2, base).conj().T @ _frame(2, np.array([0, epsilon, epsilon, 0]))
            missing = ((relative - relative.conj().T) / 2)[1, 0]
            self.assertAlmostEqual(missing.real / epsilon ** 2, .5, delta=epsilon)

    def test_explicit_width_and_inverse_call_precision_reservations(self):
        # No floating exponentials at the large endpoint dimensions.
        for n in (16, 17, 24, 32, 64, 128):
            size = 1 << n
            s = (size + n - 1) // n
            capacity = size + n + 7
            self.assertLessEqual(s + n + 12, capacity)
            self.assertLessEqual(3 * s + n + 12, capacity)
            self.assertLessEqual(2 * s + n + 13, capacity)
            for precision in (3 * s, size):
                q = precision + 20
                self.assertLessEqual(q + 2, capacity)
                self.assertGreaterEqual(q + 1, 2 * n)
            self.assertLess(Fraction(172, 1024), 1)  # controlled midpoint error / d²
        # Each filter contains C,C†,C and D,D,D†,D†. The error budget
        # charges six middle banks and eight phase-bank appearances.
        self.assertLess(Fraction((6 + 8) * 130, 1 << 20), 1)
        self.assertLess(Fraction(12 * 130, 1 << 20), 1)  # 2+sqrt(3)+8 < 12
        self.assertLess(16 + 11 + 12 * Fraction(172, 1024) + Fraction(86, 1024), 40)


if __name__ == '__main__':
    unittest.main()
