"""Finite falsification checks for the read-only three-layer obstruction.

The matrix fixture is an actual six-T, three-layer native word on three
address and two arbitrary work wires. Its deliberately non-Pauli zero
fibre prevents replacing the actual inverse by a target shortcut. Exhaustive
small binary cocycles check the independent-axis step; a CNOT counterexample
checks that the same conclusion is false for generic affine symplectic maps.
"""
from __future__ import annotations

import itertools
import unittest

import numpy as np


I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Z = np.diag([1, -1]).astype(complex)
H = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2)
S = np.diag([1, 1j])
T = np.diag([1, np.exp(1j * np.pi / 4)])


def _wire(gate, wire, width=2):
    out = np.array([[1]], dtype=complex)
    for i in reversed(range(width)):
        out = np.kron(out, gate if i == wire else I2)
    return out


def _cx(control, target, width=2):
    out = np.zeros((1 << width, 1 << width), dtype=complex)
    for value in range(1 << width):
        result = value ^ (((value >> control) & 1) << target)
        out[result, value] = 1
    return out


def _pauli(p, q, width=2):
    """Canonical X^p Z^q, retaining its literal non-Hermitian phase."""
    out = np.zeros((1 << width, 1 << width), dtype=complex)
    for value in range(1 << width):
        out[value ^ p, value] = (-1) ** ((q & value).bit_count() & 1)
    return out


def _linear(rows, value):
    return sum(((row & value).bit_count() & 1) << i
               for i, row in enumerate(rows))


# Each row is a parity of the three address bits. All controls are genuine
# address-to-work CX or CZ; expanding each CZ as H-CX-H charges three gates.
P_ROWS = [((1, 2), (4, 1)), ((5, 2), (2, 1)),
          ((3, 6), (4, 1)), ((2, 5), (5, 2))]
ALPHAS = [(1, 0, 2), (0, 3, 0), (2, 0, 1), (1, 1, 0)]
BETAS = [[(0, 1)], [(1, 2)], [(0, 2)], [(0, 1), (1, 2)]]


def _p(j, x):
    return _pauli(_linear(P_ROWS[j][0], x), _linear(P_ROWS[j][1], x))


def _lambda(j, x):
    bits = [(x >> i) & 1 for i in range(3)]
    linear = sum(a * b for a, b in zip(ALPHAS[j], bits))
    quadratic = sum(bits[i] * bits[k] for i, k in BETAS[j])
    return 1j ** linear * (-1) ** quadratic


def _fixture():
    c = [_wire(H, 0) @ _cx(0, 1) @ _wire(S, 1),
         _wire(H, 1) @ _wire(S, 0),
         _cx(1, 0) @ _wire(H, 0),
         _wire(S, 1) @ _wire(H, 0) @ _cx(0, 1)]
    d = [_wire(T, 0) @ _wire(T.conj().T, 1),
         _wire(T.conj().T, 0) @ _wire(T, 1),
         _wire(T, 0) @ _wire(T, 1)]
    return c, d


def _fibre(c, d, x):
    cx = [_lambda(j, x) * _p(j, x) @ c[j] for j in range(4)]
    return cx[3] @ d[2] @ cx[2] @ d[1] @ cx[1] @ d[0] @ cx[0]


def _bands(c, d, x):
    w3 = c[3] @ d[2]
    a = w3 @ _p(2, x) @ w3.conj().T
    prefix = w3 @ c[2] @ d[1]
    b = prefix @ _p(1, x) @ prefix.conj().T
    r = c[0].conj().T @ _p(0, x) @ c[0]
    phase = np.prod([_lambda(j, x) for j in range(4)])
    return phase, _p(3, x), a, b, r


def _cocycle(j, x, y):
    ux = _linear(P_ROWS[j][0], x)
    vx = _linear(P_ROWS[j][1], x)
    uy = _linear(P_ROWS[j][0], y)
    return (-1) ** ((vx & uy).bit_count() & 1)


def _native_action(rows, address, label, width):
    p = label & ((1 << width) - 1)
    return label ^ ((_linear(rows, address) & p) << width)


def _additive(values):
    return all(values[x ^ y] == values[x] ^ values[y]
               for x in range(len(values)) for y in range(len(values)))


def _degree(values):
    coeff = list(values)
    for bit in range((len(coeff) - 1).bit_length()):
        for x in range(len(coeff)):
            if x >> bit & 1:
                coeff[x] ^= coeff[x ^ (1 << bit)]
    return max((x.bit_count() for x, value in enumerate(coeff) if value),
               default=-1)


class ReadOnlyThreeLayerRigidity(unittest.TestCase):
    def test_literal_native_normal_form_actual_inverse_and_cocycle(self):
        c, d = _fixture()
        u0 = _fibre(c, d, 0)
        # Two T/T-dagger gates in each of three physical layers; two dirty
        # work wires and no clean wire. Clifford ledger expands every CZ.
        control_gates = sum(sum(row.bit_count() for row in ur)
                            + 3 * sum(row.bit_count() for row in vr)
                            for ur, vr in P_ROWS)
        address_phase_gates = sum(sum(alpha) + 3 * len(beta)
                                  for alpha, beta in zip(ALPHAS, BETAS))
        self.assertEqual((6, 3, 10 + control_gates + address_phase_gates, 0, 2),
                         (6, 3, 75, 0, 2))
        paulies = [_pauli(p, q) for p in range(4) for q in range(4)]
        self.assertLess(max(abs(np.trace(p.conj().T @ u0)) / 4
                            for p in paulies), 1 - 1e-8)
        bands = [_bands(c, d, x) for x in range(8)]
        for x in range(8):
            phase, p, a, b, r = bands[x]
            ux = _fibre(c, d, x)
            np.testing.assert_allclose(ux, phase * p @ a @ b @ u0 @ r,
                                       atol=3e-14)
            # Retains the actual U0 dagger, which is not a target Pauli.
            extracted = (p.conj().T @ ux @ r.conj().T @ u0.conj().T) / phase
            np.testing.assert_allclose(extracted, a @ b, atol=3e-14)
            # A is exactly Clifford even though its Pauli expansion can have
            # several terms; B and U0 are not assumed to be Clifford.
            for pauli in paulies:
                image = a @ pauli @ a.conj().T
                self.assertAlmostEqual(max(abs(np.trace(q.conj().T @ image)) / 4
                                           for q in paulies), 1, places=13)
        for x in range(8):
            for y in range(8):
                ax, bx = bands[x][2:4]
                ay, by = bands[y][2:4]
                az, bz = bands[x ^ y][2:4]
                sign = _cocycle(2, y, x) * _cocycle(1, x, y)
                np.testing.assert_allclose(
                    az @ bz, sign * ay @ (ax @ bx) @ ay.conj().T @ (ay @ by),
                    atol=4e-14)

    def test_exhaustive_independent_axis_cocycles_are_additive(self):
        # Arbitrary functions m, not an assumed linear/quadratic ansatz.
        # a=3,w=1 and a=2,w=2 include dependent and inactive address rows.
        for address_bits, width, expected in [(3, 1, 120), (2, 2, 784)]:
            size = 1 << address_bits
            count = 0
            for rows in itertools.product(range(size), repeat=width):
                for tail in itertools.product(range(1 << (2 * width)),
                                              repeat=size - 1):
                    m = (0,) + tail
                    # Diagonal tests are both necessary and an efficient
                    # filter; the remaining pairs test the entire cocycle.
                    if any(_native_action(rows, x, m[x], width) != m[x]
                           for x in range(size)):
                        continue
                    if not all(m[x ^ y] == (_native_action(rows, y, m[x], width)
                                            ^ m[y])
                               for x in range(size) for y in range(size)):
                        continue
                    self.assertTrue(_additive(m))
                    count += 1
            self.assertEqual(count, expected)

    def test_generic_affine_symplectic_cocycle_can_be_quadratic(self):
        # S(y)=CX(1->3)^y1 CX(2->3)^y2, restricted to X labels.
        # Its overlapping images violate the independent physical-T axes.
        m = [0, 2, 1, 7]

        def action(y, p):
            return p ^ ((((y & 1) * (p & 1))
                         ^ (((y >> 1) & 1) * ((p >> 1) & 1))) << 2)

        for x in range(4):
            for y in range(4):
                self.assertEqual(m[x ^ y], action(y, m[x]) ^ m[y])
        self.assertFalse(_additive(m))
        self.assertEqual(_degree([(p >> 2) & 1 for p in m]), 2)

    def test_literal_address_phase_polarization_is_bilinear(self):
        phase = [np.prod([_lambda(j, x) for j in range(4)]) for x in range(8)]
        polar = [[phase[x ^ y] / (phase[x] * phase[y])
                  for y in range(8)] for x in range(8)]
        self.assertTrue(any(value in (1j, -1j) for value in phase))
        for x, y, z in itertools.product(range(8), repeat=3):
            self.assertIn(polar[x][y], (1, -1))
            self.assertEqual(polar[x ^ z][y], polar[x][y] * polar[z][y])
            self.assertEqual(polar[x][y ^ z], polar[x][y] * polar[x][z])
        # Degree-two scalar signs are genuine Clifford targets, so a claim
        # of affine scalar rigidity would already fail without any T gates.
        f = [((x & 1) * ((x >> 1) & 1)) ^ (x & 1) for x in range(8)]
        self.assertEqual(_degree(f), 2)


if __name__ == '__main__':
    unittest.main()
