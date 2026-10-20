"""Native finite checks for the zero-clean read-only four-layer theorem.

The five-band normalization retains literal address phases and the actual
zero-address inverse. The second fixture has commuting native one-T
symplectic families whose shear products are nonzero; its quadratic
Pauli defect therefore tests the image constraint, not only its linear
degeneration. These fixtures support the analytic proof, not a search
over arbitrary circuit width.
"""
from __future__ import annotations

import itertools
import unittest

import numpy as np

from tests.test_read_only_three_layer_rigidity import (
    ALPHAS, BETAS, P_ROWS, H, S, T, _cx, _degree, _fixture,
    _linear, _pauli, _wire,
)
from tests.exact_matrices import CO, CZ, native_columns
from tests.test_read_only_query_layer_fixtures import (
    _flatten, _zero_clean_four_layers,
)


ROWS = P_ROWS + [((7, 3), (6, 4))]
PHASE_ALPHAS = ALPHAS + [(1, 2, 3)]
PHASE_BETAS = BETAS + [[(0, 2)]]


def _p(j, x):
    return _pauli(_linear(ROWS[j][0], x), _linear(ROWS[j][1], x))


def _phase(j, x):
    bits = [(x >> i) & 1 for i in range(3)]
    linear = sum(a * b for a, b in zip(PHASE_ALPHAS[j], bits))
    quadratic = sum(bits[i] * bits[k] for i, k in PHASE_BETAS[j])
    return 1j ** linear * (-1) ** quadratic


def _four_layer_fixture():
    c, d = _fixture()
    c.append(_wire(H, 0) @ _cx(1, 0) @ _wire(S, 1))
    d.append(_wire(T.conj().T, 0) @ _wire(T, 1))
    return c, d


def _fibre(c, d, x):
    result = _phase(0, x) * _p(0, x) @ c[0]
    for j in range(1, 5):
        result = (_phase(j, x) * _p(j, x) @ c[j]) @ d[j - 1] @ result
    return result


def _label(matrix):
    """Return a phase-Pauli's binary label, without discarding a residual."""
    coefficients = [np.trace(_pauli(p, q).conj().T @ matrix) / 4
                    for q in range(4) for p in range(4)]
    best = int(np.argmax(np.abs(coefficients)))
    if abs(abs(coefficients[best]) - 1) > 2e-12:
        raise AssertionError("operator is not a phase times one Pauli")
    np.testing.assert_allclose(
        matrix, coefficients[best] * _pauli(best & 3, best >> 2),
        atol=2e-12)
    return best


def _action(unitary, label):
    pauli = _pauli(label & 3, label >> 2)
    return _label(unitary @ pauli @ unitary.conj().T)


def _anf(values):
    coefficients = list(values)
    for bit in range((len(values) - 1).bit_length()):
        for x in range(len(values)):
            if x >> bit & 1:
                coefficients[x] ^= coefficients[x ^ (1 << bit)]
    return coefficients


class ReadOnlyFourLayerRigidity(unittest.TestCase):
    def test_two_quadratic_outputs_share_address_with_full_dirty_return(self):
        # The retained native words run together on distinct arbitrary
        # work wires. Shared address controls are never initialized copies.
        first = _zero_clean_four_layers(a=0, b=1, t=3, d=5)
        second = _zero_clean_four_layers(a=0, b=2, t=4, d=6)
        schedule = []
        for (kind, gates), (other_kind, other_gates) in zip(first, second):
            self.assertEqual(kind, other_kind)
            schedule.append((kind, gates + other_gates))
        layers = [gates for kind, gates in schedule if kind == 'T']
        self.assertEqual(len(layers), 4)
        self.assertEqual(sum(map(len, layers)), 16)
        self.assertEqual(sum(len(gates) for kind, gates in schedule
                             if kind == 'C'), 36)
        for layer in layers:
            targets = [gate[1] for gate in layer]
            self.assertEqual(len(targets), len(set(targets)))
            self.assertTrue(all(3 <= target < 7 for target in targets))
        word = _flatten(schedule)
        for gate in word:
            if gate[0] == 'CX':
                self.assertGreaterEqual(gate[2], 3)
            elif gate[0] == 'H':
                self.assertIn(gate[1], (3, 4))
        actual = native_columns(7, word)
        inverse = [('TDG' if gate[0] == 'T' else
                    'T' if gate[0] == 'TDG' else gate[0], *gate[1:])
                   for gate in reversed(word)]
        actual_inverse = native_columns(7, inverse)
        for column in range(128):
            a, b, c = ((column >> bit) & 1 for bit in range(3))
            row = column ^ ((a * b) << 3) ^ ((a * c) << 4)
            expected = [CO if i == row else CZ for i in range(128)]
            self.assertEqual([actual[i][column] for i in range(128)],
                             expected)
            self.assertEqual([actual_inverse[i][column] for i in range(128)],
                             expected)

    def test_literal_five_band_normalization_with_actual_inverse(self):
        c, d = _four_layer_fixture()
        u0 = _fibre(c, d, 0)
        control_gates = sum(sum(row.bit_count() for row in ur)
                            + 3 * sum(row.bit_count() for row in vr)
                            for ur, vr in ROWS)
        address_phase_gates = sum(sum(alpha) + 3 * len(beta)
                                  for alpha, beta in
                                  zip(PHASE_ALPHAS, PHASE_BETAS))
        # Eight native T/T-dagger gates, four physical layers, no clean
        # and two arbitrary work wires. Every address CZ costs H-CX-H.
        self.assertEqual((8, 4, 13 + control_gates + address_phase_gates, 0, 2),
                         (8, 4, 101, 0, 2))
        self.assertLess(max(abs(np.trace(_pauli(p, q).conj().T @ u0)) / 4
                            for p in range(4) for q in range(4)), 1 - 1e-8)
        w4 = c[4] @ d[3]
        middle = []
        for x in range(8):
            prefix_a = c[3] @ d[2]
            a = prefix_a @ _p(2, x) @ prefix_a.conj().T
            prefix_b = prefix_a @ c[2] @ d[1]
            b = prefix_b @ _p(1, x) @ prefix_b.conj().T
            r = c[0].conj().T @ _p(0, x) @ c[0]
            phase = np.prod([_phase(j, x) for j in range(5)])
            ux = _fibre(c, d, x)
            np.testing.assert_allclose(
                ux, phase * _p(4, x) @ w4 @ _p(3, x) @ a @ b
                @ w4.conj().T @ u0 @ r, atol=4e-14)
            extracted = (_p(3, x).conj().T @ w4.conj().T
                         @ _p(4, x).conj().T @ ux @ r.conj().T
                         @ u0.conj().T @ w4) / phase
            np.testing.assert_allclose(extracted, a @ b, atol=4e-14)
            middle.append((a, b))
        for x, y in itertools.product(range(8), repeat=2):
            ax, bx = middle[x]
            ay, by = middle[y]
            az, bz = middle[x ^ y]
            # These are the actual Pauli multiplication factors, not a
            # common phase optimized separately for each matrix equality.
            sign_a = (-1) ** (
                (_linear(ROWS[2][1], y) & _linear(ROWS[2][0], x))
                .bit_count() & 1)
            sign_b = (-1) ** (
                (_linear(ROWS[1][1], x) & _linear(ROWS[1][0], y))
                .bit_count() & 1)
            np.testing.assert_allclose(
                az @ bz, sign_a * sign_b * ay @ (ax @ bx)
                @ ay.conj().T @ (ay @ by), atol=5e-14)

    def test_native_auxiliary_cocycle_and_quadratic_defect_image(self):
        # N has the physical Z0,Z1 axes. The inner Clifford maps its two
        # axes to X0X1 and -Y0Y1. The complete two-axis shears commute,
        # but their products are NONZERO, exercising equation (45).
        diagonal = _wire(T, 0) @ _wire(T, 1)
        inner = _cx(0, 1) @ _wire(H, 0) @ _cx(0, 1)
        n, a, p, q = [], [], [], []
        p_labels = []
        for x in range(8):
            n.append(diagonal.conj().T
                     @ _pauli(_linear((1, 1), x), _linear((2, 4), x))
                     @ diagonal)
            a.append(inner @ diagonal
                     @ _pauli(_linear((1, 1), x), _linear((2, 0), x))
                     @ diagonal.conj().T @ inner.conj().T)
            p_labels.append(_linear((2, 4), x)
                            | (_linear((1, 3), x) << 2))
            p.append(_pauli(p_labels[-1] & 3, p_labels[-1] >> 2))
            q.append(p[-1] @ n[-1])
        # Cache exact finite symplectic maps extracted from native words.
        sn = [[_action(value, label) for label in range(16)] for value in n]
        sa = [[_action(value, label) for label in range(16)] for value in a]
        saw_nonzero_product = False
        saw_quadratic_defect = False
        for y in range(8):
            for x in range(8):
                for label in range(16):
                    self.assertEqual(sa[y][sn[x][label]], sn[x][sa[y][label]])
                    k_label = sa[y][label] ^ label
                    saw_nonzero_product |= (sn[x][k_label] ^ k_label) != 0
            c_labels = [_label(nx.conj().T @ a[y] @ nx @ a[y].conj().T)
                        for nx in n]
            for x, z in itertools.product(range(8), repeat=2):
                self.assertEqual(c_labels[x ^ z], c_labels[x] ^ c_labels[z])
                self.assertEqual(c_labels[x ^ z],
                                 sn[z][c_labels[x]] ^ c_labels[z])
            defects = []
            image_k = {sa[y][label] ^ label for label in range(16)}
            for x in range(8):
                commutator = q[x].conj().T @ a[y] @ q[x] @ a[y].conj().T
                k = _label(commutator)
                self.assertEqual(k, sn[x][sa[y][p_labels[x]] ^ p_labels[x]]
                                 ^ c_labels[x])
                product_label = _label(q[x ^ y].conj().T @ q[x] @ q[y])
                baseline = sn[x][p_labels[y]] ^ p_labels[y]
                self.assertEqual(product_label, baseline)
                defect = _label(q[x ^ y].conj().T @ a[y] @ q[x]
                                @ a[y].conj().T @ q[y])
                self.assertEqual(defect, baseline ^ sn[y][k])
                defects.append(defect)
            coefficients = _anf(defects)
            for monomial, coefficient in enumerate(coefficients):
                if monomial.bit_count() > 2:
                    self.assertEqual(coefficient, 0)
                if monomial.bit_count() == 2:
                    self.assertIn(coefficient, image_k)
                    saw_quadratic_defect |= coefficient != 0
        self.assertTrue(saw_nonzero_product)
        self.assertTrue(saw_quadratic_defect)

    def test_codimension_one_affine_derivatives_exclude_cubic_degree(self):
        # Exhaust all 256 three-variable Boolean tables and all seven
        # codimension-one direction spaces; no polynomial-degree ansatz.
        for row in range(1, 8):
            directions = [y for y in range(8) if (row & y).bit_count() % 2 == 0]
            accepted = 0
            for table in range(256):
                values = [(table >> x) & 1 for x in range(8)]
                affine_derivatives = all(
                    _degree([values[x ^ y] ^ values[x] for x in range(8)]) <= 1
                    for y in directions)
                if affine_derivatives:
                    self.assertLessEqual(_degree(values), 2)
                    accepted += 1
            self.assertEqual(accepted, 128)


if __name__ == "__main__":
    unittest.main()
