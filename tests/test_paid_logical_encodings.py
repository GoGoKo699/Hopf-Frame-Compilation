"""Exact codewords and native controls; no arbitrary-angle endpoint emitter."""
from fractions import Fraction
import unittest

import numpy as np

from compiler_robust_hopf.paid_logical_encodings import (
    apply_boolean, binary_rank, decode_six, edge_pauli_witnesses,
    encode_six, encoded_pauli_word, expand_native, extract_bits,
    product_encoder, product_six_native, six_encoder_native, zeta_encoder,
)
from compiler_robust_hopf.native_coarse_fixture import (
    gate_counts, inverse_word, simulate_columns,
)
from verification.fault_tolerant.exact_arithmetic import CQ2
from tests.test_conjugate_recovery import exact_native


def pauli(size, x, z):
    p = np.zeros((size, size), dtype=complex)
    for j in range(size):
        p[j ^ x, j] = (1j)**((x & z).bit_count()) * (-1)**((j & z).bit_count())
    return p


class EncodingTests(unittest.TestCase):
    def test_complete_boolean_permutations(self):
        for m in range(1, 10):
            gates, layout = zeta_encoder(m)
            images = [apply_boolean(x, gates) for x in range(1 << m)]
            self.assertEqual(sorted(images), list(range(1 << m)))
            self.assertEqual([apply_boolean(y, inverse_word(gates)) for y in images],
                             list(range(1 << m)))
            self.assertEqual(len(layout.active), m.bit_length())

    def test_vacuum_and_one_hot_counts(self):
        for m in range(1, 65):
            gates, layout = zeta_encoder(m)
            self.assertEqual(sum(g == 'CCX' for g, _ in gates), m - m.bit_length())
            self.assertEqual(sum(g == 'CX' for g, _ in gates),
                             sum(x.bit_count() for x in range(1, m + 1)) - m.bit_length())
            for label in range(m + 1):
                value = 0 if label == 0 else 1 << (label - 1)
                image = apply_boolean(value, gates)
                self.assertEqual(extract_bits(image, layout.active), label)
                self.assertEqual(extract_bits(image, layout.remainder), 0)

    def test_boolean_span_loss(self):
        for m in range(2, 13):
            inputs = [0] + [1 << j for j in range(m)]
            rows = inputs.copy()
            ones = (1 << len(rows)) - 1
            def span():
                functions = [sum(((x >> j) & 1) << i for i, x in enumerate(rows))
                             for j in range(m)]
                return binary_rank([ones] + functions)
            self.assertEqual(span(), m + 1)
            gates, _ = zeta_encoder(m)
            for g in gates:
                old = span()
                rows[:] = [apply_boolean(x, (g,)) for x in rows]
                self.assertEqual(span(), old - (g[0] == 'CCX'))
            self.assertEqual(span(), m.bit_length() + 1)

    def test_six_formula_padding_and_crossing(self):
        gates, layout = zeta_encoder(6)
        for x in range(64):
            self.assertEqual(encode_six(x), apply_boolean(x, gates))
            self.assertEqual(decode_six(encode_six(x)), x)
        sector = [x for x in range(64) if extract_bits(encode_six(x), layout.remainder) == 0]
        self.assertEqual(sector, [0, 1, 2, 4, 8, 16, 32, 63])
        self.assertEqual(extract_bits(encode_six(63), layout.active), 7)
        self.assertEqual((encode_six(2), encode_six(3)), (2, 7))

    def test_six_native_exact_matrix(self):
        word = six_encoder_native()
        self.assertEqual(gate_counts(word)['T_count'], 15)
        self.assertEqual(len(word), 60)
        rows = exact_native([(g.lower(), *qs) for g, qs in word], 6)
        zero, one = CQ2(), CQ2(1)
        for col in range(64):
            for row in range(64):
                self.assertEqual(rows[row][col], one if row == encode_six(col) else zero)

    def test_product_branching_code_and_full_inverse(self):
        gates, layout = product_encoder((6, 6))
        support = [x | (y << 6) for x in [0, 1, 2, 4, 8, 16, 32]
                   for y in [0, 1, 2, 4, 8, 16, 32]]
        self.assertEqual(len(set(support)), 49)
        for x in support:
            if x:
                self.assertIn(x - (x & -x), support)
            y = apply_boolean(x, gates)
            self.assertEqual(extract_bits(y, layout.remainder), 0)
            self.assertNotEqual(extract_bits(y, layout.active) & 7, 7)
            self.assertNotEqual(extract_bits(y, layout.active) >> 3, 7)
        self.assertEqual(gate_counts(product_six_native(2))['T_count'], 30)
        for x in range(1 << 12):
            y = apply_boolean(x, gates)
            self.assertEqual(apply_boolean(y, inverse_word(gates)), x)

    def test_all_conjugated_single_paulis(self):
        permutation = np.array([encode_six(x) for x in range(64)])
        for bit in range(6):
            for axis, x, z in (('X', 1 << bit, 0), ('Y', 1 << bit, 1 << bit),
                               ('Z', 0, 1 << bit)):
                expected = np.zeros((64, 64), complex)
                expected[np.ix_(permutation, permutation)] = pauli(64, x, z)
                word = encoded_pauli_word(axis, bit)
                self.assertEqual(gate_counts(word)['T_count'], 0)
                np.testing.assert_allclose(simulate_columns(word, nqubits=6),
                                           expected, atol=2e-14, rtol=0)

    def test_complete_encoded_hopf_frame(self):
        rng = np.random.default_rng(10401)
        perm = np.array([encode_six(x) for x in range(64)])
        zprime = [simulate_columns(encoded_pauli_word('Z', j), nqubits=6) for j in range(6)]
        yprime = [simulate_columns(encoded_pauli_word('Y', j), nqubits=6) for j in range(6)]
        eye = np.eye(64)
        for corners in (False, True):
            target, encoded = eye.copy(), eye.astype(complex)
            edges = sorted((6 - (v & -v).bit_length(), v - (v & -v), v)
                           for v in range(1, 64))
            for depth, u, v in edges:
                theta = rng.choice([0, np.pi / 2, -np.pi, .4]) if corners else rng.uniform(-1, 1)
                c, s = np.cos(theta), np.sin(theta)
                target[[u, v], :] = np.array([[c, -s], [s, c]]) @ target[[u, v], :]
                bit = 5 - depth
                projector = eye.astype(complex)
                for j in range(6):
                    if j != bit:
                        projector = projector @ ((eye + (-1)**((u >> j) & 1) * zprime[j]) / 2)
                encoded = (eye + (c - 1)*projector - 1j*s*projector@yprime[bit]) @ encoded
            expected = np.zeros((64, 64), complex)
            expected[np.ix_(perm, perm)] = target
            np.testing.assert_allclose(encoded, expected, atol=8e-13, rtol=0)

    def test_relabelled_edge_clique(self):
        rng = np.random.default_rng(10402)
        for n in (2, 3, 4, 5):
            size = 1 << n
            for _ in range(3):
                perm = rng.permutation(size)
                edges = [(int(perm[v - (v & -v)]), int(perm[v])) for v in range(1, size)]
                witnesses = edge_pauli_witnesses(n, edges)
                self.assertEqual(len(witnesses), n)
                for j, (u, v, x, z) in enumerate(witnesses):
                    p = pauli(size, x, z)
                    h = np.zeros((size, size), complex); h[v, u], h[u, v] = 1j, -1j
                    self.assertAlmostEqual(abs(np.trace(p @ h)), 2)
                    for _, _, xx, zz in witnesses[:j]:
                        self.assertEqual(((x & zz).bit_count() + (z & xx).bit_count()) % 2, 1)

    def test_initialized_code_dimension(self):
        # Enumerate finite isotropic spans, then restrict to clean-preserving labels.
        rng = np.random.default_rng(10403)
        for total in range(2, 7):
            for clean in range(total):
                labels = [(0, 1 << j) for j in range(total)]
                for _ in range(30):
                    kind = rng.integers(3); bit = int(rng.integers(total))
                    updated = []
                    if kind == 2:
                        other = (bit + 1) % total
                    for x, z in labels:
                        if kind == 0 and ((x >> bit) & 1) != ((z >> bit) & 1):
                            x ^= 1 << bit; z ^= 1 << bit
                        elif kind == 1:
                            z ^= ((x >> bit) & 1) << bit
                        elif kind == 2:
                            x ^= ((x >> bit) & 1) << other
                            z ^= ((z >> other) & 1) << bit
                        updated.append((x, z))
                    labels = updated
                data = total - clean
                cleanmask = ((1 << clean) - 1) << data
                induced = []
                for mask in range(1 << total):
                    x = z = 0
                    for j, (xx, zz) in enumerate(labels):
                        if (mask >> j) & 1:
                            x ^= xx; z ^= zz
                    if not x & cleanmask:
                        low = (1 << data) - 1
                        induced.append((x & low) | ((z & low) << data))
                self.assertGreaterEqual(binary_rank(induced), data)

    def test_reference_columns_and_offsets(self):
        rng = np.random.default_rng(10404)
        word = six_encoder_native()
        columns = rng.normal(size=(64, 3)) + 1j*rng.normal(size=(64, 3))
        out = simulate_columns(word, columns, nqubits=6)
        np.testing.assert_allclose(simulate_columns(inverse_word(word), out, nqubits=6),
                                   columns, atol=2e-13, rtol=0)
        shifted, layout = zeta_encoder(6, 2)
        self.assertEqual(layout.active, (2, 3, 5))
        for low in range(4):
            for x in [0] + [1 << j for j in range(6)]:
                self.assertEqual(apply_boolean(low | (x << 2), shifted), low | (encode_six(x) << 2))

    def test_reject_invalid_inputs(self):
        for value in (True, 0, -1, 1.5):
            with self.assertRaises(ValueError):
                zeta_encoder(value)
        for value in (True, -1, 64):
            with self.assertRaises(ValueError):
                encode_six(value)
        for word in ((('CX', (0, 0)),), (('CCX', (0, 1)),), (('Q', (0,)),)):
            with self.assertRaises(ValueError):
                expand_native(word)
        with self.assertRaises(ValueError):
            product_encoder(())
        with self.assertRaises(ValueError):
            encoded_pauli_word('X', 6)
        with self.assertRaises(ValueError):
            edge_pauli_witnesses(3, [(0, 8)])
        with self.assertRaises(ValueError):
            extract_bits(0, (0, 0))


if __name__ == '__main__':
    unittest.main()
