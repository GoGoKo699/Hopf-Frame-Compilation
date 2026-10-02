"""Bounded native audits of the conditionally clean geometric source.

The m=2..4 PREP words use a small explicit OR network, not the asymptotic
balanced prefix emitter. Every CH and CCX is expanded into literal native
Clifford+T gates. Source masks are direct diagonal phase programs, not a
QROM implementation. OAA uses the explicitly stated diagonal reflection;
these fixtures do not certify its scalable native gate/depth ledger.

Only the active h=1 input assumes a zero core and zero prefix scratch.
The small inactive fixtures include every core/helper/target/b input column.
All comparisons retain global phases, actual adjoints, and borrowed work.
"""
from __future__ import annotations

import unittest

import numpy as np

try:
    from .test_operator_source_compiler import (
        H, X, Z, _adjoint, _apply_native_word, _expand_toffolis, _word_matrix,
    )
except ImportError:
    from test_operator_source_compiler import (
        H, X, Z, _adjoint, _apply_native_word, _expand_toffolis, _word_matrix,
    )


ATOL = 5e-11
K = X @ Z


def _native_ch(control, target):
    # V_y = (SH) T (SH)^dagger, so V_y Z V_y^dagger = H exactly.
    # Its scalar cancels between V_y and its actual adjoint.
    vy = [('SDG', target), ('H', target), ('T', target),
          ('H', target), ('S', target)]
    cz = [('H', target), ('CX', control, target), ('H', target)]
    return _adjoint(vy) + cz + vy


def _native(word):
    expanded = []
    for gate in word:
        name, *qubits = gate
        if name == 'CH':
            expanded += _native_ch(*qubits)
        elif name == 'CZ':
            control, target = qubits
            expanded += [('H', target), ('CX', control, target), ('H', target)]
        elif name == 'CCZ':
            left, right, target = qubits
            expanded += [('H', target), ('CCX', left, right, target), ('H', target)]
        elif name == 'Z':
            expanded += [('S', qubits[0]), ('S', qubits[0])]
        else:
            expanded.append(gate)
    return _expand_toffolis(expanded)


def _or_xor(left, right, output):
    """output ^= left OR right, with read-only inputs and literal phases."""
    return [('CX', left, output), ('CX', right, output),
            ('CCX', left, right, output)]


def _columns(width, indices):
    indices = list(indices)
    result = np.zeros((1 << width, len(indices)), dtype=complex)
    result[indices, np.arange(len(indices))] = 1
    return result


class _ConditionalGeometricSource:
    """Finite direct-mask source; qubit zero is the low integer bit."""

    def __init__(self, m, masks=(0, 0)):
        assert 2 <= m <= 4
        assert len(masks) == 2 and all(0 <= mask < 1 << m for mask in masks)
        self.m = m
        self.core = list(range(m))
        self.coins, self.endpoint = self.core[:-1], self.core[-1]
        self.prefix = list(range(m, 2 * m - 1))
        self.scratch = [2 * m - 1] if m == 4 else []
        self.prep_width = 2 * m - 1 + len(self.scratch)
        self.helpers = self.prefix + self.scratch
        self.target, self.b, self.h = range(self.prep_width, self.prep_width + 3)
        self.width = self.prep_width + 3
        self.masks = masks

        prefix = [('CX', self.coins[0], self.prefix[0])]
        if m >= 3:
            prefix += _or_xor(*self.coins[:2], self.prefix[1])
        interval = []
        if m == 4:
            # This interval is deliberately not a prefix and is not CH-invariant.
            interval = _or_xor(self.coins[1], self.coins[2], self.scratch[0])
            prefix += interval + _or_xor(
                self.coins[0], self.scratch[0], self.prefix[2])
        self.leaky_prefix = prefix
        self.prefix_word = prefix + _adjoint(interval)
        middle = [('CH', self.prefix[j - 1], self.coins[j])
                  for j in range(1, m - 1)]
        middle += [('X', self.endpoint), ('CX', self.prefix[-1], self.endpoint)]
        self.prep = _native(
            [('H', q) for q in self.coins] + self.prefix_word
            + middle + _adjoint(self.prefix_word))
        self.bad_prep = _native(
            [('H', q) for q in self.coins] + self.leaky_prefix
            + middle + _adjoint(self.leaky_prefix))

        phase = []
        for j, q in enumerate(self.core):
            if (masks[0] >> j) & 1:
                phase.append(('CZ', self.h, q))
            if ((masks[0] ^ masks[1]) >> j) & 1:
                phase.append(('CCZ', self.h, self.b, q))
        self.phase = _native(phase)
        self.select = _native([('CCZ', self.h, self.b, self.target),
                               ('CCX', self.h, self.b, self.target)])
        # Chronological order: H_b, V, P, V^dagger, C_hb(K), H_b.
        self.word = ([('H', self.b)] + self.prep + self.phase
                     + _adjoint(self.prep) + self.select + [('H', self.b)])
        self.amplitudes = np.array(
            [2 ** (-(j + 1) / 2) for j in range(m - 1)]
            + [2 ** (-(m - 1) / 2)])
        self.coefficients = tuple(sum(
            amplitude ** 2 * (-1) ** ((mask >> j) & 1)
            for j, amplitude in enumerate(self.amplitudes)) for mask in masks)

    def apply(self, columns, inverse=False):
        return _apply_native_word(
            self.width, _adjoint(self.word) if inverse else self.word, columns)

    def embedding(self):
        return _columns(self.width, [(1 << self.h) | (target << self.target)
                                     for target in range(2)])

    def reflect(self, columns):
        # R_h = I - 2[h=1,b=0,core=0]; it does not test any helper.
        # Zero helpers form an invariant active subspace, audited below.
        # For h=0 this is identity for every core/helper/target/b input.
        result = columns.copy()
        indices = np.arange(1 << self.width)
        zero_mask = ((1 << self.m) - 1) | (1 << self.b)
        marked = ((indices & (1 << self.h)) != 0) & ((indices & zero_mask) == 0)
        result[marked] *= -1
        return result

    def amplify(self, columns, wrong_inverse=False, global_minus=False):
        result = self.apply(columns)
        result = self.reflect(result)
        result = self.apply(result, inverse=not wrong_inverse)
        result = self.reflect(result)
        result = self.apply(result)
        if global_minus:
            return -result
        # Literal Z_h supplies the active minus while preserving inactive inputs.
        return _apply_native_word(self.width, _native([('Z', self.h)]), result)


class ConditionalGeometricSourceTests(unittest.TestCase):
    def test_native_controlled_h_has_literal_phase(self):
        expected = np.eye(4, dtype=complex)
        expected[np.ix_([1, 3], [1, 3])] = H
        word = _native_ch(0, 1)
        np.testing.assert_allclose(_word_matrix(2, word), expected, atol=ATOL, rtol=0)
        np.testing.assert_allclose(_word_matrix(2, word + _adjoint(word)),
                                   np.eye(4), atol=ATOL, rtol=0)
        self.assertEqual(sum(gate[0] in ('T', 'TDG') for gate in word), 2)

    def test_native_prefix_xor_clears_internal_work(self):
        for m in range(2, 5):
            with self.subTest(m=m):
                source = _ConditionalGeometricSource(m)
                varying = len(source.coins) + len(source.prefix)
                inputs, outputs = [], []
                for assignment in range(1 << varying):
                    coin_value = assignment & ((1 << len(source.coins)) - 1)
                    g_value = assignment >> len(source.coins)
                    basis = coin_value | (g_value << m)
                    image = basis
                    for j, q in enumerate(source.prefix, start=1):
                        if coin_value & ((1 << j) - 1):
                            image ^= 1 << q
                    inputs.append(basis)
                    outputs.append(image)
                actual = _apply_native_word(source.prep_width,
                                            _native(source.prefix_word),
                                            _columns(source.prep_width, inputs))
                np.testing.assert_allclose(actual, _columns(source.prep_width, outputs),
                                           atol=ATOL, rtol=0)

    def test_native_prep_onehot_isometry_and_actual_inverse(self):
        for m in range(2, 5):
            with self.subTest(m=m):
                source = _ConditionalGeometricSource(m)
                initial = _columns(source.prep_width, [0])
                expected = np.zeros_like(initial)
                expected[[1 << j for j in range(m)], 0] = source.amplitudes
                actual = _apply_native_word(source.prep_width, source.prep, initial)
                np.testing.assert_allclose(actual, expected, atol=ATOL, rtol=0)
                restored = _apply_native_word(source.prep_width, _adjoint(source.prep), actual)
                np.testing.assert_allclose(restored, initial, atol=ATOL, rtol=0)
                wrong = _apply_native_word(source.prep_width, source.prep, actual)
                self.assertGreater(np.linalg.norm(wrong - initial), 0.5)

    def test_interval_garbage_must_be_cleared_before_controlled_h(self):
        source = _ConditionalGeometricSource(4)
        initial = _columns(source.prep_width, [0])
        actual = _apply_native_word(source.prep_width, source.bad_prep, initial)
        scratch_one = np.flatnonzero(np.arange(1 << source.prep_width)
                                     & (1 << source.scratch[0]))
        self.assertGreater(np.linalg.norm(actual[scratch_one]), 0.25)

    def test_prep_preserves_zero_helpers_for_every_core_input(self):
        for m in range(2, 5):
            with self.subTest(m=m):
                source = _ConditionalGeometricSource(m)
                initial = _columns(source.prep_width, range(1 << m))
                actual = _apply_native_word(source.prep_width, source.prep, initial)
                np.testing.assert_allclose(actual[1 << m:], 0, atol=ATOL, rtol=0)
                np.testing.assert_allclose(actual.conj().T @ actual,
                                           np.eye(1 << m), atol=ATOL, rtol=0)
                np.testing.assert_allclose(
                    _apply_native_word(source.prep_width, _adjoint(source.prep), actual),
                    initial, atol=ATOL, rtol=0)

    def test_controlled_k_and_phase_masks_on_every_basis_input(self):
        source = _ConditionalGeometricSource(3, (5, 2))
        identity = np.eye(1 << source.width, dtype=complex)
        expected_select = np.zeros_like(identity)
        expected_phase = identity.copy()
        for basis in range(1 << source.width):
            active = (basis >> source.h) & 1
            b = (basis >> source.b) & 1
            target = (basis >> source.target) & 1
            image = basis ^ ((active & b) << source.target)
            expected_select[image, basis] = (-1) ** (active & b & target)
            parity = (basis & source.masks[b]).bit_count() & 1
            expected_phase[basis, basis] = (-1) ** (active & parity)
        np.testing.assert_allclose(_apply_native_word(source.width, source.select, identity),
                                   expected_select, atol=ATOL, rtol=0)
        np.testing.assert_allclose(_apply_native_word(source.width, source.phase, identity),
                                   expected_phase, atol=ATOL, rtol=0)

    def test_source_block_and_native_adjoint_on_active_inputs(self):
        for m, masks in ((2, (0, 1)), (3, (2, 3)), (4, (4, 2)), (4, (2, 11))):
            with self.subTest(m=m, masks=masks):
                source = _ConditionalGeometricSource(m, masks)
                embedding = source.embedding()
                actual = source.apply(embedding)
                c, s = source.coefficients
                np.testing.assert_allclose(embedding.conj().T @ actual,
                                           (c * np.eye(2) + s * K) / 2,
                                           atol=ATOL, rtol=0)
                np.testing.assert_allclose(source.apply(actual, inverse=True), embedding,
                                           atol=ATOL, rtol=0)
                helper_mask = sum(1 << q for q in source.helpers)
                nonzero_helpers = (np.arange(1 << source.width) & helper_mask) != 0
                np.testing.assert_allclose(actual[nonzero_helpers], 0, atol=ATOL, rtol=0)

    def test_inactive_source_and_amplification_on_all_dirty_columns(self):
        for m in (2, 3):
            with self.subTest(m=m):
                source = _ConditionalGeometricSource(m, ((1 << m) - 2, 1))
                # h is the top bit: all lower core/helper/target/b columns occur.
                inactive = _columns(source.width, range(1 << source.h))
                np.testing.assert_allclose(source.apply(inactive), inactive, atol=ATOL, rtol=0)
                np.testing.assert_allclose(source.amplify(inactive), inactive,
                                           atol=ATOL, rtol=0)
                wrong = source.amplify(inactive[:, :1], global_minus=True)
                np.testing.assert_allclose(wrong, -inactive[:, :1], atol=ATOL, rtol=0)

    def test_oaa_accepted_block_and_full_output_error(self):
        for m, masks in ((2, (0, 1)), (3, (2, 2)), (3, (2, 1)),
                         (4, (4, 2)), (4, (2, 11))):
            with self.subTest(m=m, masks=masks):
                source = _ConditionalGeometricSource(m, masks)
                embedding = source.embedding()
                c, s = source.coefficients
                r = np.hypot(c, s)
                self.assertGreater(r, 0)
                self.assertLessEqual(r, 1 + ATOL)
                polar = (c * np.eye(2) + s * K) / r
                gamma = r * (3 - r * r) / 2
                actual = source.amplify(embedding)
                np.testing.assert_allclose(embedding.conj().T @ actual,
                                           gamma * polar, atol=ATOL, rtol=0)
                error = actual - embedding @ polar
                np.testing.assert_allclose(error.conj().T @ error,
                                           (2 - 2 * gamma) * np.eye(2),
                                           atol=ATOL, rtol=0)
                np.testing.assert_allclose(actual.conj().T @ actual,
                                           np.eye(2), atol=ATOL, rtol=0)
                if abs(r - 1) < ATOL:
                    np.testing.assert_allclose(actual, embedding @ polar,
                                               atol=ATOL, rtol=0)

    def test_oaa_requires_actual_inverse_and_literal_active_minus(self):
        source = _ConditionalGeometricSource(4, (4, 2))
        embedding = source.embedding()
        actual = source.amplify(embedding)
        wrong = source.amplify(embedding, wrong_inverse=True)
        self.assertGreater(np.linalg.norm(wrong - actual), 0.1)
        # A coherent superposition detects replacing Z_h by a global minus.
        inactive = _columns(source.width, [1 | (1 << source.b)])
        coherent = (embedding[:, :1] + inactive) / np.sqrt(2)
        expected = (actual[:, :1] + inactive) / np.sqrt(2)
        np.testing.assert_allclose(source.amplify(coherent), expected, atol=ATOL, rtol=0)
        self.assertGreater(np.linalg.norm(source.amplify(coherent, global_minus=True)
                                         - expected), 1)


if __name__ == '__main__':
    unittest.main()
