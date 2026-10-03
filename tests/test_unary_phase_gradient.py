"""Small audits of a unary phase-gradient group, with explicit scope.

Karatsuba identities and source-bit permutations are exact integer checks.
The scalar guarded phase gadget and q=4,8 PREP use literal Clifford+T
words. Group actions are reduced operators on the invariant unary source
space; they are not a native complete-frame/QROM emitter. Their columns
retain every logical input and the actual source inverse, without resets.
The asymptotic count, workspace and T-depth claims require the proof.
"""
from __future__ import annotations

import unittest

import numpy as np

try:
    from .test_operator_source_compiler import (
        H, _adjoint, _apply_native_word, _expand_toffolis,
    )
except ImportError:
    from test_operator_source_compiler import (
        H, _adjoint, _apply_native_word, _expand_toffolis,
    )


ATOL = 8e-11


def _karatsuba(q):
    """(left mask, right mask, output mask), before cyclic reduction."""
    if q == 1:
        return [(1, 1, 1)]
    half = q // 2
    terms = []
    for a, b, v in _karatsuba(half):
        terms.append((a, b, v ^ (v << half)))
        terms.append((a << half, b << half,
                      (v << half) ^ (v << (2 * half))))
        terms.append((a ^ (a << half), b ^ (b << half), v << half))
    return terms


def _cyclic_terms(q):
    terms = []
    for a, b, v in _karatsuba(q):
        folded = 0
        for j in range(2 * q - 1):
            if (v >> j) & 1:
                folded ^= 1 << (j % q)
        terms.append((a, b, folded))
    return terms


def _bilinear(terms, a, x):
    answer = 0
    for left, right, output in terms:
        if (a & left).bit_count() % 2 and (x & right).bit_count() % 2:
            answer ^= output
    return answer


def _convolution(q, a, x):
    answer = 0
    for i in range(q):
        for j in range(q):
            if ((a >> i) & 1) and ((x >> j) & 1):
                answer ^= 1 << ((i + j) % q)
    return answer


def _reverse_program(q, a):
    return sum(((a >> j) & 1) << ((-j) % q) for j in range(q))


def _shift_pair(q, a, x, y=0, h=1):
    """Two multiply-XORs, then a guarded swap (all-input classical model)."""
    if h:
        terms = _cyclic_terms(q)
        y ^= _bilinear(terms, a, x)
        x ^= _bilinear(terms, _reverse_program(q, a), y)
        x, y = y, x
    return x, y


def _native_gradient(q):
    """Native product phases and reversible binary-to-unary map, q=4,8."""
    r = q.bit_length() - 1
    binary = list(range(r))
    unary = list(range(r, r + q))
    helper = r + q
    width = helper + (r == 3)
    word = [('H', bit) for bit in binary]
    # q=8: T,S,Z; q=4: S,Z. No floating rotation oracle is used.
    for bit in binary:
        for _ in range(1 << (bit + 3 - r)):
            word.append(('T', bit))
    for label, target in enumerate(unary):
        flips = [('X', bit) for bit in binary if not ((label >> bit) & 1)]
        word += flips
        if r == 2:
            word.append(('CCX', binary[0], binary[1], target))
        else:
            word += [('CCX', binary[0], binary[1], helper),
                     ('CCX', helper, binary[2], target),
                     ('CCX', binary[0], binary[1], helper)]
        word += list(reversed(flips))
    # Exactly erase the binary label from the one-hot word.
    for bit in binary:
        word += [('CX', unary[label], bit) for label in range(q)
                 if (label >> bit) & 1]
    return width, unary, _expand_toffolis(word)


def _binary_source(q, phase_errors=None):
    r = q.bit_length() - 1
    hadamards = np.array([[1.0]])
    for _ in range(r):
        hadamards = np.kron(H, hadamards)
    phase_errors = [0.0] * r if phase_errors is None else phase_errors
    phases = np.array([
        np.exp(1j * sum(((j >> bit) & 1) *
                       (2 * np.pi * (1 << bit) / q + phase_errors[bit])
                       for bit in range(r)))
        for j in range(q)
    ])
    return phases[:, None] * hadamards


def _target_matrix(state, ell, matrix):
    result = state.copy()
    for low in range(len(state)):
        if (low >> ell) & 1:
            continue
        high = low | (1 << ell)
        result[low] = matrix[0, 0] * state[low] + matrix[0, 1] * state[high]
        result[high] = matrix[1, 0] * state[low] + matrix[1, 1] * state[high]
    return result


def _source_matrix(state, matrix):
    return np.einsum('ij,ajc->aic', matrix, state)


def _reduced_stages(state, rows, q, h=1, wrong_sign=False):
    # B=SH maps the Z eigenbasis to the Y eigenbasis.
    basis = np.diag([1, 1j]) @ H
    result = state
    for ell, table in enumerate(rows):
        result = _target_matrix(result, ell, basis.conj().T)
        shifted = result.copy()
        for local in range(len(result)):
            if not h or local >> (ell + 1):
                continue
            a = table[local & ((1 << ell) - 1)]
            direction = -1 if ((local >> ell) & 1) else 1
            if wrong_sign:
                direction = -direction
            shifted[local] = np.roll(result[local], direction * a, axis=0)
        result = _target_matrix(shifted, ell, basis)
    return result


def _ideal_group(rows, q):
    dimension = 1 << len(rows)
    result = np.eye(dimension, dtype=complex)
    for ell, table in enumerate(rows):
        stage = np.eye(dimension, dtype=complex)
        for low in range(dimension):
            if ((low >> ell) & 1) or low >> (ell + 1):
                continue
            high = low | (1 << ell)
            angle = 2 * np.pi * table[low & ((1 << ell) - 1)] / q
            c, s = np.cos(angle), np.sin(angle)
            stage[np.ix_([low, high], [low, high])] = [[c, -s], [s, c]]
        result = stage @ result
    return result


def _embedding(g, q):
    state = np.zeros((1 << g, q, 1 << g), dtype=complex)
    state[:, 0, :] = np.eye(1 << g)
    return state


class UnaryPhaseGradientTests(unittest.TestCase):
    def test_karatsuba_cyclic_bilinear_identity_and_rank(self):
        for q in (1, 2, 4, 8):
            terms = _cyclic_terms(q)
            self.assertEqual(len(terms), 3 ** (q.bit_length() - 1))
            # Agreement on all basis pairs proves this bilinear identity.
            for i in range(q):
                for j in range(q):
                    self.assertEqual(_bilinear(terms, 1 << i, 1 << j),
                                     1 << ((i + j) % q))
            for a in (0, 1, (1 << q) - 1, (1 << q) // 3):
                for x in range(1 << q):
                    self.assertEqual(_bilinear(terms, a, x),
                                     _convolution(q, a, x))

    def test_native_guarded_trilinear_phase_and_inactive_arbitrary_work(self):
        # Controls h,a,b,c; work t=a*b, root=t*c. 28 T gates total.
        compute = _expand_toffolis([('CCX', 1, 2, 4), ('CCX', 4, 3, 5)])
        central = [('H', 5), ('CX', 0, 5), ('H', 5)]
        word = compute + central + _adjoint(compute)
        self.assertEqual(sum(g[0] in ('T', 'TDG') for g in word), 28)
        inactive = [j for j in range(64) if not (j & 1)]
        active = [1 | (j << 1) for j in range(8)]
        columns = np.eye(64, dtype=complex)[:, inactive + active]
        expected = columns.copy()
        expected[:, -1] *= -1
        np.testing.assert_allclose(_apply_native_word(6, word, columns),
                                   expected, atol=ATOL, rtol=0)
        # Without the original h, arbitrary inactive inputs are not protected.
        unguarded = compute + [('S', 5), ('S', 5)] + _adjoint(compute)
        self.assertGreater(np.linalg.norm(
            _apply_native_word(6, unguarded, columns) - expected, 2), 1.9)

    def test_in_place_shift_returns_output_for_every_source_bitstring(self):
        q = 4
        for shift in range(q):
            for sign in (-1, 1):
                a = 1 << ((sign * shift) % q)
                for x in range(1 << q):
                    self.assertEqual(_shift_pair(q, a, x),
                                     (_convolution(q, a, x), 0))
        # No validity or clean-output promise is needed on the h=0 sector.
        for a in range(1 << q):
            for x in range(1 << q):
                for y in range(1 << q):
                    self.assertEqual(_shift_pair(q, a, x, y, h=0), (x, y))
        # The reverse program is essential: two forward convolutions fail.
        a, x = 2, 1
        y = _convolution(q, a, x)
        self.assertNotEqual(x ^ _convolution(q, a, y), 0)

    def test_native_product_phases_unary_decode_and_actual_inverse(self):
        for q in (4, 8):
            width, unary, word = _native_gradient(q)
            initial = np.zeros((1 << width, 1), dtype=complex)
            initial[0, 0] = 1
            prepared = _apply_native_word(width, word, initial)
            expected = np.zeros_like(initial)
            for j, bit in enumerate(unary):
                expected[1 << bit, 0] = np.exp(2j * np.pi * j / q) / np.sqrt(q)
            np.testing.assert_allclose(prepared, expected, atol=ATOL, rtol=0)
            np.testing.assert_allclose(
                _apply_native_word(width, _adjoint(word), prepared),
                initial, atol=ATOL, rtol=0)

    def test_unequal_three_stage_rows_give_literal_rounded_frame(self):
        q, rows = 8, [[1], [3, 6], [2, 7, 1, 5]]
        source = _binary_source(q)
        initial = _embedding(3, q)
        prepared = _source_matrix(initial, source)
        acted = _reduced_stages(prepared, rows, q)
        result = _source_matrix(acted, source.conj().T)
        expected = initial.copy()
        expected[:, 0, :] = _ideal_group(rows, q)
        np.testing.assert_allclose(result, expected, atol=ATOL, rtol=0)
        wrong = _source_matrix(_reduced_stages(prepared, rows, q, wrong_sign=True),
                               source.conj().T)
        self.assertGreater(np.linalg.norm((wrong - expected).reshape(64, 8), 2), 1)
        # Every source/input column, not merely the prepared eigenstate.
        arbitrary = np.eye(64, dtype=complex).reshape(8, 8, 64)
        np.testing.assert_allclose(_reduced_stages(arbitrary, rows, q, h=0),
                                   arbitrary, atol=ATOL, rtol=0)

    def test_approximate_source_is_charged_once_per_group_without_reset(self):
        q, rows = 8, [[1], [3, 6], [2, 7, 1, 5]]
        exact, actual = _binary_source(q), _binary_source(q, [.019, -.027, .013])
        delta = np.linalg.norm(actual[:, 0] - exact[:, 0])
        initial = _embedding(3, q)
        result = _source_matrix(
            _reduced_stages(_source_matrix(initial, actual), rows, q),
            actual.conj().T)
        expected = initial.copy()
        expected[:, 0, :] = _ideal_group(rows, q)
        error = np.linalg.norm((result - expected).reshape(64, 8), 2)
        self.assertLessEqual(error, 2 * delta + ATOL)
        self.assertGreater(np.linalg.norm(result[:, 1:, :]), 0.01)
        # For identity rows actual PREP† gives exact return; the ideal inverse
        # incorrectly leaves the preparation error in the source register.
        identity_rows = [[0], [0, 0], [0, 0, 0, 0]]
        acted = _reduced_stages(_source_matrix(initial, actual), identity_rows, q)
        np.testing.assert_allclose(_source_matrix(acted, actual.conj().T), initial,
                                   atol=ATOL, rtol=0)
        wrong = _source_matrix(acted, exact.conj().T)
        self.assertAlmostEqual(np.linalg.norm((wrong - initial).reshape(64, 8), 2),
                               delta, places=10)


if __name__ == '__main__':
    unittest.main()
