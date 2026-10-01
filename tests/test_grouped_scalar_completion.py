"""Canonical grouped-scalar compatibility, with small native diagnostics.

The positive table has one shared rotation per SELECT and uses controlled
Z_sigma for actual inverses. The retained one-tail scalar has the same
inverse symmetry, so this is not a unique factor-two improvement. q=2 is
only an exact-word/error-interface diagnostic, not precision certification.
Large grouped operators are applied to batches, never made dense.
"""
from __future__ import annotations

import unittest

import numpy as np

from tests.test_conditional_suffix_compiler import _apply_axes, _star_column_dilation
from tests.test_joint_source_body import _body
from tests.test_one_clean_compiler import (
    _amplification_word, _encoded_signs, _literal_mask, _mask_bits,
    _paired_weights, _scalar_word,
)
from tests.test_operator_source_compiler import (
    ATOL, H, I2, X, Y, Z, _adjoint, _expand_toffolis, _operator_source,
    _pauli, _rotation, _scalar_block, _sign_word, _word_matrix,
)


COEFFICIENTS = np.array([[0, .25], [.125, 0], [0, .125], [.25, .125]])
PHASES = np.array([1, 1j, -1, -1j])


def _table_query(core, address, rows):
    """Exact finite bitwise query, keeping sigma out of temporary work."""
    masks = [_mask_bits(row) for row in rows]

    def query(component):
        word = []
        for row, pair in enumerate(masks):
            flips = [('X', bit) for j, bit in enumerate(address) if not (row >> j) & 1]
            word += flips
            for target in range(core):
                if not (pair[component] >> target) & 1:
                    continue
                if len(address) <= 2:
                    word.append((('X', 'CX', 'CCX')[len(address)], *address, target))
                else:
                    assert len(address) == 3
                    helper = (target + 1) % core
                    word += [('CCX', address[0], address[1], helper),
                             ('CCX', helper, address[2], target)] * 2
            word += flips[::-1]
        return word

    hadamards = [('H', bit) for bit in range(core)]
    return query(0) + hadamards + query(1) + hadamards


def _native_rotation(q, coefficients, address, predicate=()):
    # Core 0..2, scalar rejection sigma=3, arbitrary synthesis signal b=4.
    # The current logical bit is 5; label or branch bits occupy 6 and 7.
    assert q == 2
    sigma, signal, core = 3, 4, 3
    rows = [_encoded_signs(q, np.arccos(c)) for c in coefficients]
    mask = _table_query(core, address, rows)
    a = _scalar_word(q, _literal_mask(_paired_weights(q)[1]), signal)
    scalar = _scalar_word(q, mask, signal, predicate, helper=1)
    cy = [('SDG', sigma), ('CX', signal, sigma), ('S', sigma)]
    b = cy + [('SDG', signal)] + scalar + [('S', signal)] + cy
    primitive = a + b + _adjoint(a)
    return a, b, _amplification_word(primitive, signal), _body(a, b)


def _canonical_table():
    result = np.zeros((256, 256), dtype=complex)
    for label in range(4):
        for logical, coefficient in enumerate(COEFFICIENTS[label]):
            start = label * 64 + logical * 32
            result[start:start + 32, start:start + 32] = np.kron(
                I2, np.kron(_rotation(np.arccos(coefficient)), np.eye(8)))
    return result


class _GroupedCompletion:
    """h, mode, label, sigma, atom, logical, borrowed signal, core axes."""

    shape = (2, 2, 4, 2, 2, 2, 2, 8)
    size = int(np.prod(shape))
    data = 64  # h, logical, arbitrary signal, and every dirty-core input.

    def __init__(self, table):
        self.table = table
        self.atom = _word_matrix(2, [('CX', 0, 1), ('H', 0)])
        self.label_h = np.kron(H, H)
        self.program_calls = 0

    @staticmethod
    def direction_z(selected):
        selected = selected.copy()
        selected[2:, 1] *= -1
        return selected

    def select(self, selected, inverse=False):
        result = selected.copy()
        if inverse:
            result *= PHASES.conj().reshape(4, 1, 1, 1, 1, 1, 1)
            result[2] = _apply_axes(result[2], self.atom, (1, 2))
        else:
            result[0] = _apply_axes(result[0], self.atom, (1, 2))
        result = self.direction_z(result)
        # One table application over the full label/current-logical address.
        self.program_calls += 1
        result = _apply_axes(result, self.table.conj().T if inverse else self.table,
                             (0, 3, 4, 1, 5))
        result = self.direction_z(result)
        if inverse:
            result[0] = _apply_axes(result[0], self.atom.conj().T, (1, 2))
        else:
            result[2] = _apply_axes(result[2], self.atom.conj().T, (1, 2))
            result *= PHASES.reshape(4, 1, 1, 1, 1, 1, 1)
        return result

    def apply(self, columns, inverse=False):
        state = columns.reshape(self.shape + (-1,)).copy()
        active = _apply_axes(state[1], H, (0,))
        selected = _apply_axes(active[1], self.label_h, (0,))
        active[1] = _apply_axes(self.select(selected, inverse), self.label_h, (0,))
        state[1] = _apply_axes(active, H, (0,))
        return state.reshape(self.size, -1)

    def reflection(self, columns):
        state = columns.reshape(self.shape + (-1,)).copy()
        # Borrowed b and dirty core are EXCLUDED; sigma and all private work
        # are included. The h=0 sector has exact identity reflection.
        state[1, 0, 0, 0, 0] *= -1
        return state.reshape(self.size, -1)

    def amplify(self, columns):
        result = self.apply(columns)
        result = self.reflection(result)
        result = self.apply(result, inverse=True)
        result = self.reflection(result)
        result = self.apply(result)
        state = result.reshape(self.shape + (-1,))
        state[1] *= -1  # Literal Z_h, rather than a dropped global phase.
        return state.reshape(self.size, -1)

    def embedding(self):
        columns = np.zeros((self.size, self.data), dtype=complex)
        columns[:32, :32] = np.eye(32)
        columns[self.size // 2:self.size // 2 + 32, 32:] = np.eye(32)
        return columns

    def parity(self, columns):
        state = columns.reshape(self.shape + (-1,))
        return _apply_axes(state, np.diag([1, 1, 1, -1]), (3, 6)).reshape(self.size, -1)

    def source(self, columns, a, inverse=False):
        state = columns.reshape(self.shape + (-1,))
        return _apply_axes(state, a.conj().T if inverse else a, (6, 7)).reshape(self.size, -1)

    def sigma_z(self, columns):
        state = columns.reshape(self.shape + (-1,))
        return _apply_axes(state, Z, (3,)).reshape(self.size, -1)

    def residual(self):
        result = np.zeros((2, 2), dtype=complex)
        for label in range(4):
            base = self.atom[:2, :2] if label % 2 == 0 else I2
            term = np.diag(COEFFICIENTS[label]) @ base
            if label >= 2:
                term = term.conj().T
            result += PHASES[label] * term / 4
        return np.kron(result, np.eye(16))


class GroupedScalarCompletionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.a_word, cls.b_word, cls.inner_word, cls.body_word = _native_rotation(
            2, COEFFICIENTS.reshape(-1), [5, 6, 7])
        cls.native = _word_matrix(8, _expand_toffolis(cls.inner_word))
        cls.body = _word_matrix(8, _expand_toffolis(cls.body_word))
        cls.a_small = _word_matrix(4, _scalar_word(
            2, _literal_mask(_paired_weights(2)[1]), 3))

    def test_canonical_atoms_and_shared_actual_inverse_select(self):
        canonical = _canonical_table()
        block = _GroupedCompletion(canonical)
        np.testing.assert_allclose(block.atom, _star_column_dilation(1, 'zero')[0],
                                   atol=ATOL, rtol=0)
        # Complete SELECT columns include occupied scalar/atom flags and all
        # arbitrary signal/core inputs, not merely the accepted coefficient.
        selected = np.eye(512, dtype=complex).reshape(4, 2, 2, 2, 2, 8, 512)
        actual = block.select(selected)
        expected = selected.copy()
        for label in range(4):
            term = expected[label]
            scalar = canonical[label * 64:(label + 1) * 64,
                               label * 64:(label + 1) * 64]
            if label >= 2:
                term = _apply_axes(term, scalar.conj().T, (2, 3, 0, 4))
                if label % 2 == 0:
                    term = _apply_axes(term, block.atom.conj().T, (1, 2))
            else:
                if label % 2 == 0:
                    term = _apply_axes(term, block.atom, (1, 2))
                term = _apply_axes(term, scalar, (2, 3, 0, 4))
            expected[label] = PHASES[label] * term
        np.testing.assert_allclose(actual, expected, atol=ATOL, rtol=0)
        np.testing.assert_allclose(block.select(actual, inverse=True), selected,
                                   atol=ATOL, rtol=0)
        for coefficient in (0, .125, .25):
            rotation = _rotation(np.arccos(coefficient))
            self.assertAlmostEqual(rotation[0, 0], coefficient, delta=ATOL)
            np.testing.assert_allclose(Z @ rotation @ Z, rotation.conj().T,
                                       atol=ATOL, rtol=0)
        self.assertGreater(np.linalg.norm(_rotation(np.pi / 2) - I2, 2), 1)
        # The retained ONE-TAIL scalar has the same merge. This optimization
        # must also be applied to the comparison compiler.
        source = _operator_source(5)[0]
        for coefficient in (0, .125, .25):
            scalar = _scalar_block(source, _sign_word(round(8 * (1 - coefficient)), 4))[0]
            sigma_z = np.kron(Z, np.eye(32))
            np.testing.assert_allclose(sigma_z @ scalar @ sigma_z, scalar.conj().T,
                                       atol=ATOL, rtol=0)

    def test_native_positive_table_symmetry_and_exact_branch_controls(self):
        sigma_z = _pauli(8, {3: Z})
        np.testing.assert_allclose(sigma_z @ self.native @ sigma_z,
                                   self.native.conj().T, atol=ATOL, rtol=0)
        parity = _word_matrix(8, [('H', 3), ('CX', 4, 3), ('H', 3)])
        corrected = parity @ self.native @ parity
        canonical = _canonical_table()
        # q=2 is coarse, but this independent target comparison fixes the
        # intended +arccos(c) orientation on all borrowed-signal columns.
        positive_error = np.linalg.norm(corrected - canonical, 2)
        wrong_error = np.linalg.norm(corrected - canonical.conj().T, 2)
        self.assertLess(positive_error, wrong_error)
        a, b, word, _ = _native_rotation(2, [0, .25], [5], predicate=(6, 7))
        actual = _word_matrix(8, _expand_toffolis(word))
        inactive = [basis for basis in range(256) if (basis >> 6) & 3 != 3]
        np.testing.assert_allclose(actual[:, inactive], np.eye(256)[:, inactive],
                                   atol=ATOL, rtol=0)
        # Sigma is never a lookup/center helper. Both branch controls and
        # the current logical query are actual native gate operations.
        np.testing.assert_allclose(sigma_z @ actual @ sigma_z, actual.conj().T,
                                   atol=ATOL, rtol=0)

    def test_group_amplification_borrowed_signal_reflection_and_program_ledger(self):
        canonical = _GroupedCompletion(_canonical_table())
        native = _GroupedCompletion(self.native)
        body = _GroupedCompletion(self.body)
        embedding = canonical.embedding()
        qj = canonical.apply(embedding)
        accepted = embedding.conj().T @ qj
        expected = np.eye(64, dtype=complex)
        expected[32:, 32:] = (np.eye(32) + canonical.residual()) / 2
        np.testing.assert_allclose(accepted, expected, atol=ATOL, rtol=0)
        np.testing.assert_allclose(canonical.apply(qj, inverse=True), embedding,
                                   atol=ATOL, rtol=0)
        canonical.program_calls = native.program_calls = body.program_calls = 0
        ideal_group = canonical.amplify(embedding)
        normalized = expected.copy()
        half_block = expected[32:, 32:]
        normalized[32:, 32:] = 3 * half_block - 4 * half_block @ half_block.conj().T @ half_block
        np.testing.assert_allclose(embedding.conj().T @ ideal_group, normalized,
                                   atol=ATOL, rtol=0)
        left, _, right = np.linalg.svd(2 * half_block)
        polar_target = np.eye(64, dtype=complex)
        polar_target[32:, 32:] = left @ right
        zeta = np.linalg.norm(2 * half_block - polar_target[32:, 32:], 2)
        self.assertLessEqual(np.linalg.norm(ideal_group - embedding @ polar_target, 2),
                             4 * zeta + ATOL)
        uncorrected = native.amplify(native.parity(embedding))
        corrected = native.parity(uncorrected)
        shared = body.source(body.parity(embedding), self.a_small)
        shared = body.amplify(shared)
        shared = body.parity(body.source(shared, self.a_small, inverse=True))
        np.testing.assert_allclose(shared, corrected, atol=ATOL, rtol=0)
        self.assertEqual((canonical.program_calls, native.program_calls, body.program_calls), (3, 3, 3))
        # One table in Q, three appearances after Q R Q† R Q. A fresh
        # full-operator perturbation bound is needed for these inner rotations.
        parity = _word_matrix(8, [('H', 3), ('CX', 4, 3), ('H', 3)])
        delta = np.linalg.norm(parity @ self.native @ parity - _canonical_table(), 2)
        self.assertLessEqual(np.linalg.norm(corrected - ideal_group, 2), 3 * delta + ATOL)
        np.testing.assert_allclose(corrected.conj().T @ corrected, np.eye(64), atol=ATOL, rtol=0)
        np.testing.assert_allclose(corrected[:, :32], embedding[:, :32], atol=ATOL, rtol=0)
        # Commutant audit of the non-scalar pieces, on occupied sigma and
        # arbitrary private work. The atom/coarse words never use sigma as a
        # temporary buffer, while the diagonal reflection includes sigma.
        rng = np.random.default_rng(913)
        probe = rng.normal(size=(canonical.size, 3)) + 1j * rng.normal(size=(canonical.size, 3))
        for operation in (canonical.reflection, canonical.parity):
            np.testing.assert_allclose(operation(canonical.sigma_z(probe)),
                                       canonical.sigma_z(operation(probe)), atol=ATOL, rtol=0)
        atom = lambda x: _apply_axes(x.reshape(canonical.shape + (-1,)),
                                     canonical.atom, (4, 5)).reshape(canonical.size, -1)
        coarse_matrix = _word_matrix(1, [('T', 0), ('H', 0), ('S', 0)])
        coarse = lambda x: _apply_axes(x.reshape(canonical.shape + (-1,)),
                                       coarse_matrix, (5,)).reshape(canonical.size, -1)
        def preparations(columns):
            state = columns.reshape(canonical.shape + (-1,)).copy()
            state[1] = _apply_axes(state[1], H, (0,))
            state[1, 1] = _apply_axes(state[1, 1], canonical.label_h, (0,))
            return state.reshape(canonical.size, -1)

        def phase_and_direction(columns):
            state = columns.reshape(canonical.shape + (-1,)).copy()
            state[1, 1] = canonical.direction_z(state[1, 1])
            state[1, 1] *= PHASES.reshape(4, 1, 1, 1, 1, 1, 1)
            return state.reshape(canonical.size, -1)

        for operation in (atom, coarse, preparations, phase_and_direction):
            np.testing.assert_allclose(operation(canonical.sigma_z(probe)),
                                       canonical.sigma_z(operation(probe)), atol=ATOL, rtol=0)
        # Parse the actual native words and separate the query/center T gates
        # from source loaders. Program appearances come from the executed
        # outer circuit above, not from the number of labels in the table.
        def t_count(word):
            return sum(gate[0] in ('T', 'TDG') for gate in _expand_toffolis(word))

        loader_t_per_source = 8  # Two q=2 paired loaders, four T gates each.
        self.assertEqual(t_count(self.a_word), 3 * loader_t_per_source)
        query_t = t_count(self.b_word) - t_count(self.a_word)
        literal_t = native.program_calls * t_count(self.inner_word)
        shared_t = body.program_calls * t_count(self.body_word) + 2 * t_count(self.a_word)
        self.assertEqual(literal_t - 5 * native.program_calls * query_t,
                         135 * loader_t_per_source)
        self.assertEqual(shared_t - 5 * body.program_calls * query_t,
                         87 * loader_t_per_source)
        # Hoisting the loader and simplifying fixed masks must also be done
        # in the baseline; these appearance counts are not a leading-q gain.


if __name__ == '__main__':
    unittest.main()
