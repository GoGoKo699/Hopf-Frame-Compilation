"""Finite audits of the conditional-suffix grouping interface.

The scalar source below is the explicit native word from the operator-source
tests. Vectorized application keeps every dirty input column while avoiding
large dense circuit matrices. These fixtures check block composition,
failure flags, coherent inactive sectors, star forests, and sparse support;
they do not prove asymptotic gate counts or a complete synthesis theorem.
The m=3 native fixtures exercise the general scalar interface rather than
the theorem's target precision. A separate reconstruction fixture checks
the small-residual coefficient range.
"""
from __future__ import annotations

import unittest
import math

import numpy as np

try:
    from .test_operator_source_compiler import (
        H, _operator_source, _scalar_block, _sign_word, _word_matrix,
    )
except ImportError:
    from test_operator_source_compiler import (
        H, _operator_source, _scalar_block, _sign_word, _word_matrix,
    )


ATOL = 4e-11


def _apply_axes(state, matrix, axes):
    """Apply an ordinary matrix on named tensor axes, preserving batch axes."""
    axes = tuple(axes)
    order = axes + tuple(axis for axis in range(state.ndim) if axis not in axes)
    arranged = np.transpose(state, order)
    shape = arranged.shape
    result = matrix @ arranged.reshape(matrix.shape[1], -1)
    return np.transpose(result.reshape(shape), np.argsort(order))


def _matrix_unit_dilation(u, v):
    """Flag,target ordering; predicate first, then the endpoint XOR."""
    matrix = np.zeros((4, 4), dtype=complex)
    for flag in range(2):
        for target in range(2):
            new_flag = flag ^ (target != v)
            new_target = target ^ (u ^ v)
            matrix[2 * new_flag + new_target, 2 * flag + target] = 1
    return matrix


class _SmallResidualBlock:
    """Q encodes (I+E)/2 using a mode, four labels, and separate flags.

    Axes are mode, label, scalar flag, matrix-unit flag, local target, dirty
    core, and an arbitrary column batch. The five first scratch bits can
    model a conditionally initialized logical suffix. No state is prepared
    on the three-qubit dirty core.
    """

    def __init__(self, atoms):
        assert len(atoms) == 4
        self.atoms = atoms
        self.core_size = 8
        self.data_size = 2 * self.core_size
        self.scratch_size = 2 * 4 * 2 * 2
        self.size = self.scratch_size * self.data_size
        self.shape = (2, 4, 2, 2, 2, self.core_size)
        source = _operator_source(3)[0]
        self.scalars = []
        for coefficient, _, _, _ in atoms:
            k = int(round(2 * (1 - coefficient)))
            assert abs(coefficient - (1 - k / 2)) < ATOL
            self.scalars.append(_scalar_block(source, _sign_word(k, 2))[0])
        self.label_h = np.kron(H, H)

    def apply(self, columns, inverse=False):
        state = columns.reshape(self.shape + (-1,)).copy()
        state = _apply_axes(state, H, (0,))
        selected = _apply_axes(state[1], self.label_h, (0,))
        for label, (_, phase, u, v) in enumerate(self.atoms):
            term = selected[label]
            scalar = self.scalars[label]
            matrix_unit = _matrix_unit_dilation(u, v)
            if inverse:
                term = np.conj(phase) * _apply_axes(
                    term, matrix_unit.conj().T, (1, 2))
                term = _apply_axes(term, scalar.conj().T, (0, 3))
            else:
                term = _apply_axes(term, scalar, (0, 3))
                term = phase * _apply_axes(term, matrix_unit, (1, 2))
            selected[label] = term
        state[1] = _apply_axes(selected, self.label_h, (0,))
        state = _apply_axes(state, H, (0,))
        return state.reshape(self.size, -1)

    def embedding(self):
        result = np.zeros((self.size, self.data_size), dtype=complex)
        result[:self.data_size] = np.eye(self.data_size)
        return result

    def amplify(self, columns):
        # -Q R Q^dagger R Q with R = I - 2 J J^dagger.
        result = self.apply(columns)
        result[:self.data_size] *= -1
        result = self.apply(result, inverse=True)
        result[:self.data_size] *= -1
        return -self.apply(result)

    def expected_residual(self):
        residual = np.zeros((2, 2), dtype=complex)
        for coefficient, phase, u, v in self.atoms:
            residual[u, v] += coefficient * phase / 4
        return np.kron(residual, np.eye(self.core_size))


def _small_frame(height, word_factory):
    """Complete binary frame in addressed shallow-to-deep gate order."""
    size = 1 << height
    frame = np.eye(size, dtype=complex)
    for depth in range(height):
        stride = 1 << (height - depth - 1)
        for prefix in range(1 << depth):
            rows = [prefix * (stride << 1), prefix * (stride << 1) + stride]
            frame[rows] = word_factory(depth, prefix) @ frame[rows]
    return frame


def _column_support(marker, height):
    if marker == 0:
        return set(range(1 << height))
    lowbit = marker & -marker
    subtree_size = 2 * lowbit
    start = marker - marker % subtree_size
    return set(range(start, start + subtree_size))


def _group_ledger(n, ratio=64, base=256):
    """An illustrative fixed partition, not a synthesis-workspace certificate."""
    r = min(n, base)
    groups = []
    while r < n:
        s = min(n - r, r // ratio)
        # Exact ceil(log2(4*((2*s-1)*2**s+2))), without a large power of two.
        label_bits = 4 if s == 1 else s + (2 * s - 2).bit_length() + 2
        groups.append((r, s, label_bits))
        r += s
    return groups


def _star_column_dilation(height, kind, depth=0):
    """One flag: test the ancestor marker, then spread it in its subtree."""
    size = 1 << height
    if kind == 'diagonal':
        return np.eye(2 * size), np.ones(size, dtype=bool)
    varying = height if kind == 'zero' else height - depth
    subtree_size = 1 << varying
    marker = 0 if kind == 'zero' else subtree_size // 2
    hadamards = np.ones((1, 1))
    for _ in range(varying):
        hadamards = np.kron(hadamards, H)
    mixing = np.kron(np.eye(size // subtree_size), hadamards)
    displacement = np.eye(size)[:, np.arange(size) ^ marker]
    predicate = np.zeros((2 * size, 2 * size))
    valid = np.zeros(size, dtype=bool)
    for logical in range(size):
        local = logical % subtree_size
        valid[logical] = local != 0 and local != marker
        for flag in range(2):
            output_flag = flag ^ (local != marker)
            predicate[output_flag * size + logical, flag * size + logical] = 1
    dilation = np.kron(np.eye(2), mixing @ displacement) @ predicate
    return dilation, valid


def _star_labels(height):
    labels = [('diagonal', 0, 0, phase) for phase in range(4)]
    for kind, depth in [('zero', 0)] + [('node', j) for j in range(height)]:
        for orientation in (0, 1):
            labels.extend((kind, depth, orientation, phase) for phase in range(4))
    count = 1 << (len(labels) - 1).bit_length()
    return labels + [('padding', 0, 0, 0)] * (count - len(labels))


class _SmallStarBlock:
    """Actual column-star/adjoint SELECT, with a native m=3 dirty source.

    The coefficient lookup is addressed by the current logical register.
    The test directly represents its coherent table action; it does not
    assert a native query gate count or the final asymptotic synthesis cost.
    """

    def __init__(self, height):
        self.height = height
        self.logical_size = 1 << height
        self.core_size = 8
        self.data_size = self.logical_size * self.core_size
        self.labels = _star_labels(height)
        self.count = len(self.labels)
        self.shape = (2, self.count, 2, 2, self.logical_size, self.core_size)
        self.size = math.prod(self.shape)
        source = _operator_source(3)[0]
        self.scalars = {
            c: _scalar_block(source, _sign_word(int(2 * (1 - c)), 2))[0]
            for c in (0, 0.5, 1)
        }
        self.dilations = {}
        self.coefficients = []
        for index, (kind, depth, _, _) in enumerate(self.labels):
            if kind == 'padding':
                values = np.zeros(self.logical_size)
            else:
                dilation, valid = _star_column_dilation(height, kind, depth)
                self.dilations[kind, depth] = dilation
                # Vary coefficients with both logical address and label.
                values = np.array([0.5 if (z + index) % 3 else 0
                                   for z in range(self.logical_size)]) * valid
            self.coefficients.append(values)
        self.label_h = np.ones((1, 1))
        for _ in range(self.count.bit_length() - 1):
            self.label_h = np.kron(self.label_h, H)

    def _scalar(self, state, coefficients, inverse=False):
        result = state.copy()
        for logical, coefficient in enumerate(coefficients):
            circuit = self.scalars[coefficient]
            if inverse:
                circuit = circuit.conj().T
            result[:, :, logical] = _apply_axes(
                result[:, :, logical], circuit, (0, 2))
        return result

    def term(self, state, index, inverse=False):
        kind, depth, orientation, phase = self.labels[index]
        coefficient = self.coefficients[index]
        if kind == 'padding':
            # A zero scalar block encodes padding; identity would be wrong.
            return self._scalar(state, coefficient, inverse)
        dilation = self.dilations[kind, depth]
        # T = S_c D. Reverse stars use the actual T^dagger; the external
        # phase i**phase is retained, rather than inadvertently conjugated.
        adjoint = bool(orientation) ^ bool(inverse)
        if adjoint:
            state = self._scalar(state, coefficient, inverse=True)
            state = _apply_axes(state, dilation.conj().T, (1, 2))
        else:
            state = _apply_axes(state, dilation, (1, 2))
            state = self._scalar(state, coefficient)
        return (1j ** (-phase if inverse else phase)) * state

    def apply(self, columns, inverse=False):
        state = columns.reshape(self.shape + (-1,)).copy()
        state = _apply_axes(state, H, (0,))
        selected = _apply_axes(state[1], self.label_h, (0,))
        for index in range(self.count):
            selected[index] = self.term(selected[index], index, inverse)
        state[1] = _apply_axes(selected, self.label_h, (0,))
        return _apply_axes(state, H, (0,)).reshape(self.size, -1)

    def residual(self, perturbation=None):
        result = np.zeros((self.logical_size, self.logical_size), dtype=complex)
        for index, (kind, depth, orientation, phase) in enumerate(self.labels):
            if kind == 'padding':
                continue
            coefficients = self.coefficients[index].copy()
            if perturbation is not None:
                coefficients += perturbation[index]
            accepted = self.dilations[kind, depth][:self.logical_size,
                                                        :self.logical_size]
            term = np.diag(coefficients) @ accepted
            if orientation:
                term = term.conj().T
            result += 1j ** phase * term / self.count
        return np.kron(result, np.eye(self.core_size))

    def embedding(self):
        result = np.zeros((self.size, self.data_size), dtype=complex)
        result[:self.data_size] = np.eye(self.data_size)
        return result

    def amplify(self, columns):
        result = self.apply(columns)
        result[:self.data_size] *= -1
        result = self.apply(result, inverse=True)
        result[:self.data_size] *= -1
        return -self.apply(result)


def _iterated_group_ledger(n, ratio=8, base=64):
    """Avoid constructing 2**n or even an overshooting candidate group size."""
    r = min(n, base)
    groups = []
    while r < n:
        remaining = n - r
        exponent = r // ratio
        s = remaining if exponent >= remaining.bit_length() else 1 << exponent
        label_bits = (8 * s + 11).bit_length()
        groups.append((r, s, label_bits))
        r += s
    return groups


def _coarse_interpreter_matrix(programs, streamed):
    """Literal two-bit symbols; low bits are target, address, and h."""
    slots = len(programs[0])
    word_bits = 2 if streamed else 2 * slots
    width = 3 + word_bits
    dimension = 1 << width
    unitary = np.eye(dimension, dtype=complex)
    alphabet = [np.eye(2), H, _word_matrix(1, [('T', 0)]),
                _word_matrix(1, [('TDG', 0)])]

    def load(matrix, slot=None):
        images = []
        for basis in range(dimension):
            address, active = (basis >> 1) & 1, (basis >> 2) & 1
            symbol = (programs[address][slot] if streamed else
                      sum(code << (2 * j) for j, code in enumerate(programs[address])))
            images.append(basis ^ ((symbol << 3) if active else 0))
        return matrix[np.argsort(images)]

    def interpret(matrix, slot):
        offset = 3 if streamed else 3 + 2 * slot
        result = matrix.copy()
        for low in range(0, dimension, 2):
            if not ((low >> 2) & 1):
                continue
            code = (low >> offset) & 3
            result[[low, low + 1]] = alphabet[code] @ result[[low, low + 1]]
        return result

    if not streamed:
        unitary = load(unitary)
    for slot in range(slots):
        if streamed:
            unitary = load(unitary, slot)
        unitary = interpret(unitary, slot)
        if streamed:
            unitary = load(unitary, slot)
    return unitary if streamed else load(unitary)


class ConditionalSuffixCompilerTests(unittest.TestCase):
    def test_controlled_t_dirty_phase_echo_and_helper_return(self):
        # Qubit 0 stores the preserved predicate F; 1 and 2 are arbitrary
        # dirty helpers z,w. CCZ is expanded as H-CCX-H, so this is the
        # chronological X/CX/CCX/H/T/S word, including its literal phase.
        word = [
            ('CX', 0, 1), ('T', 1), ('CX', 0, 1), ('TDG', 1),
            ('CCX', 0, 1, 2), ('S', 2), ('CCX', 0, 1, 2), ('SDG', 2),
            ('H', 2), ('CCX', 0, 1, 2), ('H', 2),
        ]
        actual = _word_matrix(3, word)
        omega = np.exp(1j * np.pi / 4)
        expected = np.diag([omega ** (basis & 1) for basis in range(8)])
        np.testing.assert_allclose(actual, expected, atol=ATOL, rtol=0)
        for predicate in (0, 1):
            for z in (0, 1):
                for w in (0, 1):
                    phase = (omega ** (predicate - 2 * z * predicate)
                             * 1j ** (z * predicate - 2 * w * z * predicate)
                             * (-1) ** (w * z * predicate))
                    self.assertAlmostEqual(phase, omega ** predicate, delta=ATOL)
        # The full diagonal equality checks every helper input, and hence
        # coherent helper superpositions and external entanglement as well.

    def test_matrix_unit_predicate_and_actual_inverse(self):
        for u in range(2):
            for v in range(2):
                dilation = _matrix_unit_dilation(u, v)
                expected = np.zeros((2, 2))
                expected[u, v] = 1
                np.testing.assert_allclose(dilation[:2, :2], expected, atol=0)
                np.testing.assert_allclose(
                    dilation.conj().T @ dilation, np.eye(4), atol=0)

    def test_native_scalar_and_separate_failure_flags(self):
        source = _operator_source(3)[0]
        scalar = _scalar_block(source, _sign_word(1, 2))[0]
        core_size = source.shape[0]
        data_size = 2 * core_size
        embedding = np.zeros((4 * data_size, data_size), dtype=complex)
        embedding[:data_size] = np.eye(data_size)
        state = embedding.reshape(2, 2, 2, core_size, data_size)
        state = _apply_axes(state, scalar, (0, 3))
        state = _apply_axes(state, _matrix_unit_dilation(1, 0), (1, 2))
        atom = np.array([[0, 0], [1, 0]])
        expected = 0.5 * np.kron(atom, np.eye(core_size))
        np.testing.assert_allclose(state[0, 0].reshape(data_size, data_size),
                                   expected, atol=ATOL, rtol=0)

        # Reusing one failure flag allows a rejected scalar path to return.
        shared = np.zeros((2 * data_size, data_size), dtype=complex)
        shared[:data_size] = np.eye(data_size)
        shared = shared.reshape(2, 2, core_size, data_size)
        shared = _apply_axes(shared, scalar, (0, 2))
        shared = _apply_axes(shared, _matrix_unit_dilation(1, 0), (0, 1))
        discrepancy = shared[0].reshape(data_size, data_size) - expected
        self.assertGreater(np.linalg.norm(discrepancy, ord=2), 0.5)

    def test_padded_lcu_preserves_literal_complex_phases(self):
        atoms = [(0.5, 1, 0, 0), (1, -1, 0, 1),
                 (0.5, 1j, 1, 0), (0, -1j, 1, 1)]
        block = _SmallResidualBlock(atoms)
        embedding = block.embedding()
        actual = block.apply(embedding)
        expected = (np.eye(block.data_size) + block.expected_residual()) / 2
        np.testing.assert_allclose(actual[:block.data_size], expected,
                                   atol=ATOL, rtol=0)
        np.testing.assert_allclose(block.apply(actual, inverse=True), embedding,
                                   atol=ATOL, rtol=0)

    def test_oblivious_amplification_includes_every_dirty_input(self):
        atoms = [(0, 1, 0, 0), (0.5, -1, 0, 1),
                 (0.5, 1, 1, 0), (0, 1, 1, 1)]
        block = _SmallResidualBlock(atoms)
        embedding = block.embedding()
        residual = block.expected_residual()
        approximate = np.eye(block.data_size) + residual
        target = approximate / np.sqrt(1 + (1 / 8) ** 2)
        np.testing.assert_allclose(target.conj().T @ target,
                                   np.eye(block.data_size), atol=ATOL, rtol=0)
        b = block.apply(embedding)[:block.data_size]
        amplified = block.amplify(embedding)
        np.testing.assert_allclose(amplified[:block.data_size],
                                   3 * b - 4 * b @ b.conj().T @ b,
                                   atol=ATOL, rtol=0)
        zeta = np.linalg.norm(2 * b - target, ord=2)
        full_error = np.linalg.norm(amplified - embedding @ target, ord=2)
        self.assertLessEqual(full_error, 4 * zeta + ATOL)
        self.assertGreater(np.linalg.norm(amplified[block.data_size:]), 1e-3)

    def test_conditional_suffix_uncompute_keeps_leakage_in_full_error(self):
        atoms = [(0, 1, 0, 0), (0.5, -1, 0, 1),
                 (0.5, 1, 1, 0), (0, 1, 1, 1)]
        block = _SmallResidualBlock(atoms)
        data = block.data_size
        # The external flag h starts zero; every suffix/core input column
        # is included. Compute h=[suffix=0], control A on h, then uncompute.
        columns = np.zeros((2, block.size, block.size), dtype=complex)
        columns[0] = np.eye(block.size)
        columns[:, :data] = columns[::-1, :data]
        columns[1] = block.amplify(columns[1])
        columns[:, :data] = columns[::-1, :data]
        np.testing.assert_allclose(columns[0, :, data:],
                                   np.eye(block.size)[:, data:], atol=0, rtol=0)
        np.testing.assert_allclose(columns[1, :, data:], 0, atol=0, rtol=0)

        target = (np.eye(data) + block.expected_residual()) / np.sqrt(65 / 64)
        ideal = np.zeros_like(columns)
        ideal[0] = np.eye(block.size)
        ideal[0, :data, :data] = target
        difference = (columns - ideal).reshape(2 * block.size, block.size)
        # All nonzero error columns are active; their spectral norm covers
        # arbitrary superpositions and reference-entangled dirty inputs.
        full_error = np.linalg.norm(difference[:, :data], ord=2)
        local_error = np.linalg.norm(
            block.amplify(block.embedding()) - block.embedding() @ target, ord=2)
        self.assertAlmostEqual(full_error, local_error, delta=ATOL)
        self.assertGreater(np.linalg.norm(columns[1]), 1e-3)

    def test_sparse_group_support_with_exact_complex_coarse_words(self):
        native_words = [
            [('H', 0), ('T', 0), ('H', 0)],
            [('T', 0), ('H', 0), ('S', 0)],
            [('H', 0), ('TDG', 0), ('H', 0), ('T', 0)],
        ]
        def target_word(depth, prefix):
            angle = 0.19 + 0.13 * depth + 0.07 * prefix
            return np.array([[np.cos(angle), -np.sin(angle)],
                             [np.sin(angle), np.cos(angle)]])

        for height in (1, 2, 3):
            size = 1 << height
            target = _small_frame(height, target_word)
            coarse = _small_frame(
                height, lambda depth, prefix: _word_matrix(
                    1, native_words[(depth + prefix) % len(native_words)]))
            support = [_column_support(marker, height) for marker in range(size)]
            permitted = np.array([[bool(left & right) for right in support]
                                  for left in support])
            for frame in (target, coarse):
                np.testing.assert_allclose(frame.conj().T @ frame, np.eye(size),
                                           atol=ATOL, rtol=0)
                for marker in range(size):
                    outside = sorted(set(range(size)) - support[marker])
                    np.testing.assert_allclose(frame[outside, marker], 0,
                                               atol=ATOL, rtol=0)
            residual = coarse.conj().T @ target - np.eye(size)
            np.testing.assert_allclose(residual[~permitted], 0, atol=ATOL, rtol=0)
            self.assertEqual(int(permitted.sum()), (2 * height - 1) * size + 2)
            self.assertGreater(np.linalg.norm(coarse.imag), 0.1)

    def test_geometric_group_integer_dirty_and_error_ledgers(self):
        # The proof chooses C sufficiently large for the exact one-qubit
        # synthesis primitive. This finite test isolates the integer dirty,
        # precision, and lookup ledger from that unspecified fixed constant.
        ratio, base = 64, 256
        for n in (1, 2, 64, 255, 256, 257, 511, 1024, 10_000, 1_000_000):
            groups = _group_ledger(n, ratio, base)
            self.assertEqual(sum(s for _, s, _ in groups) + min(n, base), n)
            self.assertLessEqual(len(groups), 2 * ratio * n.bit_length())
            self.assertLessEqual(sum(r for r, _, _ in groups), 2 * ratio * n)
            self.assertLessEqual(sum(r // 4 for r, _, _ in groups), ratio * n)
            normalized_lookup = 0.0
            normalized_error = 0.0
            previous_end = min(n, base)
            for r, s, label_bits in groups:
                self.assertEqual(r, previous_end)
                self.assertLessEqual(s, r // ratio)
                previous_end = r + s
                if s <= 100:
                    padded_entries = 4 * ((2 * s - 1) * (1 << s) + 2)
                    self.assertEqual(label_bits, (padded_entries - 1).bit_length())
                address_bits = n - r - s + label_bits + 2
                for precision in (6, 31, 1000):
                    source_width = precision + r // 4 + 8
                    self.assertLessEqual(source_width + address_bits + 8,
                                         precision + n + 7)
                # Q/N = 2**(label_bits-s-r), with the unchanged prefix used
                # directly as an address. There is no extra copied prefix.
                normalized_lookup += math.ldexp(1.0, label_bits - s - r)
                normalized_error += 10 * math.ldexp(1.0, -(r // 4) - 8)
            self.assertLessEqual(normalized_lookup, 32 * (base + 1) * 2.0 ** -base)
            self.assertLessEqual(normalized_error, 80 * 2.0 ** (-8 - base // 4))
            baseline_error = 5 * np.sqrt(2) / 8 * (1 - 2.0 ** -min(n, base))
            self.assertLess(baseline_error + normalized_error, 1)

    def test_star_forests_partition_support_and_reconstruct_complex_entries(self):
        rng = np.random.default_rng(20261001)
        for height in (1, 2, 3, 4):
            size = 1 << height
            support = [_column_support(marker, height) for marker in range(size)]
            permitted = np.array([[bool(left & right) for right in support]
                                  for left in support])
            multiplicity = np.eye(size, dtype=int)
            forests = {}
            for kind, depth in [('zero', 0)] + [('node', j) for j in range(height)]:
                dilation, valid = _star_column_dilation(height, kind, depth)
                accepted = dilation[:size, :size]
                nonzeros = abs(np.diag(valid) @ accepted) > ATOL
                multiplicity += nonzeros.astype(int) + nonzeros.T.astype(int)
                forests[kind, depth] = accepted, valid
                np.testing.assert_allclose(dilation.conj().T @ dilation,
                                           np.eye(2 * size), atol=ATOL, rtol=0)
                self.assertLessEqual(np.linalg.norm(accepted, ord=2), 1 + ATOL)
            np.testing.assert_array_equal(multiplicity, permitted.astype(int))
            np.testing.assert_array_equal(np.diag(multiplicity), 1)
            np.testing.assert_array_equal(multiplicity[0], 1)
            # Missing the zero-at-block-start rule adds ancestor pairs a
            # second time; for example height=2, subtree [2,4), marker=3.
            if height >= 2:
                _, valid = forests['node', height - 1]
                self.assertFalse(valid[2])

            entries = (rng.uniform(-1, 1, (size, size))
                       + 1j * rng.uniform(-1, 1, (size, size))) * permitted
            entries /= 64 * (height + 1) * size
            labels = _star_labels(height)
            count = len(labels)
            reconstructed = np.zeros_like(entries)
            for kind, depth, orientation, phase in labels:
                if kind == 'padding':
                    continue
                if kind == 'diagonal':
                    values = np.diag(entries)
                    accepted = np.eye(size)
                    scale = count
                else:
                    accepted, valid = forests[kind, depth]
                    ancestors = np.argmax(abs(accepted), axis=1)
                    rows = np.arange(size)
                    values = (entries[ancestors, rows] if orientation
                              else entries[rows, ancestors]) * valid
                    subtree_size = size if kind == 'zero' else 1 << (height - depth)
                    scale = count * np.sqrt(subtree_size)
                split = (values.real, values.imag, -values.real, -values.imag)[phase]
                coefficients = scale * np.maximum(split, 0)
                self.assertLessEqual(float(coefficients.max()), 0.25)
                term = np.diag(coefficients) @ accepted
                if orientation:
                    term = term.conj().T
                reconstructed += 1j ** phase * term / count
            np.testing.assert_allclose(reconstructed, entries, atol=ATOL, rtol=0)

    def test_native_logical_address_star_terms_all_dirty_columns(self):
        block = _SmallStarBlock(2)
        data, size = block.data_size, block.logical_size
        embedding = np.zeros((4 * data, data), dtype=complex)
        embedding[:data] = np.eye(data)
        state = embedding.reshape(2, 2, size, block.core_size, data)
        for index, (kind, depth, orientation, phase) in enumerate(block.labels):
            actual = block.term(state.copy(), index)
            if kind == 'padding':
                expected = np.zeros((size, size), dtype=complex)
            else:
                accepted = block.dilations[kind, depth][:size, :size]
                expected = np.diag(block.coefficients[index]) @ accepted
                if orientation:
                    expected = expected.conj().T
                expected = expected * 1j ** phase
            np.testing.assert_allclose(actual[0, 0].reshape(data, data),
                                       np.kron(expected, np.eye(block.core_size)),
                                       atol=ATOL, rtol=0)
            np.testing.assert_allclose(block.term(actual, index, inverse=True), state,
                                       atol=ATOL, rtol=0)

        # A shared failure flag does not implement the column-star product.
        index = next(i for i, label in enumerate(block.labels)
                     if label == ('zero', 0, 0, 0))
        shared = np.zeros((2 * data, data), dtype=complex)
        shared[:data] = np.eye(data)
        shared = shared.reshape(2, size, block.core_size, data)
        shared = _apply_axes(shared, block.dilations['zero', 0], (0, 1))
        for logical, coefficient in enumerate(block.coefficients[index]):
            shared[:, logical] = _apply_axes(shared[:, logical],
                                             block.scalars[coefficient], (0, 1))
        expected = (np.diag(block.coefficients[index])
                    @ block.dilations['zero', 0][:size, :size])
        discrepancy = shared[0].reshape(data, data) - np.kron(
            expected, np.eye(block.core_size))
        self.assertGreater(np.linalg.norm(discrepancy, ord=2), 0.5)

    def test_star_lcu_rounding_has_no_subtree_dimension_amplification(self):
        block = _SmallStarBlock(2)
        embedding = block.embedding()
        qj = block.apply(embedding)
        b = qj[:block.data_size]
        expected = (np.eye(block.data_size) + block.residual()) / 2
        np.testing.assert_allclose(b, expected, atol=ATOL, rtol=0)
        np.testing.assert_allclose(block.apply(qj, inverse=True), embedding,
                                   atol=ATOL, rtol=0)

        rounding_error = 0.031
        perturbation = np.array([
            [rounding_error * (-1) ** (label + 2 * logical + logical // 2)
             for logical in range(block.logical_size)]
            for label in range(block.count)
        ])
        desired = np.eye(block.data_size) + block.residual(perturbation)
        self.assertLessEqual(np.linalg.norm(2 * b - desired, ord=2),
                             rounding_error + ATOL)
        # The same bound is per forest, including large uniform subtrees.
        for height in (1, 2, 3, 4):
            dimension = 1 << height
            errors = rounding_error * np.array([(-1) ** z for z in range(dimension)])
            for kind, depth in [('zero', 0)] + [('node', j) for j in range(height)]:
                dilation, _ = _star_column_dilation(height, kind, depth)
                accepted = dilation[:dimension, :dimension]
                discrepancy = np.diag(errors) @ accepted
                self.assertLessEqual(np.linalg.norm(discrepancy, ord=2),
                                     rounding_error + ATOL)

        # Audit full rejected-space error for the new SELECT, not just its
        # projected block. The polar factor supplies an independent nearby
        # unitary target for this finite interface fixture.
        approximate = np.eye(block.data_size) + block.residual()
        left, _, right = np.linalg.svd(approximate)
        target = left @ right
        amplified = block.amplify(embedding)
        zeta = np.linalg.norm(2 * b - target, ord=2)
        error = np.linalg.norm(amplified - embedding @ target, ord=2)
        self.assertLessEqual(error, 4 * zeta + ATOL)
        np.testing.assert_allclose(amplified[:block.data_size],
                                   3 * b - 4 * b @ b.conj().T @ b,
                                   atol=ATOL, rtol=0)

    def test_iterated_group_integer_ledger_without_exponential_allocations(self):
        # These illustrative constants test rounding and bookkeeping only;
        # they do not certify the proof's unknown fixed circuit constants.
        ratio, base = 8, 64
        for n in (1, 63, 64, 65, 100, 1000, 10**6, 10**30, 1 << 4096):
            groups = _iterated_group_ledger(n, ratio, base)
            self.assertEqual(sum(s for _, s, _ in groups) + min(n, base), n)
            lookup_ratio = 0.0
            normalized_error = 0.0
            previous_end = min(n, base)
            for r, s, label_bits in groups:
                self.assertEqual(r, previous_end)
                previous_end = r + s
                self.assertLessEqual(s.bit_length() - 1, r // ratio)
                self.assertEqual(label_bits, (8 * s + 11).bit_length())
                # Prefix plus current logical address, label, h, and mode.
                address_bits = n - r + label_bits + 2
                source_extra = r // 4 + 8
                self.assertLessEqual(source_extra + address_bits + 8, n + 7)
                lookup_exponent = label_bits - r
                if lookup_exponent > -1075:
                    lookup_ratio += math.ldexp(1.0, lookup_exponent)
                error_exponent = -(r // 4) - 8
                if error_exponent > -1075:
                    normalized_error += 10 * math.ldexp(1.0, error_exponent)
            self.assertLess(lookup_ratio, 1e-12)
            self.assertLessEqual(normalized_error, 80 * 2.0 ** (-8 - base // 4))
            baseline_error = 5 * np.sqrt(2) / 8 * (1 - 2.0 ** -min(n, base))
            self.assertLess(baseline_error + normalized_error, 1)
            self.assertLessEqual(sum(r // 4 for r, _, _ in groups), n)
            self.assertLessEqual(len(groups), 2 + 2 * n.bit_length().bit_length())

    def test_streamed_coarse_symbols_match_stored_literal_words(self):
        programs = [[1, 2, 1], [2, 1, 3]]
        expected = np.zeros((4, 4), dtype=complex)
        expected[:2, :2] = _word_matrix(1, [('H', 0), ('T', 0), ('H', 0)])
        expected[2:, 2:] = _word_matrix(1, [('T', 0), ('H', 0), ('TDG', 0)])
        accepted = []
        for streamed in (False, True):
            unitary = _coarse_interpreter_matrix(programs, streamed)
            active = np.arange(4, 8)  # h=1, every address/target, zero program.
            columns = unitary[:, active]
            target = np.zeros_like(columns)
            target[active] = expected
            np.testing.assert_allclose(columns, target, atol=ATOL, rtol=0)
            accepted.append(unitary[np.ix_(active, active)])
            # Inactive program bits may be arbitrary. A zero inactive query
            # alone is insufficient: the interpreter itself is h-controlled.
            inactive = [basis for basis in range(unitary.shape[0])
                        if not ((basis >> 2) & 1)]
            np.testing.assert_allclose(unitary[:, inactive],
                                       np.eye(unitary.shape[0])[:, inactive],
                                       atol=0, rtol=0)
        np.testing.assert_allclose(accepted[0], accepted[1], atol=ATOL, rtol=0)


if __name__ == '__main__':
    unittest.main()
