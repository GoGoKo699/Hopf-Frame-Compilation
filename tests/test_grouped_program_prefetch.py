"""Finite audits of conditional program prefetch and grouped source use.

The mask emitter uses literal X/CX/CCX/CZ words and exact basis propagation;
the five-leaf AND is separately expanded to native Clifford+T. The small
group fixtures retain native PREP auxiliaries and all source/branch leakage,
but represent the already-audited mask and reflection as diagonal operators.
Program loading is an explicit one-address-bit XOR word. Sparse program
labels retain coherent external addresses without allocating 2**w columns.
These checks are not a general lookup emitter or an asymptotic depth proof.
"""
from __future__ import annotations

import unittest

import numpy as np

try:
    from .test_conditional_geometric_source import _ConditionalGeometricSource, K
    from .test_operator_source_compiler import (
        _adjoint, _apply_native_word, _expand_toffolis,
    )
except ImportError:
    from test_conditional_geometric_source import _ConditionalGeometricSource, K
    from test_operator_source_compiler import (
        _adjoint, _apply_native_word, _expand_toffolis,
    )


ATOL = 8e-11


def _basis_word(word, inputs):
    """Exact phase/permutation propagation, including words wider than 64 bits."""
    outputs = np.array(list(inputs), dtype=object)
    phases = np.ones(len(outputs), dtype=np.int8)
    for name, *qubits in word:
        enabled = np.ones(len(outputs), dtype=bool)
        controls = qubits if name == 'CZ' else qubits[:-1]
        for q in controls:
            enabled &= np.asarray((outputs >> q) & 1, dtype=bool)
        if name == 'CZ':
            phases[enabled] *= -1
        else:
            assert name in ('X', 'CX', 'CCX')
            outputs[enabled] ^= 1 << qubits[-1]
    return outputs, phases


def _five_and(leaves, nodes):
    """Four Toffolis in three disjoint-gate rounds."""
    return [[('CCX', leaves[0], leaves[1], nodes[0]),
             ('CCX', leaves[2], leaves[3], nodes[1])],
            [('CCX', nodes[0], nodes[1], nodes[2])],
            [('CCX', nodes[2], leaves[4], nodes[3])]]


class _MaskWord:
    def __init__(self, m, ell):
        self.m, self.ell = m, ell
        self.h, self.u, self.b = 0, 1, 2
        self.core = list(range(3, 3 + m))
        self.selectors = list(range(3 + m, 3 + m + (1 << ell)))
        start = self.selectors[-1] + 1
        self.terms = [(p, beta, j) for p in range(1 << ell)
                      for beta in range(2) for j in range(m)]
        self.program = list(range(start, start + len(self.terms)))
        self.work_start = self.program[-1] + 1
        copies, rounds, roots = [], [[], [], []], []
        for term, (p, beta, j) in enumerate(self.terms):
            work = list(range(self.work_start + 9 * term,
                              self.work_start + 9 * (term + 1)))
            leaves, nodes = work[:5], work[5:]
            originals = [self.selectors[p], self.u, self.b,
                         self.program[term], self.core[j]]
            copies += [('CX', original, copy) for original, copy in zip(originals, leaves)]
            if beta == 0:
                copies.append(('X', leaves[2]))
            for level, gates in enumerate(_five_and(leaves, nodes)):
                rounds[level] += gates
            roots.append(nodes[-1])
        self.rounds = rounds
        self.compute = copies + [gate for level in rounds for gate in level]
        self.word = (self.compute + [('CZ', self.h, root) for root in roots]
                     + _adjoint(self.compute))
        self.width = self.work_start + 9 * len(self.terms)

    def sign(self, basis):
        exponent = 0
        for term, (p, beta, j) in enumerate(self.terms):
            exponent ^= (((basis >> self.selectors[p]) & 1)
                         & ((basis >> self.u) & 1)
                         & (((basis >> self.b) & 1) == beta)
                         & ((basis >> self.program[term]) & 1)
                         & ((basis >> self.core[j]) & 1))
        return (-1) ** (((basis >> self.h) & 1) & exponent)


def _program_offset(m, ell, p, beta):
    return 2 * m * ((1 << ell) - 1) + m * (2 * p + beta)


def _pack_program(m, rows):
    return sum(mask << _program_offset(m, ell, p, beta)
               for ell, layer in enumerate(rows)
               for p, masks in enumerate(layer) for beta, mask in enumerate(masks))


def _loader_word(program_bits, tables):
    # Label qubits: x=0,h=1, followed by the program. No hidden helper.
    word = []
    for x, table in enumerate(tables):
        negate = [('X', 0)] if x == 0 else []
        word += negate + [('CCX', 1, 0, 2 + j) for j in range(program_bits)
                          if (table >> j) & 1] + negate
    return word


def _apply_label_word(word, state):
    """Sparse coherent labels; each value is a data-column batch."""
    labels = list(state)
    images, phases = _basis_word(word, labels)
    return {int(image): phase * state[label]
            for label, image, phase in zip(labels, images, phases)}


class _ReducedGroup:
    """Actual native PREP, with diagonal program mask/reflection represented directly.

    The private mask/selector work is omitted on its proved invariant zero
    subspace. Native preparation auxiliaries remain explicit. Internal u
    and prefix selectors are evaluated from the unchanged local control bits.
    """

    def __init__(self, m, g):
        assert 1 <= g <= 2 and 2 <= m <= 3
        self.m, self.g = m, g
        self.source = _ConditionalGeometricSource(m)
        self.local = list(range(self.source.prep_width, self.source.prep_width + g))
        self.b = self.local[-1] + 1
        self.width = self.b + 1
        self.indices = np.arange(1 << self.width)
        self.core = self.indices & ((1 << m) - 1)
        self.locals = (self.indices >> self.local[0]) & ((1 << g) - 1)
        self.branch = (self.indices >> self.b) & 1
        self.weights = self.source.amplitudes ** 2

    def masks(self, program, ell, p):
        return tuple((program >> _program_offset(self.m, ell, p, beta))
                     & ((1 << self.m) - 1) for beta in range(2))

    def coefficients(self, masks):
        return tuple(sum(weight * (-1) ** ((mask >> j) & 1)
                         for j, weight in enumerate(self.weights)) for mask in masks)

    def apply_q(self, columns, program, ell, h, inverse=False):
        enable = bool(h) & ((self.locals >> (ell + 1)) == 0)
        prefix = self.locals & ((1 << ell) - 1)
        masks = np.array([self.masks(program, ell, int(p))[int(beta)]
                          for p, beta in zip(prefix, self.branch)])
        parity = np.array([(int(core) & int(mask)).bit_count() & 1
                           for core, mask in zip(self.core, masks)])
        phase = (-1.0) ** (enable & parity)
        select = enable & self.branch
        target = (self.indices >> self.local[ell]) & 1
        kphase = (-1.0) ** (select & target)
        permutation = self.indices ^ (select.astype(int) << self.local[ell])

        def apply_select(state, adjoint=False):
            if adjoint:
                return state[permutation] * kphase[:, None]
            return (state * kphase[:, None])[permutation]

        result = _apply_native_word(self.width, [('H', self.b)], columns)
        if inverse:
            result = apply_select(result, adjoint=True)
        result = _apply_native_word(self.width, self.source.prep, result)
        result *= phase[:, None]
        result = _apply_native_word(self.width, _adjoint(self.source.prep), result)
        if not inverse:
            result = apply_select(result)
        return _apply_native_word(self.width, [('H', self.b)], result)

    def layer(self, columns, program, ell, h=1, wrong_minus=False):
        enable = bool(h) & ((self.locals >> (ell + 1)) == 0)
        reflection = np.ones(len(self.indices))
        reflection[enable & (self.branch == 0) & (self.core == 0)] = -1
        result = self.apply_q(columns, program, ell, h)
        result *= reflection[:, None]
        result = self.apply_q(result, program, ell, h, inverse=True)
        result *= reflection[:, None]
        result = self.apply_q(result, program, ell, h)
        minus = np.full(len(self.indices), -1.0 if h else 1.0) if wrong_minus else (-1.0) ** enable
        return result * minus[:, None]

    def apply(self, columns, program, h=1):
        result = columns.copy()
        for ell in range(self.g):
            result = self.layer(result, program, ell, h)
        return result

    def embedding(self):
        result = np.zeros((1 << self.width, 1 << self.g), dtype=complex)
        result[np.arange(1 << self.g) << self.local[0], np.arange(1 << self.g)] = 1
        return result

    def ideal(self, program):
        result = np.eye(1 << self.g, dtype=complex)
        error_bound = 0.0
        for ell in range(self.g):
            layer = np.eye(1 << self.g, dtype=complex)
            worst = 0.0
            for p in range(1 << ell):
                c, s = self.coefficients(self.masks(program, ell, p))
                radius = np.hypot(c, s)
                assert radius > 0
                rows = [p, p | (1 << ell)]
                layer[np.ix_(rows, rows)] = (c * np.eye(2) + s * K) / radius
                worst = max(worst, abs(radius - 1))
            result = layer @ result
            error_bound += 4 * worst
        return result, error_bound


class GroupedProgramPrefetchTests(unittest.TestCase):
    def test_five_input_tree_native_phase_count_and_rounds(self):
        rounds = _five_and(list(range(5)), list(range(5, 9)))
        word = [gate for level in rounds for gate in level]
        for level in rounds:
            support = [q for gate in level for q in gate[1:]]
            self.assertEqual(len(support), len(set(support)))
        native = _expand_toffolis(word)
        self.assertEqual(sum(gate[0] in ('T', 'TDG') for gate in native), 28)
        inputs = np.zeros((512, 32), dtype=complex)
        inputs[np.arange(32), np.arange(32)] = 1
        images, phases = _basis_word(word, range(32))
        expected = np.zeros_like(inputs)
        expected[np.asarray(images, dtype=int), np.arange(32)] = phases
        actual = _apply_native_word(9, native, inputs)
        np.testing.assert_allclose(actual, expected, atol=ATOL, rtol=0)
        np.testing.assert_allclose(_apply_native_word(9, _adjoint(native), actual),
                                   inputs, atol=ATOL, rtol=0)
        self.assertEqual(len(rounds) * 4 * 2, 24)

    def test_one_term_inactive_identity_on_every_dirty_work_input(self):
        # h=0; five originals and all nine work bits are arbitrary (2**14 columns).
        leaves, nodes = list(range(6, 11)), list(range(11, 15))
        compute = [('CX', j + 1, leaves[j]) for j in range(5)]
        compute += [gate for level in _five_and(leaves, nodes) for gate in level]
        word = compute + [('CZ', 0, nodes[-1])] + _adjoint(compute)
        inputs = list(range(0, 1 << 15, 2))
        images, phases = _basis_word(word, inputs)
        np.testing.assert_array_equal(images, inputs)
        np.testing.assert_array_equal(phases, 1)
        # Moving the central guard onto a dirty copy is not an inactive identity.
        wrong = compute + [('CZ', leaves[0], nodes[-1])] + _adjoint(compute)
        _, wrong_phases = _basis_word(wrong, inputs)
        self.assertTrue(np.any(wrong_phases == -1))

    def test_complete_mask_zero_work_all_data_and_readonly_program(self):
        for m, ell in ((2, 0), (2, 1)):
            with self.subTest(m=m, ell=ell):
                mask = _MaskWord(m, ell)
                # Enumerate all original selector/program/core/u/b values, h=1.
                inputs = list(range(1, 1 << mask.work_start, 2))
                images, phases = _basis_word(mask.word, inputs)
                np.testing.assert_array_equal(images, inputs)
                np.testing.assert_array_equal(phases, [mask.sign(basis) for basis in inputs])
                for gate in mask.word:
                    if gate[0] in ('X', 'CX', 'CCX'):
                        self.assertNotIn(gate[-1], mask.program)
                for level in mask.rounds:
                    support = [q for gate in level for q in gate[1:]]
                    self.assertEqual(len(support), len(set(support)))
                # Omitting the inverse leaves a real work-return failure.
                computed, _ = _basis_word(mask.compute, inputs)
                self.assertTrue(any(int(basis) >> mask.work_start for basis in computed))

    def test_xor_preload_actual_inverse_on_every_small_input(self):
        word = _loader_word(4, (0b0110, 0b1001))
        inputs = list(range(64))
        loaded, phases = _basis_word(word, inputs)
        expected = [basis ^ (((0b0110, 0b1001)[basis & 1] << 2)
                              if (basis >> 1) & 1 else 0) for basis in inputs]
        np.testing.assert_array_equal(loaded, expected)
        np.testing.assert_array_equal(phases, 1)
        restored, phases = _basis_word(_adjoint(word), loaded)
        np.testing.assert_array_equal(restored, inputs)
        np.testing.assert_array_equal(phases, 1)

    def test_group_native_prep_composition_and_complete_error(self):
        fixtures = ((2, 1, [[(0, 1)]]),
                    (2, 2, [[(1, 0)], [(0, 1), (3, 1)]]),
                    (3, 2, [[(0, 2)], [(2, 0), (7, 2)]]))
        for m, g, rows in fixtures:
            with self.subTest(m=m, g=g):
                group = _ReducedGroup(m, g)
                program = _pack_program(m, rows)
                embedding = group.embedding()
                actual = group.apply(embedding, program)
                ideal, bound = group.ideal(program)
                error = np.linalg.norm(actual - embedding @ ideal, ord=2)
                self.assertLessEqual(error, bound + ATOL)
                np.testing.assert_allclose(actual.conj().T @ actual,
                                           np.eye(1 << g), atol=ATOL, rtol=0)
                if m == 2:
                    np.testing.assert_allclose(actual, embedding @ ideal, atol=ATOL, rtol=0)
                else:
                    self.assertGreater(error, 0.01)
                    # Native PREP auxiliaries still return exactly despite source leakage.
                    helper_mask = sum(1 << q for q in group.source.helpers)
                    np.testing.assert_allclose(actual[(group.indices & helper_mask) != 0],
                                               0, atol=ATOL, rtol=0)

    def test_inner_enable_and_literal_minus_on_leaked_inputs(self):
        group = _ReducedGroup(3, 2)
        program = _pack_program(3, [[(0, 2)], [(2, 0), (7, 2)]])
        helper_mask = sum(1 << q for q in group.source.helpers)
        # Every core, b, and first target input; t1=1 makes the first layer inactive.
        indices = group.indices[((group.locals >> 1) == 1)
                                & ((group.indices & helper_mask) == 0)]
        columns = np.eye(1 << group.width, dtype=complex)[:, indices]
        np.testing.assert_allclose(group.layer(columns, program, 0), columns,
                                   atol=ATOL, rtol=0)
        np.testing.assert_allclose(group.layer(columns, program, 0, wrong_minus=True),
                                   -columns, atol=ATOL, rtol=0)
        # Outer inactivity includes arbitrary native PREP auxiliaries as well.
        all_columns = np.eye(1 << group.width, dtype=complex)
        np.testing.assert_allclose(group.apply(all_columns, program, h=0), all_columns,
                                   atol=ATOL, rtol=0)

    def test_coherent_prefetch_unload_with_arbitrary_core_and_branch(self):
        group = _ReducedGroup(3, 2)
        tables = tuple(_pack_program(3, rows) for rows in
                       ([[ (0, 2)], [(2, 0), (7, 2)]],
                        [[ (2, 0)], [(2, 7), (0, 2)]]))
        program_bits = 2 * group.m * ((1 << group.g) - 1)
        loader = _loader_word(program_bits, tables)
        helper_mask = sum(1 << q for q in group.source.helpers)
        indices = group.indices[(group.indices & helper_mask) == 0]
        columns = np.eye(1 << group.width, dtype=complex)[:, indices]
        state = {}
        for x in range(2):
            # Every core/b/local basis input, with a coherent, column-dependent
            # external-prefix phase. All native PREP auxiliaries remain zero.
            data = columns * (1j ** (x * np.arange(len(indices)))) / np.sqrt(2)
            state[x | 2] = data  # Coherent x, h=1, W=0; core and b are arbitrary.
        loaded = _apply_label_word(loader, state)
        for ell in range(group.g):
            loaded = {label: group.layer(data, label >> 2, ell, (label >> 1) & 1)
                      for label, data in loaded.items()}
            self.assertEqual(set(loaded), {x | 2 | (tables[x] << 2) for x in range(2)})
        restored = _apply_label_word(_adjoint(loader), loaded)
        self.assertEqual(set(restored), set(state))
        for x in range(2):
            np.testing.assert_allclose(restored[x | 2], group.apply(state[x | 2], tables[x]),
                                       atol=ATOL, rtol=0)
        self.assertTrue(all(label >> 2 for label in loaded))  # Missing unload leaves W nonzero.
        np.testing.assert_allclose(sum(data.conj().T @ data for data in restored.values()),
                                   np.eye(len(indices)), atol=ATOL, rtol=0)

    def test_simultaneous_reservation_and_mask_tree_counts(self):
        for m in range(2, 33):
            for g in range(1, m + 1):
                base = 1 << g
                program = 2 * m * (base - 1)
                total = program + 7 * m + g * base + (2 * g + 1) + 9 * m * base
                self.assertEqual(total, 11 * m * base + g * base + 5 * m + 2 * g + 1)
                self.assertLessEqual(total, 16 * m * base)
        for m, ell in ((2, 0), (2, 1), (3, 1)):
            mask = _MaskWord(m, ell)
            terms = 2 * m * (1 << ell)
            self.assertEqual(mask.width - mask.work_start, 9 * terms)
            self.assertEqual(sum(gate[0] == 'CCX' for gate in mask.word), 8 * terms)

    def test_global_group_schedule_integer_interfaces(self):
        # Finite interface checks only; asymptotic count sums remain analytic.
        starts_checked = 0
        for n in (1, 32, 2048, 8192, 16384):
            for precision in (6, 8, 16):
                with self.subTest(n=n, precision=precision):
                    cap = (8 * n - 1).bit_length()  # Exact ceil(log2(8n)).
                    m = precision + 4 + cap
                    cutoff = min(n, 256 * m)
                    remaining = n
                    while remaining > cutoff:
                        # floor(log2(k/(64m))) without floating logarithms.
                        g = (remaining // (64 * m)).bit_length() - 1
                        self.assertGreaterEqual(g, 2)
                        self.assertLessEqual(g, m)
                        self.assertLessEqual(16 * m * (1 << g), remaining - g)
                        program = 2 * m * ((1 << g) - 1)
                        self.assertLessEqual(32 * program, remaining)
                        self.assertGreater(remaining - g + 1, cap)
                        # A conservative integer certificate for
                        # w(n+2)^3 <= 2**(k/4), without constructing N=2**n.
                        polynomial = program * (n + 2) ** 3
                        self.assertLessEqual(4 * polynomial.bit_length(), remaining)
                        remaining -= g
                        starts_checked += 1
                    self.assertLessEqual(remaining, cutoff)
        self.assertGreater(starts_checked, 0)


if __name__ == '__main__':
    unittest.main()
