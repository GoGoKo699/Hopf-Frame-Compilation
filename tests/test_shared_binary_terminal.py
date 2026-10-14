"""Bounded boundary audits for a shared binary terminal source.

The pinned shift-kernel receipts already check generic capped-source,
history, overflow, and scalar-budget identities. This supplement adds a
two-stage addressed 16-type kernel, native guarded shift/history words,
and an occupied-inactive source wrapper. Exact sparse arithmetic retains
every kernel failure amplitude; only tiny native primitives use matrices.
M=4, J=2 and p=1/2 test local algebra, not production precision, a complete
frame theorem, or a scalable table/arithmetic emitter. The algebra fixture
supplies coefficient-table actions and folds the local suffix predicate
into atom rejection; native queries and their unload are covered elsewhere.
No receipt is edited.
"""
from __future__ import annotations

from fractions import Fraction as F
import unittest

import numpy as np

try:
    from . import exact_matrices as exact
    from .test_conditional_geometric_source import _native
    from .test_grouped_terminal_depth import _mcx, _phase_s
    from .test_operator_source_compiler import _adjoint, _apply_native_word, _word_matrix
except ImportError:
    import exact_matrices as exact
    from test_conditional_geometric_source import _native
    from test_grouped_terminal_depth import _mcx, _phase_s
    from test_operator_source_compiler import _adjoint, _apply_native_word, _word_matrix

# exact_matrices imports these pinned helpers with a temporary import path.
from gray_geometric import emit_chain, emit_decoder


ATOL = 2e-10
ZERO, ONE, IMAG = exact.CZ, exact.CO, exact.CI
HALF_ROOT = exact.CQ2(exact.Q2(0, F(1, 2)))
PHASES = (ONE, -ONE, IMAG, -IMAG)
SOURCE, DATA, MODE, TYPE, DIGIT, UNDER, ATOM, ZERO_DIGIT, COUNTER = range(9)
PRIVATE = (MODE, TYPE, DIGIT, UNDER, ATOM, ZERO_DIGIT)
LABEL_BITS = ((MODE, 0), *((TYPE, bit) for bit in range(4)), (DIGIT, 0))

# Keys are preserved prefix, then constant (u,v,phase). The rotating local
# target is never a coefficient-table address. Missing entries are zero.
TABLES = (
    {0: {(1, 0, 2): (0, 1), (0, 1, 1): (1, 0), (0, 0, 0): (0, 1)}},
    {0: {(0, 1, 3): (0, 1), (1, 0, 0): (1, 1)},
     1: {(1, 0, 1): (0, 1), (0, 1, 2): (1, 0), (1, 1, 3): (1, 1)}},
)
COARSE_WORDS = ({0: [('H', 0), ('T', 0)]},
                {0: [('H', 0)], 1: [('TDG', 0), ('H', 0)]})
COARSE = tuple({prefix: exact.native_columns(1, word) for prefix, word in rows.items()}
               for rows in COARSE_WORDS)


def _accumulate(state, key, value):
    value = state.get(key, ZERO) + value
    if value == ZERO:
        state.pop(key, None)
    else:
        state[key] = value


def _replace(basis, field, value):
    result = list(basis)
    result[field] = value
    return tuple(result)


def _hadamard(state, field, bit):
    result = {}
    for basis, amplitude in state.items():
        old = (basis[field] >> bit) & 1
        low = basis[field] & ~(1 << bit)
        _accumulate(result, _replace(basis, field, low), HALF_ROOT * amplitude)
        _accumulate(result, _replace(basis, field, low | (1 << bit)),
                    (-1 if old else 1) * HALF_ROOT * amplitude)
    return result


def _select(state, stage, inverse=False):
    result = {}
    for basis, amplitude in state.items():
        if basis[MODE]:
            _accumulate(result, basis, amplitude)
            continue
        kind = basis[TYPE]
        u, v, phase = kind & 1, (kind >> 1) & 1, kind >> 2
        prefix = basis[DATA] & ((1 << stage) - 1)
        bit = TABLES[stage][prefix].get((u, v, phase), (0, 0))[basis[DIGIT]]
        if not bit:
            _accumulate(result, _replace(basis, ZERO_DIGIT, basis[ZERO_DIGIT] ^ 1), amplitude)
            continue
        output = list(basis)
        displacement = 2 * basis[DIGIT]
        if inverse:
            # Actual inverse: restore the source index BEFORE testing underflow.
            output[SOURCE] = (output[SOURCE] + displacement) % 4
            output[UNDER] ^= int(output[SOURCE] < displacement)
            output[DATA] ^= (u ^ v) << stage
            mismatch = ((output[DATA] >> stage) & 1) != v or (output[DATA] >> (stage + 1)) != 0
            output[ATOM] ^= int(mismatch)
            amplitude *= PHASES[phase].conj()
        else:
            mismatch = ((output[DATA] >> stage) & 1) != v or (output[DATA] >> (stage + 1)) != 0
            output[ATOM] ^= int(mismatch)
            output[DATA] ^= (u ^ v) << stage
            output[UNDER] ^= int(output[SOURCE] < displacement)
            output[SOURCE] = (output[SOURCE] - displacement) % 4
            amplitude *= PHASES[phase]
        _accumulate(result, tuple(output), amplitude)
    return result


def _kernel(state, stage, inverse=False):
    for field, bit in LABEL_BITS:
        state = _hadamard(state, field, bit)
    state = _select(state, stage, inverse)
    for field, bit in reversed(LABEL_BITS):
        state = _hadamard(state, field, bit)
    return state


def _coarse(state, stage, inverse=False):
    result = {}
    for basis, amplitude in state.items():
        if basis[DATA] >> (stage + 1):
            _accumulate(result, basis, amplitude)
            continue
        prefix = basis[DATA] & ((1 << stage) - 1)
        matrix = COARSE[stage][prefix]
        if inverse:
            matrix = exact.adjoint(matrix)
        old = (basis[DATA] >> stage) & 1
        for new in (0, 1):
            _accumulate(result, _replace(basis, DATA, basis[DATA] ^ ((old ^ new) << stage)),
                        matrix[new][old] * amplitude)
    return result


def _history(state, inverse=False):
    # One bit suffices for exactly one tick between two stages; no wrap occurs
    # from counter zero. The modular action is defined on every other input.
    return {_replace(basis, COUNTER,
                     (basis[COUNTER] + (-1 if inverse else 1)
                      * int(any(basis[field] for field in PRIVATE))) % 2): amplitude
            for basis, amplitude in state.items()}


def _stream(state, *, history=True, inverse=False):
    if inverse:
        state = _kernel(_coarse(state, 1, True), 1, True)
        if history:
            state = _history(state, True)
        return _kernel(_coarse(state, 0, True), 0, True)
    state = _coarse(_kernel(state, 0), 0)
    if history:
        state = _history(state)
    return _coarse(_kernel(state, 1), 1)


def _basis(source, data):
    return (source, data, 0, 0, 0, 0, 0, 0, 0)


def _accepted(state):
    return {basis: value for basis, value in state.items()
            if not basis[COUNTER] and not any(basis[field] for field in PRIVATE)}


def _expected_stage(stage):
    # Independent matrix-unit construction of the F=0 block, with Q=16,J=2.
    contraction = exact.scale(exact.eye(16), F(1, 2))
    for prefix, row in TABLES[stage].items():
        for (u, v, phase), digits in row.items():
            column_data, row_data = prefix | (v << stage), prefix | (u << stage)
            for digit, enabled in enumerate(digits):
                if enabled:
                    for source in range(2 * digit, 4):
                        contraction[4 * (source - 2 * digit) + row_data][4 * source + column_data] += (
                            F(1, 64) * PHASES[phase])
    coarse = exact.zeros(16)
    for column in range(16):
        for basis, value in _coarse({_basis(column // 4, column % 4): ONE}, stage).items():
            coarse[4 * basis[SOURCE] + basis[DATA]][column] = value
    return exact.multiply(coarse, contraction)


def _source_prepare(state, inverse=False):
    # Ordered two-level H word; the native Gray implementation is checked
    # separately below and in the immutable source receipt.
    pairs = ((0, 1), (1, 2), (2, 3))
    for left, right in reversed(pairs) if inverse else pairs:
        result = {}
        for basis, amplitude in state.items():
            if basis[SOURCE] not in (left, right):
                _accumulate(result, basis, amplitude)
            else:
                _accumulate(result, _replace(basis, SOURCE, left), HALF_ROOT * amplitude)
                _accumulate(result, _replace(basis, SOURCE, right),
                            (-1 if basis[SOURCE] == right else 1) * HALF_ROOT * amplitude)
        state = result
    return state


def _guarded_shift(high, digit, underflow, activity, borrow):
    predicate = ([('X', high)] + _mcx([activity, digit, high], underflow, borrow)
                 + [('X', high)])
    subtraction = _mcx([activity, digit], high, borrow)  # Subtract 2k modulo four.
    return predicate + subtraction, predicate, subtraction


def _history_word(activity, private, counter, borrow):
    # INC_H followed by DEC_[H and F=0]. On H=0 both complete arithmetic
    # words are identity even if the private fields/counter/borrow are dirty.
    def increment(controls):
        return [gate for bit in reversed(range(len(counter)))
                for gate in _mcx([*controls, *counter[:bit]], counter[bit], borrow)]
    unconditional = increment([activity])
    negative = [('X', bit) for bit in private]
    zero = negative + increment([activity, *private]) + _adjoint(negative)
    return unconditional + _adjoint(zero)


def _native_source():
    # This is an UNCONDITIONAL full unitary on data 0,1, root 2, scratch 3.
    # Only its active initialized column assumes root/scratch start at zero.
    raw = [('X', 2)] + [(gate.op, *gate.qubits)
                       for gate in emit_chain(2) + emit_decoder(2)] + [('X', 2)]
    return _native(raw)


class SharedBinaryTerminalTests(unittest.TestCase):
    def test_addressed_two_stage_accepted_product_keeps_source_and_literal_phases(self):
        stages = [_expected_stage(stage) for stage in range(2)]
        expected = exact.multiply(stages[1], stages[0])
        actual = exact.zeros(16)
        evolved_columns = []
        for column in range(16):
            initial = {_basis(column // 4, column % 4): ONE}
            evolved = _stream(initial)
            evolved_columns.append(evolved)
            for basis, value in _accepted(evolved).items():
                actual[4 * basis[SOURCE] + basis[DATA]][column] = value
        self.assertEqual(actual, expected)
        # The first exact coarse H/T creates a coherent earlier target; that
        # target selects different coefficient rows at the second stage.
        self.assertNotEqual(TABLES[1][0], TABLES[1][1])
        self.assertTrue(any(value.imag != exact.Q2() for row in actual for value in row))

        g = (HALF_ROOT, ONE * F(1, 2), HALF_ROOT * F(1, 2), HALF_ROOT * F(1, 2))
        embedding = exact.zeros(16, 4)
        for source, value in enumerate(g):
            for data in range(4):
                embedding[4 * source + data][data] = value
        final_expected = exact.multiply(exact.adjoint(embedding), exact.multiply(expected, embedding))
        final = exact.zeros(4)
        for data in range(4):
            evolved = {}
            for source, weight in enumerate(g):
                for basis, amplitude in evolved_columns[4 * source + data].items():
                    _accumulate(evolved, basis, weight * amplitude)
            # Only now is the source inverse applied and its zero output tested.
            for basis, value in _accepted(_source_prepare(evolved, inverse=True)).items():
                if basis[SOURCE] == 0:
                    final[basis[DATA]][data] = value
        self.assertEqual(final, final_expected)
        per_stage = [exact.multiply(exact.adjoint(embedding), exact.multiply(stage, embedding))
                     for stage in stages]
        self.assertNotEqual(final, exact.multiply(per_stage[1], per_stage[0]))
        self.assertEqual(_history({_basis(3, 0): ONE}), {_basis(3, 0): ONE})
        initial = {_basis(0, 0): ONE}
        self.assertNotEqual(_accepted(_stream(initial, history=False)),
                            _accepted(evolved_columns[0]))
        # An occupied source/private input exercises the actual complete inverse.
        occupied = (3, 2, 1, 13, 1, 1, 0, 1, 1)
        self.assertEqual(_stream(_stream({occupied: IMAG}), inverse=True), {occupied: IMAG})

    def test_native_guarded_shift_history_full_inputs_and_actual_inverse(self):
        # source bits 0,1 (little endian), digit 2, underflow 3, H 4,
        # arbitrary phase helper 5 and borrowed arithmetic wire 6.
        shift, predicate, subtraction = _guarded_shift(1, 2, 3, 4, 6)
        phase = _phase_s([4, 2], 5, 6)
        word = shift + phase
        images, factors = [], []
        for value in range(128):
            source, digit, activity = value & 3, (value >> 2) & 1, (value >> 4) & 1
            image = value
            if activity:
                image ^= int(source < 2 * digit) << 3
                image = (image & ~3) | ((source - 2 * digit) % 4)
            images.append(image)
            factors.append(1j ** (activity * digit))
        expected = np.eye(128, dtype=complex)[:, images] * np.asarray(factors)[None, :]
        actual = _word_matrix(7, word)
        np.testing.assert_allclose(actual, expected, atol=ATOL, rtol=0)
        np.testing.assert_allclose(_word_matrix(7, _adjoint(word)), expected.conj().T,
                                   atol=ATOL, rtol=0)
        wrong_inverse = _adjoint(phase) + predicate + subtraction
        witness = (1 << 2) | (1 << 4)  # j=0,k=1,H=1 underflows.
        wrong = _apply_native_word(7, wrong_inverse, actual[:, [witness]])
        self.assertGreater(np.linalg.norm(wrong - np.eye(128)[:, [witness]]), 1)

        # H 0, two private bits 1,2, two history bits 3,4, arbitrary borrow 5.
        history = _history_word(0, [1, 2], [3, 4], 5)
        images = []
        for value in range(64):
            count = (value >> 3) & 3
            tick = (value & 1) and bool(value & 6)
            images.append((value & ~(3 << 3)) | (((count + tick) % 4) << 3))
        expected = np.eye(64, dtype=complex)[:, images]
        np.testing.assert_allclose(_word_matrix(6, history), expected, atol=ATOL, rtol=0)
        np.testing.assert_allclose(_word_matrix(6, _adjoint(history)), expected.conj().T,
                                   atol=ATOL, rtol=0)

    def test_occupied_inactive_native_source_wrapper_uses_actual_inverse(self):
        # Source data/root/scratch 0..3; H 4; digit 5; underflow 6;
        # arbitrary phase helper 7; history 8,9; returned dirty borrow 10.
        preparation = _native_source()
        shift, _, _ = _guarded_shift(0, 5, 6, 4, 10)  # Gray data bit 0 is high.
        phase = _phase_s([4, 5], 7, 10)
        first = [('H', 5)] + shift + phase + [('H', 5)]
        second = [('H', 5)] + _adjoint(phase) + shift + [('H', 5)]
        history = _history_word(4, [5, 6], [8, 9], 10)
        body = first + history + second
        word = preparation + body + _adjoint(preparation)
        rng = np.random.default_rng(7461)
        columns = rng.normal(size=(2048, 2)).astype(complex)
        columns += 1j * rng.normal(size=columns.shape)
        columns[(np.arange(2048) & (1 << 4)) != 0] = 0
        columns /= np.linalg.norm(columns, axis=0)
        actual = _apply_native_word(11, word, columns)
        np.testing.assert_allclose(actual, columns, atol=ATOL, rtol=0)
        np.testing.assert_allclose(_apply_native_word(11, _adjoint(word), actual),
                                   columns, atol=ATOL, rtol=0)
        # A raw Gray root enable of zero is NOT an occupied-scratch guard.
        raw = _native([(gate.op, *gate.qubits) for gate in emit_chain(2)])
        witness = np.eye(16, dtype=complex)[:, [1 << 3]]  # root=0, scratch=1.
        self.assertGreater(np.linalg.norm(_apply_native_word(4, raw, witness) - witness), .1)
        # The completed source wrapper must use its actual inverse, not PREP again.
        wrong = _apply_native_word(11, preparation + body + preparation, columns)
        self.assertGreater(np.linalg.norm(wrong - columns), .1)


if __name__ == '__main__':
    unittest.main()
