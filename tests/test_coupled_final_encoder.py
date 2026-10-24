"""Exact fixtures for a noncommuting final encoder, not a full source.

The analytic kernel is checked against literal native gates on every input
of two small instances, retaining dirty tags and actual inverse phases.
"""
from fractions import Fraction
import unittest

from tests.exact_matrices import (
    CI, CO, CQ2, CZ, Q2, Gate, adjoint, apply_word, multiply, native_columns, zeros,
)
from tests.test_two_layer_diagonal_rigidity import (
    _primitive_inverse, _reduce_adjacent_s,
)
from tests.test_weak_pauli_orbit_compression import _binary_rank


def _encoder(p):
    """G=H_d T(X_d)^p product T(X_ai) T(X_d X_ai Z_ri)."""
    word = [('H', 0)] + [('H', 1 + i) for i in range(p)]
    basis = [gate for i in range(p)
             for gate in [('CX', 0, 1 + p + i),
                          ('CX', 1 + i, 1 + p + i)]]
    word += basis
    word += [('S', 0)] * ((p // 2) % 4)
    if p % 2:
        word += [('T', 0)]
    word += [('T', 1 + i) for i in range(p)]
    word += [('T', 1 + p + i) for i in range(p)]
    word += list(reversed(basis))
    word += [('H', 1 + i) for i in reversed(range(p))]
    return word


def _phase8(exponent):
    omega = CQ2(Q2(0, Fraction(1, 2)), Q2(0, Fraction(1, 2)))
    result = CO
    for _ in range(exponent % 8):
        result *= omega
    return result


def _normalization(p):
    result = CO
    half_root = CQ2(Q2(0, Fraction(1, 2)), Q2(0, 0))
    for _ in range(p):
        result *= half_root
    return result


def _kernel(p, width, residual):
    size, address_mask = 1 << width, (1 << p) - 1
    expected = zeros(size)
    for column in range(size):
        d = column & 1
        address = (column >> 1) & address_mask
        dirty = (column >> (1 + p)) & address_mask
        for output_address in range(1 << p):
            intermediate, phase = 0, 0
            for i in range(p):
                right = ((dirty >> i) & 1) == d
                selected = output_address if right else address
                intermediate |= ((selected >> i) & 1) << i
                change = ((address ^ output_address) >> i) & 1
                phase += change if right else -change
            row = ((column & ~((address_mask << 1) | 1))
                   | (output_address << 1) | (1 - d))
            exponent = (p * (1 if d == 0 else -1) + 2 * phase
                        + 4 * residual(intermediate, column))
            expected[row][column] = _normalization(p) * _phase8(exponent)
    return expected


def _full_bank_wall(p):
    offset, size = p + 1, 1 << p
    return ([('CX', offset + j, 0) for j in range(size)]
            + [('CX', offset + j, 1 + i)
               for j in range(size) for i in range(p) if j >> i & 1])


def _transversal_bank_t(p):
    hadamards = [('H', p + 1 + j) for j in range(1 << p)]
    return hadamards + [('T', p + 1 + j) for j in range(1 << p)] + hadamards


def _full_bank_encoder(p):
    wall = _full_bank_wall(p)
    return wall + _transversal_bank_t(p) + list(reversed(wall)) + [('H', 0)]


def _classical_action(word, column):
    for name, source, target in word:
        assert name == 'CX'
        if column >> source & 1:
            column ^= 1 << target
    return column


def _full_bank_charge(p, column, rewritten):
    size, offset = 1 << p, p + 1
    address = (column >> 1) & (size - 1)
    bank, probe = column >> offset, column & 1
    if rewritten:
        syndrome = 0
        for j in range(size):
            if bank >> j & 1:
                syndrome ^= j
        return probe ^ (bank.bit_count() & 1) ^ ((bank >> (address ^ syndrome)) & 1)
    return probe ^ ((bank >> address) & 1)


def _native_state(word, state):
    return apply_word([Gate(name, tuple(qubits)) for name, *qubits in word], state)


def _add_state(total, state, factor):
    for row, value in state.items():
        total[row] = total.get(row, CZ) + factor * value
        if total[row] == CZ:
            del total[row]


def _full_bank_expansion_column(p, column):
    """Exact Fourier/Pauli formula, independent of the native emitter."""
    size, offset = 1 << p, p + 1
    address, probe = (column >> 1) & (size - 1), column & 1
    result = {}
    half_root = _normalization(1)
    for u in range(size):
        base = sum(1 << k for k in range(size) if (u & k).bit_count() % 2 == 0)
        for j in range(size):
            mask, state = base ^ (1 << j), {column: CO}
            for k in range(size):
                if not (mask >> k & 1):
                    continue
                changed = {}
                for row, value in state.items():
                    sign = -1 if row >> (offset + k) & 1 else 1
                    _add_state(changed, {row: value}, half_root * sign)
                    _add_state(changed, {row ^ (1 << (offset + k)): value},
                               -CI * half_root * sign)
                state = changed
            exponent = probe + (u & address).bit_count() + (u & j).bit_count()
            _add_state(result, state, CQ2(Q2(Fraction((-1) ** exponent, size))))
    return result


class CoupledFinalEncoderTests(unittest.TestCase):
    def test_complete_kernel_and_actual_inverse(self):
        # p=2: all four coefficient-bank bits and both dirty tags are inputs.
        # Wire order d,a,b,r_a,r_b,v,w,t; residual is av+bw+abt.
        def lookup_residual(address, column):
            a, b = address & 1, (address >> 1) & 1
            v, w, t = [(column >> bit) & 1 for bit in (5, 6, 7)]
            return a * v ^ b * w ^ a * b * t

        # p=3: a cubic address sign challenges the hybrid-address dependence.
        cases = ((2, 8, lookup_residual),
                 (3, 7, lambda address, _: int(address == 7)))
        for p, width, residual in cases:
            with self.subTest(p=p):
                word = _encoder(p)
                source = native_columns(width, word)
                inverse = _reduce_adjacent_s(_primitive_inverse(word))
                self.assertEqual(native_columns(width, inverse), adjoint(source))
                target = zeros(1 << width)
                for column in range(1 << width):
                    address = (column >> 1) & ((1 << p) - 1)
                    sign = (column & 1) ^ residual(address, column)
                    target[column][column] = -CO if sign else CO
                actual = multiply(multiply(source, target), adjoint(source))
                expected = _kernel(p, width, residual)
                self.assertEqual(actual, expected)
                for column in range(1 << width):
                    entries = [row[column] for row in actual
                               if row[column] != CZ]
                    self.assertEqual(len(entries), 1 << p)
                    self.assertTrue(all(value * value.conj()
                                        == _normalization(2 * p)
                                        for value in entries))

    def test_physical_axis_independence_and_common_native_ledger(self):
        for p in range(1, 9):
            with self.subTest(p=p):
                width = 1 + 2 * p
                axes = [1] if p % 2 else []
                axes += [1 << (1 + i) for i in range(p)]
                axes += [1 ^ (1 << (1 + i))
                         ^ (1 << (width + 1 + p + i)) for i in range(p)]
                self.assertEqual(_binary_rank(axes), len(axes))
                k = (p // 2) % 4
                word = _encoder(p)
                inverse = _reduce_adjacent_s(_primitive_inverse(word))
                for native, expected_cliffords in (
                        (word, 6 * p + 1 + k),
                        (inverse, 6 * p + 1 + ((-k) % 4))):
                    odd = [index for index, gate in enumerate(native)
                           if gate[0] in ('T', 'TDG')]
                    self.assertEqual(len(odd), 2 * p + p % 2)
                    self.assertEqual(odd, list(range(odd[0], odd[-1] + 1)))
                    self.assertEqual(len({native[index][1] for index in odd}),
                                     len(odd))
                    cliffords = len(native) - len(odd)
                    self.assertEqual(cliffords, expected_cliffords)
                    self.assertLessEqual(cliffords, 6 * p + 4)
                    self.assertEqual({wire for gate in native for wire in gate[1:]},
                                     set(range(width)))

    def test_full_bank_charge_rewrite_and_native_ledger(self):
        for p in (1, 2, 3):
            width = 1 + p + (1 << p)
            wall = _full_bank_wall(p)
            for column in range(1 << width):
                changed = _classical_action(wall, column)
                self.assertEqual(_classical_action(wall, changed), column)
                self.assertEqual(_full_bank_charge(p, changed, False),
                                 _full_bank_charge(p, column, True))
        for p in range(1, 7):
            size = 1 << p
            word = _full_bank_encoder(p)
            for native in (word, _primitive_inverse(word)):
                odd = [index for index, gate in enumerate(native)
                       if gate[0] in ('T', 'TDG')]
                self.assertEqual(len(odd), size)
                self.assertEqual(odd, list(range(odd[0], odd[-1] + 1)))
                self.assertEqual(len({native[index][1] for index in odd}), size)
                self.assertEqual(len(native) - len(odd), p * size + 4 * size + 1)
                self.assertEqual({wire for gate in native for wire in gate[1:]},
                                 set(range(1 + p + size)))

    def test_full_bank_actual_inverse_and_complete_pauli_expansion(self):
        for p in (1, 2):
            width = 1 + p + (1 << p)
            source = _full_bank_encoder(p)
            actual = native_columns(width, source)
            self.assertEqual(native_columns(width, _primitive_inverse(source)),
                             adjoint(actual))
            # The complete operator identity has already retained every
            # bank variable; these few columns challenge its signs and factors.
            transversal = _transversal_bank_t(p)
            for column in (0, 1, (1 << width) - 1):
                state = _native_state(_primitive_inverse(transversal), {column: CO})
                state = {row: (-value if _full_bank_charge(p, row, True) else value)
                         for row, value in state.items()}
                state = _native_state(transversal, state)
                self.assertEqual(state, _full_bank_expansion_column(p, column))


if __name__ == '__main__':
    unittest.main()
