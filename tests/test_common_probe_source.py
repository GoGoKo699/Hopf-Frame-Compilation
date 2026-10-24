"""Exact checks for common-probe pulses and their complete lookup defects.

These fixtures challenge native phases, fiber orientation, and norm witnesses.
The all-size proofs and their restricted circuit grammars live in
research/depth/COUPLED_FINAL_ENCODER.md, Section 8.
"""
from collections import Counter
from fractions import Fraction
import unittest

from tests.exact_matrices import (
    CI, CO, CZ, CQ2, Q2, Gate, adjoint, apply_word, multiply, native_columns,
)
from tests.test_joint_lookup_pauli_source import _bands, _inverse


ROOT_HALF = CQ2(Q2(0, Fraction(1, 2)))
OMEGA = ROOT_HALF * (CO + CI)
PHASES = [CO]
for _ in range(7):
    PHASES.append(PHASES[-1] * OMEGA)


def _native(word, state):
    out = apply_word([Gate(g[0], tuple(g[1:])) for g in word], state)
    return {key: value for key, value in out.items() if value != CZ}


def _parts(p):
    m, d = 1 << p, p + (1 << p)
    bank = range(p, p + m)
    wall = [('CX', p + j, d) for j in range(m)]
    wall += [('CX', p + j, i) for j in range(m) for i in range(p)
             if (j >> i) & 1]
    h = [('H', wire) for wire in bank]
    a = wall + h + [('T', wire) for wire in bank] + h + list(reversed(wall))
    parity = [('CX', d, p + j) for j in range(m)]
    parity += [('CX', i, p + j) for j in range(m) for i in range(p)
               if (j >> i) & 1]
    y_basis = [('H', d), ('S', d)]
    b = (_inverse(y_basis) + parity + [('T', wire) for wire in bank]
         + list(reversed(parity)) + y_basis)
    parity = [('CX', p + j, p) for j in range(1, m)]
    c = ([('S', wire) for wire in bank] + parity + [('S', p)] * 3
         + list(reversed(parity)))
    e = b + a + c + _inverse(b) + _inverse(a) + _inverse(c)
    return wall, a, b, c, e


def _syndrome(bank, p):
    result = 0
    for row in range(1 << p):
        if (bank >> row) & 1:
            result ^= row
    return result


def _m_map(bank, p):
    syndrome, parity = _syndrome(bank, p), bank.bit_count() & 1
    return bank ^ sum((parity ^ ((syndrome & j).bit_count() & 1)) << j
                      for j in range(1 << p))


def _d_phase(bank, p, address):
    m = 1 << p
    tag = sum(((address & j).bit_count() & 1) << j for j in range(m))
    return PHASES[((_m_map(bank, p) ^ tag).bit_count() - m // 2) % 8]


def _diagonal(state, phase):
    return {key: value * phase(key) for key, value in state.items()
            if value != CZ}


def _fiber(state, p, address, eta, reverse=False):
    m = 1 << p
    t = ([('H', j) for j in range(m)] + [('T', j) for j in range(m)]
         + [('H', j) for j in range(m)])
    operations = [('D', eta), ('T', 1), ('C', 1),
                  ('D', -eta), ('T', -1), ('C', -1)]
    if reverse:
        operations = [(kind, -sign) for kind, sign in reversed(operations)]
    for kind, sign in operations:
        if kind == 'T':
            state = _native(t if sign == 1 else _inverse(t), state)
        elif kind == 'C':
            state = _diagonal(state, lambda y:
                              -CO if (y.bit_count() * (y.bit_count() - 1) // 2) & 1
                              else CO)
        else:
            state = _diagonal(state, lambda y:
                              _d_phase(y, p, address) if sign == 1
                              else _d_phase(y, p, address).conj())
    return state


def _singlets(p):
    m, shift = 1 << p, 3
    state = {0: ROOT_HALF, shift: ROOT_HALF}
    for j in range(m):
        k = j ^ shift
        if j > k:
            continue
        updated = {}
        for column, value in state.items():
            updated[column | (1 << (p + k))] = value * ROOT_HALF
            updated[column | (1 << (p + j))] = -value * ROOT_HALF
        state = updated
    return state


class CommonProbeSourceTests(unittest.TestCase):
    def test_native_columns_actual_inverse_and_complete_fiber_orientation(self):
        p, m, d = 1, 2, 3
        wall, a, b, c, word = _parts(p)
        actual = native_columns(d + 1, word)
        self.assertEqual(native_columns(d + 1, _inverse(word)), adjoint(actual))
        translation = []
        for wire in (0, p + 1):
            translation += [('H', wire), ('S', wire), ('S', wire), ('H', wire)]
        symmetry = native_columns(d + 1, translation)
        for pulse in (a, b):
            matrix = native_columns(d + 1, pulse)
            self.assertEqual(multiply(matrix, symmetry), multiply(symmetry, matrix))
        graph = native_columns(d + 1, c)
        self.assertNotEqual(multiply(graph, symmetry), multiply(symmetry, graph))
        echo = wall + _inverse(word) + [('S', d)] * 2 + word + list(reversed(wall))
        for address in range(m):
            for bank in range(1 << m):
                column = address | (bank << p)
                # Input minus-Y, output plus-Y: the off-diagonal block is
                # E_+ Z_parity E_-†, not its adjoint.
                state = {column: ROOT_HALF, column | (1 << d): -CI * ROOT_HALF}
                observed = _native(echo, state)
                fiber = _fiber({bank: CO}, p, address, -1, reverse=True)
                fiber = _diagonal(fiber, lambda y: -CO if y.bit_count() & 1 else CO)
                fiber = _fiber(fiber, p, address, 1)
                expected = {}
                for y, value in fiber.items():
                    row = address | (y << p)
                    expected[row] = value * ROOT_HALF
                    expected[row | (1 << d)] = CI * value * ROOT_HALF
                self.assertEqual(observed, {key: value for key, value in expected.items()
                                            if value != CZ})

    def test_complete_native_singlet_witness_and_literal_phases(self):
        for p in (2, 3):
            m, d = 1 << p, p + (1 << p)
            _, a, b, c, word = _parts(p)
            state = _singlets(p)
            self.assertEqual(sum((value * value.conj() for value in state.values()), CZ), CO)
            native_phase = CI if p == 2 else -CO
            for pulse in (a, b):
                self.assertEqual(_native(pulse, state),
                                 {key: value * native_phase for key, value in state.items()})
            c_phase = -CO if (m // 2 * (m // 2 - 1) // 2) & 1 else CO
            self.assertEqual(_native(c, state),
                             {key: value * c_phase for key, value in state.items()})
            self.assertEqual(_native(word, state), state)
            self.assertEqual(_native(_inverse(word), state), state)
            observed = _native(_inverse(word) + [('S', d)] * 2 + word, state)
            self.assertEqual(observed, state)
            target = {key: -value if (key >> (p + (key & (m - 1)))) & 1 else value
                      for key, value in state.items()}
            self.assertEqual(sum((state[key].conj() * value for key, value in target.items()), CZ), CZ)
            self.assertEqual(sum(((observed[key] - value) * (observed[key] - value).conj()
                                  for key, value in target.items()), CZ), 2 * CO)

    def test_complete_p3_square_and_both_local_spectral_factors(self):
        p, m = 3, 8
        for address in range(m):
            mask = sum((1 ^ (((address ^ (m - 1)) & j).bit_count() & 1)) << j
                       for j in range(m))
            for bank in range(1 << m):
                self.assertEqual(_m_map(_m_map(bank, p), p), bank)
                phase = _d_phase(bank, p, address)
                expected = PHASES[(2 * bank.bit_count() + 4 * ((mask & bank).bit_count() & 1)) % 8]
                self.assertEqual(phase * phase, expected)
        t = [('H', 0), ('T', 0), ('H', 0)]
        for orientation in (1, 3):
            r = [('S', 0)] * orientation
            g = native_columns(1, r + _inverse(t) + r + t)
            squared = multiply(g, g)
            self.assertEqual((squared[0][0] + squared[1][1]) / 2,
                             CQ2(Q2(Fraction(1, 4), Fraction(1, 2))))
            self.assertEqual(squared[0][0] * squared[1][1]
                             - squared[0][1] * squared[1][0], CO)

    def test_distinct_targets_and_native_ledgers(self):
        for p in (1, 2, 3, 6):
            m = 1 << p
            _, a, b, c, word = _parts(p)
            for circuit in (word, _inverse(word)):
                bands = _bands(circuit)
                self.assertEqual(list(map(len, bands)), [m] * 4)
                self.assertTrue(all(len({gate[1] for gate in band}) == m for band in bands))
                counts = Counter(gate[0] for gate in circuit)
                self.assertEqual(counts['T'] + counts['TDG'], 4 * m)
                self.assertLessEqual(counts['H'] + counts['S'] + counts['CX'],
                                     4 * p * m + 20 * m + 12)
            self.assertEqual(sum(g[0] != 'T' for g in a), p * m + 4 * m)
            self.assertEqual(sum(g[0] != 'T' for g in b), p * m + 2 * m + 6)
            self.assertEqual(len(c), 3 * m + 1)
            self.assertEqual(len(_inverse(c)), 5 * m - 1)


if __name__ == '__main__':
    unittest.main()
