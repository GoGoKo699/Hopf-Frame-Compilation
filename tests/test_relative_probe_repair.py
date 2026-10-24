"""Exact fixtures for the corrected relative probe and complete lookup defects.

The native words, actual adjoints and complete witnesses are checked here.
The all-size fiber and displacement proofs live in COUPLED_FINAL_ENCODER §9.
"""
from collections import Counter
from fractions import Fraction
from math import comb
import unittest

from tests.exact_matrices import CI, CO, CZ, adjoint, native_columns
from tests.test_common_probe_source import (
    PHASES, ROOT_HALF, _d_phase, _diagonal, _native, _parts, _syndrome,
)
from tests.test_joint_lookup_pauli_source import _bands, _inverse


def _quadratic_mask(p):
    """Paired-coordinate quadratic, plus the last coordinate when p is odd."""
    selected = []
    for row in range(1 << p):
        value = sum(((row >> i) & 1) * ((row >> (i + 1)) & 1)
                    for i in range(0, p - 1, 2)) & 1
        if p & 1:
            value ^= (row >> (p - 1)) & 1
        if value:
            selected.append(row)
    return selected


def _repair_parts(p, selected=()):
    m = 1 << p
    wall, a, _, _, _ = _parts(p)
    gather = [('CX', p + j, p) for j in range(1, m)]
    spread = [('CX', p, p + j) for j in range(1, m)]
    parity_flip = gather + spread + gather
    tags = [('CX', i, p + j) for j in range(m) for i in range(p)
            if (j >> i) & 1]
    exposure = parity_flip + tags
    b = exposure + [('T', p + j) for j in range(m)] + _inverse(exposure)
    c = (wall + [('H', p + j) for j in selected] + _inverse(wall)
         if selected else [])
    word = b + a + c + _inverse(b) + _inverse(a) + _inverse(c)
    return wall, a, b, c, word


def _shared_fiber(state, p, address, selected=(), reverse=False):
    m = 1 << p
    t = ([('H', j) for j in range(m)] + [('T', j) for j in range(m)]
         + [('H', j) for j in range(m)])
    c = [('H', j) for j in selected]
    operations = [('D', 1), ('T', 1), ('C', 1),
                  ('D', -1), ('T', -1), ('C', -1)]
    if reverse:
        operations = [(kind, -sign) for kind, sign in reversed(operations)]
    for kind, sign in operations:
        if kind == 'T':
            state = _native(t if sign == 1 else _inverse(t), state)
        elif kind == 'C':
            state = _native(c, state)
        else:
            state = _diagonal(state, lambda y:
                              _d_phase(y, p, address) if sign == 1
                              else _d_phase(y, p, address).conj())
    return state


def _relative_image(state, p, address, selected=()):
    state = _shared_fiber(state, p, address, selected, reverse=True)
    state = _diagonal(state, lambda y: -CO if y.bit_count() & 1 else CO)
    return _shared_fiber(state, p, address, selected)


def _probe_state(bank_state, p, address, eta):
    d = p + (1 << p)
    state = {}
    for bank, value in bank_state.items():
        row = address | (bank << p)
        state[row] = value * ROOT_HALF
        state[row | (1 << d)] = eta * CI * value * ROOT_HALF
    return {key: value for key, value in state.items() if value != CZ}


def _q_sign(bank, p, address):
    selected = address ^ _syndrome(bank, p)
    return -1 if (bank.bit_count() + ((bank >> selected) & 1)) & 1 else 1


def _binomial_tail(n, k):
    return Fraction(sum(comb(n, j) for j in range(k, n + 1)), 1 << n)


def _cubic_face_phase(bank, p):
    exponent = sum(((bank >> j) & 1) * ((bank >> (j ^ 1)) & 1)
                   * ((bank >> (j ^ 2)) & 1) for j in range(1 << p))
    return -CO if exponent & 1 else CO


def _ccz_parity_word(a, b, c):
    """Literal seven-T, ten-CNOT phase polynomial; no initialized work."""
    return [('T', a), ('T', b), ('T', c),
            ('CX', a, b), ('TDG', b), ('CX', a, b),
            ('CX', a, c), ('TDG', c), ('CX', a, c),
            ('CX', b, c), ('TDG', c), ('CX', b, c),
            ('CX', a, c), ('CX', b, c), ('T', c),
            ('CX', b, c), ('CX', a, c)]


def _cubic_face_word(p, offset=0):
    word = []
    for local in range(4):
        faces = [_ccz_parity_word(offset + start + local,
                                 offset + start + (local ^ 1),
                                 offset + start + (local ^ 2))
                 for start in range(0, 1 << p, 4)]
        # Corresponding gates on distinct faces have disjoint wires.
        for position in range(len(faces[0])):
            word.extend(face[position] for face in faces)
    return word


def _cubic_shared_fiber(state, p, address, reverse=False):
    m = 1 << p
    h = [('H', j) for j in range(m)]
    t = h + [('T', j) for j in range(m)] + h
    operations = [('D', 1), ('T', 1), ('C', 1),
                  ('D', -1), ('T', -1), ('C', -1)]
    if reverse:
        operations = [(kind, -sign) for kind, sign in reversed(operations)]
    for kind, sign in operations:
        if kind == 'T':
            state = _native(t if sign == 1 else _inverse(t), state)
        elif kind == 'D':
            state = _diagonal(state, lambda y:
                              _d_phase(y, p, address) if sign == 1
                              else _d_phase(y, p, address).conj())
        elif sign == 1:
            state = _native(h, state)
            state = _diagonal(state, lambda y: _cubic_face_phase(y, p))
        else:
            state = _diagonal(state, lambda y: _cubic_face_phase(y, p))
            state = _native(h, state)
    return state


def _cubic_relative_image(state, p, address):
    state = _cubic_shared_fiber(state, p, address, reverse=True)
    state = _diagonal(state, lambda y: -CO if y.bit_count() & 1 else CO)
    return _cubic_shared_fiber(state, p, address)


class RelativeProbeRepairTests(unittest.TestCase):
    def test_cubic_hadamard_source_native_columns_and_pauli_collapse(self):
        face = _cubic_face_word(2)
        for bank in range(16):
            expected = {bank: _cubic_face_phase(bank, 2)}
            self.assertEqual(_native(face, {bank: CO}), expected)
            self.assertEqual(_native(_inverse(face), expected), {bank: CO})
            translated = _cubic_face_phase(bank ^ 15, 2)
            self.assertEqual(translated * _cubic_face_phase(bank, 2),
                             -CO if bank.bit_count() & 1 else CO)
        for p in (2, 3):
            m, d, address = 1 << p, p + (1 << p), (1 << p) - 1
            wall, a, b, _, _ = _repair_parts(p)
            c = (wall + [('H', p + j) for j in range(m)]
                 + _cubic_face_word(p, offset=p) + _inverse(wall))
            word = b + a + c + _inverse(b) + _inverse(a) + _inverse(c)
            for circuit in (word, _inverse(word)):
                self.assertEqual(sum(g[0] in ('T', 'TDG') for g in circuit), 18 * m)
                self.assertEqual(sum(g[0] not in ('T', 'TDG') for g in circuit),
                                 6 * p * m + 46 * m - 12)
            relative = _cubic_relative_image({0: CO}, p, address)
            if p == 2:
                self.assertEqual(relative[0], Fraction(1, 8) * CO)
            else:
                self.assertEqual(relative, {(1 << m) - 1: CO})
                for bank in (1, 0xAA):
                    self.assertEqual(
                        _cubic_relative_image({bank: CO}, p, address),
                        {bank ^ ((1 << m) - 1):
                         -CO if bank.bit_count() & 1 else CO})
            initial = _probe_state({0: CO}, p, address, -1)
            echo = wall + _inverse(word) + [('S', d)] * 2 + word + _inverse(wall)
            self.assertEqual(_native(echo, initial),
                             _probe_state(relative, p, address, 1))

    def test_native_source_inverse_middle_phase_and_relative_involution(self):
        p, d = 1, 3
        _, _, _, _, word = _repair_parts(p)
        self.assertEqual(native_columns(d + 1, _inverse(word)),
                         adjoint(native_columns(d + 1, word)))
        # Literal middle scalar is e^(im pi/8), including at p=1 and p=2.
        for p in (1, 2, 3):
            m, d = 1 << p, p + (1 << p)
            _, _, b, _, _ = _repair_parts(p)
            for address, bank, probe in ((0, 0, 0), (1, 1, 1), (m - 1, 5 % (1 << m), 1)):
                column = address | (bank << p) | (probe << d)
                s = address ^ _syndrome(bank, p)
                phase = PHASES[(m // 2) % 8] * _d_phase(bank, p, s)
                self.assertEqual(_native(b, {column: CO}), {column: phase})
        p, m = 3, 8
        d = p + m
        wall, _, _, _, word = _repair_parts(p)
        echo = wall + _inverse(word) + [('S', d)] * 2 + word + _inverse(wall)
        for address, bank in ((0, 0), (1, 1)):
            initial = _probe_state({bank: CO}, p, address, -1)
            relative = _relative_image({bank: CO}, p, address)
            self.assertEqual(_native(echo, initial),
                             _probe_state(relative, p, address, 1))
            self.assertEqual(_relative_image(relative, p, address), {bank: CO})
            self.assertEqual(len(relative), 1 << m)

    def test_nonempty_hadamard_mask_native_complete_witness(self):
        cases = ((1, [0]), (2, _quadratic_mask(2)), (3, _quadratic_mask(3)))
        for p, selected in cases:
            m, d = 1 << p, p + (1 << p)
            q = selected[0]
            bank_state = {0: ROOT_HALF, 1 << q: ROOT_HALF * PHASES[7]}
            self.assertEqual(sum((v * v.conj() for v in bank_state.values()), CZ), CO)
            self.assertTrue(all(_q_sign(y, p, 0) == 1 for y in bank_state))
            relative = _relative_image(bank_state, p, 0, selected)
            self.assertEqual(sum((bank_state.get(y, CZ).conj() * v
                                  for y, v in relative.items()), CZ), CZ)
            wall, _, _, _, word = _repair_parts(p, selected)
            self.assertEqual(_native(_inverse(word), _native(word, {0: CO})), {0: CO})
            echo = wall + _inverse(word) + [('S', d)] * 2 + word + _inverse(wall)
            actual = _native(echo, _probe_state(bank_state, p, 0, -1))
            self.assertEqual(actual, _probe_state(relative, p, 0, 1))
            target = _probe_state(bank_state, p, 0, 1)
            self.assertEqual(sum((actual.get(key, CZ).conj() * value
                                  for key, value in target.items()), CZ), CZ)
            difference = {key: actual.get(key, CZ) - target.get(key, CZ)
                          for key in set(actual) | set(target)}
            self.assertEqual(sum((v * v.conj() for v in difference.values()), CZ), 2 * CO)

    def test_walsh_masks_and_four_rational_displacement_cases(self):
        p, m = 3, 8
        full_mask = (1 << m) - 1
        rows = []
        for t in range(m):
            h_t = sum(((t & j).bit_count() & 1) << j for j in range(m))
            rows.extend((a, t, full_mask ^ h_t ^ (1 << a)) for a in range(m))
        self.assertEqual(len({mask for _, _, mask in rows}), m * m)
        self.assertEqual(Counter(mask.bit_count() for _, _, mask in rows),
                         Counter({m - 1: m, m // 2 - 1: m * (m - 1) // 2,
                                  m // 2 + 1: m * (m - 1) // 2}))
        for address in (0, 1):
            for bank in range(1 << m):
                total = sum((-1) ** (((t & (address ^ a)).bit_count()
                                     + (mask & bank).bit_count()) & 1)
                            for a, t, mask in rows)
                self.assertEqual(total, m * _q_sign(bank, p, address))
        for m, k in ((8, 4), (16, 7), (32, 12), (64, 24)):
            p_h = _binomial_tail(m, k)
            p_r = (Fraction(1, m) * _binomial_tail(m - 1, k)
                   + Fraction(m - 1, 2 * m)
                   * (_binomial_tail(m // 2 - 1, k)
                      + _binomial_tail(m // 2 + 1, k)))
            self.assertGreater(p_h, Fraction(3, 5))
            self.assertLess(p_r, Fraction(3, 20))
            if m == 8:
                self.assertEqual((p_h, p_r), (Fraction(163, 256), Fraction(37, 256)))
            # An independent exact comparison of sqrt(p_h)-sqrt(p_r)>1/3.
            gap = p_h - p_r - Fraction(1, 9)
            self.assertGreater(gap, 0)
            self.assertGreater(gap * gap, Fraction(4, 9) * p_r)

    def test_distinct_native_targets_and_simultaneous_ledgers(self):
        for p in (1, 2, 3, 6):
            m = 1 << p
            for selected in ([], _quadratic_mask(p)):
                _, a, b, c, word = _repair_parts(p, selected)
                expected = (6 * p * m + 24 * m + 2 * len(selected) - 12
                            if selected else 4 * p * m + 20 * m - 12)
                for circuit in (word, _inverse(word)):
                    bands = _bands(circuit)
                    self.assertEqual(list(map(len, bands)), [m] * 4)
                    self.assertTrue(all(len({g[1] for g in band}) == m for band in bands))
                    self.assertEqual(sum(g[0] in ('T', 'TDG') for g in circuit), 4 * m)
                    self.assertEqual(sum(g[0] not in ('T', 'TDG') for g in circuit), expected)
                self.assertEqual(sum(g[0] != 'T' for g in a), p * m + 4 * m)
                self.assertEqual(sum(g[0] != 'T' for g in b), p * m + 6 * m - 6)
                if selected:
                    self.assertEqual(len(c), p * m + 2 * m + len(selected))


if __name__ == '__main__':
    unittest.main()
