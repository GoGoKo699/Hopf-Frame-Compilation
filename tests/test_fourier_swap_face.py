"""Exact nonmonomial face, native return and held-address four-corner screen.

This fixed target has top coefficient one. No all-size source or variable
top coefficient is inferred from its complete operator identity.
"""
from collections import Counter
import unittest

from tests.exact_matrices import CO, CZ, adjoint, eye, multiply, native_columns
from tests.test_coherent_feedback_sources import _shared_ccz
from tests.test_common_probe_source import _native
from tests.test_dirty_cubic_phase_compiler import _gray_phase_batch
from tests.test_joint_address_helper_source import _integer_word
from tests.test_joint_lookup_pauli_source import _bands, _inverse


X1, X2, X3, X4, ROOT, U, V = range(7)
MASKS = (7, 8)


def _bit(value, wire):
    return (value >> wire) & 1


def _parts():
    bell = [('H', U), ('CX', U, V)]
    first = ([('H', ROOT)] + _shared_ccz(ROOT, [(U, V)])
             + [('H', ROOT)] + bell + [('H', V)])
    d2 = _gray_phase_batch([(U, X1, X2), (V, X3, X4)], MASKS, [1, 1], 3)
    middle = [('H', V)] + _inverse(bell) + [('H', ROOT)]
    d3 = _gray_phase_batch([(ROOT, U, V)], MASKS[:1], [1], 3)
    return first, d2, middle, d3, [('H', ROOT)]


def _ideal_probe():
    _, _, middle, _, final = _parts()
    swap = [('CX', U, V), ('CX', V, U), ('CX', U, V)]
    k = [('Z', ROOT), ('H', V)] + swap + [('H', V)]
    d2 = [('CCZ', U, X1, X2), ('CCZ', V, X3, X4)]
    outer = d2 + middle + [('CCZ', ROOT, U, V)] + final
    return list(reversed(outer)) + k + outer


def _target_sign(value):
    a = _bit(value, X1) * _bit(value, X2)
    b = _bit(value, X3) * _bit(value, X4)
    return (-1) ** (_bit(value, ROOT) + a * _bit(value, V)
                   + b * _bit(value, U) + a * b)


def _direct_words():
    def toffoli(left, right, target):
        return [('H', target)] + _shared_ccz(target, [(left, right)]) + [('H', target)]
    m = toffoli(U, V, ROOT)
    left, right = toffoli(X1, X2, U), toffoli(X3, X4, V)
    # Their identically shaped native words occupy disjoint triples.
    t = [gate for pair in zip(left, right) for gate in pair]
    cz = [('H', V), ('CX', U, V), ('H', V)]
    return m + t + m, cz + t + cz + t + [('S', ROOT)] * 2


class FourierSwapFaceTests(unittest.TestCase):
    def test_complete_seven_wire_identity_and_actual_inverse(self):
        word = _ideal_probe()
        for basis in range(1 << 7):
            for order in (word, list(reversed(word))):
                output, denominator = _integer_word(order, basis)
                self.assertEqual(output, {basis: _target_sign(basis) * denominator})

    def test_nonmonomial_balanced_clifford_and_native_first_source(self):
        rename = {ROOT: 0, U: 1, V: 2}
        first = [(gate[0], *(rename[q] for q in gate[1:])) for gate in _parts()[0]]
        source = native_columns(3, first)
        inverse = native_columns(3, _inverse(first))
        swap = [('CX', 1, 2), ('CX', 2, 1), ('CX', 1, 2)]
        target = native_columns(3, [('Z', 0), ('H', 2)] + swap + [('H', 2)])
        probe = native_columns(3, [('Z', 0)])
        self.assertEqual(multiply(multiply(source, probe), inverse), target)
        self.assertEqual(inverse, adjoint(source))
        self.assertEqual(multiply(inverse, source), eye(8))
        self.assertEqual(target, adjoint(target))
        self.assertEqual(sum((target[i][i] for i in range(8)), CZ), CZ)
        self.assertTrue(all(sum(value != CZ for value in row) == 4 for row in target))
        pair = native_columns(2, [('H', 1), ('CX', 0, 1), ('CX', 1, 0),
                                  ('CX', 0, 1), ('H', 1)])
        for mask in (1, 2, 3):
            z = native_columns(2, [('Z', wire) for wire in range(2) if mask >> wire & 1])
            image = multiply(multiply(pair, z), pair)
            self.assertTrue(any(image[i][j] != CZ for i in range(4) for j in range(4) if i != j))

    def test_literal_cubic_batches_return_every_arbitrary_mask(self):
        _, d2, _, d3, _ = _parts()
        for basis in range(1 << 9):
            a = _bit(basis, X1) * _bit(basis, X2)
            b = _bit(basis, X3) * _bit(basis, X4)
            sign = (-1) ** (a * _bit(basis, U) + b * _bit(basis, V))
            for word in (d2, _inverse(d2)):
                self.assertEqual(_native(word, {basis: CO}), {basis: sign * CO})
        rename = {ROOT: 0, U: 1, V: 2, MASKS[0]: 3}
        phase = [(gate[0], *(rename[q] for q in gate[1:])) for gate in d3]
        for basis in range(16):
            sign = -1 if basis & 7 == 7 else 1
            for word in (phase, _inverse(phase)):
                self.assertEqual(_native(word, {basis: CO}), {basis: sign * CO})

    def test_common_forward_inverse_native_ledger(self):
        source = [gate for part in _parts() for gate in part]
        for word in (source, _inverse(source)):
            count = Counter(gate[0] for gate in word)
            self.assertEqual(count['T'] + count['TDG'], 31)
            self.assertEqual(len(word), 75)
            self.assertEqual(len(_bands(word)), 20)
            for band in _bands(word):
                self.assertEqual(len({gate[1] for gate in band}), len(band))
        probe = _inverse(source) + [('S', ROOT)] * 2 + source
        self.assertEqual(len(probe), 152)
        self.assertEqual(sum(g[0] in ('T', 'TDG') for g in probe), 62)
        self.assertEqual(len(_bands(probe)), 40)

    def test_direct_classical_normal_form_and_complete_native_phase(self):
        direct, phase = _direct_words()
        self.assertEqual((len(direct), sum(g[0] in ('T', 'TDG') for g in direct),
                          len(_bands(direct))), (76, 28, 12))
        self.assertEqual((len(phase), sum(g[0] in ('T', 'TDG') for g in phase),
                          len(_bands(phase))), (84, 28, 8))
        self.assertTrue(all(q < 7 for word in (direct, phase) for gate in word for q in gate[1:]))
        for word in (direct, phase):
            for band in _bands(word):
                self.assertEqual(len({gate[1] for gate in band}), len(band))
        for basis in range(128):
            a = _bit(basis, X1) * _bit(basis, X2)
            b = _bit(basis, X3) * _bit(basis, X4)
            root_change = a * _bit(basis, V) ^ b * _bit(basis, U) ^ a * b
            row = basis ^ (a << U) ^ (b << V) ^ (root_change << ROOT)
            for word in (direct, _inverse(direct)):
                self.assertEqual(_native(word, {basis: CO}), {row: CO})
            for word in (phase, _inverse(phase)):
                self.assertEqual(_native(word, {basis: CO}), {basis: _target_sign(basis) * CO})

    def test_four_corner_scalar_and_variable_coefficient_contrast(self):
        swap = [('CX', 1, 2), ('CX', 2, 1), ('CX', 1, 2)]
        words = [([('Z', 0), ('H', 2)] + swap + [('H', 2)], -1),
                 ([('Z', 0), ('H', 1), ('H', 2)], 1)]
        p = native_columns(3, [('Z', 1)])
        r = native_columns(3, [('Z', 2)])
        for word, sign in words:
            k = native_columns(3, word)
            result = eye(8)
            for a, b in ((0, 0), (1, 0), (1, 1), (0, 1)):
                d = eye(8)
                if a:
                    d = multiply(p, d)
                if b:
                    d = multiply(r, d)
                result = multiply(result, multiply(multiply(d, k), d))
            self.assertEqual(result, [[sign * entry for entry in row] for row in eye(8)])
        for basis in range(16):
            root, c1, c2, d = (_bit(basis, wire) for wire in range(4))
            value = 1
            for a, b in ((0, 0), (1, 0), (1, 1), (0, 1)):
                value *= (-1) ** (root + a * c1 + b * c2 + a * b * d)
            self.assertEqual(value, (-1) ** d)


if __name__ == '__main__':
    unittest.main()
