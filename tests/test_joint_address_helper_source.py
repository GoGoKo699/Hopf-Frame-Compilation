"""Literal joint-mixing trial, its Pauli-support gap, and finite-face repair.

All-size shallow synthesis is not inferred from this fixed four-address word.
"""
from collections import Counter
from fractions import Fraction
import unittest

from tests.exact_matrices import CI, CO
from tests.test_coherent_feedback_sources import _shared_ccz
from tests.test_common_probe_source import ROOT_HALF, _native
from tests.test_dirty_cubic_phase_compiler import _gray_phase_batch
from tests.test_joint_lookup_pauli_source import _bands, _inverse


X1, X2, X3, X4, ROOT, TOP, C1, C2, U1, U2 = range(10)
MASKS = (10, 11)


def _bit(value, wire):
    return (value >> wire) & 1


def _swap(left, right):
    return [('CX', left, right), ('CX', right, left), ('CX', left, right)]


def _bell(u, c):
    return [('H', u), ('CX', u, c)]


def _walls():
    bell = _bell(U1, C1) + _bell(U2, C2)
    middle = _inverse(bell) + [('H', X2), ('H', X4)]
    final = _swap(X2, C1) + _swap(X4, C2) + [('H', C1), ('H', C2)]
    return bell, middle, final


def _cubic(terms):
    return _gray_phase_batch(terms, MASKS[:len(terms)], [1] * len(terms), 3)


def _toffoli_batch(terms):
    targets = sorted({term[2] for term in terms})
    wall = [('H', target) for target in targets]
    return wall + _cubic(terms) + list(reversed(wall))


def _words():
    bell, middle, final = _walls()
    first = [('H', ROOT)] + _shared_ccz(ROOT, [(U1, C1), (U2, C2)])
    first += [('H', ROOT)] + bell
    source = first + _cubic([(U1, X1, X2), (U2, X3, X4)])
    source += middle + _cubic([(TOP, X2, X4)]) + final
    probe = _inverse(source) + [('S', ROOT), ('S', ROOT)] + source
    background = [('H', X2), ('CX', U1, X2), ('H', X2),
                  ('H', X4), ('CX', U2, X4), ('H', X4)]
    u = _toffoli_batch([(X3, X4, U1), (X1, X2, U2)])
    w = _toffoli_batch([(TOP, U1, C1), (TOP, U2, C2)])
    correction = u + w + _inverse(u) + _inverse(w)
    return first, source, probe, correction, probe + background + correction


def _target_sign(value):
    a = _bit(value, X1) * _bit(value, X2)
    b = _bit(value, X3) * _bit(value, X4)
    exponent = (_bit(value, ROOT) + a * _bit(value, C1) + b * _bit(value, C2)
                + a * b * _bit(value, TOP))
    return (-1) ** exponent


def _trial_action(value):
    a = _bit(value, X1) * _bit(value, X2)
    b = _bit(value, X3) * _bit(value, X4)
    d = _bit(value, TOP)
    output = value ^ ((d * b) << C1) ^ ((d * a) << C2)
    spectator = _bit(value, U1) * _bit(value, X2) + _bit(value, U2) * _bit(value, X4)
    return output, _target_sign(value) * (-1) ** spectator


def _integer_word(word, value):
    """Exact ideal columns with the common Hadamard denominator omitted."""
    state, hadamards = {value: 1}, 0
    for gate in word:
        result = {}
        name, wires = gate[0], gate[1:]
        for basis, amplitude in state.items():
            if name == 'H':
                wire, = wires
                zero = basis & ~(1 << wire)
                for output, sign in ((zero, 1), (zero | (1 << wire), (-1) ** _bit(basis, wire))):
                    result[output] = result.get(output, 0) + sign * amplitude
            else:
                output, sign = basis, 1
                if name == 'CX':
                    control, target = wires
                    output ^= _bit(basis, control) << target
                elif name == 'SWAP':
                    left, right = wires
                    if _bit(basis, left) != _bit(basis, right):
                        output ^= (1 << left) | (1 << right)
                elif name == 'CCZ':
                    sign = (-1) ** (_bit(basis, wires[0]) * _bit(basis, wires[1]) * _bit(basis, wires[2]))
                elif name == 'Z':
                    sign = (-1) ** _bit(basis, wires[0])
                else:
                    raise AssertionError(name)
                result[output] = result.get(output, 0) + sign * amplitude
        state = {basis: amplitude for basis, amplitude in result.items() if amplitude}
        hadamards += name == 'H'
    assert hadamards % 2 == 0
    return state, 1 << (hadamards // 2)


def _ideal_probe():
    _, middle, final = _walls()
    d2 = [('CCZ', U1, X1, X2), ('CCZ', U2, X3, X4)]
    d3 = [('CCZ', TOP, X2, X4)]
    k = [('Z', ROOT), ('SWAP', U1, C1), ('SWAP', U2, C2)]
    return _inverse(final) + d3 + _inverse(middle) + d2 + k + d2 + middle + d3 + final


class JointAddressHelperSourceTests(unittest.TestCase):
    def test_old_shifted_cubic_has_wrong_complete_negative_rank(self):
        actual = sum(_bit(value, 4) * _bit(value, 5) for value in range(64))
        target = sum((value & 15) == 15 for value in range(64))
        self.assertEqual((actual, target), (16, 4))
        fiber = sum((1 ^ ((_bit(x, 0) * _bit(x, 1)))) * _bit(x, 2) * _bit(x, 3)
                    for x in range(16))
        self.assertEqual(fiber, 3)

    def test_changed_joint_word_on_every_ideal_column(self):
        word = _ideal_probe()
        for basis in range(1 << 10):
            actual, denominator = _integer_word(word, basis)
            output, sign = _trial_action(basis)
            self.assertEqual(actual, {output: sign * denominator})
            returned, reverse_sign = _trial_action(output)
            self.assertEqual((returned, sign * reverse_sign), (basis, 1))

    def test_two_exact_noncommuting_pauli_coefficients(self):
        # Hermitian P=Zr Zc1 Zu1 and R=Zr Xc1 Zc2 Zu2.
        p_x, p_z = 0, (1 << ROOT) | (1 << C1) | (1 << U1)
        r_x, r_z = 1 << C1, (1 << ROOT) | (1 << C2) | (1 << U2)
        coefficients = []
        for x_mask, z_mask in ((p_x, p_z), (r_x, r_z)):
            trace = 0
            for basis in range(1 << 10):
                output, sign = _trial_action(basis)
                if output ^ x_mask == basis:
                    trace += sign * (-1) ** ((output & z_mask).bit_count() & 1)
            coefficients.append(Fraction(trace, 1 << 10))
        self.assertEqual(coefficients, [Fraction(1, 16), Fraction(1, 16)])
        self.assertEqual(((p_x & r_z).bit_count() + (p_z & r_x).bit_count()) & 1, 1)

    def test_native_balanced_source_and_both_native_ledgers(self):
        first, source, _, correction, repaired = _words()
        occupied = (ROOT, U1, C1, U2, C2)
        for bits in range(32):
            basis = sum(((bits >> i) & 1) << wire for i, wire in enumerate(occupied))
            state = _native(_inverse(first), {basis: CO})
            state = {label: amplitude * (-1) ** _bit(label, ROOT) for label, amplitude in state.items()}
            actual = _native(first, state)
            expected = basis
            for left, right in ((U1, C1), (U2, C2)):
                if _bit(expected, left) != _bit(expected, right):
                    expected ^= (1 << left) | (1 << right)
            self.assertEqual(actual, {expected: CO * (-1) ** _bit(basis, ROOT)})
        for word, count, depth, gates in ((source, 36, 20, 103), (correction, 64, 32, 144),
                                          (repaired, 136, 72, 358)):
            for native in (word, _inverse(word)):
                counts = Counter(gate[0] for gate in native)
                self.assertLessEqual(counts['T'] + counts['TDG'], count)
                self.assertLessEqual(len(_bands(native)), depth)
                self.assertLessEqual(len(native), gates)
                self.assertTrue(all(0 <= wire < 12 for gate in native for wire in gate[1:]))

    def test_returned_complete_face_repair_and_literal_native_witnesses(self):
        # All 1,024 basis columns on the ten data/helper wires verify the correction.
        for basis in range(1 << 10):
            output, sign = _trial_action(basis)
            spectator = _bit(output, U1) * _bit(output, X2) + _bit(output, U2) * _bit(output, X4)
            sign *= (-1) ** spectator
            for controls, target in (((X3, X4), U1), ((X1, X2), U2),
                                     ((TOP, U1), C1), ((TOP, U2), C2),
                                     ((X3, X4), U1), ((X1, X2), U2),
                                     ((TOP, U1), C1), ((TOP, U2), C2)):
                output ^= (_bit(output, controls[0]) * _bit(output, controls[1])) << target
            self.assertEqual((output, sign), (basis, _target_sign(basis)))
        _, source, probe, _, repaired = _words()
        columns = (0, 15, 15 | (1 << TOP), 3 | (1 << C1) | (1 << U2),
                   12 | (1 << TOP) | (1 << C2) | (1 << U1), 1023)
        for basis in columns:
            for masks in (0, 1, 2, 3):
                value = basis | (masks << 10)
                expected, sign = _trial_action(value)
                self.assertEqual(_native(probe, {value: CO}), {expected: CO * sign})
                self.assertEqual(_native(repaired, {value: CO}), {value: CO * _target_sign(value)})
                output = _native(source, {value: CO})
                self.assertEqual(_native(_inverse(source), output), {value: CO})
        state = {15 | (1 << TOP): ROOT_HALF, 3 | (1 << C1) | (3 << 10): CI * ROOT_HALF}
        expected = {basis: amplitude * _target_sign(basis) for basis, amplitude in state.items()}
        self.assertEqual(_native(repaired, state), expected)
        self.assertEqual(_native(_inverse(repaired), expected), state)


if __name__ == '__main__':
    unittest.main()
