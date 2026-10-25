"""Exact complete-action checks for affine/Fourier relative-probe dressings.

The all-size transition-support and exact monomial arguments are in
COUPLED_FINAL_ENCODER, Section 11. These fixtures include actual native
echoes, both affine/Fourier orders, literal phases, and arbitrary inputs.
"""
from collections import Counter
from fractions import Fraction
import unittest

from tests.exact_matrices import CO, CZ
from tests.test_common_probe_source import _d_phase, _diagonal, _native
from tests.test_joint_lookup_pauli_source import _bands, _inverse
from tests.test_relative_probe_repair import (
    _ccz_parity_word, _probe_state, _q_sign, _repair_parts,
)


def _affine_word(m):
    """Nontrivial invertible linear mixing and a nonzero affine offset."""
    return ([('CX', j, j + 1) for j in range(m - 1)]
            + [('CX', m - 1, 0)]
            + [('H', 1), ('S', 1), ('S', 1), ('H', 1)])


def _dressing(p, order):
    m = 1 << p
    l = _affine_word(m)
    h = [('H', j) for j in range(m)]
    # Native lists are chronological; the labels denote operator order.
    if order == 'HL':
        return l + h
    if order == 'LH':
        return h + l
    raise ValueError(order)


def _x_rotation(p):
    m = 1 << p
    h = [('H', j) for j in range(m)]
    return h + [('T', j) for j in range(m)] + h


def _parity(state):
    return _diagonal(state, lambda y: -CO if y.bit_count() & 1 else CO)


def _actual_transformed(state, p, address, dressing):
    """D† C (T Z_parity T†) C† D, retaining native scalar cancellation."""
    t = _x_rotation(p)
    state = _diagonal(state, lambda y: _d_phase(y, p, address))
    state = _native(_inverse(dressing) + _inverse(t), state)
    state = _parity(state)
    state = _native(t + dressing, state)
    return _diagonal(state, lambda y: _d_phase(y, p, address).conj())


def _target_transformed(state, p, address, dressing, include_t=True):
    """T C Q_s C† T†; omitting T tests the proved X-commutation."""
    t = _x_rotation(p) if include_t else []
    state = _native(_inverse(t) + _inverse(dressing), state)
    state = _diagonal(state, lambda y: _q_sign(y, p, address) * CO)
    return _native(dressing + t, state)


def _inner(left, right):
    return sum((value.conj() * right.get(key, CZ)
                for key, value in left.items()), CZ)


def _native_source(p, dressing):
    wall, a, b, _, _ = _repair_parts(p)
    shifted = [(name, *(p + wire for wire in wires))
               for name, *wires in dressing]
    c = wall + shifted + _inverse(wall)
    word = b + a + c + _inverse(b) + _inverse(a) + _inverse(c)
    return wall, word


def _walsh(values):
    out = list(values)
    stride = 1
    while stride < len(out):
        for start in range(0, len(out), 2 * stride):
            for j in range(stride):
                left, right = out[start + j], out[start + j + stride]
                out[start + j] = left + right
                out[start + j + stride] = left - right
        stride *= 2
    return [Fraction(value, len(out)) for value in out]


class AffineFourierRepairTests(unittest.TestCase):
    def test_changed_correlated_x_phases_and_diagonal_middle_keep_the_gap(self):
        p, m, address, bank = 2, 4, 1, 3
        h = [('H', j) for j in range(m)]
        x_phases = ([('T', 0), ('S', 2), ('CX', 0, 1), ('T', 1), ('CX', 0, 1)]
                    + _ccz_parity_word(1, 2, 3))
        t = h + x_phases + h
        middle = ([('T', 0), ('CX', 1, 3), ('TDG', 3), ('CX', 1, 3)]
                  + _ccz_parity_word(0, 2, 3))
        for order in ('HL', 'LH'):
            with self.subTest(order=order):
                dressing = _dressing(p, order)
                actual = _native(middle + _inverse(dressing) + _inverse(t), {bank: CO})
                actual = _parity(actual)
                actual = _native(t + dressing + _inverse(middle), actual)
                target = _native(_inverse(t) + _inverse(dressing), {bank: CO})
                target = _diagonal(target, lambda y: _q_sign(y, p, address) * CO)
                target = _native(dressing + t, target)
                self.assertEqual(len(actual), 1)
                self.assertNotIn(next(iter(actual)), target)
                self.assertEqual(_inner(actual, actual), CO)
                self.assertEqual(_inner(target, target), CO)
                self.assertEqual(_inner(actual, target), CZ)

    def test_both_orders_have_disjoint_complete_transition_support(self):
        for p, columns in ((2, (0, 1, 7)), (3, (0, 37))):
            for order in ('HL', 'LH'):
                dressing = _dressing(p, order)
                displacement = None
                for address in (0, (1 << p) - 1):
                    for bank in columns:
                        with self.subTest(p=p, order=order, address=address, bank=bank):
                            actual = _actual_transformed({bank: CO}, p, address, dressing)
                            target = _target_transformed({bank: CO}, p, address, dressing)
                            self.assertEqual(target, _target_transformed(
                                {bank: CO}, p, address, dressing, include_t=False))
                            self.assertEqual(len(actual), 1)
                            output = next(iter(actual))
                            if displacement is None:
                                displacement = bank ^ output
                            self.assertEqual(bank ^ output, displacement)
                            self.assertNotIn(output, target)
                            self.assertEqual(_inner(actual, actual), CO)
                            self.assertEqual(_inner(target, target), CO)
                            self.assertEqual(_inner(actual, target), CZ)

    def test_p3_actual_native_echo_on_the_coherent_witness(self):
        p, m, address = 3, 8, 5
        d = p + m
        for order in ('HL', 'LH'):
            with self.subTest(order=order):
                dressing = _dressing(p, order)
                # psi=C†T†|0>; this is a full input witness, not free source work.
                psi = _native(_inverse(_x_rotation(p)) + _inverse(dressing), {0: CO})
                self.assertEqual(_inner(psi, psi), CO)
                wall, source = _native_source(p, dressing)
                initial = _probe_state(psi, p, address, -1)
                physical = _native(wall, initial)
                self.assertEqual(_native(_inverse(source), _native(source, physical)), physical)
                echo = wall + _inverse(source) + [('S', d)] * 2 + source + _inverse(wall)
                observed = _native(echo, initial)
                transformed = _actual_transformed({0: CO}, p, address, dressing)
                relative = _native(_inverse(_x_rotation(p)) + _inverse(dressing), transformed)
                self.assertEqual(observed, _probe_state(relative, p, address, 1))
                target = _probe_state(_diagonal(
                    psi, lambda y: _q_sign(y, p, address) * CO), p, address, 1)
                self.assertEqual(_inner(observed, observed), CO)
                self.assertEqual(_inner(target, target), CO)
                self.assertEqual(_inner(observed, target), CZ)

    def test_native_bands_and_paid_affine_clifford_ledger(self):
        for p in (2, 3, 4, 6):
            m = 1 << p
            g_l = len(_affine_word(m))
            for order in ('HL', 'LH'):
                _, word = _native_source(p, _dressing(p, order))
                for native in (word, _inverse(word)):
                    counts = Counter(g[0] for g in native)
                    self.assertEqual(counts['T'] + counts['TDG'], 4 * m)
                    self.assertEqual(counts['H'] + counts['S'] + counts['CX'],
                                     6 * p * m + 26 * m - 12 + 2 * g_l)
                    bands = _bands(native)
                    self.assertEqual(list(map(len, bands)), [m] * 4)
                    self.assertTrue(all(len({g[1] for g in band}) == m for band in bands))

    def test_exact_monomial_screen_noise_identity_and_flat_actual_columns(self):
        p, m, address = 2, 4, 1
        # Nonlinear diagonal phases are allowed in the monomial screen;
        # these S/T phases cancel in C Q C† but remain in actual amplitudes.
        dressing = (_affine_word(m) + [('S', 0), ('T', 2)]
                    + _ccz_parity_word(0, 2, 3))
        inverse_l = _inverse(_affine_word(m))
        truth = []
        for bank in range(1 << m):
            back = _native(inverse_l, {bank: CO})
            self.assertEqual(len(back), 1)
            truth.append(_q_sign(next(iter(back)), p, address))
        coeff = _walsh(truth)
        noise_norm = sum(c * c * Fraction(1, 1 << mask.bit_count())
                         for mask, c in enumerate(coeff))
        diagonal_norm = CZ
        for bank in range(1 << m):
            actual = _actual_transformed({bank: CO}, p, address, dressing)
            self.assertEqual(len(actual), 1 << m)
            self.assertTrue(all(value * value.conj() == Fraction(1, 1 << m) * CO
                                for value in actual.values()))
            target = _target_transformed({bank: CO}, p, address, dressing)
            entry = target.get(bank, CZ)
            diagonal_norm += entry * entry.conj()
        self.assertEqual(diagonal_norm, (1 << m) * noise_norm * CO)
        self.assertGreater(noise_norm, Fraction(1, 1 << m))
        self.assertEqual(sum(c * c for c in coeff), 1)
        self.assertNotEqual(abs(coeff[-1]), 1)

    def test_nonlinear_parity_transport_fails_the_walsh_rectangle_identity(self):
        p, m, address = 2, 4, 1
        # This explicit mathematical permutation passes the first necessary
        # condition C Q C†=Z_parity. No cheap native implementation is assumed.
        sigma = {}
        for sign in (1, -1):
            source = [y for y in range(1 << m) if _q_sign(y, p, address) == sign]
            destination = [y for y in range(1 << m)
                           if (-1 if y.bit_count() & 1 else 1) == sign]
            self.assertEqual(len(source), len(destination))
            sigma.update(zip(source, destination))
        inverse = {target: source for source, target in sigma.items()}
        self.assertEqual(len(inverse), 1 << m)
        for y in range(1 << m):
            self.assertEqual(_q_sign(inverse[y], p, address),
                             -1 if y.bit_count() & 1 else 1)
        mismatches = []
        for x in range(1 << m):
            gx = inverse[x] ^ inverse[0]
            for y in range(1 << m):
                gy = inverse[y] ^ inverse[0]
                # These are the exact signs of the phased-Walsh rectangle;
                # arbitrary diagonal row/column factors cancel identically.
                permuted = -1 if (gx & gy).bit_count() & 1 else 1
                original = -1 if (x & y).bit_count() & 1 else 1
                if permuted != original:
                    mismatches.append((x, y))
        self.assertTrue(mismatches)


if __name__ == '__main__':
    unittest.main()
