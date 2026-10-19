"""Exact checks for the external-center middle correction.

The nine-wire source returns its one clean wire on every column, but its
occupied diagonal has non-Boolean phases.  The separate six-axis example
tests the intrinsic obstruction to putting the corrected middle operation
back into one native band.  A small complex-alignment check distinguishes
an identity-mode equality of the whole residual from its separate factors.
None of these checks asserts a growing-degree query.
"""
from collections import Counter
from itertools import combinations
import unittest

from tests.exact_matrices import (
    CO, CI, Gate, adjoint, apply_word, multiply, native_columns,
)
from tests.test_bounded_relation_native_matching import _rank
from tests.test_two_layer_galois_return import _conjugate, _inverse, _symplectic, _tau


def _native(word):
    return [Gate(name, tuple(qubits)) for name, *qubits in word]


def _linear_word(rows):
    """CNOT word for the specified invertible binary coordinate map."""
    reduced = list(rows)
    elimination = []
    for bit in range(len(rows)):
        pivot = next(index for index in range(bit, len(rows))
                     if reduced[index] >> bit & 1)
        if pivot != bit:
            for control, target in ((pivot, bit), (bit, pivot), (pivot, bit)):
                reduced[target] ^= reduced[control]
                elimination.append(('CX', control, target))
        for index in range(len(rows)):
            if index != bit and reduced[index] >> bit & 1:
                reduced[index] ^= reduced[bit]
                elimination.append(('CX', bit, index))
    assert reduced == [1 << bit for bit in range(len(rows))]
    return list(reversed(elimination))


def _complete_rows(rows, width):
    result = list(rows)
    for bit in range(width):
        if _rank(result + [1 << bit]) > len(result):
            result.append(1 << bit)
    assert len(result) == width
    return result


def _source():
    # Physical order a,c,z1,z2,z3,z4,r,q,d; a alone is initialized.
    width = 9
    a, c, r, q, d = 0, 1, 6, 7, 8
    z = tuple(range(2, 6))
    am, cm, rm, qm, dm = [1 << bit for bit in (a, c, r, q, d)]
    zm = [1 << bit for bit in z]
    w, h = sum(zm), am ^ cm
    p, y = h ^ rm, h ^ qm
    active_a = [am] + [am ^ axis for axis in zm]
    active_c = [cm] + [cm ^ axis for axis in zm]
    negative = [am ^ w, cm ^ w]
    positive = [y ^ axis for axis in zm] + [y, y ^ p ^ w]

    # C0=S_a S_w CZ_(c,w), with one shared parity computation.
    parity = [('CX', bit, z[-1]) for bit in z[:-1]]
    c0 = ([('S', a)] + parity + [('S', z[-1]), ('H', z[-1]),
          ('CX', c, z[-1]), ('H', z[-1])] + _inverse(parity))
    a_rows = active_a + [cm, am ^ cm ^ rm, dm]
    c_rows = active_c + [am, am ^ cm ^ rm, dm ^ cm ^ rm]
    changes = (
        [('CX', c, r), ('H', c)] + _linear_word(_complete_rows(a_rows, width)),
        _linear_word(_complete_rows(negative + positive, width)),
        _inverse(c0) + [('CX', a, r), ('H', a)]
        + _linear_word(_complete_rows(c_rows, width)),
    )
    phases = [[('TDG' if band == 1 and bit == 6 else 'T', bit)
               for bit in range(8)] for band in range(3)]
    bands = [change + phase + _inverse(change)
             for change, phase in zip(changes, phases)]
    ka = [('H', a), ('CX', a, c)]
    kc = ka + c0
    word = ka + bands[0] + bands[1] + _inverse(bands[2]) + _inverse(kc)

    # Binary labels retain the actual C0-conjugated C-silent axes.
    zp = lambda value: value << width
    qa = cm ^ rm
    qc = am ^ rm ^ zp(am)  # Y_a X_r; its literal sign is positive.
    labels = ([zp(axis) for axis in active_a] + [qa, qa ^ zp(p), zp(dm)],
              [zp(axis) for axis in negative + positive],
              [zp(axis) for axis in active_c] + [qc, qc ^ zp(p), zp(dm ^ p)])
    return word, bands, labels, changes, phases


def _occupied_phase(column):
    t, r, q, d = [(column >> bit) & 1 for bit in (1, 6, 7, 8)]
    z = [(column >> bit) & 1 for bit in range(2, 6)]
    w = sum(z) % 2
    beta = sum(z[i] * z[j] for i, j in combinations(range(4), 2)) % 2
    control, p, y = beta ^ 1 ^ w, t ^ r, t ^ q
    return (sum(z) + w + 4 * q * control
            + 2 * p * (d - w - y + 2 * w * y)) % 8


class ExternalCenterMiddleCorrectionTests(unittest.TestCase):
    def test_complex_alignment_and_whole_identity_mode_hadamard(self):
        # C=S0 CZ01 is an actual complex diagonal Clifford. The residual
        # R=T1 T_(Z1 Z2) is identity on mode0 but remains non-Clifford.
        c_word = [('S', 0), ('H', 1), ('CX', 0, 1), ('H', 1)]
        r_word = [('T', 1), ('CX', 1, 2), ('T', 2), ('CX', 1, 2)]
        d_word = c_word + r_word
        c, residual, diagonal = [native_columns(3, word)
                                 for word in (c_word, r_word, d_word)]
        self.assertEqual(diagonal, multiply(c, residual))
        self.assertEqual(multiply(adjoint(c), diagonal), residual)
        # The entire commuting Pauli family is intertwined, with literal
        # signs. The i entry in the X0 image rules out a real aligner:
        # conjugating the real X0 matrix by any real matrix stays real.
        x_image = _tau(1 | (1 << 3) | (1 << 4), 3)  # Y0 Z1.
        self.assertEqual(_conjugate(c, _tau(1, 3)), x_image)
        self.assertEqual(x_image[1][0], CI)
        for mask in range(8):
            label = ((mask & 1) | ((mask & 2) << 3) | ((mask & 4) << 3))
            pauli = _tau(label, 3)
            self.assertEqual(_conjugate(diagonal, pauli), _conjugate(c, pauli))
            self.assertEqual(multiply(residual, pauli), multiply(pauli, residual))

        h = native_columns(3, [('H', 0)])
        self.assertEqual(_conjugate(h, residual), residual)
        # Peeling R=W F does not make either factor identity on mode0.
        # Only their complete product may use that identity-mode equality.
        w_word, f_word = [('T', 0)] + r_word, [('TDG', 0)]
        w, f = [native_columns(3, word) for word in (w_word, f_word)]
        self.assertEqual(multiply(w, f), residual)
        changed_w, changed_f = _conjugate(h, w), _conjugate(h, f)
        self.assertNotEqual(changed_w, w)
        self.assertNotEqual(changed_f, f)
        self.assertNotEqual(multiply(changed_w, f), residual)
        self.assertNotEqual(multiply(w, changed_f), residual)
        self.assertEqual(multiply(changed_w, changed_f), residual)
        for word, matrix in ((c_word, c), (r_word, residual), (d_word, diagonal),
                             (w_word, w), (f_word, f)):
            self.assertEqual(native_columns(3, _inverse(word)), adjoint(matrix))

    def test_actual_clean_source_every_column_and_reversed_native_word(self):
        word, bands, labels, _, phases = _source()
        for axes, phase in zip(labels, phases):
            self.assertEqual(_rank(axes), 8)
            self.assertEqual(len({gate[1] for gate in phase}), 8)
            for left in axes:
                for right in axes:
                    self.assertEqual(_symplectic(left, right, 9), 0)
        self.assertEqual(sum(gate[0] in ('T', 'TDG') for gate in word), 24)
        omega = native_columns(1, [('T', 0)])[1][1]
        roots = [CO]
        for _ in range(7):
            roots.append(roots[-1] * omega)
        native, inverse = _native(word), _native(_inverse(word))
        for column in range(1 << 9):
            clean = column & 1
            p = ((column >> 1) ^ (column >> 6)) & 1
            flip = clean * (1 ^ p)
            row = column ^ ((1 << 1) | (1 << 6) if flip else 0)
            # V=diag(G,G L), L=-i X_t X_r on p=0 and I on p=1.
            amplitude = roots[_occupied_phase(row)] * (-CI if flip else CO)
            self.assertEqual(apply_word(native, {column: CO}), {row: amplitude})
            self.assertEqual(apply_word(inverse, {row: CO}),
                             {column: amplitude.conj()})
            self.assertEqual(row & 1, clean)
        self.assertEqual(_occupied_phase(0), 0)
        self.assertEqual(_occupied_phase((1 << 2) | (1 << 3)), 2)
        self.assertEqual(roots[2], CI)  # No global scalar makes both columns real.

        # Three actual Clifford-conjugated simultaneous bands, including
        # the displayed orientations; the inverse is independently emitted.
        counts = Counter(gate[0] for gate in word)
        inverse_counts = Counter(gate[0] for gate in _inverse(word))
        self.assertEqual(counts,
                         {'H': 12, 'CX': 205, 'T': 15, 'TDG': 9, 'S': 2, 'SDG': 4})
        self.assertEqual(inverse_counts,
                         {'H': 12, 'CX': 205, 'T': 9, 'TDG': 15, 'S': 4, 'SDG': 2})
        # Raw SDG->SSS expansion, without post-expansion cancellations.
        for ledger, cliffords in ((counts, 231), (inverse_counts, 227)):
            self.assertEqual(ledger['H'] + ledger['CX'] + ledger['S']
                             + 3 * ledger['SDG'], cliffords)
        self.assertTrue(all(sum(g[0] in ('T', 'TDG') for g in band) == 8
                            for band in bands))

    def test_characteristic_defect_and_quadratic_middle_correction(self):
        width, all_bits = 6, 63
        # The actual characteristic transvection sends every Z_i to the
        # product of the other five. It is not a native-axis permutation.
        rows = [all_bits ^ (1 << bit) for bit in range(width)]
        change = _native(_linear_word(rows))
        inverse = _native(_inverse(_linear_word(rows)))
        # One fixed original-native product v=e0 solves Nv=p.
        fixed_product = 1
        fixed_image = sum(((row & fixed_product).bit_count() & 1) << bit
                          for bit, row in enumerate(rows))
        self.assertEqual(fixed_image ^ fixed_product, all_bits)
        middle = _native([('T', bit) for bit in range(width)])
        omega = native_columns(1, [('T', 0)])[1][1]
        roots = [CO]
        for _ in range(7):
            roots.append(roots[-1] * omega)

        def control(state):
            r = [((state >> bit) ^ (state >> 5)) & 1 for bit in range(5)]
            return 1 ^ (sum(r[i] * r[j] for i, j in combinations(range(5), 2)) % 2)

        def correction(state):
            return control(state) * (state & 1)

        correction_derivative = []
        for state in range(64):
            parity = state.bit_count() & 1
            image = state ^ (all_bits if parity else 0)
            self.assertEqual(apply_word(change, {state: CO}), {image: CO})
            self.assertEqual(apply_word(inverse, {image: CO}), {state: CO})
            self.assertEqual(apply_word(middle, {state: CO}),
                             {state: roots[state.bit_count() % 8]})
            self.assertEqual((state.bit_count() - image.bit_count()) % 8,
                             4 * parity * control(state))
            self.assertEqual(control(state ^ all_bits), control(state))
            corrected = (state.bit_count() + 4 * correction(state)) % 8
            self.assertEqual(2 * corrected % 8, 2 * state.bit_count() % 8)
            derivative = correction(state) ^ correction(state ^ all_bits)
            self.assertEqual(derivative, control(state))
            correction_derivative.append(derivative)
        # The derivative is genuinely quadratic, not an affine sign
        # available from the square-preserving native one-band alternatives.
        self.assertEqual(correction_derivative[0] ^ correction_derivative[1]
                         ^ correction_derivative[2] ^ correction_derivative[3], 1)
        polar_columns = []
        for i in range(width):
            polar_columns.append(sum(
                (control(0) ^ control(1 << i) ^ control(1 << j)
                 ^ control((1 << i) ^ (1 << j))) << j for j in range(width)))
        self.assertEqual(_rank(polar_columns), 4)
        self.assertEqual(_rank(polar_columns + [all_bits]), 5)
        support = {0}
        for axis in polar_columns:
            support |= {value ^ axis for value in support}
        self.assertEqual(len(support), 16)
        # In the quotient by the quadratic support, the six original
        # native axes occupy exactly two nonzero cosets, both oddly often.
        quotient = Counter(min((1 << bit) ^ value for value in support)
                           for bit in range(width))
        self.assertEqual(sorted(quotient.values()), [1, 5])
        self.assertNotIn(0, quotient)


if __name__ == '__main__':
    unittest.main()
