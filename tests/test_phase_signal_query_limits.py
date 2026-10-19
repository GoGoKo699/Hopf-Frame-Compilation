"""Exact native checks of coherent frame readout and a dirty signal echo.

The frame fixture keeps every key arbitrary, including all coherent Z
keys via exact phase coefficients. The signal fixture checks every column
of a five-qubit native word. Neither fixture supplies a faster source or
an unrestricted depth lower bound.
"""
from itertools import combinations
import unittest

import numpy as np

from tests.exact_matrices import CI, CO, CZ, adjoint, eye, multiply, native_columns
from tests.test_majorana_walsh_query import _inverse


def _cz(control, target):
    return [('H', target), ('CX', control, target), ('H', target)]


def _x(target):
    return [('H', target), ('S', target), ('S', target), ('H', target)]


def _parity_phase(wires, inverse=False):
    compute = [('CX', control, wires[-1]) for control in wires[:-1]]
    return compute + [('TDG' if inverse else 'T', wires[-1])] + list(reversed(compute))


def _ccz(wires):
    return [gate for size in (1, 2, 3) for subset in combinations(wires, size)
            for gate in _parity_phase(subset, inverse=size == 2)]


def _counts(word):
    return (sum(gate[0] in ('T', 'TDG') for gate in word),
            sum(gate[0] not in ('T', 'TDG') for gate in word))


def _phase_query_word():
    # x occupies wires 0,1 and Y occupies 2,3,4,5. The ANF is
    # Y0 + x0(Y0+Y1) + x1(Y0+Y2) + x0*x1(Y0+Y1+Y2+Y3).
    word = [('S', 2), ('S', 2)]
    for address, target in ((0, 3), (1, 4)):
        word += [('CX', 2, target)] + _cz(address, target) + [('CX', 2, target)]
    parity = [('CX', source, 5) for source in (2, 3, 4)]
    return word + parity + _ccz((0, 1, 5)) + list(reversed(parity))


def _monomial_table(matrix):
    roots = (CO, CI, -CO, -CI)
    images, phases = [], []
    for column in range(len(matrix)):
        supported = [(row, matrix[row][column]) for row in range(len(matrix))
                     if matrix[row][column] != CZ]
        assert len(supported) == 1
        row, coefficient = supported[0]
        images.append(row)
        phases.append(2 * roots.index(coefficient))
    return np.array(images, dtype=np.int64), np.array(phases, dtype=np.int64)


def _frame_sectors(events, d_tables, p_tables):
    """All x,Y,a,b,V sectors; each c,d bit remains a symbolic phase variable.

    Only the certified three-wire P/P-dagger blocks use c,d, linearly in
    the phase. A zero final coefficient mask therefore checks all their
    coherent superpositions, not just a selected classical assignment.
    """
    sectors = np.arange(1 << 16, dtype=np.int64)
    data, x_keys, bank = sectors & 63, (sectors >> 6) & 63, sectors >> 12
    phases, z_coefficients = np.zeros_like(sectors), np.zeros_like(sectors)
    for event in events:
        if event in ('D', 'Dinv'):
            images, diagonal = d_tables[event]
            phases += diagonal[data]
            data = images[data]
        elif event in ('P', 'Pinv'):
            images, phase = p_tables[event]
            for wire in range(6):
                local = ((data >> wire) & 1) | (((x_keys >> wire) & 1) << 1)
                low, high = phase[local], phase[local | 4]
                coefficient = (high - low) % 8
                assert np.all((coefficient == 0) | (coefficient == 4))
                phases += low
                z_coefficients ^= (coefficient // 4) << wire
                incoming = (data >> wire) & 1
                data ^= (incoming ^ (images[local] & 1)) << wire
        elif event == 'W':
            data ^= x_keys & 3  # x <- x+a
        elif event == 'L':
            x_keys ^= bank << 2  # b <- b+V
        elif event == 'Skey':
            phases += 2 * ((x_keys >> 2) & 1)
        elif event == 'Skey_inv':
            phases += 6 * ((x_keys >> 2) & 1)  # actual positive-S inverse
        else:
            raise ValueError(event)
    return sectors, data, x_keys, bank, phases % 8, z_coefficients


class PhaseSignalQueryLimitsTests(unittest.TestCase):
    def test_native_coherent_frame_extraction_returns_every_arbitrary_key(self):
        d_word = _phase_query_word()
        d = native_columns(6, d_word)
        d_inverse = native_columns(6, _inverse(d_word))
        self.assertEqual(d_inverse, adjoint(d))
        d_tables = {'D': _monomial_table(d), 'Dinv': _monomial_table(d_inverse)}
        for orientation in (d_word, _inverse(d_word)):
            self.assertEqual(_counts(orientation), (7, 28))
        for column in range(64):
            x, y = column & 3, column >> 2
            self.assertEqual(d_tables['D'][0][column], column)
            self.assertEqual(d_tables['D'][1][column], 4 * ((y >> x) & 1))

        # Native local frame: data wire 0, X key 1, Z key 2.
        local_p_word = _cz(2, 0) + [('CX', 1, 0)]
        local_p = native_columns(3, local_p_word)
        local_p_inverse = native_columns(3, _inverse(local_p_word))
        self.assertEqual(local_p_inverse, adjoint(local_p))
        p_tables = {'P': _monomial_table(local_p),
                    'Pinv': _monomial_table(local_p_inverse)}
        for kind, (images, phases) in p_tables.items():
            for column in range(8):
                data, a, z = column & 1, (column >> 1) & 1, column >> 2
                self.assertEqual(images[column], column ^ a)
                self.assertEqual(phases[column], 4 * z * (data if kind == 'P' else data ^ a))

        # Compile the actual 18-wire F source and its 22-wire wrapper.
        # All pieces were evaluated above on their full native columns;
        # their monomial composition below avoids a 2^22 square matrix.
        p_word = [gate for wire in range(6) for gate in _cz(12 + wire, wire)]
        p_word += [('CX', 6 + wire, wire) for wire in range(6)]
        f_word = _inverse(d_word) + p_word + d_word
        w_word = [('CX', 6 + wire, wire) for wire in range(2)]
        l_word = [('CX', 18 + wire, 8 + wire) for wire in range(4)]
        h_word = [('H', 18 + wire) for wire in range(4)]
        wrapper = (h_word + w_word + l_word + _inverse(f_word) + p_word
                   + l_word + _inverse(p_word) + f_word + w_word + h_word)
        for orientation in (f_word, _inverse(f_word)):
            self.assertEqual(_counts(orientation), (14, 80))
        for orientation in (wrapper, _inverse(wrapper)):
            self.assertEqual(_counts(orientation), (28, 228))
            self.assertTrue(all(gate[0] in ('H', 'S', 'CX', 'T', 'TDG')
                                for gate in orientation))
            self.assertLess(max(wire for gate in orientation for wire in gate[1:]), 22)

        f, f_inverse = ['Dinv', 'P', 'D'], ['Dinv', 'Pinv', 'D']
        for events, inverse in ((f, False), (f_inverse, True)):
            sectors, data, keys, bank, phases, z_coefficients = _frame_sectors(
                events, d_tables, p_tables)
            original, original_keys = sectors & 63, (sectors >> 6) & 63
            np.testing.assert_array_equal(data, original ^ original_keys)
            np.testing.assert_array_equal(keys, original_keys)
            np.testing.assert_array_equal(bank, sectors >> 12)
            np.testing.assert_array_equal(phases, (d_tables['D'][1][original]
                                                  + d_tables['D'][1][data]) % 8)
            np.testing.assert_array_equal(z_coefficients, data if inverse else original)

        events = ['W', 'L'] + f_inverse + ['P', 'L', 'Pinv'] + f + ['W']
        sectors, data, keys, bank, phases, z_coefficients = _frame_sectors(
            events, d_tables, p_tables)
        np.testing.assert_array_equal(data, sectors & 63)
        np.testing.assert_array_equal(keys, (sectors >> 6) & 63)
        np.testing.assert_array_equal(bank, sectors >> 12)
        expected = 4 * ((bank >> (sectors & 3)) & 1)
        np.testing.assert_array_equal(phases, expected)
        np.testing.assert_array_equal(z_coefficients, np.zeros_like(sectors))
        # H_V endpoints turn this complete phase identity into the exact
        # indicator on every output column. No key is fixed or initialized.

        # A phase depending on a measurement key is not a free gauge for
        # this coherent oracle. Keep S and its literal three-S inverse.
        gauged_f = f_word + [('S', 8)]
        self.assertEqual(_counts(gauged_f), (14, 81))
        self.assertEqual(_counts(_inverse(gauged_f)), (14, 83))
        gauged_events = (['W', 'L', 'Skey_inv'] + f_inverse
                         + ['P', 'L', 'Pinv'] + f + ['Skey', 'W'])
        _, gauged_data, gauged_keys, _, gauged_phases, gauge_z = _frame_sectors(
            gauged_events, d_tables, p_tables)
        b0, v0 = (keys >> 2) & 1, bank & 1
        np.testing.assert_array_equal(gauged_data, data)
        np.testing.assert_array_equal(gauged_keys, keys)
        np.testing.assert_array_equal(gauge_z, np.zeros_like(sectors))
        np.testing.assert_array_equal((gauged_phases - expected) % 8,
                                      (2 * b0 - 2 * (b0 ^ v0)) % 8)
        self.assertTrue(np.all((gauged_phases - expected)[v0 == 0] % 8 == 0))
        self.assertEqual((CI - CO).abs2(), 2)
        self.assertEqual((-CI - CO).abs2(), 2)

    def test_native_two_signal_echo_is_clean_correct_but_has_dirty_error_two(self):
        # p=3, ell=2: only pi and pi/2 controlled phases occur. CS is the
        # exact three-T parity word, not a supplied fine-angle rotation.
        signal = (3, 4)
        r_word = [('H', wire) for wire in signal]
        for address in range(3):
            r_word += _cz(address, signal[0])
            r_word += [('T', address), ('T', signal[1]), ('CX', address, signal[1]),
                       ('TDG', signal[1]), ('CX', address, signal[1])]
        r_word += [('H', wire) for wire in signal]
        local_s = _x(0) + _x(1) + _cz(0, 1) + _x(1) + _x(0)
        s_word = [(name, *(wire + 3 for wire in wires)) for name, *wires in local_s]
        echo_word = s_word + r_word + s_word + _inverse(r_word)  # R† S R S
        for orientation in (r_word, _inverse(r_word)):
            self.assertEqual(_counts(orientation), (9, 19))
        for orientation in (s_word, _inverse(s_word)):
            self.assertEqual(_counts(orientation), (0, 19))
        for orientation in (echo_word, _inverse(echo_word)):
            self.assertEqual(_counts(orientation), (18, 76))
            self.assertTrue(all(gate[0] in ('H', 'S', 'CX', 'T', 'TDG')
                                for gate in orientation))
        r = native_columns(5, r_word)
        r_inverse = native_columns(5, _inverse(r_word))
        self.assertEqual(r_inverse, adjoint(r))
        s = native_columns(2, local_s)
        self.assertEqual(s, [[(-CO if row == 0 else CO) if row == column else CZ
                              for column in range(4)] for row in range(4)])
        echo = native_columns(5, echo_word)
        self.assertEqual(native_columns(5, _inverse(echo_word)), adjoint(echo))
        for row in range(32):
            for column in range(32):
                if row & 7 != column & 7:
                    self.assertEqual(r[row][column], CZ)
                    self.assertEqual(echo[row][column], CZ)
        for x in range(8):
            weight = x.bit_count()
            block = [[echo[x | (row << 3)][x | (column << 3)]
                      for column in range(4)] for row in range(4)]
            r_block = [[r[x | (row << 3)][x | (column << 3)]
                        for column in range(4)] for row in range(4)]
            expected_r_word = [('H', 0), ('H', 1)]
            expected_r_word += [('S', 0)] * (2 * weight) + [('S', 1)] * weight
            expected_r_word += [('H', 0), ('H', 1)]
            self.assertEqual(r_block, native_columns(2, expected_r_word))
            self.assertEqual(block, multiply(adjoint(r_block), multiply(s, multiply(r_block, s))))
            if not x:
                self.assertEqual(block, eye(4))
                continue
            self.assertEqual(r_block[0][0], CZ)  # exact zero-state orthogonality
            self.assertEqual([row[0] for row in block], [-CO, CZ, CZ, CZ])
            self.assertEqual(block, adjoint(block))
            self.assertEqual(multiply(block, block), eye(4))
            self.assertEqual(sum((block[j][j] for j in range(4)), CZ), CZ)
            # A Hermitian involution of dimension four and trace zero has
            # two +1 and two -1 eigenvalues. This explicit +1 column also
            # certifies distance 2 from -I on an allowed dirty signal.
            witness = 2 if weight % 2 else 1
            self.assertEqual([row[witness] for row in block],
                             [CO if row == witness else CZ for row in range(4)])
            self.assertEqual(sum(((block[row][witness] + (CO if row == witness else CZ)).abs2()
                                  for row in range(4)), 0), 4)


if __name__ == '__main__':
    unittest.main()
