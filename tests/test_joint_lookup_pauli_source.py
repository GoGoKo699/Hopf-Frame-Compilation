"""Literal three-band extraction for a complete two-address bank.

The source uses relative phases and no initialized or borrowed work.
Its actual Pauli echo, rather than its basis permutation alone, is the
complete indicator contract. These finite identities are not an all-p
compiler or an optimal-depth result.
"""
from collections import Counter
import unittest

from tests.exact_matrices import CO, CZ, adjoint, native_columns
from tests.test_two_layer_diagonal_rigidity import (
    _primitive_inverse, _reduce_adjacent_s,
)
from tests.test_weak_pauli_orbit_compression import (
    _binary_rank, _input_mask, _target_mask,
)


def _inverse(word):
    return _reduce_adjacent_s(_primitive_inverse(word))


def _entry_encoding(a, b, d, t):
    """exp[-i pi Y_d (1-Z_t)/8], with both physical axes present."""
    assert len({a, b, d, t}) == 4
    return ([('S', d)] * 3
            + [('H', d), ('CX', d, t), ('T', d), ('TDG', t),
               ('CX', d, t), ('H', d), ('S', d)])


def _middle_band(a, b, d, t):
    """q=d-(b+d)+(a+b+d)-(a+d+t), where each sum is XOR."""
    basis = [('CX', d, b), ('CX', b, a), ('CX', a, t),
             ('CX', b, t), ('CX', d, t)]
    return (basis + [('T', d), ('TDG', b), ('T', a), ('TDG', t)]
            + list(reversed(basis)))


def _cubic_source():
    a, b, d, t = range(4)
    # The final word is F_0 dagger, with axes X_d Z_a and X_d Z_a Z_t.
    final = [('H', d), ('CX', a, d), ('CX', d, t),
             ('T', d), ('TDG', t), ('CX', d, t),
             ('CX', a, d), ('H', d)]
    return _entry_encoding(a, b, d, t) + _middle_band(a, b, d, t) + final


def _parity_basis(rows):
    """Emit an invertible CNOT basis, using only the specified wires."""
    width = len(rows)
    work = list(rows)
    reduction = []

    def add(source, target):
        work[target] ^= work[source]
        reduction.append(('CX', source, target))

    for column in range(width):
        pivot = next(row for row in range(column, width)
                     if (work[row] >> column) & 1)
        if pivot != column:
            add(pivot, column)
            add(column, pivot)
            add(pivot, column)
        for row in range(width):
            if row != column and (work[row] >> column) & 1:
                add(column, row)
    assert work == [1 << bit for bit in range(width)]
    return list(reversed(reduction))


def _exit_rows():
    # Wire order a,b,d,t,v,w. Each row includes d; H_d changes Z_d to X_d.
    a, b, d, t, v, w = [1 << bit for bit in range(6)]
    return [d ^ v, d ^ b, d ^ w, d ^ a ^ v, d ^ b ^ w, d ^ a ^ t]


def _exit_encoding():
    a, _, d, _, _, _ = range(6)
    basis = _parity_basis(_exit_rows())
    phase = [('TDG', 0), ('TDG', 1), ('TDG', 2),
             ('T', 3), ('T', 4), ('T', 5)]
    # Clifford coefficients 2-2 Z_a in s. Opposite scalar phases cancel.
    clifford = ([('H', d), ('S', d), ('CX', a, d)]
                + [('S', d)] * 3 + [('CX', a, d), ('H', d)])
    return (clifford + [('H', d)] + basis + phase
            + list(reversed(basis)) + [('H', d)])


def _normal_form_source():
    return (_entry_encoding(0, 1, 2, 3) + _middle_band(0, 1, 2, 3)
            + _inverse(_exit_encoding()))


def _bank_wall():
    _, _, d, t, v, w = range(6)
    return [('CX', d, v), ('CX', w, t), ('CX', d, w), ('CX', v, t)]


def _joint_lookup_source():
    # Physical bank order is Y0=d, Y1=v, Y2=w, Y3=t.
    hadamards = [('H', wire) for wire in (2, 4, 5, 3)]
    bank = _bank_wall()
    return (hadamards + bank + _normal_form_source()
            + list(reversed(bank)) + hadamards)


def _nested_coefficient_source():
    """Four bands transport Z_d to (-1)^(d+bt(a+cv)) on all six inputs."""
    a, b, d, t, c, v = range(6)
    # This Clifford simultaneously sends X_d Z_a to Z_d and X_a Z_d
    # to Z_a. The following five parity rows are independent.
    diagonalizer = [('H', d), ('CX', a, d), ('H', a)]
    basis = [('CX', d, t), ('CX', a, c), ('CX', a, v)]
    last = (diagonalizer + basis
            + [('T', d), ('TDG', t), ('TDG', a), ('T', c), ('T', v)]
            # The last rotation is the product of the a,c,v parity rows.
            + [('CX', a, v), ('CX', c, v), ('TDG', v),
               ('CX', c, v), ('CX', a, v)]
            + list(reversed(basis)) + list(reversed(diagonalizer)))
    return _entry_encoding(a, b, d, t) + _middle_band(a, b, d, t) + last


def _mobius_bank_word(p):
    """Fixed binary bank change Y -> its complete lookup ANF coefficients."""
    return [('CX', p + (row ^ (1 << bit)), p + row)
            for bit in range(p) for row in range(1 << p)
            if (row >> bit) & 1]


def _nested_coefficient_lift(p):
    assert p >= 3
    wall = _mobius_bank_word(p)
    # b,c,t are the first three addresses; d,a,v are coefficients 0,5,7.
    mapping = (p + 5, 0, p, 2, 1, p + 7)
    source = [(gate[0], *(mapping[wire] for wire in gate[1:]))
              for gate in _nested_coefficient_source()]
    return wall + source + list(reversed(wall))


def _two_coordinate_bank_basis(p):
    assert p >= 2
    size = 1 << p
    columns = [_target_mask(p, address) for address in range(4)]
    assert _binary_rank(columns) == 4
    for bit in range(size):
        if _binary_rank(columns + [1 << bit]) > len(columns):
            columns.append(1 << bit)
    assert len(columns) == size
    rows = [sum(((column >> row) & 1) << index
                for index, column in enumerate(columns))
            for row in range(size)]
    return columns, _parity_basis(rows)


def _two_coordinate_source(p):
    _, wall = _two_coordinate_bank_basis(p)
    wall = [('CX', p + source, p + target) for _, source, target in wall]
    # Input bank coordinates 0,1,2,3 represent the four independent masks.
    mapping = (0, 1, p, p + 3, p + 1, p + 2)
    source = [(gate[0], *(mapping[wire] for wire in gate[1:]))
              for gate in _joint_lookup_source()]
    return list(reversed(wall)) + source + wall


def _bands(word):
    result, active = [], []
    for gate in word:
        if gate[0] in ('T', 'TDG'):
            active.append(gate)
        elif active:
            result.append(active)
            active = []
    if active:
        result.append(active)
    return result


def _rows_after(width, word):
    rows = [1 << bit for bit in range(width)]
    for name, source, target in word:
        assert name == 'CX'
        rows[target] ^= rows[source]
    return rows


class JointLookupPauliSourceTests(unittest.TestCase):
    def test_cubic_source_all_columns_and_literal_actual_inverse(self):
        word = _cubic_source()
        actual = native_columns(4, word)
        for column in range(16):
            row = column ^ (4 if column & 11 == 11 else 0)
            support = [index for index in range(16)
                       if actual[index][column] != CZ]
            self.assertEqual(support, [row])
            phase = actual[row][column]
            self.assertEqual(phase * phase.conj(), CO)
        # Relative phases are retained, not mistaken for a literal C3X.
        self.assertNotEqual(actual[3][3], CO)
        self.assertEqual(native_columns(4, _inverse(word)), adjoint(actual))

    def test_joint_lookup_complete_native_echo_and_actual_inverse(self):
        normal = native_columns(6, _normal_form_source())
        for column in range(64):
            a, b, _, t, v, w = [(column >> bit) & 1 for bit in range(6)]
            row = column ^ ((a * v ^ b * w ^ a * b * t) << 2)
            support = [index for index in range(64)
                       if normal[index][column] != CZ]
            self.assertEqual(support, [row])
            phase = normal[row][column]
            self.assertEqual(phase * phase.conj(), CO)

        source = _joint_lookup_source()
        inverse = _inverse(source)
        actual = native_columns(6, source)
        self.assertEqual(native_columns(6, inverse), adjoint(actual))
        # Chronological E dagger, X_Y0, E realizes E X_Y0 E dagger.
        middle = [('H', 2), ('S', 2), ('S', 2), ('H', 2)]
        echo = native_columns(6, inverse + middle + source)
        bank_wires = (2, 4, 5, 3)
        for column in range(64):
            row = column ^ (1 << bank_wires[column & 3])
            for output in range(64):
                self.assertEqual(echo[output][column], CO if output == row else CZ)
        # Full matrix equality includes coherent bank/address/reference inputs.

    def test_physical_independence_complete_phase_identity_and_ledgers(self):
        basis = _parity_basis(_exit_rows())
        self.assertEqual(len(basis), 20)
        self.assertEqual(_rows_after(6, basis), _exit_rows())
        self.assertEqual(_rows_after(6, basis + list(reversed(basis))),
                         [1 << bit for bit in range(6)])
        for column in range(64):
            a, b, _, t, v, w = [(column >> bit) & 1 for bit in range(6)]
            za, zb, zt, zv, zw = [1 - 2 * bit for bit in (a, b, t, v, w)]
            polynomial = 4 * a * v + 4 * b * w + 4 * a * t - 2 * t
            axes = 2 - 2 * za - zv - zb - zw + za * zv + zb * zw + za * zt
            self.assertEqual((polynomial - axes) % 8, 0)

        for word, targets, cliffords in (
                (_cubic_source(), [2, 4, 2], 24),
                (_normal_form_source(), [2, 4, 6], 68),
                (_joint_lookup_source(), [2, 4, 6], 84)):
            for native in (word, _inverse(word)):
                bands = _bands(native)
                self.assertEqual(sorted(map(len, bands)), sorted(targets))
                self.assertEqual(len(bands), 3)
                for band in bands:
                    self.assertEqual(len({gate[1] for gate in band}), len(band))
                counts = Counter(gate[0] for gate in native)
                self.assertEqual(counts['T'] + counts['TDG'], sum(targets))
                self.assertEqual(counts['H'] + counts['S'] + counts['CX'], cliffords)
        # The only physical wires are the six specified inputs.
        self.assertEqual({wire for gate in _joint_lookup_source() for wire in gate[1:]},
                         set(range(6)))

    def test_four_band_nested_coefficient_complete_native_action(self):
        word = _nested_coefficient_source()
        inverse = _inverse(word)
        actual = native_columns(6, word)
        for column in range(64):
            a, b, _, t, c, v = [(column >> bit) & 1 for bit in range(6)]
            row = column ^ (c * v) ^ ((a * b * t) << 2)
            support = [index for index in range(64)
                       if actual[index][column] != CZ]
            self.assertEqual(support, [row])
            phase = actual[row][column]
            self.assertEqual(phase * phase.conj(), CO)
        self.assertEqual(native_columns(6, inverse), adjoint(actual))
        echo = native_columns(6, inverse + [('S', 2), ('S', 2)] + word)
        for column in range(64):
            a, b, d, t, c, v = [(column >> bit) & 1 for bit in range(6)]
            sign = -CO if d ^ (a * b * t) ^ (b * c * t * v) else CO
            for row in range(64):
                self.assertEqual(echo[row][column], sign if row == column else CZ)
        for native in (word, inverse):
            bands = _bands(native)
            self.assertEqual(sorted(map(len, bands)), [1, 2, 4, 5])
            for band in bands:
                self.assertEqual(len({gate[1] for gate in band}), len(band))
            counts = Counter(gate[0] for gate in native)
            self.assertEqual(counts['T'] + counts['TDG'], 12)
            self.assertEqual(counts['H'] + counts['S'] + counts['CX'], 34)
        self.assertEqual({wire for gate in word for wire in gate[1:]}, set(range(6)))

    def test_nested_coefficient_lift_and_missing_independent_coefficients(self):
        for p in (3, 4, 6):
            size = 1 << p
            wall = _mobius_bank_word(p)
            self.assertEqual(len(wall), p * size // 2)
            bank_word = [('CX', a - p, b - p) for _, a, b in wall]
            rows = _rows_after(size, bank_word)
            self.assertEqual(rows[0], 1)
            self.assertEqual(rows[5], sum(1 << row for row in (0, 1, 4, 5)))
            self.assertEqual(rows[7], 255)
            source = _nested_coefficient_lift(p)
            counts = Counter(gate[0] for gate in source)
            self.assertEqual(counts['T'] + counts['TDG'], 12)
            self.assertEqual(counts['H'] + counts['S'] + counts['CX'], 34 + p * size)
            self.assertEqual([len(band) for band in _bands(source)], [2, 4, 5, 1])

        p, size = 3, 8
        source = _nested_coefficient_lift(p)
        inputs = [0, 5 | (1 << p), 7 | (1 << (p + 7)), 1 | (1 << (p + 1))]
        echo = native_columns(p + size,
                              _inverse(source) + [('S', p), ('S', p)] + source,
                              inputs)
        for index, column in enumerate(inputs):
            b, c, t = [(column >> bit) & 1 for bit in range(3)]
            bank = column >> p
            a = (bank & 51).bit_count() & 1
            v = bank.bit_count() & 1
            sign = -CO if (bank & 1) ^ (b * t * a) ^ (b * c * t * v) else CO
            for row in range(1 << (p + size)):
                self.assertEqual(echo[row][index], sign if row == column else CZ)
        # The omitted linear-b coefficient is independent of all paid terms.
        self.assertEqual(echo[inputs[-1]][-1], CO)
        self.assertEqual((inputs[-1] >> (p + 1)) & 1, 1)

    def test_parameterized_two_coordinate_lift_and_active_spectator(self):
        for p in (2, 3, 5):
            size = 1 << p
            columns, wall = _two_coordinate_bank_basis(p)
            rows = _rows_after(size, wall)
            emitted_columns = [sum(((rows[row] >> bit) & 1) << row
                                   for row in range(size))
                               for bit in range(size)]
            self.assertEqual(emitted_columns, columns)
            self.assertEqual(columns[:4], [_target_mask(p, address)
                                           for address in range(4)])
            self.assertLessEqual(len(wall), size * size + 3 * size)
            source = _two_coordinate_source(p)
            self.assertEqual([len(band) for band in _bands(source)], [2, 4, 6])
            counts = Counter(gate[0] for gate in source)
            self.assertEqual(counts['T'] + counts['TDG'], 12)
            self.assertLessEqual(counts['H'] + counts['S'] + counts['CX'],
                                 84 + 2 * (size * size + 3 * size))

        # These exact p=3 columns challenge the emitted lift with both values
        # of the unprocessed address bit and several arbitrary bank values.
        p, size = 3, 8
        source = _two_coordinate_source(p)
        pauli = [gate for bit in range(size) if (_input_mask(p) >> bit) & 1
                 for gate in [('H', p + bit), ('S', p + bit),
                              ('S', p + bit), ('H', p + bit)]]
        inputs = [0, 1 | (17 << p), 2 | (85 << p), 3 | (255 << p),
                  4 | (68 << p), 7 | (170 << p)]
        actual = native_columns(p + size, _inverse(source) + pauli + source,
                                inputs)
        for index, column in enumerate(inputs):
            output = column ^ (_target_mask(p, column & 3) << p)
            for row in range(1 << (p + size)):
                self.assertEqual(actual[row][index], CO if row == output else CZ)
        self.assertNotEqual(_target_mask(p, 7), _target_mask(p, 3))


if __name__ == '__main__':
    unittest.main()
