"""Exact weak single-coordinate sources and their complete-orbit limitation.

The native fixture retains every address, bank and arbitrary-helper column,
literal phase and actual inverse. The all-size failure and prefix-rank
requirement are analytic statements in WEAK_PAULI_ORBIT_COMPRESSION.md;
these binary finite checks are not unrestricted native-depth bounds.
"""
from itertools import product
import unittest

from tests.exact_matrices import CO, adjoint, native_columns, zeros
from tests.test_bent_query_baseline import _compact_toffoli
from tests.test_phase_signal_query_limits import _counts
from tests.test_t_depth import _inverse, _word


def _signs(p):
    return tuple((1,) + tail for tail in product((1, -1), repeat=p-1))


def _input_mask(p):
    return sum(1 << (2*j) for j in range(1 << (p-1)))


def _target_mask(p, address):
    powers = (1, 2, 3)  # A^r(1,0), encoded as first + 2*second.
    return sum(powers[sum(sign*((address >> i) & 1)
                          for i, sign in enumerate(row)) % 3] << (2*j)
               for j, row in enumerate(_signs(p)))


def _displacement(p, coordinate):
    return _input_mask(p) ^ _target_mask(p, 1 << coordinate)


def _representative(p, coordinate, vector):
    return vector ^ (_displacement(p, coordinate)
                     if (vector ^ (vector >> 1)) & 1 else 0)


def _clifford_wall(p, coordinate):
    first, second = p, p+1
    gates = []
    for j, row in enumerate(_signs(p)[1:], start=1):
        gates += [('CX', first, p+2*j), ('CX', second, p+2*j+1)]
        if row[coordinate] == -1:
            gates.append(('CX', second, p+2*j))
    return gates


def _native_source(p, coordinate):
    assert p >= 1 and 0 <= coordinate < p
    first, second, dirty = p, p+1, p+(1 << p)
    wall = _clifford_wall(p, coordinate)
    fredkin = ([('C', [('CX', first, second)])]
               + _compact_toffoli(coordinate, second, first, dirty)
               + [('C', [('CX', first, second)])])
    return dirty+1, [('C', wall)] + fredkin + [('C', wall)]


def _source_columns(p, coordinate):
    bank_size, width = 1 << p, p+(1 << p)+1
    result = zeros(1 << width, 1 << width)
    for column in range(1 << width):
        bank = (column >> p) & ((1 << bank_size)-1)
        changed = (_representative(p, coordinate, bank)
                   if column >> coordinate & 1 else bank)
        result[column ^ ((bank ^ changed) << p)][column] = CO
    return result


def _echo_columns(p, coordinate):
    width = p+(1 << p)+1
    result = zeros(1 << width, 1 << width)
    for column in range(1 << width):
        mask = _target_mask(p, column & (1 << coordinate))
        result[column ^ (mask << p)][column] = CO
    return result


def _binary_rank(vectors):
    pivots = {}
    for vector in vectors:
        while vector:
            pivot = vector.bit_length()-1
            if pivot not in pivots:
                pivots[pivot] = vector
                break
            vector ^= pivots[pivot]
    return len(pivots)


class WeakPauliOrbitCompressionTests(unittest.TestCase):
    def test_native_source_echo_all_columns_and_actual_inverse(self):
        # p=2, i=1 exercises both sign orientations, an untouched address
        # bit, every bank column, and both values of the returned helper.
        for p, coordinate in ((1, 0), (2, 1)):
            width, source = _native_source(p, coordinate)
            actual = native_columns(width, _word(source))
            self.assertEqual(actual, _source_columns(p, coordinate))
            self.assertEqual(native_columns(width, _word(_inverse(source))),
                             adjoint(actual))
            pauli = [('C', [gate for j in range(1 << (p-1))
                            for gate in [('H', p+2*j), ('S', p+2*j),
                                         ('S', p+2*j), ('H', p+2*j)]])]
            # Chronology E^-1,P0,E has operator E P0 E^-1.
            echo = _inverse(source) + pauli + source
            self.assertEqual(native_columns(width, _word(echo)),
                             _echo_columns(p, coordinate))
            self.assertEqual(native_columns(width, _word(_inverse(echo))),
                             adjoint(_echo_columns(p, coordinate)))

    def test_native_ledger_and_physical_bands(self):
        for p in range(1, 8):
            q = 1 << (p-1)
            for coordinate in range(p):
                width, source = _native_source(p, coordinate)
                negative = 0 if coordinate == 0 else q//2
                cliffs = 4*q+16+2*negative
                self.assertEqual(width, p+2*q+1)
                self.assertEqual(len(_clifford_wall(p, coordinate)),
                                 2*(q-1)+negative)
                for schedule in (source, _inverse(source)):
                    self.assertEqual(_counts(_word(schedule)), (8, cliffs))
                    bands = [gates for kind, gates in schedule if kind == 'T']
                    self.assertEqual(len(bands), 2)
                    for band in bands:
                        self.assertEqual(len(band), 4)
                        self.assertEqual(len({gate[1] for gate in band}), 4)
                    self.assertTrue(all(0 <= wire < width
                                        for gate in _word(schedule)
                                        for wire in gate[1:]))

    def test_every_small_address_and_all_size_failure_witness(self):
        for p in range(1, 8):
            v = _input_mask(p)
            for coordinate in range(p):
                w = _target_mask(p, 1 << coordinate)
                self.assertEqual(_representative(p, coordinate, v), w)
                self.assertEqual(_representative(p, coordinate, w), v)
                displacement = [_representative(p, coordinate, 1 << j)
                                ^ (1 << j) for j in range(1 << p)]
                self.assertEqual(_binary_rank(displacement), 1)
            for address in range(1 << p):
                actual, affine = v, v
                for coordinate in range(p):
                    if address >> coordinate & 1:
                        actual = _representative(p, coordinate, actual)
                        affine ^= _displacement(p, coordinate)
                self.assertEqual(actual, affine)
                if address.bit_count() <= 1:
                    self.assertEqual(actual, _target_mask(p, address))
            if p >= 2:
                actual = _representative(p, 1, _representative(p, 0, v))
                expected = sum((1 if row[1] == 1 else 0) << (2*j)
                               for j, row in enumerate(_signs(p)))
                self.assertEqual(actual, expected)
                target = _target_mask(p, 3)
                self.assertNotEqual(actual, target)
                # A negative character of their nonzero difference supplies
                # an X-eigenvector with eigenvalues of opposite signs.
                difference = actual ^ target
                witness = difference & -difference
                self.assertEqual((difference & witness).bit_count() % 2, 1)

    def test_prefix_differences_have_the_required_independent_rank(self):
        for p in range(1, 8):
            orbit = [_target_mask(p, address) for address in range(1 << p)]
            self.assertEqual(_binary_rank(orbit), 1 << p)
            total = 0
            for coordinate in range(p):
                previous = 1 << coordinate
                differences = [orbit[y] ^ orbit[y+previous]
                               for y in range(previous)]
                self.assertEqual(_binary_rank(differences), previous)
                total += previous
            self.assertEqual(total, (1 << p)-1)


if __name__ == '__main__':
    unittest.main()
