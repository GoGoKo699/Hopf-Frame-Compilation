"""Exact paid selector-tag fixtures and the scoped orbit-pair obstruction.

All logical bank/selector columns and native dirty-helper columns are kept.
The singlet tests concern a single covariant local-tag echo, not arbitrary
interleavings of different matchings or a native lookup-depth lower bound.
"""
from fractions import Fraction
from itertools import product
import unittest

from tests.exact_matrices import (
    CI, CO, CZ, CQ2, Q2, add, adjoint, eye, multiply, native_columns, scale,
)
from tests.test_co_moving_lookup_stage import _Matching, _ccz_7t
from tests.test_majorana_walsh_query import _inverse


_HALF = Fraction(1, 2)
_ROOT_HALF = CQ2(Q2(0, _HALF))
_I = eye(2)
_H = native_columns(1, [('H', 0)])
_SH = native_columns(1, [('H', 0), ('S', 0)])
_PHASES = (CO, CI, -CO, -CI)


def _accumulate(vector, key, value):
    vector[key] = vector.get(key, CZ) + value
    if vector[key] == CZ:
        del vector[key]


def _local(vector, row, matrix):
    result = {}
    for (x, bank, selectors), amplitude in vector.items():
        old = (selectors >> row) & 1
        for new in range(2):
            if matrix[new][old] != CZ:
                key = x, bank, selectors ^ ((old ^ new) << row)
                _accumulate(result, key, amplitude * matrix[new][old])
    return result


def _joint(vector, matching):
    return {(matching.images[x], matching.relabel(bank), matching.relabel(selectors)): value
            for (x, bank, selectors), value in vector.items()}


def _stage(vector, matching, inverse=False):
    result = {}
    for (x, bank, selectors), amplitude in vector.items():
        xx, yy, ss, borrowed, phase = matching.stage_literal(
            (x, bank, selectors, 1, 0), inverse=inverse)
        assert borrowed == 1
        _accumulate(result, (xx, yy, ss), amplitude * _PHASES[phase])
    return result


def _tag(vector, matching, matrices, inverse=False):
    """The same V_d on all rows at matching-control distance d."""
    matrices = tuple(adjoint(matrix) for matrix in matrices) if inverse else matrices
    result = {}
    for key, amplitude in vector.items():
        x, _, _ = key
        term = {key: amplitude}
        for row in range(matching.m):
            distance = sum(((x ^ row) >> control) & 1 for control in matching.controls)
            if matrices[distance] != _I:
                term = _local(term, row, matrices[distance])
        for output, value in term.items():
            _accumulate(result, output, value)
    return result


def _echo(vector, matching, matrices, inverse=False):
    if inverse:
        # Actual reverse chronology of B,R,B^dagger,R^dagger.
        vector = _stage(vector, matching)
        vector = _tag(vector, matching, matrices)
        vector = _stage(vector, matching, inverse=True)
        return _tag(vector, matching, matrices, inverse=True)
    vector = _tag(vector, matching, matrices)
    vector = _stage(vector, matching)
    vector = _tag(vector, matching, matrices, inverse=True)
    return _stage(vector, matching, inverse=True)


def _sector_factors(matching, matrices):
    factors = {}
    for sigma in (1, -1):
        factors[sigma] = []
        for distance, matrix in enumerate(matrices):
            q = [[CO, CZ], [CZ, _PHASES[distance % 4] * sigma]]
            factors[sigma].append(multiply(multiply(multiply(
                adjoint(q), adjoint(matrix)), q), matrix))
    return factors


def _sector_echo(vector, matching, factors):
    result = {}
    for sigma in (1, -1):
        branch = {}
        for key, amplitude in vector.items():
            x, bank, _ = key
            term = {key: amplitude}
            for row in range(matching.m):
                if (bank >> row) & 1:
                    distance = sum(((x ^ row) >> control) & 1
                                   for control in matching.controls)
                    term = _local(term, row, factors[sigma][distance])
            for output, value in term.items():
                _accumulate(branch, output, value)
        for key, value in branch.items():
            _accumulate(result, key, value * _HALF)
        for key, value in _joint(branch, matching).items():
            _accumulate(result, key, value * _HALF * sigma)
    return result


def _basis_change(selector):
    # Chronology for V=S H T H S^dagger, with literal positive-S native gates.
    return [('S', selector)] * 3 + [('H', selector), ('T', selector),
                                    ('H', selector), ('S', selector)]


def _cz(control, target):
    return [('H', target), ('CX', control, target), ('H', target)]


def _hadamard_tag_row(address, selector, constant):
    v = _basis_change(selector)
    return (_inverse(v) + _cz(address, selector) + v
            + ([('H', selector)] if constant else []))


def _controlled_tag_row(signal, address, selector, dirty1, dirty2, constant):
    v = _basis_change(selector)
    corner = _ccz_7t(dirty1, dirty2, selector)
    rectangle = (corner + [('CX', signal, dirty1)]
                 + corner + [('CX', address, dirty2)]
                 + corner + [('CX', signal, dirty1)]
                 + corner + [('CX', address, dirty2)])
    return (_inverse(v) + rectangle
            + (_cz(signal, selector) if constant else []) + v)


def _parallel_rows(rows):
    """One T per row per band; all shared wires are Clifford-only controls."""
    positions = [0] * len(rows)
    word, bands = [], []
    while any(position < len(row) for position, row in zip(positions, rows)):
        band = []
        for index, row in enumerate(rows):
            while positions[index] < len(row) and row[positions[index]][0] not in ('T', 'TDG'):
                word.append(row[positions[index]])
                positions[index] += 1
            if positions[index] < len(row):
                band.append(row[positions[index]])
                positions[index] += 1
        if band:
            bands.append(tuple(band))
            word.extend(band)
    return word, bands


def _expected_hadamards(width, columns, selector_controls):
    """Literal controlled-H action, with untouched address/signal/helpers."""
    matrix = [[CZ for _ in columns] for _ in range(1 << width)]
    for index, column in enumerate(columns):
        state = {column: CO}
        for selector, enabled in selector_controls(column):
            if not enabled:
                continue
            result = {}
            for key, value in state.items():
                old = (key >> selector) & 1
                for new in range(2):
                    _accumulate(result, key ^ ((old ^ new) << selector), value * _H[new][old])
            state = result
        for row, value in state.items():
            matrix[row][index] = value
    return matrix


class CoMovingSelectorReadoutTests(unittest.TestCase):
    def test_local_characteristic_polynomial_and_sixth_roots(self):
        factors = _sector_factors(_Matching(2, ((0, 1),)), (_I, _H))
        for sigma in (1, -1):
            self.assertEqual(factors[sigma][0], _I)
            g = factors[sigma][1]
            self.assertEqual(g[0][0] * g[1][1] - g[0][1] * g[1][0], CO)
            self.assertEqual(add(add(multiply(g, g), scale(g, -CO)), _I), scale(_I, CZ))
            cube = multiply(multiply(g, g), g)
            self.assertEqual(cube, scale(_I, -CO))
            self.assertEqual(multiply(cube, cube), _I)
            self.assertEqual(add(g, adjoint(g)), _I)
        # Odd distances3,5,... change the sign of the sector's S exponent;
        # a two-control fixture alone would never exercise this sign.
        matching = _Matching(6, ((0, 1), (2, 3), (4, 5)))
        many = _sector_factors(matching, (_I, _H, _I, _H))
        for sigma, distance in product((1, -1), range(4)):
            expected = (_I if distance % 2 == 0 else
                        factors[sigma * (-1)**((distance - 1) // 2)][1])
            self.assertEqual(many[sigma][distance], expected)

    def test_complete_two_bit_echo_sector_identity(self):
        for pairs in (((0, 1),), ((1, 0),)):
            matching = _Matching(2, pairs)
            factors = _sector_factors(matching, (_I, _H))
            for x, bank, selectors in product(range(4), range(16), range(16)):
                vector = {(x, bank, selectors): CO}
                self.assertEqual(_echo(vector, matching, (_I, _H)),
                                 _sector_echo(vector, matching, factors))
                # The local cube is -I on each active mismatch row.  Since
                # sector products commute with the joint involution, these
                # powered factors give the complete global cube/sixth power.
                for power in (3, 6):
                    powered = {sigma: [(_I if distance == 0 or power == 6
                                        else scale(_I, -CO)) for distance in range(2)]
                               for sigma in (1, -1)}
                    _, _, total, _ = matching.moments(x, bank, selectors)
                    expected = {(x, bank, selectors): CO * (-1 if power == 3 and total % 2 else 1)}
                    self.assertEqual(_sector_echo(vector, matching, powered), expected)

    def test_singleton_gram_and_literal_cube(self):
        for pairs in (((0, 1),), ((1, 0),)):
            matching = _Matching(2, pairs)
            control = pairs[0][0]
            for x, row, selectors in product(range(4), range(4), range(16)):
                vector = {(x, 1 << row, selectors): CO}
                forward = _echo(vector, matching, (_I, _H))
                reverse = _echo(vector, matching, (_I, _H), inverse=True)
                real = {}
                for result in (forward, reverse):
                    for key, value in result.items():
                        _accumulate(real, key, value * _HALF)
                mismatch = ((x ^ row) >> control) & 1
                self.assertEqual(real, {next(iter(vector)): CO * (_HALF if mismatch else 1)})
            # An independent literal product check includes multirow banks,
            # not only the powered-factor consequence above.
            for x, bank, selectors in product(range(4), (0, 1, 3, 10, 15), (0, 5, 15)):
                vector = {(x, bank, selectors): CO}
                cube = vector
                for _ in range(3):
                    cube = _echo(cube, matching, (_I, _H))
                _, _, total, _ = matching.moments(x, bank, selectors)
                self.assertEqual(cube, {next(iter(vector)): CO * (-1)**total})

    def test_orbit_pair_singlet_all_inactive_columns_and_noninvolutive_tags(self):
        for pairs in (((0, 1),), ((1, 0),)):
            matching = _Matching(2, pairs)
            u = 1 << pairs[0][0]
            v = matching.images[u]
            bank = (1 << u) | (1 << v)
            inactive = [row for row in range(4) if row not in (u, v)]
            marked, unmarked = 0, 0
            for matrix, x, other in product((_H, _SH), range(4), range(4)):
                selectors = sum(((other >> bit) & 1) << row for bit, row in enumerate(inactive))
                vector = {(x, bank, selectors | (1 << u)): _ROOT_HALF,
                          (x, bank, selectors | (1 << v)): -_ROOT_HALF}
                self.assertEqual(_echo(vector, matching, (_I, matrix)), vector)
                self.assertEqual(_echo(vector, matching, (_I, matrix), inverse=True), vector)
                if (bank >> x) & 1:
                    marked += 1
                else:
                    unmarked += 1
            self.assertEqual((marked, unmarked), (16, 16))

    def test_two_pair_matching_distance_tags(self):
        matching = _Matching(4, ((0, 1), (2, 3)))
        u, v = 5, 15
        bank = (1 << u) | (1 << v)
        # A noninvolutive V_2 on four rows, with V_0=V_1=I.  This probes
        # distance-dependent tags beyond the one-control H specialization.
        matrices = (_I, _I, _SH)
        factors = _sector_factors(matching, matrices)
        for x, selectors in product((0, u, v), (0, 1, (1 << 1) | (1 << 7))):
            vector = {(x, bank, selectors | (1 << u)): _ROOT_HALF,
                      (x, bank, selectors | (1 << v)): -_ROOT_HALF}
            self.assertEqual(_echo(vector, matching, matrices), vector)
            self.assertEqual(_sector_echo(vector, matching, factors), vector)
        # The parity-mismatch H tag used in the sixth-root theorem itself.
        # Keep a nonsymmetric inactive-selector input: no inactive L-sector
        # restriction is needed for the pair-singlet cancellation.
        matrices = (_I, _H, _I)
        factors = _sector_factors(matching, matrices)
        vector = {(1, bank, (1 << 1) | (1 << u)): _ROOT_HALF,
                  (1, bank, (1 << 1) | (1 << v)): -_ROOT_HALF}
        self.assertEqual(_echo(vector, matching, matrices), vector)
        self.assertEqual(_sector_echo(vector, matching, factors), vector)

    def test_native_hadamard_tags_and_affine_parity_all_columns(self):
        v = _basis_change(0)
        self.assertEqual(native_columns(1, _inverse(v) + [('S', 0)] * 2 + v), _H)
        for p, support, alpha in ((1, (0,), 0), (2, (0,), 0), (2, (0, 1), 1)):
            m = 1 << p
            pivot = support[0]
            parity = [('CX', wire, pivot) for wire in support[1:]]
            rows = [_hadamard_tag_row(pivot, p + row,
                                     alpha ^ (sum((row >> bit) & 1 for bit in support) % 2))
                    for row in range(m)]
            body, bands = _parallel_rows(rows)
            word = parity + body + list(reversed(parity))
            self.assertEqual(len(bands), 2)
            for band in bands:
                self.assertEqual(len({gate[1] for gate in band}), m)
                self.assertTrue(all(gate[1] >= p for gate in band))
            for direction in (word, _inverse(word)):
                self.assertEqual(sum(gate[0] in ('T', 'TDG') for gate in direction), 2 * m)
                self.assertEqual(sum(gate[0] not in ('T', 'TDG') for gate in direction),
                                 31 * m // 2 + 2 * (len(support) - 1))
            expected = _expected_hadamards(
                p + m, range(1 << (p + m)),
                lambda column: ((p + row, alpha ^
                                 (sum(((column ^ row) >> bit) & 1 for bit in support) % 2))
                                for row in range(m)))
            actual = native_columns(p + m, word)
            self.assertEqual(actual, expected)
            self.assertEqual(native_columns(p + m, _inverse(word)), adjoint(actual))

    def test_native_controlled_tags_all_dirty_columns_and_parallel_bands(self):
        for constant in (0, 1):
            word = _controlled_tag_row(0, 1, 2, 3, 4, constant)
            self.assertEqual(sum(gate[0] in ('T', 'TDG') for gate in word), 30)
            self.assertEqual(sum(gate[0] not in ('T', 'TDG') for gate in word), 40 + 3 * constant)
            self.assertTrue(all(gate[1] not in (0, 1) for gate in word if gate[0] in ('T', 'TDG')))
            self.assertTrue(all(gate[0] == 'CX' and gate[1] in (0, 1) and gate[2] >= 2
                                for gate in word if set(gate[1:]) & {0, 1}))
            expected = _expected_hadamards(
                5, range(32), lambda column: ((2, (column & 1) & (((column >> 1) & 1) ^ constant)),))
            actual = native_columns(5, word)
            self.assertEqual(actual, expected)
            self.assertEqual(native_columns(5, _inverse(word)), adjoint(actual))
        # Signal0/address1 are shared; selectors2/3 and helpers4..7 are
        # arbitrary.  Every T band acts on distinct private row wires.
        rows = [_controlled_tag_row(0, 1, 2 + row, 4 + 2 * row, 5 + 2 * row, row)
                for row in range(2)]
        word, bands = _parallel_rows(rows)
        self.assertEqual(len(bands), 30)
        for band in bands:
            self.assertEqual(len({gate[1] for gate in band}), 2)
            self.assertTrue(all(gate[1] not in (0, 1) for gate in band))
        for direction in (word, _inverse(word)):
            self.assertEqual(sum(gate[0] in ('T', 'TDG') for gate in direction), 60)
            self.assertEqual(sum(gate[0] not in ('T', 'TDG') for gate in direction), 83)
        actual = native_columns(8, word)
        expected = _expected_hadamards(
            8, range(256), lambda column: ((2 + row, (column & 1) & (((column >> 1) & 1) ^ row))
                                          for row in range(2)))
        self.assertEqual(actual, expected)
        self.assertEqual(native_columns(8, _inverse(word)), adjoint(actual))


if __name__ == '__main__':
    unittest.main()
