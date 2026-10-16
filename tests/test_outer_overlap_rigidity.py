"""Exact bounded audits for outer-overlap rigidity and its scope.

These fixtures check signs, native layers, occupied phases and retained
counterexamples. They do not establish an infinite-width operator gap.
"""
from fractions import Fraction
from itertools import combinations, product
import unittest

import numpy as np

from tests.test_one_clean_inner_spread import (
    _ONE, _ZERO, _INV_SQRT2, _q2, _add, _mul, _scale,
    _pauli_product, _rotate, _conjugate, _multiply,
    _counterexample_expansions, _native_counterexample,
)
from tests.test_operator_source_compiler import _adjoint
from tests.test_native_representation_limits import _Cyclotomic


def _label(width, wire, letter):
    chars = ['I'] * width
    chars[width - 1 - wire] = letter
    return ''.join(chars)


def _rank(rows):
    pivots = {}
    for row in rows:
        while row:
            pivot = row.bit_length() - 1
            if pivot in pivots:
                row ^= pivots[pivot]
            else:
                pivots[pivot] = row
                break
    return len(pivots)


def _binary(label):
    width = len(label)
    x = z = 0
    for wire, letter in enumerate(reversed(label)):
        if letter in 'XY':
            x |= 1 << wire
        if letter in 'ZY':
            z |= 1 << wire
    return x | (z << width)


def _pair(left, right):
    return _pauli_product(left, right)[1] % 2


def _independent_isotropic(axes):
    return (_rank([_binary(axis) for axis in axes]) == len(axes)
            and all(_pair(a, b) == 0 for a in axes for b in axes))


def _cx_label(label, control, target):
    width = len(label)
    result, exponent = 'I' * width, 0
    for wire, letter in enumerate(reversed(label)):
        if letter == 'Y':
            exponent += 1  # Y = i X Z.
        for factor in ('X', 'Z'):
            if not (letter == factor or letter == 'Y'):
                continue
            image = _label(width, wire, factor)
            if factor == 'X' and wire == control:
                image, extra = _pauli_product(image, _label(width, target, 'X'))
                exponent += extra
            if factor == 'Z' and wire == target:
                image, extra = _pauli_product(image, _label(width, control, 'Z'))
                exponent += extra
            result, extra = _pauli_product(result, image)
            exponent += extra
    exponent %= 4
    assert exponent in (0, 2)
    return result, 1 if exponent == 0 else -1


def _clifford_label(label, gate):
    name, wire = gate[:2]
    if name == 'CX':
        return _cx_label(label, wire, gate[2])
    position = len(label) - 1 - wire
    letter = label[position]
    tables = {
        'H': {'I': ('I', 1), 'X': ('Z', 1), 'Y': ('Y', -1), 'Z': ('X', 1)},
        'S': {'I': ('I', 1), 'X': ('Y', 1), 'Y': ('X', -1), 'Z': ('Z', 1)},
        'SDG': {'I': ('I', 1), 'X': ('Y', -1), 'Y': ('X', 1), 'Z': ('Z', 1)},
        'X': {'I': ('I', 1), 'X': ('X', 1), 'Y': ('Y', -1), 'Z': ('Z', -1)},
        'Z': {'I': ('I', 1), 'X': ('X', -1), 'Y': ('Y', -1), 'Z': ('Z', 1)},
    }
    mapped, sign = tables[name][letter]
    return label[:position] + mapped + label[position + 1:], sign


def _word_conjugate(expansion, word):
    for gate in word:
        if gate[0] in ('T', 'TDG'):
            width = len(next(iter(expansion)))
            expansion = _rotate(expansion, _label(width, gate[1], 'Z'),
                                1 if gate[0] == 'T' else -1)
        else:
            result = {}
            for label, value in expansion.items():
                mapped, sign = _clifford_label(label, gate)
                result[mapped] = _add(result.get(mapped, _ZERO), _scale(value, sign))
            expansion = {label: value for label, value in result.items() if value != _ZERO}
    return expansion


def _is_signed_pauli(expansion):
    return len(expansion) == 1 and next(iter(expansion.values())) in (_ONE, _q2(-1))


def _is_clifford_word(width, word):
    return all(_is_signed_pauli(_word_conjugate({_label(width, q, p): _ONE}, word))
               for q in range(width) for p in ('X', 'Z'))


def _signed_permutes(axes, word):
    labels = []
    for axis in axes:
        image = _word_conjugate({axis: _ONE}, word)
        if not _is_signed_pauli(image):
            return False
        labels.append(next(iter(image)))
    return set(labels) == set(axes) and len(labels) == len(set(labels))


def _linear_z_word(columns, width):
    """Synthesize exact Z-label map by CNOT row elimination."""
    rows = [sum(((columns[col] >> row) & 1) << col for col in range(width))
            for row in range(width)]
    operations = []
    for pivot in range(width):
        if not (rows[pivot] >> pivot) & 1:
            other = next(row for row in range(pivot + 1, width)
                         if (rows[row] >> pivot) & 1)
            for dst, src in ((pivot, other), (other, pivot), (pivot, other)):
                rows[dst] ^= rows[src]
                operations.append(('CX', dst, src))
        for row in range(width):
            if row != pivot and (rows[row] >> pivot) & 1:
                rows[row] ^= rows[pivot]
                operations.append(('CX', row, pivot))
    assert rows == [1 << row for row in range(width)]
    return list(reversed(operations))


def _tensor(left, right):
    return {a + b: _mul(x, y) for a, x in left.items() for b, y in right.items()}


def _tensor_power(expansion, count):
    result = {'': _ONE}
    for _ in range(count):
        result = _tensor(result, expansion)
    return result


def _real_product(left, right):
    result = _multiply(left, right)
    assert all(imag == _ZERO for real, imag in result.values())
    return {label: real for label, (real, _) in result.items() if real != _ZERO}


def _projector(width, axes, signs):
    result = {'I' * width: _ONE}
    for axis, sign in zip(axes, signs):
        result = _real_product(result, {'I' * width: _q2(Fraction(1, 2)),
                                        axis: _q2(Fraction(sign, 2))})
    return result


def _sandwich(expansion, axis, sign):
    result = {}
    for label, coefficient in expansion.items():
        if _pair(axis, label):
            continue
        mapped, exponent = _pauli_product(axis, label)
        assert exponent in (0, 2)
        half = _scale(coefficient, Fraction(1, 2))
        result[label] = _add(result.get(label, _ZERO), half)
        result[mapped] = _add(result.get(mapped, _ZERO),
                              _scale(half, sign * (1 if exponent == 0 else -1)))
    return {label: value for label, value in result.items() if value != _ZERO}


def _cancel_adjacent(word):
    result = []
    for gate in word:
        if result and result[-1] == _adjoint([gate])[0]:
            result.pop()
        else:
            result.append(gate)
    return result


def _lift_word(word, offset):
    return [(gate[0], *(q + offset for q in gate[1:])) for gate in word]


def _parallel_r77(count):
    _, _, layers = _native_counterexample()
    clean = 1
    common = [('CX', 2 * block + 1, clean) for block in range(1, count)]
    words = []
    physical_layers = []
    for initial, selected in (([('H', 1), ('S', 1), ('H', 0)], layers[:2]),
                              ([('H', 1), ('S', 1)], layers[2:])):
        word = list(common)
        for block in range(count):
            word += _lift_word(initial, 2 * block)
        actual_layers = []
        for routing, layer in selected:
            routes = sum((_lift_word(routing, 2 * block) for block in range(count)), [])
            gates = sum((_lift_word(layer, 2 * block) for block in range(count)), [])
            word += _adjoint(routes) + gates + routes
            actual_layers.append(gates)
        words.append(word)
        physical_layers.append(actual_layers)
    return words, physical_layers, common


def _exact_word_matrix(width, word):
    """Literal native word, reusing the retained Q(sqrt(2),i) ring."""
    q = _Cyclotomic
    size = 1 << width
    matrix = np.array([[q(int(i == j)) for j in range(size)] for i in range(size)],
                      dtype=object)
    halfroot = q(0, Fraction(1, 2))
    omega = q(0, Fraction(1, 2), 0, Fraction(1, 2))
    phases = {'Z': q(-1), 'S': q(0, 0, 1), 'SDG': q(0, 0, -1),
              'T': omega, 'TDG': omega.conjugate()}
    for gate in word:
        name = gate[0]
        if name in ('X', 'CX'):
            images = [basis ^ (1 << gate[-1])
                      if all((basis >> control) & 1 for control in gate[1:-1])
                      else basis for basis in range(size)]
            matrix = matrix[np.argsort(images)]
        elif name == 'H':
            mask = 1 << gate[1]
            for low in range(size):
                if low & mask:
                    continue
                high = low | mask
                row0, row1 = matrix[low].copy(), matrix[high].copy()
                matrix[low] = [(a + b) * halfroot for a, b in zip(row0, row1)]
                matrix[high] = [(a - b) * halfroot for a, b in zip(row0, row1)]
        else:
            for row in range(size):
                if row >> gate[1] & 1:
                    matrix[row] = [value * phases[name] for value in matrix[row]]
    return matrix


def _power(value, exponent):
    if exponent < 0:
        return _power(value.conjugate(), -exponent)
    result = _Cyclotomic(1)
    for _ in range(exponent):
        result *= value
    return result


class OuterOverlapRigidityTests(unittest.TestCase):
    def test_exact_clifford_conjugator_and_linear_synthesis(self):
        # Small generator identities independently calibrate the routing helper.
        self.assertEqual(_cx_label('XI', 1, 0), ('XX', 1))
        self.assertEqual(_cx_label('IZ', 1, 0), ('ZZ', 1))
        self.assertEqual(_cx_label('XZ', 1, 0), ('YY', -1))
        self.assertEqual(_clifford_label('Y', ('H', 0)), ('Y', -1))
        self.assertEqual(_clifford_label('Y', ('SDG', 0)), ('X', 1))
        for width, columns in ((3, [7, 3, 5]), (4, [14, 13, 11, 7])):
            word = _linear_z_word(columns, width)
            for wire, column in enumerate(columns):
                expected = ''.join('Z' if column >> q & 1 else 'I'
                                   for q in reversed(range(width)))
                self.assertEqual(_word_conjugate({_label(width, wire, 'Z'): _ONE}, word),
                                 {expected: _ONE})

    def test_signed_permutation_inverses_and_inactive_qubits(self):
        width = 3
        axes = ['IIZ', 'IZI']
        swap = [('CX', 0, 1), ('CX', 1, 0), ('CX', 0, 1)]
        cases = [[], [('X', 0)], swap, swap + [('X', 0)], [('S', 1)],
                 [('H', 2)], [('CX', 0, 2)], [('CX', 2, 0)],
                 [('CX', 0, 1)], [('H', 0)]]
        for signs in product((1, -1), repeat=2):
            layer = [('T' if sign == 1 else 'TDG', wire)
                     for wire, sign in enumerate(signs)]
            for clifford in cases:
                with self.subTest(signs=signs, clifford=clifford):
                    conjugated = _adjoint(layer) + clifford + layer
                    self.assertEqual(_is_clifford_word(width, conjugated),
                                     _signed_permutes(axes, clifford))
        # Changing the common coordinates preserves the statement, including signs.
        routing = [('H', 0), ('S', 0), ('CX', 0, 2), ('H', 1)]
        routed_axes = [next(iter(_word_conjugate({axis: _ONE}, routing))) for axis in axes]
        self.assertTrue(_independent_isotropic(routed_axes))
        for clifford in cases:
            routed_clifford = _adjoint(routing) + clifford + routing
            routed_layer = _adjoint(routing) + [('T', 0), ('TDG', 1)] + routing
            conjugated = _adjoint(routed_layer) + routed_clifford + routed_layer
            self.assertEqual(_is_clifford_word(width, conjugated),
                             _signed_permutes(routed_axes, routed_clifford))

    def test_pair_overlap_conditions_do_not_replace_triple_condition(self):
        # A=I+all-ones at size four is orthogonal and invertible over F2,
        # but its distinct triple overlaps are odd; it is not a permutation.
        columns = [14, 13, 11, 7]
        self.assertEqual(_rank(columns), 4)
        self.assertTrue(all(column.bit_count() % 2 == 1 for column in columns))
        self.assertTrue(all((columns[i] & columns[j]).bit_count() % 2 == 0
                            for i in range(4) for j in range(i)))
        self.assertEqual((columns[0] & columns[1] & columns[2]).bit_count() % 2, 1)
        clifford = _linear_z_word(columns, 4)
        layer = [('T', wire) for wire in range(4)]
        self.assertFalse(_signed_permutes([_label(4, wire, 'Z') for wire in range(4)],
                                         clifford))
        self.assertFalse(_is_clifford_word(4, _adjoint(layer) + clifford + layer))

    def test_transverse_nonclifford_midpoint_literal_routing_and_gauge(self):
        # Input Z axes: wire2 -> XXZ, wire1 -> ZZZ, wire0 -> YYZ.
        # The X0 corrects the sign B(Z2 Z1 Z0)=-YYZ after the Bell route.
        route = [('X', 0)] + _linear_z_word([7, 3, 5], 3) + [('H', 2), ('CX', 2, 1)]
        axes_by_wire = {2: 'XXZ', 1: 'ZZZ', 0: 'YYZ'}
        for wire, axis in axes_by_wire.items():
            self.assertEqual(_word_conjugate({_label(3, wire, 'Z'): _ONE}, route),
                             {axis: _ONE})
        self.assertTrue(_independent_isotropic(list(axes_by_wire.values())))
        layer = [('T', 0), ('T', 1), ('T', 2)]
        dword = _adjoint(route) + layer + route
        core = _word_conjugate({'IIX': _ONE}, dword)
        expected = {pair + last: _scale(_INV_SQRT2, Fraction(1, 2))
                    for pair in ('II', 'XX', 'YY', 'ZZ') for last in ('X', 'Y')}
        self.assertEqual(core, expected)  # SWAP_12 tensor (X+Y)_3/sqrt(2).
        self.assertEqual(core, _conjugate({'IIX': _ONE},
                                        [('XXZ', 1), ('ZZZ', 1), ('YYZ', 1)]))
        a = [('H', 0)] + dword + [('T', 2)]
        b = [('H', 0)] + dword + [('TDG', 1)]
        midpoint = _word_conjugate({'IIZ': _ONE}, a)
        self.assertEqual(midpoint, _word_conjugate({'IIZ': _ONE}, b))
        midpoint_word = _adjoint(a) + [('Z', 0)] + a
        self.assertFalse(_is_clifford_word(3, midpoint_word))
        self.assertEqual(_rank([_binary('ZII'), _binary('IZI')]), 2)
        gauged_a, gauged_b = a + [('TDG', 2)], b + [('TDG', 2)]
        self.assertEqual(_word_conjugate({'IIZ': _ONE}, gauged_a), core)
        self.assertEqual(_word_conjugate({'IIZ': _ONE}, gauged_b), core)
        self.assertTrue(_is_clifford_word(3, _adjoint(gauged_a) + [('Z', 0)] + gauged_a))
        # Gauge leaves the complete circuit unchanged and returns the sole clean wire.
        original = a + _adjoint(b)
        gauged = gauged_a + _adjoint(gauged_b)
        # Exact word cancellation also certifies the common scalar phase.
        self.assertEqual(_cancel_adjacent(original), _cancel_adjacent(gauged))
        for wire in range(3):
            for letter in ('X', 'Z'):
                pauli = {_label(3, wire, letter): _ONE}
                self.assertEqual(_word_conjugate(pauli, original),
                                 _word_conjugate(pauli, gauged))
        self.assertEqual(_word_conjugate({'IIZ': _ONE}, original), {'IIZ': _ONE})
        # The last two TDG gates are one layer, on distinct physical wires.
        self.assertEqual(len({gate[1] for gate in gauged_b[-2:]}), 2)

    def test_r77_literal_layer_routes_are_exact(self):
        _, _, layers = _native_counterexample()
        expected = [(('YY', 1), ('ZX', 0)), (('ZY', 1), ('IY', 0)),
                    (('ZX', 0),), (('ZZ', 1), ('IZ', 0))]
        for (route, layer), images in zip(layers, expected):
            self.assertEqual(len(layer), len({gate[1] for gate in layer}))
            axes = []
            for label, wire in images:
                self.assertEqual(_word_conjugate({_label(2, wire, 'Z'): _ONE}, route),
                                 {label: _ONE})
                axes.append(label)
            self.assertTrue(_independent_isotropic(axes))

    def test_r77_tensor_parity_one_clean_shared_image_and_rank(self):
        base_s, other_s, _, _ = _counterexample_expansions()
        self.assertEqual(base_s, other_s)
        for count in (1, 2, 3):
            width = 2 * count
            words, layers, common = _parallel_r77(count)
            clean = _label(width, 1, 'Z')
            center_axes = [_label(width, 2 * block + 1, 'Z') for block in range(count)]
            parity = 'ZI' * count
            self.assertEqual(_word_conjugate({clean: _ONE}, common), {parity: _ONE})
            expected = _tensor_power(base_s, count)
            for word in words:
                self.assertEqual(_word_conjugate({clean: _ONE}, word), expected)
            self.assertEqual(_word_conjugate({clean: _ONE}, words[0] + _adjoint(words[1])),
                             {clean: _ONE})
            for pair in layers:
                for layer in pair:
                    self.assertEqual(len(layer), len({gate[1] for gate in layer}))
                    self.assertTrue(all(gate[0] in ('T', 'TDG') for gate in layer))
            a_axes = []
            b_axes = []
            for block in range(count):
                for local in ('IY', 'ZY'):
                    a_axes.append('II' * (count - 1 - block) + local + 'II' * block)
                for local in ('IZ', 'ZZ'):
                    b_axes.append('II' * (count - 1 - block) + local + 'II' * block)
            self.assertTrue(_independent_isotropic(a_axes))
            self.assertTrue(_independent_isotropic(b_axes))
            dim_intersection = (_rank([_binary(x) for x in a_axes])
                                + _rank([_binary(x) for x in b_axes])
                                - _rank([_binary(x) for x in a_axes + b_axes]))
            crossing_rows = [sum(_pair(a, b) << j for j, b in enumerate(b_axes))
                             for a in a_axes]
            self.assertEqual(dim_intersection, count)
            self.assertEqual(_rank(crossing_rows), count)
            # The common S flips, rather than fixes, every center generator.
            sword = _adjoint(words[0]) + [('Z', 1)] + words[0]
            for axis in center_axes:
                self.assertEqual(_word_conjugate({axis: _ONE}, sword), {axis: _q2(-1)})

    def test_r77_center_mixing_exact_small_parity_families(self):
        words, _, _ = _parallel_r77(1)
        a_observable = _adjoint(words[0]) + [('Z', 0)] + words[0]
        b_observable = _adjoint(words[1]) + [('Z', 0)] + words[1]
        self.assertEqual(_word_conjugate({'ZI': _ONE}, a_observable),
                         {'XY': _INV_SQRT2, 'YI': _scale(_INV_SQRT2, -1)})
        self.assertEqual(_word_conjugate({'ZI': _ONE}, b_observable), {'ZI': _ONE})
        for count in (1, 2):
            width = 2 * count
            words, _, _ = _parallel_r77(count)
            axes = [_label(width, 2 * block + 1, 'Z') for block in range(count)]
            data_parity_word = [('Z', 2 * block) for block in range(count)]
            aword = _adjoint(words[0]) + data_parity_word + words[0]
            for signs_z in product((1, -1), repeat=count):
                image = _word_conjugate(_projector(width, axes, signs_z), aword)
                for signs_w in product((1, -1), repeat=count):
                    compressed = image
                    for axis, sign in zip(axes, signs_w):
                        compressed = _sandwich(compressed, axis, sign)
                    expected = {label: _scale(value, Fraction(1, 1 << count))
                                for label, value in _projector(width, axes, signs_w).items()}
                    self.assertEqual(compressed, expected)

    def test_r77_displayed_complexity_has_pauli_preserving_resplit(self):
        # Growing displayed intersection/matching does not obstruct a different
        # native representation. This literal three-layer word fixes clean Z.
        a, b, layers = _native_counterexample()
        route_zx = layers[2][0]
        alternative = ([('X', 0), ('S', 0)]
                       + _adjoint(route_zx) + [('TDG', 0)] + route_zx
                       + [('TDG', 0)] + [('H', 0), ('T', 0), ('H', 0)])
        omega = _Cyclotomic(0, Fraction(1, 2), 0, Fraction(1, 2))
        original_matrix = _exact_word_matrix(2, a + _adjoint(b))
        alternative_matrix = _exact_word_matrix(2, alternative)
        for index in np.ndindex(original_matrix.shape):
            self.assertEqual(original_matrix[index], omega * alternative_matrix[index])
        self.assertEqual(_word_conjugate({'ZI': _ONE}, alternative), {'ZI': _ONE})
        self.assertTrue(all(_pair(axis, 'ZI') == 0 for axis in ('ZX', 'IZ', 'IX')))
        self.assertEqual(sum(gate[0] in ('T', 'TDG') for gate in alternative), 3)

    def test_exact_return_does_not_imply_same_space_phase_polynomial_descent(self):
        # Three data wires and one initialized clean wire. The actual native
        # middle word is omega times R_XXXX R_YYYY; that phase is retained.
        q = _Cyclotomic
        omega = q(0, Fraction(1, 2), 0, Fraction(1, 2))
        routing = [('CX', 2, 1), ('CX', 3, 1), ('CX', 0, 1), ('H', 0),
                   ('CX', 0, 1), ('CX', 0, 2), ('CX', 0, 3)]
        self.assertEqual(_word_conjugate({'IIIZ': _ONE}, routing), {'XXXX': _ONE})
        self.assertEqual(_word_conjugate({'IIZI': _ONE}, routing), {'YYYY': _ONE})
        self.assertTrue(_independent_isotropic(['XXXX', 'YYYY']))
        cat = [('H', 3), ('CX', 3, 0), ('CX', 3, 1), ('CX', 3, 2)]
        self.assertEqual(_word_conjugate({'ZIII': _ONE}, cat), {'XXXX': _ONE})
        layers = [[('T', wire) for wire in range(4)], [('T', 0), ('T', 1)],
                  [('TDG', wire) for wire in range(4)]]
        word = (cat + layers[0] + _adjoint(routing) + layers[1] + routing
                + layers[2] + _adjoint(cat))
        for layer in layers:
            self.assertEqual(len(layer), len({gate[1] for gate in layer}))
        self.assertEqual(_word_conjugate({'ZIII': _ONE}, word), {'ZIII': _ONE})
        matrix = _exact_word_matrix(4, word)
        determinant = q(1)
        for source in range(8):
            bits = [(source >> wire) & 1 for wire in range(3)]
            exponent = 2 - sum(bits) + 2 * bits[0] * bits[1] * bits[2]
            expected = _power(omega, exponent)
            determinant *= expected
            for target in range(16):
                self.assertEqual(matrix[target, source], expected if target == source else q())
        self.assertEqual(determinant, q(0, 0, -1))
        # Check the determinant obstruction on the generating phase monomials.
        for width in (3, 4):
            exponents = {1 << width}  # A constant eighth-root phase.
            for degree in (1, 2, 3):
                for indices in combinations(range(width), degree):
                    exponent = sum((1 << (degree - 1))
                                   * all((bits >> index) & 1 for index in indices)
                                   for bits in range(1 << width))
                    exponents.add(exponent)
            residues = {exponent % 8 for exponent in exponents}
            self.assertEqual(residues, {0, 4} if width == 3 else {0})
        self.assertEqual(determinant * determinant, q(-1))  # One appended dirty.
        self.assertEqual(_power(determinant, 4), q(1))  # Two remove this obstruction only.
        for phase in range(8):
            self.assertEqual(_power(_power(omega, phase), 8), q(1))
            self.assertEqual(_power(_power(omega, phase), 16), q(1))


if __name__ == '__main__':
    unittest.main(verbosity=2)
