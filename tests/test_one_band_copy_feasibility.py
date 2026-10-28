"""Exact finite checks of the complete one-band copy feasibility reduction.

The projected Pauli-transfer fixed space is the full compatible Hermitian
seed space. All-output dual equations alone can retain a wrong carrier
action. These fixtures are not an all-size seed/copy construction.
"""
from fractions import Fraction
import unittest

from tests.exact_matrices import (
    Q2, add, adjoint, multiply, native_columns, scale, zeros,
)
from tests.test_field_shear_intertwiner import _pencil_columns
from tests.test_joint_lookup_pauli_source import _bands
from tests.test_native_copy_interaction import _one_band_copy
from tests.test_pauli_orbit_copy import _one_address_copy, _pauli


_ZERO, _ONE = Q2(0), Q2(1)
_ISQRT2 = Q2(0, Fraction(1, 2))


def _inverse(word):
    result = []
    for name, *wires in reversed(word):
        if name == 'S':
            result += [('S', *wires)] * 3
        elif name in ('T', 'TDG'):
            result.append(('TDG' if name == 'T' else 'T', *wires))
        else:
            result.append((name, *wires))
    return result


def _add_to(row, index, value):
    result = row.get(index, _ZERO) + value
    if result != _ZERO:
        row[index] = result
    else:
        row.pop(index, None)


def _gate_image(label, value, gate, width):
    """Literal conjugation on tau(x,z)=i**(x.z) X**x Z**z."""
    name, *wires = gate
    mask = (1 << width) - 1
    x, z = label & mask, label >> width
    wire = wires[0]
    xb, zb = (x >> wire) & 1, (z >> wire) & 1
    if name in ('T', 'TDG'):
        if not xb:
            return {label: value}
        sign = -1 if zb ^ (name == 'TDG') else 1
        return {label: value * _ISQRT2,
                label ^ (1 << (width + wire)): sign * value * _ISQRT2}
    sign = 1
    if name == 'H':
        sign = -1 if xb and zb else 1
        if xb != zb:
            x ^= 1 << wire
            z ^= 1 << wire
    elif name == 'S':
        sign = -1 if xb and zb else 1
        if xb:
            z ^= 1 << wire
    elif name == 'Z':
        sign = -1 if xb else 1
    elif name == 'X':
        sign = -1 if zb else 1
    elif name == 'Y':
        sign = -1 if xb ^ zb else 1
    elif name == 'CX':
        control, target = wires
        before = (x & z).bit_count()
        if (x >> control) & 1:
            x ^= 1 << target
        if (z >> target) & 1:
            z ^= 1 << control
        exponent = (before - (x & z).bit_count()) % 4
        assert exponent in (0, 2)
        sign = 1 if exponent == 0 else -1
    else:
        raise ValueError(gate)
    return {x | (z << width): sign * value}


def _conjugate(terms, word, width):
    for gate in word:
        result = {}
        for label, value in terms.items():
            for other, coefficient in _gate_image(label, value, gate, width).items():
                _add_to(result, other, coefficient)
        terms = result
    return terms


def _lift(label, carrier, width):
    return (label & ((1 << carrier) - 1)) | ((label >> carrier) << width)


def _fixed_point_rows(carrier, outputs, copy, queries):
    """Carrier-only rows of the complete copy fixed-point equations.

    Retain the physical output Pauli X_address after conjugation by C,
    then pull the carrier coefficient back by Q_address. Equality to the
    seed saturates the Hilbert--Schmidt projection norm; discarded output
    components must therefore vanish. No clean-input restriction is used.
    """
    width = carrier + outputs
    carrier_mask = (1 << carrier) - 1
    for address, query in enumerate(queries):
        rows = {}
        for column in range(1 << (2 * carrier)):
            source = _conjugate({column: _ONE}, query, carrier)
            source = {_lift(label, carrier, width): value
                      for label, value in source.items()}
            image = _conjugate(source, copy, width)
            kept = {}
            for label, value in image.items():
                x, z = label & ((1 << width) - 1), label >> width
                if x >> carrier == 1 << address and z >> carrier == 0:
                    target = (x & carrier_mask) | ((z & carrier_mask) << carrier)
                    _add_to(kept, target, value)
            kept = _conjugate(kept, _inverse(query), carrier)
            _add_to(kept, column, -_ONE)
            for row, value in kept.items():
                _add_to(rows.setdefault(row, {}), column, value)
        yield from rows.values()


def _kernel(rows, size):
    """Sparse exact elimination over Q(sqrt(2)), without rank tolerances."""
    pivots = {}
    for original in rows:
        row = dict(original)
        while row:
            pivot = min(row)
            value = row[pivot]
            if pivot not in pivots:
                pivots[pivot] = {column: entry / value
                                 for column, entry in row.items()}
                break
            for column, entry in pivots[pivot].items():
                _add_to(row, column, -value * entry)
    basis = []
    for free in range(size):
        if free in pivots:
            continue
        vector = {free: _ONE}
        for pivot, row in sorted(pivots.items(), reverse=True):
            value = -sum((entry * vector.get(column, _ZERO)
                          for column, entry in row.items() if column != pivot),
                         _ZERO)
            if value != _ZERO:
                vector[pivot] = value
        basis.append(vector)
    return basis


def _complete_seed_space(carrier, outputs, copy, queries):
    return _kernel(_fixed_point_rows(carrier, outputs, copy, queries),
                   1 << (2 * carrier))


def _pauli_product(left, right, width):
    mask = (1 << width) - 1
    x, z = left & mask, left >> width
    u, v = right & mask, right >> width
    phase = ((x & z).bit_count() + (u & v).bit_count()
             - ((x ^ u) & (z ^ v)).bit_count() + 2 * (z & u).bit_count()) % 4
    return (x ^ u) | ((z ^ v) << width), phase


def _tableau(word, width):
    images = []
    for index in range(2 * width):
        image = _conjugate({1 << index: _ONE}, word, width)
        assert len(image) == 1
        (label, value), = image.items()
        assert value in (_ONE, -_ONE)
        images.append((label, 0 if value == _ONE else 2))
    return images


def _tableau_image(label, tableau, width):
    mask = (1 << width) - 1
    phase = ((label & mask) & (label >> width)).bit_count() % 4
    result = 0
    for index, (image, sign) in enumerate(tableau):
        if (label >> index) & 1:
            result, extra = _pauli_product(result, image, width)
            phase = (phase + sign + extra) % 4
    assert phase in (0, 2)
    return result, 0 if phase == 0 else 1


def _dual_seed_space(carrier, outputs, copy, queries):
    """Necessary all-output X/Z equations as signed permutation components."""
    width, size = carrier + outputs, 1 << (2 * carrier)
    parents, parity, bad = list(range(size)), [0] * size, [False] * size

    def find(label):
        if parents[label] != label:
            root, sign = find(parents[label])
            parity[label] ^= sign
            parents[label] = root
        return parents[label], parity[label]

    def require(left, right, sign):
        a, sa = find(left)
        b, sb = find(right)
        if a == b:
            bad[a] |= (sa ^ sb) != sign
        else:
            parents[b] = a
            parity[b] = sa ^ sb ^ sign
            bad[a] |= bad[b]

    for address, query in enumerate(queries):
        for output in range(outputs):
            for axis in ('X', 'Z'):
                word = (query + copy + [(axis, carrier + output)]
                        + _inverse(copy) + _inverse(query))
                tableau = _tableau(word, width)
                wanted = int(axis == 'Z' and output == address)
                for label in range(size):
                    target, sign = _tableau_image(
                        _lift(label, carrier, width), tableau, width)
                    x, z = target & ((1 << width) - 1), target >> width
                    if x >> carrier or z >> carrier:
                        root, _ = find(label)
                        bad[root] = True
                    else:
                        require(label, x | (z << carrier), wanted ^ sign)
    groups = {}
    for label in range(size):
        root, sign = find(label)
        if not bad[root]:
            groups.setdefault(root, {})[label] = -_ONE if sign else _ONE
    return list(groups.values())


def _matrix(terms, width):
    result = zeros(1 << width)
    for label, value in terms.items():
        result = add(result, scale(_pauli(width, label), value))
    return result


def _product_copy(size):
    """Reschedule disjoint pair components into one literal physical band."""
    prefix, phases, suffix = [], [], []
    for carrier in range(size):
        word = _one_band_copy(carrier, size + carrier)
        first = next(index for index, gate in enumerate(word)
                     if gate[0] in ('T', 'TDG'))
        assert all(gate[0] in ('T', 'TDG') for gate in word[first:first + 2])
        prefix += word[:first]
        phases += word[first:first + 2]
        suffix += word[first + 2:]
    return prefix + phases + suffix


class OneBandCopyFeasibilityTests(unittest.TestCase):
    def test_signed_pauli_transfer_matches_every_native_small_column(self):
        words = (
            (1, [('H', 0), ('T', 0), ('S', 0), ('TDG', 0), ('H', 0)]),
            (2, [('H', 0), ('S', 1), ('CX', 0, 1), ('T', 0),
                 ('TDG', 1), ('H', 1), ('CX', 1, 0)]),
            (3, _one_address_copy(0, 1, 2)),
        )
        for width, word in words:
            unitary = native_columns(width, word)
            for label in range(1 << (2 * width)):
                actual = multiply(multiply(unitary, _pauli(width, label)),
                                  adjoint(unitary))
                self.assertEqual(actual, _matrix(
                    _conjugate({label: _ONE}, word, width), width))

    def test_known_one_address_copy_has_the_complete_hadamard_seed_space(self):
        copy, queries = _one_address_copy(0, 1, 2), [[], [('Z', 0)]]
        expected = [{1: _ONE, 2: _ONE}]
        self.assertEqual(_complete_seed_space(1, 2, copy, queries), expected)
        self.assertEqual(_dual_seed_space(1, 2, copy, queries), expected)
        self.assertEqual(_matrix({1: _ISQRT2, 2: _ISQRT2}, 1),
                         native_columns(1, [('H', 0)]))
        self.assertEqual(native_columns(3, _inverse(copy)),
                         adjoint(native_columns(3, copy)))
        self.assertEqual([len(band) for band in _bands(copy)], [2])

    def test_all_output_duals_can_miss_a_nontrivial_carrier_action(self):
        copy, queries = _one_address_copy(0, 1, 2), [[], [('Z', 0)]]
        changed = copy + [('Z', 0)]
        self.assertEqual(_dual_seed_space(1, 2, changed, queries),
                         _dual_seed_space(1, 2, copy, queries))
        self.assertEqual(_complete_seed_space(1, 2, changed, queries), [])
        first, second = native_columns(3, copy), native_columns(3, changed)
        for output in (1, 2):
            for label in (1 << output, 1 << (3 + output)):
                probe = _pauli(3, label)
                self.assertEqual(multiply(multiply(adjoint(first), probe), first),
                                 multiply(multiply(adjoint(second), probe), second))
        seed = native_columns(3, [('H', 0)])
        wrong = multiply(multiply(second, seed), adjoint(second))
        target = multiply(seed, _pauli(3, 1 << 1))
        defect = add(wrong, scale(target, -1))
        self.assertEqual(multiply(adjoint(defect), defect),
                         scale(_pauli(3, 0), 2))

    def test_arbitrary_carrier_helper_retains_all_four_pauli_directions(self):
        copy, queries = _one_address_copy(0, 2, 3), [[], [('Z', 0)]]
        basis = _complete_seed_space(2, 2, copy, queries)
        expected = [{1 ^ helper: _ONE, 4 ^ helper: _ONE}
                    for helper in (0, 2, 8, 10)]
        self.assertCountEqual(basis, expected)
        self.assertCountEqual(_dual_seed_space(2, 2, copy, queries), expected)

    def test_four_address_product_copy_has_no_complete_nonzero_seed(self):
        copy = _product_copy(4)
        bands = _bands(copy)
        self.assertEqual([len(band) for band in bands], [8])
        self.assertEqual({gate[1] for gate in bands[0]}, set(range(8)))
        queries = [
            [('CX', source, 2 + target)
             for source, column in enumerate(_pencil_columns(2, address))
             for target in range(2) if (column >> target) & 1]
            for address in range(4)
        ]
        self.assertEqual(_complete_seed_space(4, 4, copy, queries), [])
        self.assertEqual(_dual_seed_space(4, 4, copy, queries), [])


if __name__ == '__main__':
    unittest.main()
