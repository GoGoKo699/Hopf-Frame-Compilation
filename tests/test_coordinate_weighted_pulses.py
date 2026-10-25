"""Exact checks for coordinate signatures in paid bipartite pulses.

Weights may differ between coordinates, but each coordinate has one weight
across its edges in a pulse. The all-size collision and packing bounds are
analytic pulse-word results, not unrestricted native T-depth lower bounds.
"""
from collections import Counter
import unittest

from tests.exact_matrices import CI, CO, CZ
from tests.test_bipartite_shear_mixing import (
    _edge_word, _lookup_sign, _quarter_rotation,
)
from tests.test_common_probe_source import ROOT_HALF, _native
from tests.test_dirty_cubic_phase_compiler import _dirty_cubic_phase
from tests.test_joint_lookup_pauli_source import _bands, _inverse


def _even_edge(p, subset, coordinate, signed, half_weight):
    """The literal Clifford part of k=2*l+e quarter-angle units."""
    a, b = p + subset, p + (subset | (1 << coordinate))
    power = half_weight % 4
    word = [('S', b)] * ((-power if signed else power) % 4)
    if power % 2:
        word.append(('CX', b, a))
        if signed:
            word += [('H', b), ('CX', coordinate, b), ('H', b)]
    return word


def _weighted_layer(p, parity, signed, weights):
    """Compile arbitrary integer quarter weights with one odd-edge batch.

    Return the literal word and total width. The even part is Clifford;
    the remaining subset uses the retained dirty cubic compiler once.
    """
    assert p >= 2 and parity in (0, 1) and len(weights) == p
    assert all(isinstance(weight, int) for weight in weights)
    size = 1 << p
    linear = [0] * (p + size)
    quadratic, cubic, even = [], [], []
    for subset in range(size):
        if subset.bit_count() % 2 != parity:
            continue
        for coordinate, weight in enumerate(weights):
            if (subset >> coordinate) & 1:
                continue
            a, b = p + subset, p + (subset | (1 << coordinate))
            half_weight, odd = divmod(weight % 8, 2)
            even += _even_edge(p, subset, coordinate, signed, half_weight)
            if not odd:
                continue
            if signed:
                linear[b] -= 1
                quadratic += [(coordinate, b, 1), (a, b, 1)]
                cubic.append((coordinate, a, b))
            else:
                linear[b] += 1
                quadratic.append((a, b, -1))
    phase, width = _dirty_cubic_phase(linear, quadratic, cubic)
    wall = [('H', p + subset) for subset in range(size)
            if subset.bit_count() % 2 == parity]
    return even + wall + phase + list(reversed(wall)), width


def _independent_layer(p, parity, signed, weights):
    """Separate two-rotation edges, without using the compiler polynomial."""
    word = []
    for subset in range(1 << p):
        if subset.bit_count() % 2 != parity:
            continue
        for coordinate, weight in enumerate(weights):
            if not (subset >> coordinate) & 1:
                word += _edge_word(p, subset, coordinate, signed) * (weight % 8)
    return word


def _embed_face(virtual, kind, p=3, pair=(0, 2)):
    """Normalized paired face; payload starts after the whole logical bank."""
    i, j = pair
    address = 1 << i
    result = {}
    for (probe, payload), value in virtual.items():
        base = probe | (1 << (1 << i))
        if kind == 'A':
            key = address | (base << p) | (payload << (p + (1 << p)))
            result[key] = value
        else:
            top = 1 << ((1 << i) | (1 << j))
            first = base | top
            second = probe | (1 << (1 << j)) | top
            prefix = address | (1 << j)
            for bank, sign in ((first, 1), (second, -1)):
                key = prefix | (bank << p) | (payload << (p + (1 << p)))
                result[key] = sign * ROOT_HALF * value
    return result


def _rotate(virtual, quarters):
    for _ in range(quarters % 8):
        virtual = _quarter_rotation(virtual)
    return virtual


def _norm_squared(state):
    return sum((value.conj() * value for value in state.values()), CZ)


class CoordinateWeightedPulseTests(unittest.TestCase):
    def test_even_clifford_split_retains_every_scalar_and_input_column(self):
        # One edge has three arbitrary logical inputs x,a,b. Testing the whole
        # period includes zero, negative representatives, and signed controls.
        for signed in (False, True):
            primitive = _edge_word(1, 0, 0, signed)
            for weight in range(8):
                even, odd = divmod(weight, 2)
                split = _even_edge(1, 0, 0, signed, even) + primitive * odd
                for basis in range(8):
                    actual = _native(split, {basis: CO})
                    self.assertEqual(actual, _native(primitive * weight, {basis: CO}))
                    self.assertEqual(_native(_inverse(split), actual), {basis: CO})
            self.assertEqual(_native(primitive * 8, {7: CO}), {7: CO})

    def test_weighted_compiler_returns_dirty_work_and_coherent_reference(self):
        weights = (9, -2)  # An odd residue and an independently paid even one.
        for parity in (0, 1):
            for signed in (False, True):
                word, width = _weighted_layer(2, parity, signed, weights)
                separate = _independent_layer(2, parity, signed, weights)
                pool = (1 << (width - 6)) - 1
                for index, principal in enumerate((0, 3, 10, 23, 42, 63)):
                    helper = (pool, pool // 3, 0)[index % 3]
                    basis = principal | (helper << 6)
                    actual = _native(word, {basis: CO})
                    self.assertEqual(actual, _native(separate, {basis: CO}))
                    self.assertEqual(_native(_inverse(word), actual), {basis: CO})
                # The last bit is an untouched reference, correlated with the
                # logical and occupied-helper branches; all relative phases stay.
                state = {3: ROOT_HALF,
                         42 | (pool << 6) | (1 << width): CI * ROOT_HALF}
                actual = _native(word, state)
                self.assertEqual(actual, _native(separate, state))
                self.assertEqual(_native(_inverse(word), actual), state)

    def test_all_quarter_weights_keep_the_original_nonclifford_budgets(self):
        for p in (2, 3, 5):
            size, edges = 1 << p, p * (1 << p) // 4
            for weights in ([7] * p, [2 * (i % 4) for i in range(p)],
                            [3 * i - 5 for i in range(p)]):
                for parity in (0, 1):
                    for signed in (False, True):
                        word, width = _weighted_layer(p, parity, signed, weights)
                        self.assertLessEqual(width - p - size, 5 * edges)
                        for native in (word, _inverse(word)):
                            counts = Counter(gate[0] for gate in native)
                            self.assertLessEqual(counts['T'] + counts['TDG'],
                                                 size + 32 * edges)
                            self.assertLessEqual(len(native), 6 * size + 71 * edges)
                            bands = _bands(native)
                            self.assertLessEqual(len(bands), 13)
                            for band in bands:
                                self.assertEqual(len({gate[1] for gate in band}), len(band))
                        if all(weight % 2 == 0 for weight in weights):
                            self.assertFalse(_bands(word))
                            self.assertEqual(width, p + size)

    def test_equal_coordinate_signatures_ignore_other_direction_weights(self):
        # Directions zero and two agree modulo eight; direction one does not.
        # The latter's entire bank support is outside the invariant face.
        virtual = {(0, 0): ROOT_HALF, (1, 3): CI * ROOT_HALF}
        for parity, signed, weights in ((0, False, (1, 2, 9)),
                                       (1, True, (3, -2, -5)),
                                       (0, True, (2, 7, 10)),
                                       (1, False, (2, 3, 10))):
            word = _independent_layer(3, parity, signed, weights)
            for orientation, native in ((1, word), (-1, _inverse(word))):
                expected = _rotate(virtual, orientation * weights[0]) if parity == 0 else virtual
                for kind in ('A', 'B'):
                    self.assertEqual(_native(native, _embed_face(virtual, kind)),
                                     _embed_face(expected, kind))

    def test_a_paid_colliding_word_saturates_the_square_root_two_gap(self):
        specifications = ((0, False, (1, 2, 9)), (1, True, (3, -2, -5)),
                          (0, True, (2, 7, 10)), (1, False, (2, 3, 10)))
        pieces = [_weighted_layer(3, parity, signed, weights)
                  for parity, signed, weights in specifications]
        word = [gate for piece, _ in pieces for gate in piece]
        width = max(width for _, width in pieces)
        occupied_with_reference = (1 << (width - 11 + 1)) - 1
        virtual = {(0, 0): ROOT_HALF, (0, occupied_with_reference): CI * ROOT_HALF}
        # Both faces see e^(3*pi*i*X/4), whose Z image is -Y. Their targets
        # are opposite Z observables. Each normalized input has error sqrt(2).
        echo = _inverse(word) + [('S', 3), ('S', 3)] + word
        for kind, sign in (('A', -1), ('B', 1)):
            state = _embed_face(virtual, kind)
            actual = _native(echo, state)
            expected_virtual = {(1, payload): -CI * value
                                for (_, payload), value in virtual.items()}
            self.assertEqual(actual, _embed_face(expected_virtual, kind))
            target = _lookup_sign(state, 3)
            self.assertEqual(target, {basis: sign * value for basis, value in state.items()})
            difference = {basis: actual.get(basis, CZ) - target.get(basis, CZ)
                          for basis in actual.keys() | target.keys()}
            self.assertEqual(_norm_squared(state), CO)
            self.assertEqual(_norm_squared(difference), 2 * CO)
            self.assertEqual(_native(word + _inverse(word), state), state)

    def test_a_selective_clifford_pulse_breaks_the_singlet_premise(self):
        # Only c->d, in the odd parity class, has angle pi/2. This costs no
        # T gates and maps the middle singlet to another Bell state. Thus a
        # bound on total pulse count is not a bound on non-Clifford pulses.
        virtual = {(0, 0): CO}
        state = _embed_face(virtual, 'B', p=2, pair=(0, 1))
        clifford = _even_edge(2, 2, 0, False, 1)
        expected = {3 | ((0b1110) << 2): CI * ROOT_HALF,
                    3 | ((0b1000) << 2): -CI * ROOT_HALF}
        self.assertFalse(_bands(clifford))
        actual = _native(clifford, state)
        self.assertEqual(actual, expected)
        self.assertNotEqual(actual, state)
        self.assertEqual(_native(_edge_word(2, 2, 0, False) * 2, state), expected)
        self.assertEqual(_native(_inverse(clifford), actual), state)


if __name__ == '__main__':
    unittest.main()
