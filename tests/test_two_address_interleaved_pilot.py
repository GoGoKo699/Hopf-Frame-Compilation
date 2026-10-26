"""A complete two-address controlled-H pilot with a paid Clifford decoder.

Every phase and dirty mask is retained. This finite source validates the
interleaved family contract; it is not a uniform shallow lookup compiler.
"""
from itertools import product
import unittest

from tests.exact_matrices import CO, CZ, adjoint, native_columns
from tests.test_common_probe_source import _native
from tests.test_complete_observable_bases import _controlled_h
from tests.test_dirty_cubic_phase_compiler import _gray_phase_batch
from tests.test_joint_lookup_pauli_source import _bands, _inverse


ADDRESS = (8, 9)
MASKS = (10, 11, 12)
CANONICAL = ((0, 1, 0), (0, 0, 1), (0, 2, 1), (0, 1, 2))


def _multiply(left, right):
    """Literal i^k X(x) Z(z), with all X factors before all Z factors."""
    k, x, z = left
    ell, a, b = right
    return ((k + ell + 2 * (z & a).bit_count()) % 4, x ^ a, z ^ b)


def _conjugate(pauli, word):
    k, x, z = pauli
    for name, *wires in word:
        if name == 'H':
            wire, = wires
            a, b = (x >> wire) & 1, (z >> wire) & 1
            k = (k + 2 * a * b) % 4
            x ^= (a ^ b) << wire
            z ^= (a ^ b) << wire
        elif name == 'S':
            wire, = wires
            a = (x >> wire) & 1
            k = (k + a) % 4
            z ^= a << wire
        elif name == 'CX':
            control, target = wires
            x ^= ((x >> control) & 1) << target
            z ^= ((z >> target) & 1) << control
        else:
            raise AssertionError(name)
    return k, x, z


def _decoder():
    # Chronological signed symplectic reduction E; D is its actual adjoint.
    reduction = [
        ('S', 0), ('S', 0), ('S', 0), ('S', 4), ('CX', 0, 4),
        ('H', 1), ('S', 4), ('H', 4), ('S', 5), ('H', 5),
        ('CX', 1, 0), ('CX', 4, 0), ('CX', 5, 0),
        ('H', 0), ('S', 0), ('H', 0), ('H', 1), ('H', 5),
        ('S', 6), ('CX', 1, 2), ('CX', 1, 5), ('CX', 1, 6),
        ('H', 3), ('S', 6), ('H', 6), ('S', 7), ('H', 7),
        ('CX', 3, 1), ('CX', 6, 1), ('CX', 7, 1),
        ('H', 1), ('S', 1), ('S', 1), ('S', 1), ('H', 1),
    ]
    return _inverse(reduction)


def _parallel_controlled_h(control, targets):
    words = [_controlled_h(control, target) for target in targets]
    partitions = []
    for word in words:
        chunks = [[]]
        previous_t = False
        for gate in word:
            is_t = gate[0] in ('T', 'TDG')
            if chunks[-1] and is_t != previous_t:
                chunks.append([])
            chunks[-1].append(gate)
            previous_t = is_t
        partitions.append(chunks)
    assert all(len(chunks) == 5 for chunks in partitions)
    # Whole Clifford chunks keep S/S-dagger inverses adjacent and paid.
    return [gate for i in range(5) for chunks in partitions for gate in chunks[i]]


def _w():
    c = [('CX', 1, 0), ('H', 1)]
    second = _inverse(c) + _controlled_h(ADDRESS[1], 0) + c
    return second + _parallel_controlled_h(ADDRESS[0], (0, 1))


def _query_parts():
    terms = ((ADDRESS[0], 1, 4), (ADDRESS[1], 2, 4), (ADDRESS[0], 3, 6))
    phase = _gray_phase_batch(terms, MASKS, [1, 1, 1], 3)
    wall = [('H', 4), ('H', 6)]
    query = [('CX', j, 4 + j) for j in range(4)] + wall + phase + list(reversed(wall))
    return query, phase, terms


def _bell_wall():
    return [('H', j) for j in range(4)] + [('CX', j, 4 + j) for j in range(4)]


def _source():
    query, _, _ = _query_parts()
    return [('H', 0)] + _w() + _decoder() + query + _bell_wall()


def _query_fiber(address):
    word = [('CX', j, 4 + j) for j in range(4)]
    for subset in (1, 2, 3):
        bit = (subset & -subset).bit_length() - 1
        if address >> bit & 1:
            word.append(('CX', subset, 4 + (subset ^ (1 << bit))))
    return word


def _target_paulis():
    result = []
    for address in range(4):
        v = sum(1 << subset for subset in range(4) if subset & address == subset)
        b = _query_fiber(address) + _bell_wall()
        result.append(_conjugate((0, 0, v << 4), _inverse(b)))
    return tuple(result)


class TwoAddressInterleavedPilotTests(unittest.TestCase):
    def test_literal_decoder_maps_the_whole_noncommuting_family(self):
        targets = _target_paulis()
        self.assertEqual(targets, ((0, 17, 17), (0, 35, 49),
                                   (0, 69, 81), (0, 191, 241)))
        decoder = _decoder()
        self.assertEqual(tuple(_conjugate(p, decoder) for p in CANONICAL), targets)
        frame = (targets[0], targets[1], _multiply(targets[1], targets[2]),
                 _multiply(targets[0], targets[3]))
        self.assertEqual(tuple(_conjugate(p, _inverse(decoder)) for p in frame),
                         ((0, 1, 0), (0, 0, 1), (0, 2, 0), (0, 0, 2)))
        for p, target in zip(CANONICAL, targets):
            self.assertEqual(_conjugate(target, _inverse(decoder)), p)
        # This checks all products, not just four individually balanced spectra.
        for powers in product((0, 1), repeat=4):
            left = right = (0, 0, 0)
            for include, p, target in zip(powers, CANONICAL, targets):
                if include:
                    left, right = _multiply(left, p), _multiply(right, target)
            self.assertEqual(_conjugate(left, decoder), right)
            if any(powers):
                self.assertNotEqual((right[1], right[2]), (0, 0))

    def test_controlled_h_word_all_columns_and_actual_inverse(self):
        rename = {0: 0, 1: 1, ADDRESS[0]: 2, ADDRESS[1]: 3}
        word = [(g[0], *(rename[q] for q in g[1:])) for g in _w()]
        matrix = native_columns(4, word)
        expected = [[CZ for _ in range(16)] for _ in range(16)]
        c = [('CX', 1, 0), ('H', 1)]
        for address in range(4):
            fiber = (_inverse(c) + [('H', 0)] + c) if address & 2 else []
            if address & 1:
                fiber += [('H', 0), ('H', 1)]
            block = native_columns(2, fiber)
            for row, column in product(range(4), repeat=2):
                expected[row | (address << 2)][column | (address << 2)] = block[row][column]
        self.assertEqual(matrix, expected)
        self.assertEqual(native_columns(4, _inverse(word)), adjoint(expected))

    def test_complete_probe_on_all_four_address_fibers(self):
        c = [('CX', 1, 0), ('H', 1)]
        for address in range(4):
            w = (_inverse(c) + [('H', 0)] + c) if address & 2 else []
            if address & 1:
                w += [('H', 0), ('H', 1)]
            fiber = [('H', 0)] + w + _decoder() + _query_fiber(address) + _bell_wall()
            v = sum(1 << subset for subset in range(4) if subset & address == subset)
            target = (0, 0, v << 4)
            self.assertEqual(_conjugate((0, 0, 1), fiber), target)
            self.assertEqual(_conjugate(target, _inverse(fiber)), (0, 0, 1))

    def test_all_active_cubic_inputs_and_arbitrary_masks_return(self):
        _, phase, terms = _query_parts()
        active = sorted(set(q for term in terms for q in term) | set(MASKS))
        self.assertEqual(len(active), 10)
        for packed in range(1 << len(active)):
            basis = sum(((packed >> i) & 1) << wire for i, wire in enumerate(active))
            exponent = sum(all(basis >> wire & 1 for wire in term) for term in terms)
            expected = {basis: (-1) ** exponent * CO}
            for word in (phase, _inverse(phase)):
                self.assertEqual(_native(word, {basis: CO}), expected)

    def test_paid_common_native_ledger(self):
        decoder, w, source = _decoder(), _w(), _source()
        self.assertEqual((len(decoder), len(_inverse(decoder))), (45, 35))
        self.assertEqual((len(w), len(_inverse(w))), (55, 55))
        self.assertEqual((len(source), len(_inverse(source))), (165, 155))
        for word in (source, _inverse(source)):
            self.assertEqual(sum(g[0] in ('T', 'TDG') for g in word), 30)
            self.assertEqual(len(_bands(word)), 12)
            for band in _bands(word):
                self.assertEqual(len({g[1] for g in band}), len(band))
        probe = _inverse(source) + [('S', 0), ('S', 0)] + source
        self.assertEqual((len(probe), len(_inverse(probe))), (322, 322))
        self.assertEqual(sum(g[0] in ('T', 'TDG') for g in probe), 60)
        self.assertEqual(len(_bands(probe)), 24)
        self.assertEqual(max(q for g in source for q in g[1:]), 12)


if __name__ == '__main__':
    unittest.main()
