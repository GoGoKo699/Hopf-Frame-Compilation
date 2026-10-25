"""Exact native checks for globally shared-angle bipartite edge layers.

The dirty cubic compiler pays each layer. Two invariant, coherent input
sectors nevertheless have identical source action and opposite lookup
signs. These tests support that specified grammar screen, not a native
depth lower bound or a prohibition on other bank/address dynamics.
"""
import unittest

from tests.exact_matrices import CI, CO, CZ
from tests.test_common_probe_source import ROOT_HALF, _native
from tests.test_dirty_cubic_phase_compiler import _parity_edge_layer
from tests.test_joint_lookup_pauli_source import _inverse


def _edge_word(p, subset, coordinate, signed):
    """Independent literal two-rotation implementation of one edge."""
    a, b = p + subset, p + (subset | (1 << coordinate))
    wall = [('H', a)]
    parity = [('CX', a, b)]
    if signed:
        parity = [('CX', coordinate, a)] + parity
        band = [('T', a), ('TDG', b)]
    else:
        band = [('TDG', a), ('T', b)]
    return wall + parity + band + list(reversed(parity)) + wall


def _separate_edge_layer(p, parity, signed):
    return [gate for subset in range(1 << p)
            if subset.bit_count() % 2 == parity
            for coordinate in range(p) if not (subset >> coordinate) & 1
            for gate in _edge_word(p, subset, coordinate, signed)]


def _face_state(virtual, kind):
    """Embed (probe, arbitrary helper/reference payload) in the p=2 face.

    A: address (1,0), d=0, b=1,c=0.
    B: address (1,1), d=1, singlet on b,c; omit its common sqrt(2).
    Logical wires are x0,x1,a,b,c,d. Payload begins at wire six.
    """
    result = {}
    for (probe, payload), value in virtual.items():
        if kind == 'A':
            result[1 | ((probe | 2) << 2) | (payload << 6)] = value
        else:
            result[3 | ((probe | 10) << 2) | (payload << 6)] = value
            result[3 | ((probe | 12) << 2) | (payload << 6)] = -value
    return result


def _quarter_rotation(virtual, orientation=1):
    """exp(orientation*i*pi*X/4), retaining its literal scalar phase."""
    result = {}
    for (probe, payload), value in virtual.items():
        for bit, amplitude in ((probe, ROOT_HALF * value),
                               (probe ^ 1, orientation * CI * ROOT_HALF * value)):
            key = bit, payload
            result[key] = result.get(key, CZ) + amplitude
    return {key: value for key, value in result.items() if value != CZ}


def _lookup_sign(state, p):
    result = {}
    for basis, value in state.items():
        address = basis & ((1 << p) - 1)
        bank = basis >> p
        parity = sum((bank >> subset) & 1 for subset in range(1 << p)
                     if subset & address == subset) & 1
        result[basis] = -value if parity else value
    return result


def _calibrated_word(first_parity):
    """Matrix word E_x O_x O E, or O_x E_x E O, in chronological order."""
    specifications = [(first_parity, False), (1 - first_parity, False),
                      (1 - first_parity, True), (first_parity, True)]
    pieces = [_parity_edge_layer(2, parity, signed)
              for parity, signed in specifications]
    return [gate for word, _ in pieces for gate in word], max(width for _, width in pieces)


class BipartiteShearMixingTests(unittest.TestCase):
    def test_paid_layers_match_independent_native_edges_and_return_dirty_work(self):
        # Sparse columns retain occupied helpers without allocating a full matrix.
        columns = (0, 1, 2, 3, 7, 10, 17, 23, 29, 42, 55, 63)
        for parity in (0, 1):
            for signed in (False, True):
                word, width = _parity_edge_layer(2, parity, signed)
                separate = _separate_edge_layer(2, parity, signed)
                dirty = (1 << (width - 6)) - 1
                for index, principal in enumerate(columns):
                    helper = dirty if index % 3 == 0 else dirty // 3 if index % 3 == 1 else 0
                    basis = principal | (helper << 6)
                    expected = _native(separate, {basis: CO})
                    actual = _native(word, {basis: CO})
                    self.assertEqual(actual, expected)
                    self.assertEqual(_native(_inverse(word), actual), {basis: CO})

    def test_paired_sectors_agree_for_every_generator_and_actual_inverse(self):
        # Ten helper bits cover either signed p=2 layer; bit sixteen is an
        # untouched reference correlated with the probe and the occupied pool.
        virtual = {(0, 0): CO, (1, (1 << 11) - 1): CI}
        for parity in (0, 1):
            for signed in (False, True):
                word, _ = _parity_edge_layer(2, parity, signed)
                for orientation, native in ((1, word), (-1, _inverse(word))):
                    expected = (_quarter_rotation(virtual, orientation)
                                if parity == 0 else virtual)
                    for kind in ('A', 'B'):
                        self.assertEqual(_native(native, _face_state(virtual, kind)),
                                         _face_state(expected, kind))

    def test_calibrated_native_words_have_maximal_complete_probe_error(self):
        for first_parity in (0, 1):
            word, width = _calibrated_word(first_parity)
            # Probe zero, an arbitrary occupied helper pool, and a coherent
            # helper/reference branch. The target fixes this whole witness.
            virtual = {(0, 0): CO, (0, (1 << (width - 6 + 1)) - 1): CI}
            witness = _face_state(virtual, 'B')
            self.assertEqual(_lookup_sign(witness, 2), witness)
            echo = _inverse(word) + [('S', 2), ('S', 2)] + word
            self.assertEqual(_native(echo, witness),
                             {basis: -value for basis, value in witness.items()})
            self.assertEqual(_native(word + _inverse(word), witness), witness)
            # Normalizing gives norm error exactly two, not an average-column gap.

    def test_zero_address_calibration_preserves_arbitrary_bank_columns(self):
        for first_parity in (0, 1):
            word, width = _calibrated_word(first_parity)
            dirty = ((1 << (width - 6)) - 1) << 6
            for bank in (0, 1, 3, 5, 10, 15):
                basis = (bank << 2) | dirty
                self.assertEqual(_native(word, {basis: CO}), {basis: CO})

    def test_nonuniform_edge_weights_break_the_paired_sector_hypothesis(self):
        # Keeping only a->b does not annihilate the singlet's individual b
        # response. The all-word screen therefore cannot include arbitrary
        # coordinate-dependent weights or changed bank interlayers.
        single_edge = _edge_word(2, 0, 0, False)
        virtual = {(0, 0): CO}
        self.assertEqual(_native(single_edge, _face_state(virtual, 'A')),
                         _face_state(_quarter_rotation(virtual), 'A'))
        self.assertNotEqual(_native(single_edge, _face_state(virtual, 'B')),
                            _face_state(_quarter_rotation(virtual), 'B'))

    def test_two_coordinate_faces_embed_at_all_sizes_and_targets_are_opposite(self):
        for p in range(2, 9):
            for i, j in {(0, 1), (0, p - 1), (p - 2, p - 1)}:
                face = {0, 1 << i, 1 << j, (1 << i) | (1 << j)}
                for parity in (0, 1):
                    edges = [(subset, subset | (1 << coordinate), coordinate)
                             for subset in range(1 << p)
                             if subset.bit_count() % 2 == parity
                             for coordinate in range(p)
                             if not (subset >> coordinate) & 1]
                    self.assertEqual(len(edges), p * (1 << p) // 4)
                    surviving = [(lower, upper, bit) for lower, upper, bit in edges
                                 if upper in face]
                    self.assertEqual(len(surviving), 2)
                    self.assertTrue(all(lower in face and bit in (i, j)
                                        for lower, _, bit in surviving))
                    self.assertTrue(all(lower != ((1 << i) | (1 << j))
                                        for lower, _, _ in surviving))
                    # Other controls remain zero: no source gate can create a
                    # coefficient outside the face from an active upper row.
                    self.assertTrue(all(upper not in face for lower, upper, _ in edges
                                        if lower not in face))
                for probe in (0, 1):
                    bank_a = probe | (1 << (1 << i))
                    state_a = {(1 << i) | (bank_a << p): CO}
                    self.assertEqual(_lookup_sign(state_a, p),
                                     {basis: (-value if not probe else value)
                                      for basis, value in state_a.items()})
                    top = 1 << ((1 << i) | (1 << j))
                    bank_b = probe | top | (1 << (1 << i))
                    bank_c = probe | top | (1 << (1 << j))
                    address = (1 << i) | (1 << j)
                    state_b = {address | (bank_b << p): CO,
                               address | (bank_c << p): -CO}
                    self.assertEqual(_lookup_sign(state_b, p),
                                     {basis: (-value if probe else value)
                                      for basis, value in state_b.items()})


if __name__ == '__main__':
    unittest.main()
