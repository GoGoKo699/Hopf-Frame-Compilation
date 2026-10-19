"""Exact small checks for complete-indicator unitary stabilizer nullity.

The source-count proof is analytic. These fixtures keep all native
columns and inverse phases, including arbitrary dirty Pauli directions.
The cheap shared band is deliberately not an indicator implementation.
"""
from functools import lru_cache
import unittest

from tests.exact_matrices import CO, CZ, adjoint, native_columns
from tests.test_majorana_walsh_query import _inverse
from tests.test_native_indicator_sharing import _affine_rank, _propagate
from tests.test_phase_signal_query_limits import _cz, _phase_query_word
from tests.test_three_layer_native_mechanisms import _native_t_bands


def _affine_mask(values):
    """Return the linear mask when a complete binary table is affine."""
    width = (len(values) - 1).bit_length()
    constant = values[0]
    mask = sum((values[1 << bit] ^ constant) << bit for bit in range(width))
    if all(value == (constant ^ ((mask & index).bit_count() & 1))
           for index, value in enumerate(values)):
        return mask
    return None


@lru_cache(None)
def _native_target(p, dirty):
    if p == 1:
        word = [('S', 1)] * 2 + _cz(0, 1) + _cz(0, 2)
    elif p == 2:
        word = _phase_query_word()
    else:
        raise ValueError('only the small literal native fixtures are used')
    width = p + (1 << p) + dirty
    actual = native_columns(width, word)
    phases = []
    for column in range(1 << width):
        for row in range(1 << width):
            if row != column:
                assert actual[row][column] == CZ
        assert actual[column][column] in (CO, -CO)
        phases.append(int(actual[column][column] == -CO))
    return word, actual, tuple(phases)


def _preserved_x_masks(phases):
    return [mask for mask in range(len(phases))
            if _affine_mask([phase ^ phases[column ^ mask]
                             for column, phase in enumerate(phases)]) is not None]


class IndicatorPauliNullityTests(unittest.TestCase):
    def test_all_native_target_paulis_include_both_dirty_directions(self):
        for p, dirty in ((1, 0), (2, 0), (2, 1)):
            with self.subTest(p=p, dirty=dirty):
                m, width = 1 << p, p + (1 << p) + dirty
                word, actual, phases = _native_target(p, dirty)
                self.assertEqual(native_columns(width, _inverse(word)), adjoint(actual))
                for column, phase in enumerate(phases):
                    address = column & (m - 1)
                    bank = (column >> p) & ((1 << m) - 1)
                    self.assertEqual(phase, (bank >> address) & 1)
                preserved_x = _preserved_x_masks(phases)
                if p == 1:
                    self.assertEqual(preserved_x, list(range(1 << width)))
                    continue
                expected_x = []
                for mask in range(1 << width):
                    address_mask = mask & (m - 1)
                    bank_mask = (mask >> p) & ((1 << m) - 1)
                    bank_table = [(bank_mask >> x) & 1 for x in range(m)]
                    if address_mask == 0 and _affine_mask(bank_table) is not None:
                        expected_x.append(mask)
                self.assertEqual(preserved_x, expected_x)
                # Every Z direction is independent; both X and Z on the
                # arbitrary dirty wire contribute to the actual subgroup.
                all_labels = [x | (z << width) for x in preserved_x
                              for z in range(1 << width)]
                dimension = m + 2 * p + 1 + 2 * dirty
                self.assertEqual(_affine_rank(all_labels), dimension)
                self.assertEqual(len(all_labels), 1 << dimension)
                self.assertEqual(2 * width - dimension, m - 1)

    def test_nonpauli_native_images_have_the_claimed_signed_gap(self):
        # A diagonal target's image of X^a retains that exact permutation.
        # An input Z mask only adds another Pauli character, so these
        # comparisons cover its gap for every input Pauli Z mask as well.
        _, actual, phases = _native_target(2, 0)
        width, size = 6, 64
        preserved = set(_preserved_x_masks(phases))
        saw_imaginary_equality = False
        for x_mask in range(size):
            if x_mask in preserved:
                continue
            derivative = [phase ^ phases[column ^ x_mask]
                          for column, phase in enumerate(phases)]
            self.assertIsNone(_affine_mask(derivative))
            for z_mask in range(size):
                # The actual Hermitian Pauli with these masks has the
                # literal coefficient i^|x&z| (-1)^(z.column).
                scalar = (x_mask & z_mask).bit_count() % 4
                for sign in (0, 1):
                    squares = []
                    for column in range(size):
                        phase_exponent = (scalar + 2 * sign
                                          + 2 * ((z_mask & column).bit_count() & 1)) % 4
                        real_part = (1, 0, -1, 0)[phase_exponent]
                        squares.append(2 - 2 * (-1 if derivative[column] else 1)
                                       * real_part)
                    self.assertGreaterEqual(max(squares), 2)
                    if scalar % 2:
                        self.assertEqual(set(squares), {2})
                        saw_imaginary_equality = True
            # A different candidate X mask has an orthogonal destination
            # on every column, independent of all literal phases.
            column, different = 0, x_mask ^ 1
            self.assertNotEqual(column ^ x_mask, column ^ different)
            self.assertIn(actual[column][column], (CO, -CO))
        self.assertTrue(saw_imaginary_equality)

    def test_one_native_layer_saturates_incidence_but_misses_joint_action(self):
        p, m, width = 2, 4, 6
        encode = [('CX', p + u, j) for u in range(m) for j in range(p)]
        word = encode + [('T', j) for j in range(p)] + _inverse(encode)
        actual = native_columns(width, word)
        omega = native_columns(1, [('T', 0)])[1][1]
        roots = [CO]
        for _ in range(7):
            roots.append(roots[-1] * omega)
        for column in range(1 << width):
            address, bank = column & (m - 1), column >> p
            parity = bank.bit_count() & 1
            exponent = sum(((address >> j) & 1) ^ parity for j in range(p))
            for row in range(1 << width):
                self.assertEqual(actual[row][column], roots[exponent] if row == column else CZ)
        self.assertEqual(native_columns(width, _inverse(word)), adjoint(actual))
        self.assertEqual(sum(gate[0] in ('T', 'TDG') for gate in word), p)
        self.assertEqual(sum(gate[0] not in ('T', 'TDG') for gate in word), 2 * p * m)
        self.assertEqual([len(band) for band in _native_t_bands(word)], [p])

        probes = [{1 << (p + u): CO} for u in range(m)]
        incidence = []
        for gate in word:
            active = 0
            for u in range(m):
                before = probes[u]
                probes[u] = _propagate(before, width, gate)
                if gate[0] in ('T', 'TDG') and probes[u] != before:
                    active += 1
            if gate[0] in ('T', 'TDG'):
                incidence.append(active)
        self.assertEqual(incidence, [m] * p)
        self.assertTrue(all(_affine_rank(probe) == p for probe in probes))

        # This joint Pauli commutes with every P_j axis. The complete
        # target instead gives it a non-Pauli image (nonzero address X).
        joint = ((1 << p) - 1) | (1 << p)
        propagated = {joint: CO}
        for gate in word:
            propagated = _propagate(propagated, width, gate)
        self.assertEqual(propagated, {joint: CO})
        _, target, phases = _native_target(p, 0)
        self.assertNotIn(joint, _preserved_x_masks(phases))
        wrong_column = (1 << p) | (1 << (p + 1))
        self.assertEqual(actual[wrong_column][wrong_column], CO)
        self.assertEqual(target[wrong_column][wrong_column], -CO)


if __name__ == '__main__':
    unittest.main()
