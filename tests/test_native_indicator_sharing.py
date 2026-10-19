"""Exact finite checks of the native indicator sharing requirement.

The support calculation retains literal native phases, arbitrary dirty
inputs and temporary address motion. A seven-T joint bank fixture lies
below the separate-row bound; no unrestricted depth bound is inferred.
"""
from functools import lru_cache
import unittest

from tests.exact_matrices import CO, CZ, add, adjoint, eye, multiply, native_columns, scale
from tests.test_majorana_walsh_query import _inverse
from tests.test_phase_signal_query_limits import _phase_query_word, _x
from tests.test_two_layer_galois_return import _coefficients, _conjugate, _tau


@lru_cache(None)
def _local_images(name):
    width = 2 if name == 'CX' else 1
    gate = native_columns(width, [(name, *range(width))])
    return tuple(tuple((label, coefficient) for label, coefficient in enumerate(
        _coefficients(_conjugate(gate, _tau(source, width)), width))
                       if coefficient != CZ)
                 for source in range(1 << (2 * width)))


def _propagate(poly, width, gate):
    """Use exact one/two-wire native matrices, including every Pauli sign."""
    name, *wires = gate
    local_width = len(wires)
    clear = sum((1 << wire) | (1 << (width + wire)) for wire in wires)
    result = {}
    for label, coefficient in poly.items():
        local = sum(((label >> wire) & 1) << index
                    for index, wire in enumerate(wires))
        local |= sum(((label >> (width + wire)) & 1) << (local_width + index)
                     for index, wire in enumerate(wires))
        for image, factor in _local_images(name)[local]:
            output = label & ~clear
            for index, wire in enumerate(wires):
                output |= ((image >> index) & 1) << wire
                output |= ((image >> (local_width + index)) & 1) << (width + wire)
            result[output] = result.get(output, CZ) + coefficient * factor
    return {label: coefficient for label, coefficient in result.items()
            if coefficient != CZ}


def _affine_rank(support):
    origin = next(iter(support))
    basis = {}
    for label in support:
        vector = label ^ origin
        while vector:
            pivot = vector.bit_length() - 1
            if pivot not in basis:
                basis[pivot] = vector
                break
            vector ^= basis[pivot]
    return len(basis)


def _operator(poly, width):
    result = scale(eye(1 << width), CZ)
    for label, coefficient in poly.items():
        result = add(result, scale(_tau(label, width), coefficient))
    return result


class NativeIndicatorSharingTests(unittest.TestCase):
    def test_exact_native_support_with_address_motion_and_arbitrary_dirty_wire(self):
        # Address 0,1; bank probe 2; arbitrary dirty wire 3. The word is
        # intentionally nonmonomial and moves the address internally.
        width, probe = 4, 1 << 2
        word = [('T', 3), ('H', 0), ('CX', 2, 0), ('T', 0),
                ('CX', 0, 3), ('TDG', 3), ('H', 3), ('H', 1),
                ('CX', 3, 1), ('T', 1), ('CX', 1, 2), ('S', 2)]
        poly, prefix, active = {probe: CO}, eye(16), 0
        ranks, inactive_t = [], 0
        for gate in word:
            before = poly
            poly = _propagate(poly, width, gate)
            if gate[0] in ('T', 'TDG'):
                z = _tau(1 << (width + gate[1]), width)
                actual_before = _operator(before, width)
                commutes = multiply(z, actual_before) == multiply(actual_before, z)
                self.assertEqual(commutes, poly == before)
                active += not commutes
                inactive_t += commutes
            prefix = multiply(native_columns(width, [gate]), prefix)
            self.assertEqual(_operator(poly, width),
                             _conjugate(prefix, _tau(probe, width)))
            ranks.append(_affine_rank(poly))
            self.assertLessEqual(ranks[-1], active)
        self.assertGreater(max(ranks), 1)
        self.assertGreater(inactive_t, 0)
        self.assertTrue(any(sum(entry != CZ for entry in row) > 1 for row in prefix))
        actual = native_columns(width, word)
        self.assertEqual(prefix, actual)
        self.assertEqual(native_columns(width, _inverse(word)), adjoint(actual))
        self.assertEqual(native_columns(width, word + _inverse(word)), eye(16))

    def test_literal_signed_translation_gap_and_one_bit_exception(self):
        # B_u=X_Y(1-2[x=u]). Since the common X_Y is unitary, the
        # signed-translation gap is the maximum diagonal difference.
        # Translation permutations below come from actual native words.
        for p in range(1, 5):
            m = 1 << p
            for displacement in range(1, m):
                word = [gate for wire in range(p) if displacement >> wire & 1
                        for gate in _x(wire)]
                translation = native_columns(p, word)
                for column in range(m):
                    self.assertEqual(translation[column ^ displacement][column], CO)
                    self.assertEqual(sum(entry != CZ for entry in
                                         (row[column] for row in translation)), 1)
                for u in range(m):
                    diagonal = [1 - 2 * (x == u) for x in range(m)]
                    shifted = [sum(translation[row][column] * diagonal[row]
                                   for row in range(m)) for column in range(m)]
                    for sign in (-1, 1):
                        differences = [shifted[x] - sign * diagonal[x] for x in range(m)]
                        if p == 1 and sign == -1:
                            self.assertEqual(set(differences), {CZ})
                        else:
                            self.assertTrue(all(value in (CZ, 2 * CO, -2 * CO)
                                                for value in differences))
                            self.assertTrue(any(value in (2 * CO, -2 * CO)
                                                for value in differences))

    def test_joint_seven_t_bank_shares_native_incidence_below_row_bound(self):
        # This retained p=2 native ANF word jointly mixes bank parities.
        # Its seven T gates are fewer than mp=8, so a universal mp bound
        # would already be false. The incidence requirement survives.
        p, m, width = 2, 4, 6
        word = _phase_query_word()
        self.assertEqual(sum(gate[0] in ('T', 'TDG') for gate in word), 7)
        self.assertLess(7, m * p)
        actual = native_columns(width, word)
        for column in range(1 << width):
            address, bank = column & (m - 1), column >> p
            expected = -CO if bank >> address & 1 else CO
            for row in range(1 << width):
                self.assertEqual(actual[row][column], expected if row == column else CZ)
        self.assertEqual(native_columns(width, _inverse(word)), adjoint(actual))

        polynomials = [{1 << (p + u): CO} for u in range(m)]
        active_counts, incidences = [0] * m, []
        for gate in word:
            incidence = 0
            for u in range(m):
                before = polynomials[u]
                polynomials[u] = _propagate(before, width, gate)
                if gate[0] in ('T', 'TDG') and polynomials[u] != before:
                    active_counts[u] += 1
                    incidence += 1
                self.assertLessEqual(_affine_rank(polynomials[u]), active_counts[u])
            if gate[0] in ('T', 'TDG'):
                incidences.append(incidence)
        self.assertTrue(all(count >= p for count in active_counts))
        self.assertGreaterEqual(sum(incidences), m * p)
        self.assertEqual(sum(incidences), sum(active_counts))
        self.assertGreater(max(incidences), 1)
        for u in range(m):
            self.assertEqual(_operator(polynomials[u], width),
                             _conjugate(actual, _tau(1 << (p + u), width)))


if __name__ == '__main__':
    unittest.main()
