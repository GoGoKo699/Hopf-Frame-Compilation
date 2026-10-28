"""Exact fixtures for the complete field decoder's Pauli-domain restriction.

The analytic proof, rather than the bounded checks, supplies the all-size
count bound. No fixture synthesizes a uniform decoder or implies a depth bound.
"""
from itertools import combinations, product
import unittest

from tests.exact_matrices import Q2
from tests.test_flat_field_bank_native import AXES, decoder_word
from tests.test_one_band_copy_feasibility import _conjugate


def _rank(values):
    pivots = {}
    for value in values:
        while value:
            bit = value.bit_length() - 1
            if bit not in pivots:
                pivots[bit] = value
                break
            value ^= pivots[bit]
    return len(pivots)


def _span(values):
    result = {0}
    for value in values:
        result |= {old ^ value for old in tuple(result)}
    return result


def _pairing(left, right, width):
    mask = (1 << width) - 1
    return (((left & mask) & (right >> width)).bit_count()
            + ((right & mask) & (left >> width)).bit_count()) % 2


def _hermitian_phase(label, width):
    return ((label & ((1 << width) - 1)) & (label >> width)).bit_count() % 4


def _product_phase(left, right, width):
    return (_hermitian_phase(left, width) + _hermitian_phase(right, width)
            - _hermitian_phase(left ^ right, width)
            + 2 * ((left >> width)
                   & (right & ((1 << width) - 1))).bit_count()) % 4


def _field_multiply(left, right, degree, polynomial):
    result = 0
    while right:
        if right & 1:
            result ^= left
        right >>= 1
        left <<= 1
        if (left >> degree) & 1:
            left ^= polynomial
    return result


def _field_image(label, address, degree, polynomial, helpers=0):
    width = 2 * degree + helpers
    mask = (1 << width) - 1
    low = (1 << degree) - 1
    a, b = label & mask, label >> width
    columns = [_field_multiply(address, 1 << j, degree, polynomial)
               for j in range(degree)]
    change_a = _field_multiply(address, a & low, degree, polynomial)
    change_b = sum(((col & ((b >> degree) & low)).bit_count() % 2) << j
                   for j, col in enumerate(columns))
    return (a ^ (change_a << degree)) | ((b ^ change_b) << width)


def _translate(mask, address, size):
    return sum(((mask >> j) & 1) << (j ^ address) for j in range(size))


def _target_image(label, address, p=3, degree=3, polynomial=0b1011,
                  helpers=0):
    size = 1 << p
    residual_width = 2 * degree + helpers
    width = size + residual_width
    mask = (1 << width) - 1
    bank_mask = (1 << size) - 1
    a, b = label & mask, label >> width
    helper_label = (a >> size) | ((b >> size) << residual_width)
    helper_image = _field_image(helper_label, address, degree, polynomial,
                                helpers)
    helper_mask = (1 << residual_width) - 1
    aa = (_translate(a & bank_mask, address, size)
          | ((helper_image & helper_mask) << size))
    bb = (_translate(b & bank_mask, address, size)
          | ((helper_image >> residual_width) << size))
    return aa | (bb << width)


class FieldDecoderPauliDomainTests(unittest.TestCase):
    def test_source_fixed_kernels_and_literal_derivatives(self):
        # Original p=3 field and a character-compatible enlarged source.
        # Basis and pair probes include every term in the bilinear phases.
        for degree, polynomial in ((4, 0b10011), (7, 0b10000011)):
            for helpers in (0, 1):
                width = 2 * degree + helpers
                basis = [1 << j for j in range(2 * width)]
                probes = ({0, (1 << (2 * width)) - 1} | set(basis)
                          | {a ^ b for a, b in combinations(basis, 2)})
                fixed = ([1 << (width + j) for j in range(degree)]
                         + [1 << (degree + j) for j in range(degree)]
                         + [1 << j for j in range(2 * degree, width)]
                         + [1 << (width + j)
                            for j in range(2 * degree, width)])
                self.assertEqual(_rank(fixed), 2 * degree + 2 * helpers)
                for address in range(1, 8):
                    image = lambda v: _field_image(
                        v, address, degree, polynomial, helpers)
                    self.assertEqual(_rank([image(v) ^ v for v in basis]),
                                     2 * degree)
                    self.assertTrue(all(image(v) == v for v in fixed))
                    for label in probes:
                        target = image(label)
                        self.assertEqual(_pairing(label, target, width), 0)
                        self.assertEqual(image(target), label)
                        sign = (_hermitian_phase(label, width)
                                - _hermitian_phase(target, width)) % 4
                        self.assertIn(sign, (0, 2))
                        derivative = label ^ target
                        phase = (sign + _product_phase(target, label, width)) % 4
                        self.assertIn(phase, (0, 2))
                        for other in range(8):
                            self.assertEqual(_field_image(
                                derivative, other, degree, polynomial, helpers),
                                derivative)
                        if target == label:
                            self.assertEqual(sign, 0)

    def test_target_second_difference_space_and_literal_commutation(self):
        p, degree, helpers = 3, 3, 1
        size, residual_width = 1 << p, 2 * degree + helpers
        width = size + residual_width
        image = lambda v, x: _target_image(v, x, p, degree, 0b1011, helpers)
        basis = [1 << j for j in range(2 * width)]
        first_columns = []
        second_columns = []
        for label in basis:
            first, second = 0, 0
            for x in range(8):
                delta = image(label, x) ^ label
                first |= delta << (2 * width * x)
                for y in range(8):
                    delta2 = image(delta, y) ^ delta
                    second |= delta2 << (2 * width * (8 * x + y))
            first_columns.append(first)
            second_columns.append(second)
        self.assertEqual(2 * width - _rank(first_columns),
                         2 + 2 * degree + 2 * helpers)
        self.assertEqual(2 * width - _rank(second_columns),
                         2 * (p + 1) + 4 * degree + 2 * helpers)
        probes = basis + [a ^ b for a, b in combinations(basis, 2)]
        for label, x in product(probes, range(8)):
            target = image(label, x)
            self.assertEqual(_pairing(label, target, width), 0)
            if target == label:
                # Both factors are literal linear basis permutations.
                self.assertEqual((_hermitian_phase(label, width)
                                  - _hermitian_phase(target, width)) % 4, 0)

    def test_helper_quotient_kernel_has_a_nonzero_bank_stabilizer(self):
        p, size, degree, helpers = 3, 8, 3, 1
        width = size + 2 * degree + helpers
        affine = [mask for mask in range(1 << size)
                  if all((_translate(mask, x, size) ^ mask) in (0, 255)
                         for x in range(8))]
        self.assertEqual(len(affine), 16)
        helper_fixed = ([1 << (width + size + j) for j in range(degree)]
                        + [1 << (size + degree + j) for j in range(degree)]
                        + [1 << (width - 1), 1 << (2 * width - 1)])
        self.assertEqual(_rank(helper_fixed), 2 * degree + 2 * helpers)
        for a, b in product(affine, repeat=2):
            bank_label = a | (b << width)
            for helper in [0] + helper_fixed:
                label = bank_label ^ helper
                stabilizers = [x for x in range(1, 8)
                               if _target_image(label, x, helpers=helpers)
                               == label]
                self.assertTrue(stabilizers)
        # The p=2 boundary has two independent mask slopes and no nonzero
        # common translating stabilizer; the finite decoder is not excluded.
        self.assertFalse([x for x in range(1, 4)
                          if _translate(0b1010, x, 4) == 0b1010
                          and _translate(0b1100, x, 4) == 0b1100])

    def test_moving_helper_requires_the_quotient_argument(self):
        # This helper Pauli has query-fixed first derivatives but no nonzero
        # stabilizing query. It falsifies the tempting individual-fixedness
        # inference for the enlarged target, not the quotient theorem.
        label = 1 << 8  # X on the residual field's first u wire.
        for address in range(1, 8):
            target = _target_image(label, address)
            self.assertNotEqual(target, label)
            derivative = target ^ label
            for other in range(8):
                self.assertEqual(_target_image(derivative, other), derivative)
        # Arbitrary Clifford helper motion need not have the commuting gap.
        # H maps X to Z; ZX=iY, so a derivative need not be Hermitian.
        self.assertEqual(_pairing(1, 2, 1), 1)
        self.assertEqual(_product_phase(2, 1, 1), 1)

    def test_finite_decoder_domain_and_outer_clifford_invariants(self):
        word = decoder_word()
        domain = {}
        for label in range(256):
            expansion = _conjugate({label: Q2(1)}, word, 4)
            if len(expansion) == 1:
                target, coefficient = next(iter(expansion.items()))
                self.assertIn(coefficient, (Q2(1), Q2(-1)))
                domain[label] = (target, coefficient)
        expected = _span((1, 4, 8, 18))
        self.assertEqual(set(domain), expected)
        self.assertTrue(all(domain[v] == (v, Q2(1)) for v in expected))
        self.assertEqual(_rank(expected), 4)
        gram = [sum(_pairing(v, w, 4) << j
                    for j, w in enumerate((1, 4, 8, 18)))
                for v in (1, 4, 8, 18)]
        self.assertEqual(_rank(gram), 2)
        axes = [axis for band in AXES for axis in band]
        self.assertEqual(_rank(axes), 4)
        self.assertEqual({v for v in range(256)
                          if all(_pairing(v, a, 4) == 0 for a in axes)},
                         expected)
        # Tensor cells have an extensive exact domain even before any
        # optional fixed Clifford wall is chosen. Clifford walls preserve
        # its dimension, including walls mixing the identity helpers.
        for cells, helpers in ((2, 0), (2, 3), (4, 5)):
            size = 4 * cells
            width = size + helpers
            generators = []
            for cell in range(cells):
                for value in (1, 4, 8, 18):
                    a, b = value & 15, value >> 4
                    generators.append((a << (4 * cell))
                                      | (b << (width + 4 * cell)))
            generators += [1 << j for j in range(size, width)]
            generators += [1 << (width + j) for j in range(size, width)]
            self.assertEqual(_rank(generators), 4 * cells + 2 * helpers)
            self.assertGreater(_rank(generators), 2 + 2 * helpers)


if __name__ == '__main__':
    unittest.main()
