"""Exact recursive field decoder and paid multiplication components.

The recoupling permutations below are abstract reversible maps, not native
emitters.  Checking their actual reversed word supplies no native depth or
helper-return budget for those maps.  The separate multiplication fixture
does emit its cubic phase word and tests arbitrary returned dirty masks.
"""
from functools import lru_cache
import unittest

from tests.test_dirty_indicator_depth import _exact_phase_action
from tests.test_field_shear_intertwiner import _field_product, _translation
from tests.test_joint_lookup_pauli_source import _bands, _inverse


@lru_cache(None)
def _tower_product(level, left, right):
    """F_j=F_(j-1)[theta]/(theta^2+theta+delta), delta top basis bit."""
    if level == 1:
        return left & right
    half = 1 << (level - 2)
    mask = (1 << half) - 1
    a, b = left & mask, left >> half
    c, d = right & mask, right >> half
    mul = lambda x, y: _tower_product(level - 1, x, y)
    bd = mul(b, d)
    return (mul(a, c) ^ mul(1 << (half - 1), bd)
            ^ ((mul(a, d) ^ mul(b, c) ^ bd) << half))


@lru_cache(None)
def _tower_inverse(level, value):
    assert value
    # Bounded exact fixture only: this search is not a proposed compiler.
    return next(other for other in range(1, 1 << (1 << (level - 1)))
                if _tower_product(level, value, other) == 1)


def _tower_trace(level, value):
    return (value if level == 1 else
            _tower_trace(level - 1, value >> (1 << (level - 2))))


def _tower_power(level, value, exponent):
    result = 1
    while exponent:
        if exponent & 1:
            result = _tower_product(level, result, value)
        value = _tower_product(level, value, value)
        exponent >>= 1
    return result


def _address_field(level, address):
    return ((address & 1)
            ^ sum(((address >> bit) & 1) << (1 << (bit - 1))
                  for bit in range(1, level)))


def _tower_shear(level, address, basis):
    degree = 1 << (level - 1)
    u = basis & ((1 << degree) - 1)
    return basis ^ (_tower_product(level, _address_field(level, address), u)
                    << degree)


@lru_cache(None)
def _recoupling_permutation(level, kind, inverse=False):
    """Whole D_d or F_d map, including their entire zero-difference sector."""
    half, degree = 1 << (level - 2), 1 << (level - 1)
    mask = (1 << half) - 1
    field_mask = (1 << degree) - 1
    theta, delta = 1 << half, 1 << (half - 1)
    result = []
    for basis in range(1 << (2 * degree)):
        u, v = basis & field_mask, basis >> degree
        if kind == 'D':
            if u:
                d0, d1 = u & mask, u >> half
                norm = (_tower_product(level - 1, d0, d0)
                        ^ _tower_product(level - 1, d0, d1)
                        ^ _tower_product(level - 1, delta,
                                         _tower_product(level - 1, d1, d1)))
                uv = _tower_product(level, u, v)
                phi0 = uv >> half
                phi1 = _tower_product(level, theta, uv) >> half
                s = _tower_product(level - 1,
                                   phi0 ^ _tower_product(level - 1, d0, d1),
                                   _tower_inverse(level - 1, norm))
                v = phi1 | (s << half)
            image = u | (v << degree)
        else:
            assert kind == 'F'
            a = _tower_product(level, u, v) if u else v
            image = a | ((a ^ u) << degree)
        result.append(image)
    if inverse:
        reverse = [None] * len(result)
        for column, row in enumerate(result):
            assert reverse[row] is None
            reverse[row] = column
        return tuple(reverse)
    return tuple(result)


def _recursive_decoder_word(level, offset=0):
    if level == 1:
        return [('CX', offset + 1, offset)]
    half, degree = 1 << (level - 2), 1 << (level - 1)
    return ([('Hf', level, offset + degree), ('D', level, offset, False),
             ('Hf', level - 1, offset + degree), ('F', level, offset, False),
             ('Hf', level - 1, offset + half),
             ('Hf', level - 1, offset + degree + half)]
            + _recursive_decoder_word(level - 1, offset)
            + _recursive_decoder_word(level - 1, offset + degree))


def _reverse_decoder(word):
    return [(kind, *args[:-1], not args[-1]) if kind in ('D', 'F')
            else (kind, *args) for kind, *args in reversed(word)]


def _apply_decoder(word, state, normalization=0):
    """Exact sparse integer numerators over 2**(normalization/2).

    Each field Fourier has one common positive normalization; retaining
    it separately avoids expensive repeated quadratic-field arithmetic.
    This is still literal Q(sqrt(2)) arithmetic, without a phase quotient.
    """
    for kind, *args in word:
        result = {}
        if kind == 'Hf':
            level, offset = args
            degree = 1 << (level - 1)
            mask = ((1 << degree) - 1) << offset
            normalization += degree
            for basis, amplitude in state.items():
                source = (basis & mask) >> offset
                for target in range(1 << degree):
                    sign = (-1) ** _tower_trace(
                        level, _tower_product(level, source, target))
                    row = (basis & ~mask) | (target << offset)
                    result[row] = result.get(row, 0) + amplitude * sign
        elif kind == 'CX':
            control, target = args
            result = {basis ^ (((basis >> control) & 1) << target): value
                      for basis, value in state.items()}
        else:
            level, offset, inverse = args
            width = 1 << level
            mask = ((1 << width) - 1) << offset
            permutation = _recoupling_permutation(level, kind, inverse)
            result = {(basis & ~mask)
                      | (permutation[(basis & mask) >> offset] << offset): value
                      for basis, value in state.items()}
        state = {basis: value for basis, value in result.items() if value}
    return state, normalization


def _polynomial_remainder(value, modulus):
    while value.bit_length() >= modulus.bit_length():
        value ^= modulus << (value.bit_length() - modulus.bit_length())
    return value


def _multiplication_forms(degree, modulus):
    """Evaluation/interpolation into at most (2k-1) ceil(log2(2k-1))^2 terms."""
    count = 2 * degree - 1
    extension_degree = max(1, (count - 1).bit_length())
    extension_modulus = {1: 0b11, 2: 0b111, 3: 0b1011}[extension_degree]
    mul = lambda a, b: _field_product(a, b, extension_modulus)
    forms = []
    for point in range(count):
        powers = [1]
        for _ in range(degree - 1):
            powers.append(mul(powers[-1], point))
        evaluation = [sum(((power >> bit) & 1) << j
                          for j, power in enumerate(powers))
                      for bit in range(extension_degree)]
        polynomial, denominator = [1], 1
        for other in range(count):
            if other == point:
                continue
            updated = [0] * (len(polynomial) + 1)
            for index, coefficient in enumerate(polynomial):
                updated[index] ^= mul(coefficient, other)
                updated[index + 1] ^= coefficient
            polynomial = updated
            denominator = mul(denominator, point ^ other)
        inverse = next(value for value in range(1, 1 << extension_degree)
                       if mul(denominator, value) == 1)
        polynomial = [mul(value, inverse) for value in polynomial]
        for left in range(extension_degree):
            for right in range(extension_degree):
                coefficient = mul(1 << left, 1 << right)
                binary = sum((mul(value, coefficient) & 1) << j
                             for j, value in enumerate(polynomial))
                output = _polynomial_remainder(binary, modulus)
                if evaluation[left] and evaluation[right] and output:
                    forms.append((evaluation[left], evaluation[right], output))
    return forms, count, extension_degree


def _linear_form_cube_word(degree, forms, first_width=None):
    widths = (degree if first_width is None else first_width, degree, degree)
    offsets = (0, widths[0], widths[0] + degree)
    principal_width = sum(widths)
    helpers = list(range(principal_width, principal_width + len(forms)))
    word, old = [], 0
    for step in range(8):
        vertex = step ^ (step >> 1)
        difference = old ^ vertex
        if difference:
            index = difference.bit_length() - 1
            word += [('CX', offsets[index] + bit, helper)
                     for form, helper in zip(forms, helpers)
                     for bit in range(widths[index]) if (form[index] >> bit) & 1]
        word += [('T' if vertex.bit_count() % 2 == 0 else 'TDG', helper)
                 for helper in helpers]
        old = vertex
    word += [('CX', offsets[2] + bit, helper)
             for form, helper in zip(forms, helpers)
             for bit in range(degree) if (form[2] >> bit) & 1]
    return word, helpers


class RecursiveFieldShearDecoderTests(unittest.TestCase):
    def test_tower_trace_and_complete_recoupling_permutations(self):
        for level in (1, 2, 3):
            degree = 1 << (level - 1)
            size = 1 << degree
            self.assertEqual(_tower_trace(level, 1 << (degree - 1)), 1)
            for value in range(1, size):
                self.assertEqual(_tower_product(
                    level, value, _tower_inverse(level, value)), 1)
                self.assertTrue(any(_tower_trace(level, _tower_product(level, value, y))
                                    for y in range(size)))
            self.assertEqual(len({_address_field(level, x)
                                  for x in range(1 << level)}), 1 << level)
            if level >= 2:
                for kind in ('D', 'F'):
                    direct = _recoupling_permutation(level, kind)
                    reverse = _recoupling_permutation(level, kind, True)
                    self.assertEqual(sorted(direct), list(range(size * size)))
                    self.assertTrue(all(reverse[row] == column
                                        for column, row in enumerate(direct)))

    def test_recursive_complete_intertwining_and_actual_abstract_inverse(self):
        for level in (1, 2, 3):
            dimension = 1 << (1 << level)
            word = _recursive_decoder_word(level)
            reverse = _reverse_decoder(word)
            columns = [_apply_decoder(word, {basis: 1}) for basis in range(dimension)]
            for basis, (column, normalization) in enumerate(columns):
                restored, total = _apply_decoder(reverse, column, normalization)
                self.assertEqual(total % 2, 0)
                self.assertEqual(restored, {basis: 1 << (total // 2)})
                for address in range(1 << level):
                    target = {_translation(level, address, row): value
                              for row, value in column.items()}
                    self.assertEqual(columns[_tower_shear(level, address, basis)],
                                     (target, normalization),
                                     (level, address, basis))

    def test_conditional_cell_retains_half_degree_in_its_symplectic_action(self):
        for level in (2, 3, 4):
            half, degree = 1 << (level - 2), 1 << (level - 1)
            mask, delta = (1 << half) - 1, 1 << (half - 1)
            coefficients = []
            for u in range(1 << degree):
                d0, d1 = u & mask, u >> half
                coefficient = 0
                if u:
                    norm = (_tower_product(level - 1, d0, d0)
                            ^ _tower_product(level - 1, d0, d1)
                            ^ _tower_product(level - 1, delta,
                                             _tower_product(level - 1, d1, d1)))
                    coefficient = _tower_product(
                        level - 1, _tower_product(level - 1, delta,
                                                 _tower_product(level - 1, d1, d1)),
                        _tower_inverse(level - 1, norm))
                coefficients.append(coefficient)
                normalized = _tower_product(
                    level - 1, coefficient, _tower_inverse(level - 1, delta))
                q = 1 << half
                self.assertEqual(normalized, _tower_power(level, u, q - 1)
                                 ^ _tower_power(level, u, q * q - q))
                if level <= 3:
                    # Independently extract the XX corner of F_U H_r D_U:
                    # affine offsets cancel, H_r removes the first X half.
                    permutation = _recoupling_permutation(level, 'D')
                    origin = permutation[u] >> degree
                    for bit in range(half):
                        after_d = (permutation[u | (1 << (degree + bit))]
                                   >> degree) ^ origin
                        after_h = after_d & ~mask
                        after_f = (_tower_product(level, u, after_h)
                                   if u else after_h)
                        self.assertEqual(after_f & mask, _tower_product(
                            level - 1, coefficient, 1 << bit))
            anf = coefficients.copy()
            for bit in range(degree):
                for mask_index in range(1 << degree):
                    if mask_index >> bit & 1:
                        anf[mask_index] ^= anf[mask_index ^ (1 << bit)]
            self.assertEqual(max(index.bit_count() for index, value in enumerate(anf)
                                 if value), half)

    def test_near_linear_multiplication_tensor_and_complete_native_masks(self):
        for degree, modulus in ((2, 0b111), (4, 0b10011)):
            forms, count, extension_degree = _multiplication_forms(degree, modulus)
            self.assertLessEqual(len(forms), count * extension_degree ** 2)
            for left in range(1 << degree):
                for right in range(1 << degree):
                    output = 0
                    for a, b, c in forms:
                        if ((left & a).bit_count() & (right & b).bit_count()) & 1:
                            output ^= c
                    self.assertEqual(output, _field_product(left, right, modulus))
            if degree == 2:
                word, helpers = _linear_form_cube_word(degree, forms)
                for actual in (word, _inverse(word)):
                    self.assertEqual(len(_bands(actual)), 8)
                    self.assertEqual(sum(gate[0] in ('T', 'TDG') for gate in actual),
                                     8 * len(forms))
                    self.assertTrue(all({gate[1] for gate in band} == set(helpers)
                                        for band in _bands(actual)))
                    for basis in range(1 << (3 * degree + len(helpers))):
                        principal = basis & ((1 << (3 * degree)) - 1)
                        left = principal & 3
                        right = (principal >> degree) & 3
                        target = principal >> (2 * degree)
                        sign = (_field_product(left, right, modulus) & target).bit_count() % 2
                        self.assertEqual(_exact_phase_action(basis, actual), (basis, 4 * sign))

    def test_dense_tower_query_uses_linear_forms_and_returned_masks(self):
        for level in (2, 3):
            degree = 1 << (level - 1)
            forms = [(1 << address_bit, 1 << input_bit,
                      _tower_product(level, _address_field(level, 1 << address_bit),
                                     1 << input_bit))
                     for address_bit in range(level) for input_bit in range(degree)]
            word, helpers = _linear_form_cube_word(degree, forms, first_width=level)
            self.assertEqual(len(helpers), level * degree)
            self.assertEqual(helpers[0], level + 2 * degree)
            for actual in (word, _inverse(word)):
                self.assertEqual(sum(g[0] in ('T', 'TDG') for g in actual),
                                 8 * level * degree)
                self.assertEqual(len(_bands(actual)), 8)
                self.assertTrue(all(len(band) == len(helpers)
                                    and {gate[1] for gate in band} == set(helpers)
                                    for band in _bands(actual)))
                # Exhaust every physical column at p=2. The p=3 fixture
                # exhausts all principal columns with two opposite masks;
                # its algebraic linear-form identity is also checked below.
                masks = (range(1 << len(helpers)) if level == 2
                         else (0, (1 << len(helpers)) - 1))
                for address in range(1 << level):
                    for u in range(1 << degree):
                        product = _tower_product(level, _address_field(level, address), u)
                        for v in range(1 << degree):
                            parity = (product & v).bit_count() % 2
                            predicted = sum(((address & a).bit_count() % 2)
                                            * ((u & b).bit_count() % 2)
                                            * ((v & c).bit_count() % 2)
                                            for a, b, c in forms) % 2
                            self.assertEqual(predicted, parity)
                            principal = address | (u << level) | (v << (level + degree))
                            for mask in masks:
                                basis = principal | (mask << (level + 2 * degree))
                                self.assertEqual(_exact_phase_action(basis, actual),
                                                 (basis, 4 * parity))


if __name__ == '__main__':
    unittest.main()
