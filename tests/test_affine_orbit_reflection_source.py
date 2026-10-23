"""Exact guards for the affine orbit reflection and its scoped completions.

The full p=1 native columns include every address/bank/private-helper
input.  The p=2 field fixture checks complete Boolean carrier action;
neither the finite repaired pair nor the degree guard proves a general
native-depth upper or lower bound.
"""
import unittest

from tests.exact_matrices import CO, CZ, adjoint, native_columns
from tests.test_corank_fourier_bell_source import _kernel
from tests.test_half_bilinear_gram_source import _binary_image, _private_toffoli
from tests.test_majorana_walsh_query import _inverse
from tests.test_phase_signal_query_limits import _counts


def _x(wire):
    return [('H', wire), ('S', wire), ('S', wire), ('H', wire)]


def _private_ccz(a, b, c, phase):
    # Removing the two target Hadamards leaves the exact diagonal phase.
    return _private_toffoli(a, b, c, phase)[1:-1]


def _private_c3x(a, b, c, target, dirty, phase):
    toggle = _private_toffoli(a, b, dirty, phase)
    sign = _private_ccz(dirty, c, target, phase)
    return ([('H', target)] + toggle + sign + toggle + sign
            + [('H', target)])


def _p1_words():
    # x,z0,z1,t0,t1,d,g.  Both d and g are arbitrary private inputs.
    source = (_private_toffoli(0, 2, 3, 6) + [('CX', 2, 4)]
              + _x(1) + _x(3))
    normalizer = ([('CX', 1, 3)] + _private_toffoli(1, 2, 4, 6)
                  + _private_c3x(1, 0, 2, 3, 5, 6))
    local = [('S', 1), ('S', 1), ('H', 3), ('CX', 1, 3), ('H', 3)]
    repaired = local + _private_ccz(0, 1, 2, 5) + _private_ccz(0, 1, 4, 6)
    return source, normalizer, local, repaired


def _parity(value):
    return value.bit_count() & 1


def _s(rows, z, t):
    return z ^ 1, t ^ _binary_image(rows, z) ^ 1


def _c(rows, z, t):
    return z, t ^ ((_binary_image(rows, z) ^ 1) if z & 1 else 0)


def _l(rows, z, t):
    return z, t ^ (_binary_image(rows, z) if z & 1 else 0)


def _degree(values):
    coefficients = list(values)
    for bit in range((len(values) - 1).bit_length()):
        for mask in range(len(values)):
            if mask >> bit & 1:
                coefficients[mask] ^= coefficients[mask ^ (1 << bit)]
    return max((mask.bit_count() for mask, value in enumerate(coefficients)
                if value), default=0)


def _consistent(rows, rhs, columns):
    pivots = {}
    for row, value in zip(rows, rhs):
        augmented = row | (value << columns)
        while augmented & ((1 << columns) - 1):
            pivot = (augmented & ((1 << columns) - 1)).bit_length() - 1
            if pivot in pivots:
                augmented ^= pivots[pivot]
            else:
                pivots[pivot] = augmented
                break
        else:
            if augmented >> columns:
                return False
    return True


class AffineOrbitReflectionSourceTests(unittest.TestCase):
    def test_native_affine_reflection_and_cubic_normalizer_return_all_helpers(self):
        source, normalizer, _, _ = _p1_words()
        gauge = normalizer[1:]  # Remove only CX(z0,t0), not a hidden inverse.
        self.assertEqual(_counts(source), (8, 19))
        self.assertEqual(_counts(normalizer), (40, 49))
        self.assertEqual(_counts(gauge), (40, 48))
        for word, action in ((source, _s), (normalizer, _c), (gauge, _l)):
            actual = native_columns(7, word)
            self.assertEqual(native_columns(7, _inverse(word)), adjoint(actual))
            for column in range(128):
                x = column & 1
                z, t = column >> 1 & 3, column >> 3 & 3
                rows, _, _ = _kernel(1, x)
                new_z, new_t = action(rows, z, t)
                output = x | (new_z << 1) | (new_t << 3) | (column & 96)
                for row in range(128):
                    self.assertEqual(actual[row][column], CO if row == output else CZ)
        # The conjugacy is not restricted to diagonal or fixed-address R.
        # This Hermitian native R moves the address and has a tilted axis.
        r_word = [('H', 0), ('TDG', 2), ('H', 2), ('T', 2)]
        s_zero = _x(1) + _x(3)
        changed_r = gauge + r_word + gauge
        left = s_zero + changed_r + s_zero + changed_r
        right = gauge + source + r_word + source + r_word + gauge
        self.assertEqual(native_columns(7, left), native_columns(7, right))

    def test_native_completed_p1_pair_and_actual_reverse_on_complete_columns(self):
        source, _, local, repaired = _p1_words()
        self.assertEqual(_counts(local), (0, 5))
        self.assertEqual(_counts(repaired), (16, 21))
        word = source + repaired + source + repaired
        self.assertEqual(_counts(word), (48, 80))
        actual = native_columns(7, word)
        self.assertEqual(native_columns(7, _inverse(word)), adjoint(actual))
        for column in range(128):
            x = column & 1
            z, t = column >> 1 & 3, column >> 3 & 3
            charge = _parity((z ^ t) & (1 | (x << 1)))
            sign = -CO if charge else CO
            for row in range(128):
                self.assertEqual(actual[row][column], sign if row == column else CZ)

    def test_field_orbits_cochains_and_both_all_size_trial_witnesses(self):
        for p in (1, 2):
            for x in range(1 << p):
                rows, vector, n = _kernel(p, x)
                images = [_binary_image(rows, z) for z in range(1 << n)]
                self.assertEqual(images[1], 0)
                self.assertTrue(all(not _parity(vector & image) for image in images))
                for z in range(1 << n):
                    for t in range(1 << n):
                        a = images[z] ^ 1
                        zz, tt = z ^ 1, t ^ a
                        self.assertNotEqual((z, t), (zz, tt))
                        self.assertEqual((zz ^ 1, tt ^ images[zz] ^ 1), (z, t))
                        q = _parity(vector & (z ^ t))
                        self.assertEqual(_parity(vector & (zz ^ tt)), q)
                        # L is a charge-preserving involution and L S0 L=S.
                        lt = t ^ (images[z] if z & 1 else 0)
                        self.assertEqual(lt ^ (images[z] if z & 1 else 0), t)
                        self.assertEqual(_parity(vector & (z ^ lt)), q)
                        lst = lt ^ 1 ^ (images[z ^ 1] if (z ^ 1) & 1 else 0)
                        self.assertEqual((z ^ 1, lst), (zz, tt))
                        u = t ^ (a if z & 1 else 0)
                        self.assertEqual(tt ^ (a if zz & 1 else 0), u)
                        self.assertEqual((u ^ a) ^ (a if zz & 1 else 0), t)
                        self.assertEqual(_parity(vector & (z ^ u)), q ^ (z & 1))
                        # A nonconstant section function tests the free invariant.
                        h = _parity((z & ~1) & 6) ^ ((u & 1) & (u >> 1 & 1))
                        f = ((z & 1) & q) ^ h
                        ff = ((zz & 1) & q) ^ h
                        self.assertEqual(f ^ ff, q)
                        local = (z & 1) ^ ((z & 1) & (t & 1))
                        local_s = (zz & 1) ^ ((zz & 1) & (tt & 1))
                        expected = ((z & 1) ^ (t & 1)
                                    ^ ((1 ^ (z & 1)) & (images[z] & 1)))
                        self.assertEqual(local ^ local_s, expected)
                        if not x:
                            self.assertEqual(expected, q)
                        scalar = (z & 1) & int(u == 0)
                        scalar_s = (zz & 1) & int(u == 0)
                        self.assertEqual(scalar ^ scalar_s, int(u == 0))
                # Scalar-zero trial fails at z=0,t=e for every address.
                witness_z, witness_t = 0, 1
                _, witness_u = _c(rows, witness_z, witness_t)
                self.assertEqual(int(witness_u == 0), 0)
                self.assertEqual(_parity(vector & (witness_z ^ witness_t)), 1)
                if x:
                    tail = 1 << (n - p)
                    self.assertEqual(_parity(vector & tail), 1)
                    self.assertEqual(tail & 1, 0)  # False acceptance.
                    self.assertEqual(_parity(vector & (tail ^ 1)), 0)
                    self.assertEqual((tail ^ 1) & 1, 1)  # False rejection.

    def test_sharp_degree_on_the_affine_restriction_at_three_small_sizes(self):
        # Retain x, lambda=z0, t0, and the first bit of the last t block.
        # On this restriction S flips lambda and t0, and the last v bit
        # is [x != 0].  Test linear inconsistency for every lower-degree
        # phase, not just a few chosen section functions.
        for p in (1, 2, 3):
            variables = p + 3
            direction = (1 << p) | (1 << (p + 1))
            size = 1 << variables
            monomials = [mask for mask in range(size) if mask.bit_count() <= p + 1]
            rows, rhs, canonical = [], [], []
            for value in range(size):
                x = value & ((1 << p) - 1)
                lam, t0, tail = (value >> bit & 1 for bit in range(p, p + 3))
                q = lam ^ t0 ^ (bool(x) & tail)
                image = value ^ direction
                rows.append(sum((int((value & mask) == mask)
                                 ^ int((image & mask) == mask)) << j
                                for j, mask in enumerate(monomials)))
                rhs.append(q)
                canonical.append(lam & q)
            self.assertFalse(_consistent(rows, rhs, len(monomials)))
            self.assertEqual(_degree(canonical), p + 2)
            self.assertTrue(all(canonical[value] ^ canonical[value ^ direction] == rhs[value]
                                for value in range(size)))


if __name__ == '__main__':
    unittest.main()
