"""Exact checks of a specified regrouping of coefficient-shear half words.

The completed coordinate maps remain valid sources. Their physical halves
cannot be regrouped by the tested palindrome: a coherent two-coordinate
face has maximal probe-image error. This is not a native-depth lower bound.
"""
import unittest

from tests.exact_matrices import CI, CO
from tests.test_common_probe_source import _native
from tests.test_joint_lookup_pauli_source import _bands, _inverse


def _half_words(p, coordinate):
    """Literal F_i and R_i halves, with their own Clifford walls.

    Addresses occupy 0..p-1 and coefficient A_S occupies p+S.
    F_i=exp(i*pi*sum_S X_{A_S} n_{A_{S+i}}/4).
    R_i replaces the angle by -Z_{x_i} times that angle.
    Each word has one physical T band. The complete shear is R_i F_i.
    """
    size = 1 << p
    left = [p + subset for subset in range(size)
            if not (subset >> coordinate) & 1]
    right = [wire + (1 << coordinate) for wire in left]
    hadamards = [('H', wire) for wire in left]
    pairs = [('CX', a, b) for a, b in zip(left, right)]
    controls = [('CX', coordinate, wire) for wire in left]
    first_band = [gate for a, b in zip(left, right)
                  for gate in [('TDG', a), ('T', b)]]
    second_band = [gate for a, b in zip(left, right)
                   for gate in [('T', a), ('TDG', b)]]
    first = (hadamards + pairs + first_band
             + list(reversed(pairs)) + hadamards)
    conditional = (hadamards + controls + pairs + second_band
                   + list(reversed(pairs)) + list(reversed(controls))
                   + hadamards)
    return first, conditional


def _completed_word(p, coordinate):
    first, conditional = _half_words(p, coordinate)
    return first + conditional


def _palindrome(p):
    halves = [_half_words(p, coordinate) for coordinate in range(p)]
    return ([gate for first, _ in halves for gate in first]
            + [gate for _, conditional in reversed(halves)
               for gate in conditional])


def _completed_column(p, coordinate, basis):
    """Complete monomial action, including the literal relative i phase."""
    active = (basis >> coordinate) & 1
    bank = basis >> p
    moved, exponent = bank, 0
    for subset in range(1 << p):
        if (subset >> coordinate) & 1:
            continue
        upper = (bank >> (subset | (1 << coordinate))) & 1
        exponent += active * upper
        moved ^= (active * upper) << subset
    phase = (CO, CI, -CO, -CI)[exponent % 4]
    return (basis & ((1 << p) - 1)) | (moved << p), phase


def _symplectic_columns(width, word):
    """Binary Pauli action of literal H/S/CX gates; signs are irrelevant here."""
    columns = []
    for basis in range(2 * width):
        x = (1 << basis) if basis < width else 0
        z = (1 << (basis - width)) if basis >= width else 0
        for gate, *wires in word:
            if gate == 'H':
                bit = wires[0]
                if ((x ^ z) >> bit) & 1:
                    x ^= 1 << bit
                    z ^= 1 << bit
            elif gate == 'S':
                z ^= ((x >> wires[0]) & 1) << wires[0]
            elif gate == 'CX':
                control, target = wires
                x ^= ((x >> control) & 1) << target
                z ^= ((z >> target) & 1) << control
            else:
                raise ValueError(gate)
        columns.append(x | (z << width))
    return tuple(columns)


def _binary_product(left, right):
    columns = []
    for column in right:
        value = 0
        while column:
            low = column & -column
            value ^= left[low.bit_length() - 1]
            column ^= low
        columns.append(value)
    return tuple(columns)


def _difference_from_identity(columns):
    return tuple(column ^ (1 << index) for index, column in enumerate(columns))


def _paid_half_crossing(p, i, j):
    """K_ij=R_i† L_j† R_i L_j, with four native bands and no helpers."""
    li, _ = _half_words(p, i)
    lj, _ = _half_words(p, j)
    bottoms = {p + subset for subset in range(1 << p)
               if not (subset & ((1 << i) | (1 << j)))}
    # Chronological L_j,L_i,L_j†,L_i†. The two positive and two
    # negative rotations on each bottom's X axis are central and cancel,
    # including the literal scalar factors of the physical T gates.
    raw = lj + li + _inverse(lj) + _inverse(li)
    core = [gate for gate in raw
            if not (gate[0] in ('T', 'TDG') and gate[1] in bottoms)]
    x = [('H', i), ('S', i), ('S', i), ('H', i)]
    controlled = []
    for bottom in sorted(bottoms):
        b, c = bottom + (1 << i), bottom + (1 << j)
        controlled += [('CX', i, b), ('H', c), ('CX', i, c), ('H', c)]
    negative_control = x + controlled + x
    return negative_control + core + negative_control, core, raw


class NativeShearFusionTests(unittest.TestCase):
    def test_completed_coordinate_all_columns_phases_and_actual_adjoint(self):
        p = 2
        for coordinate in range(p):
            word = _completed_word(p, coordinate)
            for basis in range(1 << (p + (1 << p))):
                row, phase = _completed_column(p, coordinate, basis)
                actual = _native(word, {basis: CO})
                self.assertEqual(actual, {row: phase})
                self.assertEqual(_native(_inverse(word), actual), {basis: CO})
        # Equal binary permutations do not permit replacing the actual inverse.
        basis = 1 | (8 << p)
        word = _completed_word(p, 0)
        self.assertEqual(_native(word + word, {basis: CO}), {basis: -CO})

    def test_completed_coordinate_exchange_retains_its_quadratic_phase(self):
        p = 2
        first, second = _completed_word(p, 0), _completed_word(p, 1)
        for basis in range(1 << (p + (1 << p))):
            x0, x1 = basis & 1, (basis >> 1) & 1
            bank = basis >> p
            b, c, d = (bank >> 1) & 1, (bank >> 2) & 1, (bank >> 3) & 1
            sign = -1 if x0 * x1 * d * (b ^ c) else 1
            forward = _native(first + second, {basis: CO})
            reversed_order = _native(second + first, {basis: CO})
            self.assertEqual(forward,
                             {row: sign * value for row, value in reversed_order.items()})

    def test_palindrome_zero_address_cancels_with_literal_phases(self):
        p = 2
        word = _palindrome(p)
        for bank in range(1 << (1 << p)):
            basis = bank << p
            self.assertEqual(_native(word, {basis: CO}), {basis: CO})

    def test_one_active_address_has_maximal_coherent_probe_error(self):
        p = 2
        word = _palindrome(p)
        echo = _inverse(word) + [('S', p), ('S', p)] + word
        # x=(0,1), d=A_{0,1}=1. On this entire block the actual probe
        # image is -Y_a Z_b Y_c; the required image is Z_a Z_c.
        for bank in range(8, 16):
            basis = 2 | (bank << p)
            a, b, c = bank & 1, (bank >> 1) & 1, (bank >> 2) & 1
            row = basis ^ (1 << p) ^ (1 << (p + 2))
            sign = -1 if a ^ b ^ c else 1
            self.assertEqual(_native(echo, {basis: CO}), {row: sign * CO})
        # Omit the common normalization sqrt(2): Bell-minus on a,c,
        # with b=0,d=1, is fixed by the target and negated by the source.
        witness = {2 | (8 << p): CO, 2 | (13 << p): -CO}
        target = {basis: (-value if (((basis >> p) ^ (basis >> (p + 2))) & 1)
                          else value) for basis, value in witness.items()}
        self.assertEqual(target, witness)
        self.assertEqual(_native(echo, witness),
                         {basis: -value for basis, value in witness.items()})
        self.assertEqual(_native(word + _inverse(word), witness), witness)
        # After normalization, actual-minus-target is -2 times this vector.

    def test_face_embedding_and_independent_native_half_bands(self):
        for p in range(2, 9):
            for i, j in {(0, 1), (0, p - 1), (p - 2, p - 1)}:
                face = {0, 1 << i, 1 << j, (1 << i) | (1 << j)}
                for coordinate in range(p):
                    pairs = [(subset, subset | (1 << coordinate))
                             for subset in range(1 << p)
                             if not (subset >> coordinate) & 1]
                    surviving = [(a, b) for a, b in pairs if b in face]
                    if coordinate not in (i, j):
                        self.assertEqual(surviving, [])
                    else:
                        self.assertEqual(len(surviving), 2)
                        self.assertTrue(all(a in face for a, _ in surviving))
                    # All other rotations have a zero upper control on the
                    # whole face subspace. Surviving targets stay in the face;
                    # the top coefficient remains fixed. No large matrices.
                    self.assertTrue(all(a != ((1 << i) | (1 << j))
                                        for a, _ in surviving))
            for coordinate in range(p):
                for half in _half_words(p, coordinate):
                    bands = _bands(half)
                    self.assertEqual(len(bands), 1)
                    self.assertEqual(len(bands[0]), 1 << p)
                    self.assertEqual(len({gate[1] for gate in bands[0]}), 1 << p)

    def test_no_common_clifford_diagonalizer_with_identity_helpers(self):
        for p in range(2, 7):
            size = 1 << p
            i, j = 0, p - 1
            words = [[('CX', subset | (1 << coordinate), subset)
                      for subset in range(size)
                      if not (subset >> coordinate) & 1]
                     for coordinate in (i, j)]
            for helpers in (0, 2):
                width = size + helpers
                identity = tuple(1 << bit for bit in range(2 * width))
                zero = (0,) * (2 * width)
                first, second = [_symplectic_columns(width, word) for word in words]
                self.assertEqual(_binary_product(first, second),
                                 _binary_product(second, first))
                ni, nj = map(_difference_from_identity, (first, second))
                self.assertEqual(_binary_product(ni, ni), zero)
                self.assertEqual(_binary_product(nj, nj), zero)
                cross = _binary_product(ni, nj)
                # On X labels this is N_i N_j e_{ij}=e_empty.
                self.assertEqual(cross[(1 << i) | (1 << j)], 1)
                for helper in range(size, width):
                    self.assertEqual(cross[helper], 0)
                    self.assertEqual(cross[width + helper], 0)

                basis = [('H', 0), ('S', 1), ('CX', 0, 1), ('H', 2)]
                if helpers:
                    basis += [('CX', size, 0), ('H', size + 1), ('CX', 2, size + 1)]
                change = _symplectic_columns(width, basis)
                # S and S-dagger have the same binary symplectic action.
                reverse = _symplectic_columns(width, list(reversed(basis)))
                self.assertEqual(_binary_product(change, reverse), identity)
                changed = [_binary_product(change, _binary_product(f, reverse))
                           for f in (first, second)]
                changed_cross = _binary_product(
                    *map(_difference_from_identity, changed))
                self.assertEqual(changed_cross,
                                 _binary_product(change, _binary_product(cross, reverse)))
                self.assertNotEqual(changed_cross, zero)

                # Two explicit diagonal Cliffords built from native S and CZ
                # decompositions have square-zero cross differences instead.
                diagonal_words = [[('S', 0)], [('S', width - 1)]]
                for target in range(1, width):
                    for index, control in ((0, target - 1), (1, 0)):
                        diagonal_words[index] += [('H', target), ('CX', control, target),
                                                  ('H', target)]
                differences = [_difference_from_identity(_symplectic_columns(width, word))
                               for word in diagonal_words]
                self.assertTrue(all(any(difference) for difference in differences))
                self.assertEqual(_binary_product(*differences), zero)
                self.assertEqual(_binary_product(*reversed(differences)), zero)

    def test_paid_half_crossing_all_columns_and_actual_native_inverse(self):
        p, i, j = 2, 0, 1
        _, ri = _half_words(p, i)
        lj, _ = _half_words(p, j)
        target = lj + ri + _inverse(lj) + _inverse(ri)
        paid, core, raw = _paid_half_crossing(p, i, j)
        for basis in range(1 << (p + (1 << p))):
            state = {basis: CO}
            self.assertEqual(_native(core, state), _native(raw, state))
            actual = _native(paid, state)
            self.assertEqual(actual, _native(target, state))
            self.assertEqual(_native(_inverse(paid), actual), state)
            # Literal interchange L_j R_i=R_i L_j K_ij†.
            self.assertEqual(_native(ri + lj, state),
                             _native(_inverse(paid) + lj + ri, state))

    def test_paid_half_crossing_common_forward_inverse_ledgers(self):
        for p in range(2, 7):
            size, cells = 1 << p, 1 << (p - 2)
            paid, core, raw = _paid_half_crossing(p, 0, p - 1)
            self.assertEqual(len(raw) - len(core), 4 * cells)
            self.assertEqual(len(core), 44 * cells)
            for word in (paid, _inverse(paid)):
                t_count = sum(gate[0] in ('T', 'TDG') for gate in word)
                self.assertEqual(t_count, 3 * size)
                self.assertEqual(len(word), 13 * size + 16)
                self.assertEqual(len(word) - t_count, 10 * size + 16)
                bands = _bands(word)
                self.assertEqual(list(map(len, bands)), [3 * cells] * 4)
                self.assertTrue(all(len({gate[1] for gate in band}) == len(band)
                                    for band in bands))
                self.assertTrue(all(0 <= wire < p + size
                                    for gate in word for wire in gate[1:]))


if __name__ == '__main__':
    unittest.main()
