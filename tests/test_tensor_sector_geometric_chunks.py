"""Bounded checks for tensor-sector promotion of geometric operator chunks.

Menu normalization is checked in exact rational pairs. Small dense operator
fixtures check literal phases, all payload/helper columns, occupied keys and
a noncommuting pair of seed bases. The proof handles the all-size statements.
"""
from fractions import Fraction
import math
import unittest

import numpy as np


IDENTITY = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.diag([1, -1]).astype(complex)
H = (X + Z) / math.sqrt(2)


def _kron(factors):
    matrix = np.ones((1, 1), dtype=complex)
    for factor in factors:
        matrix = np.kron(matrix, factor)
    return matrix


def _majoranas(wires):
    return [
        _kron([Z] * wire + [axis] + [IDENTITY] * (wires - wire - 1))
        for wire in range(wires)
        for axis in (X, Y)
    ]


def _geometric_powers(length):
    return [min(position + 1, length - 1) for position in range(length)]


def _exact_clifford_menu(length):
    """Atoms are a sign and one/two Majorana axes, or signed identity."""
    powers = _geometric_powers(length)
    rational = {
        index: Fraction(1, 2 ** (power // 2))
        for index, power in enumerate(powers)
        if power % 2 == 0
    }
    irrational = {
        index: Fraction(1, 2 ** (power // 2))
        for index, power in enumerate(powers)
        if power % 2 == 1
    }
    atoms = [(weight / 2, (index,), 1) for index, weight in rational.items()]
    partners = {index: weight for index, weight in irrational.items() if index}
    atoms.extend((weight / 2, (0, index), 1) for index, weight in partners.items())
    leftover = irrational[0] - sum(partners.values())
    atoms.extend((leftover / 4, (0, 1), sign) for sign in (-1, 1))
    padding = (1 - sum(weight for weight, _, _ in atoms)) / 2
    atoms.extend((padding, (), sign) for sign in (-1, 1))
    return atoms


def _weighted_geometric(axes):
    return sum(
        (2 ** (-power / 2)) * axis
        for power, axis in zip(_geometric_powers(len(axes)), axes)
    )


def _conjugate(word, operator):
    return word @ operator @ word.conj().T


def _exact_priority(local_keys):
    """Return the unique first passing subgroup's selected slot."""
    for subgroup, keys in enumerate(local_keys):
        predicted = 0
        for slot, key in enumerate(keys):
            if key == 0:
                predicted ^= slot
        if keys[predicted] == 0:
            return subgroup, predicted
    return None


def _local_partial_rotation(bad, angle, axis):
    whole = np.eye(len(bad), dtype=complex)
    rotation = math.cos(angle) * IDENTITY + 1j * math.sin(angle) * axis
    return np.kron(bad, IDENTITY) + np.kron(whole - bad, rotation)


class TensorSectorGeometricChunkTests(unittest.TestCase):
    def test_dyadic_clifford_menu_has_literal_coefficients_and_unit_weight(self):
        for length in range(2, 23):
            with self.subTest(length=length):
                atoms = _exact_clifford_menu(length)
                self.assertEqual(sum(weight for weight, _, _ in atoms), 1)
                self.assertLessEqual(len(atoms), length + 4)
                rational = [Fraction(0)] * length
                irrational = [Fraction(0)] * length
                scalar = Fraction(0)
                for weight, axes, sign in atoms:
                    self.assertGreaterEqual(weight, 0)
                    denominator = weight.denominator
                    self.assertEqual(denominator & (denominator - 1), 0)
                    self.assertLessEqual(denominator.bit_length(), length + 4)
                    if len(axes) == 0:
                        scalar += sign * weight
                    elif len(axes) == 1:
                        rational[axes[0]] += sign * weight
                    else:
                        irrational[axes[0]] += weight
                        irrational[axes[1]] += sign * weight
                self.assertEqual(scalar, 0)
                for index, power in enumerate(_geometric_powers(length)):
                    coefficient = Fraction(1, 2 ** (power // 2 + 1))
                    if power % 2:
                        self.assertEqual(rational[index], 0)
                        self.assertEqual(irrational[index], coefficient)
                    else:
                        self.assertEqual(rational[index], coefficient)
                        self.assertEqual(irrational[index], 0)

    def test_majorana_chunk_preserves_outside_axes_and_composes(self):
        axes = _majoranas(4)
        original = axes[:6]
        spare, auxiliary = axes[6:]
        dimension = len(spare)
        comparison = np.eye(dimension, dtype=complex)
        for first, last in ((0, 2), (2, 3), (3, 5)):
            active = original[first:last + 1]
            chunk = _weighted_geometric(active)
            gp = _weighted_geometric([spare] + active)
            gpp = _weighted_geometric([auxiliary, spare] + active)
            np.testing.assert_allclose(gp, (spare + chunk) / math.sqrt(2), atol=2e-15)
            np.testing.assert_allclose(_conjugate(gpp, auxiliary), gp, atol=2e-15)
            plane = (np.eye(dimension) + spare @ active[0]) / math.sqrt(2)
            word = gp @ spare @ plane
            np.testing.assert_allclose(_conjugate(word, active[0]), chunk, atol=2e-15)
            for position, axis in enumerate(original):
                if not first <= position <= last:
                    np.testing.assert_allclose(_conjugate(word, axis), axis, atol=2e-15)
            comparison = word @ comparison
        target = _weighted_geometric(original)
        image = _conjugate(comparison, original[0])
        np.testing.assert_allclose(image, target, atol=3e-15)
        # The scalar of the promoted unitary cancels in the loader conjugation.
        np.testing.assert_allclose(
            _conjugate(1j * comparison, original[0]), target, atol=3e-15
        )

    def test_guarded_center_certifies_all_columns_with_occupied_other_slots(self):
        # Qubit order from low to high is main, slot-0, key-0, slot-1, key-1.
        # On a zero key, each slot word acts by H. On a nonzero key it acts by Y.
        slot_word = np.block(
            [[H, np.zeros_like(H)], [np.zeros_like(Y), Y]]
        )
        count = 32
        v_all = np.kron(np.kron(slot_word, slot_word), IDENTITY)
        selection = np.zeros((count, count), dtype=complex)
        central = np.eye(count, dtype=complex)
        expected = np.eye(count, dtype=complex)
        angle = math.pi / 4
        for column in range(count):
            keys = ((column >> 2) & 1, (column >> 4) & 1)
            slot = next((index for index, key in enumerate(keys) if key == 0), None)
            if slot is None:
                selection[column, column] = 1
                continue
            bit = (1, 3)[slot]
            different = (column ^ (column >> bit)) & 1
            swapped_column = column ^ (different << bit) ^ different
            selection[swapped_column, column] = 1
            eigenvalue = 1 - 2 * ((column >> bit) & 1)
            central[column, column] = np.exp(1j * angle * eigenvalue)
            expected[column, column] = math.cos(angle)
            expected[column ^ 1, column] = 1j * math.sin(angle)
        actual = selection.conj().T @ v_all.conj().T @ central @ v_all @ selection
        np.testing.assert_allclose(actual, expected, atol=1e-15)
        # Include raw subgroup collisions and priority across distinct groups.
        self.assertEqual(_exact_priority(((0, 0), (1, 0))), (0, 1))
        self.assertEqual(_exact_priority(((1, 1), (1, 0))), (1, 1))
        self.assertIsNone(_exact_priority(((1, 0, 0, 1), (1, 1, 1, 1))))

    def test_tensor_bad_cover_noncommuting_seed_bases(self):
        # One fixed key block's bad sector occupies one Walsh basis vector.
        bad_z = np.diag([1, 0, 0, 0]).astype(complex)
        walsh = np.kron(H, H)
        bad_x = _conjugate(walsh, bad_z)
        for blocks in (1, 2, 3):
            with self.subTest(blocks=blocks):
                bad_first = _kron([bad_z] * blocks)
                bad_second = _kron([bad_x] * blocks)
                dimension = len(bad_first)
                overlap = math.sqrt(np.trace(bad_first @ bad_second).real)
                # Both projectors are rank one, so this is the exact norm.
                self.assertAlmostEqual(overlap, 2 ** (-blocks), places=14)
                a_first = _local_partial_rotation(bad_first, math.pi / 2, X)
                a_second = _local_partial_rotation(bad_second, math.pi / 4, X)
                mirror = np.kron(bad_first, IDENTITY) + np.kron(
                    np.eye(dimension) - bad_first, Z
                )
                word = a_first @ a_second @ mirror @ a_second @ mirror.conj().T
                target = np.kron(np.eye(dimension), 1j * X)
                # The difference vanishes outside the span of the two bad
                # lines. Compress to that seed span tensor the payload.
                first_line = np.eye(dimension)[:, 0]
                second_line = np.ones(dimension) / math.sqrt(dimension)
                perpendicular = second_line - overlap * first_line
                perpendicular /= np.linalg.norm(perpendicular)
                support = np.kron(np.column_stack((first_line, perpendicular)), IDENTITY)
                delta = word - target
                compressed = support.conj().T @ delta @ support
                restored = support @ compressed @ support.conj().T
                np.testing.assert_allclose(delta, restored, atol=1e-14)
                self.assertLessEqual(np.linalg.norm(compressed, ord=2), 4 * overlap + 1e-14)
                np.testing.assert_allclose(word.conj().T + target, delta.conj().T, atol=1e-14)

    def test_capped_reference_ray_returns_the_arbitrary_dirty_port(self):
        for levels in (1, 2, 3):
            dimension = 2 ** levels
            amplitudes = 2.0 ** (-np.arange(dimension) / 2)
            source = amplitudes / np.linalg.norm(amplitudes)
            ideal = amplitudes.copy()
            ideal[-1] *= math.sqrt(2)
            ideal /= np.linalg.norm(ideal)
            p = source[-1] ** 2
            q = ideal[-1] ** 2
            c, d = 1 - 2 * p, 2 * (1 - p)
            cosine = (c * c + d * d - q / p) / (2 * c * d)
            phase = math.acos(cosine)
            reflection = np.eye(dimension) - 2 * np.outer(source, source)
            prepared = source.astype(complex)
            prepared[-1] *= np.exp(1j * phase)
            output = reflection @ prepared
            phase_global = -np.angle(output[0])
            phase_relative = np.angle(output[0]) - np.angle(output[-1])
            # The two Z sectors of the same arbitrary port carry conjugate operators.
            for sign in (-1, 1):
                sector = source.astype(complex)
                sector[-1] *= np.exp(1j * sign * phase)
                sector = reflection @ sector
                sector[-1] *= np.exp(1j * sign * phase_relative)
                sector *= np.exp(1j * sign * phase_global)
                np.testing.assert_allclose(sector, ideal, atol=2e-15)


if __name__ == "__main__":
    unittest.main()
