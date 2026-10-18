"""Bounded evidence for the general affine diagonal-source-word theorem.

These checks retain every work column. Arbitrary rational source tables here
are mathematical fixtures, not claims of free native source implementations.
"""

from cmath import exp
from fractions import Fraction
from math import pi, sin, sqrt
import unittest

from tests.test_returned_diagonal_echo_degree import _phase_echo_column


def _linear(columns, value):
    result = 0
    for bit, column in enumerate(columns):
        if value >> bit & 1:
            result ^= column
    return result


def _parity(value):
    return value.bit_count() % 2


def _cube(function, base, directions):
    """Signed additive difference, retaining multiplicity of repeated points."""
    total = 0
    for subset in range(1 << len(directions)):
        point = base
        for bit, direction in enumerate(directions):
            if subset >> bit & 1:
                point ^= direction
        total += (-1) ** (len(directions) - subset.bit_count()) * function(point)
    return total


class DiagonalSourceWordTests(unittest.TestCase):
    def test_joint_cubes_annihilate_distinct_affinely_transported_sources(self):
        # Work width 3, address width 4. All listed work blocks are invertible;
        # arbitrary successive prefix blocks can be connected by affine
        # monomial-Clifford intervals. Denominator 308 retains quarters and
        # odd factors 7 and 11, so this is not Boolean-only phase arithmetic.
        work_mask = 7
        denominator = 308
        blocks = ((1, 2, 4), (3, 2, 4), (2, 4, 1), (1, 3, 5))
        maps = ((1, 2, 3, 4), (0, 6, 1, 1), (7, 0, 4, 2), (2, 2, 0, 3))
        frames = ((1, 2, 4, 8, 3, 5), (3, 5, 6, 0, 3, 5), (1,) * 6)

        def clifford_phase(point):
            # S phases on affine joint parities; CZ includes address-work
            # cross terms. The constant is a permitted fixed global phase.
            return (denominator // 2 + denominator // 4 * _parity((point ^ 1) & 0b1010101)
                    + denominator // 4 * 3 * _parity(point & 0b1100010)
                    + denominator // 2 * _parity(point & 0b0010101)
                    * _parity((point ^ 8) & 0b1110010))

        for calls in range(7):
            order = max(calls, 3)
            sources = []
            for j in range(calls):
                block, address_map = blocks[j % 4], maps[j % 4]
                inverse = {_linear(block, z): z for z in range(8)}
                self.assertEqual(len(inverse), 8)
                table = [((j + 2) * z**3 + (j + 5) * z**2
                          + 17 * (z ^ (j + 1))) % denominator for z in range(8)]
                sources.append((block, address_map, (3 * j + 1) & 7,
                                inverse, table))

            def argument(source, point):
                block, address_map, offset, _, _ = source
                return (_linear(block, point & work_mask)
                        ^ _linear(address_map, point >> 3) ^ offset)

            def phase(point):
                return clifford_phase(point) + sum(source[4][argument(source, point)]
                                                     for source in sources)

            for frame in frames:
                address_directions = frame[:order]
                joint_directions = []
                for j, address_direction in enumerate(address_directions):
                    shift = (sources[j][3][_linear(sources[j][1], address_direction)]
                             if j < calls else 0)
                    joint_directions.append((address_direction << 3) ^ shift)
                for point in range(128):
                    for j, source in enumerate(sources):
                        self.assertEqual(argument(source, point),
                                         argument(source, point ^ joint_directions[j]))
                # The q=4 fixture covers all physical bases. Other q values
                # use a bounded varied sample; the analytic proof is general.
                bases = range(128) if calls == 4 else (0, 1, 7, 8, 19, 57, 94, 127)
                for point in bases:
                    self.assertEqual(_cube(phase, point, joint_directions) % denominator, 0)

    def test_literal_clifford_phase_needs_degree_two_allowance(self):
        # Quarter parity is an S phase, including repeated translations.
        # Half a cross-address product is a literal CZ with no source call.
        phase = lambda point: _parity(point & 0b1011) + 2 * ((point >> 1 & 1)
                                                               * (point >> 3 & 1))
        for directions in ((1, 2, 8), (1, 1, 1), (3, 5, 6), (0, 1, 8)):
            for point in range(16):
                self.assertEqual(_cube(phase, point, directions) % 4, 0)
        cz = lambda point: Fraction((point & 1) * (point >> 1 & 1), 2)
        self.assertEqual(_cube(cz, 0, (1, 2)) % 1, Fraction(1, 2))
        self.assertEqual(_cube(cz, 0, (1, 2, 1)) % 1, 0)

    def test_cube_angle_threshold_and_moved_column_distance(self):
        # This saturates the abstract error-product inequality only; no
        # realizability or sharpness of the architecture's radius is claimed.
        for order in range(3, 9):
            radius = 2 * sin(pi / 2 ** (order + 1))
            product = 1 + 0j
            for corner in range(1 << order):
                sign = (-1) ** (order - corner.bit_count())
                error = exp(1j * sign * pi / (1 << order))
                self.assertAlmostEqual(abs(error - 1), radius, places=13)
                product *= error if sign == 1 else error.conjugate()
            self.assertAlmostEqual(abs(product + 1), 0, places=13)
            self.assertLess(radius, sqrt(2))
        # An affine permutation that moves even one column has this full
        # operator error lower bound, whatever the column's literal phase.
        for angle in (0, pi / 7, pi / 2, pi):
            actual = (0j, exp(1j * angle))
            for target_sign in (-1, 1):
                distance = sqrt(abs(actual[0] - target_sign)**2 + abs(actual[1])**2)
                self.assertAlmostEqual(distance, sqrt(2), places=14)

    def test_positive_error_echo_in_the_improved_radius(self):
        source = []
        source_error = 0.0
        for z in range(64):
            p = [(z >> i) & 1 for i in range(3)]
            q = [(z >> (i + 3)) & 1 for i in range(3)]
            cubic = (p[0] * q[1] * q[2] ^ p[1] * q[0] * q[2]
                     ^ p[2] * q[0] * q[1])
            perturbation = Fraction(((z * z + 3 * z) >> 1) % 7, 1000)
            source.append(Fraction(cubic, 2) + perturbation)
            source_error = max(source_error, abs(exp(2j * pi * float(perturbation)) - 1))
        complete_error = 0.0
        for x in range(8):
            target = (-1) ** int(x == 7)
            for z in range(64):
                returned, phase = _phase_echo_column(source, x, x << 3, z)
                self.assertEqual(returned, z)
                complete_error = max(complete_error,
                                     abs(exp(2j * pi * float(phase)) - target))
        self.assertGreater(complete_error, 2 / 275)
        self.assertLess(complete_error, 2 * sin(pi / 32))
        self.assertLessEqual(complete_error, 4 * source_error + 1e-12)


if __name__ == '__main__':
    unittest.main()
