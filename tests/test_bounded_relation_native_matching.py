"""Exact finite checks for bounded-relation native axis matching.

The proof supplies the arbitrary-width estimates. These fixtures retain
dependent axes, literal signed T/T-dagger phases and actual inverses. A
rank-profile symmetry is deliberately distinguished from Clifford passage.
"""
from itertools import product
import unittest

from tests.exact_matrices import CZ, CO, adjoint, native_columns
from tests.test_two_layer_galois_return import (
    _clifford_columns, _conjugate, _inverse, _product, _tau,
)


def _rank(columns):
    pivots = {}
    for column in columns:
        while column:
            bit = column.bit_length() - 1
            if bit not in pivots:
                pivots[bit] = column
                break
            column ^= pivots[bit]
    return len(pivots)


def _apply(columns, vector):
    image = 0
    for bit, column in enumerate(columns):
        if vector >> bit & 1:
            image ^= column
    return image


def _dot(left, right):
    return (left & right).bit_count() & 1


def _relation_rank_formula(axes, probe):
    active = [axis for axis in axes if _dot(axis, probe)]
    relations = [v for v in range(1 << len(active))
                 if _apply(active, v) == 0]
    radical = [v for v in relations
               if all(_dot(v, other) == 0 for other in relations)]
    relation_dimension = len(relations).bit_length() - 1
    radical_dimension = len(radical).bit_length() - 1
    return len(active) - relation_dimension - radical_dimension


def _phase(axes, width, signs=0, daggers=0):
    # A negative signed axis gives the complementary parity, including
    # its literal scalar phase. T-dagger negates that integer exponent.
    return [sum((-1 if daggers >> j & 1 else 1)
                * (_dot(axis, state) ^ ((signs >> j) & 1))
                for j, axis in enumerate(axes)) % 8
            for state in range(1 << width)]


def _phase_shear(phase, probe, width):
    """Recover the exact symplectic shear of D X_probe D-dagger."""
    derivative = [(phase[state ^ probe] - phase[state]) % 8
                  for state in range(1 << width)]
    columns = []
    for bit in range(width):
        difference = [(derivative[state ^ (1 << bit)] - derivative[state]) % 8
                      for state in range(1 << width)]
        column = 0
        for other in range(width):
            coefficient = (difference[1 << other] - difference[0]) % 8
            if coefficient not in (0, 4):
                raise ValueError('The image is not a Clifford phase')
            column |= (coefficient // 4) << other
        if any(value != (difference[0] + 4 * _dot(column, state)) % 8
               for state, value in enumerate(difference)):
            raise ValueError('The image is not a Clifford phase')
        columns.append(column)
    return columns


def _native_phase_word(axes, signs=0, daggers=0):
    word = []
    for index, axis in enumerate(axes):
        target = axis.bit_length() - 1
        parity = [('CX', bit, target) for bit in range(target)
                  if axis >> bit & 1]
        flip = ([('H', target), ('S', target), ('S', target), ('H', target)]
                if signs >> index & 1 else [])
        word += parity + flip
        word += [('TDG' if daggers >> index & 1 else 'T', target)]
        word += flip + list(reversed(parity))
    return word


class BoundedRelationNativeMatchingTests(unittest.TestCase):
    def test_relation_radical_formula_with_literal_signed_native_phases(self):
        width = 3
        # All distinct nonempty Z-axis lists in this three-dimensional
        # isotropic space, including kernels with nontrivial radicals.
        for subset in range(1, 128):
            axes = [axis for axis in range(1, 8) if subset >> (axis - 1) & 1]
            relation_rank = len(axes) - _rank(axes)
            phase = _phase(axes, width)
            for probe in range(8):
                rank = _rank(_phase_shear(phase, probe, width))
                self.assertEqual(rank, _relation_rank_formula(axes, probe))
                active_count = sum(_dot(axis, probe) for axis in axes)
                self.assertLessEqual(active_count - rank, 2 * relation_rank)
                self.assertLessEqual(rank, active_count)

        axes = (1, 2, 4, 7)
        self.assertEqual(_relation_rank_formula(axes, 7), 2)
        self.assertEqual(sum(_dot(axis, 7) for axis in axes), 4)
        # Exhaust every axis orientation and T/T-dagger choice separately.
        for signs, daggers in product(range(16), repeat=2):
            phase = _phase(axes, width, signs, daggers)
            for probe in range(8):
                self.assertEqual(_rank(_phase_shear(phase, probe, width)),
                                 _relation_rank_formula(axes, probe))

        omega = native_columns(1, [('T', 0)])[1][1]
        roots = [CO]
        for _ in range(7):
            roots.append(roots[-1] * omega)
        for signs, daggers in ((0, 0), (5, 10), (15, 15)):
            word = _native_phase_word(axes, signs, daggers)
            source = native_columns(width, word)
            phase = _phase(axes, width, signs, daggers)
            self.assertEqual(source, [[roots[phase[column]] if row == column else CZ
                                       for column in range(8)] for row in range(8)])
            self.assertEqual(native_columns(width, _inverse(word)), adjoint(source))
            for probe in ((1, 7) if signs == 0 else (7,)):
                image = _conjugate(source, _tau(probe, width))
                columns = _clifford_columns(image, width)
                rank = _rank([value ^ (1 << bit)
                              for bit, value in enumerate(columns)])
                self.assertEqual(rank, _relation_rank_formula(axes, probe))

    def test_relation_four_symmetries_add_axes_to_the_majority_proxy(self):
        axes = frozenset((1, 2, 4, 7))
        relation_rank = len(axes) - _rank(axes)
        self.assertEqual(relation_rank, 1)
        profile = [_relation_rank_formula(tuple(axes), probe) for probe in range(8)]
        self.assertEqual(profile, [0] + [2] * 7)
        group = [columns for columns in product(range(1, 8), repeat=3)
                 if _rank(columns) == 3]
        self.assertEqual(len(group), 168)
        memberships = {axis: 0 for axis in range(1, 8)}
        differences = []
        for columns in group:
            transformed = {_apply(columns, axis) for axis in axes}
            self.assertEqual([_relation_rank_formula(tuple(transformed), probe)
                              for probe in range(8)], profile)
            differences.append(len(axes ^ transformed))
            self.assertLessEqual(differences[-1], 16 * relation_rank ** 2)
            for axis in transformed:
                memberships[axis] += 1
        self.assertEqual(set(memberships.values()), {96})  # 96/168 = 4/7.
        self.assertEqual(max(differences), 4)
        majority = {axis for axis, count in memberships.items()
                    if 2 * count > len(group)}
        self.assertEqual(majority, set(range(1, 8)))
        self.assertEqual(len(axes ^ majority), 3)
        self.assertLessEqual(len(axes ^ majority), 32 * relation_rank ** 2)
        self.assertEqual(len(majority) - _rank(majority), 4)
        for columns in group:
            self.assertEqual({_apply(columns, axis) for axis in majority}, majority)
        # Every proxy axis remains in the original commuting Z span.
        self.assertEqual(_rank((*axes, *majority)), _rank(axes))

    def test_profile_symmetry_does_not_supply_a_clifford_passage(self):
        axes, alternate = (1, 2, 4, 7), (1, 3, 4, 6)
        axis_map = (1, 3, 4)  # Z action of the native CX(0;1).
        self.assertEqual({_apply(axis_map, axis) for axis in axes}, set(alternate))
        self.assertEqual([_relation_rank_formula(axes, v) for v in range(8)],
                         [_relation_rank_formula(alternate, v) for v in range(8)])
        phase, other_phase = _phase(axes, 3), _phase(alternate, 3)
        difference = [(left - right) % 8 for left, right in zip(phase, other_phase)]
        for state in range(8):
            x0, x1, x2 = [(state >> bit) & 1 for bit in range(3)]
            self.assertEqual(difference[state], (6 * x0 * x2 + 4 * x0 * x1 * x2) % 8)

        word, other_word = _native_phase_word(axes), _native_phase_word(alternate)
        source, other = native_columns(3, word), native_columns(3, other_word)
        inner = native_columns(3, [('CX', 0, 1)])
        self.assertEqual(_conjugate(inner, source), other)
        correction = native_columns(3, _inverse(other_word) + word)
        self.assertEqual(correction, _product(source, adjoint(other)))
        with self.assertRaises(ValueError):
            _clifford_columns(correction, 3)
        passage_word = _inverse(word) + [('CX', 0, 1)] + word
        passage = native_columns(3, passage_word)
        self.assertEqual(passage, _product(correction, inner))
        self.assertEqual(native_columns(3, _inverse(passage_word)), adjoint(passage))
        with self.assertRaises(ValueError):
            _clifford_columns(passage, 3)


if __name__ == '__main__':
    unittest.main()
