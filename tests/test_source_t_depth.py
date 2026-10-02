"""Bounded native and exact-arithmetic checks of Majorana-linear schedules.

The unrestricted Clifford+T depth problem is not tested or settled here.
Full matrices check literal phases and all dirty inputs; exact coefficient
arithmetic checks the denominator witnesses independently of floating point.
"""
from __future__ import annotations

from fractions import Fraction
import unittest

import numpy as np

from tests.test_one_clean_compiler import _paired_source_word, _paired_weights
from tests.test_operator_source_compiler import (
    ATOL, X, Y, Z, _adjoint, _optimized_source_word, _pauli, _word_matrix,
)
from tests.test_t_depth import _inverse, _word


def _pauli_s_word(factors):
    """exp(-i*pi*P/4), up to a scalar that cancels with its inverse."""
    basis = []
    for wire, factor in sorted(factors.items()):
        if factor == 'X':
            basis.append(('H', wire))
        elif factor == 'Y':
            basis += [('SDG', wire), ('H', wire)]
    support = sorted(factors)
    parity = [('CX', wire, support[-1]) for wire in support[:-1]]
    return (basis + parity + [('S', support[-1])]
            + _adjoint(parity) + _adjoint(basis))


def _even_plane(wire):
    diagonalizer = _pauli_s_word({wire: 'X', wire + 1: 'X'})
    return [('C', diagonalizer), ('T', [('TDG', wire)]),
            ('C', _adjoint(diagonalizer))]


def _paired_tail_layer(wire):
    # D=G_34 G_13 G_01 in the local six-Majorana frame.
    diagonalizer = ([('S', wire)]
                    + _pauli_s_word({wire: 'X', wire + 1: 'Y'})
                    + _pauli_s_word({wire + 1: 'X', wire + 2: 'X'}))
    return [('C', diagonalizer),
            ('T', [('T', wire), ('TDG', wire + 1)]),
            ('C', _adjoint(diagonalizer))]


def _paired_loader(q):
    assert q >= 2
    schedule = [('T', [('T', 0)])] + _even_plane(0)
    for wire in range(q - 1):
        schedule += _paired_tail_layer(wire)
    return schedule


def _paired_seed_word():
    # (HS)^3=omega I; keep omega^-1 in C0=omega^-1 S X literally.
    inverse_scalar = _adjoint([('S', 0), ('H', 0)] * 3)
    return inverse_scalar + [('X', 0), ('S', 0)]


def _paired_source(q):
    tail = _paired_loader(q)[1:]
    return _inverse(tail) + [('C', _paired_seed_word())] + tail


def _majoranas(width):
    return [_pauli(width, {**{k: Z for k in range(j)}, j: factor})
            for j in range(width) for factor in (X, Y)]


def _add(left, right):
    return tuple(a + b for a, b in zip(left, right))


def _minus(left, right):
    return tuple(a - b for a, b in zip(left, right))


def _divide_root_two(value):
    rational, irrational = value
    return irrational, rational / 2


def _half_power(exponent):
    if exponent % 2:
        return Fraction(0), Fraction(1, 2 ** ((exponent + 1) // 2))
    return Fraction(1, 2 ** (exponent // 2)), Fraction(0)


def _sde(value):
    rational, irrational = value
    exponent = 0
    while rational.denominator != 1 or irrational.denominator != 1:
        rational, irrational = 2 * irrational, rational
        exponent += 1
    return exponent


def _reflection_diagonal(coefficient):
    rational, irrational = coefficient
    return 2 * (rational ** 2 + 2 * irrational ** 2) - 1, 4 * rational * irrational


def _symbolic_loader(width, layers):
    vector = [(Fraction(0), Fraction(0)) for _ in range(2 * width)]
    vector[0] = (Fraction(1), Fraction(0))
    for layer_number, pairs in enumerate(layers, 1):
        indices = [index for pair in pairs for index in pair]
        assert len(indices) == len(set(indices))
        for left, right in pairs:
            first, second = vector[left], vector[right]
            vector[left] = _divide_root_two(_minus(first, second))
            vector[right] = _divide_root_two(_add(first, second))
        assert max(map(_sde, vector)) <= layer_number
    return vector


class SourceTDepthTests(unittest.TestCase):
    def assert_schedule(self, schedule, width):
        for kind, gates in schedule:
            self.assertTrue(all(0 <= q < width for gate in gates for q in gate[1:]))
            if kind == 'T':
                self.assertTrue(all(gate[0] in ('T', 'TDG') for gate in gates))
                targets = [gate[1] for gate in gates]
                self.assertEqual(len(targets), len(set(targets)))
            else:
                self.assertTrue(all(gate[0] in ('X', 'CX', 'H', 'S', 'SDG', 'Z')
                                    for gate in gates))

    def test_literal_parallel_tail_block_and_majorana_cliffords(self):
        width = 3
        gammas = _majoranas(width)
        schedule = _paired_tail_layer(0)
        self.assert_schedule(schedule, width)
        actual = _word_matrix(width, _word(schedule))
        identity = np.eye(1 << width)
        odd = -_pauli(width, {0: X, 1: Y})
        even = _pauli(width, {1: Y, 2: X})
        cosine, sine = np.cos(np.pi / 8), np.sin(np.pi / 8)
        expected = ((cosine * identity + 1j * sine * even)
                    @ (cosine * identity + 1j * sine * odd))
        np.testing.assert_allclose(actual, expected, atol=ATOL, rtol=0)
        np.testing.assert_allclose(
            _word_matrix(width, _word(_inverse(schedule))) @ actual,
            identity, atol=ATOL, rtol=0)
        # Complete Clifford blocks, not each elementary H/CX, preserve the span.
        for gates in (schedule[0][1], _even_plane(0)[0][1],
                      _paired_seed_word(), _optimized_source_word(2)):
            clifford = _word_matrix(width, gates)
            images = []
            for gamma in gammas:
                transformed = clifford @ gamma @ clifford.conj().T
                overlaps = [np.trace(other @ transformed).real / (1 << width)
                            for other in gammas]
                index = int(np.argmax(np.abs(overlaps)))
                self.assertAlmostEqual(abs(overlaps[index]), 1)
                images.append(index)
                np.testing.assert_allclose(transformed,
                                           np.sign(overlaps[index]) * gammas[index],
                                           atol=ATOL, rtol=0)
            self.assertEqual(len(set(images)), len(gammas))

    def test_paired_loader_and_absorbed_source_literal_full_matrices(self):
        for q in (2, 3, 4):
            width = q + 1
            loader, source = _paired_loader(q), _paired_source(q)
            for schedule in (loader, source):
                self.assert_schedule(schedule, width)
            actual_loader = _word_matrix(width, _word(loader))
            # Independent prior serial word: equality includes global phase.
            np.testing.assert_allclose(actual_loader,
                                       _word_matrix(width, _paired_source_word(q)),
                                       atol=ATOL, rtol=0)
            weights, _ = _paired_weights(q)
            expected_source = sum(np.sqrt(weight) * gamma
                                  for weight, gamma in zip(weights, _majoranas(width)))
            actual_source = _word_matrix(width, _word(source))
            np.testing.assert_allclose(actual_source, expected_source,
                                       atol=ATOL, rtol=0)
            np.testing.assert_allclose(
                _word_matrix(width, _word(_inverse(source))) @ actual_source,
                np.eye(1 << width), atol=ATOL, rtol=0)
            # A spectator is arbitrary input, not an initialized helper.
            if q == 2:
                np.testing.assert_allclose(_word_matrix(width + 1, _word(source)),
                                           np.kron(np.eye(2), expected_source),
                                           atol=ATOL, rtol=0)

    def test_exact_coefficients_denominators_and_schedule_resources(self):
        zero = (Fraction(0), Fraction(0))
        for m in range(2, 19):
            layers = [[(2 * j, 2 * j + 2)] for j in range(m - 1)]
            vector = _symbolic_loader(m, layers)
            expected = [zero for _ in range(2 * m)]
            for j in range(m):
                expected[2 * j] = _half_power(min(j + 1, m - 1))
            self.assertEqual(vector, expected)
            self.assertEqual(max(map(_sde, vector)), m - 1)
            witness = _reflection_diagonal(vector[2 * (m - 1)])
            self.assertEqual(witness, (Fraction(1, 2 ** (m - 2)) - 1, Fraction(0)))
            self.assertEqual(_sde(witness), 2 * m - 4)
        for q in range(2, 19):
            layers = [[(0, 1)], [(0, 2)]] + [
                [(2 * j + 1, 2 * j + 3), (2 * j + 2, 2 * j + 4)]
                for j in range(q - 1)]
            vector = _symbolic_loader(q + 1, layers)
            expected = [zero for _ in range(2 * (q + 1))]
            expected[0] = (Fraction(1, 2), Fraction(0))
            for j in range(q):
                expected[2 * j + 1] = _half_power(min(j + 2, q))
                expected[2 * j + 2] = _half_power(min(j + 3, q + 1))
            self.assertEqual(vector, expected)
            self.assertEqual(max(map(_sde, vector)), q + 1)
            witness = _reflection_diagonal(vector[2 * q])
            self.assertEqual(witness, (Fraction(1, 2 ** q) - 1, Fraction(0)))
            self.assertEqual(_sde(witness), 2 * q)
            for schedule, depth, count in ((_paired_loader(q), q + 1, 2 * q),
                                            (_paired_source(q), 2 * q, 4 * q - 2)):
                self.assert_schedule(schedule, q + 1)
                self.assertEqual(sum(kind == 'T' for kind, _ in schedule), depth)
                self.assertEqual(sum(len(gates) for kind, gates in schedule
                                     if kind == 'T'), count)
                self.assertLessEqual(sum(len(gates) for kind, gates in schedule
                                         if kind == 'C'), 80 * q)

    def test_geometric_native_planes_and_absorbed_source(self):
        for width in (2, 3, 4):
            loader = [stage for wire in range(width - 1)
                      for stage in _even_plane(wire)]
            tail = loader[3:]
            center = _optimized_source_word(2)
            source = _inverse(tail) + [('C', center)] + tail
            self.assert_schedule(source, width)
            gammas = _majoranas(width)
            expected = sum(2 ** (-min(j + 1, width - 1) / 2) * gammas[2 * j]
                           for j in range(width))
            unitary = _word_matrix(width, _word(loader))
            np.testing.assert_allclose(unitary @ gammas[0] @ unitary.conj().T,
                                       expected, atol=ATOL, rtol=0)
            np.testing.assert_allclose(_word_matrix(width, _word(source)), expected,
                                       atol=ATOL, rtol=0)
            self.assertEqual(sum(kind == 'T' for kind, _ in loader), width - 1)
            self.assertEqual(sum(kind == 'T' for kind, _ in source), 2 * width - 4)


if __name__ == '__main__':
    unittest.main()
