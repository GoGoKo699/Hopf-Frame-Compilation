"""Finite algebra audits for the clean-aware shallow-depth obstruction.

These checks cover initialized-projector normalization, actual complete-frame
witnesses, explicit shallow-depth gaps, terminal-sector rigidity, and
mod-eight phase digit identities.
They do not certify the imported polynomial-bias closure theorem or infer
uniform infinite-width statements from finite matrices.
"""
from __future__ import annotations

from decimal import Decimal, localcontext
from fractions import Fraction
from itertools import combinations, product
import unittest

import numpy as np

from compiler_robust_hopf.frames import (
    direct_real_frame, hopf_ry, in_canonical_magnitude_domain,
    real_frame_matrix,
)
from tests.test_operator_source_compiler import (
    ATOL, I2, X, Y, Z, _adjoint, _pauli, _word_matrix,
)
from tests.test_shallow_source_obstruction import _clifford_word, _pauli_basis


def _initialized_transfer(unitary, logical, dirty, clean):
    """Trace out all output work, averaging every dirty input column."""
    data_width = logical + dirty
    embedding = np.eye(1 << (data_width + clean))[:, :1 << data_width]
    isometry = unitary @ embedding
    data_paulis = np.array([np.kron(np.eye(1 << dirty), pauli)
                            for pauli in _pauli_basis(logical)])
    output_paulis = np.array([np.kron(np.eye(1 << (dirty + clean)), pauli)
                              for pauli in _pauli_basis(logical)])
    evolved = isometry @ data_paulis @ isometry.conj().T
    return np.einsum('pij,qji->pq', output_paulis, evolved) / (1 << data_width)


def _quadratic_decimal(value):
    rational, radical = value
    return (Decimal(rational.numerator) / Decimal(rational.denominator)
            + Decimal(radical.numerator) / Decimal(radical.denominator)
            * Decimal(2).sqrt())


def _f2_product(left, right):
    """Squarefree F2 polynomials, represented by sets of monomial masks."""
    result = set()
    for first in left:
        for second in right:
            result.symmetric_difference_update((first | second,))
    return result


def _phase_digits(coefficients):
    # The coefficient c means c repeated copies of its monomial bit.
    # Elementary symmetric polynomial e_(2^j) is integer-sum digit j.
    elementary = [{0}, set(), set(), set(), set()]
    for monomial, coefficient in coefficients.items():
        for _ in range(coefficient % 8):
            for degree in range(4, 0, -1):
                elementary[degree] ^= _f2_product(elementary[degree - 1],
                                                 {monomial})
    return [elementary[1], elementary[2], elementary[4]]


def _phase_events(digits):
    result = []
    for residue in range(8):
        event = {0}
        for bit, digit in enumerate(digits):
            factor = digit if (residue >> bit) & 1 else digit ^ {0}
            event = _f2_product(event, factor)
        result.append(event)
    return result


def _f2_value(polynomial, assignment):
    return sum((assignment & mask) == mask for mask in polynomial) % 2


def _degree(polynomial):
    return max((mask.bit_count() for mask in polynomial), default=0)


def _integer_mobius(values, width):
    coefficients = list(values)
    for bit in range(width):
        for mask in range(1 << width):
            if mask & (1 << bit):
                coefficients[mask] -= coefficients[mask ^ (1 << bit)]
    return {mask: value % 8 for mask, value in enumerate(coefficients)
            if value % 8}


class CleanAwareDepthObstructionTests(unittest.TestCase):
    def test_native_projector_sum_matches_all_initialized_transfer_entries(self):
        rng = np.random.default_rng(73129)
        saw_non_alphabet_sum = False
        for logical, clean in product((1, 2), repeat=2):
            dirty = 1
            width = logical + dirty + clean
            dimension = 1 << width
            logical_factors = list(product((I2, X, Y, Z), repeat=logical))
            output_paulis = np.array([
                _pauli(width, dict(enumerate(factors)))
                for factors in logical_factors])
            alphabet = np.array([0.0] + [sign * 2 ** (-j / 2)
                                for sign in (-1, 1)
                                for j in range(2 * width + 1)])
            for depth in (0, 1, 2):
                for sample in range(3):
                    word = _clifford_word(width, rng)
                    for _ in range(depth):
                        layer = [('T' if rng.integers(2) else 'TDG', wire)
                                 for wire in range(width)]
                        word += layer + _clifford_word(width, rng)
                    unitary = _word_matrix(width, word)
                    direct = _initialized_transfer(unitary, logical, dirty, clean)
                    reconstructed = np.zeros_like(direct)
                    for clean_bits in product((0, 1), repeat=clean):
                        clean_factors = {logical + dirty + wire: Z
                                         for wire, bit in enumerate(clean_bits) if bit}
                        inputs = np.array([
                            _pauli(width, {**dict(enumerate(factors)), **clean_factors})
                            for factors in logical_factors])
                        evolved = unitary @ inputs @ unitary.conj().T
                        summand = np.einsum('pij,qji->pq', output_paulis,
                                           evolved) / dimension
                        self.assertLess(np.max(np.abs(summand.imag)), ATOL)
                        distances = np.min(np.abs(summand.real[..., None] - alphabet),
                                           axis=-1)
                        self.assertLess(np.max(distances), ATOL)
                        reconstructed += summand
                    with self.subTest(logical=logical, clean=clean,
                                      depth=depth, sample=sample):
                        np.testing.assert_allclose(reconstructed, direct,
                                                   atol=ATOL, rtol=0)
                        self.assertAlmostEqual(direct[0, 0].real, 1.0)
                    sum_distances = np.min(np.abs(direct.real[..., None] - alphabet),
                                           axis=-1)
                    saw_non_alphabet_sum |= np.max(sum_distances) > 0.01
        self.assertTrue(saw_non_alphabet_sum,
                        'Fixtures must distinguish the sum alphabet from a single entry.')

    def test_robustness_uses_initialized_columns_with_dirty_and_reference(self):
        target = np.kron(I2, hopf_ry(0.271))
        embedding = np.eye(8)[:, :4]
        ideal = embedding @ target
        # The uninitialized sector is deliberately far from the target extension.
        extension = np.kron(Z, target)
        generator = _pauli(3, {0: X, 1: Z, 2: X})
        logical_paulis = _pauli_basis(1)
        expected = np.array([[np.trace(p @ hopf_ry(0.271) @ q
                                      @ hopf_ry(0.271).conj().T) / 2
                              for q in logical_paulis] for p in logical_paulis])
        for angle in (0.0, 0.013, 0.19):
            perturbation = np.cos(angle) * np.eye(8) - 1j * np.sin(angle) * generator
            actual = perturbation @ extension
            difference = actual @ embedding - ideal
            eta = np.linalg.norm(difference, ord=2)
            self.assertAlmostEqual(eta, 2 * np.sin(angle / 2), places=12)
            self.assertGreater(np.linalg.norm(actual - np.kron(I2, target), ord=2), 1.9)
            self.assertAlmostEqual(np.linalg.norm(np.kron(difference, I2), ord=2),
                                   eta, places=12)
            observed = _initialized_transfer(actual, logical=1, dirty=1, clean=1)
            self.assertLessEqual(np.max(np.abs(observed - expected)), 2 * eta + ATOL)

    def test_explicit_witness_is_the_complete_frame_at_every_tested_dimension(self):
        for logical, coefficient in product(range(1, 5), (Fraction(13, 16),
                                                         Fraction(1, 3))):
            theta = 0.5 * np.arccos(float(coefficient))
            angles = np.zeros((1 << logical) - 1)
            angles[(1 << (logical - 1)) - 1:] = theta
            self.assertTrue(in_canonical_magnitude_domain(angles))
            self.assertTrue(in_canonical_magnitude_domain(angles, complex_chart=True))
            expected = np.kron(np.eye(1 << (logical - 1)), hopf_ry(theta))
            frame = real_frame_matrix(angles)
            np.testing.assert_allclose(frame, expected, atol=ATOL, rtol=0)
            np.testing.assert_allclose(direct_real_frame(logical, angles), expected,
                                       atol=ATOL, rtol=0)
            witness = _pauli(logical, {0: Z})
            observed = np.trace(witness @ frame @ witness @ frame.conj().T) / len(frame)
            self.assertLess(abs(observed - float(coefficient)), ATOL)

    def test_two_sum_gap_has_exact_nearest_pair_and_certified_infinite_tail(self):
        with localcontext() as context:
            context.prec = 90
            tolerance = Decimal('1e-75')
            cutoff = 12
            alphabet = {(Fraction(0), Fraction(0))}
            for exponent in range(cutoff + 1):
                for sign in (-1, 1):
                    if exponent % 2:
                        alphabet.add((Fraction(0), Fraction(sign, 2 ** ((exponent + 1) // 2))))
                    else:
                        alphabet.add((Fraction(sign, 2 ** (exponent // 2)), Fraction(0)))
            sums = {(p + r, q + s) for p, q in alphabet for r, s in alphabet}
            target = Decimal(13) / 16
            nearest = min(sums, key=lambda value: abs(_quadratic_decimal(value) - target))
            self.assertEqual(nearest, (Fraction(1), Fraction(-1, 8)))
            gap = abs(_quadratic_decimal(nearest) - target)
            expected_gap = (3 - 2 * Decimal(2).sqrt()) / 16
            self.assertLess(abs(gap - expected_gap), tolerance)
            self.assertEqual((nearest[0] - Fraction(13, 16), nearest[1]),
                             (Fraction(3, 16), Fraction(-1, 8)))
            # Any omitted alphabet element has magnitude at most tail. First
            # certify the single-entry distance, then exclude every pair with
            # at least one omitted term by the reverse triangle inequality.
            tail = Decimal(2).sqrt() ** (-(cutoff + 1))
            single_gap = min(abs(_quadratic_decimal(point) - target) for point in alphabet)
            self.assertGreater(target - tail, single_gap)
            self.assertGreater(single_gap - tail, gap)
            self.assertGreater(gap / 2, Decimal(1) / 256)

    def test_native_three_layer_counterexample_to_the_two_sum_alphabet(self):
        for fanout_width in range(1, 5):
            width = fanout_width + 2
            fanout_wires = range(2, width)
            axes = [
                _pauli(width, {0: Y, 1: Z}),
                _pauli(width, {0: Y, **{wire: X for wire in fanout_wires}}),
                _pauli(width, {0: Y, 1: Z, **{wire: X for wire in fanout_wires}}),
            ]
            # First reduce X_B to X_2, then Y_0 and X_2 to Z. The final
            # three CNOTs send parity labels (110),(101),(111) to distinct
            # singleton Z labels. The matrix checks include the literal signs.
            clifford_word = ([('CX', 2, wire) for wire in range(3, width)]
                              + [('SDG', 0), ('H', 0), ('H', 2),
                                 ('CX', 1, 0), ('CX', 0, 2), ('CX', 2, 1)])
            clifford = _word_matrix(width, clifford_word)
            for wire, axis in enumerate(axes):
                np.testing.assert_allclose(clifford @ axis @ clifford.conj().T,
                                           _pauli(width, {wire: Z}), atol=ATOL, rtol=0)
            outer_layer = [('T', 0), ('T', 1), ('TDG', 2)]
            source_word = clifford_word + outer_layer + _adjoint(clifford_word)
            source = _word_matrix(width, source_word)
            pauli = _pauli(width, {0: X})
            fanout = _word_matrix(width, [('CX', 1, wire) for wire in fanout_wires])
            reflected_fanout = (_pauli(width, {0: X}) - _pauli(width, {0: Z})) @ fanout / np.sqrt(2)
            np.testing.assert_allclose(source @ pauli @ source.conj().T,
                                       reflected_fanout, atol=ATOL, rtol=0)
            middle_layer = [('T', wire) for wire in fanout_wires]
            layers = (outer_layer, middle_layer, _adjoint(outer_layer))
            for layer in layers:
                self.assertEqual(len(layer), len({gate[1] for gate in layer}))
            word = source_word + middle_layer + _adjoint(source_word)
            self.assertEqual(sum(gate[0] in ('T', 'TDG') for gate in word),
                             fanout_width + 6)
            unitary = _word_matrix(width, word)
            commutator = pauli @ unitary @ pauli @ unitary.conj().T
            coefficient = np.trace(commutator) / (1 << width)
            expected = (1 + 2 ** (-fanout_width / 2)) / 2
            self.assertLess(abs(coefficient - expected), ATOL)
            if fanout_width == 1:
                self.assertAlmostEqual(1 - coefficient.real,
                                       (1 - 1 / np.sqrt(2)) / 2, places=12)
            if fanout_width == 4:
                # Evaluate the disjoint tensor witness by trace factorization;
                # no 4096-by-4096 dense tensor matrix is needed. Combining the
                # disjoint physical layers preserves their depth of three.
                for layer in layers:
                    combined = layer + [(name, wire + width) for name, wire in layer]
                    self.assertEqual(len(combined), len({gate[1] for gate in combined}))
                tensor_coefficient = np.einsum('ii,jj->', commutator, commutator) / (1 << (2 * width))
                target = Fraction(25, 64)
                self.assertEqual(Fraction(5, 8) ** 2, target)
                self.assertLess(abs(tensor_coefficient - float(target)), ATOL)

        # A nonzero rational two-term sum cannot contain an uncancelled
        # sqrt(2) term; two irrational terms cancel only to zero. Thus exact
        # membership reduces to signed sums of two rational binary powers.
        cutoff = 12
        rational_alphabet = {Fraction(0)} | {Fraction(sign, 2 ** exponent)
                                           for sign in (-1, 1)
                                           for exponent in range(cutoff + 1)}
        pair_sums = {first + second for first in rational_alphabet
                     for second in rational_alphabet}
        self.assertNotIn(target, pair_sums)
        finite_gap = min(abs(target - value) for value in pair_sums)
        single_gap = min(abs(target - value) for value in rational_alphabet)
        tail = Fraction(1, 2 ** (cutoff + 1))
        self.assertGreater(target - tail, single_gap)
        self.assertGreater(single_gap - tail, finite_gap)
        self.assertGreater(finite_gap, 0)

    def test_weighted_phase_parity_strata_and_sharp_gap(self):
        delta = (1 - 1 / np.sqrt(2)) / 2
        width = 4
        rng = np.random.default_rng(73074)
        phases = [({}, 0), ({1: 1, 3: 6}, delta),
                  ({3: 2}, 0.25), ({7: 4}, 0.25), ({0: 1}, 2 * delta)]
        for _ in range(24):
            coefficients = {0: int(rng.integers(8))}
            for degree in (1, 2, 3):
                for support in combinations(range(width), degree):
                    mask = sum(1 << bit for bit in support)
                    coefficients[mask] = ((1 << (degree - 1))
                                          * int(rng.integers(8))) % 8
            phases.append((coefficients, None))
        for coefficients, expected_loss in phases:
            values = np.array([sum(coefficient for mask, coefficient
                                   in coefficients.items() if assignment & mask == mask) % 8
                               for assignment in range(1 << width)])
            loss = np.mean(1 - np.cos(np.pi * values / 4))
            if expected_loss is not None:
                self.assertAlmostEqual(loss, expected_loss, places=12)
            if not np.any(values):
                self.assertAlmostEqual(loss, 0)
                continue
            # Audit the three valuation strata, including the phase whose
            # real expectation attains the sharp scalar upper bound.
            for divisor, degree, support_bound in ((1, 1, 0.5),
                                                    (2, 2, 0.25),
                                                    (4, 3, 0.125)):
                parity = (values // divisor) % 2
                if np.any(parity):
                    parity_coefficients = _integer_mobius(parity.tolist(), width)
                    boolean_support = {mask for mask, value in parity_coefficients.items()
                                       if value % 2}
                    self.assertLessEqual(_degree(boolean_support), degree)
                    self.assertGreaterEqual(np.mean(parity), support_bound)
                    break
                self.assertTrue(all(coefficient % (2 * divisor) == 0
                                    for coefficient in coefficients.values()))
            self.assertGreaterEqual(loss + ATOL, delta)

    def test_native_three_layer_transfer_entries_have_isolated_unit_maximum(self):
        rng = np.random.default_rng(73075)
        upper = (1 + 1 / np.sqrt(2)) / 2
        for width in range(2, 4):
            paulis = np.array(_pauli_basis(width))
            for _ in range(3):
                word = _clifford_word(width, rng)
                for _ in range(3):
                    word += [('T' if rng.integers(2) else 'TDG', wire)
                             for wire in range(width)] + _clifford_word(width, rng)
                unitary = _word_matrix(width, word)
                evolved = unitary @ paulis @ unitary.conj().T
                transfer = np.einsum('pij,qji->pq', paulis, evolved) / (1 << width)
                self.assertLess(np.max(np.abs(transfer.imag)), ATOL)
                nonunit = transfer.real[np.abs(transfer.real - 1) > ATOL]
                self.assertLessEqual(np.max(nonunit), upper + ATOL)

    def test_zero_clean_sharp_floor_has_literal_complete_frame_attainer(self):
        theta = 0.5 * np.arccos(1 / 3)
        c0 = np.sqrt(2 - (2 + np.sqrt(2)) / np.sqrt(3))
        rho = np.sqrt((1 - 1 / np.sqrt(2)) / 4)
        self.assertLess(1 / 8, c0)
        self.assertLess(c0, rho)
        # PTM angle quantization permits alpha=k*pi/8 up to phase.
        # Check that the nearest angle is pi/4, not the half-angle pi/8.
        axial_distances = [2 * np.sin(abs((theta - k * np.pi / 8 + np.pi / 2)
                                          % np.pi - np.pi / 2) / 2)
                           for k in range(8)]
        self.assertEqual(int(np.argmin(axial_distances)), 2)
        self.assertAlmostEqual(min(axial_distances), c0, places=12)
        for logical, dirty in product(range(1, 5), range(3)):
            angles = np.zeros((1 << logical) - 1)
            angles[(1 << (logical - 1)) - 1:] = theta
            target = np.kron(np.eye(1 << dirty), real_frame_matrix(angles))
            width = logical + dirty
            # Chronological Z,H is literally HZ=Ry(pi/4), without a
            # projective-phase adjustment. Every dirty input is retained.
            attainer = _word_matrix(width, [('Z', 0), ('H', 0)])
            expected = np.kron(np.eye(1 << (width - 1)), hopf_ry(np.pi / 4))
            np.testing.assert_allclose(attainer, expected, atol=ATOL, rtol=0)
            difference = attainer - target
            self.assertAlmostEqual(np.linalg.norm(difference, ord=2), c0, places=12)
            np.testing.assert_allclose(difference.conj().T @ difference,
                                       c0 ** 2 * np.eye(1 << width), atol=ATOL, rtol=0)

    def test_one_clean_block_norms_and_occupied_sector_caveat(self):
        target = np.kron(I2, hopf_ry(0.271))
        embedding = np.eye(8)[:, :4]
        extension = np.zeros((8, 8), dtype=complex)
        extension[:4, :4] = target
        extension[4:, 4:] = _word_matrix(2, [('H', 1), ('CX', 1, 0), ('T', 0)])
        clean_z = _pauli(3, {2: Z})
        generator = _pauli(3, {0: Z, 2: X})
        for angle in (0, 0.031, 0.19):
            # These continuous perturbations audit the general block identity,
            # not membership in the at-most-three-layer native circuit class.
            perturbation = np.cos(angle) * np.eye(8) - 1j * np.sin(angle) * generator
            unitary = perturbation @ extension
            eta = np.linalg.norm(unitary @ embedding - embedding @ target, ord=2)
            b_norm = np.linalg.norm(unitary[:4, 4:], ord=2)
            c_norm = np.linalg.norm(unitary[4:, :4], ord=2)
            self.assertAlmostEqual(b_norm, c_norm, places=12)
            commutator = unitary @ clean_z - clean_z @ unitary
            self.assertAlmostEqual(np.linalg.norm(commutator, ord=2), 2 * c_norm,
                                   places=12)
            self.assertLessEqual(c_norm, eta + ATOL)
            transfer = np.trace(clean_z @ unitary @ clean_z @ unitary.conj().T) / 8
            self.assertLessEqual(1 - transfer.real, 2 * eta ** 2 + ATOL)
            self.assertAlmostEqual(1 - transfer.real,
                                   np.linalg.norm(commutator, 'fro') ** 2 / 16, places=12)

        # Exact initialized return does not imply global logical-Y
        # conservation: the occupied clean sector may act differently.
        controlled_x = _word_matrix(3, [('CX', 2, 0)])
        np.testing.assert_allclose(controlled_x @ embedding, embedding, atol=ATOL, rtol=0)
        np.testing.assert_allclose(controlled_x @ clean_z, clean_z @ controlled_x,
                                   atol=ATOL, rtol=0)
        logical_y = _pauli(3, {0: Y})
        evolved_y = controlled_x.conj().T @ logical_y @ controlled_x
        np.testing.assert_allclose(evolved_y, clean_z @ logical_y, atol=ATOL, rtol=0)
        self.assertAlmostEqual(np.linalg.norm(evolved_y - logical_y, ord=2), 2)

    def test_phase_digits_after_affine_restriction_and_mod_eight_carries(self):
        rng = np.random.default_rng(73021)
        width = 7
        original_width = 6
        for sample in range(4):
            affine = [(int(rng.integers(1 << width)), int(rng.integers(2)))
                      for _ in range(original_width)]
            if sample == 0:
                # Include a constant coordinate and dependent coordinates,
                # not only generic affine maps of nearly full rank.
                affine[0] = (0, 1)
                affine[1] = affine[2]
            phase_terms = [(0, int(rng.integers(8)))]
            for degree in (1, 2, 3):
                for support in combinations(range(original_width), degree):
                    mask = sum(1 << bit for bit in support)
                    coefficient = (1 << (degree - 1)) * int(rng.integers(8))
                    phase_terms.append((mask, coefficient))
            values = []
            for assignment in range(1 << width):
                original = sum((((mask & assignment).bit_count() % 2) ^ offset) << bit
                               for bit, (mask, offset) in enumerate(affine))
                values.append(sum(coefficient for mask, coefficient in phase_terms
                                  if original & mask == mask) % 8)
            coefficients = _integer_mobius(values, width)
            self.assertLessEqual(max(map(int.bit_count, coefficients), default=0), 3)
            for mask, coefficient in coefficients.items():
                self.assertEqual(coefficient % (1 << max(0, mask.bit_count() - 1)), 0)
            digits = _phase_digits(coefficients)
            for digit, bound in zip(digits, (3, 6, 12)):
                self.assertLessEqual(_degree(digit), bound)
            events = _phase_events(digits)
            self.assertTrue(all(_degree(event) <= 21 for event in events))
            for assignment, value in enumerate(values):
                reconstructed = sum(_f2_value(digit, assignment) << bit
                                    for bit, digit in enumerate(digits))
                with self.subTest(sample=sample, assignment=assignment):
                    self.assertEqual(reconstructed, value)
                    self.assertEqual([_f2_value(event, assignment) for event in events],
                                     [int(value == residue) for residue in range(8)])

    def test_safe_event_degree_bound_on_seven_disjoint_cubic_tokens(self):
        coefficients = {7 << (3 * token): 1 for token in range(7)}
        digits = _phase_digits(coefficients)
        self.assertEqual([_degree(digit) for digit in digits], [3, 6, 12])
        events = _phase_events(digits)
        self.assertEqual(_degree(events[7]), 21)
        # Every seven-token truth assignment is represented without enumerating
        # all 2^21 inputs. This audits the safe generic-cubic bound, not a claim
        # of sharpness for the weighted phases in the circuit reduction.
        for token_values in range(1 << 7):
            assignment = sum((6 | ((token_values >> token) & 1)) << (3 * token)
                             for token in range(7))
            value = token_values.bit_count() % 8
            self.assertEqual([_f2_value(event, assignment) for event in events],
                             [int(value == residue) for residue in range(8)])


if __name__ == '__main__':
    unittest.main()
