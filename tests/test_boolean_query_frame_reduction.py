"""Exact small native audits of the frame-to-Boolean-query reduction.

All truth tables through two address bits are checked on the full input
space. These tests establish no growing T-depth lower bound. The one
floating-point fixture tests an analytic hybrid with nonzero clean leakage.
"""
from fractions import Fraction
import unittest

import numpy as np

from compiler_robust_hopf.native_coarse_fixture import commutator_word, toffoli_word
from tests.exact_matrices import CO, CZ, CQ2, Q2, Gate, apply_word, native_columns
from tests.test_operator_source_compiler import _adjoint, _word_matrix


def _xor_word(table, address_bits, target):
    """Bounded ANF fixture; every Toffoli is expanded into native gates."""
    coefficients = list(table)
    for bit in range(address_bits):
        for mask in range(1 << address_bits):
            if mask & (1 << bit):
                coefficients[mask] ^= coefficients[mask ^ (1 << bit)]
    word = []
    for mask, coefficient in enumerate(coefficients):
        if not coefficient:
            continue
        controls = [bit for bit in range(address_bits) if mask & (1 << bit)]
        if not controls:
            word.append(('X', target))
        elif len(controls) == 1:
            word.append(('CX', controls[0], target))
        else:
            assert len(controls) == 2
            word.extend((name, *wires) for name, wires in
                        toffoli_word(*controls, target))
    return word


def _corner_word(table, address_bits, target):
    query = _xor_word(table, address_bits, target)
    # Chronological H,Q,H,Q has matrix Q H Q H = X^f Z^f.
    return [('H', target)] + query + [('H', target)] + query


def _query_from_frame(word, frame_target, query_target):
    cnot = [('CX', frame_target, query_target)]
    return _adjoint(word) + cnot + word + cnot


def _expected_query(table, address_bits, target, width):
    size = 1 << width
    expected = [[CZ for _ in range(size)] for _ in range(size)]
    for column in range(size):
        address = column & ((1 << address_bits) - 1)
        row = column ^ (table[address] << target)
        expected[row][column] = CO
    return expected


class BooleanQueryFrameReductionTests(unittest.TestCase):
    def test_every_small_table_native_full_input_and_converse(self):
        for address_bits in range(3):
            target, output, spectator = range(address_bits, address_bits + 3)
            width = spectator + 1
            for table_code in range(1 << (1 << address_bits)):
                table = [(table_code >> row) & 1
                         for row in range(1 << address_bits)]
                with self.subTest(address_bits=address_bits, table=table):
                    frame = _corner_word(table, address_bits, target)
                    native = native_columns(width, frame)
                    expected_frame = [[CZ for _ in range(1 << width)]
                                      for _ in range(1 << width)]
                    for column in range(1 << width):
                        f = table[column & ((1 << address_bits) - 1)]
                        row = column ^ (f << target)
                        sign = -1 if f and ((column >> target) & 1) else 1
                        expected_frame[row][column] = CQ2(sign)
                    self.assertEqual(native, expected_frame)
                    reduced = _query_from_frame(frame, target, output)
                    self.assertEqual(native_columns(width, reduced),
                                     _expected_query(table, address_bits, output, width))
                    # Reclassification uses output as logical, old target as dirty.
                    # Every old target/output/spectator basis input was checked.
                    self.assertTrue(all(gate[0] in
                                        {'X', 'CX', 'H', 'S', 'SDG', 'T', 'TDG'}
                                        for gate in reduced))

    def test_entangled_reference_with_arbitrary_borrowed_inputs(self):
        address_bits, target, output, width = 2, 2, 3, 5
        table = [1, 0, 1, 1]
        word = _query_from_frame(_corner_word(table, address_bits, target),
                                 target, output)
        # Maximally entangled with an untouched five-wire external reference.
        amplitude = CQ2(Q2(0, Fraction(1, 8)))  # 1/sqrt(32).
        state = {basis | (basis << width): amplitude for basis in range(1 << width)}
        native = [Gate(name, tuple(wires)) for name, *wires in word]
        actual = apply_word(native, state)
        expected = {}
        for basis in range(1 << width):
            image = basis ^ (table[basis & 3] << output)
            expected[image | (basis << width)] = amplitude
        self.assertEqual(actual, expected)

    def test_forward_substitution_leaves_the_forbidden_literal_phase(self):
        table, address_bits, target, output, width = [0, 1], 1, 1, 2, 3
        frame = _corner_word(table, address_bits, target)
        cnot = [('CX', target, output)]
        incorrect = native_columns(width, frame + cnot + frame + cnot)
        correct = _expected_query(table, address_bits, output, width)
        self.assertNotEqual(incorrect, correct)
        for column in range(1 << width):
            sign = -1 if column & 1 else 1
            self.assertEqual([row[column] for row in incorrect],
                             [sign * row[column] for row in correct])

    def test_actual_inverse_has_separate_elementary_clifford_charge(self):
        # S on initialized clean work preserves the frame contract, but its
        # literal inverse costs three S gates in the fixed H,S,CX alphabet.
        address, target, output, clean, width = 0, 1, 2, 3, 4
        frame = [('H', target), ('CX', address, target),
                 ('H', target), ('CX', address, target)]
        forward = frame + [('S', clean)]
        inverse = [('S', clean)] * 3 + list(reversed(frame))
        self.assertEqual(native_columns(width, inverse),
                         native_columns(width, _adjoint(forward)))
        self.assertEqual((len(forward), len(inverse)), (5, 7))
        self.assertLessEqual(len(inverse), 3 * len(forward))
        cnot = [('CX', target, output)]
        reduced = inverse + cnot + forward + cnot
        self.assertEqual(len(reduced), 14)
        self.assertGreater(len(reduced), 2 * len(forward) + 2)
        self.assertEqual(native_columns(width, reduced),
                         _expected_query([0, 1], 1, output, width))
        # This supplied-word audit includes every clean/dirty basis column;
        # it is not a claim that the deliberately padded source is minimal.

    def test_scalar_phase_and_explicit_controlled_table_conversions(self):
        address_bits, target, width = 2, 2, 4
        for table_code in range(16):
            table = [(table_code >> row) & 1 for row in range(4)]
            query = _xor_word(table, address_bits, target)
            scalar_word = [('Z', target)] + _adjoint(query) + [('Z', target)] + query
            scalar = native_columns(width, scalar_word)
            controlled_word = [('H', target)] + query + [('H', target)]
            controlled = native_columns(width, controlled_word)
            for column in range(1 << width):
                f = table[column & 3]
                for row in range(1 << width):
                    expected = CQ2(-1 if f else 1) if row == column else CZ
                    self.assertEqual(scalar[row][column], expected)
                    g = f * ((column >> target) & 1)
                    expected_g = CQ2(-1 if g else 1) if row == column else CZ
                    self.assertEqual(controlled[row][column], expected_g)
            recovered = [('H', target)] + controlled_word + [('H', target)]
            self.assertEqual(native_columns(width, recovered),
                             _expected_query(table, address_bits, target, width))

    def test_actual_inverse_hybrid_with_nonzero_native_clean_leakage(self):
        address_bits, target, output, clean, width = 1, 1, 2, 3, 4
        table = [0, 1]
        frame = _corner_word(table, address_bits, target)
        # A fixed small native commutator, conditionally conjugated through
        # a borrowed wire, leaks clean amplitude. It is not an exact-return V.
        perturbation = [('CX', output, clean)]
        perturbation += [(name, *wires) for name, wires in commutator_word(3, clean)]
        perturbation += [('CX', output, clean)]
        approximate_word = frame + perturbation
        v = _word_matrix(width, approximate_word)
        exact = _word_matrix(width, frame)
        initialized = np.eye(1 << width, dtype=complex)[:, :1 << clean]
        eta = np.linalg.norm(v @ initialized - exact @ initialized, 2)
        leakage = np.linalg.norm((v @ initialized)[1 << clean:], 2)
        self.assertGreater(leakage, 1e-8)
        self.assertLess(eta, 1 / 64)
        inverse_error = np.linalg.norm(v.conj().T @ initialized
                                       - exact.conj().T @ initialized, 2)
        self.assertLessEqual(inverse_error, eta + 2e-11)
        c = _word_matrix(width, [('CX', target, output)])
        actual = c @ v @ c @ v.conj().T
        desired = _word_matrix(width, _xor_word(table, address_bits, output))
        error = np.linalg.norm(actual @ initialized - desired @ initialized, 2)
        self.assertGreater(error, 1e-8)
        self.assertLessEqual(error, 2 * eta + 5e-11)
        # Literal native inverse, evaluated independently, agrees as a matrix.
        emitted = _word_matrix(width, _query_from_frame(approximate_word, target, output))
        np.testing.assert_allclose(emitted, actual, atol=5e-11, rtol=0)


if __name__ == '__main__':
    unittest.main()
