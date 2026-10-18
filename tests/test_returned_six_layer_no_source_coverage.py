"""Exact native check: full return does not force source-wide C3 coverage.

All arithmetic is in Q(sqrt(2), i).
No native word is simplified before extracting the specified source.
"""
from fractions import Fraction
import unittest

from tests.exact_matrices import (Q2, CZ, CO, CI, add, adjoint, eye,
                                  multiply, native_columns, scale)


I = eye(2)
X = [[CZ, CO], [CO, CZ]]
Y = [[CZ, -CI], [CI, CZ]]
Z = [[CO, CZ], [CZ, -CO]]
PAULIS = [X, Y, Z]
H = native_columns(1, [("H", 0)])
T = native_columns(1, [("T", 0)])
TDG = adjoint(T)


def product(*factors):
    result = eye(len(factors[0]))
    for factor in factors:
        result = multiply(result, factor)
    return result


def conjugate(unitary, operator):
    return product(unitary, operator, adjoint(unitary))


def pauli_sum(coefficients):
    result = [[CZ, CZ], [CZ, CZ]]
    for coefficient, pauli in zip(coefficients, PAULIS):
        result = add(result, scale(pauli, coefficient))
    return result


def signed_pauli(operator):
    return any(operator == scale(pauli, sign)
               for pauli in [I] + PAULIS for sign in (-1, 1))


def clifford(operator):
    return all(signed_pauli(conjugate(operator, pauli)) for pauli in PAULIS)


# The arrays are indexed in the six-layer paper's order, not application order.
D = [None, TDG, TDG, T, T, T, TDG]
C = [H, H, I, H, H, I, I]
# Chronological order C0,D1,C1,...,D6,C6, retaining all six nonempty bands.
WORD = [("H", 0), ("TDG", 0), ("H", 0), ("TDG", 0),
        ("T", 0), ("H", 0), ("T", 0), ("H", 0),
        ("T", 0), ("TDG", 0)]
INVERSE_GATE = {"H": "H", "T": "TDG", "TDG": "T"}
INVERSE_WORD = [(INVERSE_GATE[name], qubit) for name, qubit in reversed(WORD)]
SOURCE_WORD = [("T", 0), ("H", 0), ("T", 0), ("H", 0), ("T", 0)]
SOURCE_INVERSE = [(INVERSE_GATE[name], qubit)
                  for name, qubit in reversed(SOURCE_WORD)]


class ReturnedIdentityNoSourceCoverage(unittest.TestCase):
    def test_literal_full_word_all_columns_and_actual_inverse(self):
        actual_from_slots = C[0]
        for index in range(1, 7):
            actual_from_slots = product(C[index], D[index], actual_from_slots)
        self.assertEqual(actual_from_slots, I)
        self.assertEqual(native_columns(1, WORD), actual_from_slots)
        self.assertEqual(native_columns(1, INVERSE_WORD), adjoint(actual_from_slots))
        # Two additional arbitrary untouched qubits can serve as address/reference.
        self.assertEqual(native_columns(3, WORD), eye(8))
        self.assertEqual(native_columns(3, INVERSE_WORD), eye(8))
        for word in [WORD, INVERSE_WORD]:
            self.assertEqual(sum(name in ("T", "TDG") for name, _ in word), 6)
            self.assertEqual(sum(name == "H" for name, _ in word), 4)
        self.assertEqual(len(D) - 1, 6)
        self.assertTrue(all(band in (T, TDG) for band in D[1:]))

    def test_actual_normal_form_and_literal_source_inverse(self):
        source = product(C[5], D[5], C[4], D[4], C[3], D[3])
        self.assertEqual(source, native_columns(1, SOURCE_WORD))
        self.assertEqual(adjoint(source), native_columns(1, SOURCE_INVERSE))
        self.assertEqual(product(source, adjoint(source)), I)
        for word in [SOURCE_WORD, SOURCE_INVERSE]:
            self.assertEqual(sum(name in ("T", "TDG") for name, _ in word), 3)
            self.assertEqual(sum(name == "H" for name, _ in word), 2)
        # All address-normal-form P_j and lambda_j are literally identity/one.
        s = product(C[3], D[3])
        e_in = product(C[2], D[2])
        inner = product(adjoint(s), I, s, I, e_in, I, adjoint(e_in))
        w6 = product(C[6], D[6])
        outer = conjugate(product(C[5], D[5]), I)
        target_comparison = product(I, adjoint(w6), I, I, I, I, w6)
        comparison = product(adjoint(outer), target_comparison)
        self.assertEqual(inner, I)
        self.assertEqual(comparison, I)
        self.assertEqual(conjugate(source, inner), comparison)

    def test_exact_domain_zero_and_all_inverse_images_nonclifford(self):
        source = native_columns(1, SOURCE_WORD)
        alpha = Q2(Fraction(1, 2), Fraction(-1, 4))
        beta = Q2(Fraction(1, 2), Fraction(1, 4))
        half = Q2(Fraction(1, 2))
        a = Q2(0, Fraction(1, 2))
        rows = [(alpha, -beta, half), (beta, -alpha, -half), (half, half, a)]
        for index, pauli in enumerate(PAULIS):
            forward = conjugate(source, pauli)
            backward = conjugate(adjoint(source), pauli)
            self.assertEqual(forward, pauli_sum([row[index] for row in rows]))
            self.assertEqual(backward, pauli_sum(rows[index]))
            self.assertFalse(signed_pauli(forward))
            self.assertFalse(signed_pauli(backward))
            self.assertFalse(clifford(forward))
            self.assertFalse(clifford(backward))
        # These three witnesses directly disprove inverse Cliffordness.
        expected = [(alpha, -beta, -half), (-beta, alpha, -half), (-half, half, a)]
        for source_probe, witness, coefficients in zip(PAULIS, [Z, Z, X], expected):
            image = conjugate(adjoint(source), source_probe)
            second_image = conjugate(image, witness)
            self.assertEqual(second_image, pauli_sum(coefficients))
            self.assertFalse(signed_pauli(second_image))


if __name__ == "__main__":
    unittest.main(verbosity=2)
