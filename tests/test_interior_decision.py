"""Native operator fixtures for the interior-Clifford decision calibration.

All source words use literal H/S/CX/T gates and their actual inverses.
The final compression is checked exactly in Q(sqrt(2), i) on every clean
input column. The all-size singlet and spectator arguments are analytic;
these five-qubit fixtures check conventions and the exact base identity.
"""

from fractions import Fraction
from math import comb
import unittest

import numpy as np

from tests.exact_matrices import CO, CQ2, CZ, Q2, adjoint, eye, multiply, native_columns
from tests.test_operator_source_compiler import (
    ATOL, X, Y, Z, _adjoint, _pauli, _word_matrix,
)


# Least-significant wire first: clean c, A0, A1, B0, B1.
WIDTH = 5


def _z_word(qubit):
    return [('S', qubit), ('S', qubit)]


def _xx_pulse(bank):
    # exp(-i pi X_c X_bank/8), with the literal T scalar retained.
    return [('CX', 0, bank), ('H', 0), ('T', 0),
            ('H', 0), ('CX', 0, bank)]


def _yx_pulse(bank):
    return [('SDG', 0)] + _xx_pulse(bank) + [('S', 0)]


def _pulses():
    return _xx_pulse(1) + _xx_pulse(2), _yx_pulse(3) + _yx_pulse(4)


def _graph_word(table):
    return [gate for i in range(2) for j in range(2) if table[i][j]
            for gate in [('H', 3+j), ('CX', 1+i, 3+j), ('H', 3+j)]]


def _shear_word(table, reverse=False):
    return [('CX', 3+j, 1+i) if reverse else ('CX', 1+i, 3+j)
            for i in range(2) for j in range(2) if table[i][j]]


def _symmetry_breaker():
    # Chronological U0 followed by U1; SWAP and Z are expanded natively.
    return ([('SDG', 0), ('H', 0), ('H', 3), ('CX', 0, 3),
             ('CX', 0, 1), ('CX', 1, 0), ('CX', 0, 1),
             ('H', 3), ('H', 1), ('CX', 3, 4), ('CX', 1, 2)]
            + _z_word(3) + _z_word(1))


def _boundary_word(interior):
    pulse_a, pulse_b = _pulses()
    return pulse_b + interior + pulse_a


def _middle_word(table, u=0, v=0):
    pulse_a, pulse_b = _pulses()
    gauge_a = _z_word(2) if u else []
    gauge_b = _z_word(4) if v else []
    addressed_a = gauge_a + pulse_a + _adjoint(gauge_a)
    addressed_b = gauge_b + pulse_b + _adjoint(gauge_b)
    graph = _graph_word(table)
    k0 = (addressed_b + addressed_a
          + _adjoint(addressed_b) + _adjoint(addressed_a))
    kc = (graph + addressed_b + graph + addressed_a
          + graph + _adjoint(addressed_b) + graph + _adjoint(addressed_a))
    return _adjoint(k0) + kc


def _source_word(table, interior, u=0, v=0):
    boundary = _boundary_word(interior)
    return boundary + _middle_word(table, u, v) + _adjoint(boundary)


def _asymmetric_source_word(table, u=0, v=0):
    """Literal calibrated word with independent interior code boundaries.

    Operator order is A0† E B0† (KC K0†) A0 E A0† B0 A0.
    The cancelling B0 B0† pair in Fin,C Fin,0† has been removed.
    """
    pulse_a, pulse_b = _pulses()
    graph = _graph_word(table)
    return (pulse_a + pulse_b + _adjoint(pulse_a) + graph + pulse_a
            + _middle_word(table, u, v) + _adjoint(pulse_b)
            + graph + _adjoint(pulse_a))


def _singlet_subspace(pair):
    """Isometry with the selected pair singlet and every other wire arbitrary."""
    free = [qubit for qubit in range(WIDTH) if qubit not in pair]
    result = np.zeros((1 << WIDTH, 1 << len(free)), dtype=complex)
    for column in range(1 << len(free)):
        base = sum(((column >> index) & 1) << qubit
                   for index, qubit in enumerate(free))
        result[base | (1 << pair[0]), column] = 1 / np.sqrt(2)
        result[base | (1 << pair[1]), column] = -1 / np.sqrt(2)
    return result


class InteriorDecisionTests(unittest.TestCase):
    def test_asymmetric_boundaries_have_complete_repeated_row_response(self):
        # The zero-table calibration is literal on all physical inputs.
        zero_word = _asymmetric_source_word(((0, 0), (0, 0)))
        self.assertEqual(native_columns(WIDTH, zero_word), eye(1 << WIDTH))

        root_half = CQ2(Q2(Fraction(0), Fraction(1, 2)))
        selected = sorted(signal | (bank << 3) | a_bit
                          for signal in (0, 1) for bank in range(4)
                          for a_bit in (2, 4))
        position = {basis: index for index, basis in enumerate(selected)}

        def exact_response(signal, bank, mask, inverse):
            # Independent closed form: product_j (I - i Y_c X_Bj)/sqrt(2)
            # on the unnormalized A-pair singlet, with actual inverse signs.
            base = signal | (bank << 3)
            vector = {base | 2: CO, base | 4: -CO}
            for j in range(2):
                if not (mask >> j) & 1:
                    continue
                following = {}
                for basis, value in vector.items():
                    following[basis] = following.get(basis, CZ) + value * root_half
                    changed = basis ^ 1 ^ (1 << (3+j))
                    sign = (-1 if basis & 1 else 1) * (-1 if inverse else 1)
                    following[changed] = (following.get(changed, CZ)
                                          + value * root_half * CQ2(sign))
                vector = following
            return [vector.get(row, CZ) for row in range(1 << WIDTH)]

        # Odd z gives zero bias; even z gives a nontrivial signed Pauli.
        # Both address values and every signal/B input are retained.
        for row in ((1, 0), (1, 1)):
            table = (row, row)
            mask = row[0] + 2*row[1]
            for v in (0, 1):
                word = _asymmetric_source_word(table, v=v)
                self.assertEqual(sum(gate[0] in ('T', 'TDG') for gate in word), 28)
                for inverse in (False, True):
                    with self.subTest(row=row, v=v, inverse=inverse):
                        actual_word = _adjoint(word) if inverse else word
                        columns = native_columns(WIDTH, actual_word, columns=selected)
                        outputs = {}
                        for signal in (0, 1):
                            for bank in range(4):
                                base = signal | (bank << 3)
                                output = [values[position[base | 2]]
                                          - values[position[base | 4]]
                                          for values in columns]
                                self.assertEqual(output,
                                                 exact_response(signal, bank, mask, inverse))
                                outputs[signal, bank] = output
                        if inverse:
                            continue
                        for left in range(4):
                            for right in range(4):
                                bias = sum((outputs[0, left][basis].conj()
                                            * CQ2(-1 if basis & 1 else 1)
                                            * outputs[0, right][basis]
                                            for basis in range(1 << WIDTH)), CZ)
                                # The unnormalized singlet has squared norm two.
                                bias *= CQ2(Fraction(1, 2))
                                expected = (CZ if mask.bit_count() % 2
                                            else CQ2((-1)**(mask.bit_count()//2)
                                                     * (left == (right ^ mask))))
                                self.assertEqual(bias, expected)

    def test_native_decoder_retains_selected_row_and_omits_row_zero(self):
        pulse_a, pulse_b = _pulses()
        signal = _pauli(WIDTH, {0: Z})
        b_matrix = _word_matrix(WIDTH, pulse_b)
        fixed_signal = b_matrix @ signal @ b_matrix.conj().T
        tables = (((0, 1), (0, 1)),   # The selected nonzero row is (0,1).
                  ((1, 0), (0, 1)),   # Only row zero changes.
                  ((0, 1), (1, 0)))   # The selected row changes.
        signals = {}
        for table_index, table in enumerate(tables):
            graph = _graph_word(table)
            for v in (0, 1):
                gauge = _z_word(2) + (_z_word(4) if v else [])
                # Gamma=E-dagger A-dagger G A E, with literal conjugations.
                gamma_word = graph + pulse_a + gauge + _adjoint(graph + pulse_a)
                gamma = _word_matrix(WIDTH, gamma_word)
                transformed = gamma.conj().T @ fixed_signal @ gamma
                signals[table_index, v] = transformed
                if v:
                    continue
                axis = _pauli(WIDTH, {0: X, 2: X,
                                     **{3+j: Z for j in range(2) if table[1][j]}})
                for j in range(2):
                    with self.subTest(table=table_index, column=j):
                        bank_z = _pauli(WIDTH, {3+j: Z})
                        bank_axis = _pauli(WIDTH, {0: Y, 3+j: X})
                        selected = 1 - table[1][j]
                        product_axis = axis if selected else np.eye(1 << WIDTH)
                        expected = (-1j * 1j**selected
                                    * product_axis @ bank_axis @ bank_z)
                        np.testing.assert_allclose(transformed @ bank_z @ transformed,
                                                   expected, atol=ATOL)
        for v in (0, 1):
            np.testing.assert_allclose(signals[0, v], signals[1, v], atol=ATOL)
        self.assertGreater(np.linalg.norm(signals[0, 0] - signals[2, 0], 2), 1)

    def test_one_bank_singlet_is_fixed_for_graph_and_both_shears(self):
        # Rows repeat, but columns do not. B and c remain completely arbitrary.
        # At u=v=0 the selected table entry is one, so identity action has
        # the wrong decision sign. The same dark space persists at v=1.
        table = ((1, 0), (1, 0))
        dark = _singlet_subspace((1, 2))
        interiors = (_graph_word(table), _shear_word(table),
                     _shear_word(table, reverse=True))
        for index, interior in enumerate(interiors):
            for v in (0, 1):
                with self.subTest(interior=index, v=v):
                    source = _word_matrix(WIDTH, _source_word(table, interior, v=v))
                    np.testing.assert_allclose(source @ dark, dark, atol=ATOL)

    def test_symmetry_breaker_avoids_each_known_singlet_sector(self):
        table = ((1, 1), (1, 1))
        interior = _symmetry_breaker() + _graph_word(table)
        clifford = _word_matrix(WIDTH, _symmetry_breaker())
        for original, image in (
            (_pauli(WIDTH, {0: Y}), -_pauli(WIDTH, {1: X, 2: X})),
            (_pauli(WIDTH, {0: Y, 3: X}), -_pauli(WIDTH, {3: X, 4: X})),
        ):
            np.testing.assert_allclose(clifford @ original @ clifford.conj().T,
                                       image, atol=ATOL)
        boundary = _word_matrix(WIDTH, _boundary_word(interior))
        code = boundary[:, ::2]
        code_projector = code @ code.conj().T
        for pair in ((1, 2), (3, 4)):
            dark = _singlet_subspace(pair)
            # Half overlap on the full dark sector rules out any vector
            # shared with the code; this does not assert a decision gap.
            np.testing.assert_allclose(dark.conj().T @ code_projector @ dark,
                                       np.eye(8) / 2, atol=ATOL)

    def test_exact_native_repair_compression_on_all_dirty_columns(self):
        table = ((1, 1), (1, 1))
        interior = _symmetry_breaker() + _graph_word(table)
        word = _source_word(table, interior)
        self.assertTrue(all(gate[0] in {'H', 'S', 'SDG', 'CX', 'T', 'TDG'}
                            for gate in word))
        self.assertEqual(sum(gate[0] in ('T', 'TDG') for gate in word), 24)
        # Every dirty input column is retained, including rejected-signal
        # amplitudes. The word contains actual native inverse passages.
        columns = native_columns(WIDTH, word, columns=range(0, 32, 2))
        for row in range(16):
            for column in range(16):
                actual = sum((columns[basis][row].conj()
                              * CQ2(-1 if basis & 1 else 1)
                              * columns[basis][column]
                              for basis in range(32)), CZ)
                expected = CQ2(Fraction(3, 4) * (row == column)
                               - Fraction(1, 4) * (row == (column ^ 4)))
                self.assertEqual(actual, expected, (row, column))
        # This is exactly M=(3I-X_B0)/4, not just its diagonal. It has
        # positive eigenvalues 1/2 and 1 although the selected entry is one.

    def test_all_ones_middle_has_a_majority_short_spectral_arc(self):
        # Proof home: PRESCRIBED_SCALAR_TABLES.md,
        # #boundary-independent-spectral-obstruction. For C=J_m the
        # analytic commuting-axis identity gives K_C=I, so Q=K0-dagger.
        # Evaluate that literal native word on every physical column.
        pulse_a, pulse_b = _pulses()
        k0 = pulse_b + pulse_a + _adjoint(pulse_b) + _adjoint(pulse_a)
        middle = native_columns(WIDTH, _adjoint(k0))

        # Verify translation symmetry before extracting the bank-X blocks;
        # this does not assume the desired spectral multiplicities.
        for bank in range(16):
            for signal in (0, 1):
                for row in range(32):
                    self.assertEqual(middle[row][2 * bank + signal],
                                     middle[row ^ (bank << 1)][signal])
        identity = eye(2)
        identity_multiplicity = positive_arc_multiplicity = negative_arc_multiplicity = 0
        for character in range(16):
            block = [[sum((middle[2 * bank + output][signal]
                           * (-1 if (character & bank).bit_count() % 2 else 1)
                           for bank in range(16)), CZ)
                      for signal in (0, 1)] for output in (0, 1)]
            self.assertEqual(multiply(adjoint(block), block), identity)
            if block == identity:
                identity_multiplicity += 2
            else:
                # Q^2-Q+I=0 and trace one give one root exp(+i*pi/3)
                # and one exp(-i*pi/3), exactly in this two-dimensional block.
                self.assertEqual(multiply(block, block),
                                 [[block[row][column] - identity[row][column]
                                   for column in (0, 1)] for row in (0, 1)])
                self.assertEqual(block[0][0] + block[1][1], CO)
                positive_arc_multiplicity += 1
                negative_arc_multiplicity += 1
        self.assertEqual((identity_multiplicity, positive_arc_multiplicity,
                          negative_arc_multiplicity), (24, 4, 4))
        self.assertGreater(identity_multiplicity + positive_arc_multiplicity, 32 // 2)

        # Independent integer binomial sums check the residue probabilities
        # at three small sizes. Their all-size derivation is in the proof.
        for m in (2, 4, 8):
            with self.subTest(bank_size=m):
                counts = [sum(comb(m, k) for k in range(m + 1)
                              if (m // 2 - k) % 4 == residue)
                          for residue in range(4)]
                q0 = Fraction(counts[0], 1 << m)
                qodd = Fraction(counts[1] + counts[3], 1 << m)
                self.assertEqual(q0, Fraction(1, 4) + Fraction(1, 1 << (m // 2 + 1)))
                self.assertEqual(qodd, Fraction(1, 2))
                self.assertGreater(2 * q0 - q0 * q0 + qodd * qodd / 2,
                                   Fraction(1, 2))


if __name__ == '__main__':
    unittest.main()
