"""Exact native controls for disjoint nonlinear selector phases.

The borrowed bits and control are arbitrary. Small fixtures retain every
literal phase and actual inverse; larger fixtures audit the native layers
and simultaneous counts. This pays a selector-mixing component, not a
complete accepted lookup source or an asymptotic depth improvement.
"""
from fractions import Fraction
import unittest

from tests.exact_matrices import CO, CQ2, Q2, adjoint, native_columns, zeros
from tests.test_fixed_row_transport_readout import _positive_s_schedule
from tests.test_phase_signal_query_limits import _counts
from tests.test_t_depth import _ccz_batch, _inverse, _word


def _nonlinear_selector_control(q, fourier=False, extra_selectors=0):
    """Control q disjoint CCZs, returning q arbitrary borrowed bits.

    Wire zero is the signal, the next 3q wires are selector triples,
    extra_selectors further wires are spectators, and the last q wires
    are borrowed. Fourier conjugation acts on the whole selector bank.
    """
    assert q >= 1 and extra_selectors >= 0
    m = 3*q + extra_selectors
    triples = [(1+3*j, 2+3*j, 3+3*j) for j in range(q)]
    borrowed = [1+m+j for j in range(q)]
    # These singleton CCZ batches have disjoint complete wire supports,
    # so matching native layers can be combined without shared T targets.
    individual = [_ccz_batch(r, [(b, c)])
                  for (_, b, c), r in zip(triples, borrowed)]
    phase = [(individual[0][i][0],
              [gate for schedule in individual for gate in schedule[i][1]])
             for i in range(len(individual[0]))]
    walls = [('H', r) for r in borrowed]
    add = ([('C', walls)]
           + _ccz_batch(0, [(a, r) for (a, _, _), r in zip(triples, borrowed)])
           + [('C', walls)])
    # r bc + (r xor signal*a)bc = signal*abc, with both additions undone.
    schedule = phase + add + phase + add
    if fourier:
        selector_walls = [('H', wire) for wire in range(1, m+1)]
        schedule = [('C', selector_walls)] + schedule + [('C', selector_walls)]
    return 1+m+q, _positive_s_schedule(schedule)


def _controlled_phase_columns(q, columns=None):
    width = 1+4*q
    columns = list(range(1 << width)) if columns is None else list(columns)
    expected = zeros(1 << width, len(columns))
    for j, column in enumerate(columns):
        parity = sum(((column >> (1+3*t)) & 7) == 7 for t in range(q)) % 2
        expected[column][j] = -CO if column & 1 and parity else CO
    return expected


def _controlled_fourier_columns(q, columns=None):
    """Independent kernel: H^3 CCZ H^3 = I - vv^T/4, v_s=(-1)^|s|."""
    width, m = 1+4*q, 3*q
    columns = list(range(1 << width)) if columns is None else list(columns)
    expected = zeros(1 << width, len(columns))
    mask = (1 << m)-1
    for j, column in enumerate(columns):
        if not column & 1:
            expected[column][j] = CO
            continue
        incoming = (column >> 1) & mask
        fixed = column & ~(mask << 1)
        for outgoing in range(1 << m):
            coefficient = Fraction(1)
            for triple in range(q):
                left = (incoming >> (3*triple)) & 7
                right = (outgoing >> (3*triple)) & 7
                sign = -1 if (left ^ right).bit_count() % 2 else 1
                coefficient *= int(left == right) - Fraction(sign, 4)
            expected[fixed | (outgoing << 1)][j] = CQ2(Q2(coefficient))
    return expected


class NonlinearSelectorControlTests(unittest.TestCase):
    def test_one_triple_all_columns_literal_phase_and_actual_inverse(self):
        for fourier in (False, True):
            width, schedule = _nonlinear_selector_control(1, fourier)
            expected = (_controlled_fourier_columns(1) if fourier
                        else _controlled_phase_columns(1))
            actual = native_columns(width, _word(schedule))
            self.assertEqual(actual, expected)
            reverse = _positive_s_schedule(_inverse(schedule))
            self.assertEqual(native_columns(width, _word(reverse)), adjoint(actual))

    def test_two_triples_shared_signal_and_every_borrowed_assignment(self):
        # The phase fixture includes both signal values and every borrowed
        # assignment, with zero, one, or two active selector triples.
        selectors = (0, 7, 0o70, 0o77, 0o25, 0o52)
        columns = [(dirty << 7) | (bank << 1) | signal
                   for dirty in range(4) for bank in selectors for signal in (0, 1)]
        width, schedule = _nonlinear_selector_control(2)
        expected = _controlled_phase_columns(2, columns)
        for orientation in (schedule, _positive_s_schedule(_inverse(schedule))):
            self.assertEqual(native_columns(width, _word(orientation), columns), expected)

        # The dense Fourier columns also retain all four dirty assignments.
        columns = [(dirty << 7) | (bank << 1) | signal
                   for dirty in range(4)
                   for bank, signal in ((0, 1), (0o77, 1), (0o25, 0))]
        width, schedule = _nonlinear_selector_control(2, fourier=True)
        expected = _controlled_fourier_columns(2, columns)
        for orientation in (schedule, _positive_s_schedule(_inverse(schedule))):
            self.assertEqual(native_columns(width, _word(orientation), columns), expected)

    def test_native_layers_and_positive_s_inverse_pay_complete_counts(self):
        for q in (1, 2, 3, 4, 8, 16, 32, 64):
            for fourier in (False, True):
                # One unused selector per triple gives the four-row-cell
                # layout m=4q used by the source scout. Its H walls are paid.
                width, schedule = _nonlinear_selector_control(q, fourier,
                                                               extra_selectors=q)
                for orientation in (schedule, _positive_s_schedule(_inverse(schedule))):
                    self.assertEqual(sum(kind == 'T' for kind, _ in orientation), 16)
                    for kind, gates in orientation:
                        self.assertTrue(all(0 <= wire < width
                                            for gate in gates for wire in gate[1:]))
                        if kind == 'T':
                            self.assertTrue(gates)
                            self.assertTrue(all(gate[0] in ('T', 'TDG') for gate in gates))
                            targets = [gate[1] for gate in gates]
                            self.assertEqual(len(targets), len(set(targets)))
                        else:
                            self.assertTrue(all(gate[0] in ('H', 'S', 'CX')
                                                for gate in gates))
                    t_count, clifford_count = _counts(_word(orientation))
                    self.assertEqual(t_count, 26*q + 2*(q % 2))
                    self.assertLessEqual(clifford_count, 44*q+6 + (8*q if fourier else 0))


if __name__ == '__main__':
    unittest.main()
