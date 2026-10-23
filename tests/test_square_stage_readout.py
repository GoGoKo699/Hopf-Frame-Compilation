"""A complete two-bit readout and the literal affine-Pauli boundary.

The p=2 commutator is a complete operator identity, not a growing-size
source.  Fixed Clifford mixing does not leave the tested canonical-Pauli
normal form.  Freely assigning Hermitian phases would leave that grammar.
"""
from itertools import product
import unittest

from tests.exact_matrices import CI, CO, CZ, adjoint, multiply, native_columns
from tests.test_co_moving_lookup_stage import _Matching
from tests.test_majorana_walsh_query import _inverse
from tests.test_t_depth import _ccz_batch


_PHASES = (CO, CI, -CO, -CI)


def _matched_label(p, x, bank, coordinate):
    return sum(1 << u for u in range(1 << p)
               if ((bank >> u) & 1) and not (((x ^ u) >> coordinate) & 1))


def _square_tag(state, matching):
    # A_i=Q R_C^2.  Both calls are the actual complete stage, not a
    # replacement by its desired diagonal formula.
    x, bank, selectors, borrowed, phase = matching.stage_literal(
        matching.stage_literal(state))
    phase += 2 * ((bank & selectors).bit_count() % 2)
    return x, bank, selectors, borrowed, phase % 4


def _reduce_s_runs(word):
    reduced = []
    position = 0
    while position < len(word):
        gate = word[position]
        if gate[0] != 'S':
            reduced.append(gate)
            position += 1
            continue
        end = position + 1
        while end < len(word) and word[end] == gate:
            end += 1
        reduced.extend([gate] * ((end - position) % 4))
        position = end
    return reduced


def _direct_tag_schedule(p, coordinate):
    """Literal four-layer A_i, with no square-stage implementation needed."""
    m = 1 << p
    bank = tuple(range(p, p + m))
    selectors = tuple(range(p + m, p + 2 * m))
    schedule = _ccz_batch(coordinate, list(zip(bank, selectors)))
    schedule[0] = ('C', _reduce_s_runs(schedule[0][1]))
    mask = []
    for u in range(m):
        if not ((u >> coordinate) & 1):
            mask.extend([('H', selectors[u]), ('CX', bank[u], selectors[u]),
                         ('H', selectors[u])])
    schedule.append(('C', mask))
    return schedule


def _mobius_mod(values, modulus):
    coefficients = list(values)
    p = (len(values) - 1).bit_length()
    for bit in range(p):
        for mask in range(len(values)):
            if mask & (1 << bit):
                coefficients[mask] = (
                    coefficients[mask] - coefficients[mask ^ (1 << bit)]) % modulus
    return coefficients


def _affine(x, mask, constant=0):
    return ((x & mask).bit_count() % 2) ^ constant


def _literal_pauli(a, b, width):
    # X^a Z^b acts chronologically as Z^b followed by X^a.
    word = []
    for wire in range(width):
        if (b >> wire) & 1:
            word.extend([('S', wire), ('S', wire)])
    for wire in range(width):
        if (a >> wire) & 1:
            word.extend([('H', wire), ('S', wire), ('S', wire), ('H', wire)])
    return word


def _identify_pauli(matrix):
    width = (len(matrix) - 1).bit_length()
    nonzero = [(row, matrix[row][0]) for row in range(len(matrix))
               if matrix[row][0] != CZ]
    if len(nonzero) != 1:
        raise AssertionError('Not a canonical Pauli column')
    a, phase = nonzero[0]
    if phase not in _PHASES:
        raise AssertionError('Literal scalar is not an integer power of i')
    b = 0
    for wire in range(width):
        amplitude = matrix[a ^ (1 << wire)][1 << wire]
        if amplitude == -phase:
            b |= 1 << wire
        elif amplitude != phase:
            raise AssertionError('Not a Pauli sign character')
    for column, row in product(range(len(matrix)), repeat=2):
        expected = (phase * (-1) ** ((b & column).bit_count())
                    if row == (column ^ a) else CZ)
        if matrix[row][column] != expected:
            raise AssertionError('Incorrect full Pauli matrix')
    return a, b, _PHASES.index(phase)


class SquareStageReadoutTests(unittest.TestCase):
    def test_direct_native_tag_and_actual_inverse(self):
        for p in (1, 2):
            m, width = 1 << p, p + 2 * (1 << p)
            selector_inputs = range(1 << m) if p == 1 else (0, 1, 5, 10, 15)
            # Every p=1 column; every p=2 address/bank with five distinct
            # selector columns.  No principal register is a clean helper.
            columns = [x | (bank << p) | (selector << (p + m))
                       for x, bank, selector in product(
                           range(m), range(1 << m), selector_inputs)]
            for coordinate in range(p):
                word = [gate for _, gates in _direct_tag_schedule(p, coordinate)
                        for gate in gates]
                reverse = _reduce_s_runs(_inverse(word))
                for emitted in (word, reverse):
                    actual = native_columns(width, emitted, columns)
                    for index, column in enumerate(columns):
                        x = column & (m - 1)
                        bank = (column >> p) & ((1 << m) - 1)
                        selector = column >> (p + m)
                        phase = (-1) ** ((_matched_label(p, x, bank, coordinate)
                                         & selector).bit_count())
                        for row in range(1 << width):
                            self.assertEqual(actual[row][index],
                                             CO * phase if row == column else CZ)

    def test_direct_native_schedule_and_complete_component_ledger(self):
        for p in range(1, 8):
            m = 1 << p
            schedule = _direct_tag_schedule(p, 0)
            self.assertEqual(sum(kind == 'T' for kind, _ in schedule), 4)
            for kind, gates in schedule:
                if kind == 'T':
                    targets = [gate[1] for gate in gates]
                    self.assertEqual(len(targets), len(set(targets)))
                    self.assertNotIn(0, targets)  # Shared address is Clifford-only.
            word = [gate for _, gates in schedule for gate in gates]
            for emitted in (word, _reduce_s_runs(_inverse(word))):
                t_count = sum(gate[0] in ('T', 'TDG') for gate in emitted)
                clifford_count = len(emitted) - t_count
                self.assertEqual(t_count, 6 * m)
                self.assertLessEqual(clifford_count, 23 * m // 2 + 3)
                # ABAB adds four full-selector H walls; the XOR indicator
                # adds two full-output H walls.  Every principal selector
                # returns in the exact logical commutator identity.
                self.assertLessEqual(4 * clifford_count + 4 * m, 50 * m + 12)
                self.assertLessEqual(4 * clifford_count + 6 * m, 52 * m + 12)

    def test_complete_two_bit_square_tags_and_commutator(self):
        p, m = 2, 4
        stages = (_Matching(p, ((0, 1),)), _Matching(p, ((1, 0),)))
        for x, bank, selectors, borrowed in product(
                range(m), range(1 << m), range(1 << m), range(2)):
            state = x, bank, selectors, borrowed, 0
            a, b = (_matched_label(p, x, bank, coordinate) for coordinate in (0, 1))
            for coordinate, matching in enumerate(stages):
                label = (a, b)[coordinate]
                self.assertEqual(_square_tag(state, matching),
                                 (x, bank, selectors, borrowed,
                                  2 * ((label & selectors).bit_count() % 2)))
            # B_j=H_all A_j H_all=X^b, checked by exact Walsh sums below.
            # Chronology B,A,B,A includes both actual square-tag calls.
            value = x, bank, selectors ^ b, borrowed, 0
            value = _square_tag(value, stages[0])
            value = value[0], value[1], value[2] ^ b, value[3], value[4]
            value = _square_tag(value, stages[0])
            self.assertEqual(value, (x, bank, selectors, borrowed,
                                     2 * ((bank >> x) & 1)))

    def test_full_walsh_conjugation_has_no_missing_scalar(self):
        # Four selector qubits give integer Walsh numerators with
        # denominator16 after the two Hadamard walls.
        dimension = 16
        for label, row, column in product(range(dimension), repeat=3):
            numerator = sum((-1) ** (((row ^ label ^ column) & middle).bit_count())
                            for middle in range(dimension))
            self.assertEqual(numerator, dimension if row == (column ^ label) else 0)

    def test_three_bit_extension_reads_a_subcube_not_one_row(self):
        p, m = 3, 8
        stage = _Matching(p, ((0, 1),))
        banks = (0, 1, 16, 17, 3, 129, 255)
        selectors = (0, 1, 16, 37, 170, 255)
        for x, bank, selector, borrowed in product(range(m), banks, selectors, range(2)):
            a, b = (_matched_label(p, x, bank, coordinate) for coordinate in (0, 1))
            state = x, bank, selector ^ b, borrowed, 0
            state = _square_tag(state, stage)
            state = state[0], state[1], state[2] ^ b, state[3], state[4]
            state = _square_tag(state, stage)
            subcube = sum((bank >> u) & 1 for u in range(m)
                          if ((u ^ x) & 3) == 0) % 2
            self.assertEqual(state, (x, bank, selector, borrowed, 2 * subcube))
        x, bank = 0, 16
        self.assertEqual((_matched_label(p, x, bank, 0)
                          & _matched_label(p, x, bank, 1)).bit_count() % 2, 1)
        self.assertEqual((bank >> x) & 1, 0)

    def test_literal_clifford_interlayers_retain_affine_labels_and_quadratic_phase(self):
        # Actual native matrices independently check the phase convention
        # behind the analytic normal-form proof.  Constants in the affine
        # labels are retained; selectors are never initialized.
        width, addresses = 2, 16
        cliffords = (
            [('H', 0), ('S', 1), ('CX', 0, 1)],
            [('S', 0), ('H', 1), ('CX', 1, 0)],
            [('H', 0), ('H', 1), ('S', 0), ('CX', 0, 1)],
            [('S', 1), ('H', 0), ('CX', 1, 0)],
        )
        forms = (
            ((3, 0), (10, 1), (5, 1), (12, 0)),
            ((15, 1), (6, 0), (9, 0), (7, 1)),
            ((11, 0), (13, 1), (14, 1), (1, 0)),
            ((4, 1), (2, 1), (8, 0), (15, 0)),
        )
        fixed_word = [gate for block in cliffords for gate in block]
        fixed_inverse = adjoint(native_columns(width, fixed_word))
        labels_a, labels_b, phases = [], [], []
        for x in range(addresses):
            word = []
            for block, coefficients in zip(cliffords, forms):
                word.extend(block)
                bits = [_affine(x, mask, constant) for mask, constant in coefficients]
                word.extend(_literal_pauli(bits[0] | (bits[1] << 1),
                                           bits[2] | (bits[3] << 1), width))
            relative = multiply(native_columns(width, word), fixed_inverse)
            a, b, phase = _identify_pauli(relative)
            labels_a.append(a)
            labels_b.append(b)
            phases.append(phase)
        for labels in (labels_a, labels_b):
            for wire in range(width):
                coefficients = _mobius_mod([(value >> wire) & 1 for value in labels], 2)
                self.assertTrue(all(coefficient == 0 for mask, coefficient
                                    in enumerate(coefficients) if mask.bit_count() > 1))
        coefficients = _mobius_mod(phases, 4)
        self.assertTrue(all(coefficient == 0 for mask, coefficient
                            in enumerate(coefficients) if mask.bit_count() > 2))
        self.assertTrue(all(coefficient % 2 == 0 for mask, coefficient
                            in enumerate(coefficients) if mask.bit_count() == 2))

    def test_freely_hermitian_normalized_controlled_paulis_escape_the_grammar(self):
        # P(a,b)=i^(ab)X^a Z^b has an extra address-dependent phase.
        # X^a Z^b P(a,b)=(-i)^(ab).  This is deliberately NOT asserted
        # to be a free literal affine-Pauli generator.
        signs = []
        for a, b, c in product(range(2), repeat=3):
            phase = (-a * (b ^ c) + a * b + a * c) % 4
            self.assertEqual(phase, 2 * a * b * c)
        for x in range(8):
            a, b, c = ((x >> bit) & 1 for bit in range(3))
            signs.append((-a * (b ^ c) + a * b + a * c) // 2)
        self.assertEqual(_mobius_mod(signs, 2)[7], 1)


if __name__ == '__main__':
    unittest.main()
