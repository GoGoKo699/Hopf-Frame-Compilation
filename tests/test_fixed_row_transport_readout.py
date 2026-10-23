"""Exact fixed-row readout and the cost of translating its distinguished row.

The complete p=2 logical fixtures retain address, bank, selectors and a
borrowed parity bit; stage calls use the retained literal implementation
and its actual reverse.  Larger cases are structured checks, not a full
native stage emitter.  Separate native transport fixtures use at most six
wires and keep every column and literal phase.  Translation and translated
zeta identities do not establish sublogarithmic dirty lookup.
"""
from itertools import product
import unittest

from tests.exact_matrices import CO, adjoint, native_columns, zeros
from tests.test_co_moving_lookup_stage import _Matching
from tests.test_entangling_selector_readout import _binary_image, _zeta_word
from tests.test_phase_signal_query_limits import _counts
from tests.test_t_depth import _fredkin_batch, _inverse, _word
from tests.test_two_layer_diagonal_rigidity import _reduce_adjacent_s


def _reverse_logical(word):
    return [('R', not inverse) if name == 'R' else (name, inverse)
            for name, inverse in reversed(word)]


def _commutator_word(word, pauli):
    """Chronology for [A,P]=A P A^dagger P, with P a Hermitian Pauli."""
    return [pauli] + _reverse_logical(word) + [pauli] + word


def _readout_words(row):
    echo = _commutator_word([('R', False)], ('X', 0))
    pair = _commutator_word(echo, ('X', row))
    central = _commutator_word(pair, ('Z', row))
    return echo, pair, central


def _apply_logical(matching, word, state):
    for name, argument in word:
        if name == 'R':
            state = matching.stage_literal(state, inverse=argument)
        else:
            x, bank, selectors, borrowed, phase = state
            if name == 'X':
                selectors ^= 1 << argument
            else:
                assert name == 'Z'
                phase = (phase + 2 * ((selectors >> argument) & 1)) % 4
            state = x, bank, selectors, borrowed, phase
    return state


def _echo_formula(matching, state):
    x, bank, selectors, borrowed, phase = state
    marked_zero = bank & 1
    distance_zero = sum((x >> control) & 1 for control in matching.controls)
    phase = (phase + marked_zero * distance_zero * (2 * (selectors & 1) - 1)) % 4
    if marked_zero:
        x = matching.images[x]
        bank = matching.relabel(bank)
        selectors = matching.relabel(selectors)
    return x, bank, selectors, borrowed, phase


def _translate(bank, address, p):
    """Physical row translation T_x: the value in row u moves to u xor x."""
    return sum(((bank >> u) & 1) << (u ^ address) for u in range(1 << p))


def _positive_s_schedule(schedule):
    """Expand actual S adjoints and reduce only adjacent positive S powers."""
    result = []
    for kind, gates in schedule:
        expanded = []
        for name, *wires in gates:
            if name == 'SDG':
                expanded.extend([('S', *wires)] * 3)
            else:
                expanded.append((name, *wires))
        result.append((kind, _reduce_adjacent_s(expanded)))
    return result


def _translation_schedule(p, bank_count=1):
    """Emit the address-controlled butterfly, with every bank arbitrary."""
    m = 1 << p
    schedule = []
    for coordinate in range(p):
        pairs = [(p + bank * m + u, p + bank * m + (u ^ (1 << coordinate)))
                 for bank in range(bank_count) for u in range(m)
                 if not u & (1 << coordinate)]
        schedule.extend(_fredkin_batch(coordinate, pairs))
    return _positive_s_schedule(schedule)


class FixedRowTransportReadoutTests(unittest.TestCase):
    def assert_readout(self, matching, state):
        row = matching.m - 1
        moved = matching.images[row]
        self.assertNotEqual(row, moved)
        echo, pair, central = _readout_words(row)
        x, bank, selectors, borrowed, phase = state
        echoed = _apply_logical(matching, echo, state)
        self.assertEqual(echoed, _echo_formula(matching, state))
        self.assertEqual(_apply_logical(matching, _reverse_logical(echo), echoed), state)
        paired = _apply_logical(matching, pair, state)
        pair_mask = ((1 << row) | (1 << moved)) if bank & 1 else 0
        self.assertEqual(paired, (x, bank, selectors ^ pair_mask, borrowed, phase))
        result = _apply_logical(matching, central, state)
        self.assertEqual(result, (x, bank, selectors, borrowed, (phase + 2*(bank & 1)) % 4))
        self.assertEqual(_apply_logical(matching, _reverse_logical(central), result), state)

    def test_two_bit_every_address_bank_selector_and_borrowed_column(self):
        for pairs in (((0, 1),), ((1, 0),)):
            matching = _Matching(2, pairs)
            for x, bank, selectors, borrowed in product(range(4), range(16),
                                                       range(16), range(2)):
                self.assert_readout(matching, (x, bank, selectors, borrowed, 0))

    def test_structured_three_and_four_bit_complete_return_and_nonzero_phases(self):
        families = ((3, ((0, 1),)), (3, ((2, 0),)),
                    (4, ((0, 1), (2, 3))), (4, ((0, 2), (1, 3))))
        for p, pairs in families:
            matching = _Matching(p, pairs)
            m, full = matching.m, matching.mask
            even = sum(1 << u for u in range(m) if not u & 1)
            banks = {0, full, 1, full ^ 1, 3, 1 << (m-1),
                     1 | (1 << (m-1)), even, full ^ even}
            selectors = {0, full, 1, 1 << (m-1), even, full ^ even,
                         (1 << (m//2)) | 1}
            for x, bank, selector, borrowed in product(range(m), sorted(banks),
                                                      sorted(selectors), range(2)):
                self.assert_readout(matching, (x, bank, selector, borrowed, 3))

    def test_nonzero_address_distinguishes_fixed_row_word_from_lookup(self):
        for p in (2, 3, 4):
            matching = _Matching(p, ((0, 1),))
            central = _readout_words(matching.m-1)[2]
            for x in range(1, matching.m):
                for bank in (1, 1 << x):
                    state = x, bank, matching.mask, 1, 0
                    actual = _apply_logical(matching, central, state)
                    lookup_phase = 2 * ((bank >> x) & 1)
                    self.assertEqual((actual[4] - lookup_phase) % 4, 2)
                    self.assertEqual(actual[:4], state[:4])

    def test_nested_word_counts_actual_calls_and_pauli_overhead(self):
        echo, pair, central = _readout_words(3)
        for word, calls, paulis in ((echo, 2, 2), (pair, 4, 6), (central, 8, 14)):
            for orientation in (word, _reverse_logical(word)):
                self.assertEqual(sum(name == 'R' for name, _ in orientation), calls)
                self.assertEqual(sum(name != 'R' for name, _ in orientation), paulis)
                self.assertEqual(sum(name == 'R' and inverse
                                     for name, inverse in orientation), calls // 2)
        # The 928-layer majorant charges all eight actual stage calls;
        # fixed Pauli interlayers contribute no T layer.
        self.assertEqual(sum(name == 'R' for name, _ in central) * 116, 928)
        # In the literal H,S,CX alphabet, X=HSSH and Z=SS: 14 Pauli
        # symbols require the safe bound 52 elementary Clifford gates.
        self.assertEqual(sum(4 if name == 'X' else 2 if name == 'Z' else 0
                             for name, _ in central), 52)

    def test_translation_conjugates_fixed_row_phase_to_complete_lookup(self):
        for p in range(1, 6):
            m, full = 1 << p, (1 << (1 << p)) - 1
            banks = (range(1 << m) if p <= 3 else
                     [0, full, *[1 << u for u in range(m)],
                      *[full ^ (1 << u) for u in range(m)]])
            for x, bank in product(range(m), banks):
                transported = _translate(bank, x, p)
                self.assertEqual(_translate(transported, x, p), bank)
                self.assertEqual(transported & 1, (bank >> x) & 1)

    def test_translated_zeta_flip_is_an_address_selected_selector_flip(self):
        for p in range(1, 7):
            m, full = 1 << p, (1 << (1 << p)) - 1
            word = _zeta_word(p)
            selectors = (range(1 << m) if p <= 2 else
                         [0, full, *[1 << u for u in range(m)]])
            for x, selectors_in in product(range(m), selectors):
                def moved_zeta(selectors):
                    return _translate(_binary_image(word, _translate(selectors, x, p)), x, p)
                first = moved_zeta(selectors_in)
                self.assertEqual(moved_zeta(first), selectors_in)
                self.assertEqual(moved_zeta(first ^ full), selectors_in ^ (1 << x))
            # These zero-plus-basis-column checks determine the complete
            # affine binary map at larger sizes; no phase or native depth
            # claim follows merely from its logical identity.

    def test_native_transport_and_phase_conjugation_all_small_columns(self):
        for p, bank_count in ((1, 1), (1, 2), (2, 1)):
            m, width = 1 << p, p + bank_count * (1 << p)
            schedule = _translation_schedule(p, bank_count)
            word = _word(schedule)
            reverse = _word(_positive_s_schedule(_inverse(schedule)))
            expected = zeros(1 << width)
            phase_expected = zeros(1 << width)
            for column in range(1 << width):
                x, image = column & ((1 << p)-1), column & ((1 << p)-1)
                for bank in range(bank_count):
                    bits = (column >> (p + bank*m)) & ((1 << m)-1)
                    image |= _translate(bits, x, p) << (p + bank*m)
                expected[image][column] = CO
                phase_expected[column][column] = (-CO if (column >> (p+x)) & 1 else CO)
            actual = native_columns(width, word)
            self.assertEqual(actual, expected)
            self.assertEqual(native_columns(width, reverse), adjoint(actual))
            self.assertEqual(native_columns(width, word + [('S', p), ('S', p)] + reverse),
                             phase_expected)

    def test_transport_schedules_pay_every_coordinate_and_share_no_t_target(self):
        for p, bank_count in product(range(1, 9), (1, 2)):
            m = 1 << p
            width, pairs = p + bank_count*m, bank_count*m//2
            schedule = _translation_schedule(p, bank_count)
            for orientation in (schedule, _positive_s_schedule(_inverse(schedule))):
                self.assertEqual(sum(kind == 'T' for kind, _ in orientation), 4*p)
                for kind, gates in orientation:
                    self.assertTrue(all(0 <= wire < width for gate in gates for wire in gate[1:]))
                    if kind == 'T':
                        targets = [gate[1] for gate in gates]
                        self.assertEqual(len(targets), len(set(targets)))
                        self.assertTrue(all(gate[0] in ('T', 'TDG') for gate in gates))
                        if pairs % 2 == 0:
                            self.assertTrue(all(wire >= p for wire in targets))
                    else:
                        self.assertTrue(all(gate[0] in ('H', 'S', 'CX') for gate in gates))
                t_count, clifford_count = _counts(_word(orientation))
                self.assertEqual(t_count, p * (6*pairs + pairs % 2))
                self.assertLessEqual(clifford_count, p * (14*pairs + 3))


if __name__ == '__main__':
    unittest.main()
