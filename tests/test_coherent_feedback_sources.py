"""Exact checks for paid HRT key mixers and complete feedback interfaces.

The encoded-Hadamard circuits are native components, not an implementation
of the residual probe or an improved general complete-frame depth bound.
"""
from collections import Counter
from functools import lru_cache
import unittest

from tests.exact_matrices import CI, CO, CZ
from tests.test_common_probe_source import ROOT_HALF, _native
from tests.test_complete_observable_bases import (
    _hadamard_numerator, _hrt_unit, _permutation, _quadratic_inverse,
)
from tests.test_joint_lookup_pauli_source import _bands, _inverse


def _shared_ccz(control, pairs):
    """Literal four-band batch, with reduced positive-S root phase."""
    if not pairs:
        return []
    wires = [control] + [wire for pair in pairs for wire in pair]
    assert len(set(wires)) == len(wires)
    count = len(pairs)
    first = [('T', wire) for pair in pairs for wire in pair]
    if count & 1:
        first.append(('T', control))
    shared = [('CX', control, wire) for pair in pairs for wire in pair]
    paired = [('CX', left, right) for left, right in pairs]
    triple = [('CX', control, right) for _, right in pairs] + paired
    return ([('S', control)] * ((count // 2) % 4) + first
            + shared + [('TDG', wire) for pair in pairs for wire in pair]
            + list(reversed(shared)) + paired
            + [('TDG', right) for _, right in pairs]
            + list(reversed(paired)) + triple
            + [('T', right) for _, right in pairs]
            + list(reversed(triple)))


def _encoded_shears(r, subset):
    """Physical subset S occupies wire S-1; no constant or helper wire."""
    assert 0 < subset < (1 << r)
    root = subset - 1
    bipartitions = [(left - 1, (subset ^ left) - 1)
                   for left in range(1, subset)
                   if left & subset == left and left < (subset ^ left)]
    matching = [(left - 1, (subset | left) - 1)
                for left in range(1, 1 << r) if not left & subset]
    b = ([('H', root)] + _shared_ccz(root, bipartitions) + [('H', root)]
         if bipartitions else [])
    walls = [('H', right) for _, right in matching]
    c = walls + _shared_ccz(root, matching) + list(reversed(walls))
    return b, c


def _encoded_h(r, subset):
    b, c = _encoded_shears(r, subset)
    return c + b + [('H', subset - 1)] + b + c


@lru_cache(maxsize=None)
def _hrt_tables(r):
    width = (1 << r) - 1
    forward = tuple((_hrt_unit(key << 1, r) ^ 1) >> 1
                    for key in range(1 << width))
    inverse = [None] * len(forward)
    for key, output in enumerate(forward):
        assert inverse[output] is None
        inverse[output] = key
    return forward, tuple(inverse)


def _encoded_h_column(r, subset, basis):
    """Independent F H_S F-dagger column from the product definition of F."""
    width = (1 << r) - 1
    forward, inverse = _hrt_tables(r)
    carrier, rest = basis & ((1 << width) - 1), basis >> width
    key = inverse[carrier]
    mask = 1 << (subset - 1)
    zero, one = key & ~mask, key | mask
    return {forward[zero] | (rest << width): ROOT_HALF,
            forward[one] | (rest << width): -ROOT_HALF if key & mask else ROOT_HALF}


def _encoded_h_ideal(r, subset, state):
    answer = {}
    for basis, amplitude in state.items():
        for row, value in _encoded_h_column(r, subset, basis).items():
            answer[row] = answer.get(row, CZ) + amplitude * value
    return {row: value for row, value in answer.items() if value != CZ}


def _feedback_tokens(r, calibrated, both_labels=False):
    width = (1 << r) - 1
    c, a, b = tuple(range(width)), tuple(range(width, 2 * width)), tuple(range(2 * width, 3 * width))
    word = ([('B', 1)] + [('H', wire) for wire in a]) if calibrated else []
    word += [('H', wire) for wire in a]
    word += [('CX', source, target) for source, target in zip(b, c)]
    word += [('Q', 1, 0)] + [('H', wire) for wire in b]
    # (u,y,t) -> (t,u,y), retaining every old dirty register.
    word += [('SWAP', left, right) for left, right in zip(c, b)]
    word += [('SWAP', left, right) for left, right in zip(a, b)]
    word += [('B', 0), ('H', c[-1]), ('B', 0)]
    if both_labels:
        word += [('CX', source, target) for source, target in zip(a, c)]
    word += [('Q', 0, 1)]
    return word


def _feedback_integer(state, word, r):
    """Exact sparse numerator; each literal H contributes sqrt(2)."""
    width = (1 << r) - 1
    mask = (1 << width) - 1
    top = width

    def bank(basis, index):
        return (basis >> (index * width)) & mask

    for name, *args in word:
        if name == 'H':
            state = _hadamard_numerator(state, args[0])
        elif name == 'CX':
            source, target = args
            state = _permutation(state, lambda value: value ^ (((value >> source) & 1) << target))
        elif name == 'SWAP':
            left, right = args
            state = _permutation(state, lambda value: value ^ (
                (((value >> left) ^ (value >> right)) & 1) * ((1 << left) | (1 << right))))
        elif name == 'Q':
            source, target = args
            state = _permutation(state, lambda value: value ^ (
                (_quadratic_inverse(bank(value, source) << 1, r) >> 1) << (target * width)))
        elif name == 'B':
            index = args[0]
            state = _permutation(state, lambda value: value ^ (
                (((_quadratic_inverse(bank(value, index) << 1, r) >> 1)
                  ^ bank(value, index)) & (1 << (top - 1))) << (index * width)))
        elif name == 'Z':
            state = {value: -amplitude if (value >> args[0]) & 1 else amplitude
                     for value, amplitude in state.items()}
        elif name == 'DK':
            index = args[0]
            state = {value: -amplitude if (_quadratic_inverse(bank(value, index) << 1, r) >> top) & 1
                     else amplitude for value, amplitude in state.items()}
        else:
            raise AssertionError(name)
    return state


def _feedback_native(r, calibrated, both_labels=False):
    """Small literal circuit for the same complete feedback word."""
    width = (1 << r) - 1
    native = []
    for name, *args in _feedback_tokens(r, calibrated, both_labels):
        if name in ('H', 'CX'):
            native.append((name, *args))
        elif name == 'SWAP':
            left, right = args
            native += [('CX', left, right), ('CX', right, left), ('CX', left, right)]
        elif name == 'B':
            b, _ = _encoded_shears(r, width)
            native += [(gate[0], *(wire + args[0] * width for wire in gate[1:])) for gate in b]
        elif name == 'Q':
            source, target = args
            native += [('CX', source * width + subset - 1, target * width + subset - 1)
                       for subset in range(1, 1 << r)]
            for subset in range(1, 1 << r):
                for left in range(1, subset):
                    right = subset ^ left
                    if left & subset != left or left >= right:
                        continue
                    a, b = source * width + left - 1, source * width + right - 1
                    c = target * width + subset - 1
                    native += [('H', c)] + _shared_ccz(a, [(b, c)]) + [('H', c)]
        else:
            raise AssertionError(name)
    return native


def _both_label_feedback_target(r, basis):
    width = (1 << r) - 1
    mask = (1 << width) - 1
    c, a, b = basis & mask, (basis >> width) & mask, (basis >> (2 * width)) & mask
    v = (c ^ b ^ (_quadratic_inverse(a << 1, r) >> 1)) << 1
    reflected = sum(((v >> (width ^ subset)) & 1) << (subset - 1)
                    for subset in range(1, width))
    phase = ((b >> (width - 1)) ^ ((_quadratic_inverse(v, r) ^ v) >> width)) & 1
    moved = basis ^ reflected ^ (reflected << (2 * width))
    return moved, -1 if phase else 1


class CoherentFeedbackSourcesTests(unittest.TestCase):
    def test_shared_root_batch_literal_phase_and_actual_reverse(self):
        for count in (1, 2, 3, 8, 9):
            pairs = [(2 * index + 1, 2 * index + 2) for index in range(count)]
            word = _shared_ccz(0, pairs)
            full = (1 << (2 * count + 1)) - 1
            columns = (range(1 << (2 * count + 1)) if count <= 3
                       else (0, 1, 7, full, full ^ 1, full ^ 2))
            for basis in columns:
                parity = (basis & 1) * sum(((basis >> left) & 1)
                                           * ((basis >> right) & 1)
                                           for left, right in pairs)
                expected = {basis: -CO if parity & 1 else CO}
                self.assertEqual(_native(word, {basis: CO}), expected)
                self.assertEqual(_native(_inverse(word), {basis: CO}), expected)
            for oriented in (word, _inverse(word)):
                self.assertLessEqual(len(_bands(oriented)), 4)
                self.assertTrue(all(len({gate[1] for gate in band}) == len(band)
                                    for band in _bands(oriented)))

    def test_encoded_h_complete_native_columns_against_independent_hrt(self):
        for r in (1, 2, 3):
            for subset in range(1, 1 << r):
                word = _encoded_h(r, subset)
                for basis in range(1 << ((1 << r) - 1)):
                    expected = _encoded_h_column(r, subset, basis)
                    self.assertEqual(_native(word, {basis: CO}), expected)

    def test_encoded_h_actual_inverse_and_coherent_reference(self):
        for r in (1, 2, 3):
            width = (1 << r) - 1
            occupied = (1 << width) - 1
            for subset in range(1, 1 << r):
                word = _encoded_h(r, subset)
                # Two untouched physical work wires and a reference retain
                # their original, correlated values and literal relative i.
                state = {occupied: ROOT_HALF,
                         (occupied // 3) | (7 << width): CI * ROOT_HALF}
                actual = _native(word, state)
                self.assertEqual(actual, _encoded_h_ideal(r, subset, state))
                self.assertEqual(_native(_inverse(word), actual), state)
                self.assertEqual(_native(_inverse(word), state), actual)

    def test_encoded_h_native_depth_count_and_no_helper_allocation(self):
        for r in (1, 2, 3, 4, 7):
            for size in range(1, r + 1):
                subset = (1 << size) - 1
                a, b = (1 << (size - 1)) - 1, (1 << (r - size)) - 1
                word = _encoded_h(r, subset)
                for oriented in (word, _inverse(word)):
                    count = Counter(gate[0] for gate in oriented)
                    self.assertEqual(count['T'] + count['TDG'],
                                     12 * (a + b) + 2 * ((a & 1) + (b & 1)))
                    self.assertLessEqual(count['H'] + count['CX'] + count['S'],
                                         20 * a + 24 * b + 17)
                    self.assertLessEqual(len(_bands(oriented)), 16)
                    for band in _bands(oriented):
                        self.assertEqual(len({gate[1] for gate in band}), len(band))
                    self.assertTrue(all(0 <= wire < (1 << r) - 1
                                        for gate in oriented for wire in gate[1:]))
                if size in (1, r):
                    self.assertLessEqual(len(_bands(word)), 8)

    def test_encoded_h_changes_key_algebra_with_the_literal_pauli_image(self):
        for r in (1, 2, 3):
            width = (1 << r) - 1
            for subset in range(1, 1 << r):
                for basis in range(1 << width):
                    first = _encoded_h_column(r, subset, basis)
                    phased = {row: -value if (row >> (subset - 1)) & 1 else value
                              for row, value in first.items()}
                    actual = _encoded_h_ideal(r, subset, phased)
                    coefficients = basis << 1
                    inverse = _quadratic_inverse(coefficients, r)
                    q = ((inverse ^ coefficients) >> subset) & 1
                    moved = coefficients ^ (1 << subset)
                    for left in range(1, 1 << r):
                        if not left & subset and (coefficients >> left) & 1:
                            moved ^= 1 << (left | subset)
                    self.assertEqual(actual, {moved >> 1: -CO if q else CO})
                    self.assertNotEqual(moved >> 1, basis)

    def test_top_only_collision_distinguishes_inverse_from_forward_phase(self):
        for r in range(3, 7):
            top = (1 << r) - 1
            singletons = sum(1 << (1 << bit) for bit in range(r))
            for coefficients in (0, singletons):
                self.assertEqual((_quadratic_inverse(coefficients, r) >> top) & 1, 0)
            self.assertEqual((_hrt_unit(0, r) >> top) & 1, 0)
            self.assertEqual((_hrt_unit(singletons, r) >> top) & 1, 1)
        # At r=3 the literal top-only mixer acts as the same H on both
        # backgrounds, whereas the required forward top phases are opposite.
        r, top = 3, 7
        word = _encoded_h(r, top)
        for coefficients in (0, sum(1 << (1 << bit) for bit in range(r))):
            basis = coefficients >> 1
            self.assertEqual(_native(word, {basis: CO}),
                             {basis: ROOT_HALF, basis ^ (1 << (top - 1)): ROOT_HALF})

    def test_retained_outcome_feedback_full_uncalibrated_probe(self):
        for r in (1, 2):
            width = (1 << r) - 1
            source = _feedback_tokens(r, False)
            echo = source + [('Z', 2 * width - 1)] + list(reversed(source))
            target = ([('H', wire) for wire in range(width, 2 * width)]
                      + [('DK', 1)]
                      + [('H', wire) for wire in range(width, 2 * width)]
                      + [('Z', width - 1)])
            scale = 1 << (width + 1)
            for basis in range(1 << (3 * width)):
                actual = _feedback_integer({basis: 1}, echo, r)
                expected = _feedback_integer({basis: 1}, target, r)
                self.assertEqual(actual, {row: scale * value for row, value in expected.items()})

    def test_calibrated_feedback_complete_return_and_maximal_error_witness(self):
        for r in (1, 2):
            width = (1 << r) - 1
            source = _feedback_tokens(r, True)
            echo = source + [('Z', 2 * width - 1)] + list(reversed(source))
            scale = 1 << sum(gate[0] == 'H' for gate in source)
            for basis in range(1 << (3 * width)):
                parity = ((basis >> (width - 1)) ^ (basis >> (2 * width - 1))) & 1
                self.assertEqual(_feedback_integer({basis: 1}, echo, r),
                                 {basis: -scale if parity else scale})
            # Include arbitrary work and an external reference coherently;
            # the completed diagonal action preserves all their labels.
            mask = (1 << (3 * width)) - 1
            state = {mask: 1, (mask // 3) | (1 << (3 * width)): 2}
            expected = {basis: amplitude * scale * (-1 if (
                ((basis >> (width - 1)) ^ (basis >> (2 * width - 1))) & 1) else 1)
                        for basis, amplitude in state.items()}
            self.assertEqual(_feedback_integer(state, echo, r), expected)
        # The analytic witness lifts to every r>=2. It retains an arbitrary
        # Bell bank b; the computed probe is +1 while Xi(c,a) is -1.
        for r in range(2, 7):
            top = (1 << r) - 1
            c, a = 1 << 1, 1 << (top ^ 1)
            corrected = ((c ^ a) >> top) & 1
            xi = (_hrt_unit(c ^ _quadratic_inverse(a, r), r) >> top) & 1
            self.assertEqual(corrected, 0)
            self.assertEqual(xi, 1)

    def test_calibrated_feedback_literal_native_phases_and_actual_adjoint(self):
        for r in (1, 2):
            width = (1 << r) - 1
            source = _feedback_native(r, True)
            echo = source + [('S', 2 * width - 1)] * 2 + _inverse(source)
            if r == 1:
                states = [{basis: CO} for basis in range(8)]
            else:
                # A witness column plus a coherent dirty/reference input,
                # evaluated through actual native words, not oracle tokens.
                c, a, b = 1, 2, 5
                witness = c | (a << width) | (b << (2 * width))
                states = [{witness: CO},
                          {witness: ROOT_HALF, 511 | (1 << (3 * width)): CI * ROOT_HALF}]
            for state in states:
                expected = {basis: -amplitude if (
                    ((basis >> (width - 1)) ^ (basis >> (2 * width - 1))) & 1) else amplitude
                            for basis, amplitude in state.items()}
                self.assertEqual(_native(echo, state), expected)
                self.assertEqual(_native(_inverse(source), _native(source, state)), state)

    def test_both_label_feedback_complete_signed_permutation_and_witness(self):
        for r in (1, 2):
            width = (1 << r) - 1
            source = _feedback_tokens(r, True, both_labels=True)
            echo = source + [('Z', 2 * width - 1)] + list(reversed(source))
            scale = 1 << sum(gate[0] == 'H' for gate in source)
            for basis in range(1 << (3 * width)):
                moved, sign = _both_label_feedback_target(r, basis)
                self.assertEqual(_feedback_integer({basis: 1}, echo, r), {moved: sign * scale})
            native_source = _feedback_native(r, True, both_labels=True)
            native_echo = (native_source + [('S', 2 * width - 1)] * 2
                           + _inverse(native_source))
            state = {1 << (2 * width - 1): ROOT_HALF,
                     (1 << (3 * width)) | 1: CI * ROOT_HALF}
            expected = {}
            for basis, amplitude in state.items():
                moved, sign = _both_label_feedback_target(r, basis)
                expected[moved] = sign * amplitude
            self.assertEqual(_native(native_echo, state), expected)
            self.assertEqual(_native(_inverse(native_source), _native(native_source, state)), state)
        for r in range(1, 7):
            width = (1 << r) - 1
            # c=b=0, a=e_top: both-label feedback is +1, Xi is -1.
            basis = 1 << (2 * width - 1)
            self.assertEqual(_both_label_feedback_target(r, basis), (basis, 1))
            top_coefficient = 1 << width
            self.assertEqual((_hrt_unit(_quadratic_inverse(top_coefficient, r), r) >> width) & 1, 1)


if __name__ == '__main__':
    unittest.main()
