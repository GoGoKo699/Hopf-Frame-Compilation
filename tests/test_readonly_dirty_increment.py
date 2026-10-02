"""Bounded emitted checks of the involution-based read-only increment.

The public literal appears only in Clifford CNOTs; reflection controls and
borrowed increment work are private and arbitrary. The slow serial MCX
fixture checks exact action, phases, and return, not the imported fast
incrementer's asymptotic depth. Native matrices have order at most 64.
"""
from __future__ import annotations

import unittest

import numpy as np

try:
    from .test_counter_dirty_indicator import _parallel
    from .test_operator_source_compiler import _adjoint, _basis_action, _word_matrix
    from .test_t_depth import _inverse, _native_schedule, _word
except ImportError:
    from test_counter_dirty_indicator import _parallel
    from test_operator_source_compiler import _adjoint, _basis_action, _word_matrix
    from test_t_depth import _inverse, _native_schedule, _word


ATOL = 5e-11


def _borrowed_mcx(controls, target, helper):
    """A small serial exact word, using one private arbitrary dirty helper."""
    assert controls and len(set(controls + [target, helper])) == len(controls) + 2
    if len(controls) == 1:
        return [('CX', controls[0], target)]
    if len(controls) == 2:
        return [('CCX', *controls, target)]
    first = [('CCX', controls[0], controls[1], helper)]
    # The second word returns its borrowed first-control bit before the next
    # first word. All these wires are private in the public-literal wrapper.
    second = _borrowed_mcx([helper] + controls[2:], target, controls[0])
    return first + second + _adjoint(first) + _adjoint(second)


def _private_increment(control, counter, helper):
    return [gate for bit in reversed(range(len(counter)))
            for gate in _borrowed_mcx([control] + counter[:bit], counter[bit], helper)]


def _parts(control, counter, dirty, helper, *, negative=False):
    assert len(set([control, dirty, helper] + counter)) == len(counter) + 3
    public = [('CX', control, wire) for wire in counter]
    toggle = [('CX', control, dirty)]
    if negative:
        public += [('X', wire) for wire in counter]
        toggle += [('X', dirty)]
    reflection = ([('CX', dirty, wire) for wire in counter]
                  + _private_increment(dirty, counter, helper))
    return public, reflection, toggle


def _readonly_increment(control, counter, dirty, helper, *, negative=False):
    public, reflection, toggle = _parts(control, counter, dirty, helper, negative=negative)
    return public + reflection + toggle + reflection + toggle


def _read(basis, wires):
    return sum((basis >> wire & 1) << bit for bit, wire in enumerate(wires))


def _replace(basis, wires, value):
    for bit, wire in enumerate(wires):
        basis = (basis & ~(1 << wire)) | ((value >> bit & 1) << wire)
    return basis


class ReadonlyDirtyIncrementTests(unittest.TestCase):
    def test_all_basis_inputs_both_polarities_and_actual_inverse(self):
        for bits in range(1, 6):
            counter, dirty, helper = list(range(1, bits + 1)), bits + 1, bits + 2
            for negative in (False, True):
                word = _readonly_increment(0, counter, dirty, helper, negative=negative)
                for gate in word:
                    if 0 in gate[1:]:
                        self.assertEqual(gate[0:2], ('CX', 0))
                        self.assertNotEqual(gate[-1], 0)
                for basis in range(1 << (bits + 3)):
                    literal = (basis & 1) ^ negative
                    expected = _replace(basis, counter,
                                        (_read(basis, counter) + literal) % (1 << bits))
                    self.assertEqual(_basis_action(basis, word), expected)
                    self.assertEqual(_basis_action(expected, _adjoint(word)), basis)

    def test_native_phases_inverse_and_returned_private_helpers(self):
        # m=3 additionally exercises the borrowed helper inside a 3-control
        # gate, while still using only six wires and a 64-by-64 matrix.
        for bits in (1, 2, 3):
            width = bits + 3
            counter, dirty, helper = list(range(1, bits + 1)), bits + 1, bits + 2
            for negative in (False, True):
                word = _readonly_increment(0, counter, dirty, helper, negative=negative)
                schedule = _native_schedule(word)
                native = _word(schedule)
                for gate in native:
                    if 0 in gate[1:]:
                        self.assertEqual(gate[0:2], ('CX', 0))
                images = [_replace(basis, counter,
                                   (_read(basis, counter) + ((basis & 1) ^ negative)) % (1 << bits))
                          for basis in range(1 << width)]
                expected = np.eye(1 << width)[:, images]
                actual = _word_matrix(width, native)
                inverse_word = _word(_inverse(schedule))
                self.assertEqual(inverse_word, _adjoint(native))
                inverse = _word_matrix(width, inverse_word)
                np.testing.assert_allclose(actual, expected, atol=ATOL, rtol=0)
                np.testing.assert_allclose(inverse, expected.conj().T, atol=ATOL, rtol=0)
                np.testing.assert_allclose(inverse @ actual, np.eye(1 << width),
                                           atol=ATOL, rtol=0)

    def test_opposite_literals_share_only_clifford_controls_in_a_native_batch(self):
        first_counter, second_counter = [1, 2], [5, 6]
        first = _readonly_increment(0, first_counter, 3, 4)
        second = _readonly_increment(0, second_counter, 7, 8, negative=True)
        batch = _parallel([(first, _native_schedule(first)),
                           (second, _native_schedule(second))])
        for kind, gates in batch[1]:
            if kind == 'T':
                targets = [gate[1] for gate in gates]
                self.assertEqual(len(targets), len(set(targets)))
                self.assertNotIn(0, targets)
            for gate in gates:
                if 0 in gate[1:]:
                    self.assertEqual(gate[0:2], ('CX', 0))
        separate_depth = sum(kind == 'T' for kind, _ in _native_schedule(first))
        self.assertEqual(sum(kind == 'T' for kind, _ in batch[1]), separate_depth)
        for basis in range(512):
            literal = basis & 1
            expected = _replace(basis, first_counter, (_read(basis, first_counter) + literal) % 4)
            expected = _replace(expected, second_counter,
                                (_read(basis, second_counter) + 1 - literal) % 4)
            self.assertEqual(_basis_action(basis, batch[0]), expected)
            self.assertEqual(_basis_action(expected, _adjoint(batch[0])), basis)

    def test_missing_mask_cancellation_and_wrong_reflection_order_are_detected(self):
        counter, dirty, helper = [1, 2], 3, 4
        public, reflection, toggle = _parts(0, counter, dirty, helper)
        correct = public + reflection + toggle + reflection + toggle
        missing = public + toggle + reflection + toggle
        basis = (1 << dirty) | (1 << counter[0])  # ell=0, d=1, U=1
        self.assertEqual(_basis_action(basis, correct), basis)
        self.assertEqual(_read(_basis_action(basis, missing), counter), 3)
        self.assertEqual(_basis_action(basis, missing) >> dirty & 1, 1)
        wrong_order = reflection + toggle + reflection + toggle + public
        basis = 1  # ell=1, U=0
        self.assertEqual(_read(_basis_action(basis, correct), counter), 1)
        self.assertEqual(_read(_basis_action(basis, wrong_order), counter), 3)
        # Removing one native phase is also visible on the complete space.
        native = _word(_native_schedule(correct))
        missing_phase = list(native)
        missing_phase.pop(next(index for index, gate in enumerate(missing_phase)
                               if gate[0] in ('T', 'TDG')))
        self.assertGreater(np.linalg.norm(_word_matrix(5, native)
                                         - _word_matrix(5, missing_phase)), 0.5)


if __name__ == '__main__':
    unittest.main()
