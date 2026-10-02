"""Bounded exact audits of a two-T-layer dirty-helper Toffoli and indicator.

The phase certificate follows the emitted CNOT/T words on all basis inputs.
Native matrices then check Hadamard-conjugated Toffoli and indicator words,
including their actual inverses and every arbitrary dirty-helper input.
No general shallow indicator construction is inferred from this base case.
"""
from __future__ import annotations

import unittest

import numpy as np

try:
    from .test_operator_source_compiler import _basis_action, _word_matrix
    from .test_t_depth import _inverse, _word
except ImportError:
    from test_operator_source_compiler import _basis_action, _word_matrix
    from test_t_depth import _inverse, _word


ATOL = 5e-11
PARITY_GROUPS = ((0, 1, 2, 4), (3, 5, 6, 7))


def _parity_basis(wires, masks):
    """Synthesize the specified invertible binary parity map using CNOTs."""
    assert len(wires) == len(masks) and len(set(wires)) == len(wires)
    assert all(0 <= mask < 1 << len(wires) for mask in masks)
    rows, elimination = list(masks), []
    for column in range(len(wires)):
        pivot = next((row for row in range(column, len(wires))
                      if rows[row] >> column & 1), None)
        assert pivot is not None
        if pivot != column:
            rows[pivot], rows[column] = rows[column], rows[pivot]
            elimination += [('CX', wires[pivot], wires[column]),
                            ('CX', wires[column], wires[pivot]),
                            ('CX', wires[pivot], wires[column])]
        for row in range(len(wires)):
            if row != column and rows[row] >> column & 1:
                rows[row] ^= rows[column]
                elimination.append(('CX', wires[column], wires[row]))
    assert rows == [1 << bit for bit in range(len(wires))]
    # Elimination maps the parity matrix to identity; reverse its literal
    # row operations to synthesize the requested parity matrix from identity.
    return list(reversed(elimination))


def _ccz_two_layers(a, b, c, helper):
    wires = [a, b, c, helper]
    assert len(set(wires)) == 4
    schedule = []
    for group in PARITY_GROUPS:
        basis = [('C', _parity_basis(wires, [8 | subset for subset in group]))]
        phases = [('TDG' if subset.bit_count() % 2 else 'T', wire)
                  for wire, subset in zip(wires, group)]
        schedule += basis + [('T', phases)] + _inverse(basis)
    return schedule


def _toffoli_two_layers(a, b, target, helper):
    hadamard = [('C', [('H', target)])]
    return hadamard + _ccz_two_layers(a, b, target, helper) + hadamard


def _two_bit_indicator(address, output):
    assert len(address) == 2 and len(output) == 4
    assert len(set(address + output)) == 6
    # Address lists remain little-endian: labels (a,b) select index 2*a+b.
    b, a = address
    affine = [('X', output[0]), ('CX', a, output[0]), ('CX', b, output[0]),
              ('CX', b, output[1]), ('CX', a, output[2])]
    fanout = [('C', [('CX', output[3], wire) for wire in output[:3]])]
    # output[0] is an arbitrary helper after the first fanout. Its exact
    # identity action inside Toffoli permits borrowing it without a seventh wire.
    return ([('C', affine)] + fanout
            + _toffoli_two_layers(a, b, output[3], output[0]) + _inverse(fanout))


def _exact_phase_action(basis, word):
    """Exact phase exponent modulo eight for an emitted CNOT/T-only word."""
    exponent = 0
    for gate in word:
        if gate[0] == 'CX':
            if basis >> gate[1] & 1:
                basis ^= 1 << gate[2]
        else:
            assert gate[0] in ('T', 'TDG')
            exponent += (1 if gate[0] == 'T' else -1) * (basis >> gate[1] & 1)
    return basis, exponent % 8


class DirtyIndicatorDepthTests(unittest.TestCase):
    def test_emitted_parity_bases_exact_ccz_phase_and_two_layer_ledger(self):
        for wires in ([0, 1, 2, 3], [2, 0, 3, 1]):
            for group in PARITY_GROUPS:
                masks = [8 | subset for subset in group]
                basis_word = _parity_basis(wires, masks)
                for basis in range(16):
                    packed = sum((basis >> wire & 1) << bit for bit, wire in enumerate(wires))
                    expected = sum(((packed & mask).bit_count() % 2) << wire
                                   for wire, mask in zip(wires, masks))
                    self.assertEqual(_basis_action(basis, basis_word), expected)
                    self.assertEqual(_basis_action(expected, list(reversed(basis_word))), basis)
            schedule = _ccz_two_layers(*wires)
            self.assertEqual(sum(kind == 'T' for kind, _ in schedule), 2)
            self.assertEqual(sum(len(gates) for kind, gates in schedule if kind == 'T'), 8)
            for kind, gates in schedule:
                self.assertTrue(all(wire in wires for gate in gates for wire in gate[1:]))
                if kind == 'T':
                    targets = [gate[1] for gate in gates]
                    self.assertEqual(len(targets), len(set(targets)))
            for basis in range(16):
                phase = 4 * all(basis >> wire & 1 for wire in wires[:3])
                self.assertEqual(_exact_phase_action(basis, _word(schedule)), (basis, phase))
                self.assertEqual(_exact_phase_action(basis, _word(_inverse(schedule))), (basis, phase))

    def test_native_toffoli_returns_arbitrary_helper_with_literal_phase(self):
        for a, b, target, helper in ((0, 1, 2, 3), (3, 0, 1, 2)):
            schedule = _toffoli_two_layers(a, b, target, helper)
            images = [basis ^ (1 << target) if (basis >> a & 1) and (basis >> b & 1)
                      else basis for basis in range(16)]
            expected = np.eye(16)[:, np.argsort(images)]
            np.testing.assert_allclose(_word_matrix(4, _word(schedule)), expected,
                                       atol=ATOL, rtol=0)
            np.testing.assert_allclose(_word_matrix(4, _word(_inverse(schedule))),
                                       expected.conj().T, atol=ATOL, rtol=0)
            # Full operator equality includes both helper values and their
            # off-diagonal coherences, not merely an initialized helper column.

    def test_native_indicator_uses_only_its_four_dirty_outputs(self):
        for address, output in (([0, 1], [2, 3, 4, 5]), ([4, 1], [0, 5, 2, 3])):
            schedule = _two_bit_indicator(address, output)
            self.assertEqual(sum(kind == 'T' for kind, _ in schedule), 2)
            self.assertEqual(sum(len(gates) for kind, gates in schedule if kind == 'T'), 8)
            self.assertEqual({wire for gate in _word(schedule) for wire in gate[1:]}, set(range(6)))
            images = []
            for basis in range(64):
                selected = sum((basis >> wire & 1) << bit for bit, wire in enumerate(address))
                images.append(basis ^ (1 << output[selected]))
            expected = np.eye(64)[:, np.argsort(images)]
            np.testing.assert_allclose(_word_matrix(6, _word(schedule)), expected,
                                       atol=ATOL, rtol=0)
            np.testing.assert_allclose(_word_matrix(6, _word(_inverse(schedule))),
                                       expected.conj().T, atol=ATOL, rtol=0)

    def test_missing_helper_phase_and_wrong_basis_inverse_are_detected(self):
        schedule = _ccz_two_layers(0, 1, 2, 3)
        # The first parity is d, although its encoded coordinate is wire zero.
        missing = [(kind, gates[1:] if index == 1 else gates)
                   for index, (kind, gates) in enumerate(schedule)]
        self.assertEqual(_exact_phase_action(8, _word(schedule)), (8, 0))
        self.assertEqual(_exact_phase_action(8, _word(missing)), (8, 7))
        # A synthesized parity basis need not be an involution. Reusing the
        # forward basis in place of its actual inverse changes basis states.
        wrong = [schedule[0], schedule[1], schedule[0]] + schedule[3:]
        self.assertTrue(any(_exact_phase_action(basis, _word(wrong))[0] != basis
                            for basis in range(16)))


if __name__ == '__main__':
    unittest.main()
