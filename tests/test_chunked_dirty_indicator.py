"""Bounded all-input audits of the chunked dirty-indicator tree.

Tree fixtures emit a slow dirty-ladder conjunction for each edge, not the
asymptotic counter schedule. Separate counter-row basis checks audit its
substitutable exact interface. Six-wire native matrices check literal
phases and shared-control scheduling; integer ledgers check the proposed
chunk interpolation. General counter depth remains an analytic import.
"""
from __future__ import annotations

from fractions import Fraction
import random
import unittest

import numpy as np

try:
    from .test_counter_dirty_indicator import _counter_fixture, _parallel
    from .test_operator_source_compiler import _adjoint, _basis_action, _word_matrix
    from .test_pipelined_dirty_sum import (
        _and_oracle, _expand, _form, _not, _private_native,
    )
    from .test_t_depth import _word
except ImportError:
    from test_counter_dirty_indicator import _counter_fixture, _parallel
    from test_operator_source_compiler import _adjoint, _basis_action, _word_matrix
    from test_pipelined_dirty_sum import (
        _and_oracle, _expand, _form, _not, _private_native,
    )
    from test_t_depth import _word


ATOL = 5e-11


def _schedule(native):
    stages, cliffords = [], []
    for gate in native:
        if gate[0] in ('T', 'TDG'):
            if cliffords:
                stages.append(('C', cliffords))
                cliffords = []
            stages.append(('T', [gate]))
        else:
            cliffords.append(gate)
    if cliffords:
        stages.append(('C', cliffords))
    return stages


def _tree_fixture(bits, chunk):
    """Little-endian chunk labels; every node and ladder helper is dirty."""
    assert 1 <= chunk <= bits
    address = list(range(bits))
    output = list(range(bits, bits + (1 << bits)))
    next_wire = output[-1] + 1
    boundaries = list(range(chunk, bits, chunk)) + [bits]
    nodes = [[next_wire]]
    next_wire += 1
    for boundary in boundaries:
        nodes.append(list(range(next_wire, next_wire + (1 << boundary))))
        next_wire += 1 << boundary
    helper_start = next_wire
    # Private per edge within a level; reuse only after the completed level.
    next_wire += (1 << bits) * (chunk + 1)
    levels, offset = [], 0
    for depth, boundary in enumerate(boundaries):
        length = boundary - offset
        rows = []
        for child_label, target in enumerate(nodes[depth + 1]):
            parent = nodes[depth][child_label % (1 << offset)]
            local_label = child_label >> offset
            literals = [_form(parent)] + [
                _form(address[offset + bit], constant=1 - (local_label >> bit & 1))
                for bit in range(length)]
            start = helper_start + child_label * (chunk + 1)
            helpers = list(range(start, start + length + 1))
            affine_word = _and_oracle(literals, target, helpers)
            rows.append(dict(word=_expand(affine_word), parent=parent,
                             target=target, helpers=helpers))
        levels.append(rows)
        offset = boundary
    forward = [gate for level in levels for row in level for gate in row['word']]
    # Algebraic F X_root F^dagger is chronological F^dagger, X_root, F.
    path_flip = _adjoint(forward) + [('X', nodes[0][0])] + forward
    copy = [('CX', leaf, target) for leaf, target in zip(nodes[-1], output)]
    echo = copy + path_flip + copy + _adjoint(path_flip)
    return dict(width=next_wire, address=address, output=output, nodes=nodes,
                boundaries=boundaries, levels=levels, forward=forward,
                path_flip=path_flip, copy=copy, echo=echo)


def _ceil_log2(value):
    assert value >= 1
    return (value - 1).bit_length()


class ChunkedDirtyIndicatorTests(unittest.TestCase):
    def test_native_rows_shared_controls_and_counter_interface(self):
        # Two edges share parent/address but have distinct targets and phase
        # buses. Every shared input is only a Clifford CNOT control.
        rows = []
        for value, target, phase in ((0, 2, 4), (1, 3, 5)):
            literal = _form(1) if value else _not(_form(1))
            affine = [('AA', _form(0), literal, target)]
            native = _private_native(affine, phase)
            self.assertTrue(all(gate[-1] not in (0, 1) for gate in native))
            rows.append((_expand(affine), _schedule(native)))
        classical, schedule = _parallel(rows)
        self.assertEqual(sum(kind == 'T' for kind, _ in schedule), 8)
        for kind, gates in schedule:
            if kind == 'T':
                targets = [gate[1] for gate in gates]
                self.assertEqual(len(targets), len(set(targets)))
                self.assertTrue(set(targets) <= {4, 5})
        images = [basis ^ ((1 << (2 + (basis >> 1 & 1))) if basis & 1 else 0)
                  for basis in range(64)]
        self.assertEqual([_basis_action(basis, classical) for basis in range(64)], images)
        native = _word(schedule)
        expected = np.eye(64)[:, images]
        np.testing.assert_allclose(_word_matrix(6, native), expected, atol=ATOL, rtol=0)
        np.testing.assert_allclose(_word_matrix(6, _adjoint(native)), expected.conj().T,
                                   atol=ATOL, rtol=0)

        # A two-address-bit edge plus parent uses one borrowed AND mask and
        # one phase bus: still only six physical wires, all 64 inputs checked.
        for pattern in range(4):
            literals = [_form(wire, constant=1 - (pattern >> bit & 1))
                        for bit, wire in enumerate((1, 2))]
            first = ('AA', _form(0), literals[0], 4)
            second = ('AA', _form(4), literals[1], 3)
            affine = [first, second, first, second]
            native = _private_native(affine, 5)
            images = [basis ^ ((1 << 3) if basis & 1 and (basis >> 1 & 3) == pattern else 0)
                      for basis in range(64)]
            expected = np.eye(64)[:, images]
            np.testing.assert_allclose(_word_matrix(6, native), expected,
                                       atol=ATOL, rtol=0)
            np.testing.assert_allclose(_word_matrix(6, _adjoint(native)), expected.conj().T,
                                       atol=ATOL, rtol=0)

        # The existing serial counter fixture supplies exactly the same row
        # oracle. Its gate order is not used as evidence for shallow depth.
        generator = random.Random(0xC01A)
        for controls in (2, 3):
            for local_label in range(1 << (controls - 1)):
                pattern = 1 | (local_label << 1)  # first literal is dirty parent=1
                fixture = _counter_fixture(controls, pattern)
                word = fixture['conjunction'][0]
                backgrounds = [0, (1 << fixture['width']) - 1,
                               generator.getrandbits(fixture['width'])]
                for background in backgrounds:
                    for address in range(1 << controls):
                        initial = (background & ~((1 << controls) - 1)) | address
                        expected = initial ^ ((1 << fixture['target']) if address == pattern else 0)
                        self.assertEqual(_basis_action(initial, word), expected)
                        self.assertEqual(_basis_action(expected, _adjoint(word)), initial)

    def test_full_tree_echo_dirty_inputs_and_inverse_orientation(self):
        generator = random.Random(0xC4A9)
        for bits in range(1, 5):
            for chunk in sorted({1, min(2, bits), bits}):
                fixture = _tree_fixture(bits, chunk)
                self.assertLess(sum(map(len, fixture['nodes'])), 2 * (1 << bits))
                for level in fixture['levels']:
                    private = [set(row['helpers'] + [row['target']]) for row in level]
                    self.assertEqual(sum(map(len, private)), len(set().union(*private)))
                    shared = set(fixture['address']) | {row['parent'] for row in level}
                    self.assertTrue(all(gate[-1] not in shared
                                        for row in level for gate in row['word']))
                if bits == 1:
                    # Exhaust every background at the smallest tree size.
                    backgrounds = [value << bits
                                   for value in range(1 << (fixture['width'] - bits))]
                else:
                    backgrounds = [0, (1 << fixture['width']) - 1]
                    backgrounds += [generator.getrandbits(fixture['width']) for _ in range(5)]
                for background in backgrounds:
                    for address in range(1 << bits):
                        initial = (background & ~((1 << bits) - 1)) | address
                        path = 1 << fixture['nodes'][0][0]
                        for boundary, nodes in zip(fixture['boundaries'], fixture['nodes'][1:]):
                            path ^= 1 << nodes[address % (1 << boundary)]
                        self.assertEqual(_basis_action(initial, fixture['path_flip']), initial ^ path)
                        expected = initial ^ (1 << fixture['output'][address])
                        self.assertEqual(_basis_action(initial, fixture['echo']), expected)
                        self.assertEqual(_basis_action(expected, _adjoint(fixture['echo'])), initial)
                if len(fixture['boundaries']) > 1:
                    wrong = (fixture['forward'] + [('X', fixture['nodes'][0][0])]
                             + _adjoint(fixture['forward']))
                    wrong_echo = fixture['copy'] + wrong + fixture['copy'] + _adjoint(wrong)
                    # Bottom-up conjugation does not transport the root to leaves.
                    self.assertEqual(_basis_action(0, wrong_echo), 0)
                    self.assertNotEqual(_basis_action(0, fixture['echo']), 0)

    def test_integer_chunk_budget_and_late_query_sums(self):
        for bits in range(1, 65):
            for chunk in range(1, bits + 1):
                boundaries = list(range(chunk, bits, chunk)) + [bits]
                node_count = 1 + sum(1 << boundary for boundary in boundaries)
                self.assertLess(node_count, 2 * (1 << bits))
                self.assertEqual(len(boundaries), (bits + chunk - 1) // chunk)

        precision = 6
        q = Fraction(1, 64)
        # Group k into blocks of 12. This is an exact convergent-series
        # certificate, not a floating-point fit to sampled resource totals.
        count_bound = 12 * (Fraction(precision + 15) / (1 - q)
                            + 12 * q / (1 - q) ** 2
                            + Fraction(27) / (1 - Fraction(1, 8)))
        for n in (64, 256, 1024, 4096):
            tail = min(n, 32 * _ceil_log2(n + 2))
            cap = _ceil_log2(8 * n)
            count_sum, clifford_sum = Fraction(0), Fraction(0)
            depth_sum = 0
            for remaining in range(1, tail + 1):
                address_bits = (n - remaining + 3) // 2  # larger balanced half
                block = min(address_bits, 1 << (remaining // 12))
                word_bits = precision + 4 + min(remaining, cap)
                overhead = (block + 2) ** 3
                # 2^(-floor(k/2)) conservatively bounds 2^(-k/2).
                count_sum += Fraction(word_bits + overhead, 1 << (remaining // 2))
                clifford_sum += Fraction(word_bits, 1 << remaining)
                self.assertLessEqual(overhead, 27 * (1 << (remaining // 2)))
                levels = (address_bits + block - 1) // block
                depth_sum += word_bits + levels * _ceil_log2(block + 2)
            self.assertLessEqual(count_sum, count_bound)
            self.assertLessEqual(clifford_sum, precision + 6)
            # sum_{t>=0} 12(t+2)/2^t=72 controls the n-dependent depth.
            depth_bound = 72 * (n + 2) + tail * (precision + 4 + cap + tail // 12 + 2)
            self.assertLessEqual(depth_sum, depth_bound)


if __name__ == '__main__':
    unittest.main()
