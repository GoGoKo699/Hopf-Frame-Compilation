"""Exact bounded audits of unequal dirty-indicator chunk schedules.

Small complete trees use literal reversible X/CX/CCX words and exact
Boolean polynomials in every input wire. Separate small native conjunction
matrices cover the interfaces composed by the (1,2) and (2,1) trees.
Their slow toy ladders are not used to infer the imported shallow-counter
depth. Large-address checks allocate only a short list of chunk lengths
and normalized rational ledgers, never an exponentially large register.
"""
from __future__ import annotations

from fractions import Fraction
from math import isqrt
import unittest

import numpy as np

try:
    from .test_operator_source_compiler import _adjoint, _word_matrix
    from .test_parallel_dirty_lookup import _expected, _symbolic
    from .test_pipelined_dirty_sum import _and_oracle, _expand, _form, _private_native
except ImportError:
    from test_operator_source_compiler import _adjoint, _word_matrix
    from test_parallel_dirty_lookup import _expected, _symbolic
    from test_pipelined_dirty_sum import _and_oracle, _expand, _form, _private_native


ATOL = 5e-11


def _ceil_log2(value):
    return (value - 1).bit_length()


def _partition(bits):
    """Test-only transcription of the analytic remaining-length rule."""
    chunks, remaining = [], bits
    while remaining > 64:
        following = _ceil_log2((remaining + 2) ** 4)
        chunks.append(remaining - following)
        remaining = following
    if remaining:
        chunks.append(remaining)
    return chunks


def _lengths():
    cases = set(range(8193))
    cases.update(range(8193, 65_537, 257))
    cases.update((65_535, 65_536, 65_537, 999_999, 1_000_000))
    # Independent fourth-root boundaries for each attainable ceiling value.
    for exponent in range(24, 80):
        boundary = isqrt(isqrt(1 << exponent)) - 2
        cases.update(r for r in (boundary - 1, boundary, boundary + 1)
                     if 0 <= r <= 1_000_000)
    return sorted(cases)


def _unequal_tree(chunks):
    bits = sum(chunks)
    assert all(length > 0 for length in chunks)
    if not bits:
        return dict(width=1, address=[], output=[0], nodes=[], boundaries=[],
                    levels=[], forward=[], path_flip=[], copy=[], echo=[('X', 0)])
    address = list(range(bits))
    output = list(range(bits, bits + (1 << bits)))
    next_wire = output[-1] + 1
    nodes, boundaries, offset = [[next_wire]], [], 0
    next_wire += 1
    for length in chunks:
        offset += length
        boundaries.append(offset)
        nodes.append(list(range(next_wire, next_wire + (1 << offset))))
        next_wire += 1 << offset
    helper_start = next_wire
    next_wire += max((1 << boundary) * (length + 1)
                     for boundary, length in zip(boundaries, chunks))
    levels, offset = [], 0
    for depth, (boundary, length) in enumerate(zip(boundaries, chunks)):
        rows = []
        for label, target in enumerate(nodes[depth + 1]):
            parent = nodes[depth][label % (1 << offset)]
            local = label >> offset
            literals = [_form(parent)] + [
                _form(address[offset + bit], constant=1 - ((local >> bit) & 1))
                for bit in range(length)]
            start = helper_start + label * (length + 1)
            helpers = list(range(start, start + length + 1))
            rows.append(dict(parent=parent, target=target, helpers=helpers,
                             word=_expand(_and_oracle(literals, target, helpers))))
        levels.append(rows)
        offset = boundary
    forward = [gate for level in levels for row in level for gate in row['word']]
    path = _adjoint(forward) + [('X', nodes[0][0])] + forward
    copy = [('CX', leaf, target) for leaf, target in zip(nodes[-1], output)]
    echo = copy + path + copy + _adjoint(path)
    return dict(width=next_wire, address=address, output=output, nodes=nodes,
                boundaries=boundaries, levels=levels, forward=forward,
                path_flip=path, copy=copy, echo=echo)


class NonuniformDirtyIndicatorTests(unittest.TestCase):
    def test_unequal_complete_trees_on_every_symbolic_dirty_input(self):
        for chunks in ((), (1,), (1, 2), (2, 1), (1, 3), (3, 1),
                       (1, 2, 1), (2, 1, 1)):
            f = _unequal_tree(chunks)
            initial = [{1 << wire} for wire in range(f['width'])]
            table = [1 << x for x in range(1 << sum(chunks))]
            expected = _expected(initial, f['address'], f['output'], table)
            actual = _symbolic(f['echo'], initial)
            self.assertEqual(actual, expected)
            self.assertEqual(_symbolic(_adjoint(f['echo']), actual), initial)
            if chunks:
                self.assertLess(sum(map(len, f['nodes'])), 2 * (1 << sum(chunks)))
            for level in f['levels']:
                private = [set(row['helpers'] + [row['target']]) for row in level]
                self.assertEqual(sum(map(len, private)), len(set().union(*private)))
                shared = set(f['address']) | {row['parent'] for row in level}
                self.assertTrue(all(gate[-1] not in shared
                                    for row in level for gate in row['word']))

    def test_native_edge_interfaces_used_in_two_unequal_orders(self):
        # Validate each substituted edge interface in the two unequal trees.
        # The matrices are at most six wires, not dense complete-tree matrices.
        for chunks in ((1, 2), (2, 1)):
            for length in chunks:
                for pattern in range(1 << length):
                    literals = [_form(bit + 1, constant=1 - ((pattern >> bit) & 1))
                                for bit in range(length)]
                    if length == 1:
                        width, target, phase = 4, 2, 3
                        affine = [('AA', _form(0), literals[0], target)]
                    else:
                        width, target, phase = 6, 3, 5
                        first = ('AA', _form(0), literals[0], 4)
                        second = ('AA', _form(4), literals[1], target)
                        affine = [first, second, first, second]
                    native = _private_native(affine, phase)
                    self.assertTrue(all(gate[-1] not in range(length + 1)
                                        for gate in native))
                    images = [basis ^ ((1 << target) if basis & 1
                                       and (basis >> 1 & ((1 << length) - 1)) == pattern
                                       else 0) for basis in range(1 << width)]
                    expected = np.eye(1 << width)[:, images]
                    np.testing.assert_allclose(_word_matrix(width, native), expected,
                                               atol=ATOL, rtol=0)
                    np.testing.assert_allclose(_word_matrix(width, _adjoint(native)),
                                               expected.conj().T, atol=ATOL, rtol=0)

    def test_wrong_root_order_and_missing_final_inverse_are_detected(self):
        for chunks in ((1, 2), (2, 1), (1, 2, 1)):
            f = _unequal_tree(chunks)
            initial = [{1 << wire} for wire in range(f['width'])]
            expected = _expected(initial, f['address'], f['output'],
                                 [1 << x for x in range(1 << sum(chunks))])
            wrong_path = f['forward'] + [('X', f['nodes'][0][0])] + _adjoint(f['forward'])
            wrong_echo = f['copy'] + wrong_path + f['copy'] + _adjoint(wrong_path)
            self.assertNotEqual(_symbolic(wrong_echo, initial), expected)
            missing = f['copy'] + f['path_flip'] + f['copy']
            result = _symbolic(missing, initial)
            self.assertNotEqual(result, expected)
            self.assertEqual(result[f['nodes'][0][0]], initial[f['nodes'][0][0]] ^ {0})

    def test_exact_ceiling_transitions_and_contracting_partitions(self):
        for bits in _lengths():
            chunks = _partition(bits)
            self.assertEqual(sum(chunks), bits)
            self.assertTrue(all(length > 0 for length in chunks))
            if bits:
                self.assertLessEqual(chunks[-1], 64)
            remaining = bits
            for length in chunks[:-1]:
                following = remaining - length
                target = (remaining + 2) ** 4
                # Independently certify the ceiling, including exact powers.
                self.assertLess(1 << (following - 1), target)
                self.assertLessEqual(target, 1 << following)
                self.assertGreater(remaining, 64)
                self.assertLessEqual(2 * following, remaining)
                if remaining >= 4096:
                    self.assertLessEqual((following + 2) ** 2, remaining + 2)
                    self.assertLessEqual(_ceil_log2(following + 2),
                                         _ceil_log2(remaining + 2) // 2)
                else:
                    self.assertLessEqual(following, 49)
                remaining = following
        self.assertEqual(_partition(0), [])
        self.assertEqual(_partition(64), [64])
        self.assertEqual(_partition(65), [40, 25])
        self.assertEqual(_partition(126), [98, 28])
        self.assertEqual(_partition(127), [98, 29])

    def test_normalized_count_pool_and_depth_have_exact_uniform_majorants(self):
        count_bound = 67 ** 3 + Fraction(4, 65)
        for bits in _lengths():
            chunks = _partition(bits)
            remaining, nonfinal, peak = bits, Fraction(0), Fraction(0)
            depth = 0
            for index, length in enumerate(chunks):
                following = remaining - length
                # Stage edges/2^bits = 2^(-following). The denominator's
                # bit length is logarithmic, even for a million address bits.
                normalized = Fraction((length + 3) ** 3, 1 << following)
                peak = max(peak, normalized)
                depth += _ceil_log2(length + 3)
                if index + 1 < len(chunks):
                    self.assertLessEqual(normalized, Fraction(2, remaining + 2))
                    nonfinal += normalized
                else:
                    self.assertLessEqual(normalized, 67 ** 3)
                remaining = following
            self.assertLessEqual(nonfinal, Fraction(4, 65))
            final = (chunks[-1] + 3) ** 3 if chunks else 0
            self.assertLessEqual(nonfinal + final, count_bound)
            self.assertLessEqual(peak, 67 ** 3)
            # Above 4096, the ceiling-log surrogate at least halves.
            # Below 4096, at most one nonfinal step and one final step
            # remain, costing at most ceil(log2 4098)+ceil(log2 67)=20.
            self.assertLessEqual(depth, 2 * _ceil_log2(bits + 2) + 20)
            if not bits:
                self.assertEqual((nonfinal, peak, depth), (0, 0, 0))


if __name__ == '__main__':
    unittest.main()
