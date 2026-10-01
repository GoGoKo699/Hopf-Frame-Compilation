"""Finite exact-algebra checks of native-tree antichain updates.

Complex noncommuting local rotations test the forest factorization and full
prefix-cylinder packing. Literal Clifford+T dirty-Fredkin words are checked on
every input, including an occupied dirty helper and an inactive selector.
Matrices have dimension at most 64; no asymptotic cost is inferred here.
"""
from __future__ import annotations

import unittest

import numpy as np

from tests.test_one_clean_compiler import _controlled_x
from tests.test_operator_source_compiler import (
    _adjoint, _controlled_swap, _expand_toffolis, _word_matrix,
)
from tests.test_tree_residual_structure import ATOL, _addressed_frame, _su2


CASES = ((1, ('',)), (2, ('0', '1')), (3, ('0', '10')),
         (4, ('0', '10')), (4, ('00', '01', '1')), (4, ('000', '01', '1')))


def _node(prefix):
    return int('1' + prefix, 2)


def _prefix(node):
    return bin(node)[3:]


def _pair(height, prefix):
    first = int(prefix or '0', 2) << (height - len(prefix))
    return [first, first | (1 << (height - len(prefix) - 1))]


def _fixture(height, marked):
    size = 1 << height
    coarse = {node: _su2(0.13 + 0.03 * node, -0.07 * node, 0.11 * (node + 1))
              for node in range(1, size)}
    target = {node: word.copy() for node, word in coarse.items()}
    for prefix in marked:
        node = _node(prefix)
        target[node] = _su2(-0.21 + 0.02 * node, 0.19 + 0.04 * node, -0.05 * node)
    return coarse, target


def _local_product(height, words, nodes):
    product = np.eye(1 << height, dtype=complex)
    instructions = []
    for node in sorted(nodes):
        rows = _pair(height, _prefix(node))
        product[rows] = words[node] @ product[rows]
        instructions.append((rows, words[node]))
    inverse = np.eye(1 << height, dtype=complex)
    for rows, word in reversed(instructions):
        inverse[rows] = word.conj().T @ inverse[rows]
    return product, inverse


def _forest_and_corrections(height, marked, coarse, target):
    descendants = [node for node in coarse
                   if any(_prefix(node).startswith(prefix) and len(_prefix(node)) > len(prefix)
                          for prefix in marked)]
    forest, forest_inverse = _local_product(height, coarse, descendants)
    corrections = {_node(prefix): target[_node(prefix)] @ coarse[_node(prefix)].conj().T
                   for prefix in marked}
    middle, _ = _local_product(height, corrections, corrections)
    return forest, forest_inverse, middle, corrections


def _packing_images(height, marked):
    result = []
    for original in range(1 << height):
        value = original
        for prefix in marked:
            depth = len(prefix)
            matches = not depth or value >> (height - depth) == int(prefix, 2)
            first_bit, last_bit = height - depth - 1, 0
            if matches and ((value >> first_bit) ^ (value >> last_bit)) & 1:
                value ^= (1 << first_bit) | (1 << last_bit)
        result.append(value)
    return np.array(result)


def _native_packing_word(height, marked, dirty, selector=None):
    word = []
    for depth in sorted({len(prefix) for prefix in marked}):
        first, last = height - depth - 1, 0
        if first == last:
            continue
        query = []
        for prefix in marked:
            if len(prefix) != depth:
                continue
            controls = [height - index - 1 for index in range(depth)]
            negatives = [('X', bit) for bit, value in zip(controls, prefix) if value == '0']
            if selector is not None:
                controls.append(selector)
            # For the three-control case the swap endpoint is borrowed only
            # inside this query and returned before its Fredkin is executed.
            toggle = _controlled_x(controls, dirty, helper=first)
            query += negatives + toggle + _adjoint(negatives)
        swap = _controlled_swap(dirty, first, last)
        word += query + swap + _adjoint(query) + swap
    return _expand_toffolis(word)


class AntichainCompilerTests(unittest.TestCase):
    def test_strict_descendant_forest_factorization_for_complex_mixed_depth_updates(self):
        for height, marked in CASES:
            with self.subTest(height=height, marked=marked):
                coarse, target = _fixture(height, marked)
                forest, inverse, middle, _ = _forest_and_corrections(height, marked, coarse, target)
                c_frame = _addressed_frame(height, coarse)
                w_frame = _addressed_frame(height, target)
                actual = forest @ middle @ inverse @ c_frame
                np.testing.assert_allclose(actual, w_frame, atol=ATOL, rtol=0)
                np.testing.assert_allclose(inverse @ forest, np.eye(1 << height),
                                           atol=ATOL, rtol=0)
                np.testing.assert_allclose(inverse, forest.conj().T, atol=ATOL, rtol=0)
                np.testing.assert_allclose(actual.conj().T @ actual, np.eye(1 << height),
                                           atol=ATOL, rtol=0)

    def test_full_prefix_cylinder_packing_is_one_last_bit_multiplexor(self):
        for height, marked in CASES:
            with self.subTest(height=height, marked=marked):
                coarse, target = _fixture(height, marked)
                forest, inverse, middle, corrections = _forest_and_corrections(
                    height, marked, coarse, target)
                images = _packing_images(height, marked)
                np.testing.assert_array_equal(np.sort(images), np.arange(1 << height))
                np.testing.assert_array_equal(images[images], np.arange(1 << height))
                packing = np.eye(1 << height)[np.argsort(images)]
                multiplexor = np.eye(1 << height, dtype=complex)
                addresses = []
                for prefix in marked:
                    address = int(prefix or '0', 2) << (height - len(prefix) - 1)
                    addresses.append(address)
                    packed_pair = [2 * address, 2 * address + 1]
                    np.testing.assert_array_equal(images[_pair(height, prefix)], packed_pair)
                    multiplexor[np.ix_(packed_pair, packed_pair)] = corrections[_node(prefix)]
                self.assertEqual(len(addresses), len(set(addresses)))
                np.testing.assert_allclose(packing @ middle @ packing.conj().T, multiplexor,
                                           atol=ATOL, rtol=0)
                # Equality above includes every unprogrammed last-bit row and
                # all states moved temporarily by a full-cylinder swap.
                actual = (forest @ packing.conj().T @ multiplexor @ packing @ inverse
                          @ _addressed_frame(height, coarse))
                np.testing.assert_allclose(actual, _addressed_frame(height, target),
                                           atol=ATOL, rtol=0)

    def test_native_dirty_fredkin_echo_preserves_the_address_and_optional_selector(self):
        for height, marked in CASES:
            for selected in (False, True):
                with self.subTest(height=height, marked=marked, selected=selected):
                    dirty = height
                    selector = height + 1 if selected else None
                    width = height + 1 + int(selected)
                    word = _native_packing_word(height, marked, dirty, selector)
                    actual = _word_matrix(width, word)
                    logical_images = _packing_images(height, marked)
                    images = []
                    for basis in range(1 << width):
                        logical = basis & ((1 << height) - 1)
                        active = selector is None or (basis >> selector) & 1
                        image = logical_images[logical] if active else logical
                        images.append((basis & ~((1 << height) - 1)) | image)
                    expected = np.eye(1 << width)[np.argsort(images)]
                    np.testing.assert_allclose(actual, expected, atol=ATOL, rtol=0)
                    inverse = _word_matrix(width, _adjoint(word))
                    np.testing.assert_allclose(inverse @ actual, np.eye(1 << width),
                                               atol=ATOL, rtol=0)
                    for depth in {len(prefix) for prefix in marked}:
                        # During each membership-controlled swap, the prefix
                        # used by its query is unchanged for every logical input.
                        same_depth = tuple(prefix for prefix in marked if len(prefix) == depth)
                        layer_images = _packing_images(height, same_depth)
                        for basis, image in enumerate(layer_images):
                            self.assertEqual(basis >> (height - depth), image >> (height - depth))

    def test_comparable_updates_and_omitted_forest_do_not_satisfy_this_word(self):
        # This is a countercheck for the stated factorization, not a claim that
        # comparable updates admit no other efficient implementation.
        height, marked = 3, ('', '0')
        coarse, target = _fixture(height, marked)
        forest, inverse, middle, _ = _forest_and_corrections(height, marked, coarse, target)
        c_frame = _addressed_frame(height, coarse)
        w_frame = _addressed_frame(height, target)
        self.assertGreater(np.linalg.norm(forest @ middle @ inverse @ c_frame - w_frame, ord=2), 0.01)
        packed_addresses = [int(prefix or '0', 2) << (height - len(prefix) - 1) for prefix in marked]
        self.assertEqual(len(set(packed_addresses)), 1)
        height, marked = 4, ('0', '10')
        coarse, target = _fixture(height, marked)
        _, _, middle, _ = _forest_and_corrections(height, marked, coarse, target)
        self.assertGreater(np.linalg.norm(
            middle @ _addressed_frame(height, coarse) - _addressed_frame(height, target), ord=2), 0.01)

    def test_dirty_interpreter_echoes_each_reflection_in_a_noncommuting_word(self):
        # Low to high: target, unchanged prefix, unchanged suffix, dirty
        # predicate bit beta, and dirty bank bit. Neither dirty bit is prepared.
        target, prefix, suffix, beta, bank, width = 0, 1, 2, 3, 4, 5
        h_gate = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2)
        z_gate = np.diag([1, -1]).astype(complex)
        t_gate = np.diag([1, np.exp(1j * np.pi / 4)])
        reflection_matrices = (h_gate, z_gate, t_gate @ h_gate @ t_gate.conj().T)
        tables = ((1, 1), (1, 0), (0, 1))

        def controlled_reflection(which):
            controlled_z = [('H', target), ('CCX', beta, bank, target), ('H', target)]
            if which == 1:
                return controlled_z
            # V=Ry(pi/8), up to a scalar, and V Z V^dagger=H. The scalar
            # cancels between the actual word and its actual inverse.
            v_word = [('SDG', target), ('H', target), ('T', target),
                      ('H', target), ('S', target)]
            controlled_h = _adjoint(v_word) + controlled_z + v_word
            if which == 0:
                return controlled_h
            return [('TDG', target)] + controlled_h + [('T', target)]

        def loader(table):
            if table == (1, 1):
                return [('X', bank)]
            if table == (1, 0):
                return [('X', prefix), ('CX', prefix, bank), ('X', prefix)]
            return [('CX', prefix, bank)]

        interpreters = []
        for which, table in enumerate(tables):
            query = loader(table)
            controlled = controlled_reflection(which)
            # This inner echo returns the dirty bank and implements the
            # prefix-table reflection controlled by beta, on all beta inputs.
            interpreters.append(query + controlled + _adjoint(query) + controlled)
        predicate = [('X', suffix), ('CX', suffix, beta), ('X', suffix)]
        word = []
        for interpreter in interpreters:
            word += predicate + interpreter + _adjoint(predicate) + interpreter
        word = _expand_toffolis(word)
        actual = _word_matrix(width, word)
        expected = np.zeros_like(actual)
        for assignment in range(1 << (width - 1)):
            first = assignment << 1
            row = (first >> prefix) & 1
            local = np.eye(2, dtype=complex)
            if not (first >> suffix) & 1:
                for table, reflection in zip(tables, reflection_matrices):
                    if table[row]:
                        local = reflection @ local
            indices = [first, first | 1]
            expected[np.ix_(indices, indices)] = local
        np.testing.assert_allclose(actual, expected, atol=ATOL, rtol=0)
        inverse = _word_matrix(width, _adjoint(word))
        np.testing.assert_allclose(inverse @ actual, np.eye(1 << width), atol=ATOL, rtol=0)

        # Echoing the whole noncommuting word instead would apply its square
        # when h=0 and beta=1; a general native word is not an involution.
        whole_word = [gate for interpreter in interpreters for gate in interpreter]
        wrong = _word_matrix(width, _expand_toffolis(
            predicate + whole_word + _adjoint(predicate) + whole_word))
        inactive_dirty_one = [basis for basis in range(1 << width)
                              if (basis >> suffix) & 1 and (basis >> beta) & 1]
        self.assertGreater(np.linalg.norm(
            (wrong - expected)[:, inactive_dirty_one], ord=2), 1)


if __name__ == '__main__':
    unittest.main()
