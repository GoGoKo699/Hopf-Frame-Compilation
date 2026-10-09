"""Complete native words for queries whose output also supplies selectors.

The 24 matrix fixtures check every input column, literal XZ mask phases,
and actual inverses. Independent integer echoes check role changes and an
occupied core helper. The all-size resource claims follow from the proofs.
"""
from __future__ import annotations

import unittest

import numpy as np

from tests.test_operator_source_compiler import (
    ATOL, _adjoint, _basis_action, _expand_toffolis, _word_matrix,
)


def _traversal(address, targets, selectors, table, first_updates):
    """A complete dirty traversal with disjoint targets and selectors."""
    assert len(address) == len(selectors) and address
    assert not set(targets).intersection(selectors)
    assert not set(address).intersection(targets + selectors)
    gates = []

    def visit(level, row):
        if level == len(address):
            gates.extend(('CX', selectors[-1], target)
                         for bit, target in enumerate(targets)
                         if (table[row] >> bit) & 1)
            return
        for value in (0, 1):
            negate = [('X', address[level])] if value == 0 else []
            if level == 0:
                update = [('CX', address[0], selectors[0])] if first_updates else []
            else:
                update = [('CCX', selectors[level - 1], address[level], selectors[level])]
            gates.extend(negate + update)
            visit(level + 1, row | (value << level))
            gates.extend(_adjoint(update) + _adjoint(negate))

    visit(0, 0)
    return gates


def _query(address, targets, selectors, table):
    first = _traversal(address, targets, selectors, table, True)
    unwanted = _traversal(address, targets, selectors, table, False)
    return first + _adjoint(unwanted)


def _self_borrowed(table):
    # Low four bits are the dirty output, high two bits the unchanged address.
    address = [4, 5]
    return (_query(address, [0, 1], [2, 3], [word & 3 for word in table])
            + _query(address, [2, 3], [0, 1], [word >> 2 for word in table]))


def _literal_native(gates):
    result = []
    for gate in _expand_toffolis(gates):
        if gate[0] == 'X':
            bit = gate[1]
            result.extend([('H', bit), ('S', bit), ('S', bit), ('H', bit)])
        elif gate[0] == 'SDG':
            result.extend([('S', gate[1])] * 3)
        else:
            result.append(gate)
    assert {gate[0] for gate in result} <= {'H', 'S', 'CX', 'T', 'TDG'}
    return result


def _expected_pauli(xs, zs):
    """Literal P=X^x Z^z: Z acts on the initial dirty output word."""
    target = np.zeros((64, 64), dtype=complex)
    for basis in range(64):
        address, core = basis >> 4, basis & 15
        image = (address << 4) | (core ^ xs[address])
        target[image, basis] = (-1) ** ((core & zs[address]).bit_count())
    return target


def _row_echo(address, targets, selector, table):
    """Independent two-address row implementation, including dirty echo."""
    gates = []
    for row, mask in enumerate(table):
        neg = [('X', address[bit]) for bit in range(2) if not (row >> bit) & 1]
        predicate = [('CCX', *address, selector)]
        writes = [('CX', selector, target) for bit, target in enumerate(targets)
                  if (mask >> bit) & 1]
        gates += neg + predicate + writes + predicate + writes + _adjoint(neg)
    return gates


class SelfBorrowedQueryTests(unittest.TestCase):
    def test_full_native_xor_and_literal_pauli_queries(self):
        rng = np.random.default_rng(20261009)
        tables = [[0] * 4, [15] * 4, [0, 3, 12, 15], [15, 2, 9, 4]]
        tables += rng.integers(0, 16, size=(8, 4)).tolist()
        hadamards = [('H', bit) for bit in range(4)]
        for index, (xs, zs) in enumerate(zip(tables, reversed(tables))):
            xor = _self_borrowed(xs)
            # Conjugate the entire query on all four core bits, then apply X.
            pauli = hadamards + _self_borrowed(zs) + hadamards + xor
            for kind, word, expected in (
                ('xor', xor, _expected_pauli(xs, [0] * 4)),
                ('pauli', pauli, _expected_pauli(xs, zs)),
            ):
                with self.subTest(index=index, kind=kind):
                    native = _literal_native(word)
                    actual = _word_matrix(6, native)
                    np.testing.assert_allclose(actual, expected, atol=ATOL, rtol=0)
                    np.testing.assert_allclose(
                        _word_matrix(6, _literal_native(_adjoint(native))), expected.conj().T,
                        atol=ATOL, rtol=0)

    def test_independent_exact_overlap_and_return(self):
        for table in ((0, 15, 6, 9), (3, 12, 5, 10), (15, 14, 13, 11)):
            word = (_row_echo((0, 1), (2, 3), 4, [value & 3 for value in table])
                    + _row_echo((0, 1), (4, 5), 2, [value >> 2 for value in table]))
            for basis in range(64):
                address, core = basis & 3, basis >> 2
                expected = address | ((core ^ table[address]) << 2)
                self.assertEqual(_basis_action(basis, word), expected)
                self.assertEqual(_basis_action(expected, _adjoint(word)), basis)

    def test_all_paired_majorana_sign_tables(self):
        for signs in range(256):
            xs = [((signs >> (2 * bit)) ^ (signs >> (2 * bit + 1))) & 1
                  for bit in range(4)]
            zs = [((signs >> (2 * bit)) & 1) ^ (sum(xs[:bit]) & 1)
                  for bit in range(4)]
            for bit in range(4):
                even_sign = zs[bit] ^ (sum(xs[:bit]) & 1)
                self.assertEqual(even_sign, (signs >> (2 * bit)) & 1)
                self.assertEqual(even_sign ^ xs[bit], (signs >> (2 * bit + 1)) & 1)

    def test_source_core_helper_is_returned_on_every_input(self):
        # Controls 0,1,2; source target 3; occupied core helper 4.
        center = [('CCX', 0, 1, 4), ('CCX', 4, 2, 3),
                  ('CCX', 0, 1, 4), ('CCX', 4, 2, 3)]
        expected = np.zeros((32, 32), dtype=complex)
        for basis in range(32):
            image = basis ^ ((1 << 3) if basis & 7 == 7 else 0)
            self.assertEqual(_basis_action(basis, center), image)
            expected[image, basis] = 1
        np.testing.assert_allclose(
            _word_matrix(5, _literal_native(center)), expected, atol=ATOL, rtol=0)


if __name__ == '__main__':
    unittest.main()
