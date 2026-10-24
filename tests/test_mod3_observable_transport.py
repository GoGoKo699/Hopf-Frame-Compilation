"""Exact finite checks of observable gauge freedom and shared cleanup.

For r modulo three, E_r=U^[r=1] L^[r!=0] sends e_0 to A^r e_0;
at r=2, (A^r)^-1 E_r=U fixes e_0. The literal palindrome
B_G C_f I C_f B_G sends (T,b,z) to (T+e_x,b+f,G^f z).
Two palindromes give E_* with two indicator calls and b shifted by [r=2],
while E_* P_0 E_*^dagger restores all work. Three palindromes and one
bare indicator give the full batch with four calls and complete return.

The larger checks are phase-free ideal-indicator oracle compositions,
not native indicator emitters. A separate tiny native fixture retains
literal phases in K=C_1 D_1 C_0 P C_0^dagger D_1^dagger C_1^dagger:
K is Clifford, so E P E^dagger=C_2 D_2 K D_2^dagger C_2^dagger has
two T bands. Finite fixtures do not establish the analytic error gap.
"""
from itertools import product
import unittest

from tests.exact_matrices import adjoint, eye, multiply, native_columns
from tests.test_mod3_zeta_transport import _apply, _basis, _block_action, _block_power
from tests.test_two_layer_galois_return import _clifford_columns, _inverse


_PAIR_COLUMNS = ((1, 2), (2, 3), (3, 2))  # I, A, L.


def _representative(bank, residues, inverse=False):
    result = 0
    for block, residue in enumerate(residues):
        pair = (bank >> (2 * block)) & 3
        if inverse and residue == 1:
            pair = _block_power(pair, 2)
        else:
            pair = _apply(_PAIR_COLUMNS[residue], pair)
        result |= pair << (2 * block)
    return result


def _tables():
    signs, _, _ = _basis(2)
    residues = tuple(tuple(sum(sign for i, sign in enumerate(sigma) if x >> i & 1) % 3
                           for sigma in signs) for x in range(4))
    a, beta, c = (tuple(sum(int(predicate(r)) << j for j, r in enumerate(row))
                       for row in residues)
                  for predicate in (lambda r: r == 1, lambda r: r != 0, lambda r: r == 2))
    return residues, a, beta, c


def _controlled_shear(bank, predicates, direction):
    for block in range(2):
        source, target = 2 * block, 2 * block + 1
        if direction == 'U':
            source, target = target, source
        bank ^= ((predicates >> block & 1) * (bank >> source & 1)) << target
    return bank


def _palindrome(table, direction):
    return [('B', direction), ('C', table), ('I', False), ('C', table), ('B', direction)]


def _adjoint(word):
    # The ideal oracle is self-adjoint, but keep actual-call orientation
    # in the word rather than treating a composite source as involutory.
    return [('I', not argument) if name == 'I' else (name, argument)
            for name, argument in reversed(word)]


def _run(word, state):
    x, query, predicates, bank, work = state
    for name, argument in word:
        if name == 'I':
            query ^= 1 << x
        elif name == 'C':
            predicates ^= _apply(argument, query)
        elif name == 'B':
            bank = _controlled_shear(bank, predicates, argument)
        else:
            assert name == 'P'
            bank ^= argument
    return x, query, predicates, bank, work


def _columns():
    # All 4096 address/query/predicate/block columns, with both values of
    # an additional arbitrary private-work wire retained as a spectator.
    return product(range(4), range(16), range(4), range(16), range(2))


class Mod3ObservableTransportTests(unittest.TestCase):
    def test_two_shear_gauge_and_complete_joint_observable(self):
        for residue, pair in product(range(3), range(4)):
            image = _apply(_PAIR_COLUMNS[residue], pair)
            relative = _block_power(image, -residue)
            expected_relative = _apply((1, 3), pair) if residue == 2 else pair
            self.assertEqual(relative, expected_relative)
            self.assertEqual(_apply(_PAIR_COLUMNS[residue], 1), _block_power(1, residue))
        self.assertNotEqual(_apply(_PAIR_COLUMNS[2], 2), _block_power(2, 2))

        residues, _, _, _ = _tables()
        p0 = 0b0101
        for x, bank in product(range(4), range(16)):
            row = residues[x]
            forward = _representative(bank, row)
            self.assertEqual(_representative(forward, row, inverse=True), bank)
            actual = _representative(_representative(bank, row, inverse=True) ^ p0, row)
            expected = bank ^ _block_action(p0, row)
            self.assertEqual(actual, expected)

    def test_two_query_source_and_actual_adjoint_echo_return_all_work(self):
        residues, a, beta, c = _tables()
        first, second = _palindrome(beta, 'L'), _palindrome(a, 'U')
        source = first + second
        echo = _adjoint(source) + [('P', 0b0101)] + source
        self.assertEqual(sum(name == 'I' for name, _ in source), 2)
        self.assertEqual([inverse for name, inverse in echo if name == 'I'],
                         [True, True, False, False])
        for state in _columns():
            x, query, predicates, bank, work = state
            middle = _run(first, state)
            self.assertEqual(middle, (x, query ^ (1 << x), predicates ^ beta[x],
                                      _controlled_shear(bank, beta[x], 'L'), work))
            actual = _run(second, middle)
            self.assertEqual(actual, (x, query, predicates ^ c[x],
                                      _representative(bank, residues[x]), work))
            self.assertEqual(_run(_adjoint(source), actual), state)
            self.assertEqual(_run(echo, state),
                             (x, query, predicates, bank ^ _block_action(0b0101, residues[x]), work))
        witness = (2, 9, 1, 6, 1)
        self.assertNotEqual(_run(source, witness)[2], witness[2])
        self.assertNotEqual(_run(source, _run(source, (1, 0, 0, 1, 1))), (1, 0, 0, 1, 1))

    def test_four_query_full_batch_returns_every_dirty_column(self):
        residues, a, beta, c = _tables()
        word = (_palindrome(c, 'U') + _palindrome(beta, 'L')
                + _palindrome(a, 'U') + [('I', False)])
        reverse = _adjoint(word)
        self.assertEqual(sum(name == 'I' for name, _ in word), 4)
        self.assertTrue(all(inverse for name, inverse in reverse if name == 'I'))
        for state in _columns():
            x, query, predicates, bank, work = state
            actual = _run(word, state)
            self.assertEqual(actual, (x, query, predicates, _block_action(bank, residues[x]), work))
            self.assertEqual(_run(reverse, actual), state)
            self.assertEqual(_run(reverse, state),
                             (x, query, predicates, _block_action(bank, [-r for r in residues[x]]), work))

    def test_literal_native_two_band_echo_folds_to_two_bands(self):
        width = 2
        c0 = [('CX', 0, 1), ('S', 0)]  # C0 X0 C0^dagger = Y0 X1.
        c1 = [('H', 0), ('CX', 1, 0), ('S', 1)]
        c2 = [('CX', 0, 1), ('H', 1), ('S', 0)]
        d1, d2 = [('T', 0)], [('TDG', 0), ('T', 1)]
        pauli = [('H', 0), ('S', 0), ('S', 0), ('H', 0)]
        source = c0 + d1 + c1 + d2 + c2
        actual_source = native_columns(width, source)
        self.assertEqual(native_columns(width, _inverse(source)), adjoint(actual_source))

        # T Y T^dagger = exp(-i*pi/4) S Y. The six-gate Clifford
        # (SH)^3 is exp(i*pi/4) I; retain its inverse's actual scalar.
        scalar = _inverse([('H', 0), ('S', 0)] * 3)
        k_without_scalar = _inverse(c1) + _inverse(c0) + pauli + c0 + [('S', 0)] + c1
        k_word = _inverse(c1) + _inverse(c0) + pauli + c0 + [('S', 0)] + scalar + c1
        central_word = _inverse(c1) + _inverse(d1) + _inverse(c0) + pauli + c0 + d1 + c1
        k = native_columns(width, k_word)
        self.assertEqual(k, native_columns(width, central_word))
        self.assertNotEqual(k, native_columns(width, k_without_scalar))
        self.assertTrue(all(gate[0] not in ('T', 'TDG') for gate in k_word))
        self.assertEqual(adjoint(k), k)
        self.assertEqual(multiply(k, k), eye(1 << width))
        self.assertEqual(len(_clifford_columns(k, width)), 2 * width)

        folded = _inverse(c2) + _inverse(d2) + k_word + d2 + c2
        full = _inverse(source) + pauli + source
        self.assertEqual(native_columns(width, folded), native_columns(width, full))
        self.assertEqual(native_columns(width, _inverse(folded)),
                         adjoint(native_columns(width, folded)))
        for band in (d1, d2, _inverse(d2)):
            self.assertEqual(len({gate[1] for gate in band}), len(band))
        self.assertEqual(sum(gate[0] in ('T', 'TDG') for gate in folded), 2 * len(d2))


if __name__ == '__main__':
    unittest.main()
