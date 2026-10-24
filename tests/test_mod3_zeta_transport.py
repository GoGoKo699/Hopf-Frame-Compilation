"""Exact finite evidence for the mod-three zeta transport reduction.

Binary basis columns certify the full linear maps through p=6; the affine
indicator echo also keeps every bank input through p=3. Native fixtures
retain literal phases and actual inverses only for the small binary basis
changes. These identities alone do not charge address-controlled blocks.
"""
from functools import lru_cache
from itertools import product
import unittest

from tests.exact_matrices import CO, adjoint, eye, native_columns, zeros
from tests.test_entangling_selector_readout import _binary_image, _zeta_word
from tests.test_fixed_row_transport_readout import _translate


_OMEGA_POWERS = (1, 2, 3)  # a+b*omega is encoded as a+2*b; omega^2=omega+1.


def _gf4_multiply(first, second):
    result = 0
    while second:
        if second & 1:
            result ^= first
        second >>= 1
        first <<= 1
        if first & 4:
            first ^= 7
    return result


def _trace(value):
    return value ^ _gf4_multiply(value, value)


def _apply(columns, vector):
    result = 0
    while vector:
        lowest = vector & -vector
        result ^= columns[lowest.bit_length() - 1]
        vector ^= lowest
    return result


def _inverse(columns):
    size = len(columns)
    rows = [sum(((column >> row) & 1) << j
                for j, column in enumerate(columns)) | (1 << (size + row))
            for row in range(size)]
    for column in range(size):
        pivot = next(row for row in range(column, size) if rows[row] >> column & 1)
        rows[column], rows[pivot] = rows[pivot], rows[column]
        for row in range(size):
            if row != column and rows[row] >> column & 1:
                rows[row] ^= rows[column]
    return tuple(sum(((rows[row] >> (size + column)) & 1) << row
                     for row in range(size)) for column in range(size))


@lru_cache(None)
def _basis(p):
    signs = tuple((1,) + tail for tail in product((1, -1), repeat=p - 1))
    columns = []
    for sigma in signs:
        eigenvector = []
        for row in range(1 << p):
            value = 1
            for coordinate, sign in enumerate(sigma):
                if row >> coordinate & 1:
                    value = _gf4_multiply(value, _OMEGA_POWERS[sign % 3])
            eigenvector.append(value)
        for scalar in (1, 2):
            columns.append(sum(_trace(_gf4_multiply(scalar, value)) << row
                               for row, value in enumerate(eigenvector)))
    columns = tuple(columns)
    return signs, columns, _inverse(columns)


def _coordinate_a(vector, p, coordinate):
    result = 0
    for low in range(1 << p):
        if not low >> coordinate & 1:
            high = low | (1 << coordinate)
            first, second = (vector >> low) & 1, (vector >> high) & 1
            result |= (second << low) | ((first ^ second) << high)
    return result


def _block_power(pair, exponent):
    for _ in range(exponent % 3):
        first, second = pair & 1, (pair >> 1) & 1
        pair = second | ((first ^ second) << 1)
    return pair


def _block_action(vector, exponents):
    return sum(_block_power((vector >> (2 * block)) & 3, exponent) << (2 * block)
               for block, exponent in enumerate(exponents))


def _moved_zeta_columns(p, address):
    word = _zeta_word(p)
    return tuple(_translate(_binary_image(word, 1 << (row ^ address)), address, p)
                 for row in range(1 << p))


def _native_binary_word(columns):
    """Gaussian row operations, reversed to implement the supplied matrix."""
    size, elimination = len(columns), []
    rows = [sum(((column >> row) & 1) << j for j, column in enumerate(columns))
            for row in range(size)]

    def add(source, target):
        rows[target] ^= rows[source]
        elimination.append(('CX', source, target))

    for column in range(size):
        pivot = next(row for row in range(column, size) if rows[row] >> column & 1)
        if pivot != column:
            add(pivot, column)
            add(column, pivot)
            add(pivot, column)
        for row in range(size):
            if row != column and rows[row] >> column & 1:
                add(column, row)
    assert rows == [1 << row for row in range(size)]
    return list(reversed(elimination))


class Mod3ZetaTransportTests(unittest.TestCase):
    def test_trace_basis_simultaneously_reduces_every_coordinate(self):
        self.assertEqual([_trace(value) for value in range(4)], [0, 0, 1, 1])
        for p in range(1, 7):
            signs, columns, inverse = _basis(p)
            for column in range(1 << p):
                basis = 1 << column
                self.assertEqual(_apply(inverse, columns[column]), basis)
                self.assertEqual(_apply(columns, inverse[column]), basis)
                block, bit = divmod(column, 2)
                for coordinate in range(p):
                    actual = _apply(inverse, _coordinate_a(columns[column], p, coordinate))
                    expected = _block_power(1 << bit, signs[block][coordinate]) << (2 * block)
                    self.assertEqual(actual, expected)

            transformed_ones = _apply(inverse, (1 << (1 << p)) - 1)
            for block, sigma in enumerate(signs):
                pair = (transformed_ones >> (2 * block)) & 3
                self.assertEqual(pair, _OMEGA_POWERS[sum(sigma) % 3])
                self.assertNotEqual(pair, 0)
                self.assertEqual((inverse[0] >> (2 * block)) & 3,
                                 _OMEGA_POWERS[-sum(sigma) % 3])

    def test_all_addresses_factor_into_mod_three_blocks_on_full_binary_basis(self):
        for p in range(1, 7):
            signs, columns, inverse = _basis(p)
            word = _zeta_word(p)
            for address in range(1 << p):
                translated = _moved_zeta_columns(p, address)
                exponents = [sum(sign for i, sign in enumerate(sigma) if address >> i & 1)
                             for sigma in signs]
                for column in range(1 << p):
                    factored = _binary_image(word, 1 << column)
                    conjugated = columns[column]
                    for coordinate in range(p):
                        if address >> coordinate & 1:
                            factored = _coordinate_a(factored, p, coordinate)
                            conjugated = _coordinate_a(conjugated, p, coordinate)
                    self.assertEqual(translated[column], factored)
                    block, bit = divmod(column, 2)
                    expected = _block_power(1 << bit, exponents[block]) << (2 * block)
                    self.assertEqual(_apply(inverse, conjugated), expected)

    def test_indicator_echo_retains_arbitrary_banks_and_is_an_actual_involution(self):
        for p in range(1, 7):
            m, full = 1 << p, (1 << (1 << p)) - 1
            signs, basis, inverse = _basis(p)
            # Zero and all basis columns determine the complete affine map;
            # the smaller cases also exhaust every arbitrary bank string.
            banks = range(1 << m) if p <= 3 else (0, full, *(1 << row for row in range(m)))
            for address in range(m):
                columns = _moved_zeta_columns(p, address)
                exponents = [sum(sign for i, sign in enumerate(sigma) if address >> i & 1)
                             for sigma in signs]
                reverse = [-exponent for exponent in exponents]
                self.assertEqual(_apply(columns, full), 1 << address)
                for bank in banks:
                    middle = _apply(columns, bank)
                    self.assertEqual(_apply(columns, middle), bank)
                    self.assertEqual(_apply(columns, middle ^ full), bank ^ (1 << address))
                    # Direct block echo B X^{S e_0} B^dagger needs no
                    # controlled B, and no whole moved-zeta implementation.
                    state = _block_action(_apply(inverse, bank), reverse)
                    echoed = _block_action(state ^ inverse[0], exponents)
                    self.assertEqual(_apply(basis, echoed), bank ^ (1 << address))

    def test_weak_pauli_transport_does_not_require_full_conjugator_approximation(self):
        for p in (1, 2):
            m, full = 1 << p, (1 << (1 << p)) - 1
            for address in range(m):
                columns = _moved_zeta_columns(p, address)

                def forward(bank):
                    return _apply(columns, bank ^ 1)  # E=M_x X_0.

                def backward(bank):
                    return _apply(columns, bank) ^ 1  # E^dagger=X_0 M_x.

                # X_0 has a minus eigenvector |0>-|1>. Its E and M_x
                # images differ by sign, attaining operator distance 2.
                zeta_minus = {_apply(columns, 0): 1, _apply(columns, 1): -1}
                changed_minus = {forward(0): 1, forward(1): -1}
                self.assertEqual(changed_minus, {key: -value for key, value in zeta_minus.items()})
                for bank in range(1 << m):
                    self.assertEqual(backward(forward(bank)), bank)
                    self.assertEqual(forward(backward(bank)), bank)
                    for signal in (0, 1):
                        actual = forward(backward(bank) ^ (full if signal else 0))
                        self.assertEqual(actual, bank ^ ((1 << address) if signal else 0))
                if address == 0:
                    # The actual inverse is necessary in this fiber:
                    # applying E twice is not the identity.
                    self.assertNotEqual(forward(forward(0)), 0)

    def test_indicator_oracle_converse_returns_every_borrowed_column(self):
        # Complete p=2 oracle composition, not a native indicator emitter.
        # In matrix order A^r=U^[r=1] L^[r!=0] U^[r=2].
        p, m, q = 2, 4, 2
        signs, _, _ = _basis(p)
        residues = [[sum(sign for i, sign in enumerate(sigma) if address >> i & 1) % 3
                     for sigma in signs] for address in range(m)]
        factors = [('U', tuple(sum((row[j] == 2) << j for j in range(q))
                               for row in residues)),
                   ('L', tuple(sum((row[j] != 0) << j for j in range(q))
                               for row in residues)),
                   ('U', tuple(sum((row[j] == 1) << j for j in range(q))
                               for row in residues))]

        def extract(address, query, predicates, table, inverse=False):
            # E=I^dagger C_f I C_f; reverse every actual factor for E^dagger.
            calls = 0
            for gate in (('I', 'C', 'I', 'C') if inverse else ('C', 'I', 'C', 'I')):
                if gate == 'I':
                    query ^= 1 << address
                    calls += 1
                else:
                    predicates ^= _apply(table, query)
            return query, predicates, calls

        def controlled(bank, predicates, direction):
            for block in range(q):
                source, target = 2 * block, 2 * block + 1
                if direction == 'U':
                    source, target = target, source
                bank ^= ((predicates >> block & 1) * (bank >> source & 1)) << target
            return bank

        for address, original_query, original_predicates, original_bank in product(
                range(m), range(1 << m), range(1 << q), range(1 << (2 * q))):
            query, predicates, bank, calls = original_query, original_predicates, original_bank, 0
            for direction, table in factors:
                bank = controlled(bank, predicates, direction)
                query, predicates, used = extract(address, query, predicates, table)
                calls += used
                bank = controlled(bank, predicates, direction)
                query, predicates, used = extract(address, query, predicates, table, inverse=True)
                calls += used
                self.assertEqual((query, predicates), (original_query, original_predicates))
            self.assertEqual(calls, 12)
            self.assertEqual(bank, _block_action(original_bank, residues[address]))

    def test_small_native_basis_change_all_columns_and_actual_inverse(self):
        for p in (1, 2):
            m = 1 << p
            _, columns, inverse = _basis(p)
            word = _native_binary_word(inverse)
            reverse = list(reversed(word))
            expected, backward = zeros(1 << m), zeros(1 << m)
            for bank in range(1 << m):
                expected[_apply(inverse, bank)][bank] = CO
                backward[_apply(columns, bank)][bank] = CO
            actual = native_columns(m, word)
            self.assertEqual(actual, expected)
            self.assertEqual(native_columns(m, reverse), backward)
            self.assertEqual(backward, adjoint(actual))
            self.assertEqual(native_columns(m, word + reverse), eye(1 << m))


if __name__ == '__main__':
    unittest.main()
