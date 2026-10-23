"""Exact small affine co-moving singleton transport and readout fixtures.

These are logical identities, not an affine native-word resource ledger.
The singleton bank is retained coherently; no sampled row or clean selector
is substituted for the arbitrary-input contract.
"""
from fractions import Fraction
from itertools import product
import unittest

from tests.exact_matrices import CI, CO, CZ, CQ2, Q2
from tests.test_co_moving_lookup_stage import _linear


_M = 4
_ROOT_HALF = CQ2(Q2(0, Fraction(1, 2)))
_I_POWERS = (CO, CI, -CO, -CI)


def _relabel(selectors, images):
    return sum(((selectors >> u) & 1) << images[u] for u in range(_M))


def _add(*terms):
    result = {}
    for scalar, state in terms:
        for key, value in state.items():
            result[key] = result.get(key, CZ) + scalar * value
    return {key: value for key, value in result.items() if value != CZ}


def _embed(state):
    """J on relative (z,q) states, retaining (x,singleton-row,selectors)."""
    result = {}
    for (z, selectors), value in state.items():
        for v in range(_M):
            translated = _relabel(selectors, [u ^ v for u in range(_M)])
            result[z ^ v, v, translated] = value * Fraction(1, 2)
    return result


def _stage(state, pairs, offset, coordinates, inverse=False):
    images = [_linear(u, pairs) ^ offset for u in range(_M)]
    result = {}
    for (x, row, selectors), value in state.items():
        h = (selectors >> row) & 1
        n = h * sum(((x ^ row) >> c) & 1 for c in coordinates)
        if inverse:
            n = -n
        key = ((images[x], images[row], _relabel(selectors, images))
               if h else (x, row, selectors))
        result[key] = value * _I_POWERS[n % 4]
    return result


def _relative_stage(state, pairs, coordinates, inverse=False):
    result = {}
    for (z, selectors), value in state.items():
        h = selectors & 1
        n = h * sum((z >> c) & 1 for c in coordinates)
        if inverse:
            n = -n
        key = ((_linear(z, pairs),
                _relabel(selectors, [_linear(u, pairs) for u in range(_M)]))
               if h else (z, selectors))
        result[key] = value * _I_POWERS[n % 4]
    return result


def _flip(state):
    return {key[:-1] + (key[-1] ^ 15,): value for key, value in state.items()}


def _hadamard(state):
    result = {}
    for key, value in state.items():
        for selectors in range(16):
            output = key[:-1] + (selectors,)
            coefficient = Fraction((-1) ** ((selectors & key[-1]).bit_count()), 4)
            result[output] = result.get(output, CZ) + value * coefficient
    return {key: value for key, value in result.items() if value != CZ}


def _inner(left, right):
    return sum((value.conj() * right.get(key, CZ)
                for key, value in left.items()), CZ)


def _accepted(state, first_offset, second_offset):
    def a(vector):
        return _stage(vector, ((0, 1),), first_offset, (0,))

    def a_prime(vector):
        return _flip(a(_flip(vector)))

    def b(vector):
        return _hadamard(_stage(_hadamard(vector), ((1, 0),), second_offset, (1,)))

    def b_prime(vector):
        vector = _flip(_hadamard(vector))
        vector = _stage(vector, ((1, 0),), second_offset, (1,))
        return _hadamard(_flip(vector))

    numerator = _add((CO, b(a(state))), (CI, b(a_prime(state))),
                     (-CI, b_prime(a(state))), (-CO, b_prime(a_prime(state))))
    return _add((CQ2(Q2(0, Fraction(1, 4))), numerator))


class AffineSingletonTransportTests(unittest.TestCase):
    def test_two_bit_affine_intertwiner_and_adjoint_all_columns(self):
        for pairs in ((), ((0, 1),), ((1, 0),)):
            unchanged = [c for c in range(2) if all(c != t for _, t in pairs)]
            for offset in range(_M):
                if _linear(offset, pairs) != offset:
                    continue
                for mask in range(1 << len(unchanged)):
                    coordinates = tuple(c for j, c in enumerate(unchanged)
                                        if (mask >> j) & 1)
                    for z, selectors, inverse in product(range(_M), range(16), (False, True)):
                        state = {(z, selectors): CO}
                        actual = _stage(_embed(state), pairs, offset, coordinates, inverse)
                        expected = _embed(_relative_stage(state, pairs, coordinates, inverse))
                        self.assertEqual(actual, expected)

    def test_isometry_target_and_uniform_selector_interlayers(self):
        seen = set()
        for z, selectors in product(range(_M), range(16)):
            state = {(z, selectors): CO}
            embedded = _embed(state)
            self.assertFalse(seen.intersection(embedded))
            seen.update(embedded)
            self.assertEqual(_inner(embedded, embedded), CO)
            target = {key: value * Fraction(1, 2) for key, value in embedded.items()
                      if key[0] != key[1]}
            expected = {} if z == 0 else _add((CO * Fraction(1, 2), embedded))
            self.assertEqual(target, expected)
            self.assertEqual(_flip(embedded), _embed(_flip(state)))
            self.assertEqual(_hadamard(embedded), _embed(_hadamard(state)))

    def test_four_path_marked_three_quarters_survives_all_affine_offsets(self):
        for first_offset, second_offset in product((0, 2), (0, 1)):
            accepted = []
            for selector0 in (0, 1):
                relative = {(0, 2 | selector0): _ROOT_HALF,
                            (0, 8 | selector0): -_ROOT_HALF}
                state = _embed(relative)
                self.assertEqual(_inner(state, state), CO)
                self.assertTrue(all(x == row for x, row, _ in state))
                output = _accepted(state, first_offset, second_offset)
                self.assertTrue(all(x == row for x, row, _ in output))
                self.assertEqual(_inner(output, output), CO * Fraction(3, 4))
                accepted.append(output)
            # The cross term vanishes, so arbitrary coherent selector zero
            # (and an external reference) has the same acceptance probability.
            self.assertEqual(_inner(accepted[0], accepted[1]), CZ)


if __name__ == '__main__':
    unittest.main()
