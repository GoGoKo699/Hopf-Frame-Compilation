"""Exact finite tests for the read-only zero predicate and simultaneous row filter."""
import unittest
from collections import Counter, deque
from itertools import permutations, product

import numpy as np


def parity(x):
    return x.bit_count() & 1


def mat_vec(columns, x):
    out = 0
    for i, c in enumerate(columns):
        if x >> i & 1:
            out ^= c
    return out


def vec_mat(row, columns):
    return sum(parity(row & c) << i for i, c in enumerate(columns))


def transvections():
    h_group = []
    for columns in permutations(range(1, 8), 3):
        if len({mat_vec(columns, x) for x in range(8)}) == 8:
            h_group.append(columns)
    counts = Counter()
    for columns in h_group:
        inverse = tuple(next(x for x in range(8) if mat_vec(columns, x) == 1 << i)
                        for i in range(3))
        for i, j in permutations(range(3), 2):
            counts[columns[i], vec_mat(1 << j, inverse)] += 1
    assert len(h_group) == 168
    assert len(counts) == 21 and set(counts.values()) == {48}
    assert all(not parity(v & row) for v, row in counts)
    return tuple(counts)


def mod(f, modulus):
    while f.bit_length() >= modulus.bit_length():
        f ^= modulus << (f.bit_length() - modulus.bit_length())
    return f


def times(u, v):
    return mat_vec(tuple(u << i for i in range(v.bit_length())), v)


def times_mod(u, v, f):
    return mod(times(u, v), f)


def evaluate(poly, a, f):
    val = 0
    for i in reversed(range(poly.bit_length())):
        val = times_mod(val, a, f) ^ ((poly >> i) & 1)
    return val


def minimal_poly(a, f):
    h = f.bit_length() - 1
    return next(t for t in range(2, 1 << (h+1)) if evaluate(t, a, f) == 0)


def nonzero_orbits(pairs):
    f, h = 0b1011, 3  # x^3+x+1
    zero = (0, 0, 0)
    for a in (2, 4, 6):
        field_table = tuple(times_mod(a, v, f) for v in range(1 << h))
        reached = {(1, 0, 0)}
        queue = deque(reached)
        while queue:
            values = queue.popleft()
            syndromes = tuple(mat_vec(values, j) for j in range(8))
            for column, row in pairs:
                added = field_table[syndromes[row]]
                image = tuple(v ^ (added if column >> i & 1 else 0)
                              for i, v in enumerate(values))
                if image not in reached:
                    reached.add(image)
                    queue.append(image)
        assert len(reached) == 511 and zero not in reached
    return 3 * 511


def independent_crt_lifts():
    f = 0b100101  # x^5+x^2+1
    h = f.bit_length()-1
    a, b = next((a,b) for a,b in permutations(range(2,1 << h),2)
                if len({minimal_poly(v, f) for v in (a,b,a^b)}) == 3
                and all(minimal_poly(v,f).bit_length() == h+1 for v in (a,b,a^b)))
    samples = (a,b,a^b)
    polys = tuple(minimal_poly(x, f) for x in samples)
    cases = 0
    for i,x in enumerate(samples):
        separator = 0b10  # t imposes kernel at origin
        for j, poly in enumerate(polys):
            if j != i:
                separator = times(separator, poly)
        for desired in range(1 << h):
            q = next(q for q in range(1 << h)
                     if evaluate(times(separator,q), x, f) == desired)
            lift = times(separator,q)
            assert evaluate(lift,0,f) == 0
            assert [evaluate(lift,u,f) for u in samples] == [
                desired if j==i else 0 for j in range(len(samples))]
            cases += 1
    return (a,b), cases


def joint_rank_profiles():
    # Diagonal projectors with several local ranks, including zero and full.
    dims = [2, 3, 4, 2]
    assignments = list(product(*(range(d) for d in dims)))
    count = 0
    for ranks in product(*(range(d+1) for d in dims)):
        single_defects = [0]*len(dims)
        for word in assignments:
            failures = [i for i,(u,r) in enumerate(zip(word,ranks)) if u >= r]
            base = not failures
            increments = [all(i==j for i in failures) and not base for j in range(len(dims))]
            assert sum(increments) <= 1
            for j, v in enumerate(increments):
                single_defects[j] += v
        assert sum(single_defects) <= len(assignments)
        count += 1
    return count



class ReadOnlyFieldEqualityTest(unittest.TestCase):
    """Finite field checks. The dimension-independent gap is analytical."""

    def test_parameter_orbits(self):
        self.assertEqual(nonzero_orbits(transvections()), 3 * 511)

    def test_independent_crt_lifts(self):
        _, number = independent_crt_lifts()
        self.assertEqual(number, 3 * 32)

    def test_joint_projector_ranks(self):
        self.assertEqual(joint_rank_profiles(), 3 * 4 * 5 * 3)

    def test_sparse_echo_with_coherent_comparator_error(self):
        """Fixed address blocks, arbitrary private targets and true/false row leakage."""
        x = np.array([[0., 1.], [1., 0.]])
        z = np.diag([1., -1.])
        theta = 0.03
        tilt = np.array([[np.cos(theta), -np.sin(theta)],
                         [np.sin(theta), np.cos(theta)]])
        error = np.linalg.norm(tilt - np.eye(2), 2)
        rows = 4
        for address in range(rows):
            for occupation in product((0, 1), repeat=rows):
                if sum(occupation) > 2:
                    continue
                bank = np.ones((1, 1))
                for row, bit in enumerate(occupation):
                    comparator = tilt @ (x if row == address else np.eye(2))
                    mask = z if bit else np.eye(2)
                    echo = comparator.T @ mask @ comparator @ mask
                    bank = np.kron(bank, echo)
                desired = (-1)**occupation[address] * np.eye(1 << rows)
                self.assertLessEqual(np.linalg.norm(bank - desired, 2),
                                     2 * sum(occupation) * error + 2e-14)
