"""Exact native and algebra checks for the balanced HRT source.

The construction uses no initialized ancilla. Small native words are
checked over Q(sqrt(2),i), including every dirty column and actual adjoints.
Larger checks use exact coefficient algebra or integer phases modulo eight;
they do not replace the arbitrary-width proof in BALANCED_ALGEBRA_PERMUTATION.
"""
from itertools import product
from random import Random
import unittest

from tests.exact_matrices import CO, CZ, adjoint, native_columns


def _inverse(word):
    names = {'H': 'H', 'CX': 'CX', 'T': 'TDG', 'TDG': 'T'}
    return [(names[g[0]], *g[1:]) for g in reversed(word)]


def _gray_grid(terms, private):
    """Terms are (Y input, X input, output); four private-target T bands."""
    bands = []
    word = []
    for step in range(4):
        band = [('T' if step % 2 == 0 else 'TDG', h) for h in private]
        bands.append(band)
        word += band
        word += [('CX', term[0] if step % 2 == 0 else term[2], h)
                 for term, h in zip(terms, private)]
    return word, bands


def _bilinear_phase(terms, private):
    forward, _ = _gray_grid(terms, private)
    linear = [('CX', term[1], h) for term, h in zip(terms, private)]
    return _inverse(forward) + linear + forward + _inverse(linear)


def _bilinear_word(terms, outputs, private):
    if not terms:
        return []
    hadamards = [('H', z) for z in outputs]
    return hadamards + _bilinear_phase(terms, private) + hadamards


def _phase_column(word, column):
    current, phase = column, 0
    for gate in word:
        if gate[0] in ('T', 'TDG'):
            phase += (1 if gate[0] == 'T' else -1) * ((current >> gate[1]) & 1)
        elif gate[0] == 'CX':
            current ^= ((current >> gate[1]) & 1) << gate[2]
        else:
            raise ValueError(gate)
    return current, phase % 8


def _ring_product(a, b):
    """Coefficient bits are indexed by subsets, including the empty subset."""
    result = 0
    while a:
        low_a = a & -a
        i = low_a.bit_length() - 1
        a ^= low_a
        rest = b
        while rest:
            low_b = rest & -rest
            j = low_b.bit_length() - 1
            rest ^= low_b
            if not i & j:
                result ^= 1 << (i | j)
    return result


def _shared_multiplier(pairs, b0, v, dirty):
    """Six bilinear calls, with the first shear concurrent with Q."""
    pairs = [list(pair) for pair in pairs]
    for pair in pairs:
        pair[0] ^= _ring_product(b0, pair[1])
        pair[1] ^= pair[0]
        pair[0] ^= _ring_product(b0, pair[1])
        pair[1] ^= pair[0]
    dirty ^= b0 ^ v ^ _ring_product(v, b0)
    for pair in pairs:
        pair[1] ^= _ring_product(pair[0], dirty)
    dirty ^= b0 ^ v ^ _ring_product(v, b0)
    for pair in pairs:
        pair[1] ^= _ring_product(pair[0], dirty)
    return [tuple(pair) for pair in pairs], dirty


def _small_native_multiplier():
    # t=0,N=1: x,y,v,d,h are all arbitrary. b0 is a structural zero.
    # The two zero products emit no gates; the two diagonal copies remain.
    diagonal = [('CX', 0, 1), ('CX', 0, 1)]
    qword = [('CX', 2, 3)]
    eword = _bilinear_word([(0, 3, 1)], [1], [4])
    return diagonal + qword + eword + _inverse(qword) + _inverse(eword)


def _hrt_f(a, r):
    result = 1
    for subset in range(1, 1 << r):
        if a >> subset & 1:
            result ^= _ring_product(result, 1 << subset)
    return result ^ 1


def _hrt_k(a, r):
    occupied = [subset for subset in range(1, 1 << r) if a >> subset & 1]
    result = a
    for position, left in enumerate(occupied):
        for right in occupied[position + 1:]:
            if not left & right:
                result ^= 1 << (left | right)
    return result


def _merge(a, r, lo, mid, hi, dirty):
    """Physical coefficient update through N copies of the proved multiplier.

    Missing right-child coefficients are structural zero terms. Physical
    coefficients whose minimum labels are outside this node are preserved.
    """
    m = 1 << (r - mid - 1)
    pairs = []
    for left_subset in range(1, 1 << (mid - lo)):
        left = left_subset << lo
        x = y = 0
        for suffix in range(m):
            subset = left | (suffix << (mid + 1))
            x |= ((a >> subset) & 1) << suffix
            y |= ((a >> (subset | (1 << mid))) & 1) << suffix
        pairs.append((x, y))
    b0 = v = 0
    for suffix in range(m):
        subset = suffix << (mid + 1)
        if subset and (subset & -subset).bit_length() - 1 < hi:
            b0 |= ((a >> subset) & 1) << suffix
        # This coefficient has minimum label mid, inside the right child.
        v |= ((a >> (subset | (1 << mid))) & 1) << suffix
    transformed, returned = _shared_multiplier(pairs, b0, v, dirty)
    assert returned == dirty
    for left_subset, (x, y) in enumerate(transformed, 1):
        left = left_subset << lo
        for suffix in range(m):
            subset = left | (suffix << (mid + 1))
            for target, value in ((subset, (x >> suffix) & 1),
                                  (subset | (1 << mid), (y >> suffix) & 1)):
                a ^= (((a >> target) & 1) ^ value) << target
    return a


def _balanced(a, r, inverse=False, dirty_seed=0):
    def visit(state, lo, hi):
        if hi - lo == 1:
            return state
        mid = (lo + hi) // 2
        dirty_size = 1 << (1 << (r - mid - 1))
        dirty = (dirty_seed + 13 * lo + 7 * mid + hi) % dirty_size
        if inverse:
            state = _merge(state, r, lo, mid, hi, dirty)
            state = visit(state, mid, hi)
            return visit(state, lo, mid)
        state = visit(state, lo, mid)
        state = visit(state, mid, hi)
        return _merge(state, r, lo, mid, hi, dirty)
    return visit(a, 0, r)


def _tree_nodes(r):
    nodes = []
    def visit(lo, hi, level):
        if hi - lo == 1:
            return
        mid = (lo + hi) // 2
        n = (1 << (mid - lo)) - 1
        t = r - mid - 1
        nodes.append((lo, mid - lo, level, n, 1 << t, 3 ** t))
        visit(lo, mid, level + 1)
        visit(mid, hi, level + 1)
    visit(0, r, 0)
    return nodes


class BalancedAlgebraPermutationTests(unittest.TestCase):
    def test_gray_relative_phase_and_physical_band_ledger(self):
        word, bands = _gray_grid([(0, 99, 1)], [2])
        self.assertEqual(sum(g[0] in ('T', 'TDG') for g in word), 4)
        self.assertEqual(sum(g[0] == 'CX' for g in word), 4)
        self.assertEqual(len(bands), 4)
        for column in range(8):
            expected = ((4 * ((column >> 2) & 1) - 2)
                        * (column & 1) * ((column >> 1) & 1)) % 8
            self.assertEqual(_phase_column(word, column), (column, expected))
            self.assertEqual(_phase_column(word + _inverse(word), column),
                             (column, 0))

    def test_sparse_bilinear_phase_on_all_dirty_columns(self):
        fixtures = [
            # Product in F2[e]/(e²), all three nonzero tensor terms.
            ([(0, 2, 4), (0, 3, 5), (1, 2, 5)], [6, 7, 8], 9),
            # The cardinality baseline may reuse an unchanged source wire
            # as Y in one term and X in another; outputs remain separate.
            ([(0, 1, 3), (0, 2, 3), (1, 2, 3)], [4, 5, 6], 7),
        ]
        for terms, private, width in fixtures:
            phase = _bilinear_phase(terms, private)
            self.assertEqual(sum(g[0] in ('T', 'TDG') for g in phase), 8 * 3)
            self.assertEqual(sum(g[0] == 'CX' for g in phase), 10 * 3)
            _, bands = _gray_grid(terms, private)
            for band in bands:
                self.assertEqual([g[1] for g in band], private)
            for column in range(1 << width):
                target = 4 * (sum(((column >> y) & 1) * ((column >> x) & 1)
                                 * ((column >> z) & 1) for y, x, z in terms) % 2)
                self.assertEqual(_phase_column(phase, column), (column, target))
                self.assertEqual(_phase_column(phase + _inverse(phase), column),
                                 (column, 0))

    def test_first_stage_shared_source_native_schedule(self):
        # One term from A1 and one from Q. The shared b0 input is Y in
        # A1 and X in Q; both outputs and both private T targets differ.
        # This is the nonempty part of the t=1, N=1 concurrent stage.
        terms = [(0, 2, 3), (1, 0, 4)]
        private = [5, 6]
        word = _bilinear_word(terms, [3, 4], private)
        actual = native_columns(7, word)
        for column in range(128):
            b0 = column & 1
            output = (column ^ ((b0 * ((column >> 2) & 1)) << 3)
                      ^ ((b0 * ((column >> 1) & 1)) << 4))
            for row in range(128):
                self.assertEqual(actual[row][column],
                                 CO if row == output else CZ)
        self.assertEqual(native_columns(7, _inverse(word)), adjoint(actual))
        self.assertEqual(sum(g[0] in ('T', 'TDG') for g in word), 16)
        self.assertEqual(sum(g[0] not in ('T', 'TDG') for g in word), 24)
        bands, active = [], []
        for gate in word:
            if gate[0] in ('T', 'TDG'):
                active.append(gate[1])
            elif active:
                bands.append(active)
                active = []
        if active:
            bands.append(active)
        self.assertEqual(len(bands), 8)
        self.assertEqual([sorted(band) for band in bands], [private] * 8)

    def test_full_native_multiplier_and_actual_adjoint(self):
        word = _small_native_multiplier()
        actual = native_columns(5, word)
        self.assertEqual(sum(g[0] in ('T', 'TDG') for g in word), 16)
        self.assertEqual(sum(g[0] not in ('T', 'TDG') for g in word), 28)
        for column in range(32):
            output = column ^ (((column & 1) * ((column >> 2) & 1)) << 1)
            for row in range(32):
                self.assertEqual(actual[row][column], CO if row == output else CZ)
        self.assertEqual(native_columns(5, _inverse(word)), adjoint(actual))

    def test_shared_dirty_multiplier_coefficient_columns(self):
        for x, y, b0, v, dirty in product(range(4), range(4), (0, 2),
                                          range(4), range(4)):
            got, returned = _shared_multiplier([(x, y)], b0, v, dirty)
            expected = [(x ^ _ring_product(b0, x),
                         y ^ _ring_product(b0, y) ^ _ring_product(v, x))]
            self.assertEqual((got, returned), (expected, dirty))
        rng = Random(20261018)
        for n in (1, 2, 3, 7):
            for _ in range(64):
                pairs = [(rng.randrange(256), rng.randrange(256)) for _ in range(n)]
                b0, v, dirty = 2 * rng.randrange(128), rng.randrange(256), rng.randrange(256)
                got, returned = _shared_multiplier(pairs, b0, v, dirty)
                expected = [(x ^ _ring_product(b0, x),
                             y ^ _ring_product(b0, y) ^ _ring_product(v, x))
                            for x, y in pairs]
                self.assertEqual((got, returned), (expected, dirty))
                self.assertEqual(_shared_multiplier(got, b0, v, returned),
                                 (pairs, dirty))

    def test_balanced_tree_full_forward_and_reversed_inverse(self):
        rng = Random(20261018)
        for r in range(1, 6):
            cases = (range(1 << ((1 << r) - 1)) if r <= 3 else
                     [rng.randrange(1 << ((1 << r) - 1)) for _ in range(128)])
            for raw in cases:
                source = raw << 1
                expected = _hrt_f(source, r)
                forward = _balanced(source, r, dirty_seed=raw)
                self.assertEqual(forward, expected)
                self.assertEqual(_balanced(forward, r, inverse=True,
                                           dirty_seed=raw ^ 173), source)
                self.assertEqual(_balanced(source, r, inverse=True,
                                           dirty_seed=raw + 29), _hrt_k(source, r))
        # A literal query uses every coefficient as an address; nonsingleton
        # inputs cannot silently become initialized dirty work.
        self.assertNotEqual(_hrt_f(1 << 7, 3), 0)
        for r in range(1, 7):
            singleton_input = sum(1 << (1 << bit) for bit in range(r))
            self.assertEqual((_hrt_f(singleton_input, r) >> ((1 << r) - 1)) & 1, 1)

    def test_symbolic_count_depth_and_simultaneous_workspace(self):
        for r in range(2, 65):
            total_t = total_g = total_work = nk_total = 0
            levels = {}
            seen_pairs = set()
            for lo, left_length, level, n, m, k in _tree_nodes(r):
                self.assertNotIn((lo, left_length), seen_pairs)
                seen_pairs.add((lo, left_length))
                terms = (4 * n + 2) * k - (2 * n + 2) * m
                self.assertEqual(terms, 2 * n * (k - m) + 2 * (k - m) + 2 * n * k)
                total_t += 8 * terms
                cliffords = 10 * terms + 10 * n * m + 8 * m - 2
                self.assertEqual(cliffords,
                                 (40 * n + 20) * k - (10 * n + 12) * m - 2)
                self.assertLessEqual(8 * terms, 48 * n * k)
                self.assertLessEqual(cliffords, 60 * n * k)
                total_g += cliffords
                nk_total += n * k
                first_pool = (n + 1) * (k - m)
                work = m + max(n * k, first_pool)
                self.assertLessEqual(work, (n + 1) * k)
                self.assertLessEqual(work, 2 * n * k)
                total_work += work
                levels[level] = levels.get(level, 0) + work
            self.assertLess(nk_total, 3 ** r)
            self.assertLess(total_t, 48 * 3 ** r)
            self.assertLess(total_g, 60 * 3 ** r)
            self.assertLess(total_work, 2 * 3 ** r)
            self.assertLess(max(levels.values()), 2 * 3 ** r)
            # Six calls occupy only five native stages; each stage is the
            # proved private-target eight-layer word, including the first
            # joint A1/Q word verified above.
            self.assertEqual(5 * 8 * len(levels),
                             40 * (r - 1).bit_length())
            self.assertEqual(len(levels), (r - 1).bit_length())


if __name__ == '__main__':
    unittest.main()
