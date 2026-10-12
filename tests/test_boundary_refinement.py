"""Scoped rank/boundary regressions, not native endpoint certificates."""
from fractions import Fraction
import itertools
import unittest
import numpy as np
from compiler_robust_hopf.boundary_refinement import (
    reflection, real_state_extension, commutator, spectral_reflection_lcu,
    rank_mode_bins, forest_matching, canonical_frame, canonical_star_word,
    leaf_deflation, recursive_leaf_frame)
from compiler_robust_hopf.conjugate_recovery import recover


def unitary(rng, n, real=False):
    a = rng.normal(size=(n, n))
    if not real:
        a = a+1j*rng.normal(size=(n, n))
    return np.linalg.qr(a)[0]


def embed(block, labels, n):
    out = np.eye(n, dtype=complex)
    out[np.ix_(labels, labels)] = block
    return out


def rotation(t):
    return np.array([[np.cos(t), -np.sin(t)], [np.sin(t), np.cos(t)]])


class RankWeightedTests(unittest.TestCase):
    def test_coarse_rank_doubles_not_preserves_residual_rank(self):
        target = np.diag([np.exp(.1j), np.exp(-.1j), 1, 1])
        coarse = np.diag([1, 1, np.exp(.03j), np.exp(-.03j)])
        self.assertEqual(np.linalg.matrix_rank(target-np.eye(4), tol=1e-12), 2)
        self.assertEqual(np.linalg.matrix_rank(coarse-np.eye(4), tol=1e-12), 2)
        self.assertEqual(np.linalg.matrix_rank(coarse.conj().T@target-np.eye(4), tol=1e-12), 4)

    def test_reflection_lcu_complete_recovery(self):
        rng = np.random.default_rng(10311)
        for m in (2, 4, 8):
            for rank in (1, m):
                frame = unitary(rng, m)[:, :rank]
                phases = rng.uniform(-.02, .02, rank)
                word, target, lam = spectral_reflection_lcu(frame, phases)
                self.assertLess(lam, 1.25)
                np.testing.assert_allclose(word.conj().T@word, np.eye(len(word)), atol=2e-12)
                np.testing.assert_allclose(word[:m, :m], target/2, atol=2e-12)
                output = recover(word, m)[:, :m]
                expected = np.vstack([target, np.zeros((len(word)-m, m))])
                np.testing.assert_allclose(output, expected, atol=3e-12)

    def test_zero_eigenphases(self):
        word, target, _ = spectral_reflection_lcu(np.eye(3), np.zeros(3))
        np.testing.assert_allclose(recover(word, 3)[:3, :3], np.eye(3), atol=1e-12)
        np.testing.assert_allclose(target, np.eye(3), atol=1e-12)

    def test_mode_bin_bounds(self):
        for n in range(10, 18):
            ambient = 1 << n
            sizes = [ambient//128, ambient//256, 2]
            ranks = [min(4, sizes[0]), min(6, sizes[1]), 2]
            bins = rank_mode_bins(ambient, sizes, ranks)
            seen = [0]*len(sizes)
            volume = 0
            for bank in bins:
                self.assertEqual(len({j for j, _, _ in bank}), len(bank))
                weight = sum(m*(q+3) for _, m, q in bank)
                self.assertLessEqual(weight, ambient//4)
                volume += weight
                for j, m, q in bank:
                    self.assertLessEqual(m*(q+3), ambient//8)
                    seen[j] += q
            self.assertEqual(seen, [min(2*r, m) for m, r in zip(sizes, ranks)])
            self.assertLessEqual(volume, 4*sum(m*r for m, r in zip(sizes, ranks)))

    def test_invalid_lcu_and_bins(self):
        with self.assertRaises(ValueError):
            spectral_reflection_lcu(np.ones((2, 2)), [0, 0])
        with self.assertRaises(ValueError):
            spectral_reflection_lcu(np.eye(4), [np.pi]*4)
        with self.assertRaises(ValueError):
            rank_mode_bins(1000, [2], [2])
        with self.assertRaises(ValueError):
            rank_mode_bins(256, [3], [2])


class BoundaryFusionTests(unittest.TestCase):
    def test_boundary_column_cancellation(self):
        rng = np.random.default_rng(10312)
        n, anchors = 12, [0, 4, 8]
        u = embed(unitary(rng, 3, True), anchors, n)
        v, s = np.eye(n), np.eye(n)
        for a in anchors:
            ids = list(range(a, a+4))
            local = unitary(rng, 4, True)
            v[np.ix_(ids, ids)] = local
            s[np.ix_(ids, ids)] = real_state_extension(local[:, 0])
        np.testing.assert_allclose(commutator(u, v), commutator(u, s), atol=2e-12)
        np.testing.assert_allclose((s.T@v)@u, u@(s.T@v), atol=2e-12)

    def test_transported_product(self):
        rng = np.random.default_rng(10313)
        u, a, b = [unitary(rng, 8, True) for _ in range(3)]
        np.testing.assert_allclose(commutator(u, a@b),
                                   commutator(u, a)@a@commutator(u, b)@a.T, atol=2e-12)
        self.assertGreater(np.linalg.norm(commutator(u, a@b)-commutator(u, a)@commutator(u, b)), .1)

    def test_collective_basis_and_growing_channels(self):
        rng = np.random.default_rng(10314)
        child_size, children, channels = 6, 3, 2
        n = child_size*children
        anchors = [j*child_size for j in range(children)]
        basis = np.eye(n)
        active = []
        for a in anchors:
            basis[a+1:a+child_size, a+1:a+child_size] = unitary(rng, child_size-1, True)
            active.extend(range(a, a+channels+1))
        product = np.eye(n)
        core = np.eye(len(active))
        for _ in range(11):
            u = embed(unitary(rng, children, True), anchors, n)
            small_v = np.eye(n)
            for a in anchors:
                vector = rng.normal(size=channels+1)
                vector /= np.linalg.norm(vector)
                small_v[a:a+channels+1, a:a+channels+1] = real_state_extension(vector)
            v = basis@small_v@basis.T
            k = commutator(u, v)
            reduced = basis.T@k@basis
            core = reduced[np.ix_(active, active)]@core
            product = k@product
        np.testing.assert_allclose(product, basis@embed(core, active, n)@basis.T, atol=3e-12)

    def test_rank_growth_and_multiport_formula(self):
        u = embed(rotation(.4), [0, 1], 5)
        v = embed(rotation(.6), [1, 2], 5)
        w = embed(rotation(-.3), [2, 3], 5)
        z = embed(rotation(.7), [3, 4], 5)
        first, second = commutator(u, v), commutator(w, z)
        self.assertEqual(np.linalg.matrix_rank(first-np.eye(5), tol=1e-10), 2)
        self.assertEqual(np.linalg.matrix_rank(second-np.eye(5), tol=1e-10), 2)
        self.assertEqual(np.linalg.matrix_rank(second@first-np.eye(5), tol=1e-10), 4)
        a = np.kron(np.eye(2), rotation(.3))
        c, s = np.cos(.7), np.sin(.7)
        upper = embed(a, range(4), 8)
        lower = np.block([[c*np.eye(4), -s*np.eye(4)], [s*np.eye(4), c*np.eye(4)]])
        expected = np.block([[c*c*np.eye(4)+s*s*a, c*s*(np.eye(4)-a)],
                             [c*s*(a.T-np.eye(4)), c*c*np.eye(4)+s*s*a.T]])
        np.testing.assert_allclose(commutator(upper, lower), expected, atol=1e-12)
        self.assertEqual(np.linalg.matrix_rank(expected-np.eye(8), tol=1e-10), 8)

    def test_stable_state_extension_corners(self):
        for i in range(5):
            for sign in (-1, 1):
                v = sign*np.eye(5)[:, i]
                s = real_state_extension(v)
                np.testing.assert_allclose(s[:, 0], v, atol=1e-14)
                np.testing.assert_allclose(s.T@s, np.eye(5), atol=1e-14)
        with self.assertRaises(ValueError):
            real_state_extension([0., 0.])


class CanonicalForestTests(unittest.TestCase):
    def test_matching_against_exhaustive_subsets(self):
        rng = np.random.default_rng(10315)
        for n in range(2, 10):
            for _ in range(5):
                parents = [-1]+[int(rng.integers(-1, j)) for j in range(1, n)]
                edges = [(p, j) for j, p in enumerate(parents) if p >= 0]
                maximum = 0
                for mask in range(1 << len(edges)):
                    selected = [edges[k] for k in range(len(edges)) if mask >> k & 1]
                    flat = [v for e in selected for v in e]
                    if len(set(flat)) == len(flat):
                        maximum = max(maximum, len(selected))
                self.assertEqual(len(forest_matching(parents)), maximum)

    def test_canonical_rank_and_star_formula(self):
        rng = np.random.default_rng(10316)
        for n in (4, 8, 16, 32):
            parents = [-1]+[j-(j & -j) if rng.random() < .8 else -1 for j in range(1, n)]
            l = np.zeros((n, n))
            for j, p in enumerate(parents):
                if p >= 0:
                    l[j, p] = rng.choice([-1, 1])*rng.uniform(.3, 1.8)
            w = canonical_frame(l)
            np.testing.assert_allclose(w, canonical_star_word(l), atol=3e-12)
            self.assertEqual(np.linalg.matrix_rank(w-np.eye(n), tol=1e-9), 2*len(forest_matching(parents)))

    def test_exact_leaf_schur_zero(self):
        f = Fraction
        inverse = [[f(0), f(-3, 2)], [f(3, 2), f(0)]]
        coupling = [[f(0), f(0)], [f(2, 5), f(-3, 7)]]
        for i, j in itertools.product(range(2), repeat=2):
            correction = sum(coupling[k][i]*inverse[k][l]*coupling[l][j]
                             for k, l in itertools.product(range(2), repeat=2))
            self.assertEqual(correction, 0)

    def test_deflation_with_signed_and_declared_zero_edges(self):
        rng = np.random.default_rng(10317)
        for m in (2, 4, 8, 16):
            mask = np.zeros((m, m), bool)
            for j in range(1, m):
                mask[j, j-(j & -j)] = True
            for zeros in (False, True):
                a = np.where(mask, rng.uniform(-1, 1, (m, m)), 0)
                d = rng.uniform(-1.3, 1.3, m)
                if zeros:
                    d[::2] = 0
                    a[::3] = 0
                l = np.block([[a, np.zeros_like(a)], [np.diag(d), np.zeros_like(a)]])
                out = leaf_deflation(a, d, mask)
                core = canonical_frame(out.compact_shift)
                full = np.eye(2*m)+out.embedding@(core-np.eye(len(core)))@out.embedding.T
                g, r = out.matching, out.basis
                np.testing.assert_allclose(g@canonical_frame(l)@g, canonical_frame(out.intermediate_shift), atol=3e-12)
                np.testing.assert_allclose(g.T@r@full@r.T@g.T, canonical_frame(l), atol=3e-12)
                np.testing.assert_allclose(r.T@out.intermediate_shift@r,
                                           out.embedding@out.compact_shift@out.embedding.T, atol=3e-12)
                expected = sum(np.sum(a[:, v]**2)/(1+d[v]**2) for v in out.parents)
                self.assertAlmostEqual(np.sum(out.compact_shift**2), expected, places=11)

    def test_full_closed_recursion(self):
        rng = np.random.default_rng(10318)
        for n in (2, 4, 8, 16, 32, 64):
            l = np.zeros((n, n))
            for j in range(1, n):
                l[j, j-(j & -j)] = rng.uniform(-.7, .7)
            np.testing.assert_allclose(recursive_leaf_frame(l), canonical_frame(l), atol=4e-12)

    def test_wrong_outer_orientation_is_detected(self):
        a = np.array([[0., 0.], [.5, 0.]])
        d = np.array([.3, -.7])
        out = leaf_deflation(a, d)
        l = np.block([[a, np.zeros_like(a)], [np.diag(d), np.zeros_like(a)]])
        self.assertGreater(np.linalg.norm(out.matching@canonical_frame(l)@out.matching.T-
                                           canonical_frame(out.intermediate_shift)), .1)

    def test_invalid_forest_contracts(self):
        with self.assertRaises(ValueError):
            forest_matching([-1, 2, 0])
        with self.assertRaises(ValueError):
            canonical_frame(np.array([[0., 1.], [0., 0.]]))
        with self.assertRaises(ValueError):
            leaf_deflation(np.array([[0., 0., 0.], [1., 0., 0.], [1., 1., 0.]]), np.ones(3))
        with self.assertRaises(ValueError):
            recursive_leaf_frame(np.zeros((3, 3)))


if __name__ == '__main__':
    unittest.main()
