"""Bounded all-input checks for the dirty catalogue terminal interpreter.

The native H/S/T words below retain literal phase and actual inverses.
These finite operator checks do not constitute a scalable native emitter
or certify the asymptotic resource allocation.
"""
from __future__ import annotations

import unittest

import numpy as np

from compiler_robust_hopf.frames import real_frame_matrix


ATOL = 5e-12
I2 = np.eye(2, dtype=complex)
H = np.array([[1, 1], [1, -1]], complex) / np.sqrt(2)
S = np.diag([1, 1j]).astype(complex)
T = np.diag([1, np.exp(1j*np.pi/4)]).astype(complex)


def bit(j, q, n):
    return (j >> (n - 1 - q)) & 1


def swap_bits(j, a, b, n):
    return j ^ (((1 << (n-1-a)) | (1 << (n-1-b)))
                if bit(j,a,n) != bit(j,b,n) else 0)


def permutation(fn, n):
    out = np.zeros((1 << n, 1 << n), complex)
    for j in range(1 << n):
        out[fn(j), j] = 1
    return out


def local(u, q, n):
    out = np.array([[1]], complex)
    for t in range(n):
        out = np.kron(out, u if t == q else I2)
    return out


def chronological(word, n):
    out = np.eye(1 << n, dtype=complex)
    for gate in word:
        out = gate @ out
    return out


def inverse_word(word):
    return [gate.conj().T for gate in reversed(word)]


class DirtyCatalogueTerminalTests(unittest.TestCase):
    @staticmethod
    def _native_fixture(shared_bank=False):
        n = 7
        p, r, data = 0, 1, 2
        bank_a, bank_b = (3,4), (5,6)
        if shared_bank:
            bank_b = bank_a
        query = permutation(lambda j: j ^ ((1 << (n-1-r)) if bit(j,p,n) else 0),n)

        def shifted_word(bank):
            shift_r = permutation(lambda j: swap_bits(j,*bank,n) if bit(j,r,n) else j,n)
            return [query, shift_r, query.conj().T, shift_r.conj().T]

        def selected_swap_word(bank):
            shifts = shifted_word(bank)
            fixed = permutation(lambda j: swap_bits(j,data,bank[0],n),n)
            return shifts + [fixed] + inverse_word(shifts)

        # Actual native H/S/T/T-dagger words, not abstract angle rotations.
        # Each list is chronological, and every inverse below is the actual
        # reversed-adjoint word. Scalar phases are retained throughout.
        a_words = [[H,T,H], [S,H,T.conj().T,H]]
        b_words = [[T,H,T], [H,S.conj().T,T,H,T.conj().T]]
        a = [chronological(w,1) for w in a_words]
        b = [chronological(w,1) for w in b_words]
        a_all = [local(g,bank_a[v],n) for v,w in enumerate(a_words) for g in w]
        b_all = [local(g,bank_b[v],n) for v,w in enumerate(b_words) for g in w]
        sa, sb = selected_swap_word(bank_a), selected_swap_word(bank_b)
        # K_A = S_A A_all S_A A_all^dagger, and similarly for B.
        ka = inverse_word(a_all) + sa + a_all + sa
        kb = inverse_word(b_all) + sb + b_all + sb
        # Full operator K_A K_B K_A^dagger K_B^dagger.
        actual = chronological(inverse_word(kb)+inverse_word(ka)+kb+ka,n)
        expected = np.zeros_like(actual)
        commutators = []
        for v in (0,1):
            u = a[v] @ b[v] @ a[v].conj().T @ b[v].conj().T
            commutators.append(u)
            projector = local(np.diag([1-v,v]).astype(complex),p,n)
            expected += projector @ local(u,data,n)
        return actual, expected, a, b, commutators, shifted_word, selected_swap_word

    @classmethod
    def setUpClass(cls):
        cls.fixture = cls._native_fixture()

    def test_native_catalogue_full_operator(self):
        n, p, data = 7, 0, 2
        bank_a, bank_b = (3,4), (5,6)
        actual, expected, a, b, commutators, shifted_word, selected_swap_word = self.fixture
        for v in (0,1):
            self.assertGreater(np.linalg.norm(a[v]@b[v]-b[v]@a[v]), .1)
        self.assertGreater(np.linalg.norm(commutators[0]-commutators[1]),.1)
        self.assertLess(np.linalg.norm(actual-expected,ord=2), ATOL)
        for bank in (bank_a,bank_b):
            desired_shift = permutation(lambda j: swap_bits(j,*bank,n) if bit(j,p,n) else j,n)
            self.assertTrue(np.array_equal(chronological(shifted_word(bank),n),desired_shift))
            desired_swap = permutation(lambda j: swap_bits(j,data,bank[bit(j,p,n)],n),n)
            self.assertTrue(np.array_equal(chronological(selected_swap_word(bank),n),desired_swap))

    def test_bell_reference_and_coherent_prefix(self):
        actual, expected = self.fixture[:2]
        # R is maximally entangled with an explicit reference. The other
        # six qubits have a generic coherent state, including both P values.
        rng = np.random.default_rng(9001)
        rest = rng.normal(size=64) + 1j*rng.normal(size=64)
        rest /= np.linalg.norm(rest)
        psi0 = np.zeros(128, complex)
        for j, amplitude in enumerate(rest):
            before, after = j >> 5, j & 31
            psi0[(before << 6) | after] = amplitude
        flip_r = permutation(lambda j: j ^ (1 << 5),7)
        psi1 = flip_r @ psi0
        state = (np.kron(psi0,[1,0]) + np.kron(psi1,[0,1])) / np.sqrt(2)
        output = np.kron(actual,I2) @ state
        target = np.kron(expected,I2) @ state
        np.testing.assert_allclose(output,target,atol=ATOL)
        # Trace out P, data, and the banks; the dirty/reference Bell pair
        # itself returns, rather than merely its individual marginals.
        partition = output.reshape((2,)*8).transpose(1,7,0,2,3,4,5,6).reshape(4,-1)
        bell = np.array([1,0,0,1],complex) / np.sqrt(2)
        np.testing.assert_allclose(partition @ partition.conj().T,
                                   np.outer(bell,bell.conj()),atol=ATOL)

    def test_reusing_one_array_leaves_nontrivial_dirty_commutator(self):
        actual, target, a, b = self._native_fixture(shared_bank=True)[:4]
        self.assertGreater(np.linalg.norm(actual-target,ord=2), .1)
        # On the shared selected bank, the uncancelled residue is
        # [A_v^dagger,B_v^dagger]. Disjoint arrays remove precisely this
        # obstruction; mere selection does not return arbitrary dirty work.
        expected = np.zeros_like(actual)
        for v in (0,1):
            logical = a[v] @ b[v] @ a[v].conj().T @ b[v].conj().T
            residue = a[v].conj().T @ b[v].conj().T @ a[v] @ b[v]
            self.assertGreater(np.linalg.norm(residue-I2), .1)
            projector = local(np.diag([1-v,v]).astype(complex),0,7)
            expected += projector @ local(logical,2,7) @ local(residue,3+v,7)
        np.testing.assert_allclose(actual,expected,atol=ATOL)

    def test_computable_commutator_factors_for_real_four_column_frame(self):
        # This bounded spectral calculation checks the algebraic factor
        # choice. It is not the native synthesis or a complexity certificate.
        target = real_frame_matrix([.21,-.37,.53]).astype(complex)
        self.assertAlmostEqual(np.linalg.det(target).real,1,delta=ATOL)
        eigenvalues, basis = np.linalg.eig(target)
        np.testing.assert_allclose(basis.conj().T @ basis,np.eye(4),atol=ATOL)
        diagonal = np.ones(4,complex)
        for j in range(1,4):
            diagonal[j] = diagonal[j-1] * eigenvalues[j]
        cycle = np.roll(np.eye(4),1,axis=0)
        a = basis @ np.diag(diagonal) @ basis.conj().T
        b = basis @ cycle @ basis.conj().T
        np.testing.assert_allclose(a.conj().T @ a,np.eye(4),atol=ATOL)
        np.testing.assert_allclose(b.conj().T @ b,np.eye(4),atol=ATOL)
        np.testing.assert_allclose(a @ b @ a.conj().T @ b.conj().T,
                                   target,atol=ATOL)

    def test_two_bit_nonlinear_program_translation_all_basis_states(self):
        n = 9
        pbits, rbits, data, bank = (0,1), (2,3), 4, (5,6,7,8)

        def read(j,bits):
            return sum(bit(j,q,n) << (len(bits)-1-k) for k,q in enumerate(bits))

        def f(j):
            p0,p1 = (bit(j,q,n) for q in pbits)
            return ((p0 ^ p1) << 1) | (p0 & p1)

        def query(j):
            v = f(j)
            for k,q in enumerate(rbits):
                if (v >> (1-k)) & 1:
                    j ^= 1 << (n-1-q)
            return j

        def shift_value(j,v):
            old = [bit(j,q,n) for q in bank]
            for t,q in enumerate(bank):
                if bit(j,q,n) != old[t ^ v]:
                    j ^= 1 << (n-1-q)
            return j

        def shift_r(j):
            return shift_value(j,read(j,rbits))

        def translate(j):
            for gate in (query,shift_r,query,shift_r):
                j = gate(j)
            return j

        def translate_inverse(j):
            for gate in (shift_r,query,shift_r,query):
                j = gate(j)
            return j

        for j in range(1 << n):
            self.assertEqual(translate(j),shift_value(j,f(j)))
            self.assertEqual(translate_inverse(translate(j)),j)
            routed = translate(j)
            routed = swap_bits(routed,data,bank[0],n)
            routed = translate_inverse(routed)
            self.assertEqual(routed,swap_bits(j,data,bank[f(j)],n))

    def test_shared_control_ccz_phase_polynomial(self):
        # Two CCZ gates sharing a, on otherwise disjoint pairs. This checks
        # the aggregated common-control phase used for parallel Fredkins.
        for a in (0,1):
            for b0 in (0,1):
                for c0 in (0,1):
                    for b1 in (0,1):
                        for c1 in (0,1):
                            exponent = 2*a
                            for b,c in ((b0,c0),(b1,c1)):
                                exponent += b+c-(a^b)-(a^c)-(b^c)+(a^b^c)
                            self.assertEqual(exponent % 8, (4*a*(b0*c0+b1*c1)) % 8)


if __name__ == "__main__":
    unittest.main()
