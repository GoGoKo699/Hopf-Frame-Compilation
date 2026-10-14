"""Complete-sector regressions; no fine compiler or endpoint simulation."""
from fractions import Fraction
import unittest
import numpy as np
from compiler_robust_hopf.conjugate_recovery import (
    conjugate_pair, recover, recovery_blocks, two_sector_source,
    lift_native_word, recovery_counterexample_words, finite_phase_word)
from verification.fault_tolerant.exact_arithmetic import Q2, CQ2


def exact_native(word, qubits):
    """Small exact matrix executor using the existing retained coefficient ring."""
    size = 1 << qubits
    zero, one, imag = CQ2(), CQ2(1), CQ2(0, 1)
    halfroot = Q2(0, Fraction(1, 2))
    phases = {'s': imag, 'sdg': -imag, 't': CQ2(halfroot, halfroot),
              'tdg': CQ2(halfroot, -halfroot)}
    rows = [[one if i == j else zero for j in range(size)] for i in range(size)]
    for op, *bits in word:
        if op == 'cx':
            c, t = bits
            rows = [rows[i ^ ((1 << t) if (i >> c) & 1 else 0)] for i in range(size)]
        elif op == 'h':
            bit = 1 << bits[0]
            for i in range(size):
                if i & bit:
                    continue
                a, b = rows[i], rows[i | bit]
                rows[i] = [(x+y)*halfroot for x, y in zip(a, b)]
                rows[i | bit] = [(x-y)*halfroot for x, y in zip(a, b)]
        elif op == 'x':
            rows = [rows[i ^ (1 << bits[0])] for i in range(size)]
        else:
            phase = phases[op]
            rows = [[phase*x for x in row] if (i >> bits[0]) & 1 else row
                    for i, row in enumerate(rows)]
    return rows


def numerical(rows):
    return np.array([[x.display() for x in row] for row in rows])


def unitary(rng, size, real=False):
    a = rng.normal(size=(size, size))
    if not real:
        a = a+1j*rng.normal(size=(size, size))
    return np.linalg.qr(a)[0]


class RecoveryTests(unittest.TestCase):
    def test_four_block_polynomial(self):
        rng = np.random.default_rng(10301)
        for size in (2, 4, 8, 16):
            for real in (False, True):
                u = unitary(rng, size, real)
                np.testing.assert_allclose(recover(u), recovery_blocks(u), atol=2e-12, rtol=0)

    def test_sectors_are_independently_specifiable(self):
        rng = np.random.default_rng(10302)
        for size in (1, 2, 4, 8):
            e, f = unitary(rng, size), unitary(rng, size)
            u = two_sector_source(e, f)
            expected = np.block([[e, np.zeros_like(e)], [np.zeros_like(e), f]])
            np.testing.assert_allclose(u.conj().T@u, np.eye(2*size), atol=2e-12, rtol=0)
            np.testing.assert_allclose(recover(u), expected, atol=2e-12, rtol=0)

    def test_balanced_two_sector_bound(self):
        rng = np.random.default_rng(10303)
        for size in (2, 4, 8):
            target = np.kron(unitary(rng, size//2, True), np.eye(2))
            base = two_sector_source(target, target)
            for scale in (1e-1, 1e-3, 1e-5):
                h = rng.normal(size=base.shape)+1j*rng.normal(size=base.shape)
                h = (h+h.conj().T)/2
                values, vectors = np.linalg.eigh(h)
                perturb = (vectors*np.exp(1j*scale*values))@vectors.conj().T
                u = perturb@base
                e0 = np.linalg.norm(u[:size, :size]-target/2, 2)
                e1 = np.linalg.norm(u[size:, size:]-target/2, 2)
                error = np.linalg.norm(recover(u)-np.kron(np.eye(2), target), 2)
                self.assertLessEqual(error, np.sqrt(136*(e0*e0+e1*e1))+1e-12)

    def test_conjugate_pair_is_multiplicative(self):
        rng = np.random.default_rng(10304)
        a, b = unitary(rng, 4), unitary(rng, 4)
        np.testing.assert_allclose(conjugate_pair(a@b), conjugate_pair(a)@conjugate_pair(b), atol=1e-12)

    def test_real_target_retains_only_mode_initialization(self):
        rng = np.random.default_rng(10305)
        target, reject = unitary(rng, 4, True), unitary(rng, 4)
        full = np.block([[target, np.zeros((4, 4))], [np.zeros((4, 4)), reject]])
        j = np.vstack([np.eye(4), np.zeros((4, 4))])
        lifted_j = np.kron(np.eye(2), j)
        np.testing.assert_allclose(conjugate_pair(full)@lifted_j,
                                   lifted_j@np.kron(np.eye(2), target), atol=1e-12)
        self.assertGreater(np.linalg.norm(conjugate_pair(full)-np.kron(np.eye(2), full)), .1)

    def test_invalid_recovery_inputs(self):
        with self.assertRaises(ValueError):
            recovery_blocks(np.eye(3))
        with self.assertRaises(ValueError):
            recover(np.eye(4), True)
        with self.assertRaises(ValueError):
            recover(np.eye(4), 4)
        with self.assertRaises(ValueError):
            two_sector_source(np.eye(2), np.zeros((2, 2)))


class NativeSignalTests(unittest.TestCase):
    def test_exact_native_lifts(self):
        for gate in ('h', 's', 'sdg', 't', 'tdg', 'x'):
            raw = numerical(exact_native([(gate, 0)], 1))
            lifted = numerical(exact_native(lift_native_word([(gate, 0)], 1), 2))
            np.testing.assert_allclose(lifted, conjugate_pair(raw), atol=1e-14, rtol=0)
        for gate in ('t', 'tdg'):
            self.assertEqual(sum(g[0] in ('t', 'tdg') for g in lift_native_word([(gate, 0)], 1)), 2)

    def test_borrowed_signal_recovery_reflects_only_the_original_mode_anchor(self):
        # Four wires: arbitrary signal r=3, initialized mode suffix f=2,
        # logical x=1, and arbitrary dirty d=0. The source is supplied as a
        # diagnostic unitary; only the coarse interlayer is a native word.
        # This checks the initialized-real-target contract, not an emitter.
        logical = np.array([[3., -4.], [4., 3.]]) / 5
        target = np.kron(logical, np.eye(2))
        reject = numerical(exact_native([('h', 0), ('s', 1), ('cx', 0, 1)], 2))
        source = two_sector_source(target, reject)
        coarse_word = [('h', 0), ('t', 1), ('cx', 1, 2), ('s', 0)]
        coarse = numerical(exact_native(coarse_word, 3))
        lifted_coarse = numerical(exact_native(lift_native_word(coarse_word, 3), 4))
        fine = conjugate_pair(coarse @ source)
        lifted_source = lifted_coarse.conj().T @ fine
        np.testing.assert_allclose(lifted_source, conjugate_pair(source), atol=2e-12, rtol=0)

        mode_anchor = np.eye(8)[:, :4]
        anchor = np.kron(np.eye(2), mode_anchor)
        reflection = 2 * anchor @ anchor.conj().T - np.eye(16)
        target_with_signal = np.kron(np.eye(2), target)
        wanted = anchor @ target_with_signal

        def recovery(word, reflected=reflection):
            return -word @ reflected @ word.conj().T @ reflected @ word

        recovered = recovery(lifted_source)
        np.testing.assert_allclose(recovered @ anchor, wanted, atol=3e-12, rtol=0)
        np.testing.assert_allclose(recovered, conjugate_pair(recover(source)), atol=3e-12, rtol=0)
        self.assertGreater(np.linalg.norm(recovered.imag, 2), .5)
        self.assertGreater(np.linalg.norm(recovered - np.kron(np.eye(2), recover(source)), 2), 1.)

        # Including r in the anchor would silently initialize a borrowed bit.
        wrong_reflection = 2 * anchor[:, :4] @ anchor[:, :4].conj().T - np.eye(16)
        self.assertGreater(np.linalg.norm(recovery(lifted_source, wrong_reflection) @ anchor
                                         - wanted, 2), 1.)
        unlifted_coarse = np.kron(np.eye(2), coarse).conj().T @ fine
        self.assertGreater(np.linalg.norm(recovery(unlifted_coarse) @ anchor - wanted, 2), .1)

        # Model a complete fine-slot error that mixes both r sectors and the
        # initialized suffix while acting on dirty d. Such leakage cannot be
        # removed before composing the actual three-call recovery word.
        x = np.array([[0., 1.], [1., 0.]])
        perturbation = np.kron(x, np.kron(x, np.kron(np.diag([1., -1.]), x)))
        coefficients = np.arange(1, 9) + 1j * np.arange(8, 0, -1)
        coefficients /= np.linalg.norm(coefficients)
        # Eight reference labels give full Schmidt support on all arbitrary
        # signal/logical/dirty inputs, with only the original suffix fixed.
        referenced_input = anchor * coefficients
        for delta in (1 / 256, 1 / 1024):
            perturbed_fine = (np.cos(delta) * np.eye(16)
                              + 1j * np.sin(delta) * perturbation) @ fine
            slot_error = np.linalg.norm(perturbed_fine - fine, 2)
            actual = recovery(lifted_coarse.conj().T @ perturbed_fine)
            np.testing.assert_allclose(actual.conj().T @ actual, np.eye(16), atol=4e-12, rtol=0)
            self.assertLessEqual(np.linalg.norm(actual - recovered, 2), 3 * slot_error + 4e-12)
            error = np.linalg.norm(actual @ anchor - wanted, 2)
            self.assertLessEqual(error, 3 * slot_error + 4e-12)
            self.assertGreater(np.linalg.norm((np.eye(16) - anchor @ anchor.conj().T)
                                             @ actual @ anchor, 2), delta / 10)
            self.assertLessEqual(np.linalg.norm(actual @ referenced_input - wanted * coefficients),
                                 error + 4e-12)
            np.testing.assert_allclose(actual.conj().T @ actual @ anchor, anchor, atol=4e-12, rtol=0)
            self.assertLessEqual(np.linalg.norm(actual.conj().T @ anchor
                                               - anchor @ target_with_signal.conj().T, 2),
                                 error + 4e-12)

    def test_exact_recovery_counterexample(self):
        source, recovery = recovery_counterexample_words()
        self.assertEqual((len(source), len(recovery)), (51, 169))
        self.assertEqual(sum(g[0] in ('t', 'tdg') for g in source), 4)
        self.assertEqual(sum(g[0] in ('t', 'tdg') for g in recovery), 12)
        self.assertLessEqual({g[0] for g in recovery}, {'h', 's', 'sdg', 't', 'tdg', 'cx'})
        u, r = exact_native(source, 4), exact_native(recovery, 4)
        zero, one, half = CQ2(), CQ2(1), CQ2(Fraction(1, 2))
        for i in range(16):
            for j in range(16):
                expected = (-one if i & 8 and i & 4 else one) if i == j else zero
                self.assertEqual(r[i][j], expected)
                self.assertEqual(u[i][j].imag, Q2())
                if i < 8 and j < 8:
                    self.assertEqual(u[i][j], half if i == j else zero)
        for i in range(16):
            for j in range(16):
                value = sum((u[k][i].conj()*u[k][j] for k in range(16)), zero)
                self.assertEqual(value, one if i == j else zero)

    def test_finite_phase_full_matrix(self):
        for q, phi in ((5, -.7), (5, .31), (6, 1.2)):
            block, v, b, z = finite_phase_word(phi, q)
            eye = np.eye(len(b))
            x = abs(z)**2
            f, g = 5-20*x+16*x*x, 1-12*x+16*x*x
            np.testing.assert_allclose(block, np.block([[z*eye, b], [-b, z.conjugate()*eye]]), atol=1e-12)
            np.testing.assert_allclose(b.conj().T, b, atol=1e-12)
            np.testing.assert_allclose(b@b, (1-x)*eye, atol=1e-12)
            np.testing.assert_allclose(v, np.block([[f*z*eye, g*b], [-g*b, f*z.conjugate()*eye]]), atol=2e-12)
            wanted = np.kron(np.diag([np.exp(1j*phi), np.exp(-1j*phi)]), eye)
            error = v-wanted
            alpha = f*z-np.exp(1j*phi)
            np.testing.assert_allclose(error.conj().T@error,
                                       (abs(alpha)**2+(1-x)*g*g)*np.eye(len(v)), atol=2e-12)
            self.assertLess(np.linalg.norm(error, 2), 30*2.**(-q))
            self.assertAlmostEqual(np.linalg.norm(error, 2), np.linalg.norm(error[:, :len(b)], 2), places=11)
            self.assertGreater(np.linalg.norm(v-np.exp(1j*phi)*np.eye(len(v)), 2), .1)

    def test_invalid_signal_inputs(self):
        for bad in (True, -1):
            with self.assertRaises(ValueError):
                lift_native_word([('h', 0)], bad)
        with self.assertRaises(ValueError):
            lift_native_word([('cx', 0, 1)], 1)
        with self.assertRaises(ValueError):
            finite_phase_word(.1, 100)
        with self.assertRaises(ValueError):
            finite_phase_word(float('nan'), 5)


if __name__ == '__main__':
    unittest.main()
