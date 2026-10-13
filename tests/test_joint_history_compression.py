"""Finite evidence for joint-history extraction and a label-preserving limit.

The supplied dense compressor is not a native emitter. These checks test
the full-unitary extraction identities using only its compression columns,
including arbitrary occupied-sector and dirty-register action. The n=5,
D=2 diagnostic has valid banks but is not the D=n-8 endpoint allocation.
"""
import unittest

import numpy as np

from tests.test_depth_balanced_rejection import (
    ATOL, HADAMARD, PAULIS, bank_modes, canonical_depth_fixture,
    cap_squared, routed_cap_gate,
)


def embedding(dimension, labels):
    return np.eye(dimension, dtype=complex)[:, list(labels)]


def random_unitary(dimension, rng):
    matrix = (rng.normal(size=(dimension, dimension))
              + 1j*rng.normal(size=(dimension, dimension)))
    return np.linalg.qr(matrix)[0]


def extraction_fixture():
    f = canonical_depth_fixture()
    size, dimension = f['size'], f['dimension']
    eye = np.eye(dimension, dtype=complex)
    success = embedding(dimension, range(size))
    cap = embedding(dimension, range(2*size, 3*size))
    last = f['stages'][-1]
    support, bank = last['support'], last['bank']
    failure = embedding(dimension, bank)
    history = f['history'] @ f['preweight'] @ success
    history[:2*size, :] = 0
    compressor = f['compression'] @ f['unwind']
    # The last unwind block stores the actual preceding logical inverse,
    # including all nontrivial success-only coarse interlayers.
    prefix = np.eye(size, dtype=complex)
    prefix[np.ix_(support, support)] = f['unwind'][np.ix_(bank, bank)].conj().T
    projector = np.zeros((size, size), dtype=complex)
    projector[support, support] = 1

    # Use the complete routed controlled-X construction, rather than a
    # freely supplied ideal cap-to-bank isometry.
    swap, _ = routed_cap_gate(f['n'], last['depth'], last['axis'], PAULIS[0])
    flag_flip = np.zeros_like(eye)
    flag_flip[np.arange(dimension) ^ (2*size), np.arange(dimension)] = 1
    return dict(f=f, size=size, dimension=dimension, success=success, cap=cap,
                failure=failure, history=history, compressor=compressor,
                prefix=prefix, projector=projector, support=support,
                swap=swap, flag_flip=flag_flip)


def extraction_word(compressor, cap, failure, swap, flag_flip):
    """Actual inverse and literal minus: all input columns are specified."""
    eye = np.eye(len(compressor), dtype=complex)
    r_cap = eye-2*cap @ cap.conj().T
    r_failure = eye-2*failure @ failure.conj().T
    inverse = compressor.conj().T
    recovery = -inverse @ r_cap @ compressor @ r_failure @ inverse
    return flag_flip @ swap @ recovery @ flag_flip, recovery


class JointHistoryCompressionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.fixture = extraction_fixture()
        f = cls.fixture
        dirty_dimension = 2
        cls.dirty_dimension = dirty_dimension
        cls.lift = {name: np.kron(f[name], np.eye(dirty_dimension))
                    for name in ('success', 'cap', 'failure', 'history', 'compressor',
                                 'prefix', 'projector', 'swap', 'flag_flip')}
        # Fix the entire cap output, but freely mix every other mode and
        # dirty value. This preserves HC=J/8 without fixing S or baseline.
        rng = np.random.default_rng(9407)
        other_modes = list(range(2*f['size']))+list(range(3*f['size'], f['dimension']))
        indices = [mode*dirty_dimension+b for mode in other_modes
                   for b in range(dirty_dimension)]
        arbitrary_extension = np.eye(f['dimension']*dirty_dimension, dtype=complex)
        arbitrary_extension[np.ix_(indices, indices)] = random_unitary(len(indices), rng)
        cls.arbitrary_extension = arbitrary_extension
        cls.compressor = arbitrary_extension @ cls.lift['compressor']

    def test_normalized_history_and_last_bank_half_prefix(self):
        f = self.fixture
        normalized = 8*f['history']
        np.testing.assert_allclose(normalized.conj().T @ normalized, np.eye(f['size']),
                                   atol=ATOL, rtol=0)
        np.testing.assert_allclose(f['compressor'] @ f['history'], f['cap']/8,
                                   atol=ATOL, rtol=0)
        np.testing.assert_allclose(f['failure'].conj().T @ normalized,
                                   f['prefix'][list(f['support']), :]/2,
                                   atol=ATOL, rtol=0)
        inactive = np.eye(f['size'])-f['projector']
        np.testing.assert_allclose(normalized @ inactive, f['cap'] @ inactive,
                                   atol=ATOL, rtol=0)
        self.assertGreater(np.linalg.norm(f['prefix'].imag), .1)

    def test_actual_inverse_recovers_full_prefix_with_arbitrary_dirty_extension(self):
        f, lift, h = self.fixture, self.lift, self.compressor
        eye = np.eye(len(h), dtype=complex)
        normalized = 8*lift['history']
        np.testing.assert_allclose(h @ lift['history'], lift['cap']/8,
                                   atol=ATOL, rtol=0)
        np.testing.assert_allclose(h.conj().T @ lift['cap'], normalized,
                                   atol=ATOL, rtol=0)
        self.assertGreater(np.linalg.norm(h-lift['compressor'], 2), 1)
        self.assertGreater(np.linalg.norm(h @ lift['success']-lift['success'], 2), 1)
        dirty_z = np.kron(np.eye(f['dimension']), PAULIS[2])
        self.assertGreater(np.linalg.norm(self.arbitrary_extension @ dirty_z
                                          - dirty_z @ self.arbitrary_extension, 2), 1)
        word, recovery = extraction_word(h, lift['cap'], lift['failure'],
                                         lift['swap'], lift['flag_flip'])
        np.testing.assert_allclose(word.conj().T @ word, eye, atol=ATOL, rtol=0)
        np.testing.assert_allclose(word @ lift['success'],
                                   lift['success'] @ lift['prefix'], atol=ATOL, rtol=0)

        # Distinct input/output reflections give the same accepted cubic,
        # but here the no-entry cap columns must also survive coherently.
        b = lift['failure'].conj().T @ h.conj().T @ lift['cap']
        np.testing.assert_allclose(lift['failure'].conj().T @ recovery @ lift['cap'],
                                   3*b-4*b @ b.conj().T @ b, atol=ATOL, rtol=0)
        inactive = np.eye(len(lift['prefix']))-lift['projector']
        np.testing.assert_allclose(recovery @ lift['cap'] @ inactive,
                                   lift['cap'] @ inactive, atol=ATOL, rtol=0)
        # Omitting the literal global minus reverses every recovered input.
        self.assertGreater(np.linalg.norm(-word @ lift['success']
                                          - lift['success'] @ lift['prefix'], 2), 1.99)

    def test_compression_column_error_alone_bounds_recovery_by_24_epsilon(self):
        f, lift = self.fixture, self.lift
        normalized = 8*lift['history']
        for theta in (1e-4, .03, .3):
            with self.subTest(theta=theta):
                perturb = np.eye(f['dimension'], dtype=complex)
                rotation = [[np.cos(theta), -np.sin(theta)],
                            [np.sin(theta), np.cos(theta)]]
                for mode in range(f['size']):
                    labels = [mode, 2*f['size']+mode]
                    perturb[np.ix_(labels, labels)] = rotation
                h = np.kron(perturb, np.eye(self.dirty_dimension)) @ self.compressor
                epsilon = np.linalg.norm(h @ lift['history']-lift['cap']/8, 2)
                x = h.conj().T @ lift['cap']
                delta = np.linalg.norm(x-normalized, 2)
                self.assertAlmostEqual(delta, 8*epsilon, delta=ATOL)
                projector_error = np.linalg.norm(x @ x.conj().T
                                                  - normalized @ normalized.conj().T, 2)
                self.assertLessEqual(projector_error, delta+ATOL)
                word, _ = extraction_word(h, lift['cap'], lift['failure'],
                                           lift['swap'], lift['flag_flip'])
                error = np.linalg.norm(word @ lift['success']
                                        - lift['success'] @ lift['prefix'], 2)
                self.assertGreater(error, epsilon)
                self.assertLessEqual(error, 24*epsilon+ATOL)
                # In particular, theta=.3 is outside epsilon<=1/32:
                # this reflection comparison needs no small-error premise.
                if theta == .3:
                    self.assertGreater(epsilon, 1/32)

    def test_routed_output_swap_covers_cap_bank_and_inactive_columns(self):
        f = self.fixture
        expected = np.eye(f['dimension'], dtype=complex)
        for mode, failure in zip(f['support'], f['f']['stages'][-1]['bank']):
            cap_mode = 2*f['size']+mode
            expected[cap_mode, cap_mode] = expected[failure, failure] = 0
            expected[cap_mode, failure] = expected[failure, cap_mode] = 1
        np.testing.assert_array_equal(f['swap'], expected)
        np.testing.assert_array_equal(f['flag_flip'] @ f['flag_flip'],
                                      np.eye(f['dimension']))
        np.testing.assert_allclose(f['flag_flip'] @ f['swap'] @ f['failure'],
                                   f['success'][:, list(f['support'])], atol=ATOL, rtol=0)

    def test_label_preserving_mixer_cannot_erase_a_coarse_prefix(self):
        f = self.fixture
        size, dimension = f['size'], f['dimension']
        depth_count, dirty_dimension = 2, 2
        # All six residual coins are identity. A success-only exact HZ
        # coarse word follows depth zero, so every depth-one bank contains
        # this same nontrivial root prefix; fine residuals remain zero.
        root = (0, size//2)
        coarse = np.eye(size, dtype=complex)
        coarse[np.ix_(root, root)] = HADAMARD @ PAULIS[2]
        entry = np.full(size, depth_count, dtype=int)
        groups = [[2*size+mode] for mode in range(size)]
        history = np.zeros((dimension, size), dtype=complex)
        for stage in f['f']['stages']:
            depth, axis, support = stage['depth'], stage['axis'], stage['support']
            amplitude = np.sqrt(float(cap_squared(depth, depth_count)))
            prefix = np.eye(size) if depth == 0 else coarse
            bank = bank_modes(f['f']['n'], depth, axis)
            history[bank, :] = amplitude*prefix[list(support), :]
            for mode, failure in zip(support, bank):
                entry[mode] = min(entry[mode], depth)
                groups[mode].append(failure)
        for mode in range(size):
            history[2*size+mode, mode] = np.sqrt(float(cap_squared(int(entry[mode]),
                                                                 depth_count)))
        np.testing.assert_allclose(history.conj().T @ history, np.eye(size)/64,
                                   atol=ATOL, rtol=0)
        history = np.kron(history, np.eye(dirty_dimension))
        cap = np.kron(f['cap'], np.eye(dirty_dimension))
        rng = np.random.default_rng(61014)
        mixer = np.eye(dimension*dirty_dimension, dtype=complex)
        groups = [[mode*dirty_dimension+b for mode in group
                   for b in range(dirty_dimension)] for group in groups]
        for group in groups:
            mixer[np.ix_(group, group)] = random_unitary(len(group), rng)
        mixed = mixer @ history
        # Compare every logical-label probability, including dirty values,
        # for one root input. Positive wrong-label weight cannot cancel
        # across orthogonal banks or through dirty entanglement.
        column = root[0]*dirty_dimension
        probabilities = []
        for group in groups:
            before = np.linalg.norm(history[group, column])**2
            after = np.linalg.norm(mixed[group, column])**2
            self.assertAlmostEqual(before, after, delta=ATOL)
            probabilities.append(before)
        wrong_label_weight = sum(value for mode, value in enumerate(probabilities)
                                 if mode != root[0])
        self.assertAlmostEqual(wrong_label_weight, 3/512, delta=ATOL)
        column_error = np.linalg.norm(mixed[:, column]-cap[:, column]/8)
        self.assertGreaterEqual(column_error+ATOL, np.sqrt(3/512))
        self.assertGreaterEqual(np.linalg.norm(mixed-cap/8, 2)+ATOL, np.sqrt(3/512))


if __name__ == '__main__':
    unittest.main()
