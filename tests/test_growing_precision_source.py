"""Bounded full-unitary diagnostics for the growing-source architecture.

The explicit scalar dilations are finite implementations of the unitary
contracts, not Clifford+T synthesis or a native growing-source emitter.
The arena below has 128 coordinates, including arbitrary dirty inputs.
"""
from fractions import Fraction
import unittest

import numpy as np

from compiler_robust_hopf.growing_precision_source import (
    growing_source_scale,
    refinement_parameters,
)
from tests.test_heterogeneous_precision_source import (
    _embed_modes,
    _mixed_chart_fixture,
)
from tests.test_operator_source_compiler import (
    _adjoint,
    _embed_addressed_blocks,
    _sign_word,
    _word_matrix,
)
from tests.test_packed_diagonal import _diagonal_word, _reflection_word


ATOL = 8e-13
X = np.array([[0., 1.], [1., 0.]])
H = np.array([[1., 1.], [1., -1.]]) / np.sqrt(2)
POLARS = (np.array([[0., 1j], [1., 0.]]), np.diag([1., 1j]))
SIZE = 128
ROOT_MODES = (0, 2, 6)
ANCHOR = np.arange(8)


def _source(beta, length):
    vector = float(beta) ** np.arange(length)
    return vector / np.linalg.norm(vector)


def _scalar(value, polar):
    """Specify both occupied flag columns and the actual dirty rejection."""
    reject = np.sqrt(1 - abs(value) ** 2)
    return np.block([[value * np.eye(2), -reject * polar.conj().T],
                     [reject * polar, value.conjugate() * np.eye(2)]])


def _embed(local, indices, size=SIZE):
    result = np.eye(size, dtype=complex)
    result[np.ix_(indices, indices)] = local
    return result


def _first(flag, source, mode, dirty):
    return np.ravel_multi_index((flag, source, mode, dirty), (2, 2, 8, 2))


def _second(flag, source, mode, dirty):
    return 64 + np.ravel_multi_index((flag, source, ROOT_MODES.index(mode), dirty), (2, 4, 3, 2))


def _rare_lcu_refinement(beta, probability, old_length, polar):
    parameters = refinement_parameters(beta, probability)
    length = 2 * old_length
    shape = (2, 2, length, 2)  # rare, separate raw flag, source, dirty
    size = int(np.prod(shape))
    raw = np.eye(size, dtype=complex)
    for source in range(length):
        indices = [np.ravel_multi_index((1, flag, source, dirty), shape)
                   for flag in range(2) for dirty in range(2)]
        raw[np.ix_(indices, indices)] = _scalar(float(parameters.scalar_coefficient), polar)
    selected_x = np.zeros_like(raw)
    for column in range(size):
        rare, flag, source, dirty = np.unravel_index(column, shape)
        output = (rare, flag, source ^ rare, dirty)
        selected_x[np.ravel_multi_index(output, shape), column] = 1
    p = float(probability)
    preparation = np.array([[np.sqrt(1 - p), -np.sqrt(p)],
                            [np.sqrt(p), np.sqrt(1 - p)]])
    preparation = np.kron(preparation, np.eye(4 * length))
    return preparation.T @ selected_x @ raw @ preparation


def _refinement_dilation(beta, probability, survival, length, modes, polar):
    """Alternative full diagnostic dilation, diagonalized on the new bit.

    The standalone rare-LCU test independently checks the actual selector
    contract. This smaller dilation keeps the combined arena fixture bounded.
    """
    shape = (2, length, modes, 2)
    size = int(np.prod(shape))
    raw = np.zeros((size, size), dtype=complex)
    for source in range(length):
        value = survival * (1 - float(probability)) * (1 + (-1) ** (source % 2) * float(beta))
        for mode in range(modes):
            indices = [np.ravel_multi_index((flag, source, mode, dirty), shape)
                       for flag in range(2) for dirty in range(2)]
            raw[np.ix_(indices, indices)] = _scalar(value, polar)
    basis = np.kron(np.eye(2), np.kron(np.eye(length // 2), np.kron(H, np.eye(2 * modes))))
    return basis @ raw @ basis


def _fresh_permutation(reuse_rejection=False):
    permutation = np.eye(SIZE)
    for flag in ((0, 1) if reuse_rejection else (0,)):
        for source in range(2):
            for mode in ROOT_MODES:
                for dirty in range(2):
                    old = _first(flag, source, mode, dirty)
                    fresh = _second(flag, 2 * source, mode, dirty)
                    permutation[[old, fresh]] = permutation[[fresh, old]]
    return permutation


def _geometric_rotation(beta):
    beta = float(beta)
    return np.array([[1., -beta], [beta, 1.]]) / np.sqrt(1 + beta * beta)


def _tag_preparation(beta_first, beta_second):
    result = np.eye(SIZE, dtype=complex)
    first_rotation = _geometric_rotation(beta_first)
    second_rotation = np.kron(_geometric_rotation(beta_second ** 2),
                              _geometric_rotation(beta_second))
    for mode in range(8):
        if mode in ROOT_MODES:
            indices = [_second(0, source, mode, dirty) for source in range(4) for dirty in range(2)]
            local = second_rotation
        else:
            indices = [_first(0, source, mode, dirty) for source in range(2) for dirty in range(2)]
            local = first_rotation
        result = _embed(np.kron(local, np.eye(2)), indices) @ result
    return result


def _pword(fixture, gammas, permutation, tag_preparation):
    """Gamma-compensated actual replay and an independent private failure."""
    word = np.eye(SIZE, dtype=complex)
    nodes = {node['detail']: node for node in fixture['nodes'].values()}
    for label in range(4):
        node = nodes.get(label)
        x = fixture['coefficients'][label] / (1.6 * gammas[label])
        loss = node['loss'] if node is not None else 0.
        y = fixture['corrections'][label] / (1.6 * gammas[4 + label] * loss) if loss else 0j
        radial = np.sqrt(abs(x) ** 2 + loss * abs(y) ** 2)
        if radial >= 1:
            raise ValueError("fixture radial attenuation must be a contraction")
        indices = [2 * label + dirty for dirty in range(2)] + [120 + 2 * label + dirty for dirty in range(2)]
        word = _embed(_scalar(radial, np.eye(2)), indices) @ word
        if y == 0:
            replay = np.kron(np.diag([x / radial, 1]), np.eye(2))
        else:
            accepted, raw = node['accepted'], node['scalar']
            phase = np.exp(2j * np.arcsin(abs(y / radial) / (2 * abs(accepted))))
            first = 1 + (phase - 1) * abs(accepted) ** 2
            second = (phase - 1) * accepted.conjugate()
            phases = np.exp(1j * np.angle([x / (radial * first), y / (radial * second)]))
            replay = (np.kron(np.diag(phases), np.eye(2)) @ raw @
                      np.kron(np.diag([phase, 1]), np.eye(2)) @ raw.conj().T)
        indices = [_first(0, 0, mode, dirty) for mode in (label, 4 + label) for dirty in range(2)]
        word = _embed(replay, indices) @ word
    return tag_preparation @ permutation @ word


def _arena_fixture():
    fixture = _mixed_chart_fixture()
    beta_first, beta_second = Fraction(1, 256), Fraction(1, 16)
    probability_first, probability_second = Fraction(1, 256), Fraction(1, 16)
    survival_first, survival_second = 63 / 64, 31 / 32
    gamma_first = survival_first * np.sqrt(float(refinement_parameters(
        beta_first, probability_first).survival_squared))
    gamma_second = survival_second * np.sqrt(float(refinement_parameters(
        beta_second, probability_second).survival_squared))
    gammas = np.full(8, gamma_first)
    gammas[list(ROOT_MODES)] *= gamma_second
    wrong_gammas = np.full(8, survival_first * (1 - float(probability_first)))
    wrong_gammas[list(ROOT_MODES)] *= survival_second * (1 - float(probability_second))

    attenuation = np.eye(SIZE, dtype=complex)
    for label, diagonal in enumerate(fixture['diagonal']):
        indices = [2 * label + dirty for dirty in range(2)] + [112 + 2 * label + dirty for dirty in range(2)]
        attenuation = _embed(_scalar(diagonal.conjugate() / 1.25, np.eye(2)), indices) @ attenuation
    refine_first = _embed(_refinement_dilation(beta_first, probability_first, survival_first,
                                              2, 8, POLARS[0]), np.arange(64))
    bottom = np.eye(16, dtype=complex)
    for vertex in (3, 2):
        node = fixture['nodes'][vertex]
        bottom = _embed_modes(8, (node['left'], node['detail'], 4 + node['detail']),
                              node['local_word']) @ bottom
    bottom = _embed(np.kron(np.eye(2), bottom), np.arange(32))
    permutation = _fresh_permutation()
    refine_second = _embed(_refinement_dilation(beta_second, probability_second, survival_second,
                                               4, 3, POLARS[1]), np.arange(64, 112))
    root = _embed(np.kron(np.eye(4), fixture['nodes'][1]['local_word']), np.arange(64, 88))
    before_transfer = bottom @ refine_first @ attenuation
    later = root @ refine_second @ permutation
    qword = later @ before_transfer
    wrong_qword = root @ refine_second @ _fresh_permutation(True) @ before_transfer
    tag_preparation = _tag_preparation(beta_first, beta_second)
    pword = _pword(fixture, gammas, permutation, tag_preparation)

    tags = np.zeros((SIZE, 16), dtype=complex)
    for mode in range(8):
        length, beta = (4, beta_second) if mode in ROOT_MODES else (2, beta_first)
        for source, amplitude in enumerate(_source(beta, length)):
            for dirty in range(2):
                row = (_second(0, source, mode, dirty) if mode in ROOT_MODES
                       else _first(0, source, mode, dirty))
                tags[row, 2 * mode + dirty] = amplitude
    return dict(fixture=fixture, gammas=gammas, wrong_gammas=wrong_gammas,
                pword=pword, qword=qword, wrong_qword=wrong_qword, tags=tags,
                permutation=permutation, tag_preparation=tag_preparation,
                before_transfer=before_transfer, later=later)


class GrowingPrecisionSourceTests(unittest.TestCase):
    def test_rare_lcu_refinement_has_the_normalized_source_factor(self):
        beta, probability, old_length = Fraction(1, 8), Fraction(1, 4), 2
        parameters = refinement_parameters(beta, probability)
        for polar in POLARS:
            word = _rare_lcu_refinement(beta, probability, old_length, polar)
            np.testing.assert_allclose(word.conj().T @ word, np.eye(len(word)), atol=ATOL, rtol=0)
            accepted_size = 4 * old_length
            desired = (1 - float(probability)) * np.kron(
                np.eye(old_length), np.kron(np.eye(2) + float(beta) * X, np.eye(2)))
            np.testing.assert_allclose(word[:accepted_size, :accepted_size], desired, atol=ATOL, rtol=0)
            embedded = np.kron(_source(beta ** 2, old_length), [1., 0.])[:, None]
            embedded = np.kron(embedded, np.eye(2))
            refined = np.kron(_source(beta, 2 * old_length)[:, None], np.eye(2))
            actual = word[:accepted_size, :accepted_size] @ embedded
            normalized = np.sqrt(float(parameters.survival_squared)) * refined
            np.testing.assert_allclose(actual, normalized, atol=ATOL, rtol=0)
            self.assertGreater(np.linalg.norm(actual - (1 - float(probability)) * refined, 2), .005)
            dirty_z = np.kron(np.eye(len(word) // 2), np.diag([1., -1.]))
            if not np.allclose(polar, np.diag(np.diag(polar))):
                self.assertGreater(np.linalg.norm(word @ dirty_z - dirty_z @ word, 2), .1)
        self.assertGreater(np.linalg.norm(POLARS[0] @ POLARS[1] - POLARS[1] @ POLARS[0], 2), 1.)
        for length in (1, 2, 5):
            old_norm = sum((beta ** (4 * j) for j in range(length)), Fraction())
            new_norm = sum((beta ** (2 * j) for j in range(2 * length)), Fraction())
            self.assertEqual(new_norm / old_norm, 1 + beta ** 2)
            self.assertEqual(parameters.survival_squared,
                             (1 - probability) ** 2 * new_norm / old_norm)

    def test_fresh_permutation_seals_occupied_failures_and_frozen_features(self):
        arena = _arena_fixture()
        permutation, later = arena['permutation'], arena['later']
        np.testing.assert_array_equal(permutation.T @ permutation, np.eye(SIZE))
        live = [_first(0, source, mode, dirty) for source in range(2)
                for mode in ROOT_MODES for dirty in range(2)]
        fresh = [_second(0, 2 * source, mode, dirty) for source in range(2)
                 for mode in ROOT_MODES for dirty in range(2)]
        np.testing.assert_array_equal(permutation[:, live], np.eye(SIZE)[:, fresh])
        np.testing.assert_array_equal(permutation[:, fresh], np.eye(SIZE)[:, live])
        sealed = list(range(32, 64)) + [_first(0, source, mode, dirty)
                                      for source in range(2) for mode in (1, 3, 4, 5, 7)
                                      for dirty in range(2)]
        np.testing.assert_allclose(later[:, sealed], np.eye(SIZE)[:, sealed], atol=ATOL, rtol=0)
        # Entangle two frozen sectors, dirty values, and an external reference.
        entangled = np.zeros((SIZE, 2), dtype=complex)
        entangled[_first(1, 1, 2, 0), 0] = 1 / np.sqrt(2)
        entangled[_first(0, 0, 5, 1), 1] = 1j / np.sqrt(2)
        np.testing.assert_allclose(later @ entangled, entangled, atol=ATOL, rtol=0)
        self.assertGreater(np.linalg.norm(_fresh_permutation(True) @ entangled - entangled), .5)

    def test_frozen_gamma_compensation_uses_the_actual_dirty_replay(self):
        arena = _arena_fixture()
        fixture, tags, gammas = arena['fixture'], arena['tags'], arena['gammas']
        root_polar = fixture['nodes'][1]['rejected'] / np.sqrt(fixture['nodes'][1]['loss'])
        child_polar = fixture['nodes'][2]['rejected'] / np.sqrt(fixture['nodes'][2]['loss'])
        self.assertGreater(np.linalg.norm(root_polar @ child_polar - child_polar @ root_polar, 2), 1.)
        np.testing.assert_allclose(tags.conj().T @ tags, np.eye(16), atol=ATOL, rtol=0)
        self.assertGreaterEqual(min(gammas), 7 / 8)
        self.assertGreater(gammas[1], gammas[2])
        for detail in (1, 2, 3):
            self.assertEqual(gammas[detail], gammas[4 + detail])
        row_source = fixture['row_word'][:, :8] @ np.kron(
            np.diag(fixture['diagonal'].conj() / 1.25), np.eye(2))
        np.testing.assert_allclose(tags.conj().T @ arena['qword'][:, ANCHOR],
                                   np.diag(np.repeat(gammas, 2)) @ row_source, atol=ATOL, rtol=0)
        expected_p = np.zeros((16, 8), dtype=complex)
        expected_p[:8] = np.kron(np.diag(fixture['coefficients'] / (1.6 * gammas[:4])), np.eye(2))
        for node in fixture['nodes'].values():
            detail = node['detail']
            expected_p[2 * (4 + detail):2 * (5 + detail), 2 * detail:2 * detail + 2] = (
                fixture['corrections'][detail] * node['rejected'] /
                (1.6 * gammas[4 + detail] * node['loss']))
        np.testing.assert_allclose(tags.conj().T @ arena['pword'][:, ANCHOR], expected_p, atol=ATOL, rtol=0)
        target = np.kron(fixture['target'], np.eye(2))
        np.testing.assert_allclose(arena['qword'][:, ANCHOR].conj().T @ arena['pword'][:, ANCHOR],
                                   target / 2, atol=ATOL, rtol=0)
        wrong_p = _pword(fixture, arena['wrong_gammas'], arena['permutation'], arena['tag_preparation'])
        self.assertGreater(np.linalg.norm(arena['qword'][:, ANCHOR].conj().T @
                                          wrong_p[:, ANCHOR] - target / 2, 2), 1e-4)
        self.assertGreater(np.linalg.norm(arena['wrong_qword'][:, ANCHOR].conj().T @
                                          arena['pword'][:, ANCHOR] - target / 2, 2), 1e-4)

    def test_original_anchor_and_actual_inverse_return_every_arena_coordinate(self):
        arena = _arena_fixture()
        for word in (arena['pword'], arena['qword'], arena['wrong_qword']):
            np.testing.assert_allclose(word.conj().T @ word, np.eye(SIZE), atol=ATOL, rtol=0)
        overlap = arena['qword'].conj().T @ arena['pword']
        reflection = np.eye(SIZE)
        reflection[ANCHOR, ANCHOR] = -1
        target = np.eye(SIZE)[:, ANCHOR] @ np.kron(arena['fixture']['target'], np.eye(2))
        recovered = (-overlap @ reflection @ overlap.conj().T @ reflection @ overlap)[:, ANCHOR]
        np.testing.assert_allclose(recovered, target, atol=3 * ATOL, rtol=0)
        dirty_reflection = np.eye(SIZE)
        dirty_reflection[ANCHOR[::2], ANCHOR[::2]] = -1
        final_source = arena['tags'][:, :8]
        final_reflection = np.eye(SIZE) - 2 * final_source @ final_source.conj().T
        for wrong in (dirty_reflection, final_reflection):
            failed = (-overlap @ wrong @ overlap.conj().T @ wrong @ overlap)[:, ANCHOR]
            self.assertGreater(np.linalg.norm(failed - target, 2), .1)
        wrong_inverse = (-overlap @ reflection @ overlap @ reflection @ overlap)[:, ANCHOR]
        self.assertGreater(np.linalg.norm(wrong_inverse - target, 2), .1)

    def test_accepted_error_gives_linear_full_recovery_error_including_rejection(self):
        # An explicit two-bank full unitary isolates the accepted-only promise.
        # Its rejection polar acts on arbitrary dirty inputs; completing only
        # the accepted polynomial would miss the leading output error.
        logical = np.array([[np.cos(.37), -np.sin(.37)],
                            [np.sin(.37), np.cos(.37)]])
        target = np.kron(logical, np.eye(2))
        polar = np.zeros((4, 4), dtype=complex)
        polar[:2, :2], polar[2:, 2:] = POLARS
        reflection = np.diag([-1.] * 4 + [1.] * 4)
        wanted = np.vstack((target, np.zeros((4, 4))))
        for delta in (0., 1 / 32, 1 / 128, 1 / 1024):
            amplitude = .5 + delta
            rejected = np.sqrt(1 - amplitude ** 2)
            word = np.block([[amplitude * target, -rejected * target @ polar.conj().T],
                             [rejected * polar, amplitude * np.eye(4)]])
            np.testing.assert_allclose(word.conj().T @ word, np.eye(8), atol=ATOL, rtol=0)
            recovered = (-word @ reflection @ word.conj().T @ reflection @ word)[:, :4]
            accepted_polynomial = 3 * amplitude - 4 * amplitude ** 3
            np.testing.assert_allclose(recovered[:4], accepted_polynomial * target, atol=ATOL, rtol=0)
            np.testing.assert_allclose(recovered[4:],
                                       rejected * (1 - 4 * amplitude ** 2) * polar, atol=ATOL, rtol=0)
            self.assertLessEqual(np.linalg.norm(recovered - wanted, 2), 12 * delta + ATOL)
            if delta:
                self.assertGreater(np.linalg.norm(recovered[4:], 2), delta)
                self.assertGreater(np.linalg.norm(recovered - wanted, 2),
                                   10 * np.linalg.norm(recovered[:4] - target, 2))

    def test_phase_dressing_reuses_leaked_flags_and_dirty_core(self):
        # Diagnostic real recovery followed by literal packed-diagonal words.
        # The addressed direct sum is supplied, not a native table emitter;
        # this checks composition, not the endpoint's precision or arena cost.
        # Wires: dirty core 0,1; logical bit 2; shared clean flags 3,4.
        logical = np.array([[np.cos(.37), -np.sin(.37)],
                            [np.sin(.37), np.cos(.37)]])
        real_target = np.kron(logical, np.eye(4))
        polar = np.zeros((8, 8), dtype=complex)
        polar[:4, :4] = np.kron(POLARS[0], X)
        polar[4:, 4:] = np.kron(POLARS[1], H)
        amplitude = .5 + 1 / 512
        rejected = np.sqrt(1 - amplitude ** 2)
        word = np.block([[amplitude * real_target,
                          -rejected * real_target @ polar.conj().T],
                         [rejected * polar, amplitude * np.eye(8)]])
        reflection = np.diag([-1.] * 8 + [1.] * 8)
        real = np.kron(np.eye(2), -word @ reflection @ word.conj().T @ reflection @ word)
        embedding = np.eye(32)[:, :8]
        np.testing.assert_allclose(real.conj().T @ real, np.eye(32), atol=ATOL, rtol=0)
        real_output = real @ embedding
        real_error = np.linalg.norm(real_output - embedding @ real_target, 2)
        self.assertLess(real_error, .01)
        projected = embedding @ embedding.conj().T @ real_output
        self.assertGreater(np.linalg.norm(real_output - projected, 2), .005)

        for phases in (((Fraction(3, 5), Fraction(4, 5)),
                        (Fraction(-5, 13), Fraction(12, 13))),
                       ((Fraction(0), Fraction(1)), (Fraction(0), Fraction(1)))):
            with self.subTest(phases=phases):
                blocks, inverses = [], []
                for pair in phases:
                    signs = [_sign_word(round((1 - value) * 4), 3) for value in pair]
                    native = _diagonal_word(2, *signs)
                    flag_reflection = _reflection_word(2)
                    minus = [('X', 0), ('Z', 0), ('X', 0), ('Z', 0)]
                    amplified = (native + flag_reflection + _adjoint(native)
                                 + flag_reflection + native + minus)
                    blocks.append(_word_matrix(4, amplified))
                    inverses.append(_word_matrix(4, _adjoint(amplified)))
                diagonal = _embed_addressed_blocks(
                    5, [0, 1, 3, 4], lambda values: blocks[values[2]])
                diagonal_inverse = _embed_addressed_blocks(
                    5, [0, 1, 3, 4], lambda values: inverses[values[2]])
                phases_matrix = np.diag([complex(float(c), float(s)) for c, s in phases])
                diagonal_target = np.kron(phases_matrix, np.eye(4))
                target = np.kron(phases_matrix @ logical, np.eye(4))
                diagonal_error = np.linalg.norm(
                    diagonal @ embedding - embedding @ diagonal_target, 2)
                self.assertLess(diagonal_error, .25)
                actual = diagonal @ real_output
                error = np.linalg.norm(actual - embedding @ target, 2)
                self.assertLessEqual(error, real_error + diagonal_error + 20 * ATOL)
                self.assertGreater(np.linalg.norm(actual - diagonal @ projected, 2), .005)

                inverse = real.conj().T @ diagonal_inverse
                np.testing.assert_allclose(inverse @ actual, embedding, atol=20 * ATOL, rtol=0)
                self.assertLessEqual(np.linalg.norm(
                    inverse @ embedding - embedding @ target.conj().T, 2), error + 20 * ATOL)
                # Eight reference labels give full Schmidt support on every
                # logical/dirty input, retaining coherent flag leakage.
                coefficients = np.arange(1, 9) + 1j * np.arange(8, 0, -1)
                coefficients /= np.linalg.norm(coefficients)
                self.assertLessEqual(np.linalg.norm(
                    (actual - embedding @ target) * coefficients), error + 20 * ATOL)
                if phases[0] != phases[1]:
                    diagonal_output = diagonal @ embedding
                    self.assertGreater(np.linalg.norm(
                        diagonal_output - embedding @ embedding.conj().T @ diagonal_output, 2), .1)
                    wrong_order = real @ diagonal @ embedding
                    self.assertGreater(np.linalg.norm(actual - wrong_order, 2), .3)
                else:
                    # A common phase is part of the literal operator contract.
                    np.testing.assert_allclose(diagonal @ embedding, 1j * embedding,
                                               atol=20 * ATOL, rtol=0)
                    self.assertGreater(np.linalg.norm(actual - embedding @ real_target, 2), 1.)

    def test_dyadic_schedule_refinement_conserves_source_precision_volume(self):
        for endpoint in (64, 257, (1 << 80) + 3):
            heights = (0, 1, 2, 3, 4, 7, 8, 15, 16, 31, 32)
            previous = None
            for height in heights:
                scale = growing_source_scale(endpoint, height)
                self.assertTrue(scale.integral_exponent)
                # M*q=4N is an exponent identity, not a huge precision state.
                self.assertEqual(scale.source_length.bit_length() - 1 + scale.exponent_log2,
                                 endpoint + 2)
                self.assertEqual(scale.rare_probability * 8 * scale.order, 1)
                if previous is not None:
                    steps = scale.refinement_steps_from(previous)
                    self.assertEqual(scale.source_length, previous.source_length << steps)
                    self.assertEqual(steps % 4, 0)
                previous = scale
        before, after = growing_source_scale(64, 7), growing_source_scale(64, 8)
        self.assertEqual(after.order, 4 * before.order)
        self.assertEqual(after.refinement_steps_from(before), 4)
        # Four actual normalized two-way refinements tensor to the longer
        # geometric vector, including every normalization factor.
        beta = Fraction(1, 2)
        vector = np.ones(1)
        for power in (8, 4, 2, 1):
            vector = np.kron(vector, _source(beta ** power, 2))
        np.testing.assert_allclose(vector, _source(beta, 16), atol=ATOL, rtol=0)


if __name__ == '__main__':
    unittest.main()
