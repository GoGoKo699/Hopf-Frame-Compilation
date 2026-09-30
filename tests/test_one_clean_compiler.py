"""Small native audits of the paired-Majorana one-clean rotation compiler.

Every source plane, routing gate, Pauli mask, query, and amplification inverse
is a literal Clifford+T word. Matrices have dimension at most 128 and include
all dirty input columns. These finite checks do not establish resource bounds.
"""
from __future__ import annotations

import unittest

import numpy as np

from tests.test_operator_source_compiler import (
    ATOL, I2, X, Y, Z, _adjoint, _dirty_mask_lookup, _expand_toffolis,
    _pauli, _rotation, _sign_word, _word_matrix,
)


def _pauli_t_word(factors, sign=1):
    """A literal exp(-i sign*pi*P/8), up to a scalar canceled in M."""
    basis = []
    for qubit, factor in sorted(factors.items()):
        if factor == 'X':
            basis.append(('H', qubit))
        elif factor == 'Y':
            basis.extend([('SDG', qubit), ('H', qubit)])
    support = sorted(factors)
    parity = [('CX', qubit, support[-1]) for qubit in support[:-1]]
    return (basis + parity + [('T' if sign == 1 else 'TDG', support[-1])]
            + _adjoint(parity) + _adjoint(basis))


def _paired_source_word(q):
    """Branch once, then spread along the two constant-support Majorana chains."""
    assert q >= 2
    word = [('T', 0)]
    word += _pauli_t_word({0: 'Y', 1: 'X'}, -1)
    for qubit in range(q - 1):
        word += _pauli_t_word({qubit: 'X', qubit + 1: 'Y'})
    for qubit in range(1, q):
        word += _pauli_t_word({qubit: 'Y', qubit + 1: 'X'}, -1)
    return word


def _paired_weights(q):
    geometric = np.array([2.0 ** (-j - 1) for j in range(q - 1)]
                         + [2.0 ** (1 - q)])
    weights = np.zeros(2 * (q + 1))
    weights[0] = 1 / 4
    weights[1:2 * q:2] = geometric / 2
    weights[2:2 * q + 1:2] = geometric / 4
    fixed = np.zeros(2 * (q + 1), dtype=int)
    fixed[2:2 * q + 1:2] = 1
    return weights, fixed


def _paired_source_data(q):
    core = q + 1
    gammas = [_pauli(core, {**{k: Z for k in range(qubit)}, qubit: factor})
              for qubit in range(core) for factor in (X, Y)]
    weights, fixed = _paired_weights(q)
    source_word = _paired_source_word(q)
    source = _word_matrix(core, _adjoint(source_word) + [('X', 0)] + source_word)
    return source, gammas, weights, fixed


def _mask_bits(signs):
    """Pauli exponents from the paired signs, including the running X prefix."""
    xmask = zmask = prefix = 0
    for qubit in range(len(signs) // 2):
        first, second = signs[2 * qubit:2 * qubit + 2]
        xbit = int(first ^ second)
        zbit = int(first ^ prefix)
        xmask |= xbit << qubit
        zmask |= zbit << qubit
        prefix ^= xbit
    return xmask, zmask


def _literal_mask(signs):
    xmask, zmask = _mask_bits(signs)
    core = len(signs) // 2
    return ([('X', qubit) for qubit in range(core) if (xmask >> qubit) & 1]
            + [('Z', qubit) for qubit in range(core) if (zmask >> qubit) & 1])


def _encoded_signs(q, theta):
    """Two rounded geometric tails, with a signed quarter-weight head."""
    beta = np.sin(np.pi / 10)
    desired_s, desired_p = beta * np.cos(theta), beta * np.sin(theta)
    plus = 3 * desired_s / 4 - desired_p / 2
    minus = desired_s / 4 + desired_p / 2
    head = 1 if plus >= 0 else -1
    values = (2 * (plus - head / 4), 4 * minus)
    signs = np.zeros(2 * (q + 1), dtype=int)
    signs[0] = int(head == -1)
    scale = 1 << (q - 1)
    for start, value in zip((1, 2), values):
        integer = int(np.clip(np.rint((1 - value) * scale / 2), 0, scale))
        signs[start:2 * q + 1:2] = _sign_word(integer, q - 1)
    return signs


def _controlled_x(controls, target, helper=None):
    if len(controls) <= 2:
        return [(('X', 'CX', 'CCX')[len(controls)], *controls, target)]
    assert len(controls) == 3 and helper is not None
    first, second, third = controls
    return [('CCX', first, second, helper), ('CCX', helper, third, target)] * 2


def _scalar_word(q, mask_word, flag, predicate=(), helper=None):
    source = _paired_source_word(q)
    inverse = _adjoint(source)
    mword = inverse + _controlled_x(list(predicate), 0, helper) + source
    center = _controlled_x([*predicate, flag], 0, helper)
    control1 = inverse + center + source
    control0 = inverse + [('X', flag)] + center + [('X', flag)] + source
    return ([('H', flag)] + control1 + _adjoint(mask_word) + mword
            + mask_word + control0 + [('H', flag)])


def _sandwich_word(q, mask_word, target, flag, phase=False, predicate=(), helper=None):
    fixed = _paired_source_data(q)[3]
    scalar_f = _scalar_word(q, _literal_mask(fixed), flag, predicate, helper)
    scalar_g = _scalar_word(q, mask_word, flag, predicate, helper)
    if phase:
        left = scalar_f
        middle = [('SDG', flag)] + scalar_g + [('S', flag)]
    else:
        controlled_z = [('H', target), ('CX', flag, target), ('H', target)]
        left = controlled_z + scalar_f + controlled_z
        middle = [('CX', flag, target)] + scalar_g + [('CX', flag, target)]
    # Chronological order; the resulting matrix is A^dagger B A.
    return left + middle + _adjoint(left)


def _amplification_word(word, flag):
    # The four factors R=-Z have the same combined literal phase as four Zs.
    reflection = [('Z', flag)]
    return (word + reflection + _adjoint(word) + reflection + word
            + reflection + _adjoint(word) + reflection + word)


def _coefficients(weights, fixed, signs):
    signed = 1 - 2 * signs
    s = float(weights @ signed)
    t = float(weights @ ((1 - 2 * fixed) * signed))
    return s, s / 2 - t


def _mask_table_word(core, address, selector, sign_rows):
    masks = [_mask_bits(signs) for signs in sign_rows]

    def query(table):
        if len(address) == 2:
            return _dirty_mask_lookup(address, list(range(core)), selector, table)
        word = []
        for row, mask in enumerate(table):
            negate = [('X', address[0])] if row == 0 else []
            word += negate + [('CX', address[0], qubit)
                              for qubit in range(core) if (mask >> qubit) & 1] + negate
        return word

    hadamards = [('H', qubit) for qubit in range(core)]
    return (query([mask[0] for mask in masks]) + hadamards
            + query([mask[1] for mask in masks]) + hadamards)


def _conditional_mask_word(signs, predicate, helper):
    """Each X/Z part has its own dirty-control echo, retaining literal phase."""
    core = len(signs) // 2
    xmask, zmask = _mask_bits(signs)
    toggle = [('CX', predicate, helper)]

    def echo(mask):
        query = [('CX', helper, qubit) for qubit in range(core)
                 if (mask >> qubit) & 1]
        return toggle + query + toggle + query

    hadamards = [('H', qubit) for qubit in range(core)]
    return echo(xmask) + hadamards + echo(zmask) + hadamards


class OneCleanCompilerTests(unittest.TestCase):
    def test_native_paired_source_and_every_small_prefix_mask(self):
        for q in (2, 3):
            core = q + 1
            source, gammas, weights, fixed = _paired_source_data(q)
            expected = sum(np.sqrt(weight) * gamma
                           for weight, gamma in zip(weights, gammas))
            np.testing.assert_allclose(source, expected, atol=ATOL, rtol=0)
            np.testing.assert_allclose(source @ source, np.eye(1 << core),
                                       atol=ATOL, rtol=0)
            self.assertAlmostEqual(float(weights @ (1 - 2 * fixed)), 1 / 2)
            # All masks include the unused final Majorana: the prefix formula
            # must implement every sign pattern, not only a source subset.
            for pattern in range(1 << (2 * core)):
                signs = np.array([(pattern >> j) & 1 for j in range(2 * core)])
                mask = _word_matrix(core, _literal_mask(signs))
                for index, gamma in enumerate(gammas):
                    np.testing.assert_allclose(mask @ gamma @ mask.conj().T,
                                               (-1) ** signs[index] * gamma,
                                               atol=ATOL, rtol=0)

    def test_native_sandwich_and_five_call_amplification_on_all_dirty_inputs(self):
        for q in (2, 3):
            core, target, flag = q + 1, q + 1, q + 2
            width = flag + 1
            source, _, weights, fixed = _paired_source_data(q)
            for theta in (0.0, np.pi / 2, -np.pi / 3, 2.1):
                with self.subTest(q=q, theta=theta):
                    signs = _encoded_signs(q, theta)
                    s, p = _coefficients(weights, fixed, signs)
                    mask = _word_matrix(core, _literal_mask(signs))
                    fixed_mask = _word_matrix(core, _literal_mask(fixed))
                    nf = fixed_mask @ source @ fixed_mask.conj().T
                    ng = mask @ source @ mask.conj().T
                    df = source @ nf - np.eye(1 << core) / 2
                    dg = source @ ng - s * np.eye(1 << core)
                    np.testing.assert_allclose(df @ dg + dg @ df,
                                               2 * p * np.eye(1 << core),
                                               atol=ATOL, rtol=0)

                    word = _sandwich_word(q, _literal_mask(signs), target, flag)
                    unitary = _word_matrix(width, word)
                    accepted_size = 1 << flag
                    block = np.kron(s * I2 + p * (X @ Z), np.eye(1 << core))
                    np.testing.assert_allclose(unitary[:accepted_size, :accepted_size],
                                               block, atol=ATOL, rtol=0)
                    np.testing.assert_allclose(unitary.conj().T @ unitary,
                                               np.eye(1 << width), atol=ATOL, rtol=0)
                    amplified = _word_matrix(width, _amplification_word(word, flag))
                    radius2 = s * s + p * p
                    polynomial = 5 - 20 * radius2 + 16 * radius2 ** 2
                    np.testing.assert_allclose(amplified[:accepted_size, :accepted_size],
                                               polynomial * block, atol=ATOL, rtol=0)
                    expected_target = np.kron(_rotation(theta), np.eye(1 << core))
                    embedding = np.eye(1 << width)[:, :accepted_size]
                    error = amplified @ embedding - embedding @ expected_target
                    # This identity checks the full dirty-space isometry error,
                    # including coherent rejected-flag leakage and literal phase.
                    error2 = 2 - 2 * polynomial * (s * np.cos(theta) + p * np.sin(theta))
                    np.testing.assert_allclose(error.conj().T @ error,
                                               error2 * np.eye(accepted_size),
                                               atol=ATOL, rtol=0)

    def test_native_prefix_mask_table_preserves_arbitrary_dirty_selector(self):
        q, core = 2, 3
        address, selector, width = [3, 4], 5, 6
        rows = [_encoded_signs(q, theta) for theta in (0.1, -0.9, 1.7, 3.0)]
        word = _expand_toffolis(_mask_table_word(core, address, selector, rows))
        actual = _word_matrix(width, word)
        expected = np.zeros_like(actual)
        for basis in range(1 << width):
            row = sum(((basis >> qubit) & 1) << j for j, qubit in enumerate(address))
            xmask, zmask = _mask_bits(rows[row])
            image = basis ^ xmask
            expected[image, basis] = (-1) ** ((image & zmask).bit_count())
        np.testing.assert_allclose(actual, expected, atol=ATOL, rtol=0)
        inverse = _word_matrix(width, _adjoint(word))
        np.testing.assert_allclose(inverse @ actual, np.eye(1 << width), atol=ATOL, rtol=0)

    def test_addressed_native_rotation_table_has_no_relative_row_phase(self):
        q, core, target, address, flag, width = 2, 3, 3, 4, 5, 6
        _, _, weights, fixed = _paired_source_data(q)
        rows = [_encoded_signs(q, theta) for theta in (-0.7, 1.8)]
        mask_word = _mask_table_word(core, [address], None, rows)
        word = _sandwich_word(q, mask_word, target, flag)
        actual = _word_matrix(width, word)
        expected = np.zeros((1 << flag, 1 << flag), dtype=complex)
        for row, signs in enumerate(rows):
            s, p = _coefficients(weights, fixed, signs)
            block = np.kron(s * I2 + p * (X @ Z), np.eye(1 << core))
            low, high = row << address, (row + 1) << address
            expected[low:high, low:high] = block
        np.testing.assert_allclose(actual[:1 << flag, :1 << flag], expected,
                                   atol=ATOL, rtol=0)

    def test_dirty_predicate_echo_is_exactly_inactive_on_both_flag_inputs(self):
        q, core, target, helper, predicate, flag, width = 2, 3, 3, 4, 5, 6, 7
        signs = _encoded_signs(q, 0.9)
        _, _, weights, fixed = _paired_source_data(q)
        s, p = _coefficients(weights, fixed, signs)
        mask_word = _conditional_mask_word(signs, predicate, helper)
        word = _sandwich_word(q, mask_word, target, flag)
        actual = _word_matrix(width, word)
        amplified = _word_matrix(width, _amplification_word(word, flag))
        inactive = [basis for basis in range(1 << width) if not (basis >> predicate) & 1]
        identity = np.eye(1 << width)
        np.testing.assert_allclose(actual[:, inactive], identity[:, inactive],
                                   atol=ATOL, rtol=0)
        np.testing.assert_allclose(amplified[:, inactive], identity[:, inactive],
                                   atol=ATOL, rtol=0)
        # Active compression includes both arbitrary helper inputs. Neither the
        # source nor its flag is substituted by an initialized-state shortcut.
        active = [basis for basis in range(1 << flag) if (basis >> predicate) & 1]
        block = np.kron(I2, np.kron(s * I2 + p * (X @ Z), np.eye(1 << core)))
        np.testing.assert_allclose(actual[np.ix_(active, active)], block,
                                   atol=ATOL, rtol=0)
        for helper_input in (0, 1):
            inputs = [basis for basis in range(1 << width)
                      if (basis >> helper) & 1 == helper_input]
            forbidden = [basis for basis in range(1 << width)
                         if (basis >> helper) & 1 != helper_input]
            np.testing.assert_allclose(actual[np.ix_(forbidden, inputs)], 0,
                                       atol=ATOL, rtol=0)

    def test_flag_phase_route_has_the_literal_scalar_phase(self):
        for q in (2, 3):
            core, flag, width = q + 1, q + 1, q + 2
            _, _, weights, fixed = _paired_source_data(q)
            for phase in (-1.2, 0.0, 0.8):
                with self.subTest(q=q, phase=phase):
                    signs = _encoded_signs(q, -phase)
                    s, p = _coefficients(weights, fixed, signs)
                    word = _sandwich_word(q, _literal_mask(signs), None, flag, phase=True)
                    actual = _word_matrix(width, word)
                    np.testing.assert_allclose(actual[:1 << core, :1 << core],
                                               (s - 1j * p) * np.eye(1 << core),
                                               atol=ATOL, rtol=0)
                    amplified = _word_matrix(width, _amplification_word(word, flag))
                    radius2 = s * s + p * p
                    polynomial = 5 - 20 * radius2 + 16 * radius2 ** 2
                    np.testing.assert_allclose(amplified[:1 << core, :1 << core],
                                               polynomial * (s - 1j * p) * np.eye(1 << core),
                                               atol=ATOL, rtol=0)

    def test_native_source_center_predicate_returns_its_borrowed_helper(self):
        q, core = 2, 3
        signs = _encoded_signs(q, -0.8)
        weights, fixed = _paired_weights(q)
        s, p = _coefficients(weights, fixed, signs)
        # Rotation: one unchanged predicate. Scalar phase: two predicates and
        # an arbitrary helper, whose C3X is expanded into exact native gates.
        for phase, target, predicate, helper, flag in (
                (False, 3, (4,), None, 5), (True, None, (3, 4), 5, 6)):
            with self.subTest(phase=phase):
                width = flag + 1
                word = _expand_toffolis(_sandwich_word(
                    q, _literal_mask(signs), target, flag, phase, predicate, helper))
                actual = _word_matrix(width, word)
                amplified = _word_matrix(width, _amplification_word(word, flag))
                inactive = [basis for basis in range(1 << width)
                            if not all((basis >> bit) & 1 for bit in predicate)]
                identity = np.eye(1 << width)
                for circuit in (actual, amplified):
                    np.testing.assert_allclose(circuit[:, inactive], identity[:, inactive],
                                               atol=ATOL, rtol=0)
                active = [basis for basis in range(1 << flag)
                          if all((basis >> bit) & 1 for bit in predicate)]
                block = ((s - 1j * p) * np.eye(len(active)) if phase else
                         np.kron(s * I2 + p * (X @ Z), np.eye(1 << core)))
                np.testing.assert_allclose(actual[np.ix_(active, active)], block,
                                           atol=ATOL, rtol=0)
                if helper is not None:
                    for helper_input in (0, 1):
                        inputs = [basis for basis in range(1 << width)
                                  if (basis >> helper) & 1 == helper_input]
                        forbidden = [basis for basis in range(1 << width)
                                     if (basis >> helper) & 1 != helper_input]
                        for circuit in (actual, amplified):
                            np.testing.assert_allclose(circuit[np.ix_(forbidden, inputs)],
                                                       0, atol=ATOL, rtol=0)

    def test_fine_precision_full_isometry_in_both_small_clifford_algebra_representations(self):
        # The Gram data of the three source reflections has a two-dimensional
        # representation. This checks fine precision without constructing a
        # large precision core; it is not a native source-circuit simulation.
        maxima = []
        for q in (12, 20):
            weights, fixed = _paired_weights(q)
            errors = []
            for theta in (0.0, np.pi / 2, -np.pi / 3, 2.1, np.pi):
                signs = _encoded_signs(q, theta)
                s, p = _coefficients(weights, fixed, signs)
                t = s / 2 - p
                c = 1 / 2
                transverse = (t - c * s) / np.sqrt(1 - c * c)
                square = 1 - s * s - transverse * transverse
                self.assertGreaterEqual(square, -ATOL)
                for chirality in (-1, 1):
                    with self.subTest(q=q, theta=theta, chirality=chirality):
                        nf = c * X + np.sqrt(1 - c * c) * Z
                        ng = s * X + transverse * Z + chirality * np.sqrt(max(0, square)) * Y
                        df, dg = X @ nf - c * I2, X @ ng - s * I2
                        left = c * np.eye(8) + np.kron(X, np.kron(Z, df))
                        middle = s * np.eye(8) + np.kron(X, np.kron(X, dg))
                        unitary = left.conj().T @ middle @ left
                        reflection = np.kron(-Z, np.eye(4))
                        amplified = (unitary @ reflection @ unitary.conj().T @ reflection
                                     @ unitary @ reflection @ unitary.conj().T @ reflection
                                     @ unitary)
                        embedding = np.eye(8)[:, :4]
                        expected = embedding @ np.kron(_rotation(theta), I2)
                        error = np.linalg.norm(amplified @ embedding - expected, ord=2)
                        self.assertLessEqual(error, 30 * 2.0 ** (-q))
                        errors.append(error)
                        phase_left = c * np.eye(4) + np.kron(X, df)
                        phase_middle = s * np.eye(4) + np.kron(Y, dg)
                        phase_unitary = phase_left.conj().T @ phase_middle @ phase_left
                        phase_reflection = np.kron(-Z, I2)
                        phase_amplified = (
                            phase_unitary @ phase_reflection @ phase_unitary.conj().T
                            @ phase_reflection @ phase_unitary @ phase_reflection
                            @ phase_unitary.conj().T @ phase_reflection @ phase_unitary)
                        phase_embedding = np.eye(4)[:, :2]
                        phase_error = np.linalg.norm(
                            phase_amplified @ phase_embedding
                            - np.exp(-1j * theta) * phase_embedding, ord=2)
                        self.assertLessEqual(phase_error, 30 * 2.0 ** (-q))
                        self.assertAlmostEqual(phase_error, error, delta=ATOL)
            maxima.append(max(errors))
        self.assertLess(maxima[1], maxima[0] / 64)


if __name__ == '__main__':
    unittest.main()
