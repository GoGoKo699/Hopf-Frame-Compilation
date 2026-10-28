"""Exact grading and Clifford-dual checks for the four-address flat bank.

This finite native pilot is not an all-size field compiler.  Literal Pauli
phases are retained in every word and every partial Clifford homomorphism.
"""
from fractions import Fraction
from itertools import permutations
import unittest

from tests.exact_matrices import CI, CO, CQ2, CZ, adjoint, native_columns
from tests.test_flat_field_bank_native import (
    AXES, ORBIT, TARGET, full_intertwiner, ledger,
)
from tests.test_field_copy_clifford_native import _fixed_query
from tests.test_joint_lookup_pauli_source import _inverse


def _pairing(left, right, width):
    mask = (1 << width) - 1
    return (((left & mask) & (right >> width)).bit_count()
            + ((left >> width) & (right & mask)).bit_count()) % 2


def _hermitian_phase(label, width):
    return ((label & ((1 << width) - 1)) & (label >> width)).bit_count()


def _ordered_phase(left, right, width):
    # P(x,z)=X^x Z^z, without an implicit quotient by scalar phase.
    return 2 * ((left >> width) & (right & ((1 << width) - 1))).bit_count()


def _multiply_pauli_expansions(left, right, width=4):
    result = {}
    for a, coefficient in left.items():
        for b, other in right.items():
            exponent = (_hermitian_phase(a, width) + _hermitian_phase(b, width)
                        - _hermitian_phase(a ^ b, width)
                        + _ordered_phase(a, b, width)) % 4
            value = coefficient * other * (CO, CI, -CO, -CI)[exponent]
            label = a ^ b
            result[label] = result.get(label, CZ) + value
    return {label: value for label, value in result.items() if value != CZ}


def _extend_clifford_map(mapping, source, target, phase, width):
    """Extend a map on the entire generated Pauli-label subgroup.

    mapping[v]=(k,w) means P(v) maps to i^k P(w), for the ordered XZ
    section.  Storing the complete small subgroup makes every dependency
    explicit, independently of a choice of elimination pivots.
    """
    phase %= 4
    if source in mapping:
        return mapping if mapping[source] == (phase, target) else None
    if any(_pairing(a, source, width) != _pairing(b, target, width)
           for a, (_, b) in mapping.items()):
        return None
    extended = dict(mapping)
    for a, (alpha, b) in mapping.items():
        image = ((alpha + phase + _ordered_phase(b, target, width)
                  - _ordered_phase(a, source, width)) % 4, b ^ target)
        label = a ^ source
        if label in extended and extended[label] != image:
            return None
        extended[label] = image
    if len({b for _, b in extended.values()}) != len(extended):
        return None
    return extended


def _signed_supports(width=4):
    result = []
    for operator in ORBIT:
        row = {}
        for label, coefficient in operator.items():
            lifted = (label & 15) | ((label >> 4) << width)
            # Every coefficient is exactly +/-1/2, with no sqrt(2) part.
            assert coefficient.b == 0 and abs(coefficient.a) == Fraction(1, 2)
            row[lifted] = int(coefficient.a < 0)
        result.append(row)
    return result


def _clifford_characters(mask, width=4):
    states = [{0: (0, 0)}]
    counts = []
    for address, support in enumerate(_signed_supports(width)):
        labels = sorted(support)
        next_states = {}
        for mapping in states:
            for permutation in permutations(labels):
                current = mapping
                for source, target in zip(labels, permutation):
                    sign = ((mask >> address) & 1) ^ support[source] ^ support[target]
                    phase = (2 * sign + _hermitian_phase(target, width)
                             - _hermitian_phase(source, width))
                    current = _extend_clifford_map(
                        current, source, target, phase, width)
                    if current is None:
                        break
                if current is not None:
                    next_states[tuple(sorted(current.items()))] = current
        states = list(next_states.values())
        counts.append(len(states))
    return states, counts


def _wire_translation(address):
    word = []
    for left in range(4):
        right = left ^ address
        if left < right:
            word += [('CX', left, right), ('CX', right, left), ('CX', left, right)]
    return word


class FlatFieldBankDualTests(unittest.TestCase):
    def test_complete_commuting_grading_and_literal_full_product(self):
        orbit = [{label: CQ2(value) for label, value in operator.items()}
                 for operator in ORBIT]
        for i, left in enumerate(orbit):
            self.assertEqual(_multiply_pauli_expansions(left, left), {0: CO})
            for right in orbit[i + 1:]:
                self.assertEqual(_multiply_pauli_expansions(left, right),
                                 _multiply_pauli_expansions(right, left))
        # A literal Pauli image with four nonzero terms certifies that
        # the seed itself is genuinely non-Clifford.
        image = _multiply_pauli_expansions(
            _multiply_pauli_expansions(orbit[0], {2: CO}), orbit[0])
        self.assertEqual(image, {37: -CO * Fraction(1, 2),
                                 39: -CO * Fraction(1, 2),
                                 35: CO * Fraction(1, 2),
                                 33: -CO * Fraction(1, 2)})
        for mask in range(1, 16):
            product = {0: CO}
            for address, operator in enumerate(orbit):
                if (mask >> address) & 1:
                    product = _multiply_pauli_expansions(product, operator)
            self.assertNotIn(0, product)
            if mask == 15:
                self.assertEqual(product, {8: CO})  # +X_v1, not -X_v1.
        # The exact simultaneous Pauli images are a full virtual bank.
        labels = [next(iter(operator)) for operator in TARGET]
        span = {0}
        for label in labels:
            span |= {value ^ label for value in tuple(span)}
        self.assertEqual(len(span), 16)

    def test_every_clifford_sign_character_and_singleton_rejection(self):
        for mask in range(16):
            states, counts = _clifford_characters(mask)
            self.assertEqual(counts[:3], [24, 8, 4])
            self.assertEqual(len(states), 0 if mask.bit_count() % 2 else 4)
            for mapping in states:
                self.assertEqual(len(mapping), 128)  # Seven independent labels.
            if mask in (1, 2, 4, 8):
                self.assertEqual(counts, [24, 8, 4, 0])

    def test_literal_phase_constraint_and_identity_helper_stability(self):
        # Identity on X and Z fixes Y. A projective-only mapper would
        # incorrectly accept the negative-Y image below.
        mapping = _extend_clifford_map({0: (0, 0)}, 1, 1, 0, 1)
        mapping = _extend_clifford_map(mapping, 2, 2, 0, 1)
        self.assertIsNotNone(_extend_clifford_map(mapping, 3, 3, 0, 1))
        self.assertIsNone(_extend_clifford_map(mapping, 3, 3, 2, 1))
        for width in (5, 7):
            states, counts = _clifford_characters(1, width)
            self.assertEqual(states, [])
            self.assertEqual(counts, [24, 8, 4, 0])
        # For arbitrary identity helpers, exact flat support matching
        # forces every supported input label to a carrier-only label.
        # The same seven-dimensional phase contradiction therefore persists.

    def test_literal_full_intertwiner_on_every_column_and_actual_inverse(self):
        word = full_intertwiner()
        forward = native_columns(4, word)
        self.assertEqual(native_columns(4, _inverse(word)), adjoint(forward))
        for address in range(4):
            actual = native_columns(4, _inverse(word) + _fixed_query(address) + word)
            self.assertEqual(actual, native_columns(4, _wire_translation(address)))
        self.assertEqual(ledger(word), (6, 3, 65, 4))
        self.assertEqual(ledger(_inverse(word)), (6, 3, 63, 4))

    def test_fixed_pauli_axis_rank_and_disjoint_cell_failure(self):
        def rank(rows):
            pivots = {}
            for value in rows:
                while value:
                    bit = value.bit_length() - 1
                    if bit not in pivots:
                        pivots[bit] = value
                        break
                    value ^= pivots[bit]
            return len(pivots)

        def relabel(label, placement, width):
            x, z = label & 15, label >> 4
            return (sum(1 << placement[j] for j in range(4) if (x >> j) & 1)
                    | sum(1 << (width + placement[j])
                          for j in range(4) if (z >> j) & 1))

        def rows_for(axes, fixed, width):
            return [sum(_pairing(axis, value, width) << j
                        for j, value in enumerate(fixed)) for axis in axes]

        axes = [axis for band in AXES for axis in band]
        fixed = [16, 32, 4, 8]  # Z_u0, Z_u1, X_v0, X_v1.
        rows = rows_for(axes, fixed, 4)
        self.assertEqual(rows, [1, 1, 2, 2, 3, 3])
        self.assertEqual(rank(rows), 2)
        kernel = [value for value in range(16)
                  if all((row & value).bit_count() % 2 == 0 for row in rows)]
        self.assertEqual(kernel, [0, 4, 8, 12])  # Exactly span{X_v0,X_v1}.

        # Two natural disjoint cells on the p=3 field bank retain four
        # independent field-fixed controls. Identity helpers add 2h to
        # both fixed spaces, so they cannot repair the rank deficit.
        for helpers in (0, 1, 3):
            width = 8 + helpers
            doubled = [relabel(axis, placement, width)
                       for placement in ((0, 1, 4, 5), (2, 3, 6, 7))
                       for axis in axes]
            fixed = ([1 << (width + j) for j in range(4)]
                     + [1 << (4 + j) for j in range(4)]
                     + [1 << j for j in range(8, width)]
                     + [1 << (width + j) for j in range(8, width)])
            self.assertEqual(rank(rows_for(doubled, fixed, width)), 4)
            self.assertEqual(len(fixed) - 4, 4 + 2 * helpers)
            self.assertGreater(len(fixed) - 4, 2 + 2 * helpers)

        # An involutive wire permutation maps every Hermitian Pauli to a
        # commuting Hermitian Pauli. A nonidentical image therefore has
        # operator distance 2; this is the robust rank filter's exact gap.
        target_fixed = set(range(256))
        for address in range(4):
            placement = [j ^ address for j in range(4)]
            for label in range(256):
                image = relabel(label, placement, 4)
                self.assertEqual(_pairing(label, image, 4), 0)
                if image != label:
                    target_fixed.discard(label)
        self.assertEqual(target_fixed, {0, 15, 240, 255})


if __name__ == '__main__':
    unittest.main()
