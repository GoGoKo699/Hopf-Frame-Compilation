"""Bounded permutation fixtures for SHORT_LABEL_SUBSET_ROUTING.md.

These exercise ideal XOR words, phases (all amplitudes are carried with
coefficient +1), guards, and actual inverses. They are not native query
synthesis or evidence for the asymptotic random-catalogue existence proof.
"""

from itertools import product
import unittest


def _matching(source, target, shifts):
    owner = {}

    def augment(x, seen):
        for shift in shifts:
            y = x ^ shift
            if y not in target or y in seen:
                continue
            seen.add(y)
            if y not in owner or augment(owner[y], seen):
                owner[y] = x
                return True
        return False

    for x in source:
        assert augment(x, set())
    return {x: y for y, x in owner.items()}


def _check_dense_shift_word_full_extension_guard_and_actual_inverse():
    source = {0, 1, 3, 7, 8}
    target = {2, 4, 5, 6, 10}
    # A complete tiny catalogue makes this fixture deterministic; the
    # polynomial-size existence argument is analytic, not inferred here.
    shifts = list(range(16))
    matching = _matching(sorted(source), target, shifts)
    catalogue = [0] + shifts + [0] * 15
    chosen = {x: 1 + (x ^ y) for x, y in matching.items()}
    back = {y: chosen[x] for x, y in matching.items()}

    def word(registers, reverse=False):
        active, x, label, reference = registers
        if not active:
            return registers
        if reverse:
            label ^= back.get(x, 0)
            x ^= catalogue[label]
            label ^= chosen.get(x, 0)
        else:
            label ^= chosen.get(x, 0)
            x ^= catalogue[label]
            label ^= back.get(x, 0)
        return active, x, label, reference

    occupied = list(product(range(2), range(16), range(32), range(3)))
    outputs = [word(item) for item in occupied]
    assert len(set(outputs)) == len(occupied)
    for original, output in zip(occupied, outputs):
        assert word(output, reverse=True) == original
        assert output[3] == original[3]
        if not original[0]:
            assert output == original
    for x in source:
        assert word((1, x, 0, 2)) == (1, matching[x], 0, 2)


def _check_staged_shear_compression_collision_and_reached_return():
    width, low_width = 7, 4
    low_mask = (1 << low_width) - 1
    source = {0, 16, 32, 64}
    # First hash isolates high words 1 and 2; high words 0 and 4 collide.
    # A second hash separates the remaining pair. The fixture genuinely
    # exercises retained stage information, rather than a single injection.
    rows = {0: (0, 0, 0, 0), 1: (1, 2, 0, 0), 2: (4, 0, 0, 0)}

    def shear(stage, high):
        return sum(((high & mask).bit_count() % 2) << bit
                   for bit, mask in enumerate(rows.get(stage, rows[0])))

    assigned = {16: 1, 32: 1, 0: 2, 64: 2}
    inverse_high = {}
    for x, stage in assigned.items():
        high, low = x >> low_width, x & low_mask
        key = stage, low ^ shear(stage, high)
        assert key not in inverse_high
        inverse_high[key] = high

    def word(registers, reverse=False):
        active, x, stage, reference = registers
        if not active:
            return registers
        high, low = x >> low_width, x & low_mask
        if reverse:
            high ^= inverse_high.get((stage, low), 0)
            low ^= shear(stage, high)
            stage ^= assigned.get((high << low_width) | low, 0)
        else:
            stage ^= assigned.get(x, 0)
            low ^= shear(stage, high)
            high ^= inverse_high.get((stage, low), 0)
        return active, (high << low_width) | low, stage, reference

    occupied = list(product(range(2), range(1 << width), range(4), range(3)))
    outputs = [word(item) for item in occupied]
    assert len(set(outputs)) == len(occupied)
    for original, output in zip(occupied, outputs):
        assert word(output, reverse=True) == original
        assert output[3] == original[3]
        if not original[0]:
            assert output == original
    codes = set()
    for x in range(1 << width):
        output = word((1, x, 0, 1))
        if x in source:
            assert output[1] >> low_width == 0
            assert output[2] in (1, 2)
            codes.add(output[1:3])
        else:
            assert output == (1, x, 0, 1)
    assert len(codes) == len(source)


class ShortLabelSubsetRoutingTests(unittest.TestCase):
    def test_dense_shift_word_full_extension_guard_and_actual_inverse(self):
        _check_dense_shift_word_full_extension_guard_and_actual_inverse()

    def test_staged_shear_compression_collision_and_reached_return(self):
        _check_staged_shear_compression_collision_and_reached_return()


if __name__ == "__main__":
    unittest.main()
