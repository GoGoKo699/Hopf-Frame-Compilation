"""Bounded emitted checks of the prefix-releasing masked-sum pipeline.

Affine-control chains are expanded into actual X/CX/CCX words for basis
checks. A six-wire native fixture uses private eight-phase CCX words.
Static event deadlines audit supports and release, not simulated timings;
the general depth and resource bounds remain analytic arguments.
"""
from __future__ import annotations

import random
import unittest

import numpy as np

try:
    from .test_counter_dirty_indicator import _add, _increment
    from .test_operator_source_compiler import _adjoint, _basis_action, _word_matrix
    from .test_readonly_dirty_increment import _read, _replace
except ImportError:
    from test_counter_dirty_indicator import _add, _increment
    from test_operator_source_compiler import _adjoint, _basis_action, _word_matrix
    from test_readonly_dirty_increment import _read, _replace


ATOL = 5e-11


def _form(*wires, constant=0):
    return constant, frozenset(wires)


def _xor(left, right):
    return left[0] ^ right[0], left[1] ^ right[1]


def _not(form):
    return form[0] ^ 1, form[1]


def _value(basis, form):
    return form[0] ^ (sum(basis >> wire & 1 for wire in form[1]) % 2)


def _expand(word):
    """Compile affine XOR/AND controls into literal classical elementary gates."""
    result = []
    for gate in word:
        if gate[0] == 'AX':
            _, form, target = gate
            assert target not in form[1]
            result += [('X', target)] * form[0]
            result += [('CX', wire, target) for wire in sorted(form[1])]
        elif gate[0] == 'AA':
            _, left, right, target = gate
            assert target not in left[1] | right[1]
            result += [('X', target)] * (left[0] & right[0])
            result += [('CX', wire, target) for wire in sorted(right[1])] * left[0]
            result += [('CX', wire, target) for wire in sorted(left[1])] * right[0]
            for a in sorted(left[1]):
                for b in sorted(right[1]):
                    result.append(('CX', a, target) if a == b else ('CCX', a, b, target))
        else:
            result.append(gate)
    return result


def _private_native(word, phase):
    """Localize affine CCX phases; raw affine supports are only CX controls."""
    result = []
    for gate in word:
        if gate[0] == 'CCX':
            gate = ('AA', _form(gate[1]), _form(gate[2]), gate[3])
        if gate[0] != 'AA':
            result += _expand([gate])
            continue
        _, left, right, target = gate
        assert phase not in left[1] | right[1] | {target}
        forms = [left, right, _form(target)]
        result.append(('H', target))
        for subset in range(8):
            parity = _form()
            for bit, form in enumerate(forms):
                if subset >> bit & 1:
                    parity = _xor(parity, form)
            compute = _expand([('AX', parity, phase)])
            result += compute + [('TDG' if subset.bit_count() % 2 else 'T', phase)]
            result += _adjoint(compute)
        result.append(('H', target))
    return result


def _and_oracle(literals, target, chain):
    assert len(chain) >= len(literals)
    if len(literals) == 1:
        return [('AX', literals[0], target)]
    root = [('AX', literals[0], chain[0])]
    rest = [('AA', _form(chain[index - 1]), literal, chain[index])
            for index, literal in enumerate(literals[1:], 1)]
    full, copy = root + rest, [('CX', chain[len(literals) - 1], target)]
    return full + copy + list(reversed(full)) + rest + copy + list(reversed(rest))


def _blocks(length):
    start = 0
    while start < length:
        size = min(max(1, start), length - start)
        yield start, size
        start += size


def _signed_block(block, predicate, dirty, helpers):
    complement = [('CX', dirty, wire) for wire in block]
    private_q = _increment(dirty, block, helpers)[0]
    return complement + _adjoint(private_q) + predicate + private_q + predicate + complement


def _pipeline(count):
    bits = count.bit_length()
    address = list(range(count))
    counters = [list(range(count + row * bits, count + (row + 1) * bits))
                for row in range(count)]
    accumulator = list(range(count + count * bits, count + (count + 1) * bits))
    next_free = accumulator[-1] + 1

    def reserve(size):
        nonlocal next_free
        wires = list(range(next_free, next_free + size))
        next_free += size
        return wires

    columns = [[counter[column] for counter in counters] for column in range(bits)]
    events, forest, masks, final_forms, survivors, heights, counts = [], [], [], {}, [], [], []
    producers = {}
    for column in range(bits):
        raw = list(columns[column])
        heights.append(len(raw))
        active, forms, node_count = raw, {wire: _form(wire) for wire in raw}, 0
        while len(active) > 2:
            following = []
            triples = len(active) // 3
            for offset in range(triples):
                a, b, c = active[3 * offset:3 * offset + 3]
                form_a = forms[a]
                parity = _xor(forms[a], forms[b])
                total = _xor(parity, forms[c])
                forms[b], forms[c] = parity, total
                forest += [('CX', a, b), ('CX', b, c)]
                following.append(c)
                node_count += 1
                upper = reserve(bits - column - 1)
                for index, wire in enumerate(upper):
                    columns[column + index + 1].append(wire)
                    masks.append((wire, column + index + 1))
                elapsed = 0
                for start, size in _blocks(len(upper)):
                    dirty = reserve(1)[0]
                    chain, q_helpers, phase = reserve(start + 2), reserve(size), reserve(1)[0]
                    prefix = [_not(_form(wire)) for wire in upper[:start]]
                    predicate = (_and_oracle([form_a, _not(parity)] + prefix, dirty, chain)
                                 + _and_oracle([parity, _not(total)] + prefix, dirty, chain))
                    block = upper[start:start + size]
                    word = _signed_block(block, predicate, dirty, q_helpers)
                    duration = 128 * (start + 1) + 64 * size
                    event = dict(column=column, offset=start, block=block, prefix=upper[:start],
                                 raw=raw, word=word, phase=phase,
                                 start=512 * column + elapsed,
                                 end=512 * column + elapsed + duration)
                    events.append(event)
                    for wire in block:
                        producers[wire] = event
                    elapsed += duration
            following += active[3 * triples:]
            active = following
        counts.append(node_count)
        survivors.append(active)
        final_forms.update(forms)
    for column, wires in enumerate(survivors):
        while len(wires) < 2:
            wire = reserve(1)[0]
            wires.append(wire)
            masks.append((wire, column))
    outputs = [[wires[index] for wires in survivors] for index in range(2)]
    events.sort(key=lambda event: (event['start'], event['column'], event['offset']))
    pipeline = [gate for event in events for gate in event['word']]
    u = _expand(pipeline) + forest
    loader = u + _add(outputs[0], accumulator)[0] + _add(outputs[1], accumulator)[0] + _adjoint(u)
    literals = []
    pattern = [index % 2 for index in range(count)]
    for control, counter, value in zip(address, counters, pattern):
        literals += _increment(control, counter, reserve(bits), negative=not value)[0]
    echo = _adjoint(loader) + literals + loader + _adjoint(literals)
    return dict(bits=bits, width=next_free, address=address, counters=counters,
                accumulator=accumulator, events=events, producers=producers,
                pipeline=_expand(pipeline), forest=forest, forms=final_forms,
                masks=masks, outputs=outputs, loader=loader, echo=echo,
                heights=heights, counts=counts, pattern=pattern)


class PipelinedDirtySumTests(unittest.TestCase):
    def test_deferred_parity_forest_and_complete_masked_echo(self):
        generator = random.Random(0xB10C)
        for count in (3, 4, 5):
            fixture = _pipeline(count)
            modulus = 1 << fixture['bits']
            backgrounds = [0, (1 << fixture['width']) - 1]
            backgrounds += [generator.getrandbits(fixture['width']) for _ in range(3)]
            for basis in backgrounds:
                raw = _basis_action(basis, fixture['pipeline'])
                retired = _basis_action(raw, fixture['forest'])
                for wire, form in fixture['forms'].items():
                    self.assertEqual(retired >> wire & 1, _value(raw, form))
                total = sum(_read(basis, counter) for counter in fixture['counters'])
                offset = sum((basis >> wire & 1) << column for wire, column in fixture['masks'])
                self.assertEqual(sum(_read(retired, output) for output in fixture['outputs']) % modulus,
                                 (total + offset) % modulus)
                loaded = _replace(basis, fixture['accumulator'],
                                  (_read(basis, fixture['accumulator']) + total + offset) % modulus)
                self.assertEqual(_basis_action(basis, fixture['loader']), loaded)
                self.assertEqual(_basis_action(loaded, _adjoint(fixture['loader'])), basis)
                for address in range(1 << count):
                    initial = _replace(basis, fixture['address'], address)
                    added = sum((address >> index & 1) == value
                                for index, value in enumerate(fixture['pattern']))
                    expected = _replace(initial, fixture['accumulator'],
                                        (_read(initial, fixture['accumulator']) + added) % modulus)
                    self.assertEqual(_basis_action(initial, fixture['echo']), expected)
                    self.assertEqual(_basis_action(expected, _adjoint(fixture['echo'])), initial)
            # The first low-column triple (1, 1, 0) has majority one.
            # Retiring it first changes its raw wires to (1, 0, 0), so a
            # subsequent evaluation of the recorded affine forms loses carry.
            witness = (1 << fixture['counters'][0][0]) | (1 << fixture['counters'][1][0])
            premature = _basis_action(witness, fixture['forest'] + fixture['pipeline'])
            self.assertNotEqual(
                sum(_read(premature, output) for output in fixture['outputs']) % modulus, 2)

    def test_emitted_blocked_increment_uses_updated_zero_prefix(self):
        # The first two singleton blocks fit in six physical wires: literal,
        # two target bits, two dirty masks, and one private phase bus.
        first_e = [('CX', 0, 3)]
        second_e = [('AA', _form(0), _not(_form(1)), 4)]

        def singleton(bit, dirty, predicate):
            q = [('CX', dirty, bit)]
            return q + q + predicate + q + predicate + q

        first, second = singleton(1, 3, first_e), singleton(2, 4, second_e)
        word = first + second
        classical = _expand(word)
        images = [_replace(basis, [1, 2], (_read(basis, [1, 2]) + (basis & 1)) % 4)
                  for basis in range(64)]
        for basis, expected in enumerate(images):
            self.assertEqual(_basis_action(basis, classical), expected)
            self.assertEqual(_basis_action(expected, _adjoint(classical)), basis)
        native = _private_native(word, 5)
        expected = np.eye(64)[:, images]
        np.testing.assert_allclose(_word_matrix(6, native), expected, atol=ATOL, rtol=0)
        np.testing.assert_allclose(_word_matrix(6, _adjoint(native)), expected.conj().T,
                                   atol=ATOL, rtol=0)
        # Once the low bit settles it occurs only as a Clifford control.
        for gate in _private_native(second, 5):
            if 1 in gate[1:]:
                self.assertEqual(gate[0:2], ('CX', 1))
        wrong_ones = first + singleton(2, 4, [('AA', _form(0), _form(1), 4)])
        self.assertEqual(_read(_basis_action(1, classical), [1, 2]), 1)
        self.assertEqual(_read(_basis_action(1, _expand(wrong_ones)), [1, 2]), 3)
        self.assertEqual(_read(_basis_action(1, _expand(second + first)), [1, 2]), 3)
        missing_phase = list(native)
        missing_phase.pop(next(index for index, gate in enumerate(missing_phase)
                               if gate[0] in ('T', 'TDG')))
        self.assertGreater(np.linalg.norm(_word_matrix(6, missing_phase) - expected), 0.5)

        # Deferred forests produce affine controls with overlapping supports.
        affine = [('AA', _form(0, 1), _form(1, 2, constant=1), 3)]
        affine_native = _private_native(affine, 4)
        affine_images = [_basis_action(basis, _expand(affine)) for basis in range(32)]
        affine_expected = np.eye(32)[:, affine_images]
        np.testing.assert_allclose(_word_matrix(5, affine_native), affine_expected,
                                   atol=ATOL, rtol=0)
        np.testing.assert_allclose(_word_matrix(5, _adjoint(affine_native)),
                                   affine_expected.conj().T, atol=ATOL, rtol=0)
        for gate in affine_native:
            for raw in (0, 1, 2):
                if raw in gate[1:]:
                    self.assertEqual(gate[0:2], ('CX', raw))

    def test_static_release_deadlines_supports_and_resource_recurrence(self):
        for length in range(1, 66):
            elapsed = 0
            for start, size in _blocks(length):
                elapsed += 128 * (start + 1) + 64 * size
                self.assertLessEqual(elapsed, 512 * (start + 1))
            self.assertLessEqual(elapsed, 512 * length)
        for count in (3, 4, 7, 15):
            fixture = _pipeline(count)
            bits = fixture['bits']
            for column, (height, compressed) in enumerate(zip(fixture['heights'], fixture['counts'])):
                self.assertEqual(height, count + sum(fixture['counts'][:column]))
                self.assertEqual(compressed, (height - 1) // 2)
                self.assertLessEqual(height * 2 ** column, count * 3 ** column)
            for event in fixture['events']:
                raw_controls = set(event['raw'] + event['prefix'])
                classical = _expand(event['word'])
                self.assertTrue(raw_controls.isdisjoint(gate[-1] for gate in classical))
                for wire in raw_controls:
                    if wire in fixture['producers']:
                        self.assertLessEqual(fixture['producers'][wire]['end'], event['start'])
                for offset, wire in enumerate(event['block'], event['offset']):
                    self.assertLessEqual(event['end'], 512 * (event['column'] + offset + 1))
                localized_depth = 8 * sum(gate[0] in ('AA', 'CCX') for gate in event['word'])
                self.assertLessEqual(localized_depth, event['end'] - event['start'])
                native = _private_native(event['word'], event['phase'])
                emitted_depth = sum(gate[0] in ('T', 'TDG') for gate in native)
                self.assertEqual(emitted_depth, localized_depth)
                self.assertLessEqual(event['start'] + emitted_depth,
                                     512 * (event['column'] + event['offset'] + 1))
                for gate in native:
                    if gate[0] in ('T', 'TDG'):
                        self.assertEqual(gate[1], event['phase'])
                    if raw_controls.intersection(gate[1:]):
                        self.assertEqual(gate[0], 'CX')
                        self.assertIn(gate[1], raw_controls)
                        self.assertNotIn(gate[-1], raw_controls)
                self.assertLessEqual(event['end'], 512 * bits)
            weighted = sum(used * (bits - column)
                           for column, used in enumerate(fixture['counts']))
            self.assertLessEqual(weighted * 2 ** bits, 3 * count * 3 ** bits)


if __name__ == '__main__':
    unittest.main()
