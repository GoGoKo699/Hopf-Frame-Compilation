"""Bounded emitted-circuit checks of an indicator built from dirty counters.

Every counter, ladder helper, and cyclic selector starts arbitrarily.
Classical words audit complete work return; small native matrices audit
literal phases and actual inverses. Native schedules retain parallel tree
rounds, leaf increments, and Fredkin matchings. No large Hilbert-space
matrix or general compiler interface is constructed here.
"""
from __future__ import annotations

import unittest

import numpy as np

try:
    from .test_operator_source_compiler import _adjoint, _basis_action, _word_matrix
    from .test_t_depth import _fredkin_batch, _inverse, _native_schedule, _word
except ImportError:
    from test_operator_source_compiler import _adjoint, _basis_action, _word_matrix
    from test_t_depth import _fredkin_batch, _inverse, _native_schedule, _word


ATOL = 5e-11


def _primitive(word):
    return word, _native_schedule(word)


def _serial(circuits):
    return ([gate for word, _ in circuits for gate in word],
            [stage for _, schedule in circuits for stage in schedule])


def _reverse(circuit):
    return _adjoint(circuit[0]), _inverse(circuit[1])


def _parallel(circuits):
    """Merge independent schedules by T layer, retaining every Clifford."""
    split = []
    for _, schedule in circuits:
        cliffords, phases = [[]], []
        for kind, gates in schedule:
            if kind == 'C':
                cliffords[-1] += gates
            else:
                assert kind == 'T'
                phases.append(gates)
                cliffords.append([])
        split.append((cliffords, phases))
    depth = max((len(phases) for _, phases in split), default=0)
    schedule = []
    for layer in range(depth + 1):
        gates = [gate for cliffords, _ in split if layer < len(cliffords)
                 for gate in cliffords[layer]]
        if gates:
            schedule.append(('C', gates))
        if layer < depth:
            schedule.append(('T', [gate for _, phases in split if layer < len(phases)
                                   for gate in phases[layer]]))
    return [gate for word, _ in circuits for gate in word], schedule


def _mcx(controls, target, helpers, *, negative_first=False):
    """Read-only controls; the first control enters only through a CNOT."""
    assert len(set(controls + [target] + helpers)) == len(controls) + 1 + len(helpers)
    count = len(controls)
    if not count:
        return _primitive([('X', target)])
    if count == 1:
        return _primitive([('CX', controls[0], target)]
                          + ([('X', target)] if negative_first else []))
    assert len(helpers) >= count
    root = [('CX', controls[0], helpers[0])]
    if negative_first:
        root.append(('X', helpers[0]))
    rest = [('CCX', helpers[j - 1], controls[j], helpers[j]) for j in range(1, count)]
    compute = root + rest
    copy = [('CX', helpers[count - 1], target)]
    first = compute + copy + _adjoint(compute)
    root_disabled = rest + copy + _adjoint(rest)
    return _primitive(first + _adjoint(root_disabled))


def _increment(control, counter, helpers, *, negative=False):
    parts = [_mcx([control] + counter[:j], counter[j], helpers, negative_first=negative)
             for j in reversed(range(1, len(counter)))]
    parts.append(_mcx([control], counter[0], helpers, negative_first=negative))
    return _serial(parts)


def _add_by_increments(source, target, helpers):
    assert len(source) == len(target)
    return _serial([_increment(wire, target[j:], helpers) for j, wire in enumerate(source)])


def _add(source, target):
    """Reduced TTK six-step adder, with the carry-output gates omitted.

    Source bits may change inside the word, but return exactly; there is
    no carry output or workspace. The source is always a private counter.
    """
    assert len(source) == len(target) and source
    assert len(set(source + target)) == len(source) + len(target)
    bits = len(source)
    word = [('CX', source[i], target[i]) for i in range(1, bits)]
    word += [('CX', source[i], source[i + 1]) for i in reversed(range(1, bits - 1))]
    word += [('CCX', target[i], source[i], source[i + 1]) for i in range(bits - 1)]
    for i in reversed(range(1, bits)):
        word += [('CX', source[i], target[i]),
                 ('CCX', target[i - 1], source[i - 1], source[i])]
    word += [('CX', source[i], source[i + 1]) for i in range(1, bits - 1)]
    word += [('CX', source[i], target[i]) for i in range(bits)]
    return _primitive(word)


def _sum_tree(counters):
    current, rounds = list(range(len(counters))), []
    while len(current) > 1:
        pairs = list(zip(current[::2], current[1::2]))
        rounds.append(_parallel([_add(counters[left], counters[right])
                                 for left, right in pairs]))
        current = [right for _, right in pairs] + (current[-1:] if len(current) % 2 else [])
    return _serial(rounds), current[0], len(rounds)


def _shift_matchings(size, step):
    return [[(wire, (-offset - wire) % size) for wire in range(size)
             if wire < (-offset - wire) % size] for offset in (0, step)]


def _cyclic_shift(counter, selector):
    assert len(selector) == 1 << len(counter)
    parts = []
    for bit, control in enumerate(counter):
        for matching in _shift_matchings(len(selector), 1 << bit):
            pairs = [(selector[left], selector[right]) for left, right in matching]
            if not pairs:
                continue
            word = [gate for left, right in pairs
                    for gate in [('CX', left, right), ('CCX', control, right, left),
                                 ('CX', left, right)]]
            parts.append((word, _fredkin_batch(control, pairs)))
    return _serial(parts)


def _counter_fixture(bits, pattern, *, target=None, first_work=None):
    assert bits >= 1 and 0 <= pattern < 1 << bits
    address = list(range(bits))
    target = bits if target is None else target
    start = bits + 1 if first_work is None else first_work
    counter_bits = bits.bit_length()  # ceil(log2(bits+1)).
    counters = [list(range(start + i * counter_bits, start + (i + 1) * counter_bits))
                for i in range(bits + 1)]
    start += (bits + 1) * counter_bits
    helpers = [list(range(start + i * counter_bits, start + (i + 1) * counter_bits))
               for i in range(bits + 1)]
    start += (bits + 1) * counter_bits
    selector = list(range(start, start + (1 << counter_bits)))
    leaves, counter = counters[:-1], counters[-1]
    tree, root, rounds = _sum_tree(leaves)
    loader = _serial([tree, _add(leaves[root], counter), _reverse(tree)])
    increment = _parallel([_increment(wire, leaves[i], helpers[i], negative=not (pattern >> i & 1))
                           for i, wire in enumerate(address)])
    add_literals = _serial([_reverse(loader), increment, loader, _reverse(increment)])
    shift = _cyclic_shift(counter, selector)
    routed = _serial([add_literals, shift, _reverse(add_literals), _reverse(shift)])
    indicator = _serial([routed, _primitive([('X', selector[0])]), _reverse(routed)])
    copy = _primitive([('CX', selector[bits], target)])
    conjunction = _serial([copy, indicator, copy, _reverse(indicator)])
    return dict(width=selector[-1] + 1, address=address, target=target,
                counter_bits=counter_bits, counters=counters, helpers=helpers,
                selector=selector, rounds=rounds, loader=loader, increment=increment,
                add_literals=add_literals, shift=shift, routed=routed,
                indicator=indicator, conjunction=conjunction)


def _read(basis, wires):
    return sum((basis >> wire & 1) << bit for bit, wire in enumerate(wires))


def _replace(basis, wires, value):
    for bit, wire in enumerate(wires):
        basis = (basis & ~(1 << wire)) | ((value >> bit & 1) << wire)
    return basis


def _rotate(basis, wires, amount):
    original = _read(basis, wires)
    rotated = sum((original >> ((position + amount) % len(wires)) & 1) << position
                  for position in range(len(wires)))
    return _replace(basis, wires, rotated)


def _resources(circuit):
    schedule = circuit[1]
    return (sum(len(gates) for kind, gates in schedule if kind == 'T'),
            sum(kind == 'T' for kind, _ in schedule))


class CounterDirtyIndicatorTests(unittest.TestCase):
    def assert_native(self, width, circuit, images):
        # Unlike XOR queries, modular addition and rotation need not be
        # involutions: column j is the image of j, not its inverse image.
        expected = np.eye(1 << width)[:, images]
        np.testing.assert_allclose(_word_matrix(width, _word(circuit[1])), expected,
                                   atol=ATOL, rtol=0)
        np.testing.assert_allclose(_word_matrix(width, _word(_reverse(circuit)[1])),
                                   expected.conj().T, atol=ATOL, rtol=0)

    def test_read_only_mcx_increment_and_add_native_all_dirty_inputs(self):
        for controls in (2, 3):
            target, helpers = controls, list(range(controls + 1, 2 * controls + 1))
            circuit = _mcx(list(range(controls)), target, helpers, negative_first=True)
            images = [basis ^ (1 << target)
                      if not (basis & 1) and all(basis >> wire & 1 for wire in range(1, controls))
                      else basis for basis in range(1 << (2 * controls + 1))]
            self.assert_native(2 * controls + 1, circuit, images)
            self.assertEqual(_resources(circuit), (28 * (controls - 1), 16 * (controls - 1)))
            self.assertTrue(all(gate[-1] not in range(controls) for gate in circuit[0]))
            self.assertTrue(all(0 not in gate[1:] for gate in circuit[0] if gate[0] == 'CCX'))
        counter = [1, 2, 3]
        increment = _increment(0, counter, [4, 5, 6], negative=True)
        images = [_replace(basis, counter, (_read(basis, counter) + 1 - (basis & 1)) % 8)
                  for basis in range(128)]
        self.assert_native(7, increment, images)
        addition = _add_by_increments([0, 1], [2, 3], [4, 5])
        images = [_replace(basis, [2, 3], (_read(basis, [0, 1]) + _read(basis, [2, 3])) % 4)
                  for basis in range(64)]
        self.assert_native(6, addition, images)

    def test_linear_modular_adder_emitted_action_native_phases_and_cost(self):
        for bits in range(1, 6):
            source, target = list(range(bits)), list(range(bits, 2 * bits))
            addition = _add(source, target)
            images = [_replace(basis, target,
                               (_read(basis, source) + _read(basis, target)) % (1 << bits))
                      for basis in range(1 << (2 * bits))]
            for basis, desired in enumerate(images):
                self.assertEqual(_basis_action(basis, addition[0]), desired)
                self.assertEqual(_basis_action(desired, _reverse(addition)[0]), basis)
            if bits <= 4:
                self.assert_native(2 * bits, addition, images)
            self.assertEqual(sum(gate[0] == 'CCX' for gate in addition[0]), 2 * bits - 2)
            self.assertEqual(sum(gate[0] == 'CX' for gate in addition[0]),
                             1 if bits == 1 else 5 * bits - 6)
            self.assertEqual(_resources(addition), (14 * (bits - 1), 8 * (bits - 1)))

    def test_two_matching_cyclic_shift_native_phase_and_orientation(self):
        for bits in (1, 2):
            counter = list(range(bits))
            selector = list(range(bits, bits + (1 << bits)))
            shift = _cyclic_shift(counter, selector)
            width = bits + len(selector)
            images = [_rotate(basis, selector, _read(basis, counter)) for basis in range(1 << width)]
            self.assert_native(width, shift, images)
            count = depth = 0
            for bit in range(bits):
                for pairs in _shift_matchings(len(selector), 1 << bit):
                    self.assertEqual(len({wire for pair in pairs for wire in pair}), 2 * len(pairs))
                    if pairs:
                        count += 6 * len(pairs) + len(pairs) % 2
                        depth += 4
            self.assertEqual(_resources(shift), (count, depth))
            self.assertLessEqual(depth, 8 * bits)
        # Direct modular echo: every dirty offset and every admitted sum.
        for modulus in (2, 4, 8):
            for offset in range(modulus):
                for added in range(modulus):
                    self.assertEqual(((offset + added) % modulus - offset) % modulus, added)

    def test_emitted_counter_echo_and_conjunction_return_all_work(self):
        for bits in (2, 3):
            for pattern in range(1 << bits):
                fixture = _counter_fixture(bits, pattern)
                width, counter, selector = fixture['width'], fixture['counters'][-1], fixture['selector']
                self.assertEqual(width - bits - 1,
                                 2 * (bits + 1) * fixture['counter_bits'] + len(selector))
                free_mask = ((1 << width) - 1) ^ ((1 << bits) - 1)
                backgrounds = (0, free_mask,
                               sum(1 << wire for wire in range(bits, width, 2)),
                               sum(1 << wire for wire in range(bits + 1, width, 3)))
                for address in range(1 << bits):
                    added = bits - (address ^ pattern).bit_count()
                    for background in backgrounds:
                        basis = background | address
                        added_basis = _replace(basis, counter, (_read(basis, counter) + added) % len(selector))
                        indicator_basis = basis ^ (1 << selector[added])
                        expected = basis ^ ((1 << fixture['target']) if address == pattern else 0)
                        for name, desired in (('add_literals', added_basis),
                                              ('routed', _rotate(basis, selector, added)),
                                              ('indicator', indicator_basis),
                                              ('conjunction', expected)):
                            circuit = fixture[name]
                            actual = _basis_action(basis, circuit[0])
                            self.assertEqual(actual, desired)
                            self.assertEqual(_basis_action(actual, _reverse(circuit)[0]), basis)

    def test_shared_address_row_batch_has_disjoint_native_t_targets(self):
        for bits in (2, 3):
            rows, counter_bits = 1 << bits, bits.bit_length()
            private_width = 2 * (bits + 1) * counter_bits + (1 << counter_bits)
            fixtures = [_counter_fixture(bits, row, target=bits + row,
                                         first_work=bits + rows + row * private_width)
                        for row in range(rows)]
            circuits = [fixture['conjunction'] for fixture in fixtures]
            batch = _parallel(circuits)
            total_width = bits + rows + rows * private_width
            self.assertEqual(fixtures[-1]['width'], total_width)
            counts = [_resources(circuit) for circuit in circuits]
            self.assertEqual(_resources(batch), (sum(count for count, _ in counts), max(depth for _, depth in counts)))
            address_set = set(range(bits))
            for circuit in circuits:
                for gate in circuit[0]:
                    self.assertNotIn(gate[-1], address_set)
                    if address_set.intersection(gate[1:]):
                        self.assertEqual(gate[0], 'CX')
                        self.assertIn(gate[1], address_set)
            for kind, gates in batch[1]:
                self.assertTrue(all(0 <= wire < total_width for gate in gates for wire in gate[1:]))
                if kind == 'T':
                    targets = [gate[1] for gate in gates]
                    self.assertEqual(len(targets), len(set(targets)))
                    self.assertTrue(address_set.isdisjoint(targets))
            for address in range(rows):
                for background in (0, sum(1 << wire for wire in range(bits, total_width, 2))):
                    basis = background | address
                    expected = basis ^ (1 << (bits + address))
                    self.assertEqual(_basis_action(basis, batch[0]), expected)

    def test_emitted_tree_parallel_depth_and_complete_resource_ledger(self):
        for bits in (1, 2, 3):
            fixture = _counter_fixture(bits, (1 << bits) - 1)
            m, rounds = fixture['counter_bits'], fixture['rounds']
            inc_count, inc_depth = 14 * m * (m - 1), 8 * m * (m - 1)
            add_count, add_depth = 14 * (m - 1), 8 * (m - 1)
            self.assertEqual(rounds, (bits - 1).bit_length())
            self.assertEqual(_resources(fixture['increment']), (bits * inc_count, inc_depth))
            self.assertEqual(_resources(fixture['loader']), ((2 * bits - 1) * add_count,
                                                            (2 * rounds + 1) * add_depth))
            loader_count, loader_depth = _resources(fixture['loader'])
            self.assertEqual(_resources(fixture['add_literals']),
                             (2 * loader_count + 2 * bits * inc_count, 2 * loader_depth + 2 * inc_depth))
            count_a, depth_a = _resources(fixture['add_literals'])
            count_b, depth_b = _resources(fixture['shift'])
            self.assertEqual(_resources(fixture['conjunction']), (8 * (count_a + count_b),
                                                                 8 * (depth_a + depth_b)))
        # Missing the second copy leaks the initially arbitrary selector bit.
        fixture = _counter_fixture(2, 3)
        selector, target = fixture['selector'], fixture['target']
        copy = _primitive([('CX', selector[2], target)])
        missing = _serial([copy, fixture['indicator'], _reverse(fixture['indicator'])])
        basis = 1 << selector[2]  # Address zero is inactive for pattern three.
        self.assertEqual(_basis_action(basis, fixture['conjunction'][0]), basis)
        self.assertEqual(_basis_action(basis, missing[0]), basis ^ (1 << target))


if __name__ == '__main__':
    unittest.main()
