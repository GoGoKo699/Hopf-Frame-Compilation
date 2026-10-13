"""Exact paid logical encodings and bounded algebra diagnostics.

Gate words are chronological ``(name, (wires, ...))`` tuples, with wire
zero the low integer bit. Encoders act on every input, without helpers;
zero remainders are promised only on the stated vacuum/one-hot code.
No arbitrary-angle core or unrestricted Hopf compiler is emitted here.
"""
from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from numbers import Integral

from .native_coarse_fixture import Gate, Word, inverse_word, toffoli_word


def _integer(value: int, name: str, minimum: int = 0) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral) or value < minimum:
        raise ValueError(f'{name} must be an integer at least {minimum}.')
    return int(value)


def _gate(name: str, wires: Sequence[int]) -> Gate:
    arities = {'X': 1, 'Z': 1, 'H': 1, 'S': 1, 'SDG': 1, 'T': 1,
               'TDG': 1, 'CX': 2, 'CZ': 2, 'CCX': 3}
    bits = tuple(_integer(q, 'wire') for q in wires)
    if name not in arities or len(bits) != arities[name] or len(set(bits)) != len(bits):
        raise ValueError('Unsupported gate, arity, or repeated wire.')
    return name, bits


@dataclass(frozen=True)
class EncodingLayout:
    active: tuple[int, ...]
    remainder: tuple[int, ...]


def zeta_encoder(m: int, offset: int = 0) -> tuple[Word, EncodingLayout]:
    """Compress vacuum and m one-hot labels into ceil(log2(m+1)) data bits.

    Wire A-1 represents nonempty binary subset A, 1 <= A <= m. The full
    Boolean permutation remains defined on unpromised inputs. There are
    exactly m-m.bit_length() Toffolis; no temporary wire is allocated.
    """
    m, offset = _integer(m, 'data width', 1), _integer(offset, 'offset')
    r = m.bit_length()
    gates: list[Gate] = []
    for bit in range(r):
        for a in range(1, m + 1):
            b = a | (1 << bit)
            if a != b and b <= m:
                gates.append(('CX', (offset + b - 1, offset + a - 1)))
    for a in sorted(range(1, m + 1), key=lambda x: (-x.bit_count(), x)):
        if a.bit_count() >= 2:
            low = a & -a
            gates.append(('CCX', (offset + low - 1, offset + (a ^ low) - 1,
                                  offset + a - 1)))
    layout = EncodingLayout(
        tuple(offset + (1 << bit) - 1 for bit in range(r)),
        tuple(offset + a - 1 for a in range(1, m + 1) if a & (a - 1)),
    )
    return tuple(gates), layout


def product_encoder(block_sizes: Sequence[int]) -> tuple[Word, EncodingLayout]:
    sizes = tuple(_integer(m, 'block width', 1) for m in block_sizes)
    if not sizes:
        raise ValueError('At least one block is required.')
    gates: list[Gate] = []
    active, remainder, offset = [], [], 0
    for m in sizes:
        word, layout = zeta_encoder(m, offset)
        gates.extend(word)
        active.extend(layout.active)
        remainder.extend(layout.remainder)
        offset += m
    return tuple(gates), EncodingLayout(tuple(active), tuple(remainder))


def apply_boolean(value: int, word: Iterable[Gate]) -> int:
    value = _integer(value, 'basis label')
    for name, wires in word:
        name, bits = _gate(name, wires)
        if name not in ('X', 'CX', 'CCX'):
            raise ValueError('Only Boolean permutation gates are accepted.')
        if all((value >> q) & 1 for q in bits[:-1]):
            value ^= 1 << bits[-1]
    return value


def extract_bits(value: int, positions: Sequence[int]) -> int:
    value = _integer(value, 'basis label')
    positions = tuple(_integer(q, 'wire') for q in positions)
    if len(set(positions)) != len(positions):
        raise ValueError('Positions must be distinct.')
    return sum(((value >> q) & 1) << j for j, q in enumerate(positions))


def expand_native(word: Iterable[Gate]) -> Word:
    """Use the retained literal Toffoli; expand X/Z/CZ, preserving phases."""
    out: list[Gate] = []
    for name, wires in word:
        name, bits = _gate(name, wires)
        if name == 'CCX':
            out.extend(toffoli_word(*bits))
        elif name == 'X':
            out.extend((('H', bits), ('S', bits), ('S', bits), ('H', bits)))
        elif name == 'Z':
            out.extend((('S', bits), ('S', bits)))
        elif name == 'CZ':
            out.extend((('H', (bits[1],)), ('CX', bits), ('H', (bits[1],))))
        else:
            out.append((name, bits))
    return tuple(out)


def _parity_phase(wires: Sequence[int], phase: str) -> Word:
    wires = tuple(_integer(q, 'wire') for q in wires)
    if not wires or len(set(wires)) != len(wires) or phase not in ('T', 'TDG', 'S', 'SDG'):
        raise ValueError('A distinct nonempty parity support and phase are required.')
    compute = tuple(('CX', (q, wires[-1])) for q in wires[:-1])
    return compute + ((phase, (wires[-1],)),) + inverse_word(compute)


def six_encoder_native(offset: int = 0) -> Word:
    """The complete six-bit zeta permutation in 60 gates and 15 T/Tdagger.

    The three quadratic clears are jointly phase-synthesized. This is an
    upper bound, not a T-optimality claim. No helper or input promise is
    used by the physical circuit.
    """
    offset = _integer(offset, 'offset')
    boolean, _ = zeta_encoder(6, offset)
    out = [g for g in boolean if g[0] == 'CX']
    a, b, c = offset, offset + 1, offset + 3
    u, v, w = offset + 2, offset + 4, offset + 5
    out.extend(('H', (q,)) for q in (u, v, w))
    out.extend(('S', (q,)) for q in (a, b, c))
    out.extend(('T', (q,)) for q in (u, v, w))
    for pair in ((a, b), (a, c), (b, c), (a, u), (b, u),
                 (a, v), (c, v), (b, w), (c, w)):
        out.extend(_parity_phase(pair, 'TDG'))
    for triple in ((a, b, u), (a, c, v), (b, c, w)):
        out.extend(_parity_phase(triple, 'T'))
    out.extend(('H', (q,)) for q in (u, v, w))
    return tuple(out)


def product_six_native(blocks: int) -> Word:
    blocks = _integer(blocks, 'number of blocks', 1)
    return tuple(g for j in range(blocks) for g in six_encoder_native(6 * j))


def encode_six(value: int) -> int:
    value = _integer(value, 'six-bit input')
    if value >= 64:
        raise ValueError('A six-bit input is required.')
    x1, x2, x3, x4, x5, x6 = [(value >> j) & 1 for j in range(6)]
    a, b, c = x1 ^ x3 ^ x5, x2 ^ x3 ^ x6, x4 ^ x5 ^ x6
    u, v, w = x3 ^ (a & b), x5 ^ (a & c), x6 ^ (b & c)
    return a | (b << 1) | (u << 2) | (c << 3) | (v << 4) | (w << 5)


def decode_six(value: int) -> int:
    value = _integer(value, 'six-bit input')
    if value >= 64:
        raise ValueError('A six-bit input is required.')
    a, b, u, c, v, w = [(value >> j) & 1 for j in range(6)]
    x3, x5, x6 = u ^ (a & b), v ^ (a & c), w ^ (b & c)
    return sum(x << j for j, x in enumerate((a ^ x3 ^ x5, b ^ x3 ^ x6,
                                           x3, c ^ x5 ^ x6, x5, x6)))


def encoded_pauli_word(axis: str, bit: int, offset: int = 0) -> Word:
    """Exact E6 P E6-dagger in H/S/CNOT, including the literal Y phase."""
    bit, offset = _integer(bit, 'bit'), _integer(offset, 'offset')
    if bit >= 6 or axis not in ('X', 'Y', 'Z'):
        raise ValueError('An X/Y/Z axis and bit in 0,...,5 are required.')
    a, b, c, u, v, w = 0, 1, 3, 2, 4, 5
    d1, d2, d3, d4, d5, d6 = [int(bit == j) for j in range(6)]
    alpha, beta, gamma = d1 ^ d3 ^ d5, d2 ^ d3 ^ d6, d4 ^ d5 ^ d6
    xword = [('CX', (control, target)) for yes, control, target in
             ((alpha, b, u), (beta, a, u), (alpha, c, v),
              (gamma, a, v), (beta, c, w), (gamma, b, w)) if yes]
    xword += [('X', (target,)) for target, yes in
              ((u, d3 ^ (alpha & beta)), (v, d5 ^ (alpha & gamma)),
               (w, d6 ^ (beta & gamma)), (a, alpha), (b, beta), (c, gamma)) if yes]
    linear = ((a, u, v), (b, u, w), (u,), (c, v, w), (v,), (w,))
    quadratic = (((a, b), (a, c)), ((a, b), (b, c)), ((a, b),),
                 ((a, c), (b, c)), ((a, c),), ((b, c),))
    zword = [('Z', (q,)) for q in linear[bit]] + [('CZ', qs) for qs in quadratic[bit]]
    if axis == 'X':
        word = xword
    elif axis == 'Z':
        word = zword
    else:
        word = zword + xword + [('X', (0,)), ('S', (0,)), ('X', (0,)), ('S', (0,))]
    return expand_native((g, tuple(offset + q for q in qs)) for g, qs in word)


def binary_rank(rows: Iterable[int]) -> int:
    pivots: dict[int, int] = {}
    for raw in rows:
        value = _integer(raw, 'binary row')
        while value:
            p = value.bit_length() - 1
            if p in pivots:
                value ^= pivots[p]
            else:
                pivots[p] = value
                break
    return len(pivots)


def edge_pauli_witnesses(n: int, edges: Iterable[tuple[int, int]]) -> tuple[tuple[int, int, int, int], ...]:
    """Independent edge displacements and an anticommuting Pauli clique.

    Returns (u,v,x,z). Binary elimination is exact. The coefficient of the
    Hermitian Pauli i**popcount(x&z) X**x Z**z in the oriented edge
    generator has magnitude 2/2**n. This is not a circuit lower-bound API.
    """
    n = _integer(n, 'logical width', 1)
    pivots, selected = {}, []
    for raw_u, raw_v in edges:
        u, v = _integer(raw_u, 'label'), _integer(raw_v, 'label')
        if max(u, v) >= 1 << n or u == v:
            raise ValueError('Distinct labels inside the logical register are required.')
        value = u ^ v
        while value:
            p = value.bit_length() - 1
            if p in pivots:
                value ^= pivots[p]
            else:
                pivots[p] = value
                selected.append((u, v))
                break
    rows = [[u ^ v, 1 << j] for j, (u, v) in enumerate(selected)]
    positions, rank = [], 0
    for bit in range(n):
        chosen = next((j for j in range(rank, len(rows)) if (rows[j][0] >> bit) & 1), None)
        if chosen is None:
            continue
        rows[rank], rows[chosen] = rows[chosen], rows[rank]
        for j in range(len(rows)):
            if j != rank and (rows[j][0] >> bit) & 1:
                rows[j][0] ^= rows[rank][0]
                rows[j][1] ^= rows[rank][1]
        positions.append(bit)
        rank += 1
    z, out = 0, []
    for i, (u, v) in enumerate(selected):
        dual = sum(((rows[j][1] >> i) & 1) << p for j, p in enumerate(positions))
        z ^= dual
        out.append((u, v, u ^ v, z))
    return tuple(out)
