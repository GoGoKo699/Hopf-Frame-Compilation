"""Finite algebra for a small-ratio reusable source and Cayley kernels.

These exact arithmetic helpers describe coefficients, source-column errors,
and conservative parameter choices. They are not a native circuit emitter.
In particular, a scalar coefficient is not a specification of its rejected
operator; complete dilations and retained failure history are still required.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from numbers import Integral


def _integer(value: int, name: str, minimum: int = 0) -> int:
    if not isinstance(value, Integral) or isinstance(value, bool) or value < minimum:
        raise ValueError(f"{name} must be an integer at least {minimum}")
    return int(value)


def _rational(value: Fraction | int, name: str) -> Fraction:
    if isinstance(value, bool) or not isinstance(value, (Fraction, Integral)):
        raise ValueError(f"{name} must be an exact rational")
    # Fraction can retain NumPy Integral numerators. Normalize both parts so
    # later exact powers cannot silently become fixed-width arithmetic.
    if isinstance(value, Fraction):
        return Fraction(int(value.numerator), int(value.denominator))
    return Fraction(int(value))


def _ceil_log2(value: int) -> int:
    return (value - 1).bit_length()


def dyadic_precision(value: Fraction | int) -> int:
    """Return the reduced denominator exponent, rejecting non-dyadic values."""
    denominator = _rational(value, "value").denominator
    if denominator & (denominator - 1):
        raise ValueError("value must be dyadic")
    return denominator.bit_length() - 1


@dataclass(frozen=True)
class SmallRatioParameters:
    """A conservative ledger; no exponentially large precision state is built."""

    height: int
    order: int
    source_length: int
    rare_probability: Fraction
    rescaling_bits: int
    suffix_bits: int

    @property
    def minimum_exponent_log2(self) -> int:
        """log2(E_min) for E_min=2*2**height/height, without forming 2**height."""
        return self.height + 1 - (self.height.bit_length() - 1)

    @property
    def cap_exponent_at_least_three_precision_lengths(self) -> bool:
        """Certify M(E_min-B)>=3N by checking 2*height*B<=2**height."""
        return _ceil_log2(2 * self.height * self.rescaling_bits) <= self.height

    def coefficient_precision_bound(self, mantissa_bits: int) -> int:
        mantissa_bits = _integer(mantissa_bits, "mantissa_bits")
        return self.order * (mantissa_bits + self.rescaling_bits)


def small_ratio_parameters(height: int, *, suffix_constant: int) -> SmallRatioParameters:
    """Set P=n, M=2n, p=2^-ceil(log2(8n)), B=ceil(log2(8P/p)).

    A power-of-two n makes the P-term uniform selector exactly dyadic, as
    in the mixed-height family. Small n remain finite algebraic diagnostics;
    the cap-exponent property explicitly reports the large-n condition.
    The caller supplies the accounting multiplier C in K=C*ceil(log2(n+2));
    this helper does not certify that a particular C suffices for a circuit.
    """
    height = _integer(height, "height", 2)
    suffix_constant = _integer(suffix_constant, "suffix_constant", 1)
    if height & (height - 1):
        raise ValueError("height must be a power of two")
    probability = Fraction(1, 1 << _ceil_log2(8 * height))
    rescaling_bits = _ceil_log2(8 * height * probability.denominator)
    return SmallRatioParameters(height, height, 2 * height, probability,
                                rescaling_bits, suffix_constant * _ceil_log2(height + 2))


@dataclass(frozen=True)
class CommonSourceParameters:
    """Legacy power-of-four ledger, without constructing N=2**height.

    The exponent spacing is Delta=2**spacing_log2=N/n**2. The caller's
    maximum multiplier and suffix constant are accounting inputs; this
    object does not certify a native circuit's clean-work allocation.
    """

    height: int
    order: int
    source_length: int
    rare_probability: Fraction
    rescaling_bits: int
    spacing_log2: int
    maximum_multiplier: int
    suffix_bits: int

    @property
    def cap_exponent_at_least_three_precision_lengths(self) -> bool:
        """M(Delta-B)>=3N follows from Delta>=4B, with M=4n**2."""
        return _ceil_log2(4 * self.rescaling_bits) <= self.spacing_log2

    def coefficient_precision_overhead(self, mantissa_bits: int) -> int:
        """q_coef <= n*Delta + this polynomial overhead (integer arithmetic)."""
        mantissa_bits = _integer(mantissa_bits, "mantissa_bits")
        return self.order * (mantissa_bits + self.maximum_multiplier * self.rescaling_bits)

    def stage_precision_at_most_twice_endpoint(self, mantissa_bits: int) -> bool:
        """Certify n*q_coef<=2N using n**2*Delta=N, never forming N."""
        overhead = self.height * self.coefficient_precision_overhead(mantissa_bits)
        return _ceil_log2(overhead) <= self.height


def common_source_parameters(
    height: int, *, suffix_constant: int, maximum_multiplier: int | None = None,
) -> CommonSourceParameters:
    """Use P=n, M=4n**2 and Delta=N/n**2 for power-of-four heights.

    The default multiplier floor(n**2/log2(n)) corresponds to E<=N/log2(n).
    Other supplied bounds, including n*sqrt(n), are allowed. The suffix
    multiplier is supplied explicitly and is not certified sufficient here.
    """
    base = small_ratio_parameters(height, suffix_constant=suffix_constant)
    logarithm = base.height.bit_length() - 1
    if logarithm % 2 or base.height < 4:
        raise ValueError("height must be a power of four at least four")
    if maximum_multiplier is None:
        maximum_multiplier = base.height ** 2 // logarithm
    maximum_multiplier = _integer(maximum_multiplier, "maximum_multiplier", 1)
    return CommonSourceParameters(
        base.height, base.order, 4 * base.height ** 2, base.rare_probability,
        base.rescaling_bits, base.height - 2 * logarithm,
        maximum_multiplier, base.suffix_bits,
    )


@dataclass(frozen=True)
class ClippedSourceParameters:
    """All-height accounting for a clipped frame; not a native emitter.

    Only the scalar selector order P is rounded up to a power of two. The
    logical frame still has height n. Delta=2**spacing_log2=N/P**2 and
    H=N+clipping_exponent_offset remain symbolic. Below the stated source
    regime the theorem uses a finite-height fallback, not this source.
    """

    height: int
    order: int
    source_length: int
    rare_probability: Fraction
    rescaling_bits: int
    spacing_log2: int
    suffix_bits: int

    @property
    def clipping_exponent_offset(self) -> int:
        """The additive offset in H=N+ceil(log2(n))+6."""
        return self.order.bit_length() - 1 + 6

    @property
    def clipping_bound_relative_to_endpoint(self) -> Fraction:
        """4n*2**(-H) / 2**(-N), without forming either tiny number."""
        return Fraction(self.height, 16 * self.order)

    @property
    def minimum_expensive_exponent_log2(self) -> int:
        """log2(2N/P); this cutoff times P is exactly 2N."""
        return self.height + 1 - (self.order.bit_length() - 1)

    @property
    def cap_exponent_at_least_three_precision_lengths(self) -> bool:
        """M(Delta-B)>=3N follows from Delta>=4B and M*Delta=4N."""
        return _ceil_log2(4 * self.rescaling_bits) <= self.spacing_log2

    @property
    def uses_finite_fallback(self) -> bool:
        return not self.cap_exponent_at_least_three_precision_lengths

    @property
    def maximum_multiplier(self) -> int:
        """floor((H-1)/Delta)=P**2 in the qualified source regime."""
        if self.uses_finite_fallback:
            raise ValueError("the common source requires the large-height regime")
        # Delta>=4B>ceil(log2(n))+5, so the clipping offset adds no
        # further multiple of Delta. No exponentially large N is needed.
        return self.order ** 2

    def retained_exponent_within_cap(self, exponent: int) -> bool:
        """Check E<=H-1, not membership in a supplied open branch bin."""
        exponent = _integer(exponent, "exponent")
        beyond_offset = exponent - (self.clipping_exponent_offset - 1)
        return beyond_offset <= 0 or _ceil_log2(beyond_offset) <= self.height

    def coefficient_precision_overhead(self, mantissa_bits: int) -> int:
        """q_coef<=N/P+this overhead, using r<Delta and degree<=P."""
        mantissa_bits = _integer(mantissa_bits, "mantissa_bits")
        return self.order * (mantissa_bits + self.maximum_multiplier * self.rescaling_bits)

    def stage_precision_at_most_twice_endpoint(self, mantissa_bits: int) -> bool:
        """Certify n*q_coef<=2N, using n/P<=1 and exponent comparisons."""
        mantissa_bits = _integer(mantissa_bits, "mantissa_bits")
        if self.uses_finite_fallback:
            return False
        overhead = self.height * self.coefficient_precision_overhead(mantissa_bits)
        return _ceil_log2(overhead) <= self.height


def clipped_source_parameters(height: int, *, suffix_constant: int) -> ClippedSourceParameters:
    """Set P=2**ceil(log2(n)), M=4P**2, p=1/(8P), B=6+2log2(P).

    Every positive integer height is accepted. ``uses_finite_fallback``
    marks heights for which the common-source inequalities do not yet hold.
    The supplied C gives K=min(n,C*ceil(log2(P+2))); it is an accounting
    multiplier, not a certificate that C provides sufficient clean work.
    """
    height = _integer(height, "height", 1)
    suffix_constant = _integer(suffix_constant, "suffix_constant", 1)
    logarithm = _ceil_log2(height)
    order = 1 << logarithm
    return ClippedSourceParameters(
        height, order, 4 * order ** 2, Fraction(1, 8 * order),
        6 + 2 * logarithm, height - 2 * logarithm,
        min(height, suffix_constant * _ceil_log2(order + 2)),
    )


def clipped_branch_decision(approximation: Fraction, threshold: Fraction) -> bool:
    """Return the finite certificate's decision q<=3*zeta/2.

    The caller must certify |q-s|<=zeta/8 for the nonnegative smaller
    branch s. Then clipping implies s<=13*zeta/8, while retaining implies
    s>11*zeta/8>zeta. This comparison never decides equality of s with a
    dyadic threshold and does not choose an exact logarithm of s.
    """
    approximation = _rational(approximation, "approximation")
    threshold = _rational(threshold, "threshold")
    if threshold <= 0:
        raise ValueError("threshold must be positive")
    return approximation <= Fraction(3, 2) * threshold


def bottom_precision_trace(retained_exponents_by_depth, suffix_bits: int) -> int:
    """Sum depth maxima in the bottom K layers, with empty maxima zero.

    Depths are supplied root first and contain only retained nodes' supplied
    bin exponents. This accounts for disjoint planes within a layer; it does
    not sum over nodes, infer bins, or certify any native implementation.
    """
    depths = tuple(tuple(_integer(exponent, "exponent") for exponent in depth)
                   for depth in retained_exponents_by_depth)
    suffix_bits = _integer(suffix_bits, "suffix_bits")
    if suffix_bits > len(depths):
        raise ValueError("suffix_bits cannot exceed the number of depths")
    return sum(max(depth, default=0) for depth in depths[len(depths) - suffix_bits:])


def exponent_metadata(exponent: int, spacing: int) -> tuple[int, int]:
    """Exact E=a*Delta+r metadata, with 0<=r<Delta; no bin rounding."""
    exponent = _integer(exponent, "exponent")
    spacing = _integer(spacing, "spacing", 1)
    return divmod(exponent, spacing)


def common_source_ratio(spacing: int, rescaling_bits: int) -> Fraction:
    """Explicit beta=2**(B-Delta) for finite diagnostics, not huge ledgers."""
    spacing = _integer(spacing, "spacing", 1)
    rescaling_bits = _integer(rescaling_bits, "rescaling_bits")
    if spacing <= rescaling_bits:
        raise ValueError("spacing must exceed rescaling_bits")
    return Fraction(1, 1 << (spacing - rescaling_bits))


def source_shift_residual_squared(beta: Fraction, length: int, shift: int) -> Fraction:
    """Exact finite-source residual, with S^r=0 for every r>=length."""
    beta = _rational(beta, "beta")
    length = _integer(length, "length", 1)
    shift = _integer(shift, "shift")
    if not 0 < beta < 1:
        raise ValueError("require 0 < beta < 1")
    if shift >= length:
        return beta ** (2 * shift)
    return beta ** (2 * length) * (1 - beta ** (2 * shift)) / (1 - beta ** (2 * length))


def cayley_pair(parameter: Fraction) -> tuple[Fraction, Fraction]:
    """Real rotation coefficients for (I+tJ)(I-tJ)^-1, J=[[0,-1],[1,0]]."""
    parameter = _rational(parameter, "parameter")
    denominator = 1 + parameter * parameter
    return (1 - parameter * parameter) / denominator, 2 * parameter / denominator


def cayley_kernel_coefficients(
    mantissa: Fraction, order: int, rare_probability: Fraction, rescaling_bits: int,
) -> tuple[Fraction, ...]:
    """Return the exact f_r; zero mantissa represents an inactive row.

    Its accepted kernel is still (1-p)I on that row. A separate outer activity
    guard, rather than setting this baseline to I, supplies idle sectors.
    """
    mantissa = _rational(mantissa, "mantissa")
    order = _integer(order, "order", 1)
    rare_probability = _rational(rare_probability, "rare_probability")
    rescaling_bits = _integer(rescaling_bits, "rescaling_bits")
    dyadic_precision(mantissa)
    dyadic_precision(rare_probability)
    if abs(mantissa) > 1 or rare_probability.numerator != 1 or not 0 < rare_probability < 1:
        raise ValueError("require |mantissa| <= 1 and reciprocal-power-of-two probability")
    ratio = mantissa / (1 << rescaling_bits)
    prefactor = 2 * order * (1 - rare_probability) / rare_probability
    coefficients = tuple(prefactor * ratio ** r for r in range(1, order + 1))
    if any(abs(coefficient) > 1 for coefficient in coefficients):
        raise ValueError("rescaling is insufficient for contraction coefficients")
    return coefficients


def kernel_source_error_bound(
    beta: Fraction, mantissa: Fraction, order: int, length: int,
    rare_probability: Fraction, rescaling_bits: int,
) -> Fraction:
    """A rational cap-plus-tail bound against (1-p) R(t), t=beta*u/2**B.

    Replacing sqrt(1-beta**(2M)) by 1-beta**2 makes the cap bound rational
    and conservative. The bound concerns the complete source column.
    """
    beta = _rational(beta, "beta")
    mantissa = _rational(mantissa, "mantissa")
    rare_probability = _rational(rare_probability, "rare_probability")
    order = _integer(order, "order", 1)
    length = _integer(length, "length", 1)
    rescaling_bits = _integer(rescaling_bits, "rescaling_bits")
    if not 0 < beta < 1 or order >= length:
        raise ValueError("require 0 < beta < 1 and order < length")
    coefficients = cayley_kernel_coefficients(mantissa, order, rare_probability,
                                              rescaling_bits)
    probability = Fraction(rare_probability)
    parameter = beta * Fraction(mantissa) / (1 << rescaling_bits)
    cap = probability * sum(map(abs, coefficients), Fraction()) * beta ** length
    cap /= order * (1 - beta * beta)
    tail = 2 * (1 - probability) * abs(parameter) ** (order + 1) / (1 - abs(parameter))
    return cap + tail


def common_source_kernel_coefficients(
    mantissa: Fraction, order: int, rare_probability: Fraction, rescaling_bits: int,
    *, multiplier: int, remainder: int,
) -> tuple[Fraction, ...]:
    """Exact coefficients for S^(a*nu) J^nu, keeping r in dyadic metadata.

    For E=a*Delta+r, each coefficient equals
    2P(1-p)/p * u**nu * 2**(-(r+a*B)*nu). The quarterturn power is
    nu, independently of the source shift a*nu.
    """
    multiplier = _integer(multiplier, "multiplier", 1)
    remainder = _integer(remainder, "remainder")
    rescaling_bits = _integer(rescaling_bits, "rescaling_bits")
    return cayley_kernel_coefficients(mantissa, order, rare_probability,
                                      remainder + multiplier * rescaling_bits)


def common_source_kernel_error_bound(
    beta: Fraction, mantissa: Fraction, order: int, length: int,
    rare_probability: Fraction, rescaling_bits: int, *, multiplier: int, remainder: int,
) -> Fraction:
    """Rational source-column cap-plus-tail bound, including saturated shifts."""
    beta = _rational(beta, "beta")
    mantissa = _rational(mantissa, "mantissa")
    rare_probability = _rational(rare_probability, "rare_probability")
    order = _integer(order, "order", 1)
    length = _integer(length, "length", 1)
    multiplier = _integer(multiplier, "multiplier", 1)
    remainder = _integer(remainder, "remainder")
    rescaling_bits = _integer(rescaling_bits, "rescaling_bits")
    if not 0 < beta < 1:
        raise ValueError("require 0 < beta < 1")
    coefficients = common_source_kernel_coefficients(
        mantissa, order, rare_probability, rescaling_bits,
        multiplier=multiplier, remainder=remainder,
    )
    probability = Fraction(rare_probability)
    parameter = beta ** multiplier * Fraction(mantissa)
    parameter /= 1 << (remainder + multiplier * rescaling_bits)
    cap = Fraction()
    for degree, coefficient in enumerate(coefficients, 1):
        shift = multiplier * degree
        residual_bound = (beta ** shift if shift >= length
                          else beta ** length / (1 - beta * beta))
        cap += abs(coefficient) * residual_bound
    cap *= probability / order
    tail = 2 * (1 - probability) * abs(parameter) ** (order + 1) / (1 - abs(parameter))
    return cap + tail


def splitter_gauge_coefficients(target, coarse, carrier, node_gauge=1):
    """Local detail coefficients for the mixed chart's expensive node.

    The caller supplies normalized real pairs (c,s), (a,b), a nonzero complex
    carrier r, and node gauge d. For B=R(a,b) diag(conj(r),1) and L=(a,b),
    these are the literal coefficients C,A in the full-column gauge identity.
    This algebraic helper neither synthesizes a splitter nor infers its phase.
    Known corners use their separate identity chart and are excluded here.
    """
    cosine, sine = target
    coarse_cosine, coarse_sine = coarse
    denominator = node_gauge * cosine * sine
    if denominator == 0 or carrier == 0 or coarse_cosine * coarse_sine == 0:
        raise ValueError("the splitter chart requires nonzero gauges and branches")
    factor = carrier.conjugate() / denominator
    coefficient = factor * coarse_cosine * coarse_sine
    correction = factor * (cosine ** 2 * coarse_sine ** 2 - sine ** 2 * coarse_cosine ** 2)
    return coefficient, correction
