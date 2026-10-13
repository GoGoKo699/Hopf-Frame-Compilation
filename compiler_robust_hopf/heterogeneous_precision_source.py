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
    return Fraction(value)


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


def source_shift_residual_squared(beta: Fraction, length: int, shift: int) -> Fraction:
    """Exact ||S^r g_beta-beta^r g_beta||^2 for the normalized finite source."""
    beta = _rational(beta, "beta")
    length = _integer(length, "length", 1)
    shift = _integer(shift, "shift")
    if not 0 < beta < 1 or shift > length:
        raise ValueError("require 0 < beta < 1 and shift <= length")
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
    order = _integer(order, "order", 1)
    length = _integer(length, "length", 1)
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
