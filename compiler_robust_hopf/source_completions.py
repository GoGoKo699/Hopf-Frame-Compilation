"""Finite joint-source matrix diagnostics, not a native source compiler.

Regularization avoids rank selection in the mathematical construction.
Floating tolerances below validate diagnostic inputs only. The native cost
of the fine source and its inverse is deliberately not asserted here.
"""
from __future__ import annotations

from dataclasses import dataclass
from numbers import Integral
import math
import numpy as np


def postorder(n: int) -> tuple[int, ...]:
    if isinstance(n, bool) or not isinstance(n, Integral) or n < 1:
        raise ValueError("n must be a positive integer.")
    children = [[] for _ in range(1 << n)]
    for j in range(1, 1 << n):
        children[j - (j & -j)].append(j)
    out = []
    def visit(vertex):
        for child in children[vertex]:
            visit(child)
        out.append(vertex)
    visit(0)
    return tuple(out)


def cholesky_in_order(matrix: np.ndarray, order: tuple[int, ...]) -> np.ndarray:
    """Upper Cholesky in a specified vertex order, then restore coordinates."""
    matrix = np.asarray(matrix, complex)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1] or not len(matrix):
        raise ValueError("Expected a nonempty square matrix.")
    if (len(order) != len(matrix) or any(isinstance(i, bool) or not isinstance(i, Integral) for i in order)
            or sorted(order) != list(range(len(matrix)))):
        raise ValueError("order must be a permutation of all coordinates.")
    if not np.all(np.isfinite(matrix)) or not np.allclose(matrix, matrix.conj().T, atol=1e-12, rtol=0):
        raise ValueError("Expected a finite Hermitian matrix.")
    upper = np.linalg.cholesky(matrix[np.ix_(order, order)]).conj().T
    inverse_order = np.argsort(order)
    return upper[np.ix_(inverse_order, inverse_order)]


@dataclass(frozen=True)
class JointCompletion:
    residual: np.ndarray
    joint_gram: np.ndarray
    joint_cholesky: np.ndarray
    deficit_f: np.ndarray
    deficit_g: np.ndarray
    source_f: np.ndarray
    source_g: np.ndarray
    unitary_f: np.ndarray
    unitary_g: np.ndarray


def regularized_completion(f: np.ndarray, g: np.ndarray, *, atol: float = 1e-10) -> JointCompletion:
    """Complete factors fg†=E/2 and ff†,gg† <= 5I/8 on exactly 4N modes.

    Callers obtain f,g from the canonical shared-mode factors divided by
    sqrt(2). Cut-rank conclusions additionally require those canonical
    tree-correlated factors, not arbitrary inputs accepted here.
    """
    f, g = np.asarray(f, complex), np.asarray(g, complex)
    if not math.isfinite(atol) or atol < 0:
        raise ValueError("atol must be finite and nonnegative.")
    if f.ndim != 2 or f.shape != g.shape or len(f) < 2 or len(f) & (len(f)-1):
        raise ValueError("Factors must have common NxM shape, with power-of-two N >= 2.")
    if not np.all(np.isfinite(f)) or not np.all(np.isfinite(g)):
        raise ValueError("Factors must be finite.")
    N = len(f)
    I = np.eye(N)
    E = 2 * f @ g.conj().T
    if np.linalg.norm(E.conj().T @ E - I, 2) > atol:
        raise ValueError("Twice the factor overlap must be unitary.")
    P, Q = f @ f.conj().T, g @ g.conj().T
    if max(np.linalg.eigvalsh(P)[-1], np.linalg.eigvalsh(Q)[-1]) > 5/8 + atol:
        raise ValueError("The required Gram contraction slack is absent.")
    Pp, Qp = P + I/8, Q + I/8
    K = np.block([[Pp, E/2], [E.conj().T/2, Qp]])
    order = postorder(N.bit_length()-1)
    doubled = tuple(k for i in order for k in (i, N+i))
    R = cholesky_in_order(K, doubled)
    CF = cholesky_in_order(I-Pp, order)
    CG = cholesky_in_order(I-Qp, order)
    zero = np.zeros((N, N), complex)
    SF = np.vstack([R[:, :N], CF, zero])
    SG = np.vstack([R[:, N:], zero, CG])
    J11 = np.vstack([zero, zero, zero, I])
    J10 = np.vstack([zero, zero, I, zero])
    HF = np.eye(4*N) - (SF-J11) @ (SF-J11).conj().T
    HG = np.eye(4*N) - (SG-J10) @ (SG-J10).conj().T
    AF = HF[:, np.arange(4*N) ^ (3*N)]
    BG = HG[:, np.arange(4*N) ^ (2*N)]
    return JointCompletion(E, K, R, CF, CG, SF, SG, AF, BG)


def normalization_two(unitary: np.ndarray, accepted_dimension: int) -> np.ndarray:
    """The actual matrix word -URU†RU; no projection/discarding of outputs."""
    U = np.asarray(unitary, complex)
    if U.ndim != 2 or U.shape[0] != U.shape[1] or not np.all(np.isfinite(U)):
        raise ValueError("Expected a finite square matrix.")
    if (isinstance(accepted_dimension, bool) or not isinstance(accepted_dimension, Integral)
            or not 1 <= accepted_dimension <= len(U)):
        raise ValueError("Invalid accepted dimension.")
    if not np.allclose(U.conj().T @ U, np.eye(len(U)), atol=1e-10, rtol=0):
        raise ValueError("Expected an actual unitary to diagnostic tolerance.")
    signs = np.ones(len(U))
    signs[:accepted_dimension] = -1
    return -((U * signs) @ U.conj().T * signs) @ U


def square_root_witness(u: float) -> np.ndarray:
    """Four-mode deficit whose principal root does not preserve cut rank one."""
    if not math.isfinite(u) or not 0 <= u <= 1/8:
        raise ValueError("Use a finite witness parameter in [0,1/8].")
    J = np.ones((4, 4))/4
    M = np.array([[1,1,0,0],[1,2,0,0],[0,0,2,1],[0,0,1,1]], float)/4
    return (np.eye(4) - 2*u/(1+u*u)*J - 2*u*u/(1+u*u)*M)/2
