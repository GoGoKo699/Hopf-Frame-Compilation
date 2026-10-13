"""Phase-orbit equations and exact rejection certificate; no endpoint claim."""
from fractions import Fraction as F
from functools import lru_cache
import json
from pathlib import Path
import unittest

import numpy as np

from compiler_robust_hopf.programmed_reflections import (
    cycle_discrepancy, hermitian_coordinates, nearest_reversible,
    phase_reflection, programmed_pair, recover_phase_gauge,
    recover_shared_pair, reversing_moment_system, shared_pair_floor,
    shared_relative_kernel,
)
from compiler_robust_hopf.spectral_reflection_chart import (
    hopf_corner_word, reflection_mixer_words,
)
from compiler_robust_hopf.native_coarse_fixture import (
    gate_counts, inverse_word, simulate_columns,
)
from tests.test_conjugate_recovery import exact_native
from verification.fault_tolerant.exact_arithmetic import CQ2, Q2


def unitary(rng, size, real=False):
    a = rng.normal(size=(size,size)).astype(complex)
    if not real:
        a += 1j*rng.normal(size=(size,size))
    return np.linalg.qr(a)[0]


def small_real_rotation(size, theta=.04):
    w = np.eye(size)
    c,s = np.cos(theta),np.sin(theta)
    for target in (size//2,3*size//4):
        source = 0 if target == size//2 else size//2
        w[[source,target],:] = np.array([[c,-s],[s,c]])@w[[source,target],:]
    return w


def exact_word(word, qubits):
    return exact_native([(g.lower(),*qs) for g,qs in word],qubits)


@lru_cache(maxsize=2)
def exact_mixers():
    return tuple(exact_word(word,3) for word in reflection_mixer_words())


def rational(value):
    if isinstance(value,CQ2):
        value = value.real
    if isinstance(value,Q2):
        if value.b:
            raise AssertionError('The specified certificate must have rational coefficients.')
        return value.a
    return F(value)


class ProgrammedReflectionTests(unittest.TestCase):
    def test_phase_gauge_complete_directed_action(self):
        rng = np.random.default_rng(10501)
        for size in (4,8,16):
            b = unitary(rng,size)
            d = np.exp(1j*rng.uniform(-3,3,size))
            k = d[:,None]*b*d.conj()[None,:]
            got = recover_phase_gauge(b,k)
            self.assertIsNotNone(got)
            np.testing.assert_allclose(got[:,None]*b*got.conj()[None,:],k,atol=2e-12,rtol=0)

    def test_disconnected_gauge_and_opposite_edges(self):
        b = np.array([[0,1,0,0],[1,0,0,0],[0,0,0,1],[0,0,1,0]],complex)
        d = np.exp(1j*np.array([.2,.6,-.1,1.3]))
        good = d[:,None]*b*d.conj()[None,:]
        self.assertIsNotNone(recover_phase_gauge(b,good))
        bad = good.copy(); bad[1,0] *= 1j
        self.assertIsNone(recover_phase_gauge(b,bad))

    def test_cycle_failure_even_with_matching_magnitudes(self):
        rng = np.random.default_rng(10502)
        c = unitary(rng,4)
        b = (c*np.array([1,1,-1,-1]))@c.conj().T
        k = b.conj()
        self.assertGreater(cycle_discrepancy(b,k,(0,1,2)),1e-7)
        self.assertIsNone(recover_phase_gauge(b,k))
        for _ in range(4):
            d = np.exp(1j*rng.uniform(-3,3,4))
            error = np.linalg.norm(k-d[:,None]*b*d.conj()[None,:],2)
            self.assertGreaterEqual(error+1e-14,cycle_discrepancy(b,k,(0,1,2)))

    def test_single_pair_reversibility_and_optimum(self):
        rng = np.random.default_rng(10503)
        for size in (4,8,16):
            j = np.diag([1]*(size//2)+[-1]*(size//2))
            c = unitary(rng,size)
            v = programmed_pair(c,j,np.exp(1j*rng.normal(size=size)))
            np.testing.assert_allclose(j@v@j,v.conj().T,atol=2e-12,rtol=0)
            w = small_real_rotation(size)
            best,floor = nearest_reversible(w,j)
            self.assertAlmostEqual(np.linalg.norm(best-w,2),floor,places=12)
            h = best@j
            np.testing.assert_allclose(h@h,np.eye(size),atol=2e-12,rtol=0)
            self.assertAlmostEqual(np.trace(h).real,np.trace(j),places=11)

    def test_all_diagonal_centers_obey_path_floor(self):
        theta=.23; w=small_real_rotation(8,theta)
        for bits in range(256):
            j=np.diag([1 if bits>>k&1 else -1 for k in range(8)])
            self.assertGreaterEqual(np.linalg.norm(w-j@w.T@j,2)+1e-14,np.sin(theta)**2)

    def test_shared_kernel_reconstructs_generated_products(self):
        rng=np.random.default_rng(10504)
        for size in (4,8,16):
            c=unitary(rng,size);j=np.diag([1]*(size//2)+[-1]*(size//2))
            u=np.exp(1j*rng.normal(size=size));v=np.exp(1j*rng.normal(size=size))
            w=phase_reflection(c,j,v)@phase_reflection(c,j,u)
            relative=v.conj()*u
            np.testing.assert_allclose(shared_relative_kernel(c,j,w)@relative,0,atol=2e-12,rtol=0)
            recovered=recover_shared_pair(c,j,w)
            self.assertIsNotNone(recovered)
            ru,rv=recovered
            np.testing.assert_allclose(phase_reflection(c,j,rv)@phase_reflection(c,j,ru),w,atol=2e-10,rtol=0)

    def test_shared_kernel_floor_bounds_every_sample(self):
        rng=np.random.default_rng(10505)
        for size in (4,8):
            c=unitary(rng,size);j=np.diag([1]*(size//2)+[-1]*(size//2))
            w=unitary(rng,size,real=True)
            bound=shared_pair_floor(c,j,w)
            for _ in range(5):
                u=np.exp(1j*rng.normal(size=size));v=np.exp(1j*rng.normal(size=size))
                candidate=phase_reflection(c,j,v)@phase_reflection(c,j,u)
                self.assertGreaterEqual(np.linalg.norm(candidate-w,2)+1e-12,bound)

    def test_real_basis_gap_and_multiple_kernel_scope(self):
        h=np.array([[1,1],[1,-1]])/np.sqrt(2)
        sh=np.diag([1,1j])@h
        c1=sh@np.diag([1,np.exp(-1j*np.pi/4)])@sh.conj().T
        c=np.kron(c1,c1);j=np.kron(np.diag([1,-1]),np.diag([1,-1]))
        b=c.conj().T@j@c
        np.testing.assert_allclose(abs(b)**2,np.full((4,4),.25),atol=1e-14,rtol=0)
        w=small_real_rotation(4,.09)
        self.assertGreaterEqual(shared_pair_floor(c,j,w)+1e-14,np.linalg.norm(w-np.eye(4),'fro')**2/8)
        self.assertIsNone(recover_shared_pair(np.eye(4),j,np.eye(4)))

    def test_linear_system_does_not_assume_real_factors(self):
        c1,c2=(simulate_columns(w,nqubits=3) for w in reflection_mixer_words())
        target=simulate_columns(hopf_corner_word(),nqubits=3).real
        j=np.diag([(-1)**(i&1) for i in range(8)])
        matrix,rhs=reversing_moment_system(target,c1,c2,j)
        self.assertEqual(matrix.shape,(80,64));self.assertEqual(rhs.shape,(80,))
        self.assertEqual(np.linalg.matrix_rank(matrix,tol=1e-9),62)
        self.assertEqual(len(hermitian_coordinates(8)),64)

    def test_native_mixers_and_easy_target_counts(self):
        for word in reflection_mixer_words():
            self.assertEqual((len(word),gate_counts(word)['T_count']),(68,28))
        word=hopf_corner_word()
        self.assertEqual((len(word),gate_counts(word)['T_count']),(97,14))
        self.assertLessEqual({g for g,_ in word},{'H','S','SDG','T','TDG','CX'})
        full=exact_word(word,3);back=exact_word(inverse_word(word),3)
        for row in range(8):
            for col in range(8):
                expected=(1 if col==0 else -1) if row==(col-1)%8 else 0
                self.assertEqual(full[row][col],CQ2(expected))
                self.assertEqual(back[row][col],full[col][row].conj())

    def test_exact_dual_certificate_all_hermitian_coordinates(self):
        cert=json.loads((Path(__file__).parent/'fixtures/programmed_reflection_dual.json').read_text())
        q=[F(x) for x in cert['q']];lam=[F(x) for x in cert['lambda']]
        c1,c2=exact_mixers();size=8;zero=CQ2()
        def moments(c,h):
            active=[(a,b,h[a][b]) for a in range(size) for b in range(size) if h[a][b]!=zero]
            return [rational(sum((c[a][k].conj()*value*c[b][k] for a,b,value in active),zero)) for k in range(size)]
        center=[[CQ2((-1)**(i&1)) if i==j else zero for j in range(size)] for i in range(size)]
        rhs=[F(0)]*64+moments(c1,center)+moments(c2,center)
        numeric_columns=[]
        for a,b,kind in hermitian_coordinates(size):
            h=[[zero for _ in range(size)] for _ in range(size)]
            if kind=='diagonal':h[a][a]=CQ2(1)
            elif kind=='real':h[a][b]=h[b][a]=CQ2(1)
            else:h[a][b],h[b][a]=CQ2(0,1),CQ2(0,-1)
            wh=[[-x for x in h[i+1]] if i<7 else list(h[0]) for i in range(size)]
            defect=[[wh[i][j]-wh[j][i].conj() for j in range(size)] for i in range(size)]
            col=[rational(defect[i][j].real) for i in range(size) for j in range(i+1,size)]
            col += [rational(defect[i][j].imag) for i in range(size) for j in range(i,size)]
            col += moments(c1,h)+moments(c2,wh)
            functional=rational(sum((q[i]*h[i][0] for i in range(size)),zero))
            self.assertEqual(sum((x*y for x,y in zip(lam,col)),F(0)),functional)
            numeric_columns.append([float(x) for x in col])
        qnorm=sum(x*x for x in q)
        self.assertEqual(qnorm,F(4378,963))
        self.assertEqual(sum(x*y for x,y in zip(lam,rhs)),qnorm)
        weight=2*sum(abs(x) for x in lam[:64])+sum(abs(x) for x in lam[72:])
        self.assertEqual(weight,F(94889,321))
        self.assertLess(qnorm,F(15,7)**2)
        self.assertEqual((qnorm-F(15,7))/weight,F(16201,1992669))
        self.assertGreater(F(cert['floor_rational']),F(1,128))
        c1f=np.array([[x.display() for x in row] for row in c1]);c2f=np.array([[x.display() for x in row] for row in c2])
        w=np.zeros((8,8));w[7,0]=1
        for j in range(1,8):w[j-1,j]=-1
        matrix,nrhs=reversing_moment_system(w,c1f,c2f,np.diag([(-1)**(i&1) for i in range(8)]))
        np.testing.assert_allclose(matrix,np.array(numeric_columns).T,atol=1e-13,rtol=0)
        np.testing.assert_allclose(nrhs,np.array(rhs,dtype=float),atol=1e-13,rtol=0)

    def test_invalid_diagnostic_inputs(self):
        for bad in (np.ones((2,3)),np.zeros((2,2)),np.array([[np.nan,0],[0,1]])):
            with self.assertRaises(ValueError):phase_reflection(bad,np.eye(2),[1,1])
        with self.assertRaises(ValueError):phase_reflection(np.eye(2),np.eye(2),[1,0])
        with self.assertRaises(ValueError):nearest_reversible(-np.eye(2),np.diag([1,-1]))
        with self.assertRaises(ValueError):recover_phase_gauge(np.eye(2),np.eye(2),tolerance=True)
        with self.assertRaises(ValueError):cycle_discrepancy(np.eye(2),np.eye(2),(0,0))
        with self.assertRaises(ValueError):hermitian_coordinates(64)


if __name__=='__main__':unittest.main()
