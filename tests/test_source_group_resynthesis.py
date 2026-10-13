"""Completed-group gauges, full recovery and endpoint-only allocation."""
from fractions import Fraction as F
import unittest
import numpy as np
from compiler_robust_hopf.source_group_resynthesis import (
    Group,fixed_core_ledger,gauge_error_identity,raw_connector_counterexample,
    rational_root_word,
)
from compiler_robust_hopf.source_resynthesis import I,X,Z,norm


def blockdiag(*matrices):
    result=np.zeros((sum(len(m) for m in matrices),)*2,complex);p=0
    for m in matrices:result[p:p+len(m),p:p+len(m)]=m;p+=len(m)
    return result


def unitary(rng,n):
    return np.linalg.qr(rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)))[0]


def perturb(rng,n,epsilon):
    a=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n));a=(a+a.conj().T)/2
    values,vectors=np.linalg.eigh(a);values/=max(abs(values))
    return (vectors*np.exp(1j*epsilon*values))@vectors.conj().T


class SourceGroupTests(unittest.TestCase):
    def test_exact_rational_full_word(self):
        for ts in ((F(0),F(0),F(0)),(F(1,2),F(-1,3),F(2,5))):
            W=rational_root_word(2,ts)
            np.testing.assert_array_equal(W.T@W,np.eye(4,dtype=int))

    def test_raw_connector_keeps_literal_sign(self):
        original,omitted,x=raw_connector_counterexample()
        np.testing.assert_array_equal(original,-x)
        np.testing.assert_array_equal(omitted,x)
        self.assertEqual(norm(original-omitted),2)

    def test_full_group_gauge_identity(self):
        rng=np.random.default_rng(10601)
        for _ in range(6):
            target=np.kron(unitary(rng,4),np.eye(2));g=np.kron(np.eye(4),unitary(rng,2))
            J=np.eye(32,8);G=blockdiag(g,unitary(rng,24))
            V=perturb(rng,32,1e-3)@blockdiag(target,unitary(rng,24))
            left,right=gauge_error_identity(V,J,G,target)
            np.testing.assert_allclose(left,right,atol=3e-12,rtol=0)
            self.assertAlmostEqual(norm(left),norm(V@J-J@target),places=11)

    def test_prefix_controlled_gauge(self):
        rng=np.random.default_rng(10602)
        target=blockdiag(*[np.kron(unitary(rng,2),I) for _ in range(2)])
        g=blockdiag(*[np.kron(I,unitary(rng,2)) for _ in range(2)])
        J=np.eye(32,8);G=blockdiag(g,unitary(rng,24))
        V=perturb(rng,32,.003)@blockdiag(target,unitary(rng,24))
        a,b=gauge_error_identity(V,J,G,target)
        np.testing.assert_allclose(a,b,atol=3e-12,rtol=0)

    def test_different_completed_gauges_need_no_connector(self):
        rng=np.random.default_rng(10603);J=np.eye(32,8)
        for stages in (2,3,5,9):
            actual=np.eye(32,dtype=complex);ideal=np.eye(8,dtype=complex);budget=0
            for j in range(stages):
                # Prefix alternates between the two logical bits.
                target=blockdiag(*[np.kron(unitary(rng,2),I) for _ in range(2)])
                g=blockdiag(*[np.kron(I,unitary(rng,2)) for _ in range(2)])
                if j%2:
                    swap=np.eye(8)[[0,1,4,5,2,3,6,7]]
                    target=swap@target@swap.T;g=swap@g@swap.T
                G=blockdiag(g,unitary(rng,24))
                V=perturb(rng,32,.001)@blockdiag(target,unitary(rng,24))
                K=G.conj().T@V@G
                budget+=norm(V@J-J@target);actual=K@actual;ideal=target@ideal
            self.assertLessEqual(norm(actual@J-J@ideal),budget+2e-12)

    def test_gauge_rejects_unproved_premises(self):
        J=np.eye(4,2)
        with self.assertRaises(ValueError):gauge_error_identity(np.eye(4),J,np.kron(X,I),I)
        with self.assertRaises(ValueError):gauge_error_identity(np.eye(4),J,np.kron(I,X),Z)
        with self.assertRaises(ValueError):gauge_error_identity(np.zeros((4,4)),J,np.eye(4),I)

    def test_actual_half_blocks_and_recovery(self):
        cases=(F(0),F(1,4096),F(-3,8192),(F(1,8192),F(-3,16384),F(5,32768)))
        for m in (16,20,24):
            for ts in cases:
                g=Group(2,ts,m);J=g.embed();out=g.apply(J)
                expected=np.kron((np.eye(g.L)+g.expected)/2,np.eye(g.core))
                np.testing.assert_allclose(out[:g.data],expected,atol=2e-12,rtol=0)
                self.assertLess(norm(g.recover(J)-J@np.kron(g.W,np.eye(g.core))),10*2.**-m+2e-12)

    def test_actual_group_gauge_including_leakage(self):
        rng=np.random.default_rng(10604)
        for ts in (F(0),F(1,4096),(F(1,8192),F(-3,16384),F(5,32768))):
            g=Group(2,ts,20);J=g.embed();B=unitary(rng,g.core)
            input_gauge=np.kron(np.eye(g.L),B)
            transformed=g.recover(J@input_gauge).reshape(-1,g.core,g.data)
            transformed=np.einsum('ij,ajk->aik',B.conj().T,transformed).reshape(g.size,g.data)
            goal=J@np.kron(g.W,np.eye(g.core))
            self.assertAlmostEqual(norm(transformed-goal),norm(g.recover(J)-goal),places=10)

    def test_actual_inverse_on_occupied_work(self):
        rng=np.random.default_rng(10605);g=Group(2,F(-3,8192),16)
        v=rng.normal(size=(g.size,3))+1j*rng.normal(size=(g.size,3));v/=np.linalg.norm(v)
        np.testing.assert_allclose(g.apply(g.apply(v),True),v,atol=3e-12,rtol=0)

    def test_fixed_core_every_literal_ledger(self):
        for n in range(3,1025):
            d=fixed_core_ledger(n)
            self.assertTrue(d['self_query_fits']);self.assertGreaterEqual(d['other_dirty'],2)
            self.assertLessEqual(d['error_over_eta'],F(15,64))
            self.assertEqual(n+2+d['core']+d['other_dirty'],d['total_wires'])

    def test_fixed_core_small_cases_and_domain(self):
        self.assertEqual(fixed_core_ledger(3)['core'],16)
        self.assertEqual(fixed_core_ledger(3)['other_dirty'],2)
        for n in (True,2,3.5):
            with self.assertRaises(ValueError):fixed_core_ledger(n)
        with self.assertRaises(ValueError):Group(4,F(0),20)
