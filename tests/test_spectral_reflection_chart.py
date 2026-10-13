"""Spectral identities, literal offsets and query scope; no fine-factor oracle."""
from functools import lru_cache
import math
import unittest
import numpy as np

from compiler_robust_hopf.spectral_reflection_chart import (
    alternative_baseline_word, borrowed_mcx, centered_cycle_vector,
    complement, contour_projector, matching_generator, signed_cycle,
    signed_cycle_word, spectral_chart,
)
from compiler_robust_hopf.native_coarse_fixture import gate_counts,inverse_word,simulate_columns
from compiler_robust_hopf.paid_logical_encodings import apply_boolean
from tests.test_conjugate_recovery import exact_native
from verification.fault_tolerant.exact_arithmetic import CQ2


def skew_exponential(k,t):
    values,vectors=np.linalg.eigh(1j*k)
    return ((vectors*np.exp(-1j*t*values))@vectors.conj().T).real


@lru_cache(maxsize=2)
def exact_cycle(inverse=False):
    word=signed_cycle_word(4)
    if inverse:word=inverse_word(word)
    return exact_native([(g.lower(),*qs) for g,qs in word],5)


class SpectralReflectionTests(unittest.TestCase):
    def test_signed_cycle_spectrum_and_reference_reflections(self):
        for size in (4,8,16,32):
            p,s=signed_cycle(size),complement(size)
            np.testing.assert_array_equal(np.linalg.matrix_power(p,size),-np.eye(size))
            np.testing.assert_array_equal(s@p@s,p.T)
            np.testing.assert_array_equal((p@s).T,p@s)
            self.assertEqual(np.trace(s),0);self.assertEqual(np.trace(p@s),0)
            for j in range(size):
                v=centered_cycle_vector(size,j)
                np.testing.assert_allclose(p@v,np.exp(1j*(2*j+1)*np.pi/size)*v,atol=2e-14,rtol=0)
                np.testing.assert_allclose(s@v,v.conj(),atol=2e-14,rtol=0)

    def test_regular_chart_complete_matrix_and_bounds(self):
        rng=np.random.default_rng(10521)
        for size in (4,8,16,32):
            gap=2*np.sin(np.pi/size)
            for kind in ('zero','dense','matching'):
                k=rng.normal(size=(size,size));k=k-k.T
                if kind=='zero':k[:]=0
                if kind=='matching':
                    k[:]=0
                    for j in range(0,size,2):k[j+1,j],k[j,j+1]=1,-1
                if np.linalg.norm(k):k/=np.linalg.norm(k,2)
                e=skew_exponential(k,gap/(40*np.sqrt(size)))
                r=spectral_chart(e);h1,h2=r['H1'],r['H2'];eye=np.eye(size)
                for h in (h1,h2):
                    np.testing.assert_allclose(h,h.T,atol=3e-11,rtol=0)
                    np.testing.assert_allclose(h@h,eye,atol=3e-11,rtol=0)
                    self.assertAlmostEqual(np.trace(h),0,places=9)
                np.testing.assert_allclose(h2@h1@r['P'].T,e,atol=3e-11,rtol=0)
                alpha=4*np.sqrt(size)*r['delta']/gap
                self.assertLessEqual(np.linalg.norm(r['A']-eye,2),alpha+2e-11)
                self.assertLessEqual(np.linalg.norm(r['O']-eye,2),2*alpha+2e-11)
                self.assertLessEqual(np.linalg.norm(h1-r['S'],2),4*alpha+2e-11)

    def test_fixed_circle_quadrature(self):
        size=8;k=matching_generator(8)[:size,:size]
        e=skew_exponential(k,1e-3);r=spectral_chart(e)
        for nodes in (4,8,12):
            bound=(4/11)**nodes/(1-(4/11)**nodes)
            for j in range(size):
                q=contour_projector(r['T'],r['reference'][j],r['gap']/3,nodes)
                self.assertLessEqual(np.linalg.norm(q-r['projectors'][j],2),bound+2e-13)

    def test_exact_dirty_boolean_echo_on_all_helper_values(self):
        count=0
        for controls in range(3,7):
            target=controls;dirty=tuple(range(controls+1,2*controls-1))
            word=borrowed_mcx(tuple(range(controls)),target,dirty)
            self.assertEqual(sum(g=='CCX' for g,_ in word),4*controls-8)
            for value in range(1<<(2*controls-1)):
                expected=value ^ ((1<<target) if value&((1<<controls)-1)==(1<<controls)-1 else 0)
                self.assertEqual(apply_boolean(value,word),expected);count+=1
        self.assertEqual(count,2720)

    def test_native_cycle_and_inverse_exact(self):
        word=signed_cycle_word(4)
        self.assertEqual((len(word),gate_counts(word)['T_count']),(142,63))
        full,back=exact_cycle(),exact_cycle(True)
        for col in range(32):
            x=col&15;dirty=col>>4
            row=dirty*16+(x+1)%16;sign=-1 if x==15 else 1
            for i in range(32):
                self.assertEqual(full[i][col],CQ2(sign if i==row else 0))
                self.assertEqual(back[i][col],full[col][i].conj())

    def test_native_cycle_scaling_and_registers(self):
        for m in range(4,33):
            word=signed_cycle_word(m)
            self.assertEqual(gate_counts(word)['T_count'],7*(2*m*m-6*m+1))
            self.assertEqual(max(q for _,qs in word for q in qs),2*m-4)
        self.assertEqual(gate_counts(signed_cycle_word(3))['T_count'],14)

    def test_clifford_alternate_baseline_exact(self):
        for m in (3,4,5):
            size=1<<m;word=alternative_baseline_word(m)
            self.assertEqual(gate_counts(word)['T_count'],0)
            expected=np.linalg.matrix_power(signed_cycle(size),size//4)@complement(size)
            exact=exact_native([(g.lower(),*qs) for g,qs in word],m)
            for i in range(size):
                for j in range(size):self.assertEqual(exact[i][j],CQ2(int(expected[i,j])))

    def test_forcing_coefficient_and_finite_sensitivity(self):
        for size in (8,16,32):
            total=2*size;p=signed_cycle(total);s=complement(total);k=matching_generator(size)
            fj=centered_cycle_vector(total,size//2-1);fk=centered_cycle_vector(total,size//2+1)
            self.assertAlmostEqual((fj.conj()@k@fk).real,-.5*np.cos(np.pi/(2*size)),places=12)
            self.assertLess(abs((fj.conj()@k@fk).imag),1e-12)
            f=np.cos(np.pi/size)*np.cos(np.pi/(2*size));d=2*np.sin(np.pi/size)
            forcing=abs(fj.conj()@(k@p@s+p@s@k)@fk.conj())
            self.assertAlmostEqual(forcing,f,places=12)
            t=1e-4/size;e=skew_exponential(k,t);r=spectral_chart(e)
            h=r['H1'];rho=np.linalg.norm(r['T']@h-h@r['T'].T,2)
            self.assertGreaterEqual((d+2*t)*np.linalg.norm(h-s,2)+t*t+rho+1e-13,f*t)
            self.assertGreaterEqual(np.linalg.norm(h-s,2)+1e-12,size*t/32)

    def test_alternative_baseline_cancels_full_target(self):
        for size in (8,16,32):
            p=signed_cycle(2*size);s=complement(2*size);j=np.linalg.matrix_power(p,size//2)@s
            k=matching_generator(size)
            np.testing.assert_allclose(j@k@j,-k,atol=2e-14,rtol=0)
            for t in (0,.03,-.4):
                e=skew_exponential(k,t);h1=p.T@j;h2=e@j
                np.testing.assert_allclose(h1,h1.T,atol=2e-12,rtol=0)
                np.testing.assert_allclose(h2,h2.T,atol=2e-12,rtol=0)
                np.testing.assert_allclose(h2@h1@p.T,e,atol=2e-12,rtol=0)
                self.assertAlmostEqual(np.linalg.norm(h2-j,2),np.linalg.norm(e-np.eye(2*size),2),places=12)

    def test_chart_width_and_query_error_exponents(self):
        for n in range(3,129):
            size=1<<n;root=math.isqrt(size);s=root+(root*root<size)+6
            self.assertGreaterEqual(2*s,4+3*(n+1))
            q=size//2+4
            self.assertLessEqual(q+(n+1)+1,size+n+6)
            self.assertEqual((n+1)+2+(size+n+6),size+2*n+9)
            if n>=8:self.assertLessEqual(s+12-size,n)

    def test_query_hybrid_with_actual_inverse_and_controls(self):
        rng=np.random.default_rng(10522);size=8;p=signed_cycle(size);k=matching_generator(8)[:8,:8]
        t=.02;e=skew_exponential(k,t);query=e@p
        controlled0=np.block([[p,np.zeros_like(p)],[np.zeros_like(p),np.eye(size)]])
        controlled1=np.block([[query,np.zeros_like(p)],[np.zeros_like(p),np.eye(size)]])
        a0=a1=np.eye(2*size,dtype=complex)
        for step in range(1,8):
            z=rng.normal(size=(2*size,2*size))+1j*rng.normal(size=(2*size,2*size));mix=np.linalg.qr(z)[0]
            q0=controlled0 if step%2 else controlled0.conj().T
            q1=controlled1 if step%2 else controlled1.conj().T
            a0=q0@mix@a0;a1=q1@mix@a1
            self.assertLessEqual(np.linalg.norm(a1-a0,2),step*abs(t)+2e-12)

    def test_invalid_chart_and_wire_inputs(self):
        with self.assertRaises(ValueError):signed_cycle(3)
        with self.assertRaises(ValueError):spectral_chart(-np.eye(4))
        with self.assertRaises(ValueError):spectral_chart(np.diag([1j,-1j,1,1]))
        with self.assertRaises(ValueError):borrowed_mcx((0,1,2),3,())
        with self.assertRaises(ValueError):borrowed_mcx((0,1),1,())
        for value in (True,2,257):
            with self.assertRaises(ValueError):signed_cycle_word(value)
        with self.assertRaises(ValueError):contour_projector(np.eye(4),1.,0.,4)


if __name__=='__main__':unittest.main()
