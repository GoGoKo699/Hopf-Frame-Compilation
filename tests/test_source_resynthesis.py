"""Source resynthesis algebra; finite checks are not endpoint certificates."""
from fractions import Fraction as F
from math import comb
import unittest
import numpy as np
from compiler_robust_hopf.source_resynthesis import (
    I,X,Y,Z,actual_table_family,chebyshev_fraction,construct_model,
    encoded_coefficient,exact_pivots,exp_axis,lie_labels,majoranas,move,norm,
    program_word,scalar_power_lower_bound,sign_mask,single_mask_alignment,
    source_data,weights,
)


class SourceResynthesisTests(unittest.TestCase):
    def test_geometric_weights_exact(self):
        for m in range(2,129):
            w=weights(m)
            self.assertEqual(sum(w),1);self.assertEqual(w[-1],w[-2])
            for p in range(m):self.assertEqual(sum(w[p:]),F(1,1<<p))

    def test_canonical_rounding_and_literal_zeros(self):
        for m in range(3,11):
            for c in (F(-1),F(0),F(1),F(1,7),F(1,4),F(-2,9)):
                f=sign_mask(c,m)
                self.assertLessEqual(abs(encoded_coefficient(f)-c),F(1,1<<(m-2)))
            self.assertEqual(sign_mask(0,m).tolist(),[1]+[0]*(m-1))
            self.assertEqual(encoded_coefficient(sign_mask(0,m)),0)

    def test_actual_table_binary_masks(self):
        for v in (3,6,10):
            for h in (1,2,5,11):
                m=v+h+3;cs,fs=actual_table_family(m,h,v)
                base=np.zeros(m,dtype=np.uint8);base[1:v+1]=1
                np.testing.assert_array_equal(fs[0],base)
                for j in range(1,h+1):
                    expected=base.copy();expected[v+j]^=1
                    np.testing.assert_array_equal(fs[j],expected)
                for c,f in zip(cs,fs):self.assertEqual(encoded_coefficient(f),c)
                rows=[[1]+[1-2*int(f[j]) for f in fs] for j in range(m)]
                self.assertEqual(len(exact_pivots(rows)),h+2)

    def test_exact_denominator_growth(self):
        for m in range(4,25):
            d=m-2;c=F((1<<d)-1,1<<d)
            for R in range(1,16):
                value=chebyshev_fraction(c,2*R)
                exponent=2*R*(m-3)+1
                self.assertEqual(value.denominator,1<<exponent)
                self.assertTrue(value.numerator&1)
                self.assertEqual(scalar_power_lower_bound(m,R),2*exponent)

    def test_scalar_second_and_third_stages(self):
        for m in (4,5,6):
            masks=np.zeros((1,m),dtype=np.uint8);masks[0,-2]=1
            _,_,_,mat,_,_=source_data(m,masks)
            M,N=mat;c=1-2**(2-m);D=(M@N-N@M)/2
            S=c*np.eye(1<<(m+1))+np.kron(D,X)
            for power in (2,3,5):
                accepted=np.linalg.matrix_power(S,power)[::2,::2]
                expected=float(chebyshev_fraction(F(c),power))*np.eye(1<<m)
                np.testing.assert_allclose(accepted,expected,atol=2e-12,rtol=0)
            self.assertGreater(abs((2*c*c-1)-c*c),.01)

    def test_full_axis_alignment_and_recompiled_powers(self):
        for m in (4,5,6):
            A,theta=single_mask_alignment(m);gs=majoranas(m);G=gs[-2]
            M=sum(float(w)**.5*g for w,g in zip(weights(m),gs))
            a=2**(-(m-1)/2);Ystar=1j*(M@G-a*np.eye(1<<m))/np.sqrt(1-a*a)
            np.testing.assert_allclose(A@G@A.conj().T,Ystar,atol=3e-12,rtol=0)
            S=exp_axis(theta,np.kron(Ystar,X));B=np.kron(A,I)
            for R in (2,3,17):
                target=B@exp_axis(R*theta,np.kron(G,X))@B.conj().T
                np.testing.assert_allclose(np.linalg.matrix_power(S,R),target,atol=1e-11,rtol=0)

    def test_class_rotation_handles_endpoints(self):
        for c in (0.,.1,.5,.9,1.):
            B=c*Z+np.sqrt(1-c*c)*X
            for half in (False,True):
                U=move(Z,B,c,half)
                expected=(Z+B)/np.sqrt(2+2*c) if half else B
                np.testing.assert_allclose(U@Z@U.conj().T,expected,atol=2e-12,rtol=0)
                np.testing.assert_allclose(U.conj().T@U,I,atol=2e-12,rtol=0)

    def test_common_basis_and_both_central_sectors(self):
        for m in (4,5,6):
            cases=[np.zeros((2,m),dtype=np.uint8),np.eye(m,dtype=np.uint8)[:2],np.eye(m,dtype=np.uint8)[:3]]
            for masks in cases:
                model=construct_model(m,masks)
                self.assertEqual(model.active,model.d//2)
                self.assertLess(model.max_error,2e-11)
                for scalar in model.compact_scalars:
                    np.testing.assert_allclose(scalar.conj().T@scalar,np.eye(len(scalar)),atol=2e-12,rtol=0)

    def test_changing_address_and_actual_inverse(self):
        for m in (4,5):
            masks=np.zeros((2,m),dtype=np.uint8);masks[0,1]=1;masks[1,-2]=1
            model=construct_model(m,masks)
            for steps in (2,3,17):
                full,compact,error,_=program_word(model,steps,seed=steps+20)
                self.assertLess(error,3e-11)
                np.testing.assert_allclose(full.conj().T@full,np.eye(len(full)),atol=3e-11,rtol=0)

    def test_signal_only_lie_dimension(self):
        for p in range(1,9):self.assertEqual(len(lie_labels(p,1)),comb(p+3,2))

    def test_address_lie_envelope_is_separate(self):
        for p in range(1,7):self.assertEqual(len(lie_labels(p,2)),(1<<(p+4))-(1<<(p%2)))

    def test_exact_rank_does_not_use_float_threshold(self):
        self.assertEqual(exact_pivots([[1,1,F(1,1<<200)],[0,0,1]]),(0,2))
        self.assertEqual(exact_pivots([[1,2],[2,4]]),(0,))

    def test_source_input_validation(self):
        for call in (lambda:weights(True),lambda:majoranas(20),lambda:sign_mask(2,5),
                     lambda:source_data(4,[[0,0,0,256]]),lambda:source_data(4,[[0,0,0,.5]]),
                     lambda:actual_table_family(6,5),lambda:exact_pivots([[1],[1,2]]),
                     lambda:scalar_power_lower_bound(3,2)):
            with self.assertRaises(ValueError):call()
