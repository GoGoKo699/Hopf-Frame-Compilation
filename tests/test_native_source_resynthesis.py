"""Exact selected-tail words and full dirty-query regression checks."""
from fractions import Fraction as F
import unittest
import numpy as np
from compiler_robust_hopf.native_source_resynthesis import (
    apply_macros_basis,inv_macros,self_query,selected_tail_word,tail_matrices,
)
from compiler_robust_hopf.native_coarse_fixture import inverse_word,simulate_columns,gate_counts
from tests.test_conjugate_recovery import exact_native
from verification.fault_tolerant.exact_arithmetic import Q2,CQ2


class NativeSourceResynthesisTests(unittest.TestCase):
    def test_selected_tail_complete_native_matrices(self):
        for a,m,indices in ((0,3,(0,)),(1,3,(0,2)),(1,4,(1,3)),(2,4,(0,1,2,3))):
            word=selected_tail_word(a,m,indices);size=1<<(a+m)
            full=simulate_columns(word,np.eye(size),nqubits=a+m)
            expected=np.zeros_like(full);tails=tail_matrices(m)
            for y,p in enumerate(indices):
                expected[y*(1<<m):(y+1)*(1<<m),y*(1<<m):(y+1)*(1<<m)]=tails[p]
            np.testing.assert_allclose(full,expected,atol=3e-12,rtol=0)
            np.testing.assert_allclose(simulate_columns(inverse_word(word),full,nqubits=a+m),np.eye(size),atol=3e-12,rtol=0)

    def test_selected_tail_native_exact_and_literal_inverse(self):
        word=selected_tail_word(1,3,(0,2))
        self.assertEqual((len(word),gate_counts(word)['T_count']),(374,12))
        actual=exact_native([(name.lower(),*bits) for name,bits in word],4)
        back=exact_native([(name.lower(),*bits) for name,bits in inverse_word(word)],4)
        zero=CQ2();coeff=[CQ2(Q2(0,F(1,2))),CQ2(F(1,2)),CQ2(F(1,2))]
        for column in range(16):
            prefix=column>>3;x=column&7;expected=[zero]*16
            for j in (range(3) if prefix==0 else (2,)):
                sign=(-1)**sum((x>>(2-k))&1 for k in range(j))
                row=(prefix<<3)|(x^(1<<(2-j)))
                expected[row]=sign*(coeff[j] if prefix==0 else CQ2(1))
            self.assertEqual([actual[row][column] for row in range(16)],expected)
            self.assertEqual([back[row][column] for row in range(16)],expected)

    def test_self_queries_return_every_dirty_input(self):
        rng=np.random.default_rng(10606)
        for a,m in ((0,3),(1,3),(2,4),(2,5),(3,6)):
            table=rng.integers(0,2,size=(1<<a,m),dtype=np.uint8)
            word=self_query(a,m,table);permutation=apply_macros_basis(word,a+m)
            for value,image in enumerate(permutation):
                address=value>>m;mask=sum(int(b)<<(m-1-j) for j,b in enumerate(table[address]))
                self.assertEqual(image,value^mask)
            inverse=apply_macros_basis(inv_macros(word),a+m)
            np.testing.assert_array_equal(inverse[permutation],np.arange(1<<(a+m)))

    def test_reference_entangled_selected_tail(self):
        rng=np.random.default_rng(10607);word=selected_tail_word(1,4,(0,3))
        reference=rng.normal(size=(32,7))+1j*rng.normal(size=(32,7));reference/=np.linalg.norm(reference)
        actual=simulate_columns(word,reference,nqubits=5)
        expected=np.vstack([tail_matrices(4)[0]@reference[:16],tail_matrices(4)[3]@reference[16:]])
        np.testing.assert_allclose(actual,expected,atol=3e-12,rtol=0)

    def test_native_input_validation(self):
        for args in ((2,3,(0,1,2,0)),(1,3,(0,3)),(1,3,(0,)),(1,3,(0,1.5))):
            with self.assertRaises(ValueError):selected_tail_word(*args)
        with self.assertRaises(ValueError):self_query(1,3,[[0,0,2],[1,0,0]])
