"""Exact signed-Pauli sector fixtures; no dense high-dimensional unitary claim."""
import unittest
from tests.exact_matrices import Q2
from tests.test_one_band_copy_feasibility import _conjugate
from tests.field_sector_tableaux import (
    BitslicedTableau, decoder_word, recoupler_word, phase_coefficients,
    sector_tableau, materialize_word, query_word, query_parameter_word,
    target_word, helper_shift, mul, trace)


class FieldSectorTableauTests(unittest.TestCase):
    def test_signed_conditional_gate_updates_against_exact_paulis(self):
        # Every two-qubit Hermitian Pauli, in every helper sector. Noncanonical
        # input batches exercise intermediate Y phases, not only X/Z supports.
        word=(('H',0),('CZC',0,1,1),('CX',1,0),('ZC',1,2),
              ('H',1),('CZC',0,1,3),('H',0),('CX',0,1))
        for start in range(0,16,4):
            table=BitslicedTableau(2,2)
            for q in range(2):
                table.x[q]=[table.ones if (label>>q)&1 else 0
                            for label in range(start,start+4)]
                table.z[q]=[table.ones if (label>>(q+2))&1 else 0
                            for label in range(start,start+4)]
            table.apply(word)
            for character in range(4):
                concrete=materialize_word(word,character)
                for column,label in enumerate(range(start,start+4)):
                    x,z,sign=table.at(character)[column]
                    self.assertEqual(_conjugate({label:Q2(1)},concrete,2),
                                     {x|(z<<2):Q2(-1 if sign else 1)})

    def test_quadratic_phase_expansion_on_every_bounded_column(self):
        for level in (2,3):
            ell=1<<(level-2);k=2*ell;mask=(1<<ell)-1
            for eta in range(1<<ell):
                for kind in ('before','after'):
                    terms=phase_coefficients(level,kind)
                    for basis in range(1<<(2*k)):
                        u=basis&((1<<k)-1);v=basis>>k
                        d0,d1=u&mask,u>>ell;a0,a1=v&mask,v>>ell
                        if kind=='before':
                            f=(mul(level,u,v)>>ell)^mul(level-1,d0,d1)
                        else:
                            f=mul(level-1,d0,a1)^mul(level-1,d1,a0)
                        expected=trace(level-1,mul(level-1,eta,f))
                        actual=0
                        for gate in terms:
                            selected=(eta&gate[-1]).bit_count()&1
                            value=1
                            for q in gate[1:-1]:value&=(basis>>q)&1
                            actual^=selected&value
                        self.assertEqual(actual,expected)

    def test_complete_signed_cocycle_and_actual_inverse(self):
        for p in (2,3,4):
            n=1<<p;h=(1<<(p-1))-1;word=decoder_word(p)
            inverse=BitslicedTableau(n,h).apply(word).apply(reversed(word))
            self.assertTrue(inverse.same(BitslicedTableau(n,h)))
            # Address generators suffice by the additive cocycle induction;
            # bit-slicing includes every helper character at once.
            for addr in (1<<j for j in range(p)):
                left=BitslicedTableau(n,h).apply(query_word(p,addr))
                left.apply(word,helper_shift(p,addr))
                right=BitslicedTableau(n,h).apply(word).apply(target_word(p,addr))
                self.assertTrue(left.same(right),(p,addr))
            self.assertTrue(sector_tableau(p).same(BitslicedTableau(n,h).apply(word)))

    def test_local_cocycle_uses_the_entire_lower_field(self):
        # This is stronger than checking the p-dimensional original address
        # subspace: every lambda in the lower field follows from its F2 basis.
        for p in (2,3,4):
            ell=1<<(p-2);n=1<<p;start=ell-1
            local=[(g[0],*g[1:-1],g[-1]>>start)
                   if g[0] in ('CZC','ZC') else g
                   for g in recoupler_word(p)]
            for lam in (1<<j for j in range(ell)):
                left=BitslicedTableau(n,ell).apply(query_parameter_word(p,lam))
                left.apply(local,lam)
                children=query_parameter_word(p-1,lam)
                children+=query_parameter_word(p-1,lam,1<<(p-1))
                right=BitslicedTableau(n,ell).apply(local).apply(children)
                self.assertTrue(left.same(right),(p,lam))
            theta=1<<ell
            left=BitslicedTableau(n,ell).apply(query_parameter_word(p,theta)).apply(local)
            right=BitslicedTableau(n,ell).apply(local).apply(target_word(p,1<<(p-1)))
            self.assertTrue(left.same(right),p)


if __name__=='__main__':unittest.main()
