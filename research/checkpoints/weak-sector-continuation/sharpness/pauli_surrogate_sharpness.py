"""Exact bounded sharpness fixture for a commuting-Pauli surrogate bank.

No SVD, floating arithmetic, initialized helper, query-degree or native-depth
claim is used. The unitary has a 1/sqrt(k) block; every observable block is
rational with denominator k and is checked with integer matrices.
"""
from pathlib import Path
import hashlib,json,platform,time,unittest
import numpy as np


def gamma_action(m,j):
    dimension=1<<m
    columns=np.arange(dimension,dtype=np.int64)
    rows=columns^(1<<j)
    prefix=(1<<j)-1
    signs=np.fromiter(((-1 if (b&prefix).bit_count()&1 else 1)
                       for b in range(dimension)),dtype=np.int64,count=dimension)
    return rows,signs


def sectors(m):
    a=np.arange(1<<(m-1),dtype=np.int64)
    parity=np.fromiter((b.bit_count()&1 for b in a.tolist()),dtype=np.int64,count=len(a))
    even=a|(parity<<(m-1))
    odd=a|((parity^1)<<(m-1))
    return even,odd


def sector_sum(m,selected):
    """Return S_even,odd for S=sum gamma_j, in data-label order on both sectors."""
    selected=tuple(selected);even,odd=sectors(m);half=len(even)
    row_of={int(b):a for a,b in enumerate(even)}
    block=np.zeros((half,half),dtype=np.int64)
    for j in selected:
        rows,signs=gamma_action(m,j)
        for column,b in enumerate(odd):block[row_of[int(rows[b])],column]+=signs[b]
    return block


def z_values(states,i):return 1-2*((states>>i)&1)


def observable_blocks(m,selected):
    """B_i = diag(Z_i|even, C_i/k), while A_i=diag(Z_i|even,Z_i|even)."""
    selected=tuple(selected);k=len(selected)
    assert k and len(set(selected))==k and all(0<=j<m for j in selected)
    even,odd=sectors(m);s=sector_sum(m,selected)
    z_even=[z_values(even,i) for i in range(m)]
    z_odd=[z_values(odd,i) for i in range(m)]
    c=[(s*z[np.newaxis,:])@s.T for z in z_odd]
    return k,s,z_even,z_odd,c


def verify_case(m,selected):
    assert 2 <= m <= 8, "This exact fixture is intentionally bounded"
    selected=tuple(selected);k,s,z_even,z_odd,c=observable_blocks(m,selected)
    half=len(s);identity=np.eye(half,dtype=np.int64)
    # Operator norms bound each intermediate entry by k^t; this also bounds
    # every unaccumulated int64 matrix-product sum well below overflow.
    assert half*max(4*k*k,k**m) < 1<<61
    assert np.array_equal(s@s.T,k*identity)
    assert np.array_equal(s.T@s,k*identity)
    # This proves the full unitary U=diag(I,S/sqrt(k)) in ordered parity bases.
    errors=[]
    for i,ci in enumerate(c):
        assert np.array_equal(ci,ci.T)
        assert np.array_equal(ci@ci,k*k*identity)
        assert np.array_equal(ci@s,k*s*z_odd[i][np.newaxis,:])
        # Every original bitstring has its complete distinct joint eigenvalue
        # pattern; the preceding identity covers all odd input columns.
        difference=k*np.diag(z_even[i])-ci
        if i in selected:
            assert np.array_equal(difference@difference,4*k*identity)
            errors.append({'observable':i,'norm_squared_numerator':4,'norm_squared_denominator':k})
        else:
            assert not np.any(difference)
            errors.append({'observable':i,'norm_squared_numerator':0,'norm_squared_denominator':1})
    for i in range(m):
        for j in range(i):assert np.array_equal(c[i]@c[j],c[j]@c[i])
    assert np.array_equal(np.prod(np.stack(z_even),axis=0),np.ones(half,dtype=np.int64))
    assert np.array_equal(np.prod(np.stack(z_odd),axis=0),-np.ones(half,dtype=np.int64))
    # Rational matrix products have denominator k^m. Sizes here keep integers safe.
    assert k**m < 1<<40
    product=identity.copy()
    for ci in c:product=product@ci
    assert np.array_equal(product,-(k**m)*identity)
    # Hence product_i A_i=I but product_i B_i=Z_copy, with no loss of a bank bit.
    return {'m':m,'selected_gamma_indices':list(selected),'k':k,
            'full_bank_columns':1<<m,'odd_sector_columns':half,
            'complete_unitarity_gram':'S_eo S_oe=S_oe S_eo=kI',
            'all_commuting_involution_blocks':'PASS',
            'all_joint_eigenvalue_columns':'PASS',
            'surrogate_relation':'product A_i=I',
            'independent_bank_product':'product B_i=Z_copy',
            'errors':errors}


class PauliSurrogateSharpnessTests(unittest.TestCase):
    def test_jordan_wigner_anticommutation_every_column(self):
        for m in (4,8):
            actions=[gamma_action(m,j) for j in range(m)]
            columns=np.arange(1<<m,dtype=np.int64)
            for j,(rows,signs) in enumerate(actions):
                self.assertTrue(np.array_equal(rows[rows],columns))
                self.assertTrue(np.array_equal(signs*signs[rows],np.ones_like(signs)))
                for i in range(j):
                    ri,si=actions[i]
                    self.assertTrue(np.array_equal(ri[rows],rows[ri]))
                    self.assertTrue(np.array_equal(signs*si[rows],-si*signs[ri]))

    def test_complete_surrogate_banks_and_exact_norms(self):
        for m in (4,8):
            for selected in (tuple(range(m)),tuple(range(1,m))):
                with self.subTest(m=m,selected=selected):verify_case(m,selected)

    def test_seed_queries_and_full_boolean_degree(self):
        for m in (4,8):
            def word(addr):
                if addr==0:return []
                if addr<m-1:return [('CX',0,addr),('CX',addr,0),('CX',0,addr)]
                return [('CX',j,0) for j in range(1,m-1)]
            def image(addr,basis):
                for _,control,target in word(addr):
                    basis^=((basis>>control)&1)<<target
                return basis
            labels=[1<<i for i in range(m-1)]+[(1<<(m-1))-1]
            self.assertEqual(__import__('functools').reduce(int.__xor__,labels,0),0)
            for addr,expected in enumerate(labels):
                for basis in range(1<<m):
                    self.assertEqual(image(addr,image(addr,basis)),basis)
                    self.assertEqual(image(addr,basis)&1,(expected&basis).bit_count()&1)
            # The top p-bit ANF coefficient is the XOR over every address.
            # On any Z_j with 1<=j<=m-2 it is z0+zj, hence not zero.
            # These individually Clifford queries have degree exactly p=log2(m).
            for j in range(1,m-1):
                top=0
                for addr in range(m):
                    label=sum(((image(addr,1<<bit)>>j)&1)<<bit for bit in range(m))
                    top^=label
                self.assertEqual(top,1^(1<<j))
            _,_,ze,_,ci=observable_blocks(m,tuple(range(1,m)))
            self.assertTrue(np.array_equal(ci[0],(m-1)*np.diag(ze[0])))



def run():
    start=time.monotonic();report=[]
    for m in (4,8):
        for selected in (range(m),range(1,m)):report.append(verify_case(m,selected))
    data={'status':'PASS','python':platform.python_version(),'numpy':np.__version__,
          'arithmetic':'Exact int64 matrices; integer bounds checked; no SVD or floating norm estimates',
          'cases':report,'seconds':time.monotonic()-start,
          'scope':'Commuting observables and a full unitary only. The optional Clifford query family has Boolean degree exactly p=log2(m), no additive group premise, and no coherent native compiler ledger.'}
    Path(__file__).with_suffix('.json').write_text(json.dumps(data,indent=2)+'\n')
    return data


if __name__=='__main__':
    import sys
    if '--test' in sys.argv:
        sys.argv.remove('--test');unittest.main()
    else:print(json.dumps(run(),indent=2))
