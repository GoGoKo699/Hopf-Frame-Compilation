"""Finite witnesses for the ordered-bank results; no endpoint emitter claim."""
import itertools
import math
import unittest
from fractions import Fraction as F

import numpy as np

from compiler_robust_hopf.ordered_block_banks import (
    Band, band_schedule, coarse_exponent, fixed_factor_schedule,
    grouped_factors, history_geometry, history_sources,
    interpolation_coefficients, ledger, modular_witness_exponent,
    ordered_product, prefix_matchings, prefix_matrix, refined_schedule,
    validate_bands,
)
from compiler_robust_hopf.native_ordered_banks import (
    Emitter, CompactEmitter, EchoProgram, apply_switches, benes,
    emit_readout, emit_prefix_bank, inverse, ref_word_matrix, simulate,
)
from compiler_robust_hopf.source_completions import (
    cholesky_in_order, normalization_two, postorder,
    regularized_completion, square_root_witness,
)


def sample_layers(n, seed=11, delta=1/256):
    rng = np.random.default_rng(seed)
    N = 1 << n
    X = np.array([[0,1],[1,0]], complex)
    Y = np.array([[0,-1j],[1j,0]], complex)
    Z = np.diag([1,-1])
    targets, coarse = [], []
    for d in range(n):
        W, Q = np.eye(N, dtype=complex), np.eye(N, dtype=complex)
        for p in range(1 << d):
            angle = rng.uniform(-math.pi, math.pi)
            if p % 4 == 0:
                angle = 0 if p % 8 == 0 else math.pi
            r = np.array([[np.cos(angle),-np.sin(angle)], [np.sin(angle),np.cos(angle)]], complex)
            axis = rng.normal(size=3)
            axis /= np.linalg.norm(axis)
            e = delta * 2.0**(d-n) / 4
            perturb = np.cos(e)*np.eye(2) + 1j*np.sin(e)*(axis[0]*X+axis[1]*Y+axis[2]*Z)
            a, b = p << (n-d), (2*p+1) << (n-d-1)
            W[np.ix_([a,b],[a,b])] = r
            Q[np.ix_([a,b],[a,b])] = perturb @ r
        targets.append(W)
        coarse.append(Q)
    return targets, coarse


def shared_factors(E):
    """Independent bounded reference construction of canonical balanced factors."""
    N = len(E)
    A, B = np.zeros((N,3*N-2), complex), np.zeros((N,3*N-2), complex)
    for i, z in enumerate(np.diag(E)):
        root = math.sqrt(abs(z))
        A[i,i], B[i,i] = (z/root if root else 0), root
    k = N
    for d in range(N.bit_length()-1):
        width = N >> d
        for start in range(0,N,width):
            left, right = slice(start,start+width//2), slice(start+width//2,start+width)
            for rows, cols in ((left,right),(right,left)):
                u, s, vh = np.linalg.svd(E[rows,cols],full_matrices=False)
                A[rows,k], B[cols,k] = u[:,0]*np.sqrt(s[0]), vh[0].conj()*np.sqrt(s[0])
                if len(s)>1:
                    assert s[1] < 1e-10
                k += 1
    return A/math.sqrt(2), B/math.sqrt(2)


class OrderedBankTests(unittest.TestCase):
    def test_validation(self):
        for n in (False,0,-1,2.5):
            with self.assertRaises(ValueError): band_schedule(n)
        with self.assertRaises(ValueError): Band(1,1)
        with self.assertRaises(ValueError): band_schedule(3,4)
        with self.assertRaises(ValueError): validate_bands(3,[Band(0,1),Band(2,3)])
        with self.assertRaises(ValueError): fixed_factor_schedule(3,0)

    def test_simple_exact_ledgers(self):
        for n in range(3,257):
            bands = band_schedule(n)
            record, N = ledger(n,bands), 1<<n
            self.assertLess(record['squared_size'],9*N)
            self.assertLessEqual(len(bands),1+(n-1).bit_length())
            self.assertLessEqual(6*sum(b.active_modes for b in bands[:-1]),N)
            self.assertLess(sum(b.squared_size for b in bands[:-1]),N)
            self.assertLessEqual(6*record['completed_modes'],19*N)

    def test_refined_exact_ledgers(self):
        for n in range(3,257):
            bands = refined_schedule(n)
            p = (((n-1).bit_length())-1).bit_length()
            record, N = ledger(n,bands),1<<n
            self.assertLess(record['squared_size']*(1<<p),(4*(1<<p)+5)*N)
            self.assertLessEqual(6*(1<<p)*record['history_modes'],12*(1<<p)*N-5*N)
            self.assertLess(record['completed_modes'],4*N)
        for n in range(3,35):
            for p in range(n+1):
                bands = band_schedule(n,p)
                self.assertLess(ledger(n,bands)['completed_modes'],4*(1<<n))

    def test_all_partitions_minimum_volume(self):
        for n in range(1,11):
            minimum = None
            for cut in range(1<<(n-1)):
                endpoints=[0]+[i for i in range(1,n) if cut>>(i-1)&1]+[n]
                bands=[Band(a,b) for a,b in zip(endpoints,endpoints[1:])]
                volume=ledger(n,bands)['squared_size']
                self.assertGreaterEqual(volume,4*((1<<n)-1))
                minimum=volume if minimum is None else min(minimum,volume)
            self.assertEqual(minimum,4*((1<<n)-1))

    def test_fixed_factor_exponent_integer_form(self):
        for n in range(1,65):
            for k in range(1,min(n,8)+1):
                bands=fixed_factor_schedule(n,k)
                # Each band exponent is at most ceil(n*2**k/(2**k-1))+1.
                threshold=(n*(1<<k)+(1<<k)-2)//((1<<k)-1)+1
                self.assertTrue(all(2*b.stop-b.start <= threshold for b in bands))
                self.assertLessEqual(len(bands),k)

    def test_full_factorization_and_support(self):
        for n in range(2,7):
            W,_=sample_layers(n)
            for bands in (band_schedule(n),band_schedule(n,min(2,n))):
                grouped=grouped_factors(W,bands)
                np.testing.assert_allclose(ordered_product(grouped),ordered_product(W),atol=2e-14)
                for band,matrix in zip(bands,grouped):
                    expected=np.eye(1<<n,dtype=complex)
                    for support in band.supports(n):
                        expected[np.ix_(support,support)]=matrix[np.ix_(support,support)]
                    np.testing.assert_allclose(expected,matrix,atol=2e-14)

    def test_nested_history_coefficients_and_sources(self):
        for n in range(2,6):
            target,coarse=sample_layers(n,seed=n)
            for bands in ([Band(d,d+1) for d in range(n)],band_schedule(n),refined_schedule(n) if n>=3 else band_schedule(n)):
                F1,G=grouped_factors(target,bands),grouped_factors(coarse,bands)
                C=interpolation_coefficients(F1,G)
                geom=history_geometry(n,bands)
                SA,SB=history_sources(C,geom)
                W,Q=ordered_product(target),ordered_product(coarse)
                np.testing.assert_allclose(sum(C),W,atol=2e-14)
                np.testing.assert_allclose(C[0],Q,atol=2e-14)
                np.testing.assert_allclose(SA.conj().T@SA,np.eye(1<<n),atol=2e-14)
                np.testing.assert_allclose(SB.conj().T@SB,np.eye(1<<n),atol=2e-14)
                np.testing.assert_allclose(SA.conj().T@SB,W/2,atol=2e-14)
                np.testing.assert_allclose(geom.unitary()[:,:1<<n],SB,atol=2e-14)
                total=sum(np.linalg.norm(f-g,2) for f,g in zip(F1,G))
                for k in range(1,len(C)):
                    outside=[i for i,h in enumerate(geom.heights) if h<k]
                    if outside: self.assertLess(np.linalg.norm(C[k][:,outside]),1e-14)
                    self.assertLessEqual(np.linalg.norm(C[k],2),total**k/math.factorial(k)+1e-14)

    def test_history_rejects_wrong_support(self):
        geom=history_geometry(3)
        C=[np.eye(8),np.zeros((8,8)),np.zeros((8,8)),np.zeros((8,8))]
        C[3][0,1]=.01
        with self.assertRaises(ValueError): history_sources(C,geom)
        with self.assertRaises(ValueError): interpolation_coefficients([],[])

    def test_prefix_copy_count(self):
        for n in range(1,16):
            match=prefix_matchings(n)
            self.assertEqual(sum(map(len,match)),(1<<n)-n-1)
            for ell,pairs in enumerate(match):
                self.assertEqual(len(pairs),(n-1-ell)*(1<<ell))
                labels=[x for a,b,p in pairs for x in (a,b)]
                self.assertEqual(len(labels),len(set(labels)))

    def test_small_angle_threshold(self):
        for n in range(8,1025):
            for bands in (band_schedule(n),refined_schedule(n)):
                self.assertLessEqual(modular_witness_exponent(n,len(bands)),(1<<n)-n)
        for n in range(1,100):
            s=coarse_exponent(n)-6
            self.assertGreaterEqual(s*s,1<<n)
            self.assertLess((s-1)**2,1<<n)

    def test_exact_parity_witness(self):
        for k in range(1,7):
            u=F(1,128)
            sine=2*u/(1+u*u)
            for subset in range(1<<k):
                degree=subset.bit_count()
                if degree>=k: continue
                average=sum(math.prod(signs)*math.prod(signs[j] for j in range(k) if subset>>j&1)
                            for signs in itertools.product((-1,1),repeat=k))
                self.assertEqual(average,0)
            amp=sum(F(math.prod(signs),1)*math.prod(sign*sine for sign in signs)
                    for signs in itertools.product((-1,1),repeat=k))/(1<<k)
            self.assertEqual(amp,sine**k)
            self.assertGreaterEqual(amp,u**k)


class SourceCompletionTests(unittest.TestCase):
    def test_completion_full_space_and_cut_ranks(self):
        for n in (2,3,4):
            target,coarse=sample_layers(n,seed=71+n,delta=1/1024)
            E=ordered_product(target)@ordered_product(coarse).conj().T
            f,g=shared_factors(E)
            c=regularized_completion(f,g)
            N=1<<n
            np.testing.assert_allclose(c.source_f.conj().T@c.source_g,E/2,atol=1e-13)
            for source,unitary in ((c.source_f,c.unitary_f),(c.source_g,c.unitary_g)):
                np.testing.assert_allclose(source.conj().T@source,np.eye(N),atol=1e-13)
                np.testing.assert_allclose(unitary.conj().T@unitary,np.eye(4*N),atol=1e-13)
                np.testing.assert_allclose(unitary[:,:N],source,atol=1e-13)
            U=c.unitary_f.conj().T@c.unitary_g
            expected=np.vstack([E,np.zeros((3*N,N))])
            np.testing.assert_allclose(normalization_two(U,N)[:,:N],expected,atol=2e-13)
            self.assertGreaterEqual(np.linalg.eigvalsh(c.joint_gram)[0],1/8-1e-12)
            self.assertLessEqual(np.linalg.eigvalsh(c.joint_gram)[-1],5/4+1e-12)
            for b in range(1,N):
                inside=list(range(b,b+(b&-b)))
                outside=[i for i in range(N) if i not in inside]
                for matrix,rc,cc,rank in ((c.joint_gram,2,2,2),(c.joint_cholesky,2,2,2),
                                         (c.deficit_f,1,1,1),(c.source_f,4,1,3),(c.unitary_f,4,4,12)):
                    for rows,cols in ((inside,outside),(outside,inside)):
                        r=[ch*N+i for ch in range(rc) for i in rows]
                        col=[ch*N+i for ch in range(cc) for i in cols]
                        singular=np.linalg.svd(matrix[np.ix_(r,col)],compute_uv=False)
                        if len(singular)>rank: self.assertLess(singular[rank],1e-11)

    def test_completion_input_validation(self):
        with self.assertRaises(ValueError): regularized_completion(np.eye(3),np.eye(3))
        with self.assertRaises(ValueError): regularized_completion(np.eye(2),np.eye(2))
        with self.assertRaises(ValueError): postorder(0)
        with self.assertRaises(ValueError): cholesky_in_order(np.eye(2),(0,0))
        with self.assertRaises(ValueError): normalization_two(np.zeros((2,2)),1)
        with self.assertRaises(ValueError): square_root_witness(-1)

    def test_common_factor_and_rejected_completion(self):
        rng=np.random.default_rng(115)
        def unitary():
            return np.linalg.qr(rng.normal(size=(8,8))+1j*rng.normal(size=(8,8)))[0]
        A,B,C=unitary(),unitary(),unitary()
        np.testing.assert_allclose((C@A).conj().T@(C@B),A.conj().T@B,atol=1e-13)
        U=A.conj().T@B
        D1=np.eye(8,dtype=complex);D2=D1.copy()
        D1[2:,2:]=np.linalg.qr(rng.normal(size=(6,6))+1j*rng.normal(size=(6,6)))[0]
        D2[2:,2:]=np.linalg.qr(rng.normal(size=(6,6))+1j*rng.normal(size=(6,6)))[0]
        np.testing.assert_allclose(normalization_two(D1@U@D2,2),
                                   D1@normalization_two(U,2)@D2,atol=1e-13)

    def test_fixed_leading_failure_overlap(self):
        # Two 2D ranges in a 3D bank necessarily intersect.
        X=math.sqrt(3/8)*np.eye(3)[:,:2]
        Y=math.sqrt(3/8)*np.eye(3)[:,1:]
        self.assertAlmostEqual(np.linalg.norm(X.conj().T@Y,2),3/8)

    def test_principal_root_rank_inflation_exact_series(self):
        # Work with F(u)=2D(u). All coefficients and convolutions are rational.
        zero=lambda: np.full((4,4),F(0),dtype=object)
        identity=zero()
        for i in range(4): identity[i,i]=F(1)
        J=np.full((4,4),F(1,4),dtype=object)
        M=np.array([[1,1,0,0],[1,2,0,0],[0,0,2,1],[0,0,1,1]],dtype=object)
        M=np.array([[F(int(x),4) for x in row] for row in M],dtype=object)
        coefficients=[identity]+[zero() for _ in range(8)]
        for k in range(1,9):
            coefficients[k]=-2*((-1)**((k-1)//2))*J if k%2 else -2*((-1)**((k-2)//2))*M
        root=[identity]
        for k in range(1,9):
            root.append((coefficients[k]-sum((root[j]@root[k-j] for j in range(1,k)),zero()))/F(2))
        for k in range(9):
            self.assertTrue(np.array_equal(sum((root[j]@root[k-j] for j in range(k+1)),zero()),coefficients[k]))
        determinant=[]
        for k in range(7):
            # sqrt(D)=sqrt(F)/sqrt(2), so the 2x2 determinant is halved.
            determinant.append(sum(root[j][0,2]*root[k-j][1,3]-root[j][0,3]*root[k-j][1,2]
                                   for j in range(k+1))/2)
        self.assertEqual(determinant[:6],[F(0)]*6)
        self.assertEqual(determinant[6],F(-1,2048))
        D=square_root_witness(.05)
        eig,U=np.linalg.eigh(D); R=(U*np.sqrt(eig))@U.T
        self.assertLess(np.linalg.svd(D[:2,2:],compute_uv=False)[1],1e-15)
        self.assertGreater(np.linalg.svd(R[:2,2:],compute_uv=False)[1],1e-11)

    def test_square_root_witness_is_an_actual_tree_residual(self):
        N=16;u=.005
        w=np.zeros(N);w[:8]=[.5,0,1/(2*math.sqrt(2)),1/(2*math.sqrt(2)),
                             1/(2*math.sqrt(2)),1/(2*math.sqrt(2)),.5,0]
        e=np.eye(N)[:,8]
        E=np.eye(N)-2*u*u/(1+u*u)*(np.outer(w,w)+np.outer(e,e))
        E+=2*u/(1+u*u)*(np.outer(e,w)-np.outer(w,e))
        f,g=shared_factors(E)
        V=np.zeros((N,4));V[0,0]=1;V[2:4,1]=1/math.sqrt(2)
        V[4:6,2]=1/math.sqrt(2);V[6,3]=1
        np.testing.assert_allclose(V.T@(np.eye(N)-f@f.conj().T)@V,
                                   square_root_witness(u),atol=2e-14)
        self.assertLess(np.linalg.norm(E-np.eye(N),2),1/(8*math.sqrt(N)))

    def test_exact_cayley_finite_difference(self):
        # Gaussian rationals as pairs of Fractions; no symbolic dependency.
        add=lambda a,b:(a[0]+b[0],a[1]+b[1])
        mul=lambda a,b:(a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0])
        def div(a,b):
            den=b[0]*b[0]+b[1]*b[1]
            return ((a[0]*b[0]+a[1]*b[1])/den,(a[1]*b[0]-a[0]*b[1])/den)
        for m in range(1,9):
            r=F(1,7);h=2*r/m;xs=[-r+j*h for j in range(m+1)]
            total=(F(0),F(0));den=(F(1),F(0));power=(F(1),F(0))
            for j,x in enumerate(xs):
                weight=(-1)**(m-j)*math.comb(m,j)
                total=add(total,mul((F(weight),F(0)),div((F(1),x),(F(1),-x))))
                den=mul(den,(F(1),-x))
            for _ in range(m): power=mul(power,(F(0),h))
            self.assertEqual(total,div(mul((F(2*math.factorial(m)),F(0)),power),den))
            for degree in range(m):
                self.assertEqual(sum((-1)**(m-j)*math.comb(m,j)*x**degree for j,x in enumerate(xs)),0)

    def test_averaged_block_loses_order(self):
        for t in (.01,.07,.3):
            r=np.array([[np.cos(t),-np.sin(t)],[np.sin(t),np.cos(t)]])
            U,V=np.eye(3),np.eye(3)
            U[:2,:2]=r;V[1:,1:]=r
            B=(U+V-np.eye(3))/2
            v=np.array([1,0,1])/math.sqrt(2);w=np.array([1,0,-1])/math.sqrt(2)
            np.testing.assert_allclose(B@v,np.cos(t)*v/2,atol=1e-15)
            np.testing.assert_allclose(B.T@v,np.cos(t)*v/2,atol=1e-15)
            self.assertAlmostEqual(w@(V@U)@v,-np.sin(t)**2/2)


class NativeOrderedBankTests(unittest.TestCase):
    def test_benes_permutations(self):
        for p in itertools.permutations(range(4)):
            self.assertEqual(apply_switches(benes(p),4),list(p))
        rng=np.random.default_rng(316)
        for w in range(1,9):
            for _ in range(3):
                p=rng.permutation(1<<w)
                self.assertEqual(apply_switches(benes(p),len(p)),list(p))
        with self.assertRaises(ValueError): benes([0,1,1,3])
        with self.assertRaises(ValueError): benes([0.,1.])

    def test_literal_toffoli_and_controlled_h(self):
        e=Emitter(3);e.ccx(0,1,2)
        actual=simulate(e.gates,np.eye(8))
        p=np.arange(8);p^=(((p&3)==3).astype(int)<<2)
        np.testing.assert_allclose(actual,np.eye(8)[p],atol=2e-15)
        self.assertEqual(sum(g[0] in ('t','tdg') for g in e.gates),7)
        e=Emitter(2);e.controlled('h',0,1)
        expected=np.eye(4);expected[np.ix_([1,3],[1,3])]=np.array([[1,1],[1,-1]])/math.sqrt(2)
        np.testing.assert_allclose(simulate(e.gates,np.eye(4)),expected,atol=2e-15)
        self.assertEqual(sum(g[0] in ('t','tdg') for g in e.gates),2)

    def test_full_dirty_boolean_interpreters(self):
        for G in ('x','h','z'):
            e=Emitter(3);table=[1,0,1,1];e.mux(2,table,G)
            D=1<<e.qubits;expected=np.eye(D,dtype=complex)
            for i in range(D):
                if i&4 or not table[i&3]: continue
                pair=[i,i|4]
                block={'x':np.array([[0,1],[1,0]]),'h':np.array([[1,1],[1,-1]])/math.sqrt(2),'z':np.diag([1,-1])}[G]
                expected[np.ix_(pair,pair)]=block
            np.testing.assert_allclose(simulate(e.gates,np.eye(D)),expected,atol=2e-13)

    def test_echo_programs_literal_symmetry(self):
        X=np.array([[0,1],[1,0]],complex)
        for seed in ((),('h','t'),('t','h','tdg','h'),('h','t','x','tdg','h')):
            program=EchoProgram(seed);D=program.matrix()
            self.assertAlmostEqual(np.linalg.det(D),1,places=13)
            np.testing.assert_allclose(X@D@X,D.conj().T,atol=2e-14)
        with self.assertRaises(ValueError): EchoProgram(('arbitrary',))
        with self.assertRaises(ValueError): emit_prefix_bank(3,[[[1]]])
        with self.assertRaises(ValueError): CompactEmitter(3).compressed_matching([(0,1)],[[1]])

    def test_native_readout_occupied_and_initialized_vectors(self):
        rng=np.random.default_rng(1001)
        for n in (2,3):
            e=emit_readout(n);geometry=history_geometry(n)
            M=4*(1<<n);D=1<<e.qubits
            state=rng.normal(size=(D,2))+1j*rng.normal(size=(D,2))
            state[(np.arange(D)%M)>=(1<<n),0]=0
            state/=np.linalg.norm(state,axis=0)
            ideal=np.einsum('ij,bjk->bik',geometry.unitary(),state.reshape(-1,M,2)).reshape(D,2)
            actual=simulate(e.gates,state)
            np.testing.assert_allclose(actual,ideal,atol=2e-11)
            self.assertLessEqual(e.qubits-(n+2),(1<<n)+n+7)
            if n==3:
                self.assertEqual(len(e.gates),44940)
                self.assertEqual(sum(g[0] in ('t','tdg') for g in e.gates),18492)

    def test_native_prefix_occupied_vectors_and_inverse(self):
        rng=np.random.default_rng(441)
        n=3
        seeds=[('h','t'),('t','h','tdg','h')]
        programs=[[EchoProgram(seeds[(d+p)%2]) for p in range(1<<d)] for d in range(n-1)]
        e=emit_prefix_bank(n,programs)
        ideal_matrix=prefix_matrix(n,[[p.matrix() for p in row] for row in programs])
        D=1<<e.qubits;M=1<<(n+1)
        state=rng.normal(size=(D,2))+1j*rng.normal(size=(D,2));state/=np.linalg.norm(state,axis=0)
        ideal=np.einsum('ij,bjk->bik',ideal_matrix,state.reshape(-1,M,2)).reshape(D,2)
        actual=simulate(e.gates,state)
        np.testing.assert_allclose(actual,ideal,atol=2e-11)
        np.testing.assert_allclose(simulate(inverse(e.gates),actual),state,atol=4e-11)
        self.assertLessEqual(e.qubits-(n+1),(1<<n)+n+7)

    def test_native_input_validation(self):
        with self.assertRaises(ValueError): Emitter(0)
        with self.assertRaises(ValueError): Emitter(3).mux(0,[1],'x')
        with self.assertRaises(ValueError): Emitter(3).gate('cx',0,0)
        with self.assertRaises(ValueError): simulate([],np.zeros(3))
        with self.assertRaises(ValueError): simulate([('foo',0)],np.zeros(2))


if __name__ == '__main__':
    unittest.main()
