"""Actual conditional-suffix half-blocks in faithful compact core coordinates.

The group coefficients and dyadic masks are calculated exactly with Fraction.
The geometric source may have many dirty qubits; its finite Clifford span is
represented faithfully on a small core for numerical full-word diagnostics.
This is not a simulation of 2**m physical core amplitudes or a gate emitter.
"""
from __future__ import annotations
from fractions import Fraction
import math
import numpy as np

from .source_resynthesis import (I,X,Y,Z,H,weights,kronall,norm,
                                  sign_mask,exact_pivots,integer)


def apply_axes(state,matrix,axes):
    axes=tuple(axes);order=axes+tuple(j for j in range(state.ndim) if j not in axes)
    a=np.transpose(state,order);shape=a.shape
    return np.transpose((matrix@a.reshape(matrix.shape[1],-1)).reshape(shape),np.argsort(order))


def rational_root_word(s,t):
    """All entries are Fractions; complete addressed outward word."""
    integer(s, 'diagnostic height', 1, 4)
    size = 1 << s
    W = np.array([[Fraction(i == j) for j in range(size)] for i in range(size)], dtype=object)
    ts = list(t) if isinstance(t,(tuple,list)) else [t]+[Fraction(0)]*(size-2)
    if len(ts) != size-1:
        raise ValueError('One rational tangent per edge is required.')
    j=0
    for depth in range(s):
        stride=1 << (s-depth-1)
        for prefix in range(1 << depth):
            u=Fraction(ts[j]); j+=1
            c=(1-u*u)/(1+u*u); v=2*u/(1+u*u)
            x=2*prefix*stride; y=x+stride
            old_x=W[x].copy(); old_y=W[y].copy()
            W[x]=c*old_x-v*old_y; W[y]=v old_x+c*old_y
    return W


def labels(s):
    ls=[('diag',0,0,z) for z in range(4)]
    for kind,depth in [('zero',0)]+[('node',j) for j in range(s)]:
        for rev in (0,1):ls.extend((kind,depth,rev,z) for z in range(4))
    K=1<<(len(ls)-1).bit_length()
    return ls+[('pad',0,0,0)]*(K-len(ls))


def component(x,phase):
    # These fixtures have exactly rational real residual entries.
    if phase==0:return max(x,Fraction(0))
    if phase==1:return max(-x,Fraction(0))
    return Fraction(0)


def group_tables(s,t,m):
    if s%2:raise ValueError('rational scaling fixture requires even s')
    W=rational_root_word(s,t);E=W-np.eye(1<<s,dtype=int);ls=labels(s);K=len(ls);L=1<<s
    tables=[]
    for kind,d,rev,ph in ls:
        row=[]
        for u in range(L):
            if kind=='pad':c=Fraction(0)
            elif kind=='diag':c=K*component(Fraction(E[u,u]),ph)
            else:
                D=L if kind=='zero' else 1<<(s-d)
                lo=(u//D)*D;v=0 if kind=='zero' else lo+D//2
                valid=(u!=0) if kind=='zero' else (u not in (lo,v))
                if not valid:c=Fraction(0)
                else:
                    value=Fraction(E[v,u] if rev else E[u,v])
                    # For odd-width types all selected entries here are zero.
                    a=component(value,ph)
                    if a and (s-d)%2 and kind!='zero':raise AssertionError('unexpected irrational nonzero scale')
                    c=K*math.isqrt(D)*a
            assert 0<=c<=Fraction(1,4)
            row.append(c)
        tables.append(row)
    masks=[[sign_mask(c,m) for c in row] for row in tables]
    unique={tuple(f) for row in masks for f in row}
    return np.array(W,dtype=complex),ls,tables,masks,list(sorted(unique))


def compact_scalars(m,mask_list):
    ws=weights(m);a=np.sqrt([float(w) for w in ws]);vectors=[a]
    vectors.extend(a*(1-2*np.array(f,dtype=int)) for f in mask_list)
    sg=[[1]+[1-2*int(f[j]) for f in mask_list] for j in range(m)]
    pivots=list(exact_pivots(sg));d=len(pivots)
    Q,_=np.linalg.qr(np.column_stack([vectors[j] for j in pivots]),mode='reduced')
    if np.dot(Q[:,0],a)<0:Q[:,0]*=-1
    b=d//2;p=d-1
    if b==0:
        paulis=[]
    else:
        paulis=[]
        for j in range(b):
            paulis.extend([kronall([Z]*j+[X]+[I]*(b-j-1)),kronall([Z]*j+[Y]+[I]*(b-j-1))])
        paulis=paulis[:p]
    dim=1<<b;out={};maxgram=0
    for f,v in zip(mask_list,vectors[1:]):
        c=sum(w*(1-2*int(bit)) for w,bit in zip(ws,f))
        coeff=Q.T@v
        assert abs(coeff[0]-float(c))<1e-12
        D=-1j*sum((coeff[j+1]*paulis[j] for j in range(p)),np.zeros((dim,dim),complex))
        S0=float(c)*np.eye(2*dim)+np.kron(X,D)
        maxgram=max(maxgram,norm(S0.conj().T@S0-np.eye(2*dim)))
        out[tuple(f)]=(S0,float(c))
    return out,b,d,maxgram


class Group:
    def __init__(self,s,t,m):
        integer(s, 'compact group height', 2, 2)
        integer(m, 'diagnostic source width', 8, 64)
        self.s=s;self.m=m
        self.W,self.labels,self.tables,self.masks,unique=group_tables(s,t,m)
        self.scalars,self.active,self.d,self.source_error=compact_scalars(m,unique)
        self.h=len(unique);self.L=1<<s;self.K=len(self.labels);self.core=1<<self.active
        self.shape=(2,self.K,2,2,self.L,self.core);self.size=math.prod(self.shape)
        self.data=self.L*self.core
        self.Hlabel=kronall([H]*(self.K.bit_length()-1))
        self.atoms={}
        for kind,d,_,_ in self.labels:
            if kind not in ('zero','node') or (kind,d) in self.atoms:continue
            varying=s if kind=='zero' else s-d;D=1<<varying
            marker=0 if kind=='zero' else D//2
            U=np.kron(np.eye(self.L//D),kronall([H]*varying))
            shift=np.eye(self.L)[:,np.arange(self.L)^marker]
            pre=np.zeros((2*self.L,2*self.L))
            for a in range(2):
                for x in range(self.L):pre[(a^(x%D!=marker))*self.L+x,a*self.L+x]=1
            self.atoms[kind,d]=np.kron(np.eye(2),U@shift)@pre
        self.expected=np.zeros((self.L,self.L),complex)
        self.phases=[1,-1,1j,-1j]
        for idx,(kind,d,rev,ph) in enumerate(self.labels):
            chat=np.array([self.scalars[tuple(f)][1] for f in self.masks[idx]])
            if kind in ('diag','pad'):term=np.diag(chat)
            else:
                atom=self.atoms[kind,d][:self.L,:self.L]
                term=atom.conj().T@np.diag(chat) if rev else np.diag(chat)@atom
            self.expected+=self.phases[ph]*term/self.K

    def scalar(self,state,idx,inverse):
        out=state.copy()
        for u,f in enumerate(self.masks[idx]):
            S0=self.scalars[tuple(f)][0]
            if inverse:S0=S0.conj().T
            out[:,:,u,:,:]=apply_axes(out[:,:,u,:,:],S0,(0,2))
        return out

    def apply(self,columns,inverse=False):
        state=np.asarray(columns).reshape(self.shape+(-1,)).copy()
        state=apply_axes(state,H,(0,));selected=apply_axes(state[1],self.Hlabel,(0,))
        for j,(kind,d,rev,ph) in enumerate(self.labels):
            term=selected[j]
            if kind in ('diag','pad'):
                term=self.scalar(term,j,inverse)
            else:
                atom=self.atoms[kind,d]
                # Forward source D, reverse D^dagger S^dagger, then invert whole word.
                effective_reverse=bool(rev)^bool(inverse)
                if effective_reverse:
                    term=self.scalar(term,j,True)
                    term=apply_axes(term,atom.conj().T,(1,2))
                else:
                    term=apply_axes(term,atom,(1,2))
                    term=self.scalar(term,j,False)
            term*=np.conj(self.phases[ph]) if inverse else self.phases[ph]
            selected[j]=term
        state[1]=apply_axes(selected,self.Hlabel,(0,));state=apply_axes(state,H,(0,))
        return state.reshape(self.size,-1)

    def embed(self):
        out=np.zeros((self.size,self.data),complex);out[:self.data]=np.eye(self.data);return out

    def recover(self,v):
        out=self.apply(v);out[:self.data]*=-1
        out=self.apply(out,True);out[:self.data]*=-1
        return -self.apply(out)


def fixed_core_ledger(n):
    """Endpoint-only, two-clean ledger. It does not improve the T bound."""
    integer(n, 'logical width', 3, 4096)
    N=1 << n; ceil_log=(n-1).bit_length(); mu=N+6+ceil_log
    remainder=N+n+7-mu
    return {'n':n,'N':N,'core':mu,'other_dirty':remainder,
            'max_address':n+2,'self_query_fits':mu >= 2*(n+2),
            'error_over_eta':Fraction(15*n,1 << (6+ceil_log)),
            'total_wires':N+2*n+9}


def gauge_error_identity(V,J,G,F):
    """Return both complete residuals; reject a violated gauge premise."""
    matrices=[np.asarray(a,dtype=complex) for a in (V,J,G,F)]
    V,J,G,F=matrices
    if (J.ndim != 2 or V.shape != (len(J),len(J)) or G.shape != V.shape
            or F.shape != (J.shape[1],J.shape[1]) or
            any(not np.isfinite(a).all() for a in matrices)):
        raise ValueError('Finite matrices with compatible dimensions are required.')
    for U in (V,G,F):
        if not np.allclose(U.conj().T@U,np.eye(len(U)),atol=1e-10,rtol=0):
            raise ValueError('Unitary circuit and target matrices are required.')
    if not np.allclose(J.conj().T@J,np.eye(J.shape[1]),atol=1e-10,rtol=0):
        raise ValueError('J must be an isometry.')
    g=J.conj().T@G@J
    if not np.allclose(G@J,J@g,atol=1e-10,rtol=0):
        raise ValueError('The gauge must preserve the initialized code.')
    if not np.allclose(g@F,F@g,atol=1e-10,rtol=0):
        raise ValueError('The input gauge must commute with the target.')
    left=G.conj().T@V@G@J-J@F
    right=G.conj().T@(V@J-J@F)@g
    return left,right


def raw_connector_counterexample():
    """Exact Gaussian-integer matrices for the norm-two sign error."""
    g0=np.kron(X,I); g1=np.kron(Z,X)
    D=g0@g1
    source=np.kron(np.eye(2),np.kron(X,D))
    prefix_x=np.kron(X,np.eye(8))
    P=np.kron(Z,I)
    gauge=np.zeros((16,16),complex)
    gauge[:8,:8]=np.eye(8); gauge[8:,8:]=np.kron(I,P)
    transported=gauge@source@gauge
    return source@prefix_x@source, transported@prefix_x@transported, prefix_x
