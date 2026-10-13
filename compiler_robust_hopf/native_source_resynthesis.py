"""Exact native selected-tail words and finite controlled-source diagnostics.

Big-endian wires. Prefix wires precede the arbitrary source core. Internal macro
words use mathematical big-endian indices; selected_tail_word returns the
repository native convention with wire zero the low integer bit. Internal output
uses H,S,Sdg,T,Tdg,CX only; X/Z/CCX below are expanded before output.
The generic fine theorem is proved in SOURCE_PROGRAM_RESYNTHESIS.md, not
inferred from these finite floating-point matrices.
"""
from __future__ import annotations
import math
from typing import Sequence
import numpy as np
from .source_resynthesis import integer, majoranas, weights
from .native_coarse_fixture import inverse_word, Word
from .paid_logical_encodings import expand_native


def source_basis_word(m):
    """Internal big-endian scalar-free loader; phases cancel in conjugations."""
    integer(m, 'source width', 2, 4096)
    return [g for j in range(m-1) for g in
            [('Sdg',j),('H',j),('H',j+1),('CX',j,j+1),('Tdg',j+1),
             ('CX',j,j+1),('H',j),('S',j),('H',j+1)]]


def inv_macros(word):
    naming = {'Sdg':'SDG','Tdg':'TDG'}
    reverse_names = {'SDG':'Sdg','TDG':'Tdg'}
    raw = tuple((naming.get(g[0],g[0]), tuple(g[1:])) for g in word)
    return [(reverse_names.get(op,op),*qs) for op,qs in inverse_word(raw)]


def expand(word):
    """Reuse the retained exact Toffoli and phase expansion, not a new copy."""
    naming = {'Sdg':'SDG','Tdg':'TDG'}
    reverse_names = {'SDG':'Sdg','TDG':'Tdg'}
    raw = tuple((naming.get(g[0],g[0]), tuple(g[1:])) for g in word)
    return [(reverse_names.get(op,op),*qs) for op,qs in expand_native(raw)]


def tcount(word):return sum(g[0] in ('T','Tdg') for g in word)


def xor_query(address,targets,selectors,table):
    """Two-pass exact dirty traversal; targets/selectors are disjoint."""
    address=list(address);targets=list(targets);selectors=list(selectors)
    a=len(address);tab=np.asarray(table)
    if not np.isin(tab,(0,1)).all():raise ValueError('binary table required')
    tab=tab.astype(np.uint8)
    if tab.shape!=(1<<a,len(targets)) or np.any(tab>1):raise ValueError('binary table shape')
    if set(address)&(set(targets)|set(selectors)) or set(targets)&set(selectors):raise ValueError('overlapping registers')
    if len(selectors)<a:raise ValueError('too few selectors')
    if not a:return [('X',t) for t,b in zip(targets,tab[0]) if b]
    def traversal(omit):
        out=[]
        def visit(level,row):
            for bit in (0,1):
                if level==0:
                    update=[] if omit else ([('X',selectors[0])] if bit==0 else [])+[('CX',address[0],selectors[0])]
                else:
                    update=([('X',address[level])] if bit==0 else [])+[('CCX',selectors[level-1],address[level],selectors[level])]+([('X',address[level])] if bit==0 else [])
                out.extend(update)
                child=2*row+bit
                if level+1==a:
                    out.extend(('CX',selectors[a-1],t) for t,b in zip(targets,tab[child]) if b)
                else:visit(level+1,child)
                out.extend(inv_macros(update))
        visit(0,0)
        return out
    return traversal(False)+inv_macros(traversal(True))


def self_query(a,m,table):
    """Whole-core XOR with two self-borrowing partial banks."""
    integer(a, 'prefix width', 0, 12);integer(m, 'source width', 2, 4096)
    tab=np.asarray(table)
    if tab.shape!=(1<<a,m) or not np.isin(tab,(0,1)).all():raise ValueError('binary table shape')
    tab=tab.astype(np.uint8)
    if m<2*a:raise ValueError('need m>=2a')
    core=list(range(a,a+m))
    if not a:return [('X',core[j]) for j,b in enumerate(tab[0]) if b]
    mid=m//2;A=core[:mid];B=core[mid:]
    return xor_query(range(a),A,B[:a],tab[:,:mid])+xor_query(range(a),B,A[:a],tab[:,mid:])


def selected_gamma(a,m,indices):
    """Exact SELECT Gamma_p, with -1 denoting identity; no clean work."""
    if len(indices)!=(1<<a):raise ValueError('wrong row count')
    xtab=np.zeros((1<<a,m),np.uint8);ztab=xtab.copy()
    for x,p in enumerate(indices):
        if p<-1 or p>=m:raise ValueError('invalid index')
        if p>=0:xtab[x,p]=1;ztab[x,:p]=1
    had=[('H',a+j) for j in range(m)]
    return had+self_query(a,m,ztab)+had+self_query(a,m,xtab)


def sign_phase(a,m,values):
    """(-1)^f(prefix), using one arbitrary core bit with exact return."""
    if m<max(a+1,2):raise ValueError('insufficient dirty data for sign echo')
    target=a;sel=list(range(a+1,a+1+a))
    q=xor_query(range(a),[target],sel,np.array(values,dtype=np.uint8)[:,None])
    return q+[('Z',target)]+inv_macros(q)+[('Z',target)]


def selected_tail(a,m,indices):
    """Exact full-space SELECT A_p, including p=0 and the final tail."""
    integer(a, 'prefix width', 0, 12);integer(m, 'source width', 2, 4096)
    if any(isinstance(p,bool) or not isinstance(p,(int,np.integer)) for p in indices):raise ValueError('integer tail indices required')
    if m<max(2*a,a+1,2):raise ValueError('source/prefix width')
    if len(indices)!=(1<<a) or any(p<0 or p>=m for p in indices):raise ValueError('tail indices')
    U=[(g[0],*(a+j for j in g[1:])) for g in source_basis_word(m)]
    invU=inv_macros(U)
    p=[j if j>0 else -1 for j in indices]
    prev=[j-1 if j>0 else -1 for j in indices]
    atzero=[0 if j==0 else -1 for j in indices]
    V=invU+selected_gamma(a,m,p)+U
    # For p>0, A_p=-V Gamma_(p-1) V. For p=0, this part is I.
    tailpart=V+selected_gamma(a,m,prev)+V+sign_phase(a,m,[j>0 for j in indices])
    zeropart=invU+selected_gamma(a,m,atzero)+U
    return expand(tailpart+zeropart)


def tail_matrices(m):
    gs=majoranas(m);a=np.sqrt([float(w) for w in weights(m)])
    return [sum((2**(p/2))*a[j]*gs[j] for j in range(p,m)) for p in range(m)]


def apply_macros_basis(word,n):
    ids=np.arange(1<<n,dtype=np.int64);out=ids.copy()
    for g in word:
        if g[0]=='X':out^=1<<(n-1-g[1])
        elif g[0]=='CX':
            _,c,t=g;out^=((out>>(n-1-c))&1)<<(n-1-t)
        elif g[0]=='CCX':
            _,c,d,t=g;out^=(((out>>(n-1-c))&1)&((out>>(n-1-d))&1))<<(n-1-t)
        else:raise ValueError('permutation gates only')
    return out


def block_diag(blocks):
    sizes=[len(b) for b in blocks];out=np.zeros((sum(sizes),sum(sizes)),complex);i=0
    for b in blocks:out[i:i+len(b),i:i+len(b)]=b;i+=len(b)
    return out



def selected_tail_word(address_bits, core_bits, indices) -> Word:
    """Exact native SELECT of normalized tails, including p=0 and final p.

    Prefix occupies the high address_bits integer bits; core occupies the
    low core_bits bits. Source coordinate j is physical wire core_bits-1-j.
    Every input is arbitrary. No initialized or additional helper is used.
    """
    total=address_bits+core_bits
    names={'Sdg':'SDG','Tdg':'TDG'}
    raw=selected_tail(address_bits,core_bits,tuple(indices))
    return tuple((names.get(g[0],g[0]),tuple(total-1-q for q in g[1:])) for g in raw)
