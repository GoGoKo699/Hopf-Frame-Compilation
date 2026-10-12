"""Exact native right-history and coarse-prefix emitters.

One shared Beneš/dirty-table implementation serves both constructions. Low
bits hold the occupied mode register; all higher helpers are arbitrary dirty
inputs. Returned gates are chronological and use only H,S,Sdg,T,Tdg,X,CX.
The fine bank, fine left source, and endpoint compiler are not emitted here.
"""
from __future__ import annotations
from dataclasses import dataclass
from numbers import Integral
from typing import Sequence
import math
import numpy as np
from .ordered_block_banks import (Band, _integer, history_geometry,
                                  prefix_matchings, prefix_matrix)

def benes(p):
    raw = tuple(p)
    if any(isinstance(x, bool) or not isinstance(x, Integral) for x in raw):
        raise ValueError('Permutation labels must be integers.')
    p=np.asarray(raw,dtype=int); m=len(p); w=m.bit_length()-1
    if sorted(p.tolist())!=list(range(m)) or m!=(1<<w): raise ValueError('not a power-of-two permutation')
    if w==0:return []
    if w==1:return [(0,[int(p[0])])]
    half=m//2; inv=np.argsort(p); col=[-1]*m
    for start in range(m):
        if col[start]>=0:continue
        col[start]=0; todo=[start]
        while todo:
            x=todo.pop()
            for y in (x^half,int(inv[int(p[x])^half])):
                if col[y]<0: col[y]=1-col[x];todo.append(y)
                elif col[y]!=1-col[x]:raise AssertionError('color conflict')
    first=col[:half]; last=[col[int(inv[v])] for v in range(half)]
    subs=[]
    for c in (0,1):
        q=[int(p[u if col[u]==c else u+half])%half for u in range(half)]
        subs.append(benes(q))
    middle=[]
    for (t,f0),(u,f1) in zip(*subs):
        assert t==u; middle.append((t,f0+f1))
    return [(w-1,first)]+middle+[(w-1,last)]

def apply_switches(stages,m):
    out=[]
    for x in range(m):
        y=x
        for t,f in stages:
            address=(y&((1<<t)-1))|((y>>(t+1))<<t)
            y^=int(f[address])<<t
        out.append(y)
    return out

class Emitter:
    def __init__(self,w):
        w = _integer(w, 'mode qubits', 1)
        self.w=w; self.r=w-1; self.ell=self.r//2; self.lam=1<<self.ell
        self.banks=list(range(w,w+self.lam));self.sel=list(range(w+self.lam,w+self.lam+self.r-self.ell))
        self.gates=[]
    @property
    def qubits(self):return self.w+self.lam+len(self.sel)
    def gate(self,op,*q):
        arity = 2 if op == 'cx' else 1
        if op not in ('h','s','sdg','t','tdg','x','cx') or len(q) != arity:
            raise ValueError('Unsupported elementary gate or arity.')
        if len(set(q)) != len(q) or any(not isinstance(x, Integral) or isinstance(x, bool)
                                       or not 0 <= x < self.qubits for x in q):
            raise ValueError('Invalid or repeated gate wires.')
        self.gates.append((op,*map(int,q)))
    def ccx(self,a,b,t):
        for op,qs in [('h',(t,)),('cx',(b,t)),('tdg',(t,)),('cx',(a,t)),('t',(t,)),('cx',(b,t)),('tdg',(t,)),('cx',(a,t)),('t',(b,)),('t',(t,)),('h',(t,)),('cx',(a,b)),('t',(a,)),('tdg',(b,)),('cx',(a,b))]:self.gate(op,*qs)
    def fredkin(self,c,a,b):self.gate('cx',a,b);self.ccx(c,b,a);self.gate('cx',a,b)
    def literal_update(self,addr,val,t,other=None):
        if not val:self.gate('x',addr)
        if other is None:self.gate('cx',addr,t)
        else:self.ccx(other,addr,t)
        if not val:self.gate('x',addr)
    def loader(self,address,f):
        low=address[:self.ell]; high=list(reversed(address[self.ell:])); k=len(high)
        if not k:
            for j,v in enumerate(f):
                if v:self.gate('x',self.banks[j])
            return
        def traversal(skip):
            def walk(d,prefix):
                for val in (0,1):
                    if d or not skip:self.literal_update(high[d],val,self.sel[d],self.sel[d-1] if d else None)
                    now=2*prefix+val
                    if d==k-1:
                        for j in range(self.lam):
                            if f[now*self.lam+j]:self.gate('cx',self.sel[-1],self.banks[j])
                    else:walk(d+1,now)
                    if d or not skip:self.literal_update(high[d],val,self.sel[d],self.sel[d-1] if d else None)
            walk(0,0)
        traversal(False)
        at=len(self.gates);traversal(True); temp=self.gates[at:];self.gates[at:]=inverse(temp)
    def router(self,address):
        for d in reversed(range(self.ell)):
            for j in range(1<<d):self.fredkin(address[d],self.banks[j],self.banks[j+(1<<d)])
    def controlled(self,G,c,t):
        if G=='x':self.gate('cx',c,t)
        elif G=='z':self.gate('h',t);self.gate('cx',c,t);self.gate('h',t)
        elif G=='h':
            # V†, CZ, V, with V=S H T H S† and V Z V†=H.
            for op in ('sdg','h','tdg','h','s'):self.gate(op,t)
            self.controlled('z',c,t)
            for op in ('sdg','h','t','h','s'):self.gate(op,t)
        else:raise ValueError(G)
    def mux(self,t,f,G='x'):
        t = _integer(t, 'target')
        f = tuple(f)
        if t >= self.w or len(f) != 1 << self.r or any(v not in (0,1) for v in f):
            raise ValueError('Boolean table or target is invalid.')
        if G not in ('x','z','h'):
            raise ValueError('Only involutions X,Z,H are supported.')
        addr=[q for q in range(self.w) if q!=t]
        at=len(self.gates);self.loader(addr,f);L=self.gates[at:]
        at=len(self.gates);self.router(addr);R=self.gates[at:]
        self.controlled(G,self.banks[0],t);self.gates.extend(inverse(R));self.gates.extend(inverse(L))
        self.gates.extend(R);self.controlled(G,self.banks[0],t);self.gates.extend(inverse(R))
    def permutation(self,p):
        for t,f in benes(p):self.mux(t,f)

def inverse(gates):
    inv={'t':'tdg','tdg':'t','s':'sdg','sdg':'s'}
    return [(inv.get(g[0],g[0]),*g[1:]) for g in reversed(gates)]

def emit_readout(n: int, bands: Sequence[Band] | None = None):
    """Compile the complete readout, not just its initialized columns."""
    geom = history_geometry(n, bands)
    M = 4 * geom.dimension
    em = Emitter(n + 2)
    for pairs in geom.matches:
        p = [-1] * M
        for j, (a, b) in enumerate(pairs):
            p[a], p[b] = 2*j, 2*j+1
        free = iter(range(2*len(pairs), M))
        for i in range(M):
            if p[i] < 0:
                p[i] = next(free)
        em.permutation(p)
        table = [int(j < len(pairs)) for j in range(M//2)]
        em.mux(0, table, 'z')
        em.mux(0, table, 'h')
        em.permutation(np.argsort(p))
    return em

def simulate(gates,v):
    v=np.asarray(v,complex).copy(); one=v.ndim==1
    if v.ndim not in (1,2) or len(v) == 0 or len(v) & (len(v)-1):
        raise ValueError('State arrays must have a positive power-of-two leading dimension.')
    if not np.all(np.isfinite(v)):
        raise ValueError('State entries must be finite.')
    if one:v=v[:,None]
    d=len(v);ids=np.arange(d);cache={};omega=np.exp(1j*np.pi/4)
    for raw_gate in gates:
        g = tuple(raw_gate)
        if g[0] not in ('h','s','sdg','t','tdg','x','cx') or len(g) != (3 if g[0]=='cx' else 2):
            raise ValueError('Invalid native gate.')
        if any(not isinstance(q, Integral) or not 0 <= q < d.bit_length()-1 for q in g[1:]):
            raise ValueError('Gate outside the state register.')
        if g not in cache:
            op=g[0];t=g[-1]; mask=(ids>>t)&1
            if op in ('x','cx'):
                yes=np.ones(d,dtype=bool) if op=='x' else ((ids>>g[1])&1).astype(bool)
                cache[g]=ids^(yes.astype(int)<<t)
            elif op=='h':cache[g]=ids[mask==0]
            else:cache[g]=ids[mask==1]
        index=cache[g];op=g[0]
        if op in ('x','cx'):v=v[index]
        elif op=='h':
            other=index|(1<<g[1]);a=v[index].copy();b=v[other].copy()
            v[index]=(a+b)/math.sqrt(2);v[other]=(a-b)/math.sqrt(2)
        else:v[index]*={'s':1j,'sdg':-1j,'t':omega,'tdg':omega.conjugate()}[op]
    return v[:,0] if one else v


H=np.array([[1,1],[1,-1]],complex)/math.sqrt(2)
X=np.array([[0,1],[1,0]],complex)
Z=np.diag([1,-1]).astype(complex)
I2=np.eye(2,dtype=complex)
OMEGA=np.exp(1j*np.pi/4)

def refs_from_native(word):
    """Chronological H,T,Tdg,X,Z -> chronological R_j, Z (literal equality)."""
    expanded=[]
    for op in word:
        if op=='x': expanded.extend(['h']+['t']*4+['h'])
        elif op=='z': expanded.extend(['t']*4)
        else: expanded.append(op)
    power=0; left=[]
    for op in reversed(expanded):
        if op=='h':left.append(power%8)
        elif op=='t':power+=1
        elif op=='tdg':power-=1
        else:raise ValueError(op)
    if power%8 not in (0,4):raise ValueError('word has determinant other than +/-1')
    if power%8==4:left.append('z')
    return list(reversed(left))

def ref_matrix(ref):
    if ref=='z':return Z
    T=np.diag([1,OMEGA**int(ref)])
    return T@H@T.conj().T

def ref_word_matrix(refs):
    ans=I2.copy()
    for g in refs:ans=ref_matrix(g)@ans
    return ans

@dataclass(frozen=True)
class EchoProgram:
    """A literal D=(X P† X P)^2 word; its echo identity is syntactic."""
    seed: tuple[str, ...]

    def __post_init__(self):
        seed = tuple(self.seed)
        if any(op not in ('h','t','tdg','x','z') for op in seed):
            raise ValueError('Seed gates must be H,T,Tdg,X,Z.')
        object.__setattr__(self, 'seed', seed)

    @property
    def reflections(self):
        inv = {'t':'tdg','tdg':'t','h':'h','x':'x','z':'z'}
        pinv = [inv[g] for g in reversed(self.seed)]
        return tuple(refs_from_native((list(self.seed)+['x']+pinv+['x'])*2))

    def matrix(self):
        return ref_word_matrix(self.reflections)


class CompactEmitter(Emitter):
    def __init__(self,w):
        super().__init__(w)
        self.z=super().qubits
    @property
    def qubits(self):return self.z+1
    def mcx(self,controls,target,borrow):
        controls=list(controls)
        if not controls:self.gate('x',target);return
        if len(controls)==1:self.gate('cx',controls[0],target);return
        if len(controls)==2:self.ccx(*controls,target);return
        a=controls[:(len(controls)+1)//2];b=controls[len(a):]
        self.mcx(a,borrow,target)
        self.mcx([borrow]+b,target,a[0])
        self.mcx(a,borrow,target)
        self.mcx([borrow]+b,target,a[0])
    def mcx_zero(self,controls,target,borrow):
        for c in controls:self.gate('x',c)
        self.mcx(controls,target,borrow)
        for c in reversed(controls):self.gate('x',c)
    def cz_many(self,controls,t):
        self.gate('h',t)
        if len(controls)==1:self.gate('cx',controls[0],t)
        elif len(controls)==2:self.ccx(*controls,t)
        else:raise ValueError('one or two controls required')
        self.gate('h',t)
    def controlled_ref(self,ref,controls,t):
        if ref=='z':self.cz_many(controls,t);return
        j=int(ref)%8
        # T^-j, V†, CCZ/CZ, V, T^j; V=S H T H S† maps Z to H.
        for _ in range(j):self.gate('tdg',t)
        for op in ('sdg','h','tdg','h','s'):self.gate(op,t)
        self.cz_many(controls,t)
        for op in ('sdg','h','t','h','s'):self.gate(op,t)
        for _ in range(j):self.gate('t',t)
    def table_ref(self,address,target,table,ref,control=None):
        address=list(address); k=len(address)
        if len(table)!=(1<<k):raise ValueError('table length mismatch')
        if not any(table):return
        if not k:
            if control is None:
                if ref=='z':self.gate('h',target);self.gate('x',target);self.gate('h',target)
                else:
                    for _ in range(int(ref)):self.gate('tdg',target)
                    self.gate('h',target)
                    for _ in range(int(ref)):self.gate('t',target)
            else:self.controlled_ref(ref,[control],target)
            return
        # Temporarily choose banks/selectors for the smaller address table.
        old=(self.ell,self.lam,self.banks,self.sel)
        self.ell=k//2;self.lam=1<<self.ell
        self.banks=old[2][:self.lam];self.sel=old[3][:k-self.ell]
        at=len(self.gates);self.loader(address,table);L=self.gates[at:]
        at=len(self.gates);self.router(address);R=self.gates[at:]
        ctrls=[self.banks[0]]+([] if control is None else [control])
        self.controlled_ref(ref,ctrls,target)
        self.gates.extend(inverse(R));self.gates.extend(inverse(L))
        self.gates.extend(R);self.controlled_ref(ref,ctrls,target)
        self.gates.extend(inverse(R))
        self.ell,self.lam,self.banks,self.sel=old
    def word_mux(self,address,target,words,control=None):
        alphabet=list(range(8))+['z']
        length=max(map(len,words),default=0)
        for pos in range(length):
            for ref in alphabet:
                table=[int(pos<len(w) and w[pos]==ref) for w in words]
                self.table_ref(address,target,table,ref,control)
    def compressed_matching(self,pairs,programs):
        pairs, programs = tuple(pairs), tuple(programs)
        if len(pairs) != len(programs) or any(not isinstance(p, EchoProgram) for p in programs):
            raise ValueError('Every pair requires a syntactically certified EchoProgram.')
        labels = [v for pair in pairs for v in pair]
        if any(len(pair) != 2 for pair in pairs) or len(set(labels)) != len(labels):
            raise ValueError('Pairs must be disjoint.')
        if any(not isinstance(v, Integral) or isinstance(v, bool) or not 0 <= v < (1<<self.w) for v in labels):
            raise ValueError('Invalid matching label.')
        words = [p.reflections for p in programs]
        M=1<<self.w;count=len(pairs)
        if not count:return
        size=1<<(count-1).bit_length();k=size.bit_length()-1
        if 2*size>M:raise ValueError('matching too large')
        p=[-1]*M
        for j,(a,b) in enumerate(pairs):p[a]=2*j;p[b]=2*j+1
        free=iter(range(2*count,M))
        for a in range(M):
            if p[a]<0:p[a]=next(free)
        self.permutation(p)
        tab=list(words)+[[] for _ in range(size-count)]
        addr=list(range(1,k+1));high=list(range(k+1,self.w))
        if not high:self.word_mux(addr,0,tab)
        else:
            # Y_z (T_h G_z T_h) G_z† Y_z, read rightmost first.
            at=len(self.gates);self.word_mux(addr,0,tab,self.z)
            G=self.gates[at:];del self.gates[at:]
            self.gate('cx',self.z,0)
            self.gates.extend(inverse(G))
            self.mcx_zero(high,self.z,0)
            self.gates.extend(G)
            self.mcx_zero(high,self.z,0)
            self.gate('cx',self.z,0)
        self.permutation(np.argsort(p))

def emit_prefix_bank(n: int, programs: Sequence[Sequence[EchoProgram]]):
    """Compile actual coarse prefixes with exact logical/dirty separation."""
    n = _integer(n, 'n', 2)
    programs = tuple(tuple(row) for row in programs)
    if len(programs) != n-1 or any(len(row) != 1<<d for d,row in enumerate(programs)):
        raise ValueError('Supply rows of 2**d EchoPrograms for d < n-1.')
    if any(not isinstance(p, EchoProgram) for row in programs for p in row):
        raise ValueError('Prefix blocks require EchoProgram certificates.')
    em=CompactEmitter(n+1)
    for ell,pairs in enumerate(prefix_matchings(n)):
        em.compressed_matching([(a,b) for a,b,p in pairs],
                               [programs[ell][p] for a,b,p in pairs])
    return em

