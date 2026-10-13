"""Bounded diagnostics for source-program resynthesis, not a fine compiler.

Big-endian wires. No initialized inputs are used. Matrices here diagnose
explicit formulae; arbitrary rotations are not claimed to be native gates.
The native costs and interface restrictions are proved in
research/endpoint/SOURCE_PROGRAM_RESYNTHESIS.md. Exact ranks use rational
arithmetic; floating QR and matrix routines diagnose identities only.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction
import math
import numpy as np
from numbers import Integral

I=np.eye(2,dtype=complex)
X=np.array([[0,1],[1,0]],complex)
Y=np.array([[0,-1j],[1j,0]],complex)
Z=np.diag([1.,-1.]).astype(complex)
H=(X+Z)/math.sqrt(2)
S=np.diag([1,1j]); T=np.diag([1,np.exp(1j*np.pi/4)])
ONE={'H':H,'X':X,'S':S,'Sdg':S.conj().T,'T':T,'Tdg':T.conj().T}

def integer(value, name, minimum=0, maximum=None):
    if (isinstance(value, bool) or not isinstance(value, Integral)
            or value < minimum or (maximum is not None and value > maximum)):
        raise ValueError(f'{name} is outside its integer domain.')
    return int(value)


def validate_masks(masks, m):
    integer(m, 'source width', 2, 4096)
    raw = np.asarray(masks)
    if (raw.ndim != 2 or raw.shape[1] != m or raw.shape[0] == 0
            or not np.isin(raw, (0, 1)).all()):
        raise ValueError('Nonempty binary mask rows of the source width are required.')
    return raw.astype(np.uint8)


def exact_pivots(rows):
    """Pivot columns over Q; never infer source rank from floating QR."""
    a = [[Fraction(x) for x in row] for row in rows]
    if not a or not a[0] or any(len(row) != len(a[0]) for row in a):
        raise ValueError('A nonempty rectangular rational matrix is required.')
    rank, pivots = 0, []
    for col in range(len(a[0])):
        pivot = next((r for r in range(rank, len(a)) if a[r][col]), None)
        if pivot is None:
            continue
        a[rank], a[pivot] = a[pivot], a[rank]
        scale = a[rank][col]
        a[rank] = [x / scale for x in a[rank]]
        for r in range(rank + 1, len(a)):
            factor = a[r][col]
            if factor:
                a[r] = [x - factor*y for x, y in zip(a[r], a[rank])]
        pivots.append(col)
        rank += 1
        if rank == len(a):
            break
    return tuple(pivots)


def sign_mask(coefficient, m):
    """Canonical rational rounding, with literal endpoint/zero encodings."""
    integer(m, 'source width', 2, 4096)
    c = Fraction(coefficient)
    if not -1 <= c <= 1:
        raise ValueError('Coefficient must be in [-1,1].')
    if c == 1:
        return np.zeros(m, dtype=np.uint8)
    if c == 0:
        result = np.zeros(m, dtype=np.uint8); result[0] = 1
        return result
    if c == -1:
        return np.ones(m, dtype=np.uint8)
    x = (1-c) * (1 << (m-2))
    k = (2*x.numerator + x.denominator) // (2*x.denominator)
    if k == 1 << (m-1):
        return np.ones(m, dtype=np.uint8)
    return np.array([(k >> (m-2-j)) & 1 for j in range(m-1)] + [0], dtype=np.uint8)


def encoded_coefficient(mask):
    bits = validate_masks([mask], len(mask))[0]
    return sum((w * (1-2*int(f)) for w, f in zip(weights(len(bits)), bits)), Fraction(0))


def actual_table_family(m, h, v=3):
    integer(v, 'leading precision', 3)
    integer(h, 'mask family size', 1)
    integer(m, 'source width', v+h+2, 4096)
    coefficients = [Fraction(1, 1 << v)]
    coefficients += [coefficients[0]-Fraction(1, 1 << (v+j)) for j in range(1,h+1)]
    return coefficients, np.array([sign_mask(c,m) for c in coefficients])


def chebyshev_fraction(c, order):
    c = Fraction(c); integer(order, 'Chebyshev order', 0)
    previous, current = Fraction(1), c
    if order == 0:
        return previous
    for _ in range(1, order):
        previous, current = current, 2*c*current-previous
    return current


def scalar_power_lower_bound(m, repetitions):
    integer(m, 'source width', 4); integer(repetitions, 'repetitions', 1)
    return 4*repetitions*(m-3)+2


def single_mask_alignment(m):
    """Ideal alignment, no claim that its arbitrary angles are native gates."""
    integer(m, 'single-mask alignment width', 4, 8)
    gs = majoranas(m); a = 2**(-(m-1)/2)
    source = sum(float(w)**0.5*g for w,g in zip(weights(m),gs))
    beta = 0.5*math.atan2(math.sqrt(1-2*a*a), -a*a)
    phi = -0.5*math.atan2(a, math.sqrt(1-2*a*a))
    return exp_axis(phi,gs[-2]) @ exp_axis(beta,source), 2*math.asin(a)


def symplectic(a,b):
    return ((a[0]&b[1]).bit_count() + (a[1]&b[0]).bit_count()) & 1


def transverse_labels(p):
    """Independent anticommuting labels, retaining odd-p central parity."""
    integer(p, 'transverse rank', 1, 12)
    width = (p+1)//2
    labels = []
    for j in range(width):
        labels.extend([(1 << j, (1 << j)-1), (1 << j, (1 << (j+1))-1)])
    return tuple(labels[:p]), width


def lie_labels(p, external_qubits=1):
    """Exact binary commutator closure; a representation diagnostic, not T cost."""
    integer(p, 'transverse rank', 1, 10)
    integer(external_qubits, 'external width', 1, 2)
    ps, _ = transverse_labels(p)
    q = external_qubits; limit = 1 << q
    local = [(x,z) for x in range(limit) for z in range(limit) if x or z]
    labels = set(local)
    labels.update(((x << q)^a, (z << q)^b) for x,z in ps for a,b in local)
    queue = list(labels)
    # Closing under the original generators suffices by Jacobi's identity.
    generators = tuple(queue)
    for current in queue:
        for other in generators:
            if symplectic(current,other):
                new = current[0]^other[0], current[1]^other[1]
                if new not in labels:
                    labels.add(new); queue.append(new)
    return frozenset(labels)


def norm(a):return float(np.linalg.norm(a,2))
def kronall(ops):
    out=np.ones((1,1),complex)
    for a in ops:out=np.kron(out,a)
    return out

def majoranas(m):
    m = integer(m, 'dense core width', 2, 8)
    return [kronall([Z]*j+[X]+[I]*(m-j-1)) for j in range(m)]

def weights(m):
    m = integer(m, 'source width', 2, 4096)
    return [Fraction(1,1 << (j+1)) for j in range(m-1)]+[Fraction(1,1 << (m-1))]

def source_data(m,masks):
    masks=validate_masks(masks,m)
    gs=majoranas(m);a=np.sqrt(np.array([float(t) for t in weights(m)]))
    v=[a]+[a*(1-2*f.astype(int)) for f in masks]
    matrices=[sum(c*g for c,g in zip(row,gs)) for row in v]
    groups={}
    for j in range(m):groups.setdefault(tuple(int(f[j]) for f in masks),[]).append(j)
    groups=list(groups.values())
    signmatrix=[[1]+[1-2*int(masks[t,j]) for t in range(len(masks))] for j in range(m)]
    return gs,a,v,matrices,groups,len(exact_pivots(signmatrix))

def exp_axis(t,A):return math.cos(t)*np.eye(len(A))+1j*math.sin(t)*A

def rotation_parameters(c,half):
    if c < -1e-10 or c > 1+1e-10:raise ValueError('overlap must be in [0,1]')
    c=float(np.clip(c,0,1));t=math.sqrt((1+c)/2)
    r=1/(4*t*math.sqrt((1+t)/2)) if half else 1/(2*t)
    beta=math.asin(r);zeta=math.atan2(c*r,math.sqrt(1-r*r))
    return math.pi/4-zeta/2,beta,-math.pi/4-zeta/2

def move(A,B,c,half=False):
    """Full move A -> B, or half move A -> normalized(A+B).

The returned product uses three exponentials of A/B, with a continuous
formula at c=1. For Clifford vectors it is an even plane rotation.
"""
    alpha,beta,gamma=rotation_parameters(c,half)
    return exp_axis(alpha,A)@exp_axis(beta,B)@exp_axis(gamma,A)

def class_basis(gs,a,groups):
    """Unitary B_class maps each first-coordinate Gamma to its class vector."""
    m=len(gs);size=len(gs[0]);U=np.eye(size,dtype=complex);errors=[];details=[]
    for ids in groups:
        p=ids[0];mass=float(np.sum(a[ids]**2));v=np.zeros(m);v[ids]=a[ids]/math.sqrt(mass)
        V=sum(c*g for c,g in zip(v,gs))
        tailvec=np.zeros(m);tailvec[p:]=a[p:]/math.sqrt(float(np.sum(a[p:]**2)))
        A=sum(c*g for c,g in zip(tailvec,gs))
        signs=-np.ones(m);signs[ids]=1
        N=sum(c*g for c,g in zip(tailvec*signs,gs))
        c=float(tailvec@(tailvec*signs))
        G=move(A,N,c,half=True)
        errors.append(norm(G@A@G.conj().T-V))
        # This *whole* second rotor is supported in the class Clifford span.
        C=move(gs[p],V,float(v[p]),half=False)
        errors.append(norm(C@gs[p]@C.conj().T-V))
        for j in range(m):
            if j not in ids:errors.append(norm(C@gs[j]@C.conj().T-gs[j]))
        U=C@U
        details.append({'indices':ids,'first':p,'mass':mass,'tail_overlap':c,
                        'endpoint_overlap':float(v[p])})
    for ids in groups:
        v=np.zeros(m);v[ids]=a[ids]/np.linalg.norm(a[ids])
        errors.append(norm(U@gs[ids[0]]@U.conj().T-sum(c*g for c,g in zip(v,gs))))
    return U,max(errors,default=0),details

def pauli_apply(word,n,columns=None):
    dim=1<<n;out=np.eye(dim,dtype=complex) if columns is None else np.array(columns,complex,copy=True)
    was_vector=out.ndim==1
    if was_vector:out=out[:,None]
    ids=np.arange(dim)
    for gate in word:
        if gate[0]=='CX':
            _,c,t=gate
            perm=ids ^ (((ids>>(n-1-c))&1) << (n-1-t));out=out[perm]
        else:
            name,q=gate;mat=ONE[name];bit=1<<(n-1-q)
            lo=ids[(ids&bit)==0];hi=lo|bit
            a=out[lo].copy();b=out[hi].copy()
            out[lo]=mat[0,0]*a+mat[0,1]*b
            out[hi]=mat[1,0]*a+mat[1,1]*b
    return out[:,0] if was_vector else out

def compress_paulis(m,generators):
    """Binary symplectic elimination with an emitted H/S/CNOT circuit.

Rows are binary [x|z]; phase signs are carried by the actual circuit, not
used in pivot selection. Returns circuit and number of occupied output bits.
"""
    rows=[np.array(r,dtype=np.uint8).copy() for r in generators]
    originals=[r.copy() for r in rows];word=[];slot=0
    def gate(g):
        word.append(g)
        for r in rows+originals:
            if g[0]=='H':
                q=g[1];r[q],r[m+q]=r[m+q],r[q]
            elif g[0] in ('S','Sdg'):
                q=g[1];r[m+q]^=r[q]
            elif g[0]=='CX':
                _,c,t=g;r[t]^=r[c];r[m+c]^=r[m+t]
            else:raise ValueError(g)
    def swap(a,b):
        if a!=b:
            gate(('CX',a,b));gate(('CX',b,a));gate(('CX',a,b))
    def to_z(index):
        A=rows[index]
        active=[j for j in range(slot,m) if A[j] or A[m+j]]
        if not active:raise AssertionError('zero pivot')
        for j in active:
            if A[j] and A[m+j]:gate(('Sdg',j));gate(('H',j))
            elif A[j]:gate(('H',j))
        pivot=next(j for j in range(slot,m) if A[m+j])
        swap(pivot,slot)
        for j in range(slot+1,m):
            if A[m+j]:gate(('CX',j,slot))
    while rows:
        pair=None
        for i in range(len(rows)):
            for j in range(i+1,len(rows)):
                omega=(int(rows[i][:m]@rows[j][m:])+int(rows[i][m:]@rows[j][:m]))%2
                if omega:pair=(i,j);break
            if pair is not None:break
        if pair is None:
            to_z(0);A=rows[0].copy()
            for r in rows[1:]:
                if r[m+slot]:r^=A
            rows=rows[1:];slot+=1
        else:
            i,j=pair
            rows[0],rows[i]=rows[i],rows[0]
            if j==0:j=i
            rows[1],rows[j]=rows[j],rows[1]
            to_z(0)
            B=rows[1]
            assert B[slot]
            if B[m+slot]:gate(('Sdg',slot))
            for k in range(slot+1,m):
                x,z=int(B[k]),int(B[m+k])
                if (x,z)==(1,0):gate(('CX',slot,k))
                elif (x,z)==(0,1):
                    gate(('H',k));gate(('CX',slot,k));gate(('H',k))
                elif (x,z)==(1,1):
                    gate(('Sdg',k));gate(('CX',slot,k));gate(('S',k))
            A=rows[0].copy();B=rows[1].copy()
            for r in rows[2:]:
                if r[slot]:r^=B
                if r[m+slot]:r^=A
            rows=rows[2:];slot+=1
        rows=[r for r in rows if np.any(r)]
        assert all(not np.any(r[:slot]) and not np.any(r[m:m+slot]) for r in rows)
    assert all(not np.any(r[slot:m]) and not np.any(r[m+slot:]) for r in originals)
    return word,slot

def coefficient_alignment(vectors):
    """Orthogonal coefficient rotation through QR Givens; finite diagnostic.

Exact rank and independence in the proof are determined from the integer
sign matrix, not from this floating QR.
"""
    mat=np.column_stack(vectors)
    q=mat.shape[0]
    Q,_=np.linalg.qr(mat,mode='complete')
    if np.linalg.det(Q)<0:Q[:,-1]*=-1
    work=Q.copy();rotations=[]
    for col in range(q-1):
        for row in range(q-1,col,-1):
            a=float(work[col,col]);b=float(work[row,col])
            if abs(b)<2e-15:continue
            r=math.hypot(a,b);c=a/r;s=b/r
            old=work[[col,row],:].copy()
            work[col,:]=c*old[0]+s*old[1]
            work[row,:]=-s*old[0]+c*old[1]
            rotations.append((col,row,math.atan2(s,c)))
    # Zero subdiagonal entries may leave paired negative diagonal pivots.
    negative=[j for j in range(q) if work[j,j]<0]
    assert len(negative)%2==0
    for i,j in zip(negative[::2],negative[1::2]):
        work[i,:]*=-1;work[j,:]*=-1
        rotations.append((i,j,math.pi))
    assert norm(work-np.eye(q))<1e-10
    return rotations,Q.T

@dataclass
class Model:
    m:int
    masks:np.ndarray
    d:int
    q:int
    active:int
    decoder:np.ndarray
    scalars:list[np.ndarray]
    compact_scalars:list[np.ndarray]
    max_error:float
    class_details:list[dict]
    clifford_word:list[tuple]

def construct_model(m,masks):
    integer(m, 'dense model width', 2, 8)
    gs,a,vs,Ms,groups,d=source_data(m,masks)
    # A_class takes chosen Gamma coordinates to disjoint class vectors.
    A_class,err,details=class_basis(gs,a,groups)
    reps=[g[0] for g in groups];q=len(reps);dim=1<<m
    Cperm=np.eye(dim,dtype=complex);images=list(range(m));signs=np.ones(m)
    for i,orig in enumerate(reps):
        p=images[orig]
        if p==i:continue
        R=(np.eye(dim)+gs[i]@gs[p])/math.sqrt(2)
        Cperm=R@Cperm
        for j in range(m):
            if images[j]==p:images[j]=i
            elif images[j]==i:images[j]=p;signs[j]*=-1
    # Columns of these q-dimensional vectors are exact algebraic in proof.
    smallvecs=[]
    for v in vs:
        sm=np.zeros(q)
        for i,ids in enumerate(groups):sm[i]=signs[reps[i]]*float(np.dot(v[ids],a[ids]))/np.linalg.norm(a[ids])
        smallvecs.append(sm)
    # Select an independent generating list by exact sign rows.
    Sgn=[[1]+[1-2*int(f[ids[0]]) for f in masks] for ids in groups]
    pivots=list(exact_pivots(Sgn))
    rotations,O=coefficient_alignment([smallvecs[i] for i in pivots])
    Crot=np.eye(dim,dtype=complex)
    for i,j,theta in rotations:
        R=math.cos(theta/2)*np.eye(dim)+math.sin(theta/2)*(gs[i]@gs[j])
        Crot=R@Crot
    Cvec=Crot@Cperm@A_class.conj().T
    for v,M in zip(smallvecs,Ms):
        small=O@v
        expected=sum(small[j]*gs[j] for j in range(q))
        err=max(err,norm(Cvec@M@Cvec.conj().T-expected),float(np.linalg.norm(small[d:])))
    # The EVEN algebra is generated by i Gamma_0 Gamma_j, j=1,...,d-1.
    labels=[]
    for j in range(1,d):
        r=np.zeros(2*m,dtype=np.uint8);r[0]=r[j]=1;r[m:m+j]=1;labels.append(r)
    word,active=compress_paulis(m,labels)
    assert active==d//2
    Ccl=pauli_apply(word,m);C=Ccl@Cvec
    scalars=[];compact=[];spec=1<<(m-active)
    for M in Ms[1:]:
        c=float(np.trace(Ms[0]@M).real/dim)
        D=(Ms[0]@M-M@Ms[0])/2
        Xc=C@D@C.conj().T
        shaped=Xc.reshape(1<<active,spec,1<<active,spec)
        Dc=np.trace(shaped,axis1=1,axis2=3)/spec
        err=max(err,norm(Xc-np.kron(Dc,np.eye(spec))))
        scalars.append(c*np.eye(2*dim)+np.kron(D,X))
        compact.append(c*np.eye(2*(1<<active))+np.kron(Dc,X))
    return Model(m,np.array(masks),d,q,active,C,scalars,compact,err,details,word)

def addressed_scalar(model,compact=False,swap=False):
    """Core, address, signal. One address chooses between the first two masks."""
    if len(model.scalars)<2:raise ValueError('two masks required')
    src=model.compact_scalars if compact else model.scalars
    a,b=(src[1],src[0]) if swap else (src[0],src[1])
    dc=a.shape[0]//2;out=np.zeros((4*dc,4*dc),complex)
    for address,S0 in enumerate((a,b)):
        ids=np.array([4*i+2*address+signal for i in range(dc) for signal in range(2)])
        out[np.ix_(ids,ids)]=S0
    return out

def embed_compact(V,m,active,logical=2):
    """Embed [active core, logical] between active core and core spectators."""
    ds=1<<(m-active);dl=1<<logical;da=1<<active
    out=np.kron(V,np.eye(ds)).reshape(da,dl,ds,da,dl,ds).transpose(0,2,1,3,5,4)
    return out.reshape((da*ds*dl,)*2)

def program_word(model,steps,seed=0):
    """Address changes include H, S and signal-to-address CNOT; no fixed row."""
    rng=np.random.default_rng(seed);dim=1<<(model.m+2);smalldim=1<<(model.active+2)
    V=np.eye(dim,dtype=complex);Vsmall=np.eye(smalldim,dtype=complex)
    sources=[addressed_scalar(model,swap=t) for t in (False,True)]
    smallsources=[addressed_scalar(model,compact=True,swap=t) for t in (False,True)]
    cliffwords=[[('H',0)],[('S',0)],[('H',1)],[('S',1)],[('CX',1,0)],[('CX',0,1)]]
    logical=[pauli_apply(w,2) for w in cliffwords];choices=[]
    for j in range(steps):
        choice=int(rng.integers(2));inv=bool(rng.integers(2));g=int(rng.integers(len(logical)))
        S0=sources[choice].conj().T if inv else sources[choice]
        S1=smallsources[choice].conj().T if inv else smallsources[choice]
        V=np.kron(np.eye(1<<model.m),logical[g])@S0@V
        Vsmall=np.kron(np.eye(1<<model.active),logical[g])@S1@Vsmall
        choices.append((choice,inv,g))
    Cbig=np.kron(model.decoder,np.eye(4))
    expected=Cbig.conj().T@embed_compact(Vsmall,model.m,model.active)@Cbig
    return V,Vsmall,norm(V-expected),choices

