"""Exact bit-sliced Clifford tableaux for every augmented-helper character.

A bit position denotes a complete helper-X character. This bounded evaluator
is not a native emitter or a scalable synthesis algorithm.
Signed Hermitian Pauli images are exact; tableau equality alone quotients a
unitary's global scalar, so literal scalar claims use the analytic cell proof
and independent bounded matrix calibrations.
"""
from functools import lru_cache, reduce
from tests.test_recursive_field_shear_decoder import (_tower_product as mul,
    _tower_trace as trace,_address_field)
from tests.test_projective_field_recoupler import fourier_word,linear_word

@lru_cache(None)
def phase_coefficients(level,kind):
    ell=1<<(level-2); k=2*ell; width=2*k
    def field_value(basis):
        mask=(1<<ell)-1;u=basis&((1<<k)-1);v=basis>>k
        d0,d1=u&mask,u>>ell;a0,a1=v&mask,v>>ell
        if kind=='before':return (mul(level,u,v)>>ell)^mul(level-1,d0,d1)
        return mul(level-1,d0,a1)^mul(level-1,d1,a0)
    linear=[field_value(1<<i) for i in range(width)]
    out=[]
    for i,c in enumerate(linear):
        if c:
            mask=sum(trace(level-1,mul(level-1,1<<b,c))<<b for b in range(ell))
            out.append(('ZC',i,mask))
    for i in range(width):
        for j in range(i+1,width):
            c=field_value((1<<i)|(1<<j))^linear[i]^linear[j]
            if c:
                mask=sum(trace(level-1,mul(level-1,1<<b,c))<<b for b in range(ell))
                out.append(('CZC',i,j,mask))
    return tuple(out)

@lru_cache(None)
def recoupler_word(level,offset=0):
    ell=1<<(level-2);k=2*ell;theta=1<<ell
    def phase(kind):
        return [(g[0],*(q+offset for q in g[1:-1]),g[-1]<<(ell-1)) for g in phase_coefficients(level,kind)]
    htheta=linear_word([mul(level,theta,1<<j) for j in range(k)],range(offset+k,offset+2*k))+fourier_word(level,offset+k)
    wall=[g for j in range(k) for g in [('CX',offset+j,offset+k+j),('CX',offset+k+j,offset+j)]]
    return tuple(fourier_word(level,offset+k)+phase('before')+htheta+phase('after')+wall+fourier_word(level-1,offset+ell)+fourier_word(level-1,offset+k+ell))

@lru_cache(None)
def decoder_word(level,offset=0):
    if level==1:return (('CX',offset+1,offset),)
    return recoupler_word(level,offset)+decoder_word(level-1,offset)+decoder_word(level-1,offset+(1<<(level-1)))

class BitslicedTableau:
    def __init__(self,n,h):
        if n < 1 or not 0 <= h <= 15:
            raise ValueError("Bounded fixture requires n >= 1 and 0 <= helper bits <= 15")
        self.n,self.h=n,h;self.sectors=1<<h;self.ones=(1<<self.sectors)-1
        self.x=[[self.ones if c==q else 0 for c in range(2*n)] for q in range(n)]
        self.z=[[self.ones if c==n+q else 0 for c in range(2*n)] for q in range(n)]
        self.sign=[0]*(2*n);self.masks={0:0}
        self.variables=[]
        for bit in range(h):
            step=1<<bit;mask=0
            # Repeated blocks of step zeros and step ones.
            block=(1<<step)-1
            for start in range(step,self.sectors,2*step):mask|=block<<start
            self.variables.append(mask)
    def condition(self,mask,shift=0):
        if mask not in self.masks:
            value=0
            for bit in range(self.h):
                if mask>>bit&1:value^=self.variables[bit]
            self.masks[mask]=value
        return self.masks[mask] ^ (self.ones if (mask&shift).bit_count()%2 else 0)
    def apply(self,word,shift=0):
        cols=range(2*self.n)
        for gate in word:
            kind=gate[0]
            if kind=='H':
                q=gate[1];x,z=self.x[q],self.z[q]
                for c in cols:self.sign[c]^=x[c]&z[c]
                self.x[q],self.z[q]=z,x
            elif kind=='CX':
                a,b=gate[1:];xa,xb,za,zb=self.x[a],self.x[b],self.z[a],self.z[b]
                for c in cols:
                    self.sign[c]^=xa[c]&zb[c]&(xb[c]^za[c]^self.ones)
                    xb[c]^=xa[c];za[c]^=zb[c]
            elif kind=='CZC':
                a,b,mask=gate[1:];cond=self.condition(mask,shift)
                xa,xb,za,zb=self.x[a],self.x[b],self.z[a],self.z[b]
                for c in cols:
                    self.sign[c]^=cond&xa[c]&xb[c]&(za[c]^zb[c])
                    za[c]^=cond&xb[c];zb[c]^=cond&xa[c]
            elif kind=='ZC':
                q,mask=gate[1:];cond=self.condition(mask,shift)
                for c in cols:self.sign[c]^=cond&self.x[q][c]
            else:raise ValueError(gate)
        return self
    def same(self,other):return self.x==other.x and self.z==other.z and self.sign==other.sign
    def at(self,r):
        return tuple((sum(((self.x[q][c]>>r)&1)<<q for q in range(self.n)),sum(((self.z[q][c]>>r)&1)<<q for q in range(self.n)),(self.sign[c]>>r)&1) for c in range(2*self.n))
    def anf(self,value):
        for b,mask in enumerate(self.variables):value^=(value<<(1<<b))&mask
        return value
    def degree_profile(self):
        degree_masks=[0]*(self.h+1)
        for i in range(self.sectors):degree_masks[i.bit_count()]|=1<<i
        profile={'symplectic':{},'sign':{}};witness={}
        for kind,arrays in [('symplectic',self.x+self.z),('sign',[self.sign])]:
            for row,arr in enumerate(arrays):
                for col,value in enumerate(arr):
                    a=self.anf(value)
                    degree=next((d for d in range(self.h,-1,-1) if a&degree_masks[d]),-1)
                    profile[kind][degree]=profile[kind].get(degree,0)+1
                    if kind not in witness or degree>witness[kind]['degree']:
                        monomial=(a&degree_masks[degree]) if degree>=0 else 0
                        low=monomial&-monomial
                        witness[kind]={'degree':degree,'row':row,'column':col,'monomial_helper_mask':low.bit_length()-1 if low else None}
        return {'histograms':profile,'witnesses':witness}

def query_parameter_word(level, alpha, offset=0):
    """Full-field shear parameter; the address subspace is a separate restriction."""
    degree = 1 << (level-1)
    return [('CX', offset+j, offset+degree+i)
            for j in range(degree) for i in range(degree)
            if mul(level, alpha, 1 << j) >> i & 1]

def query_word(level,addr):
    return query_parameter_word(level, _address_field(level,addr))

def target_word(level,addr):
    n=1<<level;out=[]
    for i in range(n):
        j=i^addr
        if i<j:out.extend([('CX',i,j),('CX',j,i),('CX',i,j)])
    return out

def helper_shift(level,addr):
    return sum(_address_field(j-1,addr&((1<<(j-1))-1))<<((1<<(j-2))-1) for j in range(2,level+1))


@lru_cache(None)
def sector_tableau(level):
    """Return shared read-only exact tableaux for all helper-X characters."""
    return BitslicedTableau(1 << level, (1 << (level-1))-1).apply(decoder_word(level))

def materialize_word(word, character):
    """Replace conditional phase gates by literal Clifford gates in one sector."""
    out=[]
    for gate in word:
        if gate[0] == 'CZC':
            _,a,b,mask=gate
            if (mask & character).bit_count() & 1:
                out.extend([('H',b),('CX',a,b),('H',b)])
        elif gate[0] == 'ZC':
            _,q,mask=gate
            if (mask & character).bit_count() & 1:
                out.extend([('S',q),('S',q)])
        else: out.append(gate)
    return out

def xor_translate(table, value, shift):
    """Truth-table translation f(r) -> f(r+shift), without enumerating r."""
    for bit,ones in enumerate(table.variables):
        if shift >> bit & 1:
            zeros=table.ones ^ ones; step=1 << bit
            value=((value & zeros) << step) | ((value & ones) >> step)
    return value

def active_derivative(table, value, level):
    for bit in range(level-1):
        value ^= xor_translate(table,value,helper_shift(level,1 << bit))
    return value

def top_derivative_identity(level):
    """Verify C_r N_top=(J directsum J)M_r in every helper coset, with signs separate.

    C_r is the (p-1)-fold active-helper derivative of the symplectic tableau.
    This exact binary identity witnesses degree >=p-1; it makes no T-depth claim.
    """
    t=sector_tableau(level); n=t.n; rows=t.x+t.z
    derivative=[[active_derivative(t,v,level) for v in row] for row in rows]
    q=BitslicedTableau(n,0).apply(query_word(level,1 << (level-1)))
    qrows=q.x+q.z
    nx=[[qrows[i][j] ^ (i==j) for j in range(2*n)] for i in range(2*n)]
    sums=[]
    for block in (t.x,t.z):
        v=[]
        for col in range(2*n):
            value=0
            for row in block:value ^= row[col]
            v.append(value)
        sums.append(v)
    for row in range(2*n):
        expected=sums[row//n]
        for col in range(2*n):
            actual=0
            for k in range(2*n):
                if nx[k][col]:actual ^= derivative[row][k]
            assert actual == expected[col],(level,row,col)
    # Each all-ones output block remains nonzero in every helper character.
    assert all(reduce(int.__or__,row,0)==t.ones for row in sums)
    return True
