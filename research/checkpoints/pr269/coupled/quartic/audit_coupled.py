"""Independent algebraic and full-column audit of the coupled scratch recoupler."""
import sys,json
from pathlib import Path
import numpy as np
REPO = next(p for p in Path(__file__).resolve().parents if (p / 'WORKSPACE.md').is_file() and (p / 'tests').is_dir())
BASE=Path(__file__).resolve().parent.parent
sys.path.insert(0,str(REPO));sys.path.insert(0,str(BASE))
import coupled_projective as emit
from tests.test_recursive_field_shear_decoder import _tower_product as mul,_tower_inverse as inv,_tower_trace as trace,_tower_shear
from tests.test_field_shear_intertwiner import _translation
from tests.test_flat_field_bank_native import full_intertwiner
from tests.exact_matrices import native_columns
from tests.test_joint_lookup_pauli_source import _inverse,_bands
from private_quartic_batch import displacement_batch_word

def pauli(label):
    x,z=label&15,label>>4
    out=np.zeros((16,16),complex)
    for j in range(16):out[j^x,j]=(1,1j,-1,-1j)[(x&z).bit_count()%4]*(-1 if (z&j).bit_count()%2 else 1)
    return out

def sp(a,b):return (((a&15)&(b>>4)).bit_count()+((b&15)&(a>>4)).bit_count())%2

def reflection(a,b,p):
    assert sp(a,b)==0 and a!=b and sp(a,p)==sp(b,p)==1
    A,B,P=pauli(a),pauli(b),pauli(p)
    J=(P-1j*A@P-1j*B@P-A@B@P)/2
    assert np.array_equal(J.conj().T,J) and np.array_equal(J@J,np.eye(16))
    return J

def field_cell(u):
    assert u
    d0,d1=u&3,u>>2
    nu=mul(2,d0,d0)^mul(2,d0,d1)^mul(2,2,mul(2,d1,d1))
    out=np.zeros((16,16))
    for v in range(16):
        uv=mul(3,u,v)
        rhat=mul(3,4,uv)>>2
        s=mul(2,uv>>2,inv(2,nu))
        for r in range(4):
            a=mul(3,u,r|(s<<2))
            out[a,v]=(-1)**trace(2,mul(2,rhat,r))/2
    assert np.array_equal(out.T@out,np.eye(16))
    return out

def field_shift(u):
    if not u:return 0
    d0,d1=u&3,u>>2
    nu=mul(2,d0,d0)^mul(2,d0,d1)^mul(2,2,mul(2,d1,d1))
    s=mul(2,mul(2,d0,d1),inv(2,nu))
    return mul(3,u,s<<2)

def as_gaussian(ring):
    a,d=ring
    assert not np.any(a[1]) and not np.any(a[3])
    return (a[0]+1j*a[2])/(1<<d)

def fourier(level,offset,width=8):
    degree=1<<(level-1);size=1<<width;mask=((1<<degree)-1)<<offset
    out=np.zeros((size,size))
    norm=1<<(degree//2)
    assert degree%2==0
    for col in range(size):
        source=(col&mask)>>offset
        for target in range(1<<degree):
            row=(col&~mask)|(target<<offset)
            out[row,col]=(-1)**trace(level,mul(level,source,target))/norm
    return out

def run():
    JA=reflection(19,50,32)
    JB64=reflection(156,228,64);JB112=reflection(156,228,112)
    JC44=reflection(89,247,44);JC87=reflection(89,247,87)
    C1=field_cell(1)
    # Verify the changedprobe44 and commonA32 against independentfieldmatrices.
    for u,R in ((4,JB112@JA),(5,JA@JB64),(6,JA@JC44),(7,-JC87@JA)):
        assert np.array_equal(R,field_cell(u)@C1.T),u
    expected=np.zeros((256,256),complex)
    phase_cases=[]
    for u in range(16):
        f={d:int((mul(3,inv(3,d),u)>>2)==0) for d in(4,5,6,7)}
        g=int(bool(u&12))
        K=C1.astype(complex)
        if f[5]:K=(-1j*JB64)@K
        if f[6]:K=(-1j*JC44)@K
        if g:K=(1j*JA)@K
        if f[4]:K=(-1j*JB112)@K
        if f[7]:K=(1j*JC87)@K # relative(-i) and literalminus
        inds=np.arange(16)*16+u
        expected[np.ix_(inds,inds)]=K
        if u:assert np.array_equal(K,field_cell(u)),u
        assert np.array_equal(K.conj().T@K,np.eye(16))
        phase_cases.append((u,list(f.values()),g))
    actual=as_gaussian(emit.ring_matrix(emit.coupled_linear_word()))
    assert np.array_equal(actual,expected)
    assert np.array_equal(as_gaussian(emit.ring_matrix(_inverse(emit.coupled_linear_word()))),expected.conj().T)
    # Independentcompleteeight-wirelogicalrecouplerassembly.
    sh=np.zeros((256,256));summation=np.zeros((256,256))
    for col in range(256):
        u,v=col&15,col>>4
        sh[u|((v^field_shift(u))<<4),col]=1
        summation[v|((v^u)<<4),col]=1
    R=fourier(2,6)@fourier(2,2)@summation@sh@expected@fourier(3,4)
    assert np.array_equal(R,as_gaussian(emit.logical_ring()))
    # IndependentCQ2nativeevaluatorforbottom4-wireword; inspectirrationalparts.
    exact=native_columns(4,full_intertwiner())
    E2=np.zeros((16,16),complex)
    for i,row in enumerate(exact):
        for j,q in enumerate(row):
            assert q.real.b==q.imag.b==0
            E2[i,j]=float(q.real.a)+1j*float(q.imag.a)
    E=np.kron(E2,E2)@R
    assert np.array_equal(E.conj().T@E,np.eye(256))
    for addr in range(8):
        source=[_tower_shear(3,addr,j) for j in range(256)]
        target=[_translation(3,addr,j) for j in range(256)]
        assert np.array_equal(E[:,source],E[target,:]),addr
    # Thefull14-wirewordhasTgatesondistinctwiresateveryclaimedband.
    ledgers={}
    for name,word in [('recoupler',emit.recoupler_word()),('decoder',emit.decoder_word())]:
        for orientation,w in [('forward',word),('inverse',_inverse(word))]:
            for band in _bands(w):assert len({g[1] for g in band})==len(band)
            ledgers[name+'_'+orientation]=emit.ledger(w)
    report={'verdict':'PASS','paired_field_identities':4,'linear_full_columns':256,
        'linear_actual_inverse_columns':256,'complete_intertwiner_columns':256*8,
        'independent_assembly':'directfieldcellandPauliGaussianmatrices;CQ2bottomword;allvaluesexactdyadic',
        'helper_return':'3072exact2x2sectorchecksinthequarticbatchreceipt;tensoridentityonall6helpers',
        'ledgers':ledgers,'phase_cases':phase_cases}
    with open(BASE/'quartic'/'coupled_independent_audit.json','w') as f:json.dump(report,f,indent=2)
    print(json.dumps(report,indent=2))
if __name__=='__main__':run()
