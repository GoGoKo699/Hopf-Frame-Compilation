"""Scratch: returned private-mask quartic batching from the retained R91 word.

No repository edits. The phase helpers are arbitrary and returned.
"""
import sys,json
from pathlib import Path
from collections import Counter
REPO = next(p for p in Path(__file__).resolve().parents if (p / 'WORKSPACE.md').is_file() and (p / 'tests').is_dir())
sys.path.insert(0,str(REPO))
from tests.test_projective_field_recoupler import mul,inv,_raw_linear_cell,ledger,ring_matrix,equal_ring
from tests.test_joint_lookup_pauli_source import _inverse,_bands
import numpy as np

def affine_r91_batch(pairs,helpers):
    """Ftilde(a,b;h): four T targets h, only logical parity CNOT controls.

    A form is (linear_mask,constant). Ftilde differs from (-i)^(ab)Z_h^(ab)
    by scalar i iff both constants are1; these scalars cancel in the commutator.
    """
    word=[]
    # Parity vertices 0,a,a+b,b, thenreturn0.
    vertices=((0,0),(1,0),(1,1),(0,1))
    changes=(0,1,0,1)
    for step,(ia,ib) in enumerate(vertices):
        for (a,b),h in zip(pairs,helpers):
            parity=(ia+ib+ia*a[1]+ib*b[1])&1
            word.append(('TDG' if parity else 'T',h))
        for pair,h in zip(pairs,helpers):
            mask=pair[changes[step]][0]
            word += [('CX',j,h) for j in range(mask.bit_length()) if mask>>j&1]
    return word

def quartic_batch(forms,helpers):
    assert len(forms)==len(helpers) and len(set(helpers))==len(helpers)
    for row in forms:
        assert all((mask>>h)&1==0 for mask,const in row for h in helpers)
    u=affine_r91_batch([(r[0],r[1]) for r in forms],helpers)
    h=[('H',j) for j in helpers]
    v=h+affine_r91_batch([(r[2],r[3]) for r in forms],helpers)+h
    return u+v+_inverse(u)+_inverse(v)

def displacement_forms():
    out=[]
    for direction in (4,5,6,7):
        columns=[mul(3,inv(3,direction),1<<j) for j in range(4)]
        rows=[sum(((column>>i)&1)<<j for j,column in enumerate(columns)) for i in range(4)]
        shift=_raw_linear_cell(direction)[1]
        if not shift:continue
        for scalarbit in (0,1):
            vector=mul(3,1<<scalarbit,shift)
            out.append(((rows[2],1),(rows[3],1),(rows[scalarbit],0),(vector<<4,0)))
    assert len(out)==6
    return out

def displacement_batch_word():
    """CompleteXORdisplacement onU0..3,V4..7,helpers8..13 returned."""
    forms=displacement_forms();helpers=list(range(8,14))
    wall=[('H',j) for j in range(4,8)]
    word=wall+quartic_batch(forms,helpers)+wall
    return word,{'forms':forms,'helpers':helpers,'ledger':ledger(word),
                 'inverse_ledger':ledger(_inverse(word))}

def evalform(form,logical):return ((form[0]&logical).bit_count()&1)^form[1]

def helper_matrix(word,logical,helper):
    """Exact one-helper word at fixedlogical; otherhelperoperations omitted."""
    local=[]
    for name,*wires in word:
        if name=='CX':
            c,t=wires
            assert c<8 and t>=8
            if t==helper and logical>>c&1:
                # NativeX asHSSH.
                local += [('H',0),('S',0),('S',0),('H',0)]
        else:
            assert len(wires)==1 and wires[0]>=8
            if wires[0]==helper:local.append((name,0))
    return ring_matrix(local,1)

def run():
    # All four independent forms/all32columns; retainactualinverse.
    formal=[((1,0),(2,0),(4,0),(8,0))]
    word=quartic_batch(formal,[4])
    expected=np.zeros((4,32,32),dtype=np.int64)
    for x in range(32):expected[0,x,x]=-1 if x&15==15 else 1
    for emitted in(word,_inverse(word)):
        assert equal_ring(ring_matrix(emitted,5),(expected,0))
    print('FORMAL',ledger(word),flush=True)
    forms=displacement_forms();helpers=list(range(8,14))
    batch=quartic_batch(forms,helpers)
    # Everyprincipalbitstring, bothvaluesofeachhelper, fulltensoridentity.
    checks=0
    for logical in range(256):
        for row,helper in zip(forms,helpers):
            phase=-1 if all(evalform(f,logical) for f in row) else 1
            wanted=np.zeros((4,2,2),dtype=np.int64);wanted[0]=phase*np.eye(2,dtype=np.int64)
            for emitted in(batch,_inverse(batch)):
                assert equal_ring(helper_matrix(emitted,logical,helper),(wanted,0))
                checks+=1
    # ExpectedXORdisplacementagreesfor all16U and16FourierV.
    for logical in range(256):
        u=logical&15;v=logical>>4
        expected_shift=_raw_linear_cell(u)[1] if u else 0
        actual=sum(all(evalform(f,logical) for f in row) for row in forms)%2
        assert actual==((expected_shift&v).bit_count()%2)
    t_bands=_bands(batch)
    assert len(t_bands)==16
    assert all(len({g[1] for g in band})==len(band) for band in t_bands)
    assert all(g[1]>=8 for band in t_bands for g in band)
    print('SIXTERM',ledger(batch),'matrix_checks',checks,'forms',forms,flush=True)
    with open(Path(__file__).with_name('batch_receipt.json'),'w') as f:
        json.dump({'forms':forms,'ledger':ledger(batch),'matrix_checks':checks,'private_helpers':6,'inverse_ledger':ledger(_inverse(batch))},f,indent=2)
    return forms,batch
if __name__=='__main__':run()
