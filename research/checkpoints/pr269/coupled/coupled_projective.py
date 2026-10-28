"""Scratch coupled p3 recoupler: 136T/36bands with six returned dirty helpers.

The original zero-extra-helper PR269 word remains a distinct allocation.
This file does not modify that frozen repository construction.
"""
import sys,json
from pathlib import Path
import numpy as np
REPO = next(p for p in Path(__file__).resolve().parents if (p / 'WORKSPACE.md').is_file() and (p / 'tests').is_dir())
sys.path.insert(0,str(REPO))
from tests.test_projective_field_recoupler import (
    paired_t,pauli_diagonalizer,linear_word,_baseline_word,_raw_linear_cell,
    fourier_word,ring_matrix,equal_ring,adjoint_ring,ledger,mul,inv,
    _inverse,_native_x,_relative_toffoli,_tower_shear,_translation)
from tests.test_flat_field_bank_native import full_intertwiner
sys.path.insert(0,str(Path(__file__).parent/'quartic'))
from private_quartic_batch import displacement_batch_word
OUT=Path(__file__).parent

def relative_pair(left,right,probe,controls=(2,3)):
    tau=paired_t(left,right)
    k,pivot=pauli_diagonalizer(probe)
    pword=(_inverse(k)+[('H',pivot)]
           +_relative_toffoli(*controls,pivot)+[('H',pivot)]+k)
    return _inverse(tau)+pword+tau

def selected_pair(direction,pair,negative_phase=False):
    coords=linear_word([mul(3,inv(3,direction),1<<j) for j in range(4)],range(4))
    negative=_native_x(2)+_native_x(3)
    inner=relative_pair(*pair)
    if negative_phase:inner += [('H',3),('CX',2,3),('H',3)]
    return coords+negative+inner+negative+_inverse(coords)

def common_a():
    # g=u2 OR u3. The inner word is i^g P^g, P=Z_V1.
    # (iP)^c(iP)^d(-iP)^(cd) = i^(c OR d) P^(c OR d).
    tau=paired_t(19,50)
    cliff=[('S',2),('S',3),('H',5),('CX',2,5),('CX',3,5),('H',5)]
    relative=[('H',5)]+_relative_toffoli(2,3,5)+[('H',5)]
    return _inverse(tau)+cliff+relative+tau

def coupled_linear_word():
    return (_baseline_word()
            +selected_pair(5,(156,228,64))
            +selected_pair(6,(89,247,44))
            +common_a()
            +selected_pair(4,(156,228,112))
            +selected_pair(7,(89,247,87),negative_phase=True))

def exterior_after():
    wall=[g for j in range(4) for g in [('CX',j,4+j),('CX',4+j,j)]]
    return wall+fourier_word(2,2)+fourier_word(2,6)

def recoupler_word():
    batch,meta=displacement_batch_word()
    return fourier_word(3,4)+coupled_linear_word()+batch+exterior_after()

def decoder_word():
    bottom=[g for gate in full_intertwiner()
            for g in [gate,(gate[0],*(j+4 for j in gate[1:]))]]
    return recoupler_word()+bottom

def logical_ring():
    # The separately verified six-helper batch is exactly Shift tensor I_helpers.
    # Compose the effective eight-wire permutation with exact native principal walls.
    before=ring_matrix(fourier_word(3,4)+coupled_linear_word())
    perm=[]
    for basis in range(256):
        u=basis&15;v=basis>>4
        shift=_raw_linear_cell(u)[1] if u else 0
        perm.append(u|((v^shift)<<4))
    moved=before[0][:,perm,:]
    # Clifford suffix is evaluated as its exact row permutation/linear combination
    # by applying the existing evaluator to a composed native substitute. A simple
    # matrix product in Z[omega] retains the complete effective logical operator.
    return multiply_ring(ring_matrix(exterior_after()),(moved,before[1]))

def multiply_ring(left,right):
    a,da=left;b,db=right
    out=np.zeros_like(a)
    bound=max(1,int(np.max(abs(a))))*max(1,int(np.max(abs(b))))*a.shape[1]*4
    assert bound<1<<61
    for i in range(4):
        for j in range(4):out[(i+j)%4]+=(-1 if i+j>=4 else 1)*(a[i]@b[j])
    den=da+db
    while den and not np.any(out&1):out//=2;den-=1
    return out,den

def verify():
    # All nonzero-U linear fibers, with literal phases, all input/output columns.
    actual=ring_matrix(coupled_linear_word())
    expected=np.zeros_like(actual[0]);expected_den=1
    for u in range(1,16):
        local=[(g[0],*(j-4 for j in g[1:])) for g in _raw_linear_cell(u)[0]]
        block,den=ring_matrix(local,4)
        assert den<=1
        inds=np.arange(16)*16+u
        for power in range(4):expected[power][np.ix_(inds,inds)]=block[power]*(1<<(1-den))
    nonzero=[j for j in range(256) if j&15]
    assert equal_ring((actual[0][:,:,nonzero],actual[1]),(expected[:,:,nonzero],expected_den))
    # Native actual inverse on all eight principal inputs, before adding verified I helpers.
    assert equal_ring(ring_matrix(_inverse(coupled_linear_word())),adjoint_ring(actual))
    r=logical_ring()
    bottom=[g for gate in full_intertwiner() for g in [gate,(gate[0],*(j+4 for j in gate[1:]))]]
    e=multiply_ring(ring_matrix(bottom),r)
    identity=np.zeros_like(e[0]);identity[0]=np.eye(256,dtype=np.int64)
    assert equal_ring(multiply_ring(adjoint_ring(e),e),(identity,0))
    for addr in range(8):
        source=[_tower_shear(3,addr,j) for j in range(256)]
        target=[_translation(3,addr,j) for j in range(256)]
        assert np.array_equal(e[0][:,:,source],e[0][:,target,:])
    report={'linear':ledger(coupled_linear_word()),'recoupler':ledger(recoupler_word()),
            'decoder':ledger(decoder_word()),'inverse_recoupler':ledger(_inverse(recoupler_word())),
            'inverse_decoder':ledger(_inverse(decoder_word())),
            'verified':'all240nonzeroUfibralcolumns;all256linearactualinversecolumns;effective8wireunitarityandall8queryintertwiners;privatehelpersverifiedinquarticbatchreceipt',
            'scope':'scratchonly;uniformhigharitynonzeroPredicateandcoupledlevelsremainunpaid'}
    with open(OUT/'coupled_receipt.json','w') as f:json.dump(report,f,indent=2)
    print(json.dumps(report,indent=2))
if __name__=='__main__':verify()
