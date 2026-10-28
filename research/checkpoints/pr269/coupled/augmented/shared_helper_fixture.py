"""Exact helper-X-sector fixture for the augmented uniform recoupler identity.

This emits Clifford sector words to verify the analytic augmented construction;
its physical multiplication-XOR emitter is the retained all-size primitive.
It is not a constant-depth all-size decoder.
"""
import sys,json,itertools
from pathlib import Path
import numpy as np
REPO = next(p for p in Path(__file__).resolve().parents if (p / 'WORKSPACE.md').is_file() and (p / 'tests').is_dir())
sys.path.insert(0,str(REPO))
from tests.test_recursive_field_shear_decoder import (_tower_product as mul,
    _tower_trace as trace,_address_field,_tower_shear)
from tests.test_projective_field_recoupler import (fourier_word,linear_word,
    ring_matrix,equal_ring,adjoint_ring,_inverse)
from tests.test_field_shear_intertwiner import _translation

def phase_value(level,helper_character,basis,kind):
    ell=1<<(level-2);k=2*ell;mask=(1<<ell)-1
    u=basis&((1<<k)-1);v=basis>>k
    d0,d1=u&mask,u>>ell;a0,a1=v&mask,v>>ell
    if kind=='before':
        f=(mul(level,u,v)>>ell)^mul(level-1,d0,d1)
    else:
        assert kind=='after'
        f=mul(level-1,d0,a1)^mul(level-1,d1,a0)
    return trace(level-1,mul(level-1,helper_character,f))

def phase_word(level,helper_character,kind,offset):
    width=1<<level
    linear=[phase_value(level,helper_character,1<<j,kind) for j in range(width)]
    word=[]
    for j,c in enumerate(linear):
        if c:word += [('S',offset+j)]*2
    quadratic=[]
    for i in range(width):
        for j in range(i+1,width):
            c=phase_value(level,helper_character,(1<<i)|(1<<j),kind)^linear[i]^linear[j]
            if c:
                quadratic.append((i,j));word += [('H',offset+j),('CX',offset+i,offset+j),('H',offset+j)]
    # Complete bounded truth table verifies this exact quadratic expansion.
    for basis in range(1<<width):
        val=sum(c*((basis>>j)&1) for j,c in enumerate(linear))
        val+=sum(((basis>>i)&1)*((basis>>j)&1) for i,j in quadratic)
        assert val%2==phase_value(level,helper_character,basis,kind)
    return word

def recoupler_sector_word(level,r,offset=0):
    ell=1<<(level-2);k=2*ell
    theta=1<<ell
    htheta=(linear_word([mul(level,theta,1<<j) for j in range(k)],range(offset+k,offset+2*k))
            +fourier_word(level,offset+k))
    sumwall=[g for j in range(k) for g in [('CX',offset+j,offset+k+j),('CX',offset+k+j,offset+j)]]
    return (fourier_word(level,offset+k)
            +phase_word(level,r,'before',offset)
            +htheta+phase_word(level,r,'after',offset)
            +sumwall+fourier_word(level-1,offset+ell)
            +fourier_word(level-1,offset+k+ell))

def decoder_sector_word(level,characters,offset=0):
    if level==1:return [('CX',offset+1,offset)]
    half=1<<(level-1)
    return (recoupler_sector_word(level,characters[level],offset)
            +decoder_sector_word(level-1,characters,offset)
            +decoder_sector_word(level-1,characters,offset+half))

def run():
    reports=[]
    for level in(2,3):
        sizes=[1<<(1<<(j-2)) for j in range(2,level+1)]
        sectors=list(itertools.product(*(range(size) for size in sizes)))
        matrices={}
        for sector in sectors:
            chars={j:sector[j-2] for j in range(2,level+1)}
            word=decoder_sector_word(level,chars)
            matrix=ring_matrix(word,1<<level)
            assert equal_ring(ring_matrix(_inverse(word),1<<level),adjoint_ring(matrix))
            matrices[sector]=matrix
        dimension=1<<(1<<level);comparisons=0
        for sector in sectors:
            for addr in range(1<<level):
                shifted=tuple(sector[j-2]^_address_field(j-1,addr&((1<<(j-1))-1)) for j in range(2,level+1))
                source=[_tower_shear(level,addr,col) for col in range(dimension)]
                target=[_translation(level,addr,row) for row in range(dimension)]
                left,ld=matrices[shifted];right,rd=matrices[sector]
                assert equal_ring((left[:,:,source],ld),(right[:,target,:],rd)),(level,sector,addr)
                comparisons+=dimension
        reports.append({'p':level,'core_qubits':1<<level,'shared_moving_helper_qubits':sum(1<<(j-2) for j in range(2,level+1)),
                        'helper_X_sectors':len(sectors),'address_column_comparisons':comparisons,
                        'actual_inverse_sector_columns':len(sectors)*dimension})
    report={'checks':reports,'verdict':'PASS','scope':'exactCliffordsectorfixtureforaugmentedalgebraicidentity;retainedmultiplicationXORproofownsnativeledger;overall16(p−1)bandsnotconstant'}
    with open(Path(__file__).with_name('shared_helper_receipt.json'),'w') as f:json.dump(report,f,indent=2)
    print(json.dumps(report,indent=2))
if __name__=='__main__':run()
