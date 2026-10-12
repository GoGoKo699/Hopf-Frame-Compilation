#!/usr/bin/env python3
"""Regenerate compact exact-native readout/prefix evidence, without a fine source.

Optional gate-list export is a local artifact, not an additional canonical
copy of the circuits. Simulation uses floating arithmetic on explicitly
stated finite vectors; exact native identities justify arbitrary inputs.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from compiler_robust_hopf.native_ordered_banks import (  # noqa: E402
    EchoProgram, emit_prefix_bank, emit_readout, simulate,
)
from compiler_robust_hopf.ordered_block_banks import history_geometry, prefix_matrix  # noqa: E402


def run(n: int, output: Path | None = None) -> dict:
    if n not in (2, 3):
        raise ValueError("This bounded diagnostic supports n=2 or n=3.")
    seeds=[('h','t'),('t','h','tdg'),('h','t','h'),('t','h','t','h')]
    programs=[[EchoProgram(seeds[(d+p)%len(seeds)]) for p in range(1<<d)] for d in range(n-1)]
    readout=emit_readout(n)
    prefix=emit_prefix_bank(n,programs)
    cases=[('right_readout',readout,history_geometry(n).unitary(),n+2),
           ('coarse_prefix',prefix,prefix_matrix(n,[[p.matrix() for p in row] for row in programs]),n+1)]
    rng=np.random.default_rng(120102026)
    records=[]
    for name,emitter,ideal,w in cases:
        D=1<<emitter.qubits;M=1<<w
        vectors=rng.normal(size=(D,2))+1j*rng.normal(size=(D,2))
        if name=='right_readout': vectors[(np.arange(D)%M)>=(1<<n),0]=0
        vectors/=np.linalg.norm(vectors,axis=0)
        expected=np.einsum('ij,bjk->bik',ideal,vectors.reshape(-1,M,2)).reshape(D,2)
        error=float(np.linalg.norm(simulate(emitter.gates,vectors)-expected,2))
        if error>2e-10: raise AssertionError(f"{name} finite native check failed: {error}")
        gates=json.dumps(emitter.gates,separators=(',',':')).encode()
        record={'kind':name,'n':n,'mode_qubits':w,'dirty_helpers':emitter.qubits-w,
                'elementary_gates':len(emitter.gates),
                'T_count':sum(g[0] in ('t','tdg') for g in emitter.gates),
                'gate_list_sha256':hashlib.sha256(gates).hexdigest(),
                'vectors_checked':2,'floating_vector_error':error}
        records.append(record)
        if output is not None:
            output.mkdir(parents=True,exist_ok=True)
            (output/f'{name}_n{n}.json').write_text(json.dumps({**record,'gates':emitter.gates},indent=2)+'\n')
    return {'scope':'Exact-native component words only; finite floating vector checks, not fine source or endpoint synthesis.',
            'cases':records}


def main() -> None:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--n',type=int,choices=(2,3),default=3)
    parser.add_argument('--write-gates',type=Path,help='Optional local directory for generated gate lists.')
    args=parser.parse_args()
    print(json.dumps(run(args.n,args.write_gates),indent=2))


if __name__=='__main__':
    main()
