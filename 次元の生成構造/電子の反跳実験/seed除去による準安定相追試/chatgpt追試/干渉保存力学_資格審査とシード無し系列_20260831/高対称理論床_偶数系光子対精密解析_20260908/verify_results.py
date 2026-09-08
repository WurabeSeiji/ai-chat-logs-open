#!/usr/bin/env python3
import csv, math
from pathlib import Path
p=Path(__file__).resolve().parent
rows=list(csv.DictReader((p/'even_N_summary.csv').open(encoding='utf-8')))
exact=[int(r['N']) for r in rows if int(r['candidate_pairs'])>0]
assert exact==[8,16,24,32,40], exact
for r in rows:
    N=int(r['N']); M=int(r['M']); c=int(r['candidate_pairs']); cov=float(r['matching_coverage'])
    if N%8==0:
        assert c==N*(4*N-11)//2
        assert abs(cov-1)<1e-12
        assert int(r['component_size4'])==N//8
        assert int(r['component_size8'])==N//8
        assert int(r['component_size16'])==N*(N-4)//32
        assert int(r['degree2_vertices'])==N//2
        assert int(r['degree4_vertices'])==N
        assert int(r['degree8_vertices'])==M-3*N//2
    else:
        assert c==0 and abs(cov)<1e-15
print('VERIFY PASS: exact pair N=8,16,24,32,40; all closed-form census identities hold.')
