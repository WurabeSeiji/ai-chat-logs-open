#!/usr/bin/env python3
import csv, math, sys
from pathlib import Path
HERE=Path(__file__).resolve().parent

def spf(n):
    p=2
    while p*p<=n:
        if n%p==0:return p
        p+=1
    return n
rows=list(csv.DictReader((HERE/'odd_N_summary.csv').open(encoding='utf-8')))
assert len(rows)==19
for r in rows:
    N=int(r['N']); M=N*(N-1)//2
    assert int(r['M'])==M
    assert int(r['exact_2wave_pair_candidates'])==0
    assert int(r['max_matching_2wave_pairs'])==0
    assert float(r['2wave_pair_coverage'])==0.0
    assert int(r['unpaired_after_max_2wave'])==M
    assert r['global_exact_triplet_exists']=='True' if N%3==0 else r['global_exact_triplet_exists']=='False'
    assert r['global_exact_quad_exists']=='False'
    assert int(r['smallest_prime_factor_N'])==spf(N)
    assert abs(float(r['constructive_coverage'])-1.0)<1e-15
    assert int(r['leftover_after_constructive_factorization'])==0
    assert float(r['max_constructive_group_relative_closure'])<2e-14
    assert (r['guaranteed_3wave_regular_closure']=='True')==(N%3==0)
    assert (r['guaranteed_5wave_regular_closure']=='True')==(N%5==0)
# pair candidate CSV must have header only
with (HERE/'exact_2wave_pair_candidates.csv').open(encoding='utf-8') as f:
    lines=f.read().splitlines()
assert len(lines)==1, len(lines)
groups=list(csv.DictReader((HERE/'constructive_primitive_zero_closure_groups.csv').open(encoding='utf-8')))
byN={}
for g in groups: byN.setdefault(int(g['N']),[]).append(g)
for r in rows:
    N=int(r['N']); M=int(r['M']); p=int(r['constructive_minimal_block_size'])
    gs=byN[N]
    assert len(gs)==M//p
    assert all(int(g['group_size'])==p for g in gs)
    assert max(float(g['relative_closure']) for g in gs)<2e-14
print('PASS: all odd-N closure-factorization checks')
