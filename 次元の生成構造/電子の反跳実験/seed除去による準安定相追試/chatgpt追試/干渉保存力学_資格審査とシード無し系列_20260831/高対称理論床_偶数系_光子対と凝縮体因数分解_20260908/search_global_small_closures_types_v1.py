#!/usr/bin/env python3
import csv
from pathlib import Path
from collections import defaultdict, Counter
import numpy as np
from audit_highsym_even_photon_pairs_v2 import build_theoretical_floor
HERE=Path(__file__).resolve().parent
TOL=1e-10
ROUND=11

def ckey(x): return (round(float(x.real),ROUND),round(float(x.imag),ROUND))
def valid_mult(types, inds):
    c=Counter(inds)
    return all(c[i] <= types[i]['count'] for i in c)
def residual(types, inds):
    vals=[types[i]['value'] for i in inds]
    return abs(sum(vals))/max(sum(abs(v) for v in vals),1e-300)

def find_small(types,k):
    lookup=defaultdict(list)
    for i,t in enumerate(types): lookup[ckey(t['value'])].append(i)
    n=len(types)
    if k==3:
        best=(1e9,None)
        for a in range(n):
            for b in range(a,n):
                if not valid_mult(types,[a,b]): continue
                target=-(types[a]['value']+types[b]['value'])
                for c in lookup.get(ckey(target),[]):
                    inds=sorted([a,b,c])
                    if not valid_mult(types,inds): continue
                    r=residual(types,inds)
                    if r<best[0]: best=(r,inds)
                    if r<TOL:return best
        return best
    if k==4:
        sums=defaultdict(list)
        pairs=[]
        for a in range(n):
            for b in range(a,n):
                if not valid_mult(types,[a,b]): continue
                s=types[a]['value']+types[b]['value']
                sums[ckey(s)].append((a,b)); pairs.append((a,b,s))
        best=(1e9,None)
        for a,b,s in pairs:
            for c,d in sums.get(ckey(-s),[]):
                inds=sorted([a,b,c,d])
                if not valid_mult(types,inds):continue
                r=residual(types,inds)
                if r<best[0]:best=(r,inds)
                if r<TOL:return best
        return best

def main():
    rows=[]; exrows=[]
    for N in range(4,41,2):
        edges,d,z,lam,H,nd=build_theoretical_floor(N)
        buckets=defaultdict(list)
        for i,w in enumerate(z*z): buckets[ckey(w)].append(i)
        types=[]
        for key,idxs in buckets.items():
            v=np.mean([z[i]*z[i] for i in idxs])
            types.append({'value':v,'count':len(idxs),'indices':idxs,'dset':sorted(set(int(d[i]) for i in idxs))})
        tr,ti=find_small(types,3); qr,qi=find_small(types,4)
        row={'N':N,'M':len(z),'unique_squared_wave_types':len(types),'triplet_exact':tr<TOL,'triplet_best_relative':tr,
             'quad_exact':qr<TOL,'quad_best_relative':qr}
        rows.append(row)
        for label,inds,r in [('triplet',ti,tr),('quad',qi,qr)]:
            if inds and r<TOL:
                actual=[]; need=Counter(inds)
                for tix,cnt in need.items(): actual.extend(types[tix]['indices'][:cnt])
                exrows.append({'N':N,'kind':label,'relative_closure':r,
                               'edge_indices':';'.join(map(str,actual)),
                               'edges':';'.join(f'{edges[i][0]}-{edges[i][1]}' for i in actual),
                               'distance_classes':';'.join(map(str,[int(d[i]) for i in actual]))})
        print(N,len(types),'tri',tr<TOL,tr,'quad',qr<TOL,qr)
    with (HERE/'global_small_closure_type_search.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    with (HERE/'global_small_closure_examples.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(exrows[0]));w.writeheader();w.writerows(exrows)
if __name__=='__main__':main()
