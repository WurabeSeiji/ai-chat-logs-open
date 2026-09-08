#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Audit spontaneous two-wave zero-closure pairs in canonical make_parent and handcrafted symmetric parents.

Photon-like pair criterion:
  z_i^2 + z_j^2 = 0  <=>  z_j = +/- i z_i  (for nonzero waves)
Primary dimensionless residual:
  eps2 = |z_i^2 + z_j^2| / (|z_i|^2 + |z_j|^2)
Also records amplitude mismatch, +/-90 degree phase residual, and direct +/-i residual.
"""
from __future__ import annotations
import argparse, csv, hashlib, json, math, os, re
from pathlib import Path
import numpy as np
import networkx as nx

TOLS = (1e-12, 1e-9, 1e-6)
TOPK = 10

def wrap_pi(x):
    return (x + np.pi) % (2*np.pi) - np.pi

def edge_labels(N):
    a,b=np.triu_indices(N,k=1)
    return list(zip(a.tolist(), b.tolist()))

def analyze_state(series, N, state_name, z):
    z=np.asarray(z,dtype=np.complex128).ravel(); M=len(z); labels=edge_labels(N)
    ii,jj=np.triu_indices(M,k=1)
    zi=z[ii]; zj=z[jj]; ai=np.abs(zi); aj=np.abs(zj)
    den=ai*ai+aj*aj
    eps=np.abs(zi*zi+zj*zj)/den
    dphi=wrap_pi(np.angle(zj)-np.angle(zi))
    rp=np.abs(wrap_pi(dphi-np.pi/2)); rm=np.abs(wrap_pi(dphi+np.pi/2))
    plus=rp<=rm; pres=np.where(plus,rp,rm)
    meanamp=(ai+aj)/2
    amp_rel=np.abs(ai-aj)/meanamp
    direct=np.minimum(np.abs(zj-1j*zi),np.abs(zj+1j*zi))/np.sqrt(den)
    order=np.argsort(eps)
    strict_idx=np.flatnonzero(eps<=1e-12)
    selected=np.unique(np.concatenate([order[:TOPK],strict_idx]))
    rows=[]
    for k in selected:
        i=int(ii[k]); j=int(jj[k]); e1,e2=labels[i],labels[j]
        rows.append(dict(series=series,N=N,M=M,state=state_name,
            i=i,j=j,edge_i=f'{e1[0]}-{e1[1]}',edge_j=f'{e2[0]}-{e2[1]}',
            abs_zi=float(ai[k]),abs_zj=float(aj[k]),amp_ratio=float(aj[k]/ai[k]) if ai[k] else float('inf'),
            amp_rel_mismatch=float(amp_rel[k]),dphi_rad=float(dphi[k]),dphi_deg=float(np.degrees(dphi[k])),
            helicity='+i' if plus[k] else '-i',phase_resid_rad=float(pres[k]),phase_resid_deg=float(np.degrees(pres[k])),
            pair_closure_abs=float(np.abs(zi[k]*zi[k]+zj[k]*zj[k])),pair_closure_norm=float(eps[k]),
            plusminus_i_resid=float(direct[k])))
    strict=[r for r in rows if r['pair_closure_norm']<=1e-12]
    unique=set(); hp=hm=0
    for r in strict:
        unique|={r['i'],r['j']}; hp += r['helicity']=='+i'; hm += r['helicity']=='-i'
    denom=float(np.sum(np.abs(z)**2)); global_norm=float(abs(np.sum(z*z))/denom) if denom else float('nan')
    b=int(order[0]); bi=int(ii[b]); bj=int(jj[b]); be1,be2=labels[bi],labels[bj]
    summary=dict(series=series,N=N,M=M,state=state_name,
      global_closure_abs=float(abs(np.sum(z*z))),global_closure_norm=global_norm,pair_count=len(ii),
      best_pair_closure_norm=float(eps[b]),best_i=bi,best_j=bj,best_edge_i=f'{be1[0]}-{be1[1]}',best_edge_j=f'{be2[0]}-{be2[1]}',
      best_helicity='+i' if plus[b] else '-i',best_amp_rel_mismatch=float(amp_rel[b]),best_phase_resid_deg=float(np.degrees(pres[b])),
      count_le_1e12=int(np.sum(eps<=1e-12)),count_le_1e9=int(np.sum(eps<=1e-9)),count_le_1e6=int(np.sum(eps<=1e-6)),
      strict_unique_waves=len(unique),strict_wave_fraction=len(unique)/M,strict_plus_i=hp,strict_minus_i=hm)
    return rows, summary

def analyze_file(series, p):
    m=re.search(r'N(\d{5})',p.name); N=int(m.group(1))
    d=np.load(p,allow_pickle=False)
    out=[]; summaries=[]
    for state in ('v','Z0'):
        z=np.asarray(d[state],complex).ravel()
        rows, summary = analyze_state(series,N,state,z)
        out.extend(rows); summaries.append(summary)
    return out,summaries

def write_csv(path,rows):
    if not rows: return
    with open(path,'w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)

def sha256(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''): h.update(b)
    return h.hexdigest()

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--makeparent',required=True)
    ap.add_argument('--symmetric',required=True)
    ap.add_argument('--out',required=True)
    args=ap.parse_args()
    out=Path(args.out); out.mkdir(parents=True,exist_ok=True)
    allpairs=[]; summaries=[]
    inputs=[]
    for series,folder in [('make_parent',Path(args.makeparent)),('high_symmetry',Path(args.symmetric))]:
        files=sorted(folder.glob('parent_static_N*_makeparent_20260905.npz'))
        if len(files)!=38:
            raise RuntimeError(f'{series}: expected 38 files N=3..40, got {len(files)}')
        ns=[]
        for p in files:
            m=re.search(r'N(\d{5})',p.name); ns.append(int(m.group(1)))
            pr,sr=analyze_file(series,p); allpairs.extend(pr); summaries.extend(sr)
            inputs.append({'series':series,'file':p.name,'sha256':sha256(p)})
        if ns != list(range(3,41)):
            raise RuntimeError(f'{series}: N sequence mismatch {ns[:3]}..{ns[-3:]}')
    allpairs.sort(key=lambda r:(r['series'],r['state'],r['N'],r['pair_closure_norm']))
    summaries.sort(key=lambda r:(r['state'],r['N'],r['series']))
    top=[]; strict=[]
    groups={}
    for r in allpairs: groups.setdefault((r['series'],r['N'],r['state']),[]).append(r)
    for key,rows in groups.items():
        rows=sorted(rows,key=lambda r:r['pair_closure_norm'])
        for rank,r in enumerate(rows[:TOPK],1): top.append(dict(rank=rank,**r))
        strict.extend(r for r in rows if r['pair_closure_norm']<=1e-12)
    write_csv(out/'pair_summary_N3_N40.csv',summaries)
    write_csv(out/'pair_top10_N3_N40.csv',top)
    write_csv(out/'strict_photon_pairs_1e-12.csv',strict)
    matching_rows=[]
    for series in ('make_parent','high_symmetry'):
        for state in ('v','Z0'):
            for N in range(3,41):
                M=N*(N-1)//2
                sr=[r for r in strict if r['series']==series and r['state']==state and r['N']==N]
                G=nx.Graph(); G.add_nodes_from(range(M)); G.add_edges_from((r['i'],r['j']) for r in sr)
                matching=nx.algorithms.matching.max_weight_matching(G,maxcardinality=True)
                matched={x for e in matching for x in e}
                matching_rows.append(dict(series=series,state=state,N=N,M=M,strict_edges=len(sr),
                    max_disjoint_pairs=len(matching),matched_waves=len(matched),
                    perfect_pair_partition=(len(matched)==M),unmatched_waves=M-len(matched)))
    write_csv(out/'strict_pair_matching_summary.csv',matching_rows)
    meta={'created':'2026-09-08','criterion':'abs(zi^2+zj^2)/(abs(zi)^2+abs(zj)^2)',
          'strict_tolerance':1e-12,'other_tolerances':[1e-9,1e-6],
          'series':['make_parent','high_symmetry'],'states':['v','Z0'],'inputs':inputs}
    (out/'RUN_METADATA.json').write_text(json.dumps(meta,ensure_ascii=False,indent=2),encoding='utf-8')
    print(f'summaries={len(summaries)} pairs={len(allpairs)} strict={len(strict)}')
    # compact findings for downstream markdown generation
    for state in ('v','Z0'):
        for series in ('make_parent','high_symmetry'):
            ss=[r for r in summaries if r['state']==state and r['series']==series]
            exactNs=[r['N'] for r in ss if r['count_le_1e12']>0]
            best=min(ss,key=lambda r:r['best_pair_closure_norm'])
            print(state,series,'Ns_with_strict',exactNs,'best',best['N'],best['best_pair_closure_norm'],'count',best['count_le_1e12'])

if __name__=='__main__': main()
