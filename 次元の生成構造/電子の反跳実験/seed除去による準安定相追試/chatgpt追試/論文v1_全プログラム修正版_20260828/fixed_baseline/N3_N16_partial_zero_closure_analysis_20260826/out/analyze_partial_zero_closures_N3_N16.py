from pathlib import Path
import csv, math, itertools, json, hashlib, zipfile, shutil
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

BASE=Path(__file__).resolve().parent.parent  # PATH PATCH: was Path('/mnt/data') (sibling packages)
OUT=Path(__file__).resolve().parent/'out'  # PATH PATCH: was BASE/'N3_N16_partial_zero_closure_analysis_20260826' (would rmtree this package)
if OUT.exists(): shutil.rmtree(OUT)
OUT.mkdir()

PATHS={
3:BASE/'N3_N4_complex_simplex_complete_analysis_20260826/N3_step5000_final_edges.csv',
4:BASE/'N3_N4_complex_simplex_complete_analysis_20260826/N4_step5000_final_edges.csv',
6:BASE/'N6_N7_complex_simplex_complete_analysis_20260826/N6_step5000_final_edges.csv',
7:BASE/'N6_N7_complex_simplex_complete_analysis_20260826/N7_step5000_final_edges.csv',
8:BASE/'N8_N9_complex_simplex_complete_analysis_20260826/N8_step5000_final_edges.csv',
9:BASE/'N8_N9_complex_simplex_complete_analysis_20260826/N9_step5000_final_edges.csv',
10:BASE/'N10_N11_complex_simplex_complete_analysis_20260826/N10_step5000_final_edges.csv',
11:BASE/'N10_N11_complex_simplex_complete_analysis_20260826/N11_step5000_final_edges.csv',
12:BASE/'N12_N13_complex_simplex_complete_analysis_20260826/N12_step5000_final_edges.csv',
13:BASE/'N12_N13_complex_simplex_complete_analysis_20260826/N13_step5000_final_edges.csv',
14:BASE/'N14_N15_complex_simplex_complete_analysis_20260826/N14_step5000_final_edges.csv',
15:BASE/'N14_N15_complex_simplex_complete_analysis_20260826/N15_step5000_final_edges.csv',
16:BASE/'N16_complex_simplex_complete_analysis_20260826/N16_step5000_final_edges.csv',
}

def edges_for_n(n): return [(i,j) for i in range(1,n+1) for j in range(i+1,n+1)]

def load_n(n):
    if n==5:
        p=BASE/'N5_complex_simplex_complete_analysis_20260826/N5_all_steps_a_b_a2_b2_ab.csv'
        d=pd.read_csv(p)
        d=d[d.step==5000].copy().sort_values('edge_index')
        es=edges_for_n(n)
        d['i']=[e[0] for e in es]; d['j']=[e[1] for e in es]
        if 'z2_re' not in d: d['z2_re']=d.a2-d.b2
        if 'z2_im' not in d: d['z2_im']=2*d.ab
        return d
    return pd.read_csv(PATHS[n]).sort_values('edge_index')

def normres(vals):
    den=np.sum(np.abs(vals))
    return float(abs(np.sum(vals))/den) if den else 0.0

def round_robin_matchings(n):
    # edge-coloring: even n -> n-1 perfect matchings; odd n -> n near-perfect matchings via dummy
    dummy=None
    nn=n
    verts=list(range(1,n+1))
    if n%2:
        dummy=0; verts.append(dummy); nn=n+1
    arr=verts[:]
    rounds=[]
    for _ in range(nn-1):
        pairs=[]
        for k in range(nn//2):
            a,b=arr[k],arr[nn-1-k]
            if a!=dummy and b!=dummy: pairs.append(tuple(sorted((a,b))))
        rounds.append(pairs)
        arr=[arr[0]]+[arr[-1]]+arr[1:-1]
    return rounds

def walecki_odd_cycles(n):
    # K_(2m+1) decomposes into m Hamilton cycles using infinity + Z_(2m)
    assert n%2==1
    m=(n-1)//2; mod=2*m; inf=n
    cycles=[]
    for r in range(m):
        seq=[inf]
        # r, r-1, r+1, r-2, r+2, ... , r-m
        seq.append((r%mod)+1)
        for k in range(1,m):
            seq.append(((r-k)%mod)+1)
            seq.append(((r+k)%mod)+1)
        seq.append(((r-m)%mod)+1)
        ed=[]
        for a,b in zip(seq,seq[1:]+seq[:1]): ed.append(tuple(sorted((a,b))))
        cycles.append(ed)
    return cycles

def even_hamilton_plus_matching(n):
    # obtain Hamilton decomposition of K_(n+1), remove extra vertex => (n/2) Hamilton paths; not an edge partition into cycles.
    # For even n use standard round-robin perfect matching partition as the canonical exact edge decomposition.
    return None

def subset_search(z2, edges, size):
    best=(1e9,None,None)
    hits=[]
    for comb in itertools.combinations(range(len(z2)),size):
        vals=z2[list(comb)]
        r=normres(vals)
        if r<best[0]: best=(r,comb,np.sum(vals))
        if r<1e-6: hits.append((r,comb,np.sum(vals)))
    hits.sort(key=lambda x:x[0])
    return best,hits

def exact_cover_triples(M, triples, limit=1000):
    # triples list (res, comb, sum); find disjoint exact covers, best first. Only feasible when M divisible by 3.
    if M%3: return []
    by_edge={i:[] for i in range(M)}
    for idx,t in enumerate(triples):
        for e in t[1]: by_edge[e].append(idx)
    sols=[]
    def rec(used, chosen):
        if len(sols)>=limit: return
        if len(used)==M:
            sols.append(chosen.copy()); return
        # uncovered edge with fewest candidates compatible
        unc=[e for e in range(M) if e not in used]
        e=min(unc,key=lambda x:sum(1 for ti in by_edge[x] if not (set(triples[ti][1])&used)))
        cand=[]
        for ti in by_edge[e]:
            comb=set(triples[ti][1])
            if not comb&used: cand.append(ti)
        for ti in cand:
            comb=set(triples[ti][1]); rec(used|comb,chosen+[ti])
    rec(set(),[])
    return sols

def exact_cover_subsets(M, subsets, limit=100):
    if not subsets: return []
    size=len(subsets[0][1])
    if M%size: return []
    by_edge={i:[] for i in range(M)}
    for idx,t in enumerate(subsets):
        for e in t[1]: by_edge[e].append(idx)
    sols=[]
    def rec(used,chosen):
        if len(sols)>=limit: return
        if len(used)==M:
            sols.append(chosen.copy()); return
        unc=[e for e in range(M) if e not in used]
        e=min(unc,key=lambda x:sum(1 for ti in by_edge[x] if not (set(subsets[ti][1])&used)))
        for ti in by_edge[e]:
            comb=set(subsets[ti][1])
            if not comb&used: rec(used|comb,chosen+[ti])
    rec(set(),[])
    return sols

summary=[]; struct_rows=[]; subset_rows=[]; cover_rows=[]
all_best={}
for n in range(3,17):
    d=load_n(n); edges=[(int(i),int(j)) for i,j in zip(d.i,d.j)]
    z2=d.z2_re.to_numpy()+1j*d.z2_im.to_numpy(); M=len(z2)
    total=normres(z2)
    r2=np.abs(z2)
    # canonical round-robin edge decomposition
    rr=round_robin_matchings(n)
    emap={e:k for k,e in enumerate(edges)}
    rrres=[]
    for bi,block in enumerate(rr):
        inds=[emap[e] for e in block]; r=normres(z2[inds]); rrres.append(r)
        struct_rows.append([n,M,'round_robin_matching',bi,len(inds),r,','.join(f'{a}-{b}' for a,b in block)])
    # odd Hamilton cycle decomposition
    hcres=[]
    if n%2==1:
        cycles=walecki_odd_cycles(n)
        # validate unique cover
        flat=[e for c in cycles for e in c]
        assert len(flat)==M and len(set(flat))==M
        for bi,block in enumerate(cycles):
            inds=[emap[e] for e in block]; r=normres(z2[inds]); hcres.append(r)
            struct_rows.append([n,M,'walecki_hamilton_cycle',bi,len(inds),r,','.join(f'{a}-{b}' for a,b in block)])
    # exact small subset searches
    bests={}; hitcounts={}
    hits2=[]; hits3=[]
    for s in (2,3):
        best,hits=subset_search(z2,edges,s); bests[s]=best[0]; hitcounts[s]=len(hits)
        subset_rows.append([n,M,s,best[0],len(hits),','.join(f'{edges[k][0]}-{edges[k][1]}' for k in best[1])])
        if s==2: hits2=hits
        if s==3: hits3=hits
    # exact-cover by zero-closure pairs/triples at tol1e-6
    paircovers=exact_cover_subsets(M,hits2,limit=100) if M%2==0 else []
    for ci,sol in enumerate(paircovers[:20]):
        maxr=max(hits2[ti][0] for ti in sol); meanr=np.mean([hits2[ti][0] for ti in sol])
        blocks=['+'.join(f'{edges[k][0]}-{edges[k][1]}' for k in hits2[ti][1]) for ti in sol]
        cover_rows.append([n,M,2,ci,len(sol),maxr,meanr,' | '.join(blocks)])
    covers=[]
    if M%3==0 and hits3:
        covers=exact_cover_triples(M,hits3,limit=100)
        for ci,sol in enumerate(covers[:20]):
            maxr=max(hits3[ti][0] for ti in sol); meanr=np.mean([hits3[ti][0] for ti in sol])
            blocks=['+'.join(f'{edges[k][0]}-{edges[k][1]}' for k in hits3[ti][1]) for ti in sol]
            cover_rows.append([n,M,3,ci,len(sol),maxr,meanr,' | '.join(blocks)])
    summary.append([n,M,total,bests[2],hitcounts[2],len(paircovers),bests[3],hitcounts[3],len(covers),
                    min(rrres),np.mean(rrres),max(rrres),
                    min(hcres) if hcres else np.nan,np.mean(hcres) if hcres else np.nan,max(hcres) if hcres else np.nan])
    all_best[n]=(z2,edges)

cols=['N','M','total_closure_residual','best_pair_residual','pair_hits_lt1e-6','pair_exact_cover_count_capped100','best_triple_residual','triple_hits_lt1e-6','triple_exact_cover_count_capped100',
      'roundrobin_min','roundrobin_mean','roundrobin_max','hamilton_min','hamilton_mean','hamilton_max']
sdf=pd.DataFrame(summary,columns=cols); sdf.to_csv(OUT/'N3_N16_partial_closure_summary.csv',index=False)
pd.DataFrame(struct_rows,columns=['N','M','decomposition','block_id','block_size','closure_residual','edges']).to_csv(OUT/'N3_N16_structural_decomposition_blocks.csv',index=False)
pd.DataFrame(subset_rows,columns=['N','M','subset_size','best_closure_residual','hit_count_lt1e-6','best_edges']).to_csv(OUT/'N3_N16_best_small_subsets.csv',index=False)
pd.DataFrame(cover_rows,columns=['N','M','subset_size','cover_id','block_count','max_block_residual','mean_block_residual','blocks']).to_csv(OUT/'N3_N16_zero_triple_exact_covers.csv',index=False)

# plots
plt.figure(figsize=(9,5.5))
plt.semilogy(sdf.N,np.maximum(sdf.total_closure_residual,1e-18),marker='o',label='whole system')
plt.semilogy(sdf.N,np.maximum(sdf.best_pair_residual,1e-18),marker='o',label='best 2-edge subset')
plt.semilogy(sdf.N,np.maximum(sdf.best_triple_residual,1e-18),marker='o',label='best 3-edge subset')
plt.xlabel('N'); plt.ylabel('normalized closure residual'); plt.title('N=3..16: whole and best partial zero-closure residuals'); plt.legend(); plt.tight_layout(); plt.savefig(OUT/'N3_N16_best_partial_closure_residuals.png',dpi=180); plt.close()

plt.figure(figsize=(9,5.5))
plt.plot(sdf.N,sdf['triple_hits_lt1e-6'],marker='o',label='zero-close triples (<1e-6)')
plt.plot(sdf.N,sdf['pair_hits_lt1e-6'],marker='o',label='zero-close pairs (<1e-6)')
plt.xlabel('N'); plt.ylabel('count'); plt.title('N=3..16: number of near-exact small zero-closure subsets'); plt.legend(); plt.tight_layout(); plt.savefig(OUT/'N3_N16_small_zero_closure_counts.png',dpi=180); plt.close()

plt.figure(figsize=(9,5.5))
plt.semilogy(sdf.N,np.maximum(sdf.roundrobin_mean,1e-18),marker='o',label='round-robin matching mean')
odd=sdf[sdf.N%2==1]
plt.semilogy(odd.N,np.maximum(odd.hamilton_mean,1e-18),marker='o',label='odd-N Hamilton-cycle mean')
plt.xlabel('N'); plt.ylabel('mean normalized block residual'); plt.title('Canonical graph-factor decompositions: block closure'); plt.legend(); plt.tight_layout(); plt.savefig(OUT/'N3_N16_structural_decomposition_residuals.png',dpi=180); plt.close()

plt.figure(figsize=(9,5.5))
plt.plot(sdf.N,sdf.pair_exact_cover_count_capped100,marker='o',label='2-edge covers')
plt.plot(sdf.N,sdf.triple_exact_cover_count_capped100,marker='o',label='3-edge covers')
plt.xlabel('N'); plt.ylabel('exact covers found (cap 100)'); plt.title('Partition of all edges into near-zero subsets (<1e-6)'); plt.legend(); plt.tight_layout(); plt.savefig(OUT/'N3_N16_zero_triple_exact_cover_counts.png',dpi=180); plt.close()

# report derived strictly from calculations
lines=['# N=3〜16 部分ゼロ閉包・複数閉包分解解析','',
'## 定義','',
'全体系の関係波を `z_e` とし、複素二乗距離 `w_e=z_e^2` を使う。部分集合 B の閉包残差を', '',
'`C(B)=|sum_{e in B} w_e| / sum_{e in B}|w_e|`', '',
'とした。C=0 がその部分集合単独の二乗ゼロ閉包である。物理更新は一切変更せず、各Nの既存step=5000生データだけを再解析した。','',
'## 実施した探索','',
'- 全体系のゼロ閉包残差','- 全2辺部分集合の総当たり','- 全3辺部分集合の総当たり','- C<1e-6 の3辺閉包だけを用いた全辺exact-cover探索（最大100解）','- K_N の標準round-robin matching分解','- 奇数NについてWalecki Hamilton-cycle分解','',
'## 主要結果','']
for _,r in sdf.iterrows():
    lines.append(f"- N={int(r.N)}, M={int(r.M)}: total={r.total_closure_residual:.3e}, best2={r.best_pair_residual:.3e} (hits={int(r['pair_hits_lt1e-6'])}, covers={int(r.pair_exact_cover_count_capped100)}), best3={r.best_triple_residual:.3e} (hits={int(r['triple_hits_lt1e-6'])}), triple exact covers={int(r.triple_exact_cover_count_capped100)}")
lines += ['','## 解釈上の注意','',
'1. Mの整数因数分解だけでは部分ゼロ閉包は保証されない。必要なのは、実際の複素ベクトル `z_e^2` の部分和がゼロになること。',
'2. canonical matching / Hamilton分解で残差がゼロでないことは「別の頂点置換・別の組合せ分解にも解がない」ことの証明ではない。',
'3. 一方、2辺・3辺探索はそのサイズについて全組合せを総当たりしているので、step=5000データに対する存在判定として直接的である。',
'4. triple exact cover が見つかったNでは、全体系を複数の3辺ゼロ閉包へ実際に分割できる（tol=1e-6の数値判定）。',
'','## 次の判定軸','',
'3辺で閉じないNについては4辺以上の部分閉包探索が次段階になる。Mが大きいNでは全組合せ総当たりが急増するため、meet-in-the-middleまたはexact-cover最適化を用いるのが適切。']
(OUT/'ANALYSIS.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')

# copy exact analysis program into output
shutil.copy2(Path(__file__),OUT/'analyze_partial_zero_closures_N3_N16.py')
(OUT/'README.md').write_text('N=3..16 partial zero-closure analysis using existing step=5000 raw edge data. No physics rerun or modification. See ANALYSIS.md and CSV/PNG outputs.\n',encoding='utf-8')

# checksums
checks=[]
for p in sorted(OUT.iterdir()):
    if p.is_file() and p.name!='SHA256SUMS.txt': checks.append(f"{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.name}")
(OUT/'SHA256SUMS.txt').write_text('\n'.join(checks)+'\n',encoding='utf-8')

zip_path=BASE/'N3_N16_partial_zero_closure_analysis_20260826.zip'
with zipfile.ZipFile(zip_path,'w',zipfile.ZIP_DEFLATED) as z:
    for p in sorted(OUT.iterdir()):
        if p.is_file(): z.write(p,arcname=p.name)
print(sdf.to_string(index=False))
print('ZIP',zip_path,zip_path.stat().st_size)
