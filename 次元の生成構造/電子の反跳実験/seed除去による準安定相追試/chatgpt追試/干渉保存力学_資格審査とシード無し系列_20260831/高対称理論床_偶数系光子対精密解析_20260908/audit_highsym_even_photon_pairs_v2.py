#!/usr/bin/env python3
"""High-symmetry theoretical-floor even-N audit for exact two-wave ±i closure pairs.

Reconstructs the cyclic Fourier theoretical-floor state directly from its defining
ansatz and reduced Hermitian eigenproblem, then audits all unordered wave pairs.
No make_parent data and no saved parent vectors are required.
"""
import csv, json, math, hashlib
from pathlib import Path
from collections import Counter
import numpy as np
import networkx as nx

HERE=Path(__file__).resolve().parent
TOL=1e-10
EVEN_NS=list(range(4,41,2))

def edges_for_n(N): return [(i,j) for i in range(N) for j in range(i+1,N)]
def cyc_dist(i,j,N): return min((j-i)%N,(i-j)%N)
def build_A(edges):
    m=len(edges); A=np.zeros((m,m),float); S=[set(e) for e in edges]
    for a in range(m):
        for b in range(a+1,m):
            if S[a]&S[b]: A[a,b]=A[b,a]=1.0
    return A

def build_theoretical_floor(N):
    edges=edges_for_n(N); M=len(edges); D=N//2
    d=np.array([cyc_dist(i,j,N) for i,j in edges],int)
    theta=np.array([2*np.pi*(i+j)/N for i,j in edges])
    A=build_A(edges)
    K=A*np.sin(theta[None,:]-theta[:,None]); H=1j*K
    B=np.zeros((M,D),complex); nd=[]
    for dd in range(1,D+1):
        idx=np.where(d==dd)[0]; nd.append(len(idx))
        B[idx,dd-1]=np.exp(1j*theta[idx])/np.sqrt(len(idx))
    Q=B.conj().T@H@B
    w,V=np.linalg.eigh(Q)
    chosen=None
    for kk in np.argsort(w):
        c=V[:,kk]; p=int(np.argmax(np.abs(c))); c*=np.exp(-1j*np.angle(c[p]))
        if np.max(np.abs(c.imag))<1e-9:
            cr=c.real
            if np.all(cr>=-1e-10) or np.all(cr<=1e-10):
                if np.sum(cr)<0: cr=-cr
                chosen=(float(w[kk]),cr); break
    if chosen is None: raise RuntimeError(f'No same-sign reduced eigenvector N={N}')
    _,c=chosen
    r=np.array([c[dd-1]/np.sqrt(nd[dd-1]) for dd in d])
    z=r*np.exp(1j*theta); z/=np.linalg.norm(z)
    th=np.angle(z); K2=A*np.sin(th[None,:]-th[:,None]); H2=1j*K2
    lam=float(np.real(np.vdot(z,H2@z)/np.vdot(z,z)))
    return edges,d,z,lam,H2,nd

def pmetrics(z,a,b):
    za,zb=z[a],z[b]; den=abs(za)**2+abs(zb)**2
    cl=abs(za*za+zb*zb)/den
    ratio=zb/za; ep=abs(ratio-1j); em=abs(ratio+1j)
    sg='+' if ep<=em else '-'; re=min(ep,em)
    amp=abs(abs(za)-abs(zb))/((abs(za)+abs(zb))/2)
    dp=float(np.angle(zb*np.conj(za)))
    pe=min(abs(np.angle(np.exp(1j*(dp-np.pi/2)))),abs(np.angle(np.exp(1j*(dp+np.pi/2)))))
    return cl,re,amp,pe,sg,dp

def matching_extremes(G):
    Gmax=nx.Graph(); Gmin=nx.Graph(); Gmax.add_nodes_from(G); Gmin.add_nodes_from(G)
    for u,v,d in G.edges(data=True):
        s=d['shared']
        Gmax.add_edge(u,v,weight=1000+s); Gmin.add_edge(u,v,weight=1000-s)
    ma=nx.max_weight_matching(Gmax,maxcardinality=True,weight='weight')
    mi=nx.max_weight_matching(Gmin,maxcardinality=True,weight='weight')
    return sum(G.edges[u,v]['shared'] for u,v in mi),sum(G.edges[u,v]['shared'] for u,v in ma)

def write_csv(path, rows):
    if not rows: return
    with path.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

def main():
    summary=[]; candidates=[]; matching_rows=[]; components=[]; amplitudes=[]
    for N in EVEN_NS:
        edges,d,z,lam,H,nd=build_theoretical_floor(N); M=len(edges)
        eigres=float(np.linalg.norm(H@z-lam*z)/np.linalg.norm(z))
        G=nx.Graph(); G.add_nodes_from(range(M)); cand=[]
        for a in range(M):
            for b in range(a+1,M):
                cl,re,amp,pe,sg,dp=pmetrics(z,a,b)
                if cl<TOL and re<5*TOL:
                    e1,e2=edges[a],edges[b]; shared=len(set(e1)&set(e2))
                    G.add_edge(a,b,shared=shared)
                    rec={'N':N,'edge_index_1':a,'edge_index_2':b,
                         'edge1_i':e1[0],'edge1_j':e1[1],'edge2_i':e2[0],'edge2_j':e2[1],
                         'd1':int(d[a]),'d2':int(d[b]),'closure_residual':cl,
                         'ratio_pm_i_error':re,'amplitude_rel_error':amp,'phase_error_rad':pe,
                         'canonical_sign_e2_over_e1':sg,'phase_diff_rad':dp,'shared_vertices':shared}
                    cand.append(rec); candidates.append(rec)
        match=nx.max_weight_matching(G,maxcardinality=True)
        mp=[]
        for u,v in sorted(tuple(sorted(x)) for x in match):
            cl,re,amp,pe,sg,dp=pmetrics(z,u,v); e1,e2=edges[u],edges[v]
            rec={'N':N,'edge_index_1':u,'edge_index_2':v,
                 'edge1_i':e1[0],'edge1_j':e1[1],'edge2_i':e2[0],'edge2_j':e2[1],
                 'd1':int(d[u]),'d2':int(d[v]),'closure_residual':cl,
                 'ratio_pm_i_error':re,'amplitude_rel_error':amp,'phase_error_rad':pe,
                 'canonical_sign_e2_over_e1':sg,'phase_diff_rad':dp,'shared_vertices':len(set(e1)&set(e2))}
            mp.append(rec); matching_rows.append(rec)
        cc=Counter(len(c) for c in nx.connected_components(G)) if cand else Counter()
        deg=Counter(dict(G.degree()).values()) if cand else Counter({0:M})
        minsh,maxsh=matching_extremes(G) if cand else (0,0)
        summary.append({
            'N':N,'M':M,'D':N//2,'lambda':lam,'norm_error':abs(float(np.vdot(z,z).real)-1),
            'sum_z_abs':abs(np.sum(z)),'sum_z2_abs':abs(np.sum(z*z)),'eigen_residual':eigres,
            'candidate_pairs':len(cand),'matching_pairs':len(mp),'matching_coverage':2*len(mp)/M,
            'candidate_graph_bipartite':bool(nx.is_bipartite(G)) if cand else True,
            'component_count':sum(cc.values()),'component_size4':cc.get(4,0),'component_size8':cc.get(8,0),'component_size16':cc.get(16,0),
            'degree0_vertices':deg.get(0,0),'degree2_vertices':deg.get(2,0),'degree4_vertices':deg.get(4,0),'degree8_vertices':deg.get(8,0),
            'min_shared_vertex_pairs_in_perfect_matching':minsh,'max_shared_vertex_pairs_in_perfect_matching':maxsh,
            'canonical_plus_in_selected_matching':sum(r['canonical_sign_e2_over_e1']=='+' for r in mp),
            'canonical_minus_in_selected_matching':sum(r['canonical_sign_e2_over_e1']=='-' for r in mp),
            'selected_matching_max_closure':max((r['closure_residual'] for r in mp),default=float('nan'))
        })
        for sz,count in sorted(cc.items()): components.append({'N':N,'component_size':sz,'count':count})
        for dd in range(1,N//2+1):
            vals=np.abs(z[d==dd]); amplitudes.append({'N':N,'distance_class':dd,'edge_count':len(vals),'amplitude':float(vals.mean()),'spread':float(vals.max()-vals.min())})
    write_csv(HERE/'even_N_summary.csv',summary)
    write_csv(HERE/'exact_pair_candidates.csv',candidates)
    write_csv(HERE/'selected_maximum_matchings.csv',matching_rows)
    write_csv(HERE/'component_size_counts.csv',components)
    write_csv(HERE/'distance_class_amplitudes.csv',amplitudes)
    # class-pair counts
    c=Counter((int(r['N']),int(r['d1']),int(r['d2']),int(r['shared_vertices'])) for r in candidates)
    write_csv(HERE/'candidate_distance_class_counts.csv',[{'N':k[0],'d1':k[1],'d2':k[2],'shared_vertices':k[3],'count':v} for k,v in sorted(c.items())])
    meta={'tolerance':TOL,'even_N':EVEN_NS,'exact_pair_N':[r['N'] for r in summary if r['candidate_pairs']>0],
          'perfect_partition_N':[r['N'] for r in summary if abs(r['matching_coverage']-1)<1e-12],
          'definition':'z_ij=r_d exp(i 2pi(i+j)/N), r_d from reduced Hermitian eigenproblem Q'}
    (HERE/'RUN_METADATA.json').write_text(json.dumps(meta,ensure_ascii=False,indent=2),encoding='utf-8')

if __name__=='__main__': main()
