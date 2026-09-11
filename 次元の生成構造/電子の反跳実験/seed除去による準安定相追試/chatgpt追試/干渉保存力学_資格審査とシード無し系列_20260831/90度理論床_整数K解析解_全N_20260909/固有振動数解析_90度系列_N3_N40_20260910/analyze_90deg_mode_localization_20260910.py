import os, math, csv
import numpy as np

CANON='/mnt/data/run_N3_N40_stage123_v1.py'
OUT='/mnt/data/90deg_full_eigenspectrum_20260910'

src=open(CANON,encoding='utf-8').read()
g={'__name__':'canonical','__file__':CANON}
exec(compile(src.split('rows=[]; summaries=[]')[0],CANON,'exec'),g)
adjacency,one_step=g['adjacency'],g['one_step']

def integer_K(labels,A):
    diff=(labels[None,:]-labels[:,None])%4
    s=np.zeros_like(diff,dtype=float)
    s[diff==1]=1.0; s[diff==3]=-1.0
    return A*s

def releq_residual(z,A,N):
    w=one_step(z,A,N)
    ip=np.vdot(z,w)
    return float(np.linalg.norm(w-np.exp(1j*np.angle(ip))*z)/np.linalg.norm(z))

def build_analytic_floor(N,max_iter=500):
    ea,eb=np.triu_indices(N,k=1); A=adjacency(N)
    labels=((ea+eb)%4).astype(int)
    best=None; seen={}
    for it in range(max_iter):
        K=integer_K(labels,A)
        ew,V=np.linalg.eigh(1j*K)
        v=V[:,-1]; z=v/np.linalg.norm(v)
        res=releq_residual(z,A,N)
        if best is None or res<best[0]: best=(res,z.copy(),labels.copy())
        newlab=(np.round(np.degrees(np.angle(v))%360/90).astype(int))%4
        cand=None; agmax=-1
        for sh in range(4):
            c=(newlab+sh)%4; ag=int(np.sum(c==labels))
            if ag>agmax: agmax=ag; cand=c
        key=tuple(cand.tolist())
        if np.array_equal(cand,labels) or key in seen: break
        seen[key]=it; labels=cand
    return best[1],best[2],best[0]

def cluster_indices(vals,tol=1e-8):
    # positive eigenvalues only; return groups of original indices
    idx=[i for i,x in enumerate(vals) if x>1e-9]
    groups=[]
    for i in idx:
        x=vals[i]
        if not groups:
            groups=[[i]]
        else:
            ref=np.mean([vals[j] for j in groups[-1]])
            if abs(x-ref)<=tol*max(1.0,abs(x),abs(ref)):
                groups[-1].append(i)
            else:
                groups.append([i])
    return groups

rows=[]
for N in range(3,41):
    z,labels,res=build_analytic_floor(N)
    A=adjacency(N); K=integer_K(labels,A)
    vals,V=np.linalg.eigh(1j*K)
    ea,eb=np.triu_indices(N,k=1)
    M=len(ea)
    groups=cluster_indices(vals)
    for gi,inds in enumerate(groups,1):
        w=float(np.mean(vals[inds])); mult=len(inds)
        U=V[:,inds]
        # Basis-invariant diagonal of spectral projector P=UU^H
        pedge=np.sum(np.abs(U)**2,axis=1).real  # sum = mult
        edge_w=pedge/mult
        edge_ipr=float(np.sum(edge_w**2))
        edge_pr=float(1.0/edge_ipr)
        edge_pr_frac=float(edge_pr/M)
        edge_max=float(np.max(edge_w))
        edge_min=float(np.min(edge_w))
        edge_cv=float(np.std(edge_w)/np.mean(edge_w))
        projector_diag_max=float(np.max(pedge)) # max projection fraction of an edge basis state into subspace

        # vertex support: sum incident edge projector weights; normalize to sum 1
        pv=np.zeros(N,float)
        for e,(a,b) in enumerate(zip(ea,eb)):
            pv[a]+=pedge[e]; pv[b]+=pedge[e]
        vertex_w=pv/(2.0*mult)
        vertex_ipr=float(np.sum(vertex_w**2))
        vertex_pr=float(1.0/vertex_ipr)
        vertex_pr_frac=float(vertex_pr/N)
        vertex_max=float(np.max(vertex_w))
        vertex_min=float(np.min(vertex_w))
        vertex_cv=float(np.std(vertex_w)/np.mean(vertex_w))

        # Z4 sector support, basis-invariant
        qweights=[]
        for q in range(4):
            qweights.append(float(np.sum(pedge[labels==q])/mult))
        q_ipr=float(np.sum(np.array(qweights)**2))
        q_pr=float(1.0/q_ipr) if q_ipr>0 else 0.0

        # For non-degenerate modes, ordinary vector IPR is meaningful too.
        if mult==1:
            u=U[:,0]
            wi=np.abs(u)**2; wi/=wi.sum()
            single_edge_pr=float(1.0/np.sum(wi**2))
            vi=np.zeros(N,float)
            for e,(a,b) in enumerate(zip(ea,eb)):
                vi[a]+=wi[e]; vi[b]+=wi[e]
            vi/=vi.sum()
            single_vertex_pr=float(1.0/np.sum(vi**2))
        else:
            single_edge_pr=float('nan'); single_vertex_pr=float('nan')

        rows.append(dict(N=N,M=M,mode_family=gi,omega=w,multiplicity=mult,
                         edge_subspace_PR=edge_pr,edge_PR_fraction=edge_pr_frac,
                         edge_weight_cv=edge_cv,edge_weight_min=edge_min,edge_weight_max=edge_max,
                         max_edge_basis_projection=projector_diag_max,
                         vertex_subspace_PR=vertex_pr,vertex_PR_fraction=vertex_pr_frac,
                         vertex_weight_cv=vertex_cv,vertex_weight_min=vertex_min,vertex_weight_max=vertex_max,
                         z4_PR=q_pr,z4_q0=qweights[0],z4_q1=qweights[1],z4_q2=qweights[2],z4_q3=qweights[3],
                         single_edge_PR=single_edge_pr,single_vertex_PR=single_vertex_pr,
                         releq_residual=res))
    print('done',N,'families',len(groups),flush=True)

path=os.path.join(OUT,'mode_localization_basis_invariant_N3_N40.csv')
with open(path,'w',newline='',encoding='utf-8') as f:
    wr=csv.DictWriter(f,fieldnames=list(rows[0].keys())); wr.writeheader(); wr.writerows(rows)

# Compact findings
print('\nSAMPLE')
for r in rows:
    if r['N'] in (3,4,5,6,7,8,12,20,40):
        print(r['N'],r['mode_family'],'w',round(r['omega'],6),'mult',r['multiplicity'],
              'edgePR/M',round(r['edge_PR_fraction'],6),'vertexPR/N',round(r['vertex_PR_fraction'],6),
              'edgeCV',round(r['edge_weight_cv'],6),'vCV',round(r['vertex_weight_cv'],6),
              'maxProj',round(r['max_edge_basis_projection'],6),'z4PR',round(r['z4_PR'],6))
