
from pathlib import Path
import importlib.util, numpy as np, pandas as pd, json
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
ENGINE=HERE/"run_n_scaling_lowrank_v1_no_sigma_norm.py"

spec=importlib.util.spec_from_file_location("eng",ENGINE)
eng=importlib.util.module_from_spec(spec); spec.loader.exec_module(eng)

STEPS=5000

def centered_gram(z,n,edges):
    D2=np.zeros((n,n),complex)
    for val,(i,j) in zip(z*z,edges):
        D2[i,j]=D2[j,i]=val
    J=np.eye(n)-np.ones((n,n))/n
    return -0.5*J@D2@J

def cluster_triplets(vals,tol):
    centers=[]; members=[]
    for ridx,v in enumerate(vals):
        hit=None
        for ci,c in enumerate(centers):
            if np.max(np.abs(v-c))<tol:
                hit=ci; break
        if hit is None:
            centers.append(v.copy()); members.append([ridx])
        else:
            members[hit].append(ridx)
            centers[hit]=vals[members[hit]].mean(axis=0)
    return centers,members

def runN(N):
    M=N*(N-1)//2
    SEED=40260722+1000*N
    sys=eng.LowRankSystem(N)
    rng=np.random.default_rng(SEED)
    v,res,sig=eng.make_parent(sys,rng,iters=1200,beta=0.5,tol=1e-12,restarts=3)
    Z=v.copy()
    wp=rng.normal(size=M)
    edges=[(i,j) for i in range(N) for j in range(i+1,N)]

    Q,_=np.linalg.qr(np.column_stack([v.real,v.imag]))
    Q=Q[:,:min(2,Q.shape[1])]

    global_rows=[]; axes_rows=[]; phase_rows=[]
    for t in range(STEPS+1):
        a=Z.real.copy(); b=Z.imag.copy()
        a2=a*a; b2=b*b; ab=a*b; r2=a2+b2
        Pa=Q@(Q.T@a); Pb=Q@(Q.T@b)
        Hpar=float(np.dot(Pa,Pa)+np.dot(Pb,Pb))
        Htot=float(np.vdot(Z,Z).real)
        Hperp=max(0.0,Htot-Hpar)
        B=centered_gram(Z,N,edges)
        s=np.linalg.svd(B,compute_uv=False)
        rank=int(np.sum(s>max(float(s[0]),1.0)*1e-10))
        tak=np.sqrt(np.maximum(s,0.0))
        global_rows.append([t,Htot,Hpar,Hperp,np.sqrt(Hperp),float(abs(Z@Z)),
                            float(a2.sum()),float(b2.sum()),float(ab.sum()),
                            float(r2.min()),float(r2.max()),rank])
        axes_rows.append([t]+list(map(float,tak)))
        for k,(i,j) in enumerate(edges):
            z2=Z[k]*Z[k]
            phase_rows.append([t,k,i+1,j+1,float(a[k]),float(b[k]),float(a2[k]),float(b2[k]),
                               float(ab[k]),float(r2[k]),float(np.angle(Z[k])),
                               float(z2.real),float(z2.imag),float(np.angle(z2)/np.pi)])
        if t<STEPS:
            sys.set_state(Z)  # FIX4
            se,wp=sys.sigma_max_power(wp)
            Z=sys.linear_rotation_step(Z,se)

    glob=pd.DataFrame(global_rows,columns=["step","H_total","H_parallel","H_perp","A_perp","abs_ZtZ",
        "sum_a2","sum_b2","sum_ab","r2_min","r2_max","simplex_rank"])
    glob.to_csv(HERE/f"N{N}_global_summary.csv",index=False)
    axes=pd.DataFrame(axes_rows,columns=["step"]+[f"takagi_r{k+1}" for k in range(N)])
    axes.to_csv(HERE/f"N{N}_takagi_axes.csv",index=False)
    ph=pd.DataFrame(phase_rows,columns=["step","edge_index","i","j","a","b","a2","b2","ab","r2","theta","z2_re","z2_im","z2_phase_pi"])
    ph.to_csv(HERE/f"N{N}_all_steps_long.csv",index=False)
    gf=ph[ph.step==STEPS].copy()
    gf.to_csv(HERE/f"N{N}_step5000_final_edges.csv",index=False)

    # selected-step class counts
    sels=[0,1,10,50,100,150,200,300,500,750,1000,1500,2000,2500,3000,3500,4000,4500,5000]
    clrows=[]
    for t in sels:
        g=ph[ph.step==t].copy()
        vals=np.column_stack([g.a2/g.r2,g.b2/g.r2,g.ab/g.r2])
        for tol in [1e-2,1e-3,1e-4,1e-6,1e-8,1e-10]:
            centers,members=cluster_triplets(vals,tol)
            clrows.append([t,tol,len(centers),";".join(map(str,sorted([len(x) for x in members],reverse=True)))])
    cl=pd.DataFrame(clrows,columns=["step","tol","triplet_cluster_count","cluster_sizes_desc"])
    cl.to_csv(HERE/f"N{N}_triplet_cluster_counts.csv",index=False)

    # final classes at 1e-8 and 1e-6
    for tol in [1e-6,1e-8]:
        vals=np.column_stack([gf.a2/gf.r2,gf.b2/gf.r2,gf.ab/gf.r2])
        centers,members=cluster_triplets(vals,tol)
        out=[]
        for ci,inds in enumerate(members):
            gg=gf.iloc[inds]
            c=vals[inds].mean(axis=0)
            # circular mean of theta modulo pi
            u=np.mean(np.exp(1j*2*gg.theta.to_numpy()))
            th=((np.angle(u)/2)%np.pi)/np.pi
            out.append([ci,len(inds),c[0],c[1],c[2],th,
                        ",".join(f"{int(i)}-{int(j)}" for i,j in zip(gg.i,gg.j))])
        pd.DataFrame(out,columns=["class_id","count","a2n","b2n","abn","theta_over_pi_mod1","edges"]).to_csv(
            HERE/f"N{N}_final_classes_tol{tol:g}.csv",index=False)

    # pairwise phase differences at final
    th=gf.theta.to_numpy()
    diffs=[]
    for i in range(len(th)):
        for j in range(i+1,len(th)):
            d=abs(((th[i]-th[j]+np.pi/2)%np.pi)-np.pi/2)/np.pi
            diffs.append(d)
    diffs=np.array(diffs)
    dcent=[]; dcount=[]
    for val in sorted(diffs):
        hit=None
        for ci,c in enumerate(dcent):
            if abs(val-c)<1e-6:
                hit=ci; break
        if hit is None:
            dcent.append(val); dcount.append(1)
        else:
            dcount[hit]+=1
    pd.DataFrame({"abs_delta_theta_over_pi":dcent,"count":dcount}).to_csv(
        HERE/f"N{N}_pairwise_phase_difference_classes_tol1e-6.csv",index=False)

    # time milestones
    glob["r2_mean"]=(glob.r2_min+glob.r2_max)/2
    glob["r2_rel_spread"]=(glob.r2_max-glob.r2_min)/glob.r2_mean
    milestones=[]
    finalh=float(glob.H_perp.iloc[-1])
    if finalh>0:
        for frac in [0.5,0.9,0.95,0.99]:
            idx=np.where(glob.H_perp.to_numpy()>=frac*finalh)[0]
            milestones.append([f"H_perp >= {int(frac*100)}% final",int(glob.step.iloc[idx[0]]) if len(idx) else None])
    for thr in [1e-1,1e-2,1e-3,1e-4,1e-6,1e-8,1e-10]:
        arr=glob.r2_rel_spread.to_numpy()
        first=None
        for ii in range(len(arr)):
            if np.all(arr[ii:]<thr):
                first=int(glob.step.iloc[ii]); break
        milestones.append([f"r2 relative spread < {thr:g} forever",first])
    pd.DataFrame(milestones,columns=["metric","step"]).to_csv(HERE/f"N{N}_time_milestones.csv",index=False)

    # plots
    plt.figure(figsize=(8,5))
    plt.semilogy(glob.step,np.maximum(glob.H_perp,1e-30),label="H_perp")
    plt.semilogy(glob.step,np.maximum(glob.r2_rel_spread,1e-16),label="edge |z|^2 relative spread")
    plt.xlabel("step"); plt.ylabel("log scale"); plt.title(f"N={N}: decompactification and equalization")
    plt.legend(); plt.tight_layout(); plt.savefig(HERE/f"N{N}_inflation_equalization.png",dpi=180); plt.close()

    plt.figure(figsize=(8,5))
    for k in range(N-1):
        plt.plot(axes.step,axes[f"takagi_r{k+1}"])
    plt.xlabel("step"); plt.ylabel("Takagi axis scale"); plt.title(f"N={N}: non-null simplex axes")
    plt.tight_layout(); plt.savefig(HERE/f"N{N}_takagi_axes.png",dpi=180); plt.close()

    r2=gf.r2.to_numpy()
    plt.figure(figsize=(6,6))
    plt.scatter(gf.z2_re/r2,gf.z2_im/r2,s=40)
    for _,rr in gf.iterrows():
        plt.annotate(f"{int(rr.i)}-{int(rr.j)}",(rr.z2_re/rr.r2,rr.z2_im/rr.r2),fontsize=7)
    ang=np.linspace(0,2*np.pi,400)
    plt.plot(np.cos(ang),np.sin(ang),linewidth=0.8)
    plt.xlabel("Re(z^2)/|z|^2"); plt.ylabel("Im(z^2)/|z|^2")
    plt.title(f"N={N} step 5000: normalized complex squared distances")
    plt.axis("equal"); plt.tight_layout(); plt.savefig(HERE/f"N{N}_final_z2_phase_circle.png",dpi=180); plt.close()

    plt.figure(figsize=(8,5))
    for tol in [1e-2,1e-4,1e-6,1e-8]:
        dd=cl[cl.tol==tol]
        plt.plot(dd.step,dd.triplet_cluster_count,label=f"tol={tol:g}")
    plt.xlabel("step"); plt.ylabel("class count"); plt.title(f"N={N}: (a2,b2,ab) class count")
    plt.legend(); plt.tight_layout(); plt.savefig(HERE/f"N{N}_class_count.png",dpi=180); plt.close()

    summary={
        "N":N,"M":M,"seed":SEED,"steps":STEPS,"parent_residual":float(res),
        "final_rank":int(glob.simplex_rank.iloc[-1]),
        "final_H_perp":float(glob.H_perp.iloc[-1]),
        "final_r2_min":float(glob.r2_min.iloc[-1]),"final_r2_max":float(glob.r2_max.iloc[-1]),
        "final_sum_a2":float(glob.sum_a2.iloc[-1]),"final_sum_b2":float(glob.sum_b2.iloc[-1]),"final_sum_ab":float(glob.sum_ab.iloc[-1]),
        "final_class_count_tol1e-6":int(pd.read_csv(HERE/f"N{N}_final_classes_tol1e-06.csv").shape[0]),
        "final_class_count_tol1e-8":int(pd.read_csv(HERE/f"N{N}_final_classes_tol1e-08.csv").shape[0]),
    }
    (HERE/f"N{N}_summary.json").write_text(json.dumps(summary,indent=2),encoding="utf-8")
    return summary

print(json.dumps([runN(10),runN(11)],indent=2))
