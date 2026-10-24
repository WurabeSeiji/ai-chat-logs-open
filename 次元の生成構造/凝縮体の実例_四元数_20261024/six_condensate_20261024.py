"""
四元数 6 個（複素数 12 個）の凝縮体の実例 2026-10-24
構成（厳密）：
  三つ組 A（実部 −1/√2）：q = (−1 + n)/√2,  n ∈ { i, −½i+(√3/2)j, −½i−(√3/2)j }   （i-j 平面の正三角形）
  三つ組 B（実部 +1/√2）：q = (+1 + n)/√2,  n ∈ { ½j+(√3/2)k, −j, ½j−(√3/2)k }    （j-k 平面の正三角形）
  q = z1 + z2 j（z1 = a+bi, z2 = c+di）
判定量は condensate_example_20261024.py と同じ（閉塞、α、|k|、D、Δ3）。
あわせて「全ての対が |k| ≥ 1 になれるか」を、閉塞 Σq=0, Σq²=0 のもとでの最小角の最大化で N=4..8 について調べる。
"""
import numpy as np, itertools, json, sys
from scipy.optimize import minimize
def qmul(p, q):
    a1,b1,c1,d1 = p; a2,b2,c2,d2 = q
    return np.array([a1*a2-b1*b2-c1*c2-d1*d2, a1*b2+b1*a2+c1*d2-d1*c2,
                     a1*c2-b1*d2+c1*a2+d1*b2, a1*d2+b1*c2-c1*b2+d1*a2])
def qinv(q): return np.array([q[0],-q[1],-q[2],-q[3]])/np.dot(q,q)
r2, r3 = np.sqrt(2), np.sqrt(3)
nA = [np.array([0,1,0,0]), np.array([0,-.5,r3/2,0]), np.array([0,-.5,-r3/2,0])]
nB = [np.array([0,0,.5,r3/2]), np.array([0,0,-1,0]), np.array([0,0,.5,-r3/2])]
one = np.array([1.,0,0,0])
Q = [(-one+n)/r2 for n in nA] + [(one+n)/r2 for n in nB]
names = ["A1","A2","A3","B1","B2","B3"]

def analyze(Q):
    N=len(Q); out={}
    out["z"]=[[f"{q[0]:+.6f}{q[1]:+.6f}i", f"{q[2]:+.6f}{q[3]:+.6f}i"] for q in Q]
    zs=[complex(q[0],q[1]) for q in Q]+[complex(q[2],q[3]) for q in Q]
    out["sum_z_12"]=f"{sum(zs):.3e}"; out["sum_z2_12"]=f"{sum(z*z for z in zs):.6f}"
    out["sum_q"]=[round(x,12) for x in sum(Q)]; out["sum_q2"]=[round(x,12) for x in sum(qmul(q,q) for q in Q)]
    out["vector_part_singular_values"]=[round(x,6) for x in np.linalg.svd(np.array(Q)[:,1:],compute_uv=False)]
    out["pairs"]=[]
    for a,b in itertools.combinations(range(N),2):
        u=qmul(qinv(Q[a]),Q[b]); al=np.degrees(np.arccos(np.clip(u[0],-1,1))); k=np.tan(np.radians(al)/2)
        out["pairs"].append(dict(pair=[names[a],names[b]],alpha_deg=round(al,6),abs_k=round(k,6),D=round(1/k,6),
                                 bound=bool(k>1+1e-12)))
    out["triples"]=[]
    for a,b,c in itertools.combinations(range(N),3):
        u=qmul(qinv(Q[a]),Q[b]); v=qmul(qinv(Q[a]),Q[c]); com=qmul(u,v)-qmul(v,u)
        out["triples"].append(dict(triple=[names[a],names[b],names[c]],Delta3=round(0.25*np.dot(com,com),6)))
    return out

def maxmin_angle(N, trials=60, seed=11):
    P=list(itertools.combinations(range(N),2)); U=lambda x:x[:4*N].reshape(N,4)
    eq=lambda x:np.r_[U(x).sum(0),sum(qmul(q,q) for q in U(x)),np.sum(U(x)**2,1)-1]
    ineq=lambda x:np.array([x[-1]-U(x)[a]@U(x)[b] for a,b in P])
    rng=np.random.default_rng(seed); best=None
    for t in range(trials):
        Q0=rng.normal(size=(N,4)); Q0/=np.linalg.norm(Q0,axis=1)[:,None]
        r=minimize(lambda x:x[-1],np.r_[Q0.ravel(),.5],method="SLSQP",
                   constraints=[{"type":"eq","fun":eq},{"type":"ineq","fun":ineq}],options={"maxiter":500,"ftol":1e-13})
        if r.success and np.abs(eq(r.x)).max()<1e-9 and (best is None or r.fun<best): best=r.fun
    return float(np.degrees(np.arccos(best)))

if __name__=="__main__":
    res={"six_condensate":analyze(Q)}
    res["maxmin_angle_under_closures"]={str(N):round(maxmin_angle(N),4) for N in [4,5,6,7,8]}
    out=sys.argv[1] if len(sys.argv)>1 else "six_condensate_results_20261024.json"
    json.dump(res,open(out,"w"),ensure_ascii=False,indent=1)
    r=res["six_condensate"]
    print("z:",r["z"]); print("Σz(12)=",r["sum_z_12"],"Σz²(12)=",r["sum_z2_12"])
    print("Σq=",r["sum_q"],"Σq²=",r["sum_q2"],"sv=",r["vector_part_singular_values"])
    for p in r["pairs"]: print(p)
    print("Δ3 min/max:",min(t["Delta3"] for t in r["triples"]),max(t["Delta3"] for t in r["triples"]))
    for t in r["triples"]: print(t)
    print("max-min angle:",res["maxmin_angle_under_closures"])
