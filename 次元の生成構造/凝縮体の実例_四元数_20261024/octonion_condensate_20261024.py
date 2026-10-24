"""
八元数の凝縮体の実例 2026-10-24
八元数 o = (a, b)（Cayley–Dickson、a, b は四元数）、a = z1 + z2 j, b = z3 + z4 j → 一つの八元数 = 複素数 4 個
積 (a,b)(c,d) = (ac − d̄ b, da + b c̄)（論文 §11 の γ = −1、通常の八元数）
構成（四元数の 4 個・6 個の例と同じ族）：
  群 A（実部 −1/√2）：虚部が e1..e4 の 4 次元部分空間の正 4 単体（5 個、和 0、|v| = 1/√2）。単体は (±√5/4,±√5/4,±√5/4,−1/4) と (0,0,0,1)
  群 B（実部 +1/√2）：虚部が e4..e7 の 4 次元部分空間の正 4 単体（5 個）
  → 八元数 10 個 = 複素数 40 個。虚部は 7 方向すべてに張る。
判定量：閉塞 Σo=0, Σo²=0、相対位相の角 α（u = o_a⁻¹ o_b）、|k|、D、
        三体の非可換性 Δ3 = ¼|uv − vu|²、四体の非結合性 |(uv)w − u(vw)|（観測者 o_a から見た相対位相 u, v, w）
実行：python octonion_condensate_20261024.py         （実例の計算）
      python octonion_condensate_20261024.py --maxmin （閉塞のもとでの最小角の最大値 N=6,8,10 を追記）
"""
import numpy as np, itertools, json, sys

def qmul(p,q):
    a1,b1,c1,d1=p;a2,b2,c2,d2=q
    return np.array([a1*a2-b1*b2-c1*c2-d1*d2,a1*b2+b1*a2+c1*d2-d1*c2,a1*c2-b1*d2+c1*a2+d1*b2,a1*d2+b1*c2-c1*b2+d1*a2])
def qc(q): return np.array([q[0],-q[1],-q[2],-q[3]])
def omul(x,y):
    a,b=x[:4],x[4:]; c,d=y[:4],y[4:]
    return np.r_[qmul(a,c)-qmul(qc(d),b), qmul(d,a)+qmul(b,qc(c))]
def oconj(x): return np.r_[x[0],-x[1:]]
def oinv(x): return oconj(x)/np.dot(x,x)

def build():
    # 厳密な正 4 単体（4 次元、単位長）：(±√5/4, ±√5/4, ±√5/4, −1/4)（積が +1 の 4 通り）と (0,0,0,1)
    a=np.sqrt(5)/4
    S=np.array([[ a, a, a,-.25],[ a,-a,-a,-.25],[-a, a,-a,-.25],[-a,-a, a,-.25],[0,0,0,1.]])
    SA=S                                         # 群 A：(e1,e2,e3,e4) に S の第 1..4 成分
    SB=S[:,[3,0,1,2]][:,[1,0,2,3]]               # 群 B：(e4,e5,e6,e7) に (S1, S4, S2, S3)
    O=[]
    for v in SA:                                 # 群 A：e1..e4
        o=np.zeros(8); o[0]=-1/np.sqrt(2); o[1:5]=v/np.sqrt(2); O.append(o)
    for v in SB:                                 # 群 B：e4..e7（e4 を共有）
        o=np.zeros(8); o[0]=1/np.sqrt(2); o[4:8]=v/np.sqrt(2); O.append(o)
    return O

def analyze(O):
    N=len(O); names=[f"A{i+1}" for i in range(N//2)]+[f"B{i+1}" for i in range(N//2)]
    r={}
    zs=[]; r["z"]=[]
    for o in O:
        z=[complex(o[0],o[1]),complex(o[2],o[3]),complex(o[4],o[5]),complex(o[6],o[7])]
        zs+=z; r["z"].append([f"{w.real:+.6f}{w.imag:+.6f}i" for w in z])
    r["sum_z"]=f"{sum(zs):.3e}"; r["sum_z2"]=f"{sum(w*w for w in zs):.6f}"
    r["sum_o"]=[round(x,12)+0.0 for x in sum(O)]
    r["sum_o2"]=[round(x,12)+0.0 for x in sum(omul(o,o) for o in O)]
    r["vector_part_rank_singular_values"]=[round(x,6) for x in np.linalg.svd(np.array(O)[:,1:],compute_uv=False)]
    P=[]
    for a,b in itertools.combinations(range(N),2):
        u=omul(oinv(O[a]),O[b]); al=np.degrees(np.arccos(np.clip(u[0],-1,1))); k=np.tan(np.radians(al)/2)
        P.append(dict(pair=[names[a],names[b]],alpha_deg=round(al,6),abs_k=round(k,6),D=round(1/k,6),bound=bool(k>1+1e-12)))
    r["pairs"]=P
    d3=[]
    for a,b,c in itertools.combinations(range(N),3):
        u=omul(oinv(O[a]),O[b]); v=omul(oinv(O[a]),O[c]); com=omul(u,v)-omul(v,u); d3.append(0.25*np.dot(com,com))
    r["Delta3_min_max_count_pos"]=[round(min(d3),6),round(max(d3),6),int(sum(1 for x in d3 if x>1e-12)),len(d3)]
    assoc=[]
    for a in range(N):
        others=[i for i in range(N) if i!=a]
        for b,c,d in itertools.combinations(others,3):
            u=omul(oinv(O[a]),O[b]); v=omul(oinv(O[a]),O[c]); w=omul(oinv(O[a]),O[d])
            A=omul(omul(u,v),w)-omul(u,omul(v,w)); assoc.append(np.linalg.norm(A))
    r["associator_min_max_count_pos"]=[round(min(assoc),6),round(max(assoc),6),int(sum(1 for x in assoc if x>1e-12)),len(assoc)]
    # 群内だけ・群間だけの集計
    within=[p for p in P if p["pair"][0][0]==p["pair"][1][0]]; cross=[p for p in P if p["pair"][0][0]!=p["pair"][1][0]]
    r["within_alpha_range"]=[min(p["alpha_deg"] for p in within),max(p["alpha_deg"] for p in within)]
    r["cross_alpha_range"]=[min(p["alpha_deg"] for p in cross),max(p["alpha_deg"] for p in cross)]
    r["cross_all_bound"]=all(p["bound"] for p in cross)
    return r

if __name__=="__main__" and "--maxmin" not in sys.argv:
    O=build()
    res={"octonion_condensate_N10":analyze(O), "octonions":[[round(x,6)+0.0 for x in o] for o in O]}
    json.dump(res,open("octonion_condensate_results_20261024.json","w"),ensure_ascii=False,indent=1)
    r=res["octonion_condensate_N10"]
    for key in ["sum_z","sum_z2","sum_o","sum_o2","vector_part_rank_singular_values","within_alpha_range","cross_alpha_range","cross_all_bound","Delta3_min_max_count_pos","associator_min_max_count_pos"]:
        print(key, r[key])
    for row in r["z"]: print(row)

# ---- 閉塞のもとでの最小角の最大値（8 次元、N = 6, 8, 10）。全対が結合できるかの上限
def maxmin_angle(N, dim=8, trials=25, seed=3):
    from scipy.optimize import minimize
    P=list(itertools.combinations(range(N),2)); U=lambda x:x[:dim*N].reshape(N,dim)
    sq=lambda o:np.r_[o[0]**2-np.dot(o[1:],o[1:]),2*o[0]*o[1:]]     # 八元数の真の自乗（交代性より次元によらない形）
    eq=lambda x:np.r_[U(x).sum(0),sum(sq(o) for o in U(x)),np.sum(U(x)**2,1)-1]
    ineq=lambda x:np.array([x[-1]-U(x)[a]@U(x)[b] for a,b in P])
    rng=np.random.default_rng(seed);best=None
    for t in range(trials):
        X=rng.normal(size=(N,dim));X/=np.linalg.norm(X,axis=1)[:,None]
        r=minimize(lambda x:x[-1],np.r_[X.ravel(),.5],method="SLSQP",constraints=[{"type":"eq","fun":eq},{"type":"ineq","fun":ineq}],options={"maxiter":400,"ftol":1e-12})
        if r.success and np.abs(eq(r.x)).max()<1e-9 and (best is None or r.fun<best): best=r.fun
    return round(float(np.degrees(np.arccos(best))),4)
if __name__=="__main__" and "--maxmin" in sys.argv:
    mm={str(N):maxmin_angle(N) for N in [6,8,10]}
    res=json.load(open("octonion_condensate_results_20261024.json")); res["maxmin_angle_under_closures_8D"]=mm
    json.dump(res,open("octonion_condensate_results_20261024.json","w"),ensure_ascii=False,indent=1); print(mm)
