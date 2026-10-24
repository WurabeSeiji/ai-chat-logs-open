"""
凝縮体が存在する {z_n} の実例（四元数段）2026-10-24
q = z1 + z2 j （z1 = a+bi, z2 = c+di → q = a + b i + c j + d k）
判定量：
  閉塞    Σq = 0, Σq² = 0（真の自乗）
  相対位相 u = q_a^{-1} q_b, 角 α = arccos Re(u)（等振幅）
  結合係数 |k| = tan(α/2), 視野角距離 D = cot(α/2)
  三体判別式 Δ3 = 1/4 |uv - vu|^2 （u = q1^{-1}q2, v = q1^{-1}q3）
"""
import numpy as np, itertools, json, sys

def qmul(p, q):
    a1,b1,c1,d1 = p; a2,b2,c2,d2 = q
    return np.array([a1*a2-b1*b2-c1*c2-d1*d2,
                     a1*b2+b1*a2+c1*d2-d1*c2,
                     a1*c2-b1*d2+c1*a2+d1*b2,
                     a1*d2+b1*c2-c1*b2+d1*a2])
def qconj(q): return np.array([q[0],-q[1],-q[2],-q[3]])
def qinv(q):  return qconj(q)/np.dot(q,q)
def from_z(z1, z2): return np.array([z1.real, z1.imag, z2.real, z2.imag])

s = 1/np.sqrt(2)
def make_set(theta):
    """凝縮体候補：二組のホライズン上のペアを結合させたもの。
       m = i, n = cos(theta) i + sin(theta) j。theta=90° で m⊥n（ℍ）、theta=0 で m∥n（ℂ）"""
    n = np.array([0, np.cos(theta), np.sin(theta), 0])
    m = np.array([0, 1, 0, 0])
    one = np.array([1,0,0,0])
    return [s*(one+m), s*(one-m), s*(-one+n), s*(-one-n)]

def analyze(Q, label):
    N = len(Q)
    sumq  = sum(Q)
    sumq2 = sum(qmul(q,q) for q in Q)
    pairs = []
    for a,b in itertools.combinations(range(N),2):
        u = qmul(qinv(Q[a]), Q[b])
        alpha = np.degrees(np.arccos(np.clip(u[0],-1,1)))
        k = np.tan(np.radians(alpha)/2)
        pairs.append(dict(pair=[a+1,b+1], alpha_deg=round(alpha,6),
                          abs_k=(round(k,6) if alpha<179.999 else "inf"),
                          D=(round(1/k,6) if k>1e-12 else "inf")))
    tri = []
    for a,b,c in itertools.combinations(range(N),3):
        u = qmul(qinv(Q[a]), Q[b]); v = qmul(qinv(Q[a]), Q[c])
        com = qmul(u,v)-qmul(v,u)
        tri.append(dict(triple=[a+1,b+1,c+1], Delta3=round(0.25*np.dot(com,com),6)))
    return dict(label=label,
                q=[[round(x,6) for x in q] for q in Q],
                z=[[f"{q[0]:+.6f}{q[1]:+.6f}i", f"{q[2]:+.6f}{q[3]:+.6f}i"] for q in Q],
                sum_q=[round(x,12) for x in sumq],
                sum_q2=[round(x,12) for x in sumq2],
                pairs=pairs, triples=tri)

res = {}
res["H_condensate_theta90"] = analyze(make_set(np.pi/2), "m ⊥ n（四元数）")
res["C_degenerate_theta0"]  = analyze(make_set(0.0),    "m ∥ n（複素数）")
res["intermediate_theta45"] = analyze(make_set(np.pi/4), "中間 45°")
out = sys.argv[1] if len(sys.argv)>1 else "condensate_example_results_20261024.json"
json.dump(res, open(out,"w"), ensure_ascii=False, indent=1)
for key,r in res.items():
    print("==", key, r["label"])
    print(" z =", r["z"])
    print(" Σq =", r["sum_q"], " Σq² =", r["sum_q2"])
    for p in r["pairs"]: print("  ", p)
    for t in r["triples"]: print("  ", t)
