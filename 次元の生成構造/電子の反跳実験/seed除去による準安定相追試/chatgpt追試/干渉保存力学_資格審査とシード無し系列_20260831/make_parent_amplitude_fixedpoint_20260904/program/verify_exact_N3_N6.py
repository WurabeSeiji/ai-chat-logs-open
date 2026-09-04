#!/usr/bin/env python3
import json, math, sys
from pathlib import Path
import numpy as np
import sympy as sp

HERE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(HERE/'source_snapshot'))
import original_engine as eng
BASE=40260721; ITERS=400; BETA=0.5


def canonical_stage1(N):
    s=eng.LowRankSystem(N); rng=np.random.default_rng(BASE+1000*N)
    theta=rng.uniform(0,2*np.pi,s.m)
    for _ in range(ITERS):
        s.set_theta(theta)
        ev,EV=np.linalg.eig(s.J@s.G)
        idx=int(np.argmin(ev.imag))
        y=EV[:,idx].astype(complex)
        v=s.w(y)
        theta_new=np.angle(v)
        theta=np.angle((1-BETA)*np.exp(1j*theta)+BETA*np.exp(1j*theta_new))
    return s,v,y,float(-ev[idx].imag)


def quadrant_complex(z):
    a=float(np.angle(z))
    k=int(np.rint(a/(np.pi/2))) % 4
    return [sp.Integer(1),sp.I,sp.Integer(-1),-sp.I][k], abs(float(np.angle(np.exp(1j*(a-k*np.pi/2)))))


def exact_candidate(N):
    s,vn,yn,sig_num=canonical_stage1(N)
    exact_v=[]; phase_err=[]; amp2_exact=[]
    for z in vn:
        q,err=quadrant_complex(z)
        a2=sp.nsimplify(float(abs(z)**2), tolerance=1e-12, full=True)
        if not isinstance(a2, sp.Rational):
            raise RuntimeError(f'N={N}: non-rational amplitude square candidate {a2}')
        exact_v.append(q*sp.sqrt(a2))
        phase_err.append(err); amp2_exact.append(a2)
    v=sp.Matrix(exact_v)
    # exact phase-only W and K from exact quadrants
    ea,eb=eng.build_edges(N); M=len(ea)
    c=[]; ss=[]
    for z in exact_v:
        # q is z/|z| in {+-1,+-i}
        q=sp.simplify(z/sp.sqrt(sp.conjugate(z)*z))
        c.append(sp.re(q)); ss.append(sp.im(q))
    W=sp.zeros(M,2*N)
    for e,(a,b) in enumerate(zip(ea,eb)):
        W[e,int(a)]=c[e]; W[e,int(b)]=c[e]
        W[e,N+int(a)]=ss[e]; W[e,N+int(b)]=ss[e]
    J=sp.zeros(2*N,2*N)
    for i in range(N): J[i,N+i]=1; J[N+i,i]=-1
    G=sp.simplify(W.T*W)
    K=sp.simplify(W*J*W.T)
    iv=sp.simplify(sp.I*K*v)
    # choose eigenvalue from a nonzero component
    mu=None
    for i in range(M):
        if v[i]!=0:
            mu=sp.simplify(iv[i]/v[i]); break
    eig_res=sp.simplify(iv-mu*v)
    # sigma sign: JG y = -i sigma y corresponds iK v = sigma v for this convention
    sigma=sp.simplify(mu)
    y=sp.simplify(sp.I/sigma * J*W.T*v)
    ynorm2=sp.simplify((sp.conjugate(y).T*y)[0])
    wy_res=sp.simplify(W*y-v)
    jg_res=sp.simplify(J*G*y + sp.I*sigma*y)
    norm2=sp.simplify((sp.conjugate(v).T*v)[0])
    r2=sp.simplify(norm2/M)
    Nr2=sp.simplify(N*r2)
    return {
        'N':N,'M':M,'sigma':str(sigma),'sigma2':str(sp.simplify(sigma**2)),
        'norm2':str(norm2),'r2':str(r2),'Nr2':str(Nr2),'y_norm2':str(ynorm2),
        'max_phase_quadrant_error_rad':max(phase_err),
        'eig_residual_exact_zero':all(sp.simplify(x)==0 for x in eig_res),
        'Wy_minus_v_exact_zero':all(sp.simplify(x)==0 for x in wy_res),
        'JG_eigen_residual_exact_zero':all(sp.simplify(x)==0 for x in jg_res),
        'y_norm_exact_one':sp.simplify(ynorm2-1)==0,
        'amp2_by_edge':[str(x) for x in amp2_exact],
        'phase_quadrants':[str(sp.simplify(z/sp.sqrt(sp.conjugate(z)*z))) for z in exact_v]
    }

out=[]
for N in range(3,7):
    d=exact_candidate(N); out.append(d)
    print(N,d['sigma2'],d['norm2'],d['r2'],d['Nr2'],d['eig_residual_exact_zero'],d['Wy_minus_v_exact_zero'],d['JG_eigen_residual_exact_zero'],d['y_norm_exact_one'])
res=HERE/'results'
(res/'symbolic_exact_verification_N3_N6.json').write_text(json.dumps(out,indent=2,ensure_ascii=False),encoding='utf-8')
with open(res/'symbolic_exact_verification_N3_N6.md','w',encoding='utf-8') as f:
    f.write('# N=3..6 段階1 exact 検証\n\n')
    f.write('倍精度の canonical stage1 最終状態から、位相を最寄りの {1,i,-1,-i}、振幅二乗を有理数候補へ復元し、その候補を SymPy の厳密演算で元の W,J,G,K 方程式へ代入した。したがって以下の residual=0 は浮動小数点近似ではなく記号的恒等式。\n\n')
    f.write('|N|sigma^2|norm^2|r^2|N r^2|iKv=sigma v|Wy=v|JG y=-i sigma y|norm y=1|\n|---:|---:|---:|---:|---:|:---:|:---:|:---:|:---:|\n')
    for d in out:
        f.write(f"|{d['N']}|{d['sigma2']}|{d['norm2']}|{d['r2']}|{d['Nr2']}|{d['eig_residual_exact_zero']}|{d['Wy_minus_v_exact_zero']}|{d['JG_eigen_residual_exact_zero']}|{d['y_norm_exact_one']}|\n")
    f.write('\n注意: これは候補の exact 検証であり、固定点枝の一意性証明ではない。枝の一意性/多重性は別の seed sweep で検査する。\n')
