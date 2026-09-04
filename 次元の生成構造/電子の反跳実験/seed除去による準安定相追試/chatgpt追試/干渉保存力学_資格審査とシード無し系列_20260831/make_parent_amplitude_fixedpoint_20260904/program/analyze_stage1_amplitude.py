#!/usr/bin/env python3
import csv, json, math, sys
from pathlib import Path
import numpy as np
from sympy import nsimplify, sqrt

HERE = Path(__file__).resolve().parents[1]
SRC = HERE / 'source_snapshot'
sys.path.insert(0, str(SRC))
import original_engine as eng

NMIN, NMAX = 3, 16
ITERS = 400
BETA = 0.5
SEED_BASE = 40260721


def old_residual(sys_lr, v):
    sys_lr.set_theta(np.angle(v))
    kv = sys_lr.kmatvec(v)
    mu = float(np.real(np.conj(v) @ (1j * kv)))
    return float(np.linalg.norm(1j * kv - mu * v)), mu


def corrected_residual(sys_lr, v):
    sys_lr.set_theta(np.angle(v))
    kv = sys_lr.kmatvec(v)
    iv = 1j * kv
    mu = np.vdot(v, iv) / np.vdot(v, v)
    res = np.linalg.norm(iv - mu * v) / np.linalg.norm(v)
    return float(res), float(np.real(mu))


def stage1_trace(N):
    sys_lr = eng.LowRankSystem(N)
    rng = np.random.default_rng(SEED_BASE + 1000 * N)
    theta = rng.uniform(0.0, 2.0*np.pi, sys_lr.m)
    hist=[]
    v=None; y=None; sigma=None; chi=None
    for it in range(ITERS):
        sys_lr.set_theta(theta)
        ev, EV = np.linalg.eig(sys_lr.J @ sys_lr.G)
        idx = int(np.argmin(ev.imag))
        lam = ev[idx]
        y = EV[:, idx].astype(complex)
        v = sys_lr.w(y)
        sigma = float(-lam.imag)
        chi = float(np.real(1j * np.vdot(y, sys_lr.J @ y)))
        norm2 = float(np.vdot(v,v).real)
        theta_new=np.angle(v)
        mix=(1.0-BETA)*np.exp(1j*theta)+BETA*np.exp(1j*theta_new)
        theta=np.angle(mix)
        if it >= ITERS-20:
            # circular phase change between current output and mixed next theta
            d=np.angle(np.exp(1j*(theta-theta_new)))
            hist.append(dict(N=N,iteration=it+1,sigma=sigma,sigma2=sigma*sigma,chi=chi,
                             norm2=norm2,r2=norm2/sys_lr.m,Nr2=N*norm2/sys_lr.m,
                             phase_mix_rms=float(np.sqrt(np.mean(d*d))),
                             identity_error=float(abs(norm2-sigma*chi))))
    oldres, oldmu = old_residual(sys_lr,v)
    cres, cmu = corrected_residual(sys_lr,v)
    row=dict(N=N,M=sys_lr.m,seed=SEED_BASE+1000*N,iterations=ITERS,
             sigma=sigma,sigma2=sigma*sigma,chi=chi,norm2=float(np.vdot(v,v).real),
             r2=float(np.vdot(v,v).real/sys_lr.m),Nr2=float(N*np.vdot(v,v).real/sys_lr.m),
             old_residual=oldres,old_mu=oldmu,corrected_residual=cres,corrected_mu=cmu,
             identity_error=float(abs(np.vdot(v,v).real-sigma*chi)),
             y_norm=float(np.linalg.norm(y)))
    return row,hist

rows=[]; hist=[]
for N in range(NMIN,NMAX+1):
    r,h=stage1_trace(N); rows.append(r); hist.extend(h)
    print(f"N={N:2d} sigma2={r['sigma2']:.15g} norm2={r['norm2']:.15g} r2={r['r2']:.15g} oldres={r['old_residual']:.3e} corr={r['corrected_residual']:.3e}")

resdir=HERE/'results'
with open(resdir/'stage1_amplitude_N3_N16.csv','w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)
with open(resdir/'stage1_last20_iterations_N3_N16.csv','w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=hist[0].keys()); w.writeheader(); w.writerows(hist)

# Exact-candidate recognition only for N=3..6; this is identification, not proof.
exact={}
for r in rows[:4]:
    N=r['N']
    vals={}
    for key in ['sigma2','norm2','r2','Nr2']:
        vals[key]=str(nsimplify(r[key], tolerance=1e-12, full=True))
    vals['sigma']=str(nsimplify(r['sigma'], [sqrt(2),sqrt(3),sqrt(5),sqrt(6),sqrt(7),sqrt(14)], tolerance=1e-12, full=True))
    vals['chi']=str(nsimplify(r['chi'], [sqrt(2),sqrt(3),sqrt(5),sqrt(6),sqrt(7),sqrt(14)], tolerance=1e-12, full=True))
    exact[str(N)]=vals
with open(resdir/'exact_candidates_N3_N6.json','w') as f: json.dump(exact,f,indent=2,ensure_ascii=False)

summary=[]
summary.append('# make_parent 段階1 振幅選択の再現解析\n')
summary.append('## 目的\n')
summary.append('`original_engine.py` を変更せず、`_make_parent_phase_only` の1反復を外部診断として再現し、振幅を決める `v=W y` のノルム、`JG` 固有値、シンプレクティック偏極量を記録する。\n')
summary.append('## 固定条件\n')
summary.append(f'- N={NMIN}..{NMAX}\n- 反復 {ITERS}\n- beta={BETA}\n- seed = {SEED_BASE}+1000*N\n- `v` の正規化なし\n')
summary.append('## 恒等式の検算\n')
summary.append('各反復で `JG y=-i sigma y`, `||y||=1`, `v=W y`, `G=W^T W` より `||v||^2 = sigma * chi`, `chi=i y^† J y`。CSV の `identity_error` で確認。\n')
summary.append('## 収束残差の注意\n')
summary.append('現行 `_eigenmode_residual` は `mu=v^†(iKv)` を使い、`v` 非正規化後も分母 `v^†v` がない。比較のため現行残差と正しい Rayleigh quotient を使うスケール不変残差を両方保存した。元コードは変更していない。\n')
summary.append('## N=3..6 exact 候補\n')
for N,v in exact.items(): summary.append(f'- N={N}: {v}\n')
summary.append('## 科学的注意\n')
summary.append('`exact_candidates_N3_N6.json` は倍精度値からの代数数同定であり、記号的証明ではない。今後の証明対象は位相写像上の不変軌道と `y^†Gy` の閉形式。\n')
(resdir/'analysis_summary.md').write_text(''.join(summary),encoding='utf-8')
