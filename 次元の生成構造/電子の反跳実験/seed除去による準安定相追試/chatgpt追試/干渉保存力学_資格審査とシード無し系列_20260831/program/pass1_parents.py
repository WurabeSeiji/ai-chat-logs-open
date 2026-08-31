#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""パス1：親の生成・受け入れ検査・保存と、走行前予測の固定（走行より前に必ず完了）。
親（タグ）：
  hm_N{3..16} = 手作り等モジュラー（偶数 1-因子分解／奇数 距離クラス／N=3 Z3。局所閉塞＋等モジュラー、N=3 のみ S_i≠0）
  ne_N{3..16} = 非等モジュラー（q≥4 クラス重み付き族 a_c=r̄²(1+0.6cos(4πc/q))、q≤3 は旧フレーム多様体上の代表点）
  rb_N{5..16} = 乱数均衡親 {S_i=0, W_i=W0}（対称性なし、seed=100+N）
スケール規約：全親 ‖v‖² = M·r̄²、r̄²=1/15（スケールは時計にしか効かない）。
受け入れ：大域閉塞 |Σz²|/H < 1e-12。旧フレームの自己無撞着量は記録のみ（本パッケージの力学は新フレーム）。
新フレーム判定：res_new/r̄² < 1e-10 → equilibrium（アンカー）、それ以外 → non-equilibrium。
予測（走行前に固定）：
  equilibrium → 共回転 1 刻み線形化 G=R(β)DΦ の ρ。ρ−1>1e-3 → inflating（λ_f=2lnρ、t50=(ln0.5−ln3e-32)/λ_f）、未満 → neutral（走行中 保持）
  non-equilibrium → 即時ドリフト。1 step 相対変位の予測 disp1 = Δ·res_new。閉塞 |Σz²|/H が開く（増加する）
  閉塞：equilibrium は全走行で < 1e-10 に留まる"""
import os, sys, csv, json, math
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from common import edges, adjacency, selfconsistency
from interference_dynamics import hermitian_H, unified_interference_step, unified_readout, DELTA
import state_provider as sp
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RBAR2 = 1.0/15.0
F0 = 3e-32

def SW(N, v):
    E = edges(N); S = np.zeros(N, complex); W = np.zeros(N)
    for k, (i, j) in enumerate(E):
        S[i] += v[k]**2; S[j] += v[k]**2; W[i] += abs(v[k])**2; W[j] += abs(v[k])**2
    return S, W

def solve_balanced(N, rng, W0, iters=100):
    """v2 補完実験 pass1_parents.py と同一（コピー）。{S_i=0, W_i=W0} への Newton。"""
    E = edges(N); M = len(E)
    v = rng.standard_normal(M) + 1j*rng.standard_normal(M)
    v *= np.sqrt(N*W0/2)/np.linalg.norm(v)
    def G(x):
        vv = x[:M] + 1j*x[M:]; S, W = SW(N, vv)
        return np.r_[S.real, S.imag, W - W0]
    for it in range(iters):
        x = np.r_[v.real, v.imag]; g = G(x)
        if np.linalg.norm(g) < 1e-15: break
        J = np.zeros((3*N, 2*M)); h = 1e-7
        for j in range(2*M):
            e = np.zeros(2*M); e[j] = h
            J[:, j] = (G(x+e) - G(x-e))/(2*h)
        x = x - np.linalg.lstsq(J, g, rcond=1e-12)[0]
        v = x[:M] + 1j*x[M:]
    return v

def monodromy_new(N, v, A, h=1e-7):
    """新フレームの共回転 1 刻み線形化：G=R(β)DΦ、β=−arg⟨v,Φ(v)⟩。"""
    M = len(v)
    x = np.r_[v.real, v.imag]
    D = np.zeros((2*M, 2*M))
    for j in range(2*M):
        e = np.zeros(2*M); e[j] = h
        zp = unified_interference_step((x+e)[:M] + 1j*(x+e)[M:], A)
        zm = unified_interference_step((x-e)[:M] + 1j*(x-e)[M:], A)
        D[:, j] = np.r_[(zp - zm).real, (zp - zm).imag]/(2*h)
    beta = -float(np.angle(np.vdot(v, unified_interference_step(v, A))))
    c, s = math.cos(beta), math.sin(beta)
    R = np.block([[c*np.eye(M), -s*np.eye(M)], [s*np.eye(M), c*np.eye(M)]])
    ev = np.linalg.eigvals(R @ D)
    return float(np.abs(ev).max()), int((np.abs(ev) > 1 + 1e-9).sum())

def save(tag, N, v, design, extra):
    A = adjacency(N); E = edges(N)
    sc_old = selfconsistency(N, v, A)          # 旧フレーム量（記録のみ）
    ro = unified_readout(v, A, E)              # 新フレーム量
    r2 = ro["H_total"]/len(E)
    res_rel = ro["residual_new"]/r2
    is_eq = res_rel < 1e-10
    ok = dict(global_closure=sc_old["global_closure"] < 1e-12, norm=abs(ro["H_total"] - len(E)*RBAR2) < 1e-9)
    if is_eq:
        rho, nun = monodromy_new(N, v, A)
        lam = 2*math.log(rho) if rho > 0 else float("nan")
        kind = "eq_inflating" if rho - 1 > 1e-3 else "eq_neutral"
        pred = dict(pred_kind=kind, pred_rho_minus_1=rho-1, pred_lambda_f=lam,
                    pred_n_unstable=nun,
                    pred_t50=(float((math.log(0.5)-math.log(F0))/lam) if lam > 1e-6 else None),
                    pred_disp1=None)
    else:
        pred = dict(pred_kind="non_equilibrium", pred_rho_minus_1=None, pred_lambda_f=None,
                    pred_n_unstable=None, pred_t50=None,
                    pred_disp1=float(DELTA*ro["residual_new"]))
    dd = os.path.join(ROOT, "data", tag); os.makedirs(dd, exist_ok=True)
    np.savez_compressed(os.path.join(dd, "parent_v.npz"), v=v, edges=np.array(E),
                        design=design, mu_new=ro["mu_new"], residual_new=ro["residual_new"],
                        mu_old=sc_old["mu"], residual_old=sc_old["residual"], rbar2=r2)
    rep = dict(tag=tag, N=N, method=tag.split("_")[0], design=design, M=len(E),
               norm=float(np.linalg.norm(v)), rbar2=r2,
               amp_spread_rel=float(np.ptp(np.abs(v)**2)/(np.abs(v)**2).mean()),
               global_closure=sc_old["global_closure"], local_closure=sc_old["local_closure"],
               local_closed=bool(sc_old["local_closure"] < 1e-10),
               mu_old=sc_old["mu"], residual_old=sc_old["residual"],
               mu_new=ro["mu_new"], mu_new_over_r2=ro["mu_new"]/r2,
               residual_new_over_r2=res_rel, is_equilibrium=bool(is_eq),
               clock_pred_dphi=float(-DELTA*ro["mu_new"]),
               ok={k: bool(x) for k, x in ok.items()}, **pred, **extra)
    with open(os.path.join(dd, "parent_checks.json"), "w") as f:
        json.dump(rep, f, indent=2, ensure_ascii=False)
    print(f"{tag}: 閉塞={sc_old['global_closure']:.1e} 局所={sc_old['local_closure']:.1e} "
          f"μ_new/r²={ro['mu_new']/r2:+.6f} res_new/r²={res_rel:.1e} → {rep['pred_kind']}"
          + (f" (ρ−1={pred['pred_rho_minus_1']:.2e})" if is_eq else f" (disp1={pred['pred_disp1']:.2e})"))
    if not all(ok.values()):
        raise SystemExit(f"ABORT {tag}: {ok}")
    return rep

rows = []
for N in range(3, 17):
    M = N*(N-1)//2
    scale = math.sqrt(M*RBAR2)
    v = sp.equimodular(N); v *= scale/np.linalg.norm(v)
    rows.append(save(f"hm_N{N}", N, v,
                     "handmade_equimodular_" + ("Z3" if N == 3 else "1factor" if N % 2 == 0 else "distance_classes"), dict()))
    v, kind, col, q, step = sp.state(N); v *= scale/np.linalg.norm(v)
    rows.append(save(f"ne_N{N}", N, v,
                     "nonequimodular_" + ("class_family_k2_eps0.6" if kind == "class" else "manifold_point_oldframe"), dict(q=q)))
    if N >= 5:
        rng = np.random.default_rng(100 + N)
        W0 = 2*(M*RBAR2)/N
        v = solve_balanced(N, rng, W0)
        v *= scale/np.linalg.norm(v)
        rows.append(save(f"rb_N{N}", N, v, "random_balanced_S0_Wconst", dict(seed_rng=100 + N)))

os.makedirs(os.path.join(ROOT, "results"), exist_ok=True)
keys = [k for k in rows[0].keys() if k not in ("ok",)]
allkeys = sorted({k for r in rows for k in r.keys() if k != "ok"}, key=lambda k: (keys.index(k) if k in keys else 99, k))
with open(os.path.join(ROOT, "results", "parents_predictions.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=allkeys, extrasaction="ignore")
    w.writeheader(); w.writerows(rows)
print(f"PASS1 OK（親 {len(rows)} 個、予測は走行前に固定）")
