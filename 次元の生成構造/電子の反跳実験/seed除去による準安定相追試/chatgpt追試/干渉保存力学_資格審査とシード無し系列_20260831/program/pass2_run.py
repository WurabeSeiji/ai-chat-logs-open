#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""パス2：シード無し走行（データ収集のみ・図なし）。力学は unified_interference_step のみ（毎步一様 1 回・直書き禁止）。
Z0 = 親そのまま（シード無し。摂動源は親残差と丸めのみ）。STEPS=40000、Δ=2π/124。
記録：ユニタリ性（H_total ドリフト）・閉塞率 |ΣZ²|/H（反射チャネル）・親平面外 H⊥・重なり欠損・PR/M・振幅統計・位相進み。
状態はキー step のスナップショットと最終状態のみ保存（全 step 状態は保存しない）。"""
import os, sys, csv, json, math, re
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from common import edges, adjacency
from interference_dynamics import unified_interference_step, unified_readout, DELTA, L
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TAG = sys.argv[1]
N = int(re.search(r"N(\d+)", TAG).group(1))
STEPS = 40000
DATA = os.path.join(ROOT, "data", TAG)
KEY_STEPS = sorted(set([0, 25, 50, 75, 100, 125, 150, 200, 300, 500, 750, 1000, 1500, 2000,
                        3000, 4000, 5000, 7500, 10000, 15000, 20000, 25000, 30000, 35000, 40000]))

pz = np.load(os.path.join(DATA, "parent_v.npz"), allow_pickle=True)
v = pz["v"]
A = adjacency(N); E = edges(N); M = len(E)
p = v.real/np.linalg.norm(v.real)
q = v.imag - (v.imag @ p)*p; q /= np.linalg.norm(q)
nv2 = float(np.vdot(v, v).real)

Z = v.copy()  # SEEDLESS
rows = []; snaps = {}
h0 = float(np.vdot(Z, Z).real)
for t in range(STEPS + 1):
    d2 = Z*Z
    htot = float(np.vdot(Z, Z).real)
    closure = float(abs(d2.sum())/htot)
    Zperp = Z - p*(p @ Z) - q*(q @ Z)
    hperp = float(np.vdot(Zperp, Zperp).real)
    a2 = np.abs(Z)**2
    pr = float((a2.sum()**2)/(a2**2).sum())
    ov = float(abs(np.vdot(v, Z))**2/(nv2*htot))
    if t in KEY_STEPS:
        snaps[t] = Z.copy()
    if t == STEPS:
        rows.append([t, htot, closure, hperp, hperp/htot, pr/M,
                     float(np.abs(Z).min()), float(np.abs(Z).max()), float(np.abs(Z).std()),
                     1.0 - ov, float("nan")])
        break
    Zn = unified_interference_step(Z, A)
    dphi = float(np.angle(np.vdot(Z, Zn)))
    rows.append([t, htot, closure, hperp, hperp/htot, pr/M,
                 float(np.abs(Z).min()), float(np.abs(Z).max()), float(np.abs(Z).std()),
                 1.0 - ov, dphi])
    Z = Zn
a = np.asarray(rows, float)

headers = ["step", "H_total", "closure_frac", "H_perp", "H_perp_frac", "PR_over_M",
           "amp_min", "amp_max", "amp_std", "overlap_deficit", "dphi"]
with open(os.path.join(DATA, "timeseries.csv"), "w", newline="") as f:
    w = csv.writer(f); w.writerow(headers); w.writerows(rows)
ks = sorted(snaps.keys())
np.savez_compressed(os.path.join(DATA, "snapshots.npz"),
                    steps=np.array(ks), Z=np.array([snaps[t] for t in ks]))
np.savez_compressed(os.path.join(DATA, "final_state.npz"), Z=Z)

def growth_fit(a, lo=1e-10, hi=1e-3):
    y = a[:, 3]; mask = (y > lo) & (y < hi) & np.isfinite(y)
    if mask.sum() < 3: return None
    x = a[mask, 0]; ly = np.log(y[mask])
    slope, inter = np.polyfit(x, ly, 1)
    pred = slope*x + inter
    ssr = ((ly - pred)**2).sum(); sst = ((ly - ly.mean())**2).sum()
    return dict(slope_ln_Hperp_per_step=float(slope), intercept=float(inter),
                R2=float(1 - ssr/sst) if sst > 0 else None, n=int(mask.sum()),
                step_min=int(x.min()), step_max=int(x.max()))

ro_fin = unified_readout(Z, A, E)
r2fin = ro_fin["H_total"]/M
frac = a[:, 4]
ix = np.where(frac > 0.05)[0]
summary = dict(
    experiment=f"{TAG} interference-preserving frame, seedless, exp(-i(2pi/{L})H) STEPS={STEPS}",
    N=N, M=M, L=L, steps=STEPS,
    parent_design=str(pz["design"]), parent_mu_new=float(pz["mu_new"]),
    parent_residual_new=float(pz["residual_new"]),
    unitarity_max_rel_drift=float(np.nanmax(np.abs(a[:, 1] - h0))/h0),
    closure_max=float(np.nanmax(a[:, 2])), closure_final=float(a[-1, 2]),
    Hperp_frac_max=float(np.nanmax(frac)),
    onset_Hperp_frac_gt_0p05=(int(a[ix[0], 0]) if len(ix) else None),
    overlap_deficit_max=float(np.nanmax(a[:, 9])), overlap_deficit_final=float(a[-1, 9]),
    disp1_measured=float(math.sqrt(max(a[1, 9], 0.0)*2.0)) if STEPS >= 1 else None,
    PR_over_M_final=float(a[-1, 5]),
    amp_spread_final=float((a[-1, 7] - a[-1, 6])/max(a[-1, 7], 1e-300)),
    dphi_mean_first200=float(np.nanmean(a[:200, 10])),
    dphi_mean_last200=float(np.nanmean(a[-201:-1, 10])),
    clock_pred_dphi=float(-DELTA*float(pz["mu_new"])),
    final_mu_new=ro_fin["mu_new"], final_mu_new_over_r2=ro_fin["mu_new"]/r2fin,
    final_residual_new_over_r2=ro_fin["residual_new"]/r2fin,
    final_global_closure=ro_fin["global_closure"], final_local_closure=ro_fin["local_closure"],
    growth_fit=growth_fit(a))
with open(os.path.join(DATA, "summary.json"), "w") as f:
    json.dump(summary, f, indent=2, ensure_ascii=False)
print(f"RUN {TAG}: unitarity={summary['unitarity_max_rel_drift']:.1e} closure_max={summary['closure_max']:.1e} "
      f"Hperp_frac_max={summary['Hperp_frac_max']:.2e} ov_def_max={summary['overlap_deficit_max']:.2e} "
      f"PR/M_fin={summary['PR_over_M_final']:.3f}")
