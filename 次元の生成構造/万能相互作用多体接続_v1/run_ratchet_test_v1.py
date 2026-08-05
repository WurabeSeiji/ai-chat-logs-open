#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 本ファイルはセッション記録から回収した原本コード（2026-08-05実行・結果JSONと対応）。
# 回収時にbashラッパーとgrepフィルタのみ除去、コード本体は無変更。
"""ラチェット検証: 50/50くじ＋自己触媒→アンサンブル正味成長か"""
import importlib.util, sys, json
from pathlib import Path
import numpy as np
HERE = Path(".").resolve()
spec3 = importlib.util.spec_from_file_location("s3r", HERE/"run_genesis_v3_register_local_v1.py")
g3 = importlib.util.module_from_spec(spec3); sys.modules[spec3.name]=g3; spec3.loader.exec_module(g3)
abl, V2 = g3.abl, g3.V2
N_GRAPH, NREG, DELTA = 12, 16, 3e-2
m = N_GRAPH*(N_GRAPH-1)//2
_, v, _, _, _, _, _, Z0c, wp0 = abl.build_init(N_GRAPH, False)
Z0c = Z0c/np.linalg.norm(Z0c)
seed_edge = g3.zero_closure_state(m, np.random.default_rng(98000))
odd_ks = [k for k in range(NREG) if k%2==1]
T = 800
trajs = []
for ps in range(20):
    rng = np.random.default_rng(100000+ps)
    prof = np.zeros(NREG, complex)
    for k in odd_ks:
        prof[k] = np.exp(1j*rng.uniform(0,2*np.pi))/np.sqrt(len(odd_ks))
    C0 = np.zeros((m,NREG), complex); C0[:,2]=Z0c
    for k in range(NREG):
        if abs(prof[k])>0: C0[:,k]+=DELTA*prof[k]*seed_edge
    eng = V2(N_GRAPH, C0, wp0, vertex_on=True)
    fs = []
    for t in range(T):
        eng.step(); fs.append(eng.diagnostics()["f_seed"])
    trajs.append(fs)
trajs = np.array(trajs)
f0 = trajs[0,0]
mean_t = trajs.mean(axis=0)
frac_above = np.mean(trajs[:,-1] > trajs[:,0])
print(f"アンサンブル20本 T={T}: f0={f0:.3e}")
print(f"  平均 f_seed: t=0 {mean_t[0]:.4e} → t={T} {mean_t[-1]:.4e}（比 {mean_t[-1]/mean_t[0]:.4f}）")
print(f"  終点で初期値超え: {frac_above*100:.0f}%  最大伸び={trajs[:,-1].max()/f0:.3f}倍  最小={trajs[:,-1].min()/f0:.3f}倍")
json.dump({"T": T, "mean_ratio": float(mean_t[-1]/mean_t[0]),
           "frac_final_above_init": float(frac_above),
           "max_gain": float(trajs[:,-1].max()/f0), "min_gain": float(trajs[:,-1].min()/f0)},
          open("ratchet_test_result_v1.json","w"), indent=1)
print("saved")
