#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 本ファイルはセッション記録から回収した原本コード（2026-08-05実行・結果JSONと対応）。
# 回収時にbashラッパーとgrepフィルタのみ除去、コード本体は無変更。
"""運動学の実証: 調和閉鎖分散 ω_k = k·ω₁ をレジスタ時計として導入→パケット並進"""
import importlib.util, sys, json
from pathlib import Path
import numpy as np
HERE = Path(".").resolve()
spec3 = importlib.util.spec_from_file_location("s3kin", HERE/"run_genesis_v3_register_local_v1.py")
g3 = importlib.util.module_from_spec(spec3); sys.modules[spec3.name]=g3; spec3.loader.exec_module(g3)
abl, V2 = g3.abl, g3.V2
N_GRAPH, NREG, DELTA = 12, 16, 3e-2
m = N_GRAPH*(N_GRAPH-1)//2
_, v, _, _, _, _, _, Z0c, wp0 = abl.build_init(N_GRAPH, False)
Z0c = Z0c/np.linalg.norm(Z0c)
seed_edge = g3.zero_closure_state(m, np.random.default_rng(98000))
odd_ks = [k for k in range(NREG) if k%2==1]

class KineticEngine(V2):
    """v2＋調和閉鎖レジスタ時計: 各倍音kにe^{i k ω1}を毎step適用（分散=一様並進）"""
    def __init__(self, n, C0, wp, omega1, **kw):
        super().__init__(n, C0, wp, **kw)
        self.disp = np.exp(1j*np.arange(self.nreg)*omega1)[None,:]
    def step(self):
        super().step()
        self.C = self.C * self.disp

def centroid(C):
    ks = np.arange(C.shape[1]); Codd = C*((ks%2==1)[None,:])
    W = np.fft.ifft(Codd, axis=1)*C.shape[1]
    P = np.sum(np.abs(W)**2, axis=0)
    ph = np.exp(2j*np.pi*np.arange(C.shape[1])/C.shape[1])
    z = np.sum(P*ph)/P.sum()
    return float((np.angle(z)*C.shape[1]/(2*np.pi)) % C.shape[1]), float(P.sum()**2/np.sum(P**2))

prof = np.zeros(NREG, complex)
for k in odd_ks: prof[k] = 1.0/np.sqrt(len(odd_ks))
C0 = np.zeros((m,NREG), complex); C0[:,2]=Z0c
for k in range(NREG):
    if abs(prof[k])>0: C0[:,k]+=DELTA*prof[k]*seed_edge

omega1 = 2*np.pi/NREG * 0.05   # 群速度 0.05セル/step
eng = KineticEngine(N_GRAPH, C0, wp0, omega1, vertex_on=True)
x0, pr0 = centroid(eng.C)
print(f"初期: 重心n={x0:.2f} PR={pr0:.1f}  予言速度=0.05セル/step（調和閉鎖 ω_k=k·ω₁）")
rows=[]
for t in range(400):
    eng.step()
    if (t+1)%100==0:
        x, pr = centroid(eng.C)
        rows.append((t+1, round(x,2), round(pr,2)))
        print(f"  t={t+1}: 重心n={x:.2f}（予言 {(x0+0.05*(t+1))%NREG:.2f}） PR={pr:.2f}")
# 保存則確認
W = np.fft.ifft(eng.C, axis=1)*NREG
print(f"ノルム={float(np.sum(np.abs(W)**2)):.12f}（初期1.000...+δ²）")
json.dump({"omega1_cells_per_step":0.05, "trajectory": rows},
          open("kinetic_dispersion_demo_v1.json","w"), indent=1)
print("saved")
