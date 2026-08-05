#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 本ファイルはセッション記録から回収した原本コード（2026-08-05実行・結果JSONと対応）。
# 回収時にbashラッパーとgrepフィルタのみ除去、コード本体は無変更。
"""GENESIS 3D デモ: 3次元レジスタの粒子——位置・速度・保存則・生成局所性"""
import importlib.util, sys, json
from pathlib import Path
import numpy as np
HERE = Path(".").resolve()
spec3 = importlib.util.spec_from_file_location("s3d", HERE/"run_stage3_sharedO_v2_and_hair_v1.py")
s3 = importlib.util.module_from_spec(spec3); sys.modules[spec3.name]=s3; spec3.loader.exec_module(s3)
abl, V2, H_MAX = s3.abl, s3.VertexEngineV2, s3.s2.H_MAX

N_GRAPH = 12
D = (8, 8, 8)                      # 3次元レジスタ
NF = D[0]*D[1]*D[2]
DELTA = 3e-2
VEL = (0.05, 0.03, 0.02)           # 予言速度（セル/step、各軸）
m = N_GRAPH*(N_GRAPH-1)//2

class Engine3D(V2):
    def __init__(self, n, C3_0, wp, vel, **kw):
        super().__init__(n, C3_0.reshape(C3_0.shape[0], -1), wp, **kw)
        k1 = np.arange(D[0]); k2 = np.arange(D[1]); k3 = np.arange(D[2])
        # 分散（3軸調和閉鎖・並進演算子）: c_k ← c_k e^{-2πi k·s/N}
        ph = np.exp(-2j*np.pi*(k1[:,None,None]*vel[0]/D[0]
                              + k2[None,:,None]*vel[1]/D[1]
                              + k3[None,None,:]*vel[2]/D[2]))
        self.disp = ph.reshape(1, -1)
        self.odd1 = (np.arange(D[0]) % 2 == 1)
    def C3(self):
        return self.C.reshape(self.m, *D)
    def _readout(self):
        P = np.abs(self.C3())**2
        Pk1 = P.sum(axis=(2,3))
        Av = np.zeros((self.n, D[0]))
        np.add.at(Av, self.ia, Pk1); np.add.at(Av, self.ib, Pk1)
        Sagg = Av[self.ia]+Av[self.ib]-2*Pk1
        comb = Pk1+Sagg
        Pf = comb[:, self.odd1].sum(axis=1)
        ev = (np.arange(D[0])%2==0)&(np.arange(D[0])!=0)&(np.arange(D[0])!=D[0]//2)
        Pb = comb[:, ev].sum(axis=1)
        th = np.arctan2(np.sqrt(np.maximum(Pf,0)), np.sqrt(np.maximum(Pb,0)))
        return self.scale*np.sin(th)**2
    def _nonlinear(self):
        R = self._readout()
        if not np.any(R>0): return
        C3 = self.C3()
        W = np.fft.ifftn(C3, axes=(1,2,3))*NF
        Wf = W.reshape(self.m, -1)
        r0 = self._vertex_rate(Wf, R)
        L = float(np.max(np.abs(r0)))/max(float(np.max(np.abs(Wf))),1e-300)
        ns = max(1, int(np.ceil(L/H_MAX))); h = 1.0/ns
        for _ in range(ns):
            a = self._vertex_rate(Wf, R); b = self._vertex_rate(Wf+0.5*h*a, R)
            c = self._vertex_rate(Wf+0.5*h*b, R); d_ = self._vertex_rate(Wf+h*c, R)
            Wf = Wf + (h/6)*(a+2*b+2*c+d_)
        self.C = (np.fft.fftn(Wf.reshape(self.m,*D), axes=(1,2,3))/NF).reshape(self.m,-1)
    def step(self):
        super().step()
        self.C = self.C*self.disp

def odd_P3(C3):
    mask = (np.arange(D[0])%2==1)
    Codd = C3*mask[None,:,None,None]
    W = np.fft.ifftn(Codd, axes=(1,2,3))*NF
    return np.sum(np.abs(W)**2, axis=0)     # P(n1,n2,n3)

def centroid3(P):
    out=[]
    for ax, Nax in enumerate(D):
        marg = P.sum(axis=tuple(a for a in range(3) if a!=ax))
        z = np.sum(marg*np.exp(2j*np.pi*np.arange(Nax)/Nax))/marg.sum()
        out.append(float((np.angle(z)*Nax/(2*np.pi))%Nax))
    pr = float(P.sum()**2/np.sum(P**2))
    return out, pr

# 初期状態: ポンプ=(k1=2,0,0)にcontrol親、種=3D局在パケット
_, v, _, _, _, _, _, Z0c, wp0 = abl.build_init(N_GRAPH, False)
Z0c = Z0c/np.linalg.norm(Z0c)
spec_g3 = importlib.util.spec_from_file_location("g3z", HERE/"run_genesis_v3_register_local_v1.py")
g3 = importlib.util.module_from_spec(spec_g3); sys.modules[spec_g3.name]=g3; spec_g3.loader.exec_module(g3)
seed_edge = g3.zero_closure_state(m, np.random.default_rng(98000))

C3_0 = np.zeros((m,)+D, complex)
C3_0[:,2,0,0] = Z0c
odd1 = [k for k in range(D[0]) if k%2==1]
prof = np.zeros(D, complex)
for k1 in odd1:
    for k2 in range(D[1]):
        for k3 in range(D[2]):
            if (2*k2)%D[1]==0 and k2!=0: continue
            if (2*k3)%D[2]==0 and k3!=0: continue
            prof[k1,k2,k3] = 1.0
prof /= np.linalg.norm(prof)
C3_0 += DELTA*prof[None,:,:,:]*seed_edge[:,None,None,None]

eng = Engine3D(N_GRAPH, C3_0, wp0, VEL, vertex_on=True)
P0 = odd_P3(eng.C3()); x0, pr0 = centroid3(P0)
W0 = np.fft.ifftn(eng.C3(), axes=(1,2,3))*NF
cl0 = np.abs(np.einsum('eijk,eijk->ijk', W0, W0)); n0 = float(np.sum(np.abs(W0)**2))
print(f"3D粒子 初期位置 n⃗=({x0[0]:.2f},{x0[1]:.2f},{x0[2]:.2f}) PR_3D={pr0:.1f} 予言速度={VEL}")
T = 300
for t in range(T):
    eng.step()
    if (t+1)%100==0:
        P = odd_P3(eng.C3()); x, pr = centroid3(P)
        pred = [ (x0[a]+VEL[a]*(t+1))%D[a] for a in range(3) ]
        print(f"  t={t+1}: 実測 n⃗=({x[0]:.2f},{x[1]:.2f},{x[2]:.2f}) "
              f"予言=({pred[0]:.2f},{pred[1]:.2f},{pred[2]:.2f}) PR={pr:.1f}")
W1 = np.fft.ifftn(eng.C3(), axes=(1,2,3))*NF
cl1 = np.abs(np.einsum('eijk,eijk->ijk', W1, W1))
n1 = float(np.sum(np.abs(W1)**2))
# 分散=並進は閉塞値も並進させるので、集合として比較（ソート最大差）
cld = float(np.max(np.abs(np.sort(cl1.flatten())-np.sort(cl0.flatten()))))
print(f"保存則: 点ごと閉塞（集合比較）={cld:.2e} ノルム相対={(abs(n1-n0)/n0):.2e}")
json.dump({"x0": x0, "PR0": pr0, "VEL": VEL, "T": T,
           "closure_setdrift": cld, "norm_drift": abs(n1-n0)/n0},
          open("genesis_3d_demo_result_v1.json","w"), indent=1)
print("saved")
