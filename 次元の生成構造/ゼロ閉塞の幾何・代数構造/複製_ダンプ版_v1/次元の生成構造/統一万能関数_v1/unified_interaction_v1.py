#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""統一万能相互作用関数 v1 — 以後の全実験はこのモジュールのみを使う（木原ルール）

建築（重力論文・柱G8）: F[Ψ] = R(G[Ψ])·Ψ
  読出し G が生成子を作り、その生成する閉塞保存回転 R が状態を回す。
  IF文なし・種別分岐なし・調整定数ゼロ・宇宙の第0步から常時実行。

系譜（正本・read-only import）:
  LowRankSystem（位相差正弦生成子＋Cayley直交步・閉塞/ノルム厳密保存）
    = 自発的分裂予備実験_v1/run_n_scaling_lowrank_v1.py（abl 経由）
  → VertexEngine（媒介頂点・stage2）
  → VertexEngineV2（共有O線形部＝LowRankSystem を毎步内蔵実行・stage3）
  → 毛(η)拡張（stage3 census・和則 m*=2m_B−m_s 自動成立）
  ＝ 本モジュールの UnifiedEngine。
  二体正本（万能非弾性写像 v3・collision_step_exact）は再輸出のみ。

資格審査: run_qualification_v1.py（Q1-Q5・全通過が使用条件）。
"""
from __future__ import annotations
import importlib.util
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
MB = HERE.parent / "万能相互作用多体接続_v1"
UIM = HERE.parent / "万能非弾性写像_managed_v1"

_spec3 = importlib.util.spec_from_file_location(
    "unified_s3", MB / "run_stage3_sharedO_v2_and_hair_v1.py")
_s3 = importlib.util.module_from_spec(_spec3)
sys.modules[_spec3.name] = _s3
_spec3.loader.exec_module(_s3)

abl = _s3.abl              # LowRankSystem を含む正本（N体側）
gen3 = _s3.gen3            # 親・種の正本生成器
s2 = _s3.s2                # stage2（頂点・H_MAX）
VertexEngineV2 = _s3.VertexEngineV2

_spec_ex = importlib.util.spec_from_file_location(
    "unified_exact2body", UIM / "run_ignition_fate_exact_v3.py")
_ex = importlib.util.module_from_spec(_spec_ex)
sys.modules[_spec_ex.name] = _ex
_spec_ex.loader.exec_module(_ex)

collision_step_exact = _ex.collision_step_exact   # 二体正本（再輸出・read-only）
two_body_base = _ex.base
two_body_v1 = _ex.v1


class UnifiedEngine(VertexEngineV2):
    """統一エンジン＝共有O線形部（LowRankSystem 内蔵）＋媒介頂点＋毛(η)レジスタ。

    状態: C2 (M×Nn×Nη 複素)・M=N(N-1)/2 関係波・Nn=帯レジスタ・Nη=毛レジスタ。
    正本 stage3 main() 内 HairEngine の忠実コピー（資格審査Q1でビット同一を確認）。
    Nη=1 で正本 VertexEngineV2 に厳密退化（資格審査Q2でビット同一を確認）。
    """

    def __init__(self, n_, C2_0, wp, **kw):
        self.Nn, self.Neta = C2_0.shape[1], C2_0.shape[2]
        C0f = C2_0.reshape(C2_0.shape[0], -1)
        super().__init__(n_, C0f, wp, **kw)
        ks = np.arange(self.Nn)
        self.odd_k = (ks % 2 == 1)
        self.even_k = (ks % 2 == 0) & (ks != 0)

    def C2(self):
        return self.C.reshape(self.m, self.Nn, self.Neta)

    def _readout(self):
        P2 = np.abs(self.C2()) ** 2
        Pk = P2.sum(axis=2)
        Av = np.zeros((self.n, self.Nn))
        np.add.at(Av, self.ia, Pk)
        np.add.at(Av, self.ib, Pk)
        Sagg = Av[self.ia] + Av[self.ib] - 2 * Pk
        comb = Pk + Sagg
        Pf = comb[:, self.odd_k].sum(axis=1)
        Pb = comb[:, self.even_k].sum(axis=1)
        th = np.arctan2(np.sqrt(np.maximum(Pf, 0)), np.sqrt(np.maximum(Pb, 0)))
        return self.scale * np.sin(th) ** 2

    def _nonlinear(self):
        R = self._readout()
        if not np.any(R > 0):
            return
        C2 = self.C2()
        W = np.fft.ifft2(C2, axes=(1, 2)) * (self.Nn * self.Neta)
        Wf = W.reshape(self.m, -1)
        rate0 = self._vertex_rate(Wf, R)
        Lmax = float(np.max(np.abs(rate0))) / max(float(np.max(np.abs(Wf))), 1e-300)
        nsub = max(1, int(np.ceil(Lmax / s2.H_MAX)))
        h = 1.0 / nsub
        for _ in range(nsub):
            k1 = self._vertex_rate(Wf, R)
            k2 = self._vertex_rate(Wf + 0.5 * h * k1, R)
            k3 = self._vertex_rate(Wf + 0.5 * h * k2, R)
            k4 = self._vertex_rate(Wf + h * k3, R)
            Wf = Wf + (h / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
        W = Wf.reshape(self.m, self.Nn, self.Neta)
        self.C = (np.fft.fft2(W, axes=(1, 2)) / (self.Nn * self.Neta)
                  ).reshape(self.m, -1)


def build_standard_universe(n, delta, Nn=5, Neta=8):
    """標準宇宙の正本構成（run_genesis_v1 準拠・毛拡張・乱数規則正本準拠）。
    delta>0: 物質宇宙（シード＝gen3 親の第1セクション状態×δ・毛0）
    delta=0: 真空宇宙（完全無シード・内在床点火）
    戻り値: (engine, p2, q2)  p2,q2=誕生時親平面基底（f₂読出しの参照系）。
    """
    m = n * (n - 1) // 2
    _, v, _, _, _, _, _, Z0c, wp0 = abl.build_init(n, False)
    r2 = gen3.make_parent(n, seed=2)
    Csec = np.fft.fft(r2.relation_waves, axis=1) / n
    seed_state = Csec[:, 1] / np.linalg.norm(Csec[:, 1])
    C2_0 = np.zeros((m, Nn, Neta), complex)
    C2_0[:, 2, 0] = Z0c
    if delta > 0:
        C2_0[:, 1, 0] = delta * seed_state
    p2 = C2_0[:, 2, 0].real / np.linalg.norm(C2_0[:, 2, 0].real)
    q2 = C2_0[:, 2, 0].imag - (C2_0[:, 2, 0].imag @ p2) * p2
    nq = np.linalg.norm(q2)
    q2 = q2 / nq if nq > 1e-12 else np.zeros_like(p2)
    return UnifiedEngine(n, C2_0, wp0), p2, q2
