#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""統一万能相互作用関数 v2 — 分岐版（以後の実験はこれを使う）

v1 からの分岐（2026-08-08・木原指示）: **無駄な IF 文と閾値ガードを全廃**する。
計算不能なら NaN をそのまま伝播させる。条件判断は上位側の責務であり、
上位プログラムがアボートしても上位側で対応する。

**過去のバージョンは残す**（unified_interaction_v1.py）。公開済み論文の
再現性はそちらで確保され、本 v2 は今後の実験専用の分岐である。

--------------------------------------------------------------------------
v1 → v2 の差分（撤廃した箇所）と、それが力学に与える影響

| 箇所 | v1（撤廃前） | v2 | 力学への影響 |
|---|---|---|---|
| `_readout` | `np.maximum(Pf, 0)` / `np.maximum(Pb, 0)` | クランプ撤廃 | 数値誤差で負のパワーが出れば √ が NaN。**通常は同一** |
| `_nonlinear` | `if not np.any(R > 0): return`（頂点を丸ごと飛ばす） | 分岐撤廃・常に頂点を評価 | **★重要**: 頂点式は 0.5i(R·(cA·W−cB·W̄) + (cAR·W−cBR·W̄)) であり、**第2群は R を掛けていない**。したがって R≡0 の状況（奇数帯が厳密に空＝シード無しの真空宇宙）で v1 は頂点を完全に飛ばしていたが、v2 は R 非依存項を作用させる。**真空宇宙（δ=0）の軌道は v1 と一致しない** |
| `_nonlinear` | `max(float(np.max(np.abs(Wf))), 1e-300)` | 閾値撤廃・そのまま除算 | 状態が厳密零なら Lmax=NaN → 部分步数の算出で例外（上位で対応） |
| `build_standard_universe` | `q2 = q2/nq if nq>1e-12 else zeros` | 閾値撤廃・そのまま除算 | 親平面が退化していれば q2 が NaN |
| 同 | `if delta > 0:` でシード投入を分岐 | 分岐撤廃（δ=0 なら零を足すだけで同一） | 影響なし |
| 同 | `gen3.make_parent(n, seed=2)` の **seed 固定** | `seed` を引数化（既定 2） | 隠れた選択の明示化。**N=8 で親構成が失敗する問題**（ParentConstructionError）に対し、実験側でシードを宣言できる |

**据え置いた箇所（意図的）**: `nsub = max(1, int(np.ceil(Lmax / H_MAX)))` の
`max(1, …)`。これは部分步数が 1 以上でなければ積分そのものが定義されない
という**積分器の制御**であって、データに対する閾値ではない。H_MAX は正本
stage2 の定義値をそのまま用いる。
--------------------------------------------------------------------------

建築（重力論文・柱G8）: F[Ψ] = R(G[Ψ])·Ψ
  読出し G が生成子を作り、その生成する閉塞保存回転 R が状態を回す。

系譜（正本・read-only import）: v1 と同一（LowRankSystem → VertexEngine →
VertexEngineV2 → 毛(η)拡張）。二体正本 collision_step_exact は再輸出のみ。
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
    "unified_s3_v2", MB / "run_stage3_sharedO_v2_and_hair_v1.py")
_s3 = importlib.util.module_from_spec(_spec3)
sys.modules[_spec3.name] = _s3
_spec3.loader.exec_module(_s3)

abl = _s3.abl
gen3 = _s3.gen3
s2 = _s3.s2
VertexEngineV2 = _s3.VertexEngineV2

_spec_ex = importlib.util.spec_from_file_location(
    "unified_exact2body_v2", UIM / "run_ignition_fate_exact_v3.py")
_ex = importlib.util.module_from_spec(_spec_ex)
sys.modules[_spec_ex.name] = _ex
_spec_ex.loader.exec_module(_ex)

collision_step_exact = _ex.collision_step_exact
two_body_base = _ex.base
two_body_v1 = _ex.v1


class UnifiedEngineV2(VertexEngineV2):
    """統一エンジン v2 ＝ 共有O線形部＋媒介頂点＋毛(η)レジスタ。
    v1 との差は上表の通り——数値ガードと頂点スキップ分岐を撤廃した。"""

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
        """生成子側の読出し。負パワーのクランプを撤廃（√ が NaN を返す）。"""
        P2 = np.abs(self.C2()) ** 2
        Pk = P2.sum(axis=2)
        Av = np.zeros((self.n, self.Nn))
        np.add.at(Av, self.ia, Pk)
        np.add.at(Av, self.ib, Pk)
        Sagg = Av[self.ia] + Av[self.ib] - 2 * Pk
        comb = Pk + Sagg
        Pf = comb[:, self.odd_k].sum(axis=1)
        Pb = comb[:, self.even_k].sum(axis=1)
        with np.errstate(invalid="ignore"):
            th = np.arctan2(np.sqrt(Pf), np.sqrt(Pb))
        return self.scale * np.sin(th) ** 2

    def _nonlinear(self):
        """媒介頂点。**頂点スキップ分岐と 1e-300 ガードを撤廃**。"""
        R = self._readout()
        C2 = self.C2()
        W = np.fft.ifft2(C2, axes=(1, 2)) * (self.Nn * self.Neta)
        Wf = W.reshape(self.m, -1)
        rate0 = self._vertex_rate(Wf, R)
        with np.errstate(divide="ignore", invalid="ignore"):
            Lmax = float(np.max(np.abs(rate0))) / float(np.max(np.abs(Wf)))
        nsub = max(1, int(np.ceil(Lmax / s2.H_MAX)))   # 積分器の制御（閾値でない）
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


def build_standard_universe(n, delta, Nn=5, Neta=8, seed=2):
    """標準宇宙の構成（正本 run_genesis_v1 準拠・毛拡張）。

    真空 = control親（k=2・巻き0）／シード = 親第1セクション×δ（k=1・巻き0）。
    delta>0: 物質宇宙、delta=0: 真空宇宙（シード項が厳密に零になるだけ）。
    **seed は引数**（v1 は 2 固定。N=8 では seed=2 で親構成が失敗するため、
    実験側がシードを宣言できるようにした）。
    戻り値: (engine, p2, q2) — p2,q2 は誕生時の親平面基底。
    """
    m = n * (n - 1) // 2
    _, v, _, _, _, _, _, Z0c, wp0 = abl.build_init(n, False)
    r2 = gen3.make_parent(n, seed=seed)
    Csec = np.fft.fft(r2.relation_waves, axis=1) / n
    seed_state = Csec[:, 1] / np.linalg.norm(Csec[:, 1])
    C2_0 = np.zeros((m, Nn, Neta), complex)
    C2_0[:, 2, 0] = Z0c
    C2_0[:, 1, 0] = delta * seed_state          # δ=0 なら零（分岐不要）
    p2 = C2_0[:, 2, 0].real / np.linalg.norm(C2_0[:, 2, 0].real)
    q2 = C2_0[:, 2, 0].imag - (C2_0[:, 2, 0].imag @ p2) * p2
    with np.errstate(divide="ignore", invalid="ignore"):
        q2 = q2 / np.linalg.norm(q2)            # 退化なら NaN（代用しない）
    return UnifiedEngineV2(n, C2_0, wp0), p2, q2
