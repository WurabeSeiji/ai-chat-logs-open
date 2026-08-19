#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""統一場の万能読出し関数 v5 — 二体対読出し（流れの二文法）を追加

v4 からの分岐（2026-08-19）: **三部作第一部（二文法分解・DOI 21763995）の
流れ計器を G の正式メンバーに昇格**する。v4 までの全メンバー（座標読出しを含む）は
無改訂で動的に継承する。

--------------------------------------------------------------------------
なぜ v5 が要るか

[F2] の統一主張（重力＋クーロン＝同一相互作用の二復調）の電荷側の力観測量は、
三部作第一部が確定した「流れの二項恒等式」

    dN_B = sin²θ·(N_A − N_B)  ＋  sin2θ·Re⟨a|b⟩
           （大きさの項＝重力文法）   （重なりの項＝電荷文法）

である。符号＝共有チャネル上の相対位相（φ 半回転で流れが反転・反対称 2.2e-16）、
結合の大きさ＝振幅の積、遮蔽＝成分ゼロ／チャネル非共有、量子化＝有理ロック。
v4 までの G にはこの対読出し（チャネル別ノルム・ビン別交差スペクトル・
三項分解・相対位相）が存在しなかった（v1 の g_spin_stokes_candidate は
大域内積のみ・テストベッド適用禁止の候補扱い）。v5 がその口を作る。

--------------------------------------------------------------------------
設計（ゲージ構造の明記）

- **個体の位相はゲージ、対の相対位相が物理**。共通位相変換 (a,b)→(e^{iγ}a, e^{iγ}b)
  の下で全対読出しは不変（[F1] §2.1 の対読出し構造と同一。資格審査 A5）。
  片側の位相シフト（例 b→e^{iπ}b）は**ゲージ変換ではなく状態準備（電荷共役）**であり、
  物理（流れの向き）を変える。±電荷の表現は状態側（部分共有構成＋片側位相ラベル）が
  担い、読出しは相対位相をゲージ不変量として読む。読出しのゲージ変換で
  「同符号＝反力」を作ることはできないし、してはならない（無名性）。
- **θ は再実装しない**。θ（結合角）の正本は two_body_v1.toy.theta_from_ab であり、
  g_pair_flow は θ を引数に取る（読出し層は力学定数を持たない）。
- 「同符号＝斥力」への割当は読出しの外＝辞書（翻訳層）の仕事であり、
  三部作の予言どおり距離辞書の完成と同時に自動判定される（DL7-F2）。

資格審査: run_qualification_pair_v5.py（A1〜A5・三部作アンカー）
仕様: 統一万能関数_仕様_v2.md 末尾の v5 追記
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np

# ---- v4 全メンバーの無改訂継承（動的 import・同ディレクトリ） ----
_HERE = Path(__file__).resolve().parent
_spec4 = importlib.util.spec_from_file_location("unified_readout_v4_for_v5",
                                                _HERE / "unified_readout_v4.py")
_v4 = importlib.util.module_from_spec(_spec4)
sys.modules["unified_readout_v4_for_v5"] = _v4
_spec4.loader.exec_module(_v4)
for _name in dir(_v4):
    if _name.startswith("g_") or _name.startswith("_odd"):
        globals()[_name] = getattr(_v4, _name)


# ======================================================================
# 二体対読出し（v5 新設。三部作第一部 §2〜§6 の計器の G 化）
# ======================================================================

def g_pair_norms(a, b):
    """チャネル別ノルム（大きさの項の材料）。"""
    return {"N_A": float(np.vdot(a, a).real), "N_B": float(np.vdot(b, b).real)}


def g_pair_overlap(a, b):
    """大域内積 ⟨a|b⟩（重なりの項の材料）。共通位相に対して相対位相のみ可観測。"""
    z = complex(np.vdot(a, b))
    return {"overlap_re": z.real, "overlap_im": z.imag,
            "overlap_abs": abs(z), "overlap_phase": float(np.angle(z))}


def g_pair_cross_spectrum(a, b, nc, ne, kmax=16):
    """ビン別交差スペクトル ⟨A_j|B_k⟩ と共有チャネルの梯子分解。

    A_j, B_k は χ 方向 FFT のビン成分（η は内積に畳む）。返り値:
      ladder: {(j,k): (c_abs, phi)} で c_abs>1e-12 の対のみ
      三部作の実測: 搬送波あり→全ゼロ（ヌル定理）、なし→Δk=±2・c=0.5 厳密
    """
    Fa = np.fft.fft(a.reshape(nc, ne), axis=0)
    Fb = np.fft.fft(b.reshape(nc, ne), axis=0)
    ladder = {}
    for j in range(kmax + 1):
        for k in range(kmax + 1):
            z = complex(np.vdot(Fa[j], Fb[k])) / (nc)
            if abs(z) > 1e-12:
                ladder[f"{j},{k}"] = (abs(z), float(np.angle(z)))
    return {"ladder": ladder, "n_shared": len(ladder)}


def g_pair_flow(a, b, theta):
    """流れの二項恒等式（三部作第一部 §3）。θ は正本（toy.theta_from_ab 等）から渡す。

      dN_B_pred = sin²θ·(N_A−N_B) + sin2θ·Re⟨a|b⟩

    返り値: 大きさの項（重力文法）・重なりの項（電荷文法・符号=相対位相）・予言合計。
    """
    n = g_pair_norms(a, b)
    ov = g_pair_overlap(a, b)
    mag = float(np.sin(theta) ** 2 * (n["N_A"] - n["N_B"]))
    ovl = float(np.sin(2 * theta) * ov["overlap_re"])
    return {"flow_mag_term": mag, "flow_overlap_term": ovl,
            "dN_B_pred": mag + ovl,
            "N_A": n["N_A"], "N_B": n["N_B"],
            "overlap_re": ov["overlap_re"], "overlap_phase": ov["overlap_phase"]}


def g_pair_charge_phase(a, b, nc, ne, kmax=16):
    """電荷符号の読出し: 共有チャネルごとの相対位相 φ_ch（ゲージ不変・対量）。

    共通位相変換で不変、片側 π シフト（電荷共役）で全 φ_ch が π 反転する。
    「どちらの φ が＋電荷か」は読出しでは決めない（辞書＝翻訳層の仕事）。
    """
    cs = g_pair_cross_spectrum(a, b, nc, ne, kmax)
    phases = {ch: v[1] for ch, v in cs["ladder"].items()}
    weights = {ch: v[0] for ch, v in cs["ladder"].items()}
    tot = sum(weights.values())
    mean_phase = (float(np.angle(sum(w * np.exp(1j * phases[ch])
                                     for ch, w in weights.items())))
                  if tot > 0 else float("nan"))
    return {"phi_ch": phases, "c_ch": weights, "phi_mean": mean_phase}


def g_pair_panel(a, b, theta, nc, ne, kmax=16):
    """二体対の常時実行パネル（v5 新設メンバーの束ね）。"""
    out = {}
    out.update(g_pair_norms(a, b))
    out.update(g_pair_overlap(a, b))
    out.update(g_pair_flow(a, b, theta))
    out.update(g_pair_charge_phase(a, b, nc, ne, kmax))
    return out


def _selftest():
    rng = np.random.default_rng(5)
    nc, ne = 8, 4
    a = (rng.normal(size=nc * ne) + 1j * rng.normal(size=nc * ne))
    b = (rng.normal(size=nc * ne) + 1j * rng.normal(size=nc * ne))
    th = 0.2
    # 恒等式: プローブ回転の直接差分と一致
    a2 = a * np.cos(th) - b * np.sin(th)
    b2 = a * np.sin(th) + b * np.cos(th)
    dNb = float(np.vdot(b2, b2).real - np.vdot(b, b).real)
    fl = g_pair_flow(a, b, th)
    assert abs(dNb - fl["dN_B_pred"]) < 1e-12 * max(1.0, abs(dNb))
    # 共通位相ゲージ不変
    g = np.exp(1j * 0.7)
    f2 = g_pair_flow(g * a, g * b, th)
    assert abs(f2["overlap_re"] - fl["overlap_re"]) < 1e-12
    print("v5 selftest OK")


if __name__ == "__main__":
    _selftest()
