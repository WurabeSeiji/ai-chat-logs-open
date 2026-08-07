#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""統一場の万能読出し関数 v3 — 分岐版（以後の実験はこれを使う）

v2 からの分岐（2026-08-08・木原指示）: **無駄な IF 文と閾値ガードを全廃**する。
計算不能なら NaN をそのまま返す。条件判断（存在するか・読めたか・どこで切るか）は
**図化・分析を行う上位側の責務**であり、読出し関数の責務ではない。
上位プログラムが NaN でアボートしても、上位側で対応する。

**過去のバージョンは残す**（v1・v2）。公開済み論文の再現性はそちらで確保され、
本 v3 は今後の実験専用の分岐である。

--------------------------------------------------------------------------
v2 → v3 の差分（撤廃した箇所）

| 箇所 | v2（撤廃前） | v3 |
|---|---|---|
| g_space_history | `/ max(d2, 1e-300)` | そのまま除算（0/0→NaN） |
| g_matter_fraction | `/ max(ptot, 1e-300)` | 同上 |
| g_cell_ledger | `np.where(S2>0, …, 0.0)` | そのまま除算（NaN） |
| g_species_content | `np.where(…, 0.0)`・`if w>0 else 0.0` | そのまま除算（NaN） |
| g_position_spectrum | `if tot<=0: 零ベクトルと pr_n=0 を返す` | 分岐撤廃（NaN） |
| g_collective_residual | `ph = ip/abs(ip) if abs(ip)>0 else **1.0**` | そのまま除算（NaN）。**単位位相の捏造をやめた** |
| g_clock_phase | `if ov>0 else nan`・`coherence … else 0.0` | 分岐撤廃（NaN） |

None チェック（前時刻が存在しない）だけは残す——これは閾値ではなく
「前の状態が無い」という構造的事実であり、NaN を返す。
--------------------------------------------------------------------------

設計規約（R1–R7 を継承し、D と共通の R9 を追加）:
  R1 純関係量  R2 パラメータフリー  R3 種別分岐なし  R4 全時代同一定義
  R5 常時実行  R6 受動性
  R7 曖昧さ保存: G は選択をしない。読み値は (値, 重み) の束で返す。不在は
     真偽値でなく重み 0（または NaN）で表す。周期量は円上の値で返す。
  R9 停止条件の外部化: 閾値・判定・停止条件を持たない。計算不能なら NaN。

確定値が必要なときは selection_v1.py の選択子を宣言して適用する
（確定値 = S ∘ G(·, D(·))）。
"""
from __future__ import annotations
import numpy as np

_NAN = float("nan")


def _odd_mask(Nn):
    ks = np.arange(Nn)
    return (ks % 2 == 1)


# ---------------------------------------------------------------- 単時刻メンバー

def g_invariants(C2):
    """不変量族: 総パワー（振幅保存の監査対象）"""
    return {"P_tot": float(np.sum(np.abs(C2) ** 2))}


def g_space_history(C2, p2, q2):
    """空間形成史 f₂: 親平面（スライス2・巻き0）の面外分率。
    親平面の内容が無ければ 0/0 → NaN（v2 の 1e-300 ガードを撤廃）。"""
    Z2 = C2[:, 2, 0]
    d2 = np.float64(np.real(np.conj(Z2) @ Z2))
    Zp = Z2 - p2 * (p2 @ Z2) - q2 * (q2 @ Z2)
    with np.errstate(divide="ignore", invalid="ignore"):
        f2 = np.float64(np.real(np.conj(Zp) @ Zp)) / d2
    return {"f2": float(f2)}


def g_matter_fraction(C2):
    """物質分率 f_seed: 奇数帯（フェルミオン型）内容のパワー分率。
    総パワー 0 なら NaN（v2 の 1e-300 ガードを撤廃）。"""
    P2 = np.abs(C2) ** 2
    ptot = np.float64(P2.sum())
    podd = np.float64(P2[:, _odd_mask(C2.shape[1]), :].sum())
    with np.errstate(divide="ignore", invalid="ignore"):
        f = podd / ptot
    return {"f_seed": float(f)}


def g_cell_ledger(C2):
    """帳簿（帯 k × 巻き η）: 各セルの場の量と実効本数。

    cell_power  : セルのパワー（唯一の一次量）
    cell_pr_m   : 実効本数 PR_M（M方向の参加比・連続量）。パワー 0 のセルは
                  0/0 → **NaN**（v2 は 0 を代用していた）
    cell_support: 台集合の濃度（厳密非零判定・**個数ではない**——監査で
                  count ≡ M×(power>0) と実測され情報量ゼロ）
    """
    A2 = np.abs(C2) ** 2
    P = A2.sum(axis=0)
    S2 = (A2 ** 2).sum(axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        pr = P ** 2 / S2
    sup = (A2 > 0.0).sum(axis=0).astype(int)
    return {"cell_power": P, "cell_pr_m": pr, "cell_support": sup}


def g_species_content(C2):
    """種内容の重ね合わせ: 各関係波が同時に担うセル数（混在度）。
    パワー 0 の波は NaN（v2 は 0 を代用していた）。"""
    A2 = np.abs(C2) ** 2
    F = A2.reshape(A2.shape[0], -1)
    P = F.sum(axis=1)
    S2 = (F ** 2).sum(axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        mix = P ** 2 / S2
        mix_mean = (mix * P).sum() / P.sum()
    return {"per_wave_power": P, "per_wave_mix": mix,
            "mix_mean": float(mix_mean)}


def g_position_spectrum(C2):
    """位置読出し（束・周期型）: 双対レジスタ上の巻きモーメント族。

    奇数帯内容の双対分布 P(n) の全巻き m=1…⌊Nn/2⌋ のモーメントを返す。
    内容が無ければ全て NaN（v2 の `if tot<=0` 分岐と零代用を撤廃）。
    どの m を採るかは選択であり G は選択しない（R7）。
    """
    Nn, Neta = C2.shape[1], C2.shape[2]
    mask = _odd_mask(Nn).astype(float)
    Wo = np.fft.ifft2(C2 * mask[None, :, None], axes=(1, 2)) * (Nn * Neta)
    Pn = np.sum(np.abs(Wo) ** 2, axis=(0, 2))
    tot = np.float64(Pn.sum())
    nn = np.arange(Nn)
    ms = np.arange(1, Nn // 2 + 1)
    with np.errstate(divide="ignore", invalid="ignore"):
        z = np.array([np.sum(Pn * np.exp(2j * np.pi * m * nn / Nn)) / tot
                      for m in ms])
        pr_n = float(tot ** 2 / np.sum(Pn ** 2))
        x = (np.angle(z) * Nn / (2 * np.pi * ms)) % (Nn / ms)
    return {"pos_m": ms, "pos_x": x, "pos_weight": np.abs(z),
            "pos_modulus": Nn / ms, "content_power": float(tot), "pr_n": pr_n,
            "dual_profile": Pn}


def g_cone_components(C2):
    """閉塞錐の成分読出し（関係波ごとの局所値・M 個）。

    錐の恒等式 x²+y²+z² = t²+R²+Q²（R′²=R²+Q²）から **t² = R′² − m² − q²**。
    座標時間は残余として読む。**粒度は関係波ごと**——大域一値にすると重力
    （ゲージの目盛の不等間隔）が定義上消えるため。

    Gram はエンジン自身の二チャネル分割（帯の偶奇対 k=2j / 2j+1）で作る。
      T=S0, X=S1, Y=S2, Z=S3、R′²=T²、m²=detΓ=T²−X²−Y²−Z²（≥0 は
      Cauchy–Schwarz より自動）、q=Σ_η η_signed·P(η)、t²=R′²−m²−q²。
    本メンバーはスピン読出しではない（(X,Y,Z) をスピンとは主張しない）。
    """
    M, Nn, Neta = C2.shape
    J = Nn // 2
    ke = 2 * np.arange(J)
    ko = ke + 1
    a = C2[:, ko, :].reshape(M, -1)
    b = C2[:, ke, :].reshape(M, -1)
    pa = np.sum(np.abs(a) ** 2, axis=1)
    pb = np.sum(np.abs(b) ** 2, axis=1)
    z = np.sum(np.conj(a) * b, axis=1)
    T = pa + pb
    X = pa - pb
    Y = 2.0 * np.real(z)
    Z = 2.0 * np.imag(z)
    m2 = T ** 2 - X ** 2 - Y ** 2 - Z ** 2
    eta = np.arange(Neta)
    eta_s = ((eta + Neta // 2) % Neta) - Neta // 2
    P_eta = (np.abs(C2[:, ko, :]) ** 2 + np.abs(C2[:, ke, :]) ** 2).sum(axis=1)
    q = P_eta @ eta_s.astype(float)
    Rp2 = T ** 2
    return {"cone_T": T, "cone_X": X, "cone_Y": Y, "cone_Z": Z,
            "cone_Rp2": Rp2, "cone_m2": m2, "cone_q": q, "cone_q2": q ** 2,
            "cone_t2": Rp2 - m2 - q ** 2}


# ---------------------------------------------------------------- 二時刻メンバー

def g_collective_residual(C_flat, C_flat_prev):
    """一段残差 r（集団時計の活動度）: 大域位相を合わせた上での状態変化ノルム。
    **v2 は重なりが 0 のとき単位位相 1.0 を捏造していた**。v3 は除算をそのまま
    行い NaN を返す（読めないものを読めたことにしない）。"""
    ip = np.vdot(C_flat_prev, C_flat)
    with np.errstate(divide="ignore", invalid="ignore"):
        ph = ip / np.abs(ip)
        r = np.linalg.norm(C_flat - ph * C_flat_prev)
    return {"r": float(r)}


def g_clock_phase(c_gen, c_gen_prev):
    """物質時計（周期型・束）: 生成内容の一段位相前進。
    担体が無ければ NaN（v2 の分岐と 0 代用を撤廃）。前時刻が無い場合のみ
    None チェックを残す（閾値ではなく構造的事実）。"""
    cp = float(np.real(np.vdot(c_gen, c_gen)))
    if c_gen_prev is None:
        return {"phase": _NAN, "overlap": _NAN, "carrier_power": cp,
                "coherence": _NAN}
    zz = np.vdot(c_gen_prev, c_gen)
    n0 = np.float64(np.linalg.norm(c_gen_prev))
    n1 = np.float64(np.linalg.norm(c_gen))
    ov = np.float64(abs(zz))
    with np.errstate(divide="ignore", invalid="ignore"):
        coh = ov / (n0 * n1)
    return {"phase": float(np.angle(zz)), "overlap": float(ov),
            "carrier_power": cp, "coherence": float(coh)}


# ---------------------------------------------------------------- パネル

def g_panel(C2, p2, q2, C_flat_prev=None, c_gen_prev=None):
    """常時実行パネル（v3）: 宇宙の第0步から毎ステップ・一様ケイデンス・切替なし。
    返り値は曖昧さを保存した束（確定値が要る場合は選択層 S を宣言して適用）。"""
    out = {}
    out.update(g_invariants(C2))
    out.update(g_space_history(C2, p2, q2))
    out.update(g_matter_fraction(C2))
    out.update(g_cell_ledger(C2))
    out.update(g_species_content(C2))
    out.update(g_position_spectrum(C2))
    out.update(g_cone_components(C2))
    C_flat = C2.reshape(-1)
    if C_flat_prev is not None:
        out.update(g_collective_residual(C_flat, C_flat_prev))
    om = _odd_mask(C2.shape[1]).copy()
    om[1] = False           # 生成帯＝最低奇数帯（シード帯 k=1）を除く奇数帯
    c_gen = C2[:, om, :].reshape(-1)
    out.update(g_clock_phase(c_gen, c_gen_prev))
    out["_carry"] = {"C_flat": C_flat.copy(), "c_gen": c_gen.copy()}
    return out
