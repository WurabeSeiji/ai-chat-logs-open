#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""統一場の万能読出し関数 v4 — 座標読出し（xyz・t）を追加

v3 からの分岐（2026-08-08）: **D v2 の出力を引数に取る座標読出しを新設**する。
v3 までのメンバー（f₂・f_seed・帳簿・混在度・位置スペクトル・錐成分・一段残差・
時計位相）は無改訂でそのまま引き継ぐ。

--------------------------------------------------------------------------
なぜ v4 が要るか

xyz を読むには軸と原点と位相基準が要り、それは D が供給する。v3 までの G は
フレームを必要としない量（パワー・比・位相差）しか読んでいなかったため、
D の出力を受け取る口が API に存在しなかった。**宣言していた三層構造
`確定値 = S ∘ G(·, D(·))` が実装されていなかった。** v4 でその口を作る。

--------------------------------------------------------------------------
座標の作り方（D v2 と対）

  **移動量は位相の前進そのものである。** θ = θ₀ + ωt ⇒ X = ωt。
  静止（集団時計と同期）なら移動量 0、光なら最大。

  D が返す:
    displacement : 一段の相対位相前進 δ（＝時計との離調＝移動量）
    dir_gauge    : 方向のゲージ（局所方向を大本の基準三つ組へ射影した方向余弦）
    local_gauge  : 等方ゲージ |r| = T（スカラー1個・**規格化していない**）

  G が行う:
    この歩の変位 = δ × 方向余弦（3成分）
    位置 = 前ステップの位置 + この歩の変位

  **原点は不明である。** したがって初期値を 0 と置き、以後は差だけが物理。
  本来は積分だが、進行中の波が自身の位置を覚えていることは自明であり
  （実測: キック変位は 8000 ステップ無減衰で保持）、前ステップの位置を
  持ち回して足す離散和を計算上の便宜として採る。

  **時間方向は観測できる方向を持たない。** したがって方向ベクトルを持たず、
  ゲージ（局所の位相前進レート）だけを返す。固有時間はその累積である。

--------------------------------------------------------------------------
設計規約（R1–R7 継承・R9 継承・R10 新設）

  R7 曖昧さ保存: G は選択をしない。読み値は束で返す。不在は NaN。
  R9 停止条件の外部化: 閾値・判定を持たない。計算不能なら NaN。
     前時刻・前位置が存在しない場合の None チェックだけは残す
     （閾値ではなく「前が無い」という構造的事実）。
  **R10 スケール保存（新設）**: 比だけを返してはならない。比を返すときは
     必ず分母（スケール）を同時に返す。**規格化は R 軸の消去であり、重力の
     消去である。** 座標は位相単位で返し、等方ゲージを併せて返す——両者の
     積を取るかどうかは呼び出し側（S）の宣言に委ねる。

確定値が必要なときは selection_v1.py の選択子を宣言して適用する。
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
    **比であるため R10 により分母も返す**（v3 は比だけを返していた）。"""
    Z2 = C2[:, 2, 0]
    d2 = np.float64(np.real(np.conj(Z2) @ Z2))
    Zp = Z2 - p2 * (p2 @ Z2) - q2 * (q2 @ Z2)
    with np.errstate(divide="ignore", invalid="ignore"):
        f2 = np.float64(np.real(np.conj(Zp) @ Zp)) / d2
    return {"f2": float(f2), "f2_scale": float(d2)}


def g_matter_fraction(C2):
    """物質分率 f_seed: 奇数帯（フェルミオン型）内容のパワー分率。
    **比であるため R10 により分母も返す。**"""
    P2 = np.abs(C2) ** 2
    ptot = np.float64(P2.sum())
    podd = np.float64(P2[:, _odd_mask(C2.shape[1]), :].sum())
    with np.errstate(divide="ignore", invalid="ignore"):
        f = podd / ptot
    return {"f_seed": float(f), "f_seed_scale": float(ptot)}


def g_cell_ledger(C2):
    """帳簿（帯 k × 巻き η）: 各セルの場の量と実効本数。v3 から無改訂。"""
    A2 = np.abs(C2) ** 2
    P = A2.sum(axis=0)
    S2 = (A2 ** 2).sum(axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        pr = P ** 2 / S2
    sup = (A2 > 0.0).sum(axis=0).astype(int)
    return {"cell_power": P, "cell_pr_m": pr, "cell_support": sup}


def g_species_content(C2):
    """種内容の重ね合わせ: 各関係波が同時に担うセル数（混在度）。v3 から無改訂。"""
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
    """位置読出し（束・周期型）: 双対レジスタ上の巻きモーメント族。v3 から無改訂。

    これは**レジスタ格子上の円周位置**であり、原点を要さない（偏角）。
    公開論文の 3D デモ（centroid3_v2）は同一の式を3本の空間レジスタ軸それぞれに
    適用して (x,y,z) を得ている。本テストベッドの状態は空間レジスタ軸を1本
    （帯レジスタ Nn）しか持たないため、この経路で読めるのは x のみである。
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
    """閉塞錐の成分読出し（関係波ごとの局所値・M 個）。v3 から無改訂。

    錐の恒等式 x²+y²+z² = t²+R²+Q²（R′²=R²+Q²）から **t² = R′² − m² − q²**。
    **R′² = T² であり、|r| = T が等方ゲージである**（D v2 の local_gauge と同一）。
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
    """一段残差 r（集団時計の活動度）。v3 から無改訂。"""
    ip = np.vdot(C_flat_prev, C_flat)
    with np.errstate(divide="ignore", invalid="ignore"):
        ph = ip / np.abs(ip)
        r = np.linalg.norm(C_flat - ph * C_flat_prev)
    return {"r": float(r)}


def g_clock_phase(c_gen, c_gen_prev):
    """物質時計（周期型・束）: 生成内容の一段位相前進。v3 から無改訂。"""
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


# ------------------------------------------------- 座標読出し（D の出力を引数に取る）

def g_quadrature_xyz(dpanel, prev=None):
    """**座標読出し**: D の出力から波ごとの (x, y, z) と固有時間 τ を返す。

    dpanel : D v2 の `d_panel` が返した束。必要なのは
             `dir_gauge`（方向のゲージ＝方向余弦・(M,3)）、
             `displacement`（移動量＝一段の相対位相前進・(M,)）、
             `local_gauge`（等方ゲージ |r|・(M,)）、
             `theta_adv`（局所の一段位相前進・(M,)）。
    prev   : 前ステップの持ち回し {"pos": (M,3), "tau": (M,)}。初回は None。

    構成:
      この歩の変位   step = δ × 方向余弦        （3成分）
      位置           pos  = 前の位置 + step      （初期値 0）
      固有時間       tau  = 前の tau + Δθ_local  （初期値 0）

    **原点は不明であるため初期値を 0 と置く**（宣言された規約）。以後、絶対値に
    物理的意味はなく差だけが物理である。座標は位相単位で返し、**等方ゲージを
    併せて返す**（R10）——両者を掛けて長さにするかどうかは S の宣言に委ねる。

    **時間は方向を持たない**（観測できる方向がない）ため、方向余弦を持たず
    ゲージ（局所の位相前進）の累積だけを返す。
    """
    dg = np.asarray(dpanel["dir_gauge"])            # (M,3)
    disp = np.asarray(dpanel["displacement"])       # (M,)
    th_adv = np.asarray(dpanel["theta_adv"])        # (M,)
    step = disp[:, None] * dg                       # (M,3)

    if prev is None:                                # 前が無い＝構造的事実（R9）
        pos = np.zeros_like(step)
        tau = np.zeros_like(th_adv)
    else:
        pos = np.asarray(prev["pos"]) + step
        tau = np.asarray(prev["tau"]) + th_adv

    return {"xyz": pos, "xyz_step": step,
            "x": pos[:, 0], "y": pos[:, 1], "z": pos[:, 2],
            "tau": tau, "t_gauge": th_adv,
            "gauge_iso": np.asarray(dpanel["local_gauge"]),
            "_carry_pos": {"pos": pos.copy(), "tau": tau.copy()}}


# ---------------------------------------------------------------- パネル

def g_panel(C2, p2, q2, C_flat_prev=None, c_gen_prev=None, dpanel=None,
            prev_pos=None):
    """常時実行パネル（v4）: 宇宙の第0步から毎ステップ・一様ケイデンス・切替なし。

    dpanel を渡すと座標読出し（xyz・τ）を併せて実行する。渡さなければ v3 と
    同じ内容を返す（座標は D の出力なしには定義できないため）。
    返り値は曖昧さを保存した束（確定値が要る場合は選択層 S を宣言して適用）。
    """
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
    carry = {"C_flat": C_flat.copy(), "c_gen": c_gen.copy()}
    if dpanel is not None:
        qx = g_quadrature_xyz(dpanel, prev_pos)
        carry["pos"] = qx["_carry_pos"]
        out.update({k: v for k, v in qx.items() if k != "_carry_pos"})
    out["_carry"] = carry
    return out
