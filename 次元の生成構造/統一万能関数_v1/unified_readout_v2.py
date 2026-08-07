#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""統一場の万能読出し関数 v2 — 曖昧さ保存版（以後の全実験はこのモジュールを使う）

v1 からの設計変更（監査 2026-08-08 の結果・是正記録は §「v1 の欠陥」）:

  確定値 = S（選択層・宣言される） ∘ G（曖昧さを保存した読出し）

  v1 の G は「曖昧さの読出し」と「どれを採るかの選択」を一体に焼き込み、
  確定値をひとつだけ返していた（位置の被覆判定 if |z2|>|z1|、存在の真偽化、
  時計の取得可否の真偽化、波数の非零指示）。この設計では
  「どの選択の下での値か」を宣言できず、基準値が再現不能になる。
  v2 は G から選択を全て外に出し、**重みつきの複数読み値（束）** を返す。
  従来の確定値は selection_v1.py の選択子を適用して厳密に再現できる
  （資格審査 Q6 でビット一致を確認）。

設計規約（v1 の R1–R6 を継承し R7 を追加）:
  R1 純関係量  R2 パラメータフリー  R3 種別分岐なし  R4 全時代同一定義
  R5 常時実行  R6 受動性
  R7 曖昧さ保存: G は選択をしない。分岐・閾値・二値化・順序づけ・単一化を
     行わず、読み値は (値, 重み) の束として返す。不在は真偽値でなく
     「重み 0」で表す。周期量は円上の値として返し、直線への展開（巻数の
     累積）は呼び出し側の責務とする。縮約は関係量への集約までとし、
     波の添字・セルの添字は保持してよい（演算対象の構造だから）。

参照: 統一万能関数_仕様_v1.md／資格審査 run_qualification_readout_v2.py
"""
from __future__ import annotations
import numpy as np


def _odd_mask(Nn):
    ks = np.arange(Nn)
    return (ks % 2 == 1)


# ---------------------------------------------------------------- 単時刻メンバー

def g_invariants(C2):
    """不変量族: 総パワー（振幅保存の監査対象）"""
    return {"P_tot": float(np.sum(np.abs(C2) ** 2))}


def g_space_history(C2, p2, q2):
    """空間形成史 f₂: 親平面（スライス2・巻き0）の面外分率。
    比であり選択を含まない（v1 と同一定義）。"""
    Z2 = C2[:, 2, 0]
    d2 = float(np.real(np.conj(Z2) @ Z2))
    Zp = Z2 - p2 * (p2 @ Z2) - q2 * (q2 @ Z2)
    return {"f2": float(np.real(np.conj(Zp) @ Zp)) / max(d2, 1e-300)}


def g_matter_fraction(C2):
    """物質分率 f_seed: 奇数帯（フェルミオン型）内容のパワー分率。"""
    P2 = np.abs(C2) ** 2
    ptot = float(P2.sum())
    podd = float(P2[:, _odd_mask(C2.shape[1]), :].sum())
    return {"f_seed": podd / max(ptot, 1e-300)}


def g_cell_ledger(C2):
    """帳簿（帯 k × 巻き η）: 各セルの場の量と、そこに参加する関係波の実効本数。

    power  : セルのパワー（連続量・これが唯一の一次量）
    pr_m   : 実効本数 PR_M = (Σ_e|c|²)²/Σ_e|c|⁴（M方向の参加比・連続量）
             ——「何本の波がこのセルの内容を担っているか」の実効値。

    【v1 の欠陥・是正記録】v1 の g_wave_census は本数を
    count = #(振幅² > 0) の厳密非零判定で返していた。監査（N=12・Nn=5）で
    到達可能な全セルが例外なく count = M = 66 となり（パワー 1.0 の凝縮体も
    9.1e-22 の数値塵も同じ 66）、count ≡ M×(パワー>0 の指示関数) であって
    帳簿の台集合以外の情報を持たないことが実測された。個数の量としては
    機能しないため v2 では返さない。個数を数える量は
      (a) 実効本数 pr_m（連続量・G が返す）
      (b) 占有数 n = power / ε（ε は 1 局在ドメインあたりのパワー・
          較正が必要なので選択層 selection_v1.s_occupancy の責務）
    の二つに分離した。台集合の指示が要る実験は
    selection_v1.s_support_count を宣言して使うこと。
    """
    A2 = np.abs(C2) ** 2                      # M×Nn×Nη
    P = A2.sum(axis=0)
    S2 = (A2 ** 2).sum(axis=0)
    pr = np.where(S2 > 0.0, P ** 2 / np.where(S2 > 0.0, S2, 1.0), 0.0)
    sup = (A2 > 0.0).sum(axis=0).astype(int)   # 台集合の濃度（厳密零判定・閾値なし）
    return {"cell_power": P, "cell_pr_m": pr, "cell_support": sup}


def g_species_content(C2):
    """種内容の重ね合わせ: 各関係波が、いくつのセル（＝粒子種の内容）を
    同時に担っているか。

    per_wave_power : 波ごとの総パワー（M,）
    per_wave_mix   : 波ごとの実効セル数 PR_cell = (ΣP)²/ΣP²（M,）
                     ——1 なら単一種、大きいほど多種の重ね合わせ。
    mix_mean       : per_wave_mix のパワー重み平均（系全体の混在度）

    「一つの波は、読出しの段階では多くの粒子として見える内容の
    重ね合わせである」——この命題を、注釈でなく測定量として返す。
    """
    A2 = np.abs(C2) ** 2                      # M×Nn×Nη
    F = A2.reshape(A2.shape[0], -1)           # M×(Nn·Nη)
    P = F.sum(axis=1)
    S2 = (F ** 2).sum(axis=1)
    mix = np.where(S2 > 0.0, P ** 2 / np.where(S2 > 0.0, S2, 1.0), 0.0)
    w = P.sum()
    return {"per_wave_power": P, "per_wave_mix": mix,
            "mix_mean": float((mix * P).sum() / w) if w > 0 else 0.0}


def g_position_spectrum(C2):
    """位置読出し（束・周期型）: 双対レジスタ上の巻きモーメント族。

    奇数帯（フェルミオン型）内容の双対分布 P(n) に対し、全ての巻き
    m = 1 … ⌊Nn/2⌋ のモーメント z_m = Σ_n P(n)e^{2πi m n/Nn} / ΣP を返す。

    x_m     : m 番目の読み値 = arg(z_m)·Nn/(2πm) mod (Nn/m)   ← 円上の値
    modulus : Nn/m（その読み値の周期＝法）
    weight  : |z_m|（その読み値の重み・鋭さ）
    content_power : 奇数帯内容の総パワー（0 なら不在——真偽値は返さない）
    pr_n    : 双対占有の実効セル数（局在度・連続量）

    どの m を採るか（＝空間の被覆をいくつと読むか）は選択であり、G は
    選択しない。従来の被覆判定（|z2|>|z1|）は
    selection_v1.s_position_maxmoment で厳密に再現できる。
    """
    Nn, Neta = C2.shape[1], C2.shape[2]
    mask = _odd_mask(Nn).astype(float)
    Wo = np.fft.ifft2(C2 * mask[None, :, None], axes=(1, 2)) * (Nn * Neta)
    Pn = np.sum(np.abs(Wo) ** 2, axis=(0, 2))
    tot = float(Pn.sum())
    nn = np.arange(Nn)
    ms = np.arange(1, Nn // 2 + 1)
    if tot <= 0.0:
        z = np.zeros(len(ms), complex)
        pr_n = 0.0
    else:
        z = np.array([np.sum(Pn * np.exp(2j * np.pi * m * nn / Nn)) / tot
                      for m in ms])
        pr_n = float(tot ** 2 / np.sum(Pn ** 2))
    with np.errstate(invalid="ignore", divide="ignore"):
        x = np.where(np.abs(z) > 0.0,
                     (np.angle(z) * Nn / (2 * np.pi * ms)) % (Nn / ms), np.nan)
    return {"pos_m": ms, "pos_x": x, "pos_weight": np.abs(z),
            "pos_modulus": Nn / ms, "content_power": tot, "pr_n": pr_n,
            "dual_profile": Pn}


def g_cone_components(C2):
    """閉塞錐の成分読出し（関係波ごとの局所値・M 個）。

    錐の恒等式（[3]・空間軸3軸と固有時間の創生）:
        x² + y² + z² = t² + R² + Q²  ,  R′² = R² + Q²
        ⇒ **t² = R′² − m² − q²**（座標時間は残余として読む）

    **粒度は関係波 e ごと（M 個）である。** 宇宙全体で一つの値にしてはならない
    ——重力はこのモデルでは「ゲージ（空間・時間の目盛）が等間隔でなくなること」
    として現れるので、t を大域一値に定義すると重力の効果が定義上消える。
    局所的に異なる時計の進み方こそが測定対象であり、M 本の分布の広がりが
    重力赤方偏移に相当する。

    Gram（コヒーレンス行列）の構成——エンジン自身の二チャネル分割を使う:
      帯の偶奇対 (k_even=2j, k_odd=2j+1), j=0…⌊Nn/2⌋−1 を対にし、
      a = 奇数帯（フェルミオン型）内容、b = 偶数帯（ボゾン型）内容とする。
      これは相互作用関数 F が生成子を作るときに使う分割そのものであり、
      新たな任意構造を持ち込んでいない（R3 種別分岐なしを侵さない）。
      Nn が奇数のとき対にならない最上帯は除外する。
      S0=Σ(|a|²+|b|²)  S1=Σ(|a|²−|b|²)  S2=2ReΣā b  S3=2ImΣā b
      (T,X,Y,Z) := (S0,S1,S2,S3)

    返す量（全て (M,) 配列・パワー² の次元）:
      cone_T    : T = 対内容の総パワー（保存読出し）
      cone_X/Y/Z: Gram の Pauli 成分
      cone_Rp2  : R′² = T²
      cone_m2   : m² = detΓ = T²−X²−Y²−Z²（質量²＝非コヒーレンス）
                  Cauchy–Schwarz より恒等的に ≥0（光錐束縛が自動）
      cone_q    : q = Σ_η η_signed·P(η)（符号付き巻きの重みつき和・電荷）
                  η_signed = ((η+Nη/2) mod Nη) − Nη/2。中性レシピでは厳密 0。
      cone_q2   : q²
      cone_t2   : t² = R′² − m² − q²（座標時間の二乗・残余読出し）

    注意（正直な限界）: 本メンバーは **スピン読出しではない**。Stokes 型の
    構成を Gram 4元を作るためだけに用いており、(X,Y,Z) をスピンと解釈する
    ことは主張しない（統一エンジン上のスピン同定は TB-spin として未確定）。
    また t² ≥ 0 は中性（q=0）では detΓ の定義から自動だが、荷電種が居る場合は
    自動でなく、**そこで初めて非自明な検定になる**。
    """
    M, Nn, Neta = C2.shape
    J = Nn // 2
    ke = 2 * np.arange(J)          # 偶数帯（ボゾン型）
    ko = ke + 1                    # 奇数帯（フェルミオン型）
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
    """一段残差 r（集団時計の活動度）: 大域位相を合わせた上での状態変化ノルム。"""
    ip = np.vdot(C_flat_prev, C_flat)
    ph = ip / abs(ip) if abs(ip) > 0 else 1.0
    return {"r": float(np.linalg.norm(C_flat - ph * C_flat_prev))}


def g_clock_phase(c_gen, c_gen_prev):
    """物質時計（周期型・束）: 生成内容の一段位相前進。

    phase          : arg⟨c_prev, c⟩（円上の値・巻数の累積は呼び出し側の責務）
    overlap        : |⟨c_prev, c⟩|（重み・0 なら担体不在）
    carrier_power  : ‖c‖²（時計の担い手＝質量の量）
    coherence      : |⟨c_prev,c⟩|/(‖c_prev‖‖c‖)（読み値の鋭さ）

    取得可否の真偽判定（v1 の acquirable）は選択層
    selection_v1.s_clock_acquirable の責務。
    """
    if c_gen_prev is None:
        return {"phase": float("nan"), "overlap": 0.0,
                "carrier_power": float(np.vdot(c_gen, c_gen).real),
                "coherence": 0.0}
    zz = np.vdot(c_gen_prev, c_gen)
    n0 = float(np.linalg.norm(c_gen_prev)); n1 = float(np.linalg.norm(c_gen))
    ov = float(abs(zz))
    return {"phase": float(np.angle(zz)) if ov > 0.0 else float("nan"),
            "overlap": ov, "carrier_power": n1 * n1,
            "coherence": ov / (n0 * n1) if n0 > 0 and n1 > 0 else 0.0}


# ---------------------------------------------------------------- パネル（束ね）

def g_panel(C2, p2, q2, C_flat_prev=None, c_gen_prev=None):
    """常時実行パネル（v2）: 宇宙の第0步から毎ステップ・一様ケイデンス・切替なし。
    返り値は全て曖昧さを保存した束（確定値が要る場合は選択層を適用する）。"""
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
    om[1] = False          # 生成帯＝最低奇数帯（シード帯 k=1）を除く奇数帯
    c_gen = C2[:, om, :].reshape(-1)
    out.update(g_clock_phase(c_gen, c_gen_prev))
    out["_carry"] = {"C_flat": C_flat.copy(), "c_gen": c_gen.copy()}
    return out
