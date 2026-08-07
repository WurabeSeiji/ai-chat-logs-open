#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""統一場の万能読出し関数 v1 — 以後の全実験はこのモジュールのみを使う（木原ルール）

資格審査: run_qualification_v1.py（正本実装との数値一致・受動性・全通過が使用条件）。
各メンバーの正本出典:
  f₂（空間形成史）・一段残差 r（集団時計）= run_genesis_v1（数値一致 Q3）
  f_seed（物質分率）・census 帳簿・物質時計 ω̂ = stage3 census 系譜（TB v2 で実証）
  位置(1D) = 双対レジスタ巻きモーメント（3Dデモ centroid3 の1D版・
             二重被覆は mod Nn/2 で読む——周期表柱7）
  スピンStokes = AC-0 新設（リングエンジン4判定通過・統一エンジンでは対象
             未同定＝TB-spin 試験まで統一エンジン状態への適用禁止）

設計規約（違反はテストベッド資格の喪失）:
  R1 純関係量: 全メンバーは状態の内積・パワー・偏角のみから成る。
  R2 パラメータフリー: 調整定数・閾値・窓関数を含まない（判定閾値は
     実験スクリプト側の責務であり G の責務ではない）。
  R3 種別分岐なし: 粒子種・時代・N による場合分けを持たない（無名性）。
  R4 全時代同一定義: 宇宙誕生の第0步から終端まで同一の式。未成立の量は
     NaN でなく「曖昧な値」または不在指標として正直に返す。
  R5 常時実行: 毎ステップ・一様ケイデンスで評価される前提で設計する。
     タイミングを図った適用は無名性違反。
  R6 受動性: 状態を読むだけで書かない（力学はビット単位で不変）。

対象エンジン: HairEngine 系（状態 C2: M×Nn×Nη 複素・M=N(N-1)/2 関係波・
  Nn=帯レジスタ・Nη=毛レジスタ）。他エンジンへの移植は同名メンバーの
  再導出と対照テストを要する（エンジン間整合は登録済み未解決）。

戻り値の規約: 全メンバー dict。二時刻メンバー（ω̂ 等）は前回状態を
  呼び出し側が保持し引数で渡す（G自身は履歴を持たない＝無状態）。
"""
from __future__ import annotations
import numpy as np


# ---------------------------------------------------------------- 単時刻メンバー

def g_invariants(C2):
    """不変量族: ノルム（振幅保存の監査対象）"""
    P_tot = float(np.sum(np.abs(C2) ** 2))
    return {"P_tot": P_tot}


def g_space_history(C2, p2, q2):
    """空間形成史 f₂: 親平面（スライス2・毛0）の面外分率。
    真空期 ≈0 → インフレーションで指数増幅 → 凝縮で飽和。
    p2, q2 は誕生時の親平面基底（第0步で一度だけ構成する参照系）。"""
    Z2 = C2[:, 2, 0]
    d2 = float(np.real(np.conj(Z2) @ Z2))
    Zp = Z2 - p2 * (p2 @ Z2) - q2 * (q2 @ Z2)
    f2 = float(np.real(np.conj(Zp) @ Zp)) / max(d2, 1e-300)
    return {"f2": f2}


def _odd_mask(Nn):
    ks = np.arange(Nn)
    return (ks % 2 == 1)


def g_matter_fraction(C2):
    """物質分率 f_seed: 奇数帯（フェルミオン型）内容のパワー分率。
    帯マスクは汎用（全奇数k）——Nn=5 では {1,3} で従来定義とビット同一。"""
    P2 = np.abs(C2) ** 2
    ptot = float(P2.sum())
    podd = float(P2[:, _odd_mask(C2.shape[1]), :].sum())
    return {"f_seed": podd / max(ptot, 1e-300)}


def g_census(C2):
    """粒子census帳簿: 帯 k × 毛 η のパワー行列（何が・どれだけ居るか）。
    毛=電荷巻き（Q=η/3 の読みは周期表柱2の後段読出し）。"""
    return {"census": np.abs(C2).__pow__(2).sum(axis=0)}


def g_wave_census(C2):
    """波数census: セル（帯k×巻きη）ごとの参加関係波の勘定。
    粒子＝波・個数＝波の本数（周期表の数え上げ規定）。全て純関係量:
      count = 厳密非零の参加本数（振幅²>0 の厳密判定・閾値なし——
              荷電種の厳密0はここで厳密0本と数えられる）
      pr_m  = 実効本数（参加比 PR_M = (Σ|c|²)²/Σ|c|⁴・M方向）"""
    A2 = np.abs(C2) ** 2                       # M×Nn×Nη
    P = A2.sum(axis=0)
    count = (A2 > 0.0).sum(axis=0).astype(int)
    S2 = (A2 ** 2).sum(axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        pr = np.where(S2 > 0.0, P ** 2 / np.where(S2 > 0.0, S2, 1.0), 0.0)
    return {"wave_count": count, "wave_pr_m": pr}


def g_position_1d(C2):
    """位置読出し（1D・暫定）: 奇数帯内容の双対レジスタ重心（巻きモーメント偏角）。
    注意: フェルミオン優勢内容は対蹠二点構造（空間的二重被覆・周期表柱7）を
    持ち得るため、位置は被覆に応じ mod Nn/2 で読むこと。3D版（x,y,z）は
    3Dレジスタエンジンのみ対応（本メンバーの既知の制限）。内容が無い場合は
    存在指標 present=False を返す（NaN でなく不在の明示）。"""
    Nn, Neta = C2.shape[1], C2.shape[2]
    mask = _odd_mask(Nn).astype(float)
    Wo = np.fft.ifft2(C2 * mask[None, :, None], axes=(1, 2)) * (Nn * Neta)
    Pn = np.sum(np.abs(Wo) ** 2, axis=(0, 2))
    if Pn.sum() <= 1e-280:
        return {"x": None, "cover": None, "present": False, "pr_n": 0.0}
    pr_n = float(Pn.sum() ** 2 / np.sum(Pn ** 2))   # 双対占有の実効セル数（局在度）
    nn = np.arange(Nn)
    z1 = np.sum(Pn * np.exp(2j * np.pi * nn / Nn)) / Pn.sum()
    z2 = np.sum(Pn * np.exp(2j * np.pi * 2 * nn / Nn)) / Pn.sum()
    # 被覆判定は関係量のみ: 対蹠二点構造（奇数帯の反周期性）では z1 が相殺し
    # z2 が残る（正本 centroid3_v2 軸1 の処方・周期表柱7）。
    if abs(z2) > abs(z1):
        x = float((np.angle(z2) * Nn / (4 * np.pi)) % (Nn / 2))
        return {"x": x, "cover": 2, "present": True, "pr_n": pr_n}
    x = float((np.angle(z1) * Nn / (2 * np.pi)) % Nn)
    return {"x": x, "cover": 1, "present": True, "pr_n": pr_n}


# ---------------------------------------------------------------- 二時刻メンバー

def g_collective_residual(C_flat, C_flat_prev):
    """一段残差 r（集団時計の活動度）: 大域位相を合わせた上での状態変化ノルム。"""
    ip = np.vdot(C_flat_prev, C_flat)
    ph = ip / abs(ip) if abs(ip) > 0 else 1.0
    return {"r": float(np.linalg.norm(C_flat - ph * C_flat_prev))}


def g_matter_clock(c_gen, c_gen_prev):
    """物質時計 ω̂: 生成内容（相棒帯 k=3）の一段位相前進。
    固有時間の計測可能性そのもの——担体が無い（真空宇宙）では取得不能を
    不在として返す。安定期の期待値は柱1の集団時計 π/72（実測で確認済み）。"""
    if c_gen_prev is None:
        return {"omega": None, "acquirable": False}
    zz = np.vdot(c_gen_prev, c_gen)
    if abs(zz) <= 1e-30 or np.linalg.norm(c_gen) <= 1e-12:
        return {"omega": None, "acquirable": False}
    return {"omega": float(np.angle(zz)), "acquirable": True}


# ---------------------------------------------------------------- 未確定メンバー

def g_spin_stokes_candidate(a2, b2):
    """スピン読出し（Stokes射影族）——二チャネル (a,b) エンジン用の設計。
    リングエンジンでは検定済み（AC-0: 灰色誕生 σ1=σ2=0 厳密・灰色不変性
    1.7e-17・受動性ビット同一）。テストベッドエンジン（単チャネル関係波）
    での対応物は【未確定・スピン可読性試験 TB-spin に登録】——本メンバーを
    テストベッド状態に適用してはならない（対象構造が未同定）。"""
    S0 = np.sum(np.abs(a2) ** 2 + np.abs(b2) ** 2, axis=1)
    S1 = np.sum(np.abs(a2) ** 2 - np.abs(b2) ** 2, axis=1)
    Z = np.sum(np.conj(a2) * b2, axis=1)
    d = np.maximum(S0, 1e-300)
    return {"sigma1": S1 / d, "sigma2": 2.0 * np.real(Z) / d,
            "sigma3": 2.0 * np.imag(Z) / d}


# ---------------------------------------------------------------- パネル（束ね）

def g_panel(C2, p2, q2, C_flat_prev=None, c_gen_prev=None):
    """常時実行パネル: 全単時刻メンバー＋（履歴があれば）二時刻メンバー。
    呼び出し規約: 宇宙の第0步から毎ステップ・一様ケイデンス・切替なし。"""
    out = {}
    out.update(g_invariants(C2))
    out.update(g_space_history(C2, p2, q2))
    out.update(g_matter_fraction(C2))
    out.update(g_census(C2))
    out.update(g_wave_census(C2))
    out.update(g_position_1d(C2))
    C_flat = C2.reshape(-1)
    if C_flat_prev is not None:
        out.update(g_collective_residual(C_flat, C_flat_prev))
    om = _odd_mask(C2.shape[1]).copy()
    om[1] = False   # 生成帯=シード帯(最低奇数帯k=1)を除く奇数帯（Nn=5ではk=3と同一）
    c_gen = C2[:, om, :].reshape(-1)
    out.update(g_matter_clock(c_gen, c_gen_prev))
    out["_carry"] = {"C_flat": C_flat.copy(), "c_gen": c_gen.copy()}
    return out


# ---------------------------------------------------------------- 自己検定

def _selftest():
    """対照テスト: run_tb_genesis_universe_v2 のインライン実装と数値一致すること
    （短走行・N=4）。仕様と実測の同一性の監査。"""
    import importlib.util, sys
    from pathlib import Path
    here = Path(__file__).resolve().parent
    spec = importlib.util.spec_from_file_location(
        "ui_ref", here / "unified_interaction_v1.py")
    ui = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = ui
    spec.loader.exec_module(ui)
    eng, p2, q2 = ui.build_standard_universe(4, 1e-2)
    carry = {"C_flat": None, "c_gen": None}
    ok = True
    for t in range(50):
        eng.step()
        C2 = eng.C2()
        p = g_panel(C2, p2, q2, carry["C_flat"], carry["c_gen"])
        carry = p["_carry"]
        # インライン基準値
        Z2 = C2[:, 2, 0]
        Zp = Z2 - p2 * (p2 @ Z2) - q2 * (q2 @ Z2)
        f2_ref = float(np.real(np.conj(Zp) @ Zp)) / float(np.real(np.conj(Z2) @ Z2))
        if abs(p["f2"] - f2_ref) > 1e-15:
            ok = False
    print("selftest:", "PASS" if ok else "FAIL")
    return ok


if __name__ == "__main__":
    _selftest()
