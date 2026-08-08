#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""統一万能次元読出し関数 D v2 — 方向とゲージを返す（v1 から全面改訂）

階層:  F（力学）→ **D（次元・フレーム読出し）** → G（量の読出し）→ S（選択層）

--------------------------------------------------------------------------
v1 からの改訂理由（2026-08-08・実測に基づく）

v1 は「フレーム＝凝縮体から作った瞬時 Krylov 基底」を返し、その基底ベクトルを
**規格化して**返していた。実測で以下が判明したため全面改訂する。

 (1) v1 は方向しか返さず、**ゲージ（目盛）を規格化で捨てていた**。
     `v1 = x/‖x‖`・`d_plane /= ‖·‖`・`axis_vec /= ‖·‖`。
     規格化はスケール軸 R の消去であり、重力の消去である（規約 R10）。
 (2) v1 の返り値からは **z 成分を作れない**。法線を場の空間に戻していない。
 (3) v1 は**フレームを1本しか返さない**。位置がフレームの中にしか存在せず、
     位置の場にならない。局所フレームは M 本（関係波ごと）立つ。
 (4) 局所フレームだけでは軸の向きが波ごとにバラバラで、位置にならない。
     **大本（系全体）の基準三つ組に射影して初めて共通座標になる**
     ——その射影（方向余弦）が「方向のゲージ」である。

--------------------------------------------------------------------------
設計規約（R1–R7 を継承し、以下を改訂・追加）

  **R8'（撤回と差替え）**: v1 の「瞬時性（窓を使わない）」は撤回する。
  D は時間平均も窓も使わないが、それは規約ではなく構成の帰結である
  （二時刻メンバーの履歴は呼び出し側が `_carry` で持ち回す）。

  **R9 停止条件の外部化**（継承）: 閾値・ランク判定・打切りを持たない。
  計算不能なら NaN をそのまま返す。0 や空の代用値に置き換えない。

  **R10 スケール保存（新設）**: 読出しは比だけを返してはならない。比を返す
  ときは必ず分母（スケール）を同時に返す。**規格化は R 軸の消去であり、
  重力の消去である。** 単位方向ベクトルを返す場合は、その大きさ（ゲージ）を
  必ず別に返すこと。

--------------------------------------------------------------------------
構造（本モジュールの中核）

**局所場の単位は関係波 e**（G が返す粒度と同一）。波 e の局所場は、エンジン
自身の二チャネル分割による (a_e, b_e)＝(奇数帯＝フェルミオン型, 偶数帯＝
ボゾン型) の内容である。相互作用関数 F が生成子を作るときに使う分割そのもの
であり、新たな任意構造を持ち込まない。

各局所場について:

    代表複素数   z = Σ conj(a)·b
    ゲージ       T = Σ(|a|² + |b|²)         ← 等方スカラー1個・R′² = T²
    Stokes 3成分 X = Σ(|a|² − |b|²),  Y = 2 Re z,  Z = 2 Im z

    複素平面 = (Y, Z) 平面（代表複素数 z の実部・虚部が張る）
    直交軸   = X 方向

閉塞錐 x²+y²+z² = t²+R²+Q² に対し |r|² = T² であり、**|r| = T が等方ゲージ**。
（m² = T²−X²−Y²−Z² ≥ 0 は Cauchy–Schwarz より恒等的に成立＝光錐束縛。）

**大本の基準三つ組**は、系全体（全波の合算）の代表複素数の位相 φ_glob が作る:

    ê_X = (0,  cos φ_glob,  sin φ_glob)     ← 位相 0
    ê_Y = (0, −sin φ_glob,  cos φ_glob)     ← 位相 90°
    ê_Z = (1,  0,           0)              ← xy に直交

**方向のゲージ** = 波ごとの Stokes 単位ベクトルを、この三つ組へ射影した
方向余弦（3成分）。これがないと局所軸がバラバラのままで位置にならない。

**移動量**は位相の前進そのものである（θ = θ₀ + ωt ⇒ X = ωt）。したがって
本モジュールは一段の相対位相前進 δ = Δθ_local − Δφ_glob を返す。静止（時計と
同期）なら δ=0、光なら δ が最大。**位置はこれを積分して得る**（積分は G が
前ステップの位置を受け取って行う。原点は不明なので初期値 0・差だけが物理）。
"""
from __future__ import annotations
import numpy as np

_NAN = float("nan")


# ---------------------------------------------------------------- 子閉塞ブロック

def d_closure_blocks(C2):
    """ゼロ閉塞する子ブロックの検出（帳簿セル粒度）。v1 から無改訂。

    cell_power    : セルのパワー Σ_e|z_e|²
    cell_closure  : 閉塞残差 |Σ_e z_e²| / Σ_e|z_e|²（パワー0 のセルは NaN）
    total_closure : 全体の閉塞残差（＝物質の指標）

    凝縮体はここで閉塞残差が機械精度に落ちるセルとして現れる。
    どこで切るかは D の責務ではない（R9）——残差そのものを返す。
    """
    A2 = np.abs(C2) ** 2
    P = A2.sum(axis=0)
    B = np.abs(np.sum(C2 ** 2, axis=0))
    with np.errstate(divide="ignore", invalid="ignore"):
        R = B / P
    zf = C2.reshape(-1)
    pt = np.float64(np.real(np.vdot(zf, zf)))
    with np.errstate(divide="ignore", invalid="ignore"):
        tot = np.float64(abs(complex(np.sum(zf * zf)))) / pt
    return {"cell_power": P, "cell_closure": R, "total_closure": float(tot)}


# ---------------------------------------------------------------- 局所場の切り出し

def d_split_channels(C2):
    """二チャネル分割（エンジン自身の帯の偶奇対）。局所場の単位＝関係波 e。

    a : 奇数帯（フェルミオン型）内容  (M, J*Nη)
    b : 偶数帯（ボゾン型）内容        (M, J*Nη)
    Nn が奇数のとき対にならない最上帯は除外する（G と同一の規約）。
    """
    M, Nn, Neta = C2.shape
    J = Nn // 2
    ke = 2 * np.arange(J)
    ko = ke + 1
    a = C2[:, ko, :].reshape(M, -1)
    b = C2[:, ke, :].reshape(M, -1)
    return a, b


# ---------------------------------------------------------------- 局所フレーム

def d_local_frame(a, b):
    """局所場（二チャネル）から、代表複素数・等方ゲージ・Stokes 3成分を返す。

    軸方向は最後の軸に沿って縮約する（波ごとに1組を返す）。

    返す束:
      z        : 代表複素数 z = Σ conj(a)·b          （複素・スケールを持つ）
      gauge    : 等方ゲージ |r| = T = Σ(|a|²+|b|²)   （スカラー・**規格化しない**）
      stokes   : (…, 3) の実ベクトル (X, Y, Z)
                 X = Σ(|a|²−|b|²), Y = 2 Re z, Z = 2 Im z
      s_norm   : ‖(X,Y,Z)‖                            （スケールを持つ）
      s_hat    : 単位方向 (X,Y,Z)/‖·‖                  （‖·‖=0 なら NaN・代用しない）
      m2       : T² − X² − Y² − Z²（質量²・非コヒーレンス）
      theta    : arg(z)（複素平面内の位相・円上の値。巻数の展開は S の責務）
    """
    a = np.asarray(a); b = np.asarray(b)
    pa = np.sum(np.abs(a) ** 2, axis=-1)
    pb = np.sum(np.abs(b) ** 2, axis=-1)
    z = np.sum(np.conj(a) * b, axis=-1)
    T = pa + pb
    X = pa - pb
    Y = 2.0 * np.real(z)
    Z = 2.0 * np.imag(z)
    s = np.stack([X, Y, Z], axis=-1)
    s_norm = np.sqrt(np.sum(s ** 2, axis=-1))
    with np.errstate(divide="ignore", invalid="ignore"):
        s_hat = s / s_norm[..., None]          # 0/0 → NaN（代用しない・R9）
    m2 = T ** 2 - X ** 2 - Y ** 2 - Z ** 2
    with np.errstate(divide="ignore", invalid="ignore"):
        theta = np.angle(z)
    return {"z": z, "gauge": T, "stokes": s, "s_norm": s_norm, "s_hat": s_hat,
            "m2": m2, "theta": theta}


# ---------------------------------------------------------------- 大本の基準三つ組

def d_reference_triad(z_glob):
    """大本（系全体）の代表複素数の位相から、基準となる正規直交三つ組を作る。

    位相 0 を X 軸、位相 90° を Y 軸、複素平面に直交する方向を Z 軸とする:

        ê_X = (0,  cos φ,  sin φ)
        ê_Y = (0, −sin φ,  cos φ)
        ê_Z = (1,  0,      0)

    Stokes 空間の (Y, Z) 平面が代表複素数の住む複素平面、X が直交軸である。
    φ が読めない（z=0）場合は NaN を伝播させる（代用しない・R9）。

    返す束:
      triad   : (3, 3) 行 = ê_X, ê_Y, ê_Z
      phi     : φ_glob（円上の値）
    """
    with np.errstate(divide="ignore", invalid="ignore"):
        phi = float(np.angle(z_glob))
    if not np.isfinite(phi) or z_glob == 0:
        return {"triad": np.full((3, 3), _NAN), "phi": _NAN}
    c, s = np.cos(phi), np.sin(phi)
    triad = np.array([[0.0, c, s],
                      [0.0, -s, c],
                      [1.0, 0.0, 0.0]])
    return {"triad": triad, "phi": phi}


def d_direction_gauge(s_hat, triad):
    """方向のゲージ＝局所方向を大本の基準三つ組へ射影した方向余弦（3成分）。

    s_hat : (…, 3) 局所の単位方向
    triad : (3, 3) 大本の基準三つ組（行が ê_X, ê_Y, ê_Z）
    返す  : (…, 3) 方向余弦 (cX, cY, cZ)

    これがないと M 本の局所軸が各々勝手な向きを向いたままで、座標にならない。
    """
    return np.asarray(s_hat) @ np.asarray(triad).T


# ---------------------------------------------------------------- 移動量（位相前進）

def d_phase_advance(z, z_prev):
    """一段の位相前進 Δθ = arg(z · conj(z_prev))。円上の値（展開は S の責務）。

    z_prev が None（初回）なら NaN を返す（0 で代用しない・R9）。
    """
    if z_prev is None:
        return np.full(np.shape(z), _NAN)
    w = np.asarray(z) * np.conj(np.asarray(z_prev))
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.angle(w)


def d_displacement(theta_adv_local, phi_adv_glob):
    """移動量 δ ＝ 局所の位相前進 − 大本の位相前進（＝時計との離調）。

    θ = θ₀ + ωt ⇒ 移動量は位相の前進そのもの。
    静止（時計と同期）なら δ = 0、光なら δ が最大。
    """
    return np.asarray(theta_adv_local) - phi_adv_glob


# ---------------------------------------------------------------- 時計位相（参照系）

def d_clock_ref(Z, p2, q2):
    """集団時計の位相基準 φ（親平面の回転位相）＝復調の位相原点。v1 から無改訂。
    phi_weight が 0 のとき phi は意味を持たない（判定は呼び出し側・R9）。"""
    c = complex(p2 @ Z) + 1j * complex(q2 @ Z)
    return {"phi": float(np.angle(c)), "phi_weight": float(abs(c))}


# ---------------------------------------------------------------- 常時実行パネル

def d_panel(C2, carry=None):
    """常時実行パネル（D）。第0步から毎ステップ・一様ケイデンスで呼ぶ。

    carry : 前ステップの持ち回し {"z_local": …, "z_glob": …}。初回は None。
            **D は無状態**であり、履歴は呼び出し側が持ち回す（R5/R9）。

    返す束（すべて関係波ごと (M,) または (M,3)。大本は接尾辞 _glob）:
      local_z / local_gauge / local_stokes / local_s_hat / local_m2 / local_theta
      glob_z  / glob_gauge  / glob_stokes  / glob_s_hat  / glob_m2  / glob_theta
      triad          : (3,3) 大本の基準三つ組（位相0・位相90°・直交）
      phi_glob       : 大本の代表複素数の位相
      dir_gauge      : (M,3) 方向のゲージ＝方向余弦
      theta_adv      : (M,) 局所の一段位相前進（初回 NaN）
      phi_adv_glob   : 大本の一段位相前進（初回 NaN）
      displacement   : (M,) 移動量 δ（初回 NaN）。位置はこれを積分して得る
      _carry         : 次ステップへ渡す持ち回し
    """
    a, b = d_split_channels(C2)
    loc = d_local_frame(a, b)
    # 大本＝系全体（全波を1つの局所場として合算）
    glob = d_local_frame(a.reshape(1, -1), b.reshape(1, -1))
    z_glob = complex(np.asarray(glob["z"]).reshape(-1)[0])
    ref = d_reference_triad(z_glob)

    z_prev = carry.get("z_local") if carry else None
    zg_prev = carry.get("z_glob") if carry else None
    th_adv = d_phase_advance(loc["z"], z_prev)
    ph_adv = d_phase_advance(np.array([z_glob]), None if zg_prev is None
                             else np.array([zg_prev]))
    ph_adv = float(np.asarray(ph_adv).reshape(-1)[0])

    out = {
        "local_z": loc["z"], "local_gauge": loc["gauge"],
        "local_stokes": loc["stokes"], "local_s_norm": loc["s_norm"],
        "local_s_hat": loc["s_hat"], "local_m2": loc["m2"],
        "local_theta": loc["theta"],
        "glob_z": z_glob,
        "glob_gauge": float(np.asarray(glob["gauge"]).reshape(-1)[0]),
        "glob_stokes": np.asarray(glob["stokes"]).reshape(3),
        "glob_s_hat": np.asarray(glob["s_hat"]).reshape(3),
        "glob_m2": float(np.asarray(glob["m2"]).reshape(-1)[0]),
        "glob_theta": float(np.asarray(glob["theta"]).reshape(-1)[0]),
        "triad": ref["triad"], "phi_glob": ref["phi"],
        "dir_gauge": d_direction_gauge(loc["s_hat"], ref["triad"]),
        "theta_adv": th_adv, "phi_adv_glob": ph_adv,
        "displacement": d_displacement(th_adv, ph_adv),
        "_carry": {"z_local": np.asarray(loc["z"]).copy(), "z_glob": z_glob},
    }
    out.update(d_closure_blocks(C2))
    return out
