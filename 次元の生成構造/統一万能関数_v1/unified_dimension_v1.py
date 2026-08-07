#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""統一万能次元読出し関数 D v1 — フレーム（凝縮体）を読む

階層:  F（力学）→ **D（次元・フレーム読出し）** → G（量の読出し）→ S（選択層）
       確定値 = S ∘ G(·, D(·))

なぜ D が要るか: 座標は軸と原点と位相基準なしには読めない。G がこれまで
破綻せずに済んだのは、読んでいたのがパワー・比・位相差というフレーム非依存量
だけだったからである。xyz を読む段になって依存が露出した。フレームを供給する
のは凝縮体自身であり（原点＝親平面、軸＝回転平面と回転軸、位相基準＝集団時計）、
背景座標も絶対原点も要らない。

**凝縮体の定義（代数的）**: Σxₙ²=0 の内部で ΣAₙ²=0 を満たす子ブロック。
実測（2026-08-08）: 真空ポンプは |zᵀz|/‖z‖²=4.4e-15 で厳密ゼロ閉塞、物質シードは
0.318 で非ゼロ閉塞、全体の閉塞欠損は物質量に厳密比例（相対1.4e-10）。
すなわち**凝縮体＝ゼロ閉塞ブロック、物質＝閉塞の欠損**である。

--------------------------------------------------------------------------
設計規約: unified_readout_v2 の R1–R7 を継承し、以下を追加する。

  **R8 瞬時性**: D は窓（時間平均・移動SVD）を使わない。フレームは毎ステップ、
  その瞬間の状態と生成子だけから構成する。窓幅は調整パラメータであり R2
  （パラメータフリー）を破るため。生成子は正本 LowRankSystem.kmatvec を
  read-only で用い、状態には書き込まない（受動）。

  **R9 停止条件の外部化**: D は閾値・ランク判定・打切り・停止条件を一切
  持たない。計算が不能な場合は **NaN / None をそのまま返す**（0 や空の
  代用値に置き換えない・IF で先回りしない）。ランクや打切り次数のような
  「どこで切るか」は上位のループ（実験スクリプト）が宣言して行う。
  そのために D は**連続量の素材**を返す——逐次直交化の残差ノルム列、
  射影生成子そのもの、平面ごとの重みの列。

  【是正記録 2026-08-08】初版の D は (a) Gram–Schmidt に 1e-300 の閾値を
  埋め、(b) 値が取れないとき 0 の代用値を返す分岐を持ち、(c) Krylov 反復の
  停止条件を内部に持ち、(d) その閾値の産物である rank で分岐までしていた。
  これは統一G v1 で犯した「読出しに選択を焼き込む」失敗の再演であり、
  どの閾値の下での次元かを宣言できないため基準値が再現不能になる。全廃した。

  唯一の例外: 線形代数ルーチン（固有分解）は NaN 入力で例外を送出するため、
  例外を捕捉して NaN を返す。これは閾値ではなく「計算不能→NaN」の写像である。
--------------------------------------------------------------------------

瞬時フレームの構成: 実表現 x=[Re Z⊥; Im Z⊥] ∈ ℝ^{2M} に対し、生成子の作用
A（＝K を実部・虚部に同一に作用させたもの）で
  v₁ = x（位置）, v₂ = A x（速度）, v₃ = A² x（加速度）
を作り、逐次直交化して占有3次元部分空間 E を張る。**位置と速度の張る面が
回転平面、加速度が第3方向**を与える。生成子をこの3次元へ射影した反対称
Ω = EᵀAE の回転軸 â と、回転平面の法線（E 基底の第3軸）の一致度
|n̂·â| = |Ω₁₂|/‖ω‖ が第3次元の確定度である。

**次元の結晶化の読出しは平面の縮退度 n_eff（実効平面数）である**（単体テスト
2026-08-08 の実測: N=4 は n_eff=1.938 で2枚が拮抗＝一意な平面が選べず第3軸が
定義できない／N≥5 は 1.09–1.47 で1枚が卓越）。瞬時 align は N=4 の非結晶化を
検出できず（0.606）、二時刻の持続性も全 N で 1.0000 となり判別しない——
いずれも実測で確認した限界であり、記録として残す。
"""
from __future__ import annotations
import numpy as np

_NAN = float("nan")


# ---------------------------------------------------------------- 子閉塞ブロック

def d_closure_blocks(C2):
    """ゼロ閉塞する子ブロックの検出（帳簿セル粒度）。

    cell_power    : セルのパワー Σ_e|z_e|²
    cell_closure  : 閉塞残差 |Σ_e z_e²| / Σ_e|z_e|²
                    （パワー 0 のセルは 0/0 → **NaN**。0 で代用しない）
    total_closure : 全体の閉塞残差（＝物質の指標。0 なら物質を生まない真空）

    凝縮体はここで閉塞残差が機械精度に落ちるセルとして現れる。**どこで切るかは
    D の責務ではない**（R9）——残差そのものを返す。
    """
    A2 = np.abs(C2) ** 2
    P = A2.sum(axis=0)                                   # Nn×Nη
    B = np.abs(np.sum(C2 ** 2, axis=0))                  # |Σ_e z_e²|
    with np.errstate(divide="ignore", invalid="ignore"):
        R = B / P
    zf = C2.reshape(-1)
    # Python スカラの割り算は 0/0 で例外になるため numpy 型で割る（NaN を返す・R9）
    pt = np.float64(np.real(np.vdot(zf, zf)))
    with np.errstate(divide="ignore", invalid="ignore"):
        tot = np.float64(abs(complex(np.sum(zf * zf)))) / pt
    return {"cell_power": P, "cell_closure": R, "total_closure": float(tot)}


# ---------------------------------------------------------------- 瞬時フレーム

def _to_real(z):
    return np.concatenate([z.real, z.imag])


def _from_real(x):
    h = len(x) // 2
    return x[:h] + 1j * x[h:]


def _apply_A(x, kmatvec):
    """実表現への生成子の作用（K を実部・虚部へ同一に作用させる）"""
    return _to_real(kmatvec(_from_real(x)))


def _orthonormalize(vs):
    """逐次直交化。**ベクトルを一本も捨てない**（ランク判定をしない・R9）。
    戻り: E（列が正規直交ベクトル・退化列は NaN）、resid（各段の残差ノルム）。
    残差ノルムが「どこまで独立か」の連続量であり、ランクは呼び出し側が決める。"""
    E, resid = [], []
    for v in vs:
        w = np.asarray(v, float).copy()
        for e in E:
            with np.errstate(invalid="ignore"):
                w = w - (e @ w) * e
        r = float(np.linalg.norm(w))
        resid.append(r)
        with np.errstate(divide="ignore", invalid="ignore"):
            E.append(w / r)                       # r=0 → NaN（代用しない）
    return np.array(E).T, np.array(resid)


def d_frame(Z, kmatvec, p2=None, q2=None):
    """瞬時フレーム（位置・速度・加速度が張る占有3次元から構成）。

    Z      : 関係波ベクトル（複素 M）
    kmatvec: 正本 LowRankSystem.kmatvec（生成子の作用・read-only）
    p2, q2 : 親平面基底（与えれば Z⊥ を取る＝誕生時の参照系を除く）

    返す束（計算不能な成分は NaN・**ランク判定は行わない**）:
      weight    : ‖Z⊥‖²（フレームの存在重み。0 なら不在）
      resid     : 逐次直交化の残差ノルム3本（位置・速度・加速度の独立性の連続量。
                  **ランクはこの列から呼び出し側が決める**）
      d_plane   : 回転平面の複素方向（x,y 復調の軸）
      axis_vec  : 回転軸の複素方向
      Omega     : 射影生成子 Ω = EᵀAE（3×3 反対称・生の素材）
      omega_gen : ‖ω‖（回転の大きさ）
      align     : |n̂·â| = |Ω₁₂|/‖ω‖（第3次元の確定度・補助量）
    """
    z = np.asarray(Z, complex)
    if p2 is not None and q2 is not None:
        z = z - p2 * (p2 @ z) - q2 * (q2 @ z)
    w = float(np.real(np.vdot(z, z)))
    x = _to_real(z)
    with np.errstate(divide="ignore", invalid="ignore"):
        v1 = x / np.linalg.norm(x)
    v2 = _apply_A(v1, kmatvec)
    v3 = _apply_A(v2, kmatvec)
    E, resid = _orthonormalize([v1, v2, v3])
    try:
        AE = np.stack([_apply_A(E[:, j], kmatvec) for j in range(E.shape[1])],
                      axis=1)
        Om = E.T @ AE
        Om = 0.5 * (Om - Om.T)
        axis = np.array([Om[2, 1], Om[0, 2], Om[1, 0]])
        om = float(np.linalg.norm(axis))
        with np.errstate(divide="ignore", invalid="ignore"):
            ah = axis / om
        align = float(abs(ah[2]))
        d_plane = _from_real(E[:, 0]) + 1j * _from_real(E[:, 1])
        with np.errstate(divide="ignore", invalid="ignore"):
            d_plane = d_plane / np.linalg.norm(d_plane)
            axis_vec = _from_real(E @ ah)
            axis_vec = axis_vec / np.linalg.norm(axis_vec)
    except Exception:                       # 計算不能 → NaN（R9・閾値ではない）
        M = len(z)
        Om = np.full((3, 3), _NAN)
        om, align = _NAN, _NAN
        d_plane = np.full(M, _NAN, complex)
        axis_vec = np.full(M, _NAN, complex)
    return {"weight": w, "resid": resid, "d_plane": d_plane,
            "axis_vec": axis_vec, "Omega": Om, "omega_gen": om, "align": align}


def d_plane_ladder(Z, kmatvec, order, p2=None, q2=None):
    """回転平面の梯子（Krylov 部分空間上で射影生成子を対角化して得る）。

    **order は呼び出し側が宣言する打切り次数**（＝求める平面の枚数。D は
    停止条件を持たない・R9）。実験仕様に値を明記すること。

    返す束（計算不能なら NaN）:
      krylov_resid : 逐次直交化の残差ノルム列（長さ 2·order。**どこまで
                     独立かの連続量**——数値ランクは呼び出し側が決める）
      freqs        : 各平面の回転周波数（cycles/step・重み降順）
      weights      : 各平面への状態の射影パワー
      n_eff        : 実効平面数 ＝ 重みの参加比 PR（閾値なし・連続量）
                     **次元の結晶化の読出しはこの量である**
      Omega        : 射影生成子（生の素材）
    """
    z = np.asarray(Z, complex)
    if p2 is not None and q2 is not None:
        z = z - p2 * (p2 @ z) - q2 * (q2 @ z)
    x = _to_real(z)
    vs, cur = [], x
    for _ in range(2 * int(order)):
        vs.append(cur)
        cur = _apply_A(cur, kmatvec)
    E, kres = _orthonormalize(vs)
    d = E.shape[1]
    try:
        AE = np.stack([_apply_A(E[:, j], kmatvec) for j in range(d)], axis=1)
        Om = E.T @ AE
        Om = 0.5 * (Om - Om.T)
        ev, V = np.linalg.eig(Om)
        coef = V.conj().T @ (E.T @ x)
        im = np.imag(ev)
        sel = im > 0                       # 共役対の正の側（構造の選択・閾値でない）
        freqs = im[sel] / (2 * np.pi)
        wts = np.abs(coef[sel]) ** 2
        idx = np.argsort(-wts)
        freqs, wts = freqs[idx], wts[idx]
        with np.errstate(divide="ignore", invalid="ignore"):
            n_eff = float(wts.sum() ** 2 / (wts ** 2).sum())
    except Exception:                      # 計算不能 → NaN（R9）
        Om = np.full((d, d), _NAN)
        freqs = np.full(0, _NAN)
        wts = np.full(0, _NAN)
        n_eff = _NAN
    return {"krylov_resid": kres, "freqs": freqs, "weights": wts,
            "n_eff": n_eff, "Omega": Om}


def d_frame_persistence(frame, frame_prev):
    """フレームの持続性（二時刻）。

    【実測による限界の記録 2026-08-08】次元の結晶化の判別には**使えない**——
    N=4（非結晶化）でも N≥5 でも一様に 1.0000 を返した。結晶化の読出しは
    d_plane_ladder の n_eff（平面の縮退度）である。本メンバーは記録として残す。
    前時刻が無い場合は NaN（0 で代用しない・R9）。
    """
    if frame_prev is None:
        return {"axis_persist": _NAN, "plane_persist": _NAN}
    with np.errstate(divide="ignore", invalid="ignore"):
        a0, a1 = frame_prev["axis_vec"], frame["axis_vec"]
        d0, d1 = frame_prev["d_plane"], frame["d_plane"]
        ap = abs(np.vdot(a0, a1)) / (np.linalg.norm(a0) * np.linalg.norm(a1))
        pp = abs(np.vdot(d0, d1)) / (np.linalg.norm(d0) * np.linalg.norm(d1))
    return {"axis_persist": float(ap), "plane_persist": float(pp)}


def d_clock_ref(Z, p2, q2):
    """集団時計の位相基準 φ（親平面の回転位相）＝復調の位相原点。
    量子光学のホモダイン検波の局部発振器に当たるものを、宇宙自身が供給する。
    phi_weight が 0 のとき phi は意味を持たない（判定は呼び出し側・R9）。"""
    c = complex(p2 @ Z) + 1j * complex(q2 @ Z)
    return {"phi": float(np.angle(c)), "phi_weight": float(abs(c))}


def d_gauge(frame):
    """局所ゲージ（四脚場）: フレームが各関係波に割り当てる方向の成分と、
    その非一様度。一様な真空では滑らかに揃い、質量が偏在すると割り当てが
    ばらつく——**この非一様度が重力（ゲージの目盛の不等間隔）**である。

    col_norm : 各辺 e における局所軸の大きさ |d_plane[e]|²+|axis_vec[e]|²
    nonunif  : 非一様度 ＝ col_norm の変動係数（std/mean・0 割は NaN）
    """
    col = np.abs(frame["d_plane"]) ** 2 + np.abs(frame["axis_vec"]) ** 2
    with np.errstate(divide="ignore", invalid="ignore"):
        nonunif = float(col.std() / col.mean())
    return {"col_norm": col, "nonunif": nonunif}


# ---------------------------------------------------------------- パネル

def d_panel(C2, kmatvec, p2, q2, cell, order, frame_prev=None):
    """常時実行パネル（D）。第0步から毎ステップ・一様ケイデンスで呼ぶ。

    cell  : どのセルの内容からフレームを作るかの**宣言**（選択。実験仕様に明記）
    order : Krylov 打切り次数の**宣言**（停止条件は呼び出し側の責務・R9）
    frame_prev : 前ステップの frame（持続性メンバー用・呼び出し側が持ち回す）

    図で用いる量はすべて本パネルの出力から取る（実験側で独自計算をしない）。
    """
    Z = C2[:, cell[0], cell[1]]
    out = {}
    out.update(d_closure_blocks(C2))
    fr = d_frame(Z, kmatvec, p2, q2)
    out.update({f"frame_{k}": v for k, v in fr.items()})
    out.update({f"ladder_{k}": v for k, v in
                d_plane_ladder(Z, kmatvec, order, p2, q2).items()})
    out.update({f"clock_{k}": v for k, v in d_clock_ref(Z, p2, q2).items()})
    out.update({f"gauge_{k}": v for k, v in d_gauge(fr).items()})
    out.update({f"pers_{k}": v for k, v in
                d_frame_persistence(fr, frame_prev).items()})
    out["_frame"] = fr
    return out
