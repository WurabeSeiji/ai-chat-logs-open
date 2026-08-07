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

設計規約: R1–R7（unified_readout_v2 参照）を継承し、

  **R8 瞬時性**: D は窓（時間平均・移動SVD）を使わない。フレームは毎ステップ、
  その瞬間の状態と生成子だけから構成する。窓幅は調整パラメータであり R2
  （パラメータフリー）を破るため。生成子は正本 LowRankSystem.kmatvec を
  read-only で用い、状態には書き込まない（受動）。

瞬時フレームの構成: 実表現 x=[Re Z⊥; Im Z⊥] ∈ ℝ^{2M} に対し、生成子の作用
A（＝K を実部・虚部に同一に作用させたもの）で
  v₁ = x（位置）, v₂ = A x（速度）, v₃ = A² x（加速度）
を作り、Gram–Schmidt で占有3次元部分空間 E を張る。**位置と速度の張る面が
回転平面、加速度が第3方向**を与える。生成子をこの3次元へ射影した反対称
Ω = EᵀAE の回転軸 â と、回転平面の法線 n̂（＝E基底の第3軸）の一致度
|n̂·â| = |Ω₁₂|/‖ω‖ が**第3次元がどれだけ確定しているか**の読出しになる。

平面の梯子は Krylov 部分空間 {x, Ax, …, A^{2p−1}x} 上で A を対角化して得る。
p は**計算上の打切り次数**（物理閾値ではない）であり、値を記録し、p を上げても
結果が安定であることを審査で確認する。
"""
from __future__ import annotations
import numpy as np


# ---------------------------------------------------------------- 子閉塞ブロック

def d_closure_blocks(C2):
    """ゼロ閉塞する子ブロックの探索（帳簿セル粒度）。

    cell_power    : セルのパワー
    cell_closure  : 閉塞残差 |Σ_e z_e²| / Σ_e|z_e|²（0 なら厳密なゼロ閉塞ブロック）
    total_closure : 全体の閉塞残差（＝物質の指標。0 なら物質を生まない真空）

    凝縮体はここで閉塞残差が機械精度に落ちるセルとして現れる。閾値判定は
    しない（残差そのものを返す——どこで切るかは選択層の責務）。
    """
    M, Nn, Neta = C2.shape
    P = np.zeros((Nn, Neta))
    R = np.zeros((Nn, Neta))
    for k in range(Nn):
        for e in range(Neta):
            z = C2[:, k, e]
            p = float(np.real(np.vdot(z, z)))
            P[k, e] = p
            R[k, e] = abs(complex(np.sum(z * z))) / p if p > 0 else 0.0
    zf = C2.reshape(-1)
    pt = float(np.real(np.vdot(zf, zf)))
    return {"cell_power": P, "cell_closure": R,
            "total_closure": abs(complex(np.sum(zf * zf))) / pt if pt > 0 else 0.0}


# ---------------------------------------------------------------- 瞬時フレーム

def _to_real(z):
    return np.concatenate([z.real, z.imag])


def _from_real(x):
    h = len(x) // 2
    return x[:h] + 1j * x[h:]


def _apply_A(x, kmatvec):
    """実表現への生成子の作用（K を実部・虚部へ同一に作用させる）"""
    return _to_real(kmatvec(_from_real(x)))


def _gram_schmidt(vs, tol=0.0):
    """数値ランクを保った直交化（tol=0: 厳密に非零な成分のみ採用）"""
    E = []
    for v in vs:
        w = v.copy()
        for e in E:
            w = w - (e @ w) * e
        nrm = float(np.linalg.norm(w))
        if nrm > tol and nrm > 1e-300:
            E.append(w / nrm)
    return np.array(E).T if E else np.zeros((len(vs[0]), 0))


def d_frame(Z, kmatvec, p2=None, q2=None):
    """瞬時フレーム（位置・速度・加速度が張る占有3次元から構成）。

    Z      : 関係波ベクトル（複素 M）——凝縮体セルの内容を渡す
    kmatvec: 正本 LowRankSystem.kmatvec（生成子の作用・read-only）
    p2, q2 : 親平面基底（与えれば Z⊥ を取る＝誕生時の参照系を除く）

    返す束:
      d_plane   : 回転平面の複素方向（M）——x,y 復調の軸
      axis_vec  : 回転軸（実 2M 表現の 3 次元係数を M へ戻した複素方向）
      align     : |n̂·â| ＝ 第3次元の確定度（1 に近いほど鋭い・0 なら非結晶）
      omega_gen : 射影生成子の回転の大きさ ‖ω‖
      rank      : 占有部分空間の実次元（3 未満なら次元が立っていない）
      weight    : ‖Z⊥‖²（フレームの存在重み。0 なら不在——真偽値は返さない）
    """
    z = Z.astype(complex)
    if p2 is not None and q2 is not None:
        z = z - p2 * (p2 @ z) - q2 * (q2 @ z)
    w = float(np.real(np.vdot(z, z)))
    M = len(z)
    empty = {"d_plane": np.zeros(M, complex), "axis_vec": np.zeros(M, complex),
             "align": 0.0, "omega_gen": 0.0, "rank": 0, "weight": w}
    if w <= 0.0:
        return empty
    x = _to_real(z)
    v1 = x / np.linalg.norm(x)
    v2 = _apply_A(v1, kmatvec)
    v3 = _apply_A(v2, kmatvec)
    E = _gram_schmidt([v1, v2, v3])
    if E.shape[1] < 3:
        empty["rank"] = int(E.shape[1])
        return empty
    AE = np.stack([_apply_A(E[:, j], kmatvec) for j in range(3)], axis=1)
    Om = E.T @ AE
    Om = 0.5 * (Om - Om.T)                       # 反対称成分（射影生成子）
    axis = np.array([Om[2, 1], Om[0, 2], Om[1, 0]])
    om = float(np.linalg.norm(axis))
    if om <= 0.0:
        empty["rank"] = 3
        return empty
    ah = axis / om
    # 回転平面 = {位置, 速度} が張る面 → その法線は E 基底の第3軸 (0,0,1)
    align = float(abs(ah[2]))
    d_plane = _from_real(E[:, 0]) + 1j * _from_real(E[:, 1])
    d_plane = d_plane / np.linalg.norm(d_plane)
    axis_vec = _from_real(E @ ah)
    nv = np.linalg.norm(axis_vec)
    axis_vec = axis_vec / nv if nv > 0 else axis_vec
    return {"d_plane": d_plane, "axis_vec": axis_vec, "align": align,
            "omega_gen": om, "rank": 3, "weight": w}


def d_plane_ladder(Z, kmatvec, p2=None, q2=None, order=6):
    """回転平面の梯子（Krylov 部分空間上で生成子を対角化して得る）。

    order : 計算上の打切り次数（＝求める平面の枚数。物理閾値ではない）
    返す束:
      freqs      : 各平面の回転周波数（射影生成子の固有値の虚部・降順の重み順）
      weights    : 各平面への状態の射影パワー（重み）
      n_eff      : 実効平面数 ＝ 重みの参加比 PR（閾値なし・連続量）
      krylov_dim : 実際に張れた次元（数値ランク）
    """
    z = Z.astype(complex)
    if p2 is not None and q2 is not None:
        z = z - p2 * (p2 @ z) - q2 * (q2 @ z)
    if float(np.real(np.vdot(z, z))) <= 0.0:
        return {"freqs": np.zeros(0), "weights": np.zeros(0), "n_eff": 0.0,
                "krylov_dim": 0}
    x = _to_real(z)
    vs, cur = [], x / np.linalg.norm(x)
    for _ in range(2 * order):
        vs.append(cur)
        cur = _apply_A(cur, kmatvec)
        nn = np.linalg.norm(cur)
        if nn <= 1e-300:
            break
        cur = cur / nn
    E = _gram_schmidt(vs)
    d = E.shape[1]
    if d < 2:
        return {"freqs": np.zeros(0), "weights": np.zeros(0), "n_eff": 0.0,
                "krylov_dim": int(d)}
    AE = np.stack([_apply_A(E[:, j], kmatvec) for j in range(d)], axis=1)
    Om = E.T @ AE
    Om = 0.5 * (Om - Om.T)
    ev, V = np.linalg.eig(Om)                    # 反対称 → 純虚固有値の共役対
    coef = V.conj().T @ (E.T @ x)                # 各固有方向への射影
    im = np.imag(ev)
    pos = im > 0
    freqs, wts = [], []
    for j in np.where(pos)[0]:
        freqs.append(float(im[j] / (2 * np.pi)))     # cycles/step 換算
        wts.append(float(abs(coef[j]) ** 2))
    if not wts:
        return {"freqs": np.zeros(0), "weights": np.zeros(0), "n_eff": 0.0,
                "krylov_dim": int(d)}
    wts = np.array(wts)
    freqs = np.array(freqs)
    idx = np.argsort(-wts)
    wts, freqs = wts[idx], freqs[idx]
    s1, s2 = wts.sum(), (wts ** 2).sum()
    n_eff = float(s1 ** 2 / s2) if s2 > 0 else 0.0
    return {"freqs": freqs, "weights": wts, "n_eff": n_eff, "krylov_dim": int(d)}


def d_frame_persistence(frame, frame_prev):
    """フレームの持続性（二時刻メンバー）——**次元が結晶化しているかの読出し**。

    【設計の是正記録 2026-08-08】瞬時フレーム単独では次元の結晶化を測れない
    ことが資格審査 Q16 で判明した（N=4 は第3次元が結晶化しないはずだが、
    瞬時版は 0.527 を返した）。理由: {位置, 速度, 加速度} は一般に必ず3次元を
    張るので、瞬間のスナップショットは常に rank 3 を返す。**結晶化とは同じ
    方向が時間的に選ばれ続けること＝持続性**であり、瞬間の性質ではない。
    窓（移動平均・移動SVD）を導入せずにこれを測るため、連続する二時刻の
    フレームの重なりを返す（時計メンバーと同じ二時刻の作法・R8 を侵さない）。

    axis_persist  : |⟨â(τ−1), â(τ)⟩|（回転軸の持続・1 なら完全に同じ方向）
    plane_persist : |⟨d(τ−1), d(τ)⟩|（回転平面の持続）
    重みは frame["weight"]。前時刻が無い場合は 0（不在＝重み0で表す）。
    """
    if frame_prev is None or frame["rank"] < 3 or frame_prev.get("rank", 0) < 3:
        return {"axis_persist": 0.0, "plane_persist": 0.0}
    a0, a1 = frame_prev["axis_vec"], frame["axis_vec"]
    d0, d1 = frame_prev["d_plane"], frame["d_plane"]
    na = np.linalg.norm(a0) * np.linalg.norm(a1)
    nd = np.linalg.norm(d0) * np.linalg.norm(d1)
    return {"axis_persist": float(abs(np.vdot(a0, a1)) / na) if na > 0 else 0.0,
            "plane_persist": float(abs(np.vdot(d0, d1)) / nd) if nd > 0 else 0.0}


def d_clock_ref(Z, p2, q2):
    """集団時計の位相基準 φ（親平面の回転位相）＝復調の位相原点。
    量子光学のホモダイン検波の局部発振器に当たるものを、宇宙自身が供給する。"""
    a = complex(p2 @ Z)
    b = complex(q2 @ Z)
    c = a + 1j * b
    return {"phi": float(np.angle(c)) if abs(c) > 0 else float("nan"),
            "phi_weight": float(abs(c))}


def d_gauge(frame, kmatvec):
    """局所ゲージ（四脚場）: フレームが各関係波に割り当てる方向の成分と、
    その非一様度。一様な真空では滑らかに揃い、質量が偏在すると割り当てが
    ばらつく——**この非一様度が重力（ゲージの目盛の不等間隔）**である。

    col_norm : 各辺 e における局所軸の大きさ |d_plane[e]|²+|axis[e]|²
    nonunif  : 非一様度 ＝ col_norm の変動係数（std/mean・連続量）
    """
    d, a = frame["d_plane"], frame["axis_vec"]
    col = np.abs(d) ** 2 + np.abs(a) ** 2
    mu = float(col.mean())
    return {"col_norm": col,
            "nonunif": float(col.std() / mu) if mu > 0 else 0.0}


# ---------------------------------------------------------------- パネル

def d_panel(C2, kmatvec, p2, q2, cell=(2, 0), order=6):
    """常時実行パネル（D）。第0步から毎ステップ・一様ケイデンスで呼ぶ。
    cell は「どのセルの内容からフレームを作るか」の宣言（選択）であり、
    実験仕様に明記すること。既定は凝縮体セル（帯2・巻き0）。"""
    Z = C2[:, cell[0], cell[1]]
    out = {}
    out.update(d_closure_blocks(C2))
    fr = d_frame(Z, kmatvec, p2, q2)
    out.update({f"frame_{k}": v for k, v in fr.items()})
    out.update({f"ladder_{k}": v for k, v in
                d_plane_ladder(Z, kmatvec, p2, q2, order).items()})
    out.update({f"clock_{k}": v for k, v in d_clock_ref(Z, p2, q2).items()})
    out.update({f"gauge_{k}": v for k, v in d_gauge(fr, kmatvec).items()})
    out["_frame"] = fr
    return out
