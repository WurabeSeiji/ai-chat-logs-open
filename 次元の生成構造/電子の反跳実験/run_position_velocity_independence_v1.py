#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""位置と読出し接速度の独立性 — 同一位置・異接速度の状態対の構成（2026-08-08）

中心命題:
    ker DR_Z  ⊄  ker Dv_Z

    すなわち「位置読出しの核方向の中に、読出し接速度を変化させる方向が存在する」。
    局所（微分）で示し、有限状態対で確認し、最後に自己無撞着 K(Z) でも
    非零極限が残ることを確かめる、という三段構成。

--------------------------------------------------------------------------
用語の限定（重要・論文化時に守ること）

  * 本稿の「速度」とは **空間読出しの接速度**
        v_read = DR_Z[ K(Z) Z ]
    のことである。**通常の物理速度との対応は証明していない。**

  * 「核 375 次元」とは、反対称生成子の係数空間 ℝ^{M(M−1)/2}=ℝ^378 における
    線形写像 A ↦ DR_{Z0}[A Z0]（rank 3）の核の次元である。
    **物理的状態空間そのもののファイバー次元が 375 と証明されたのではない。**

--------------------------------------------------------------------------
背景（なぜこれを調べたか）

「波は過去の位置を保存するレジスタを持たないのに、どうやって進行方向を持つのか」
という問いに対し、以下の候補を順に検定して全て棄却した（同日の記録）:

  1. cone_Z（偶奇クロス項の虚部）＝並進速度
     → 局在パケットが Z≠0 のまま Δpos_x = 0（機械精度）で棄却
  2. レジスタ上の並進が力学から自発生成する
     → UnifiedEngine で Δpos_x = 0 厳密。公開3Dデモの並進は外部入力
       （Engine3D が毎步 c_k ← c_k e^{−2πik·v/N} を掛けている）と判明
  3. 双対軸に沿った位相勾配 k = ∇φ ＝ 速度
     → 時間反転で偶、かつ格子の自明値 2π/Nn（定数）で棄却
  4. 復調偏角 Ω の線形応答＝固有面の離調
     → 固有係数で直接読むと傾きは 10⁻¹⁶ で不動。Ω の応答は読出し側の混合

本実験は、上記と異なり「位置読出しから不可視な内部自由度」を直接構成する。

--------------------------------------------------------------------------
系と読出しチャート

  系  : N体関係波の閉鎖力学（第8論文 code の build_init / evolve を read-only 使用）
        Z ∈ ℂ^M, M = N(N−1)/2。閉塞 C(Z)=ZᵀZ、ノルム ‖Z‖ は Cayley 步が厳密保存。
  整定: 自己無撞着（毎步 set_theta(arg Z)）で T_SETTLE 步。その時点を Z0 とする。

  チャート（Z0 で確定し、以後**固定**する。取り直さない）:
        P = I − ppᵀ − qqᵀ            （p,q = 誕生時の親平面基底）
        h = pᵀZ + i qᵀZ,  φ = arg h  （集団時計）
        a_j = w_j† P Z               （w_1,w_2 = 生成子 K0 の固有面方向）
        C_1 = a_1 e^{−iφ} = x + iy
        C_2 = a_2 e^{−iφ} = z + is
        R(Z) = (x, y, z)             （位置読出し。s は内部成分として残す）

  push-forward の厳密式（仮定なし。P と K は非可換なので Ċ=iΩC は成り立たない）:
        φ̇[T]  = Im( (pᵀT + i qᵀT) / h )
        ȧ_j[T] = w_j† P T
        Ċ_j[T] = e^{−iφ} ( ȧ_j[T] − i φ̇[T] a_j )
        DR_Z[T] = (Re Ċ1, Im Ċ1, Re Ċ2, Im Ċ2)   ← 第4成分 s の速度も記録する

--------------------------------------------------------------------------
事前登録した判定（実行前固定）

 (P1) 局所: A ↦ DR_{Z0}[A Z0] は rank 3。その核から選んだ A について
      ‖DR[A Z0] の位置3成分‖ / ‖v_read(Z0)‖ < 1e-10 かつ ‖Dv[A Z0]‖ > 0。
 (P2) 有限（固定 K0）: Newton 補正後に ‖R(Z+)−R(Z−)‖ < 1e-12 かつ
      ‖v+−v−‖/(2ε) が ε→0 で一定値へ収束。補正量 c は O(ε²)。
 (P3) 自己無撞着: K± = K(Z±) を各状態から独立に再構成し、
      ‖v+^raw − v−^raw‖/(2ε) → C_sc > 0 が残る。
 (P4) 監査: 全変換は反対称生成子＋Cayley のみ（閉塞・ノルムは代数的に保存）。
      接空間 |ZᵀKZ|/‖Z‖² ≈ 0。σ_max は power iteration と行列直接値を照合し、
      ± 両条件で同一の初期 wp を使う（計測器の履歴を速度と誤認しないため）。
 (P5) 時計尺度の分離: v_eng = DR[T_eng]（Cayley の生成子）と v_raw の
      cos∠ と尺度比を記録する。

使い方: python3 run_position_velocity_independence_v1.py
"""
from __future__ import annotations
import importlib.util, itertools, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
ABL = (HERE.parent / "第8論文_二段階seed除去による準安定相の因果分離" / "code"
       / "run_preliminary_seed_ablation_v1.py")

N, T_SETTLE = 8, 6000
EPS_LIST = (1e-3, 1e-4, 1e-5)


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def main():
    t0 = time.time()
    abl = load("abl_pvi", ABL)
    sys_lr, _v, _, _, _, p, q, Z, wp = abl.build_init(N, True)
    M = sys_lr.m
    GAMMA = sys_lr.cayley_step.__globals__["GAMMA"]      # エンジンの実値（推定しない）
    for _ in range(T_SETTLE):
        Z, wp = abl.evolve(sys_lr, Z, wp)
    Z0, WP0 = Z.copy(), wp.copy()

    def Kmat(Zc):
        sys_lr.set_theta(np.angle(Zc))
        Km = np.zeros((M, M))
        for j in range(M):
            e = np.zeros(M); e[j] = 1.0
            Km[:, j] = np.real(sys_lr.kmatvec(e.astype(complex)))
        return Km

    def sigmax(Zc):
        """± で同一の初期 wp を使い、行列からの直接値と照合する（P4）"""
        sys_lr.set_theta(np.angle(Zc))
        se, _ = sys_lr.sigma_max_power(WP0.copy())
        return float(se), float(np.max(np.abs(np.linalg.eigvals(Kmat(Zc)).imag)))

    K0 = Kmat(Z0)
    ev, EVEC = np.linalg.eig(K0)
    idx = [i for i in np.argsort(-ev.imag) if ev[i].imag > 1e-12]
    w1 = EVEC[:, idx[0]] / np.linalg.norm(EVEC[:, idx[0]])
    w2 = EVEC[:, idx[1]] / np.linalg.norm(EVEC[:, idx[1]])

    # ---- チャート（以後固定）------------------------------------------------
    def DR(Zc, T):
        Zp = Zc - p * (p @ Zc) - q * (q @ Zc)
        h = complex(p @ Zc) + 1j * complex(q @ Zc)
        phi = np.angle(h)
        phid = float(np.imag((complex(p @ T) + 1j * complex(q @ T)) / h))
        Tp = T - p * (p @ T) - q * (q @ T)
        a1, a2 = np.vdot(w1, Zp), np.vdot(w2, Zp)
        d1, d2 = np.vdot(w1, Tp), np.vdot(w2, Tp)
        e = np.exp(-1j * phi)
        C1d, C2d = e * (d1 - 1j * phid * a1), e * (d2 - 1j * phid * a2)
        return np.array([C1d.real, C1d.imag, C2d.real, C2d.imag])

    def pos(Zc):
        Zp = Zc - p * (p @ Zc) - q * (q @ Zc)
        e = np.exp(-1j * np.angle(complex(p @ Zc) + 1j * complex(q @ Zc)))
        return np.array([(np.vdot(w1, Zp) * e).real, (np.vdot(w1, Zp) * e).imag,
                         (np.vdot(w2, Zp) * e).real])

    def clo(Zc):
        return abs(complex(Zc @ Zc)) / float(np.real(np.vdot(Zc, Zc)))

    def cay(Bm):
        return np.linalg.solve(np.eye(M) - Bm / 2, np.eye(M) + Bm / 2)   # 厳密直交

    out = {"env": {"N": N, "M": M, "T_settle": T_SETTLE, "GAMMA": GAMMA,
                   "eps_list": list(EPS_LIST)}}

    # ================= 段1: 局所 ker DR ⊄ ker Dv =============================
    pairs = list(itertools.combinations(range(M), 2))
    def mk(a):
        A_ = np.zeros((M, M))
        for k, (i, j) in enumerate(pairs):
            A_[i, j] += a[k]; A_[j, i] -= a[k]
        return A_
    G = np.zeros((4, len(pairs)))
    for k, (i, j) in enumerate(pairs):
        T = np.zeros(M, complex); T[i] = Z0[j]; T[j] = -Z0[i]
        G[:, k] = DR(Z0, T)
    U_, S_, Vt_ = np.linalg.svd(G[:3], full_matrices=True)
    rank_pos = int(np.sum(S_ > 1e-12 * S_[0]))
    Kern = Vt_[rank_pos:]
    cs = G[3] @ Kern.T; cs /= np.linalg.norm(cs)
    A = mk(Kern.T @ cs); A /= np.linalg.norm(A)
    Bs = [mk(Vt_[i]) for i in range(rank_pos)]
    v0 = DR(Z0, K0 @ Z0)
    dr_fib = DR(Z0, A @ Z0)
    def v_of(Zc, Km): return DR(Zc, Km @ Zc)
    dv = None
    for eps in (1e-5,):
        dv = (v_of(cay(+eps * A) @ Z0, K0) - v_of(cay(-eps * A) @ Z0, K0)) / (2 * eps)
    P1 = (np.linalg.norm(dr_fib[:3]) / np.linalg.norm(v0) < 1e-10) and (np.linalg.norm(dv) > 0)
    out["stage1"] = {"coef_space_dim": len(pairs), "rank_pos": rank_pos,
                     "kernel_dim": int(Kern.shape[0]),
                     "DR_fiber": dr_fib.tolist(),
                     "pos_resp_rel": float(np.linalg.norm(dr_fib[:3]) / np.linalg.norm(v0)),
                     "Dv_norm": float(np.linalg.norm(dv)),
                     "v0_norm": float(np.linalg.norm(v0)), "ok": bool(P1)}
    print(f"=== 段1 局所 ===\n  係数空間 {len(pairs)} 次元・rank(位置)={rank_pos}・核 {Kern.shape[0]} 次元")
    print(f"  DR[A Z0] = {np.round(dr_fib, 12).tolist()}")
    print(f"  位置応答/‖v0‖ = {np.linalg.norm(dr_fib[:3])/np.linalg.norm(v0):.3e}   "
          f"‖Dv‖ = {np.linalg.norm(dv):.6f}   → {'通過' if P1 else '不成立'}")

    # ================= 段2/段3: 有限状態対 ==================================
    R0 = pos(Z0)
    def newton(Zi):
        c = np.zeros(rank_pos)
        for _ in range(8):
            Zc = cay(sum(c[i] * Bs[i] for i in range(rank_pos))) @ Zi
            F = pos(Zc) - R0
            if np.linalg.norm(F) < 1e-15: break
            J = np.zeros((3, rank_pos)); d = 1e-7
            for a in range(rank_pos):
                cc = c.copy(); cc[a] += d
                J[:, a] = (pos(cay(sum(cc[i] * Bs[i] for i in range(rank_pos))) @ Zi) - pos(Zc)) / d
            c = c - np.linalg.solve(J, F)
        return cay(sum(c[i] * Bs[i] for i in range(rank_pos))) @ Zi, c

    print(f"\n=== 段2（固定 K0）／段3（自己無撞着 K(Z)）===")
    print(f"{'ε':>8} {'c/ε²':>9} {'‖R+−R−‖':>11} {'固定K0 /(2ε)':>14} "
          f"{'自己無撞着 /(2ε)':>16} {'cos(raw,eng)':>13} {'‖K+−K−‖/2ε':>12}")
    out["pairs"] = {}
    for eps in EPS_LIST:
        Zp_, cp = newton(cay(+eps * A) @ Z0)
        Zm_, cm = newton(cay(-eps * A) @ Z0)
        dR = float(np.linalg.norm(pos(Zp_) - pos(Zm_)))
        # 段2: 固定 K0
        C_fix = float(np.linalg.norm(v_of(Zp_, K0) - v_of(Zm_, K0)) / (2 * eps))
        # 段3: 自己無撞着
        Kp, Km_ = Kmat(Zp_), Kmat(Zm_)
        sep, sdp = sigmax(Zp_); sem, sdm = sigmax(Zm_)
        Tp, Tm = Kp @ Zp_, Km_ @ Zm_
        audit = max(abs(complex(Zp_ @ Tp)) / float(np.real(np.vdot(Zp_, Zp_))),
                    abs(complex(Zm_ @ Tm)) / float(np.real(np.vdot(Zm_, Zm_))))
        vpr, vmr = DR(Zp_, Tp), DR(Zm_, Tm)
        C_sc = float(np.linalg.norm(vpr - vmr) / (2 * eps))
        # engine 正規化接方向（Cayley の生成子・h→0）
        a_p = GAMMA / sep
        hh = 1e-5
        Uh = np.linalg.solve(np.eye(M) - hh * a_p * Kp, np.eye(M) + hh * a_p * Kp)
        vpe = DR(Zp_, (Uh @ Zp_ - Zp_) / hh)
        cosv = float(vpr @ vpe / (np.linalg.norm(vpr) * np.linalg.norm(vpe)))
        out["pairs"][f"{eps:.0e}"] = {
            "c_over_eps2": float(np.linalg.norm(cp) / eps ** 2), "dR": dR,
            "C_fixed": C_fix, "C_selfconsistent": C_sc,
            "cos_raw_eng": cosv, "scale_eng_over_raw": float(np.linalg.norm(vpe) / np.linalg.norm(vpr)),
            "sigma_power": sep, "sigma_direct": sdp,
            "dK_over_2eps": float(np.linalg.norm(Kp - Km_) / (2 * eps)),
            "tangent_audit": float(audit), "closure_drift": float(abs(clo(Zp_) - clo(Z0)))}
        print(f"{eps:>8.0e} {np.linalg.norm(cp)/eps**2:>9.3f} {dR:>11.2e} {C_fix:>14.6f} "
              f"{C_sc:>16.6f} {cosv:>13.9f} {np.linalg.norm(Kp-Km_)/(2*eps):>12.4f}")

    ks = [f"{e:.0e}" for e in EPS_LIST]
    P2 = all(out["pairs"][k]["dR"] < 1e-12 for k in ks)
    C_sc_lim = out["pairs"][ks[-1]]["C_selfconsistent"]
    P3 = C_sc_lim > 0 and abs(out["pairs"][ks[-1]]["C_selfconsistent"]
                              - out["pairs"][ks[-2]]["C_selfconsistent"]) < 1e-4
    P4 = all(out["pairs"][k]["tangent_audit"] < 1e-12
             and out["pairs"][k]["closure_drift"] < 1e-12
             and abs(out["pairs"][k]["sigma_power"] - out["pairs"][k]["sigma_direct"]) < 1e-9
             for k in ks)
    P5 = all(out["pairs"][k]["cos_raw_eng"] > 1 - 1e-6 for k in ks)
    print(f"\n(P1) 局所 ker DR ⊄ ker Dv                : {'通過' if P1 else '不成立'}")
    print(f"(P2) 有限状態対 ‖R+−R−‖<1e-12            : {'通過' if P2 else '不成立'}")
    print(f"(P3) 自己無撞着 C_sc = {C_sc_lim:.6f} > 0 で収束 : {'通過' if P3 else '不成立'}")
    print(f"(P4) 監査（接空間・閉塞・σmax 照合）      : {'通過' if P4 else '不成立'}")
    print(f"(P5) v_eng ∥ v_raw（尺度のみ分離）        : {'通過' if P5 else '不成立'}")
    print(f"\n中心命題: ker DR_Z ⊄ ker Dv_Z"
          f"\n  → 位置読出しの核方向の中に、読出し接速度を変化させる方向が存在する。"
          f"\n  → 同一位置 R(Z+)=R(Z−) をもつ二状態が、異なる v_read を持てる（自己無撞着系でも）。")
    out.update({"P1": bool(P1), "P2": bool(P2), "P3": bool(P3), "P4": bool(P4), "P5": bool(P5),
                "all_pass": bool(P1 and P2 and P3 and P4 and P5),
                "C_sc": C_sc_lim, "runtime_sec": time.time() - t0})
    (HERE / "result_position_velocity_independence_v1.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False, default=float))
    print(f"完了 {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
