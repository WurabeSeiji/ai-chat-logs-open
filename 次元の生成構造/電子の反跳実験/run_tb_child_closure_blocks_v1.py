#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""子閉塞ブロック検定 — 凝縮体は Σxₙ²=0 の内部にある「子のゼロ閉塞塊」か

木原の指摘（2026-08-08）:
  Σxₙ²=0 自身が凝縮体の存在を予言している。全体の閉塞の中に
  ΣAₙ²=0, ΣBₙ²=0 … となる塊があれば、その総和もまた 0 で閉塞条件を満たす。
  この A, B が凝縮体のはずで、各凝縮体の中で Σxₙ² を計算すると厳密に
  ゼロ閉塞するはずである。

これは凝縮体の定義を「特異値の対が立っている塊」という経験的兆候から、
**閉塞条件の代数的帰結**へ格上げする主張である。本実験はこれを検定する。

代数的背景（本実験が確かめる仮説）:
  閉塞 Σ_e z_e² = 0 は複素**双一次**形式 zᵀz（共役なし）である。状態を
  部分空間へ分解したとき、各成分が独立に自己零 (P z)ᵀ(P z)=0 であるためには、
  その部分空間が双一次形式について**完全等方 (totally isotropic)** である
  必要がある（2次元なら dᵀd=0 かつ d₁ᵀd₂=0）。すなわち
  「子のゼロ閉塞塊」＝「完全等方部分空間」であり、凝縮体の回転平面が
  それに一致するか否かが問われる。ℂ^M の最大等方部分空間の次元は ⌊M/2⌋。

事前登録した判定（実行前固定）:
 (B1) 閉塞の所在: 全体和 / セル別 / 双対点別 のどれが厳密ゼロかを判定・記録。
      （どの粒度で閉塞が成り立っているかを先に確定させる。判定でなく事実の確定）
 (B2) 平面の自己等方性: 各平面方向 d_k について |dᵀd|/(d†d) < 1e-12
 (B3) 平面間の等方性: |d_jᵀd_k|/(|d_j||d_k|) < 1e-12
 (B4) 子の閉塞: 状態を各平面へ射影した成分 Z_k について
      |Σ_e (Z_k)_e²| / Σ_e|Z_k|² < 1e-12
 (B5) 対照: ランダム部分集合・ランダム部分空間では上記が有意に非ゼロ（>1e-3）
 (B6) 個数: 子閉塞を満たす平面の数を数え、特異値対の梯子と対応づける
不成立の場合もそのまま記録する（近似ブロックなら残差の桁を報告する）。

使い方: python3 run_tb_child_closure_blocks_v1.py
"""
from __future__ import annotations
import importlib.util, json, sys, time
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams["font.family"] = ["Hiragino Sans", "sans-serif"]
plt.rcParams["axes.unicode_minus"] = False

HERE = Path(__file__).resolve().parent
UF = HERE.parent / "統一万能関数_v1"


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


ui = load("ui_cb", UF / "unified_interaction_v1.py")
G = load("ur_cb", UF / "unified_readout_v2.py")

N, NN, NETA = 12, 16, 8
T = 4000
DELTA = 1e-2
WIN = (2000, 4000)
EV = 5                    # サンプリング間引き（フレーム抽出用）
TOL = 1e-12
CTRL_MIN = 1e-3


def build(delta, seed=2):
    m = N * (N - 1) // 2
    _, v, _, _, _, _, _, Z0c, wp0 = ui.abl.build_init(N, False)
    r2 = ui.gen3.make_parent(N, seed=seed)
    Csec = np.fft.fft(r2.relation_waves, axis=1) / N
    s0 = Csec[:, 1] / np.linalg.norm(Csec[:, 1])
    C2_0 = np.zeros((m, NN, NETA), complex)
    C2_0[:, 2, 0] = Z0c
    if delta > 0:
        C2_0[:, 1, 0] = delta * s0
    p2 = C2_0[:, 2, 0].real / np.linalg.norm(C2_0[:, 2, 0].real)
    q2 = C2_0[:, 2, 0].imag - (C2_0[:, 2, 0].imag @ p2) * p2
    q2 = q2 / np.linalg.norm(q2)
    return ui.UnifiedEngine(N, C2_0, wp0), p2, q2


def bil(z):
    """双一次形式 zᵀz（共役なし）＝ゼロ閉塞の左辺"""
    return complex(np.sum(z * z))


def rel_bil(z):
    """規格化した閉塞残差 |zᵀz| / z†z"""
    d = float(np.real(np.vdot(z, z)))
    return abs(bil(z)) / d if d > 0 else 0.0


def main():
    t0 = time.time()
    eng, p2, q2 = build(DELTA)
    samples, hist = [], {"total": [], "cell_max": [], "dual_max": []}
    for t in range(T):
        eng.step()
        C2 = eng.C2()
        # --- B1: 閉塞の所在（監査量・読出しの外）
        hist["total"].append(rel_bil(C2.reshape(-1)))
        cell = [rel_bil(C2[:, k, e]) for k in range(NN) for e in range(NETA)
                if np.any(C2[:, k, e] != 0)]
        hist["cell_max"].append(max(cell) if cell else 0.0)
        W = np.fft.ifft2(C2, axes=(1, 2)) * (NN * NETA)
        dual = [rel_bil(W[:, n, h]) for n in range(NN) for h in range(NETA)
                if np.any(np.abs(W[:, n, h]) > 1e-200)]
        hist["dual_max"].append(max(dual) if dual else 0.0)
        if t >= WIN[0] and (t - WIN[0]) % EV == 0:
            samples.append(C2[:, 2, 0].copy())     # 凝縮体セル（ポンプ k=2・巻き0）
    S = np.array(samples)
    for k in hist:
        hist[k] = np.array(hist[k])
    w = slice(WIN[0], WIN[1])
    print("=== (B1) 閉塞の所在（安定窓の中央値・|zᵀz|/z†z） ===")
    loc = {}
    for k, lab in (("total", "全体和 Σ_{e,k,η}"), ("cell_max", "セル別 Σ_e（最大）"),
                   ("dual_max", "双対点別 Σ_e（最大）")):
        med = float(np.median(hist[k][w]))
        loc[k] = med
        print(f"  {lab:26s} = {med:.3e}  {'← 厳密ゼロ' if med < TOL else ''}")

    # ---- フレーム抽出（公開論文と同じ処方: Z⊥ の実表現SVD）
    Sp = np.array([z - p2 * (p2 @ z) - q2 * (q2 @ z) for z in S])
    X = np.hstack([Sp.real, Sp.imag])
    Xc = X - X.mean(axis=0)
    _, sv, Vt = np.linalg.svd(Xc, full_matrices=False)
    M = S.shape[1]
    svr = sv / sv[0]
    print(f"\n=== 特異値の梯子（上位10・相対） ===\n  {np.round(svr[:10], 4)}")

    def cdir(k):
        u = Vt[k]
        d = u[:M] + 1j * u[M:]
        return d / np.linalg.norm(d)

    n_planes = min(6, len(sv) // 2)
    dirs = [cdir(2 * j) for j in range(n_planes)]

    # ---- B2: 平面の自己等方性
    print("\n=== (B2) 平面方向の自己等方性 |dᵀd|/(d†d) ===")
    b2 = []
    for j, d in enumerate(dirs):
        r = abs(bil(d)) / float(np.real(np.vdot(d, d)))
        b2.append(r)
        print(f"  平面{j+1}（特異値比 {svr[2*j]:.4f}）: {r:.3e} "
              f"{'← 等方' if r < TOL else ''}")

    # ---- B3: 平面間の等方性
    print("\n=== (B3) 平面間 |d_jᵀd_k|/(|d_j||d_k|) ===")
    b3 = []
    for j in range(n_planes):
        for k in range(j + 1, n_planes):
            r = abs(complex(np.sum(dirs[j] * dirs[k])))
            b3.append(r)
    print(f"  最大 {max(b3) if b3 else 0.0:.3e}")

    # ---- B4: 子の閉塞（状態を各平面へ射影した成分の閉塞）
    print("\n=== (B4) 子の閉塞 |Σ_e(Z_k)²|/Σ_e|Z_k|²（安定窓の中央値） ===")
    b4 = []
    for j, d in enumerate(dirs):
        rs = []
        for z in Sp:
            c = np.vdot(d, z)              # エルミート射影係数
            rs.append(rel_bil(c * d))
        med = float(np.median(rs))
        b4.append(med)
        print(f"  平面{j+1} 成分: {med:.3e} {'← 厳密ゼロ' if med < TOL else ''}")

    # ---- B5: 対照（ランダム部分集合・ランダム方向）
    rng = np.random.default_rng(20260808)
    zrep = Sp[len(Sp) // 2]
    ctrl_subset = []
    for _ in range(200):
        k = rng.integers(2, M)
        idx = rng.choice(M, size=k, replace=False)
        ctrl_subset.append(rel_bil(zrep[idx]))
    ctrl_dir = []
    for _ in range(200):
        d = rng.normal(size=M) + 1j * rng.normal(size=M)
        d /= np.linalg.norm(d)
        ctrl_dir.append(rel_bil(np.vdot(d, zrep) * d))
    print("\n=== (B5) 対照 ===")
    print(f"  ランダム部分集合（辺の部分集合・200回）: 中央値 "
          f"{np.median(ctrl_subset):.3e}  最小 {np.min(ctrl_subset):.3e}")
    print(f"  ランダム方向への射影成分（200回）: 中央値 "
          f"{np.median(ctrl_dir):.3e}  最小 {np.min(ctrl_dir):.3e}")

    B2 = all(r < TOL for r in b2)
    B3 = (max(b3) if b3 else 0.0) < TOL
    B4 = all(r < TOL for r in b4)
    B5 = float(np.median(ctrl_subset)) > CTRL_MIN and float(np.median(ctrl_dir)) > CTRL_MIN
    n_child = int(sum(1 for r in b4 if r < TOL))
    print(f"\n(B2) 平面の自己等方性: {'通過' if B2 else '不成立'}")
    print(f"(B3) 平面間の等方性: {'通過' if B3 else '不成立'}")
    print(f"(B4) 子の閉塞（平面成分が厳密ゼロ閉塞）: {'通過' if B4 else '不成立'}")
    print(f"(B5) 対照が非ゼロ: {'通過' if B5 else '不成立'}")
    print(f"(B6) 子閉塞を満たす平面の数 = {n_child} / {n_planes}")

    out = {"env": {"N": N, "Nn": NN, "Neta": NETA, "T": T, "delta": DELTA,
                   "window": WIN, "tol": TOL},
           "B1_closure_location": loc,
           "singular_ladder_rel": [float(x) for x in svr[:12]],
           "B2_plane_isotropy": b2, "B3_cross_isotropy_max": float(max(b3) if b3 else 0.0),
           "B4_child_closure": b4,
           "B5_control": {"subset_median": float(np.median(ctrl_subset)),
                          "subset_min": float(np.min(ctrl_subset)),
                          "dir_median": float(np.median(ctrl_dir)),
                          "dir_min": float(np.min(ctrl_dir))},
           "verdict": {"B2": bool(B2), "B3": bool(B3), "B4": bool(B4),
                       "B5": bool(B5), "n_child_blocks": n_child,
                       "n_planes_tested": n_planes},
           "runtime_sec": time.time() - t0}
    (HERE / "result_tb_child_closure_blocks_v1.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False, default=float))

    ts = np.arange(1, T + 1)
    fig, ax = plt.subplots(2, 1, figsize=(9, 7))
    ax[0].semilogy(ts, np.maximum(hist["total"], 1e-20), lw=0.7, label="全体和")
    ax[0].semilogy(ts, np.maximum(hist["cell_max"], 1e-20), lw=0.7, label="セル別（最大）")
    ax[0].semilogy(ts, np.maximum(hist["dual_max"], 1e-20), lw=0.7, label="双対点別（最大）")
    ax[0].set_ylabel("|Σxₙ²| / Σ|xₙ|²"); ax[0].set_xlabel("τ（step）")
    ax[0].legend(fontsize=8); ax[0].set_title("閉塞の所在（監査量・統一Gの外で直接計測）")
    idx = np.arange(1, len(b4) + 1)
    ax[1].semilogy(idx, np.maximum(b4, 1e-20), "o-", label="平面成分の閉塞残差（子）")
    ax[1].axhline(float(np.median(ctrl_dir)), color="red", ls="--", lw=0.9,
                  label="対照: ランダム方向の中央値")
    ax[1].set_xlabel("平面（特異値の梯子の順）"); ax[1].set_ylabel("|Σ(Z_k)²|/Σ|Z_k|²")
    ax[1].legend(fontsize=8); ax[1].set_title("子の閉塞: 凝縮体＝ゼロ閉塞する部分塊か")
    fig.tight_layout()
    fig.savefig(HERE / "fig_tb_child_closure_v1.png", dpi=130)
    plt.close(fig)
    print(f"完了 {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
