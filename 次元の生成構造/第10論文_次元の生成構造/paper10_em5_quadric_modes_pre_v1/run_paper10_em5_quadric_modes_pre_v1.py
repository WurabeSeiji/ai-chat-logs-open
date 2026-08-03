#!/usr/bin/env python3
"""E-M5：ロック次数が決める状態空間の次元と E8 判定（Stage 0-2）v1

設計書: ../E-M5設計_ロック次数と状態空間次元とE8判定.md

中心同定 H5:
    N は最小分解能であって体の数ではない。粒子閉鎖はロック次数 g
    （結晶学的メニュー g∈{3,4,6}）で特徴づけられ、比だけが物理なので
    基準セル1個（差の不在）を除いた有効分解能 N = g−1。
    関係波 M = N(N−1)/2 本、零閉塞＋射影化でコンパクトなクアドリック
    Q_{M−2}（複素次元 M−2）。
        g=3: N=2, M=1  → 空集合（状態空間なし）
        g=4: N=3, M=3  → Q₁ = CP¹（リーマン球面）
        g=6: N=5, M=10 → Q₈（複素次元8、E8 候補の住処）

固定予言（実行前に固定。事後変更禁止。反証も記録する）:
    P1（次元則）  : 制約 Φ=(‖Z‖²−1, Re Z·Z, Im Z·Z) のヤコビアンの
                    ランクは 3、ゲージ方向 iZ は核に入り、物理接空間の
                    実次元 = 2(M−2)。Z₆型で16、Z₄型で2。整数なので厳密一致のみ PASS。
                    g=3（M=1）は |Z·Z|=‖Z‖²=1 が恒等的に成り立ち空集合。
    P2（楕円性）  : 準安定閉鎖点の回転枠一段写像の接空間線形化の固有値は
                    全て単位円上（max||λ|−1| < 1e-8）。Z₆型で8本の周波数対。
    P3（球面ラダー）: Z₄型（CP¹）の非線形モード信号のスペクトルは
                    等間隔ラダー（倍音間隔の偏差 < 1e-6 [cycle/step]）。
    P4（物質誕生の階層仮説・第8論文v2ゲート）:
                    軌道レベルで形成されなかった整数比ロックが、モード座標の
                    周波数間では形成される。判定は E-M4 と同一の機械的基準
                    （比の丸め偏差 < 1e-3、比≥2 を非自明ロックと数える）。
                    不成立なら「場の円環化」予想の反証として記録。
    P5（E8指紋・条件付き）: P1・P2・P4 が全て PASS の場合のみ Stage 3
                    （作用格子の Gram 監査 → E8 判定バッテリー8項）を別途実施。

規約:
    無名性——粒子種・ロック次数による IF 分岐なし。初期データは g 次配置の
    位相のみで指定し（機械的な閉塞射影）、判定はランク・スペクトル・
    ロック判定の機械計算。
    シリーズ内再現——一段写像は E-M4 の対照密形式（エンジン自身の kmatvec
    検証形 Kd=A·sin(Δθ) の Cayley）を read-only import 相当で厳密再現し、
    import 元の SHA-256 を記録する。

パラメータ（実行前固定）:
    T_RELAX = 12000（E-M4 XMAX と同一）, T_LONG = 48000（E-M4 追試と同一）,
    励起 EPS_EXC = 1e-3（全物理モード等励起）, ヤコビアン中心差分 H = 1e-6。
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PAPER8 = HERE.parent.parent / "第8論文_二段階seed除去による準安定相の因果分離"
EM4 = PAPER8 / "paper8_axiom2_coupling_harmonics_pre_v1" / "run_paper8_axiom2_coupling_harmonics_pre_v1.py"
ABL = PAPER8 / "code" / "run_preliminary_seed_ablation_v1.py"

# read-only import（E-M4 経由でエンジンも読み込まれる）
spec = importlib.util.spec_from_file_location("em4_m5", EM4)
em4 = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = em4
spec.loader.exec_module(em4)
eng = em4.eng                     # run_n_scaling_lowrank_v1
ENGINE_FILE = Path(eng.__file__)

GAMMA = eng.GAMMA                 # tan(π/144)（エンジンと同一）
T_RELAX = 12000
T_LONG = 48000
EPS_EXC = 1e-3
H_FD = 1e-6

TOL_P2 = 1e-8
TOL_P3 = 1e-6
TOL_LOCK = 1e-3                   # E-M4 と同一


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


# ---------- 構成（無名：g 次配置の位相のみで指定） ----------

def shared_vertex_adjacency(sys_lr):
    """E-M4 dense_setup と同一の隣接構成（62-68行）。"""
    ea, eb = sys_lr.ea, sys_lr.eb
    m = sys_lr.m
    A = np.zeros((m, m))
    for i in range(m):
        share = (ea == ea[i]) | (ea == eb[i]) | (eb == ea[i]) | (eb == eb[i])
        A[i, share] = 1.0
    np.fill_diagonal(A, 0.0)
    return A


def closure_project(Z, tol=1e-14, iters=100):
    """Z·Z=0, ‖Z‖=1 への機械的射影（Newton＋正規化、乱数なし）。"""
    Z = Z / np.linalg.norm(Z)
    for _ in range(iters):
        c = complex(Z @ Z)
        if abs(c) < tol:
            break
        Z = Z - (c / 2.0) * np.conj(Z) / float(np.real(np.conj(Z) @ Z))
        Z = Z / np.linalg.norm(Z)
    return Z


def lock_config(g):
    """g 次ロック配置：セル j=1..g−1 の位相 θ_j=2πj/g、辺 (j,k) は位相差。"""
    n_cells = g - 1
    sys_lr = eng.LowRankSystem(n_cells)
    m = sys_lr.m
    theta_cell = 2.0 * np.pi * (np.arange(n_cells) + 1) / g
    Z0 = np.exp(1j * (theta_cell[sys_lr.ea] - theta_cell[sys_lr.eb])) / np.sqrt(m)
    Z0 = closure_project(Z0)
    A = shared_vertex_adjacency(sys_lr)
    return sys_lr, A, Z0, m


# ---------- 一段写像（E-M4 対照密形式 91-99行の厳密再現、第1倍音のみ） ----------

def step(A, Z):
    th = np.angle(Z)
    K = A * np.sin(th[None, :] - th[:, None])
    sigma = np.linalg.norm(K, 2)
    if sigma == 0.0:
        return Z.copy()
    gs = GAMMA / sigma
    eye = np.eye(len(Z))
    return np.linalg.solve(eye - gs * K, (eye + gs * K) @ Z)


# ---------- Stage 0: 制約ランク検査 ----------

def dphi(Z):
    """Φ=(‖Z‖²−1, Re Z·Z, Im Z·Z) の実ヤコビアン 3×2M。実座標 x=(Re Z, Im Z)。"""
    x, y = Z.real, Z.imag
    return np.vstack([
        np.concatenate([2 * x, 2 * y]),
        np.concatenate([2 * x, -2 * y]),
        np.concatenate([2 * y, 2 * x]),
    ])


def rank_check(Z):
    m = len(Z)
    D = dphi(Z)
    sv = np.linalg.svd(D, compute_uv=False)
    rank = int(np.sum(sv > 1e-8 * sv[0]))
    gauge = np.concatenate([(1j * Z).real, (1j * Z).imag])
    gauge_in_kernel = float(np.linalg.norm(D @ gauge))
    dim_phys = 2 * m - rank - 1
    return {"rank_dPhi": rank, "gauge_kernel_residual": gauge_in_kernel,
            "dim_physical_tangent": dim_phys, "expected": 2 * (m - 2)}


def tangent_basis(Z):
    """物理接空間（ker dΦ ∩ gauge⊥）の正規直交基底 2M×(2M−4)。"""
    m = len(Z)
    gauge = np.concatenate([(1j * Z).real, (1j * Z).imag])
    V = np.vstack([dphi(Z), gauge])          # 4×2M
    _, _, vt = np.linalg.svd(V, full_matrices=True)
    return vt[4:].T                           # 直交補空間


# ---------- Stage 1: 回転枠一段写像のモードスペクトル ----------

def stage1(sys_lr, A, Z0, label):
    m = len(Z0)
    Z = Z0.copy()
    for _ in range(T_RELAX):
        Z = step(A, Z)
    Zs = Z
    FZ = step(A, Zs)
    phi_star = float(np.angle(np.conj(Zs) @ FZ))
    rot = np.exp(-1j * phi_star)
    residual = float(np.linalg.norm(rot * FZ - Zs))

    def G_real(xr):
        Zi = xr[:m] + 1j * xr[m:]
        Zo = rot * step(A, Zi)
        return np.concatenate([Zo.real, Zo.imag])

    x0 = np.concatenate([Zs.real, Zs.imag])
    J = np.zeros((2 * m, 2 * m))
    for k in range(2 * m):
        e = np.zeros(2 * m); e[k] = H_FD
        J[:, k] = (G_real(x0 + e) - G_real(x0 - e)) / (2 * H_FD)

    B = tangent_basis(Zs)
    Mred = B.T @ J @ B
    lam, V = np.linalg.eig(Mred)
    order = np.argsort(-np.abs(np.angle(lam)))
    lam, V = lam[order], V[:, order]
    max_dev_circle = float(np.max(np.abs(np.abs(lam) - 1.0)))
    freqs = np.abs(np.angle(lam))            # rad/step（共役対で重複）
    pos = np.sort(freqs[np.angle(lam) > 0])  # 正側のみ＝周波数対
    print(f"  [{label}] Stage1: 残差={residual:.3e} max||λ|-1|={max_dev_circle:.3e} "
          f"周波数対の本数={len(pos)}（期待 {m-2}）")
    return {"Zs": Zs, "phi_star": phi_star, "residual": residual,
            "Mred": Mred, "lam": lam, "V": V, "B": B,
            "max_dev_circle": max_dev_circle, "freq_pairs": pos}


# ---------- Stage 2: モード座標の周波数とロック判定（E-M4 基準） ----------

def em4_lock_stats(fb):
    """E-M4 111-116行と同一の判定。"""
    if fb.size > 1:
        r = fb[:, None] / np.maximum(fb[None, :], 1e-30)
        rmax = np.maximum(r, 1 / np.maximum(r, 1e-30))
        max_dev = float(np.max(np.abs(r[r >= 1] - 1)))
        pr = np.round(rmax)
        locks = int(np.sum((pr >= 2) & (np.abs(rmax - pr) < TOL_LOCK)) // 2)
        ratio_matrix = rmax.tolist()
    else:
        max_dev, locks, ratio_matrix = 0.0, 0, []
    return max_dev, locks, ratio_matrix


def stage2(sys_lr, A, s1, label):
    Zs, phi_star, B, V = s1["Zs"], s1["phi_star"], s1["B"], s1["V"]
    m = len(Zs)
    rot = np.exp(-1j * phi_star)
    # 全物理モードを等励起（機械的・乱数なし）
    exc = B @ np.ones(B.shape[1])
    exc = exc / np.linalg.norm(exc)
    Zp = (Zs.real + 1j * Zs.imag) + EPS_EXC * (exc[:m] + 1j * exc[m:])
    Zp = closure_project(Zp)

    n_modes = B.shape[1]
    coords = np.zeros((T_LONG, n_modes), dtype=complex)
    Z = Zp.copy()
    frame = 1.0 + 0.0j
    for t in range(T_LONG):
        Z = step(A, Z)
        frame *= rot
        dz = frame * Z - Zs
        y = B.T @ np.concatenate([dz.real, dz.imag])
        coords[t] = np.linalg.solve(V, y)

    # 共役対の正側モードのみ周波数抽出（位相のアンラップ勾配＝E-M4 観測量）
    lam = s1["lam"]
    idx_pos = [i for i in range(len(lam)) if np.angle(lam[i]) > 0]
    freqs, amps = [], []
    for i in idx_pos:
        a = coords[:, i]
        amp = float(np.mean(np.abs(a)))
        u = np.unwrap(np.angle(a[max(1, T_LONG - 3000):]))
        fr = abs(float(np.polyfit(np.arange(len(u)), u, 1)[0]))
        freqs.append(fr); amps.append(amp)
    fb = np.asarray(freqs)
    fb_valid = fb[fb > 1e-8]
    max_dev, locks, ratio_matrix = em4_lock_stats(fb_valid)
    print(f"  [{label}] Stage2: モード周波数={np.array2string(fb, precision=6)} "
          f"最大比ずれ={max_dev:.3e} 非自明ロック={locks}")
    out = {"mode_freqs_rad_per_step": fb.tolist(), "mode_mean_amps": amps,
           "freq_max_ratio_dev": max_dev, "nontrivial_locks": locks,
           "ratio_matrix": ratio_matrix,
           "linear_freq_pairs": s1["freq_pairs"].tolist()}

    # Z₄型（単一モード）のみ：倍音ラダーの等間隔検査（P3）
    if n_modes == 2:
        sig = coords[:, idx_pos[0]].real
        w = np.hanning(T_LONG)
        F = np.abs(np.fft.rfft(sig * w))
        fgrid = np.fft.rfftfreq(T_LONG)     # cycle/step
        peaks = []
        for k in range(2, len(F) - 2):
            if F[k] > F[k - 1] and F[k] > F[k + 1] and F[k] > 1e-9 * F.max():
                # 二次補間
                d = 0.5 * (F[k - 1] - F[k + 1]) / (F[k - 1] - 2 * F[k] + F[k + 1])
                peaks.append((fgrid[k] + d * (fgrid[1] - fgrid[0]), F[k]))
        peaks.sort(key=lambda x: -x[1])
        top = sorted([p[0] for p in peaks[:3]])
        if len(top) >= 3:
            spacings = np.diff(top)
            ladder_dev = float(abs(spacings[1] - spacings[0]))
        elif len(top) == 2:
            ladder_dev = float(abs(top[1] - 2 * top[0]))
        else:
            ladder_dev = float("nan")
        out["ladder_peaks_cyc_per_step"] = top
        out["ladder_spacing_dev"] = ladder_dev
        print(f"  [{label}] P3 ラダー: ピーク={top} 間隔偏差={ladder_dev:.3e}")
    return out


def main() -> None:
    t0 = time.time()
    print("E-M5 Stage 0-2 実行")
    print(f"  import: EM4  sha256={sha256(EM4)[:16]}…")
    print(f"  import: ABL  sha256={sha256(ABL)[:16]}…")
    print(f"  import: ENG  sha256={sha256(ENGINE_FILE)[:16]}…")

    results = {"imports": {"em4": sha256(EM4), "abl": sha256(ABL), "engine": sha256(ENGINE_FILE)},
               "params": {"GAMMA": GAMMA, "T_RELAX": T_RELAX, "T_LONG": T_LONG,
                           "EPS_EXC": EPS_EXC, "H_FD": H_FD,
                           "TOL_P2": TOL_P2, "TOL_P3": TOL_P3, "TOL_LOCK": TOL_LOCK}}

    # ---- g=3（M=1）: 空集合の検査 ----
    print("\n[g=3] N=2, M=1（空集合の検査）")
    zs = np.exp(1j * np.linspace(0, 2 * np.pi, 7, endpoint=False))
    min_closure = float(min(abs(z * z) for z in zs))
    empty_ok = abs(min_closure - 1.0) < 1e-12
    print(f"  |Z·Z| の最小値（単位円上サンプル）= {min_closure:.12f}（恒等的に1、状態空間なし）: "
          f"{'PASS' if empty_ok else 'FAIL'}")
    results["g3"] = {"min_closure_abs": min_closure, "empty_state_space": bool(empty_ok)}

    # ---- g=4, g=6 ----
    per_g = {}
    for g, label in ((4, "Z4"), (6, "Z6")):
        print(f"\n[g={g}] N={g-1}, M={(g-1)*(g-2)//2}")
        sys_lr, A, Z0, m = lock_config(g)
        rc0 = rank_check(Z0)
        print(f"  [{label}] Stage0(初期): rank={rc0['rank_dPhi']} "
              f"gauge残差={rc0['gauge_kernel_residual']:.2e} "
              f"物理接次元={rc0['dim_physical_tangent']}（期待 {rc0['expected']}）")
        s1 = stage1(sys_lr, A, Z0, label)
        rcs = rank_check(s1["Zs"])
        print(f"  [{label}] Stage0(準安定点): rank={rcs['rank_dPhi']} "
              f"物理接次元={rcs['dim_physical_tangent']}（期待 {rcs['expected']}）")
        s2 = stage2(sys_lr, A, s1, label)
        per_g[label] = {
            "m": m,
            "stage0_init": rc0, "stage0_meta": rcs,
            "stage1": {"phi_star": s1["phi_star"], "residual": s1["residual"],
                        "max_dev_circle": s1["max_dev_circle"],
                        "freq_pairs": s1["freq_pairs"].tolist(),
                        "n_freq_pairs": len(s1["freq_pairs"]), "expected_pairs": m - 2},
            "stage2": s2,
        }

    # ---- 判定 ----
    z4, z6 = per_g["Z4"], per_g["Z6"]
    p1 = (results["g3"]["empty_state_space"]
          and z4["stage0_init"]["dim_physical_tangent"] == 2
          and z4["stage0_meta"]["dim_physical_tangent"] == 2
          and z6["stage0_init"]["dim_physical_tangent"] == 16
          and z6["stage0_meta"]["dim_physical_tangent"] == 16)
    p2 = (z4["stage1"]["max_dev_circle"] < TOL_P2
          and z6["stage1"]["max_dev_circle"] < TOL_P2
          and z6["stage1"]["n_freq_pairs"] == 8
          and z4["stage1"]["n_freq_pairs"] == 1)
    ladder = z4["stage2"].get("ladder_spacing_dev", float("nan"))
    p3 = bool(np.isfinite(ladder) and ladder < TOL_P3)
    p4 = z6["stage2"]["nontrivial_locks"] > 0
    p5_gate = p1 and p2 and p4

    print("\n==== 判定 ====")
    print(f"P1 次元則（空/2/16）: {'PASS' if p1 else 'FAIL'}")
    print(f"P2 楕円性（max||λ|-1|<1e-8・対の本数 1/8）: {'PASS' if p2 else 'FAIL'} "
          f"(Z4:{z4['stage1']['max_dev_circle']:.2e}, Z6:{z6['stage1']['max_dev_circle']:.2e})")
    print(f"P3 球面等間隔ラダー（<1e-6）: {'PASS' if p3 else 'FAIL'} (偏差 {ladder:.3e})")
    print(f"P4 モードレベル整数ロック（v2ゲート）: {'PASS' if p4 else 'FAIL'} "
          f"(非自明ロック {z6['stage2']['nontrivial_locks']})")
    print(f"P5 ゲート（Stage3 E8バッテリーへ進むか）: {'開' if p5_gate else '閉'}")

    results["per_g"] = per_g
    results["verdicts"] = {"P1": bool(p1), "P2": bool(p2), "P3": bool(p3),
                            "P4": bool(p4), "P5_gate_open": bool(p5_gate)}
    results["runtime_sec"] = time.time() - t0
    (HERE / "paper10_em5_result_v1.json").write_text(
        json.dumps(results, ensure_ascii=False, indent=2, default=float), encoding="utf-8")
    print(f"\nsaved ({results['runtime_sec']:.0f}s)")


if __name__ == "__main__":
    main()
