"""第十三思考実験 §12.2〜12.6 の追加実験（厳密根 R_{124,23} = cos^2(23π/124) の系）。

原本（2026-07-13 予備実験、SHA256 f815320f…）の関数で初期状態・散乱係数・倍音分布を作り、
更新は原本と同じ交換作用素 U_R = [[r, t], [t, r]] を使う。
  P（局在波）= 原本の B（奇数倍音 1〜63、m=2、q=-1）、背景 A_j = 原本の A（単音、m=1、q=+1）の写し。
多チャネルではチャネルごとのノルムが 1 に保たれないため normalize を外した純粋な U_R を使う
（閉じた二波では両者が丸め誤差で一致することは [S7] で確認済み）。
各チャネルは常に x a0 + y b0 に留まるので、E12_3〜E12_5 は係数 (x, y) で計算し、全状態の計算と照合する。

E12_2 (§12.2) 閉条件 0..248 衝突：倍音成分の復元残差、倍音ごとのパワーの保存、二乗閉塞 Σψ²、
              エネルギー候補（原本の倍音分布の一次 E1・二次 E2 モーメント）を別々に記録
E12_3 (§12.3) 閉（M=1）と開：M 本の背景から毎回ラベルに依らず無作為に相手を選ぶ（M=2,4,16,64）、
              毎回新しい背景（M=∞、局在成分は厳密に R^k）。流出・再流入・局在の残存を記録
E12_4 (§12.4) R（0.01〜0.99 と有限位数根 n≤8, 124）と背景の初期位相の掃引。整数比を設定として与えない
E12_5 (§12.5) 交換しない基準波 F に対する P の自己振幅 y_P(k)（内部時計の読み）と N_eff（局在の読み）
E12_6 (§12.6) 反証条件 1〜4 の判定を summary に記す
"""
from __future__ import annotations

import csv
import importlib.util
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
SRC = HERE / "original_copy" / "run_exchange_scattering_matrix_fermionic_localization_transfer_preliminary_v1.py"
OUT = HERE / "results"
OUT.mkdir(exist_ok=True)
spec = importlib.util.spec_from_file_location("orig", SRC)
orig = importlib.util.module_from_spec(spec)
sys.modules["orig"] = orig
spec.loader.exec_module(orig)

PR = orig.Params()
R0 = math.cos(math.pi * 23 / 124) ** 2
A0 = orig.make_state(PR, 1, PR.q_A, PR.m_A, True, PR.A_A)
B0 = orig.make_state(PR, 63, PR.q_B, PR.m_B, True, PR.A_B)
NH = PR.chi_grid_n // 2 + 1
K3, S3, MS = 1000, 64, (1, 2, 4, 16, 64)


def tr(R):
    t, r, _, _ = orig.scattering_coefficients(orig.delta_from_reflection_rate(R))
    return complex(t), complex(r)


def hpow(v):
    """倍音次数 n ごとのパワー（原本の harmonic_distribution × ノルム²）。"""
    out = np.zeros(NH)
    for n, w in orig.harmonic_distribution(PR, v).items():
        out[n] = orig.norm2(v) * w
    return out


def mom(v):
    d = orig.harmonic_distribution(PR, v)
    return sum(n * w for n, w in d.items()), sum(n * n * w for n, w in d.items())


E1A, E2A = mom(A0)
E1B, E2B = mom(B0)


def write(name, head, rows, fmt):
    with open(OUT / name, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(head)
        for r in rows:
            w.writerow([x if isinstance(x, str) else (fmt % x if isinstance(x, float) else x) for x in r])


def e12_2():
    t, r = tr(R0)
    a, b = A0.copy(), B0.copy()
    FA0, FB0 = orig.chi_frequency_components(PR, A0), orig.chi_frequency_components(PR, B0)
    sc = max(np.abs(FA0).max(), np.abs(FB0).max())
    H0 = hpow(A0) + hpow(B0)
    om = math.pi + 2 * math.asin(math.sqrt(R0))
    rows, rec, hdev, cf, s2max = [], 0.0, 0.0, 0.0, [0.0, 0.0, 0.0]
    for k in range(249):
        x, y, u, v = orig.inner(A0, a), orig.inner(B0, a), orig.inner(A0, b), orig.inner(B0, b)
        FA, FB = orig.chi_frequency_components(PR, a), orig.chi_frequency_components(PR, b)
        rec = max(rec, np.abs(FA - x * FA0 - y * FB0).max() / sc, np.abs(FB - u * FA0 - v * FB0).max() / sc)
        hdev = max(hdev, float(np.abs(hpow(a) + hpow(b) - H0).max()))
        (e1a, e2a), (e1b, e2b) = mom(a), mom(b)
        s = math.sin(k * om / 2) ** 2
        cf = max(cf, abs(e2a - ((1 - s) * E2A + s * E2B)), abs(e2b - (s * E2A + (1 - s) * E2B)))
        s2 = [abs(complex(np.sum(a * a))), abs(complex(np.sum(b * b))), abs(complex(np.sum(a * a) + np.sum(b * b)))]
        s2max = [max(p, q) for p, q in zip(s2max, s2)]
        if k <= 124:
            rows.append((k, s2[2], e2a, e2b))
        a, b = r * a + t * b, t * a + r * b
    write("e12_2_state_tracking_R12423.csv", ["collision", "abs_sum_psi2_total", "E2_A", "E2_B"], rows, "%.5g")
    return {"harmonic_recon_residual_rel": rec, "harmonic_power_total_max_dev": hdev,
            "max_abs_sum_psi2_A_B_total": s2max, "E1_a0_b0": [E1A, E1B], "E2_a0_b0": [E2A, E2B],
            "max_abs_E2_minus_closed_form": cf}


def pool(R, M, K, S, seed, theta=False, hist=False):
    """P と M 本の背景。毎回、各試行で背景を一様に無作為に選んで U_R で交換する（M=1 は閉条件）。"""
    t, r = tr(R)
    rp, rt = np.random.default_rng(seed), np.random.default_rng(seed + 1)
    ix = np.arange(S)
    X = np.zeros((S, M + 1), complex)
    Y = np.zeros((S, M + 1), complex)
    Y[:, 0] = 1.0
    X[:, 1:] = np.exp(1j * rt.uniform(0, 2 * np.pi, (S, M))) if theta else 1.0
    W = np.empty((K + 1, S))
    W[0] = 1.0
    HX, HY = [X[:, 0].copy()], [Y[:, 0].copy()]
    for k in range(K):
        j = rp.integers(1, M + 1, S)
        for Z in (X, Y):
            p, q = Z[ix, 0].copy(), Z[ix, j].copy()
            Z[ix, 0], Z[ix, j] = r * p + t * q, t * p + r * q
        W[k + 1] = np.abs(Y[:, 0]) ** 2
        if hist:
            HX.append(X[:, 0].copy())
            HY.append(Y[:, 0].copy())
    return W, (np.array(HX), np.array(HY)), (X, Y)


def fresh(R, K):
    """毎回新しい背景（同位相の単音）を相手にする全開放。"""
    t, r = tr(R)
    x, y, HX, HY = 0j, 1 + 0j, [0j], [1 + 0j]
    for k in range(K):
        x, y = r * x + t, r * y
        HX.append(x)
        HY.append(y)
    return np.array(HX), np.array(HY)


def neff(x, y):
    p, q = np.abs(x) ** 2, np.abs(y) ** 2
    return (p * E1A + q * E1B) / (p + q)


def lifetime(Wm):
    """残存率 f(k)=(w-w_∞)/(1-w_∞) が 1/e を切る衝突数（線形補間）と w_∞（後半の平均）。"""
    pl = float(Wm[len(Wm) // 2:].mean())
    f = (Wm - pl) / (1 - pl)
    hit = np.nonzero(f <= 1 / math.e)[0]
    if not hit.size or hit[0] == 0:
        return None, pl
    k = int(hit[0])
    return k - 1 + float((f[k - 1] - 1 / math.e) / (f[k - 1] - f[k])), pl


def e12_3():
    t, r = tr(R0)
    ks = list(range(131)) + list(range(135, 301, 15))
    cols, info = {}, {}
    for M in MS:
        W, _, (X, Y) = pool(R0, M, K3, S3, 100 + M)
        Wm = W.mean(1)
        d = np.diff(W, axis=0)
        tau, pl = lifetime(Wm)
        info[f"M{M}"] = {"tau_e": None if M == 1 else tau, "plateau": pl, "w_at_124_mean": float(Wm[124]),
                         "max_w_k10_to_K_mean_over_seeds": float(W[10:].max(0).mean()),
                         "reinflow_total_mean": float(np.clip(d, 0, None).sum(0).mean()),
                         "outflow_total_mean": float(np.clip(-d, 0, None).sum(0).mean()),
                         "total_power_max_dev": float(np.abs((np.abs(X) ** 2 + np.abs(Y) ** 2).sum(1) - (M + 1)).max())}
        cols[f"M{M}"] = Wm
    _, HY = fresh(R0, K3)
    cols["fresh"] = np.abs(HY) ** 2
    info["fresh"] = {"tau_exact": -1 / math.log(R0), "max_abs_w_minus_R_pow_k": float(np.abs(cols["fresh"] - R0 ** np.arange(K3 + 1)).max())}
    # 全状態との照合（M=4、試行 1、120 衝突、同じ相手の系列）
    _, (hx, hy), _ = pool(R0, 4, 120, 1, 7, hist=True)
    rp, P, bg, res = np.random.default_rng(7), B0.copy(), [A0.copy() for _ in range(4)], 0.0
    for k in range(120):
        j = int(rp.integers(1, 5, 1)[0]) - 1
        P, bg[j] = r * P + t * bg[j], t * P + r * bg[j]
        res = max(res, float(np.linalg.norm(P - hx[k + 1, 0] * A0 - hy[k + 1, 0] * B0)))
    info["full_state_check_M4_120_collisions"] = res
    names = ["M1", "M2", "M4", "M16", "M64", "fresh"]
    write("e12_3_open_vs_closed_R12423.csv", ["collision", "closed_M1", "pool_M2", "pool_M4", "pool_M16", "pool_M64", "fresh"],
          [(k, *[float(cols[n][k]) for n in names]) for k in ks], "%.4g")
    return info


def e12_4():
    roots = sorted({(n, m) for n in range(3, 9) for m in range(1, n) if 2 * m < n and math.gcd(m, n) == 1} | {(124, 23)})

    def tp(R, M):
        W, _, _ = pool(R, M, 1000, 32, 4000 + M)
        return lifetime(W.mean(1))

    def mineps(R):
        t, r = tr(R)
        U, Z, best = np.array([[r, t], [t, r]]), np.eye(2, dtype=complex), 9.0
        for k in range(1, 1001):
            Z = U @ Z
            if k >= 10:
                best = min(best, math.sqrt(float(np.sum(np.abs(Z - np.eye(2)) ** 2)) / 2))
        return best

    pts = [(round(0.01 * i, 2), "grid") for i in range(1, 100)] + [(math.cos(math.pi * m / n) ** 2, f"root_{n}_{m}") for n, m in roots]
    rows, nb = [], {"M2": 0.0, "M4": 0.0}
    for R, tag in pts:
        (t2, p2), (t4, p4) = tp(R, 2), tp(R, 4)
        rows.append((tag, R, 0.5 + math.asin(math.sqrt(R)) / math.pi, mineps(R), t2, t4, p2, p4))
        if tag != "grid":
            for M, tv in (("M2", t2), ("M4", t4)):
                nbr = [tp(R + s, int(M[1:]))[0] for s in (-0.003, 0.003)]
                nb[M] = max(nb[M], abs(tv - sum(nbr) / 2) / (sum(nbr) / 2))
    ref = {"M2": 0.0, "M4": 0.0}
    for q in rows:
        if q[0] == "grid" and q[3] > 1e-10:
            for M, tv in (("M2", q[4]), ("M4", q[5])):
                nbr = [tp(q[1] + s, int(M[1:]))[0] for s in (-0.003, 0.003)]
                ref[M] = max(ref[M], abs(tv - sum(nbr) / 2) / (sum(nbr) / 2))
    write("e12_4_sweep_R.csv", ["point", "R", "omega_over_2pi", "min_eps_closed_k10_1000", "tau_e_M2", "tau_e_M4", "plateau_M2", "plateau_M4"], rows, "%.5g")
    g = [r for r in rows if r[0] == "grid" and r[3] > 1e-10]
    rt = [r for r in rows if r[0] != "grid"]
    W0, (x0, y0), _ = pool(R0, 4, 1000, 32, 77, theta=False, hist=True)
    W1, (x1, y1), _ = pool(R0, 4, 1000, 32, 77, theta=True, hist=True)
    return {"n_grid": 99, "n_grid_exact_return": 99 - len(g), "roots": [r[0] for r in rt],
            "tau_e_range_grid_M2_M4": [min(r[4] for r in g), max(r[4] for r in g), min(r[5] for r in g), max(r[5] for r in g)],
            "max_rel_tau_diff_root_vs_neighbors_pm0003": nb, "max_rel_tau_diff_nonroot_grid_vs_neighbors_pm0003": ref,
            "min_eps_closed_roots_max": max(r[3] for r in rt), "min_eps_closed_nonroot_grid_min": min(r[3] for r in g),
            "initial_phase_test_M4": {"max_abs_dw": float(np.abs(W0 - W1).max()),
                                      "max_abs_dNeff_P": float(np.abs(neff(x0, y0) - neff(x1, y1)).max())}}


def e12_5():
    t, r = tr(R0)
    _, (xc, yc), _ = pool(R0, 1, 124, 1, 11, hist=True)
    _, (xp, yp), _ = pool(R0, 4, 124, 1, 11, hist=True)
    xf, yf = fresh(R0, 124)
    xc, yc, xp, yp = xc[:, 0], yc[:, 0], xp[:, 0], yp[:, 0]
    rows = [(k, float(yc[k].real), float(yc[k].imag), float(neff(xc[k], yc[k])), float(yf[k].real), float(yf[k].imag),
             float(neff(xf[k], yf[k])), float(yp[k].real), float(yp[k].imag), float(neff(xp[k], yp[k]))) for k in range(125)]
    write("e12_5_clock_width_R12423.csv", ["collision", "yP_re_closed", "yP_im_closed", "Neff_closed", "yP_re_fresh", "yP_im_fresh",
                                           "Neff_fresh", "yP_re_pool_M4", "yP_im_pool_M4", "Neff_pool_M4"], rows, "%.4g")
    mu = r - t
    return {"closed_max_dev_from_circle_center_half_radius_half": float(np.abs(np.abs(yc - 0.5) - 0.5).max()),
            "closed_yP_at_124": [float(yc[124].real), float(yc[124].imag)],
            "fresh_max_abs_absy2_minus_R_pow_k": float(np.abs(np.abs(yf) ** 2 - R0 ** np.arange(125)).max()),
            "arg_r_deg": math.degrees(np.angle(r)), "arg_mu_deg": math.degrees(np.angle(mu)),
            "arg_mu_minus_2arg_r_wrapped_deg": math.degrees((np.angle(mu) - 2 * np.angle(r) + math.pi) % (2 * math.pi) - math.pi),
            "Neff_P_at_124_closed_fresh_pool": [float(neff(xc[124], yc[124])), float(neff(xf[124], yf[124])), float(neff(xp[124], yp[124]))]}


def main():
    t0 = time.time()
    t, r = tr(R0)
    s = {"R": R0, "R_repr": repr(R0), "omega_over_2pi": 0.5 + math.asin(math.sqrt(R0)) / math.pi, "t": [t.real, t.imag], "r": [r.real, r.imag],
         "E12_2": e12_2(), "E12_3": e12_3(), "E12_4": e12_4(), "E12_5": e12_5()}
    e3, e4 = s["E12_3"], s["E12_4"]
    nb = e4["max_rel_tau_diff_root_vs_neighbors_pm0003"]
    ref = e4["max_rel_tau_diff_nonroot_grid_vs_neighbors_pm0003"]
    sel = ("根での差は根でない格子点で測った差の範囲（M=2 で {:.1e}、M=4 で {:.1e}）に収まり、整数比への選別は現れない".format(ref["M2"], ref["M4"])
           if all(nb[M] <= ref[M] for M in nb) else "根での差が根でない格子点の範囲を超える")
    s["E12_6"] = {
        "1": f"閉条件（M=1）は k=124 で局在の重みが {e3['M1']['w_at_124_mean']:.15f} に戻り、全パワーのずれは {e3['M1']['total_power_max_dev']:.1e}。"
             "単調な流出は起きない：この系では該当しない",
        "2": f"閉条件はどの R でも散逸がなく、有限位数根で厳密に戻り（ε≤{e4['min_eps_closed_roots_max']:.1e}）、それ以外は準周期。"
             f"開条件の寿命 τ_e の有限位数根と ±0.003 の隣接点との相対差は M=2 で {nb['M2']:.1e}、M=4 で {nb['M4']:.1e}：{sel}",
        "3": "自己振幅 y_P（内部時計の読み）と N_eff（局在の読み）は独立に記録した。相対運動を変えた比較には衝突の間の伝播の規則が要るため未判定",
        "4": "質量様応答はこの系では定義していない：判定対象外"}
    s["seeds"] = {"E12_3": "100+M（照合 7）", "E12_4": "4000+M（初期位相 77）", "E12_5": 11}
    s["versions"] = {"python": sys.version.split()[0], "numpy": np.__version__}
    (OUT / "summary_section12_R12423.json").write_text(json.dumps(s, ensure_ascii=False, indent=1), encoding="utf-8")
    print(json.dumps({k: s[k] for k in ("E12_2", "E12_3", "E12_4", "E12_5")}, ensure_ascii=False))
    print("elapsed_sec", round(time.time() - t0, 1))


if __name__ == "__main__":
    main()
