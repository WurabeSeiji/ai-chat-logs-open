#!/usr/bin/env python3
"""骨格写像の定常状態（関係の固定点）の数値実験 v1  2026-11-08

骨格（一行）:
    q_i' = normalize_i( sum_j q_j <q_j, q_i> / E_i ),   <a,b> = mean_psi( conj(a) b )
自己項 j = i を含む。交差項 Y_ij = <q_i, q_j> が相互作用であり、同時に読出しである。
右から関係を掛けるので、全体の左乗 q -> g q（無名性）に対して共変。
各波のパワー E_i は保存（弱観測：放射と吸収の釣り合い）。

実験
  E1 基準配置（4波）で複素（共通軸）と四元数（4軸）を反復し、グラム固有値・関係の変化・自己読みの虚部を記録
  E2 固定点が無名性の変換（共通左乗・共通並進）を除いて固定であることの確認
  E3 自分の大きさの1/2の規則から有効階数の閾値を決める
  E4 ランダム初期配置の集団で固定点の固有値・有効階数を集計
出力: results/skeleton_fixedpoint_v1.json, results/skeleton_fixedpoint_v1.npz
"""
import json, math, platform
from pathlib import Path
import numpy as np

NG = 721                       # 弦の標本数（奇数）
M0 = 12                        # 初期パルスの倍音数
TOL = 1e-13                    # 収束判定（関係の相対変化）
MAXSTEP = 2000
SEED = 20261108
OUT = Path(__file__).resolve().parent / "results"
psi = np.arange(NG) * 2 * np.pi / NG


def qmul(a, b):
    s = a[..., 0] * b[..., 0] - np.sum(a[..., 1:] * b[..., 1:], -1)
    v = a[..., :1] * b[..., 1:] + b[..., :1] * a[..., 1:] + np.cross(a[..., 1:], b[..., 1:])
    return np.concatenate([s[..., None], v], -1)


def qconj(q):
    return q * np.array([1., -1., -1., -1.])


def ip(a, b):
    """<a,b> = mean_psi conj(a) b（四元数）"""
    return np.mean(qmul(qconj(a), b), axis=0)


def pulse(theta, phi0, axis, M=M0):
    """中心 theta、内部位相 phi0、軸 axis の局在パルス（ハン窓で振幅を漸減）。E = 1 に規格化"""
    n = np.asarray(axis, float); n = n / np.linalg.norm(n)
    m = np.arange(1, M + 1); w = np.cos(0.5 * np.pi * m / (M + 1)) ** 2
    al = phi0 + np.outer(m, psi - theta)
    re = (w[:, None] * np.cos(al)).sum(0); im = (w[:, None] * np.sin(al)).sum(0)
    q = np.zeros((NG, 4)); q[:, 0] = re; q[:, 1:] = np.outer(im, n)
    return q / math.sqrt(np.mean(re ** 2 + im ** 2))


def relations(Q):
    N = len(Q)
    return np.array([[ip(Q[i], Q[j]) for j in range(N)] for i in range(N)])   # (N,N,4)


def step(Q):
    """骨格写像の一回。返り値: 次の状態、自己読み <q_i,q_i'>/E_i の四元数"""
    N = len(Q); Y = relations(Q); E = Y[np.arange(N), np.arange(N), 0]
    out, selfread = [], []
    for i in range(N):
        acc = np.zeros_like(Q[i])
        for j in range(N):
            acc += qmul(Q[j], np.broadcast_to(Y[j, i], Q[j].shape))     # q_j <q_j, q_i>
        acc /= E[i]
        selfread.append(ip(Q[i], acc) / E[i])
        acc *= math.sqrt(E[i] / ip(acc, acc)[0])                           # E_i を保存
        out.append(acc)
    return np.array(out), np.array(selfread)


def gram(Q):
    Y = relations(Q); E = Y[:, :, 0].diagonal()
    G = Y[:, :, 0] / np.sqrt(np.outer(E, E))
    return np.sort(np.linalg.eigvalsh(G)), Y


def qgram_eigs(Y):
    """四元数エルミートなグラム行列（規格化）の固有値。2N×2N 複素表現の固有値を一つおきに取る"""
    N = Y.shape[0]; E = Y[:, :, 0].diagonal(); Yn = Y / np.sqrt(np.outer(E, E))[:, :, None]
    A = np.zeros((2 * N, 2 * N), complex)
    for i in range(N):
        for j in range(N):
            a, b, c, d = Yn[i, j]
            A[2 * i:2 * i + 2, 2 * j:2 * j + 2] = [[a + 1j * b, c + 1j * d], [-c + 1j * d, a - 1j * b]]
    return np.sort(np.linalg.eigvalsh(A))[::2]


PEAK_FRAC = 0.05   # 共通密度の山の検出下限（最大値に対する比）


def profile(Q):
    """各波の密度 |q_i(psi)|^2 が共通か、共通密度の山の位置（度）と高さ比"""
    rho = (Q ** 2).sum(-1)
    same = float(np.max(np.abs(rho - rho[0])) / np.max(rho[0]))
    r = rho.sum(0)
    pk = [i for i in range(NG) if r[i] > r[i - 1] and r[i] >= r[(i + 1) % NG] and r[i] > PEAK_FRAC * r.max()]
    return same, [float(np.degrees(psi[i])) for i in pk], [float(r[i] / r.max()) for i in pk]


def iterate(Q, record=False):
    Y0 = relations(Q); hist = []
    for k in range(1, MAXSTEP + 1):
        Qn, sr = step(Q)
        Yn = relations(Qn)
        drel = float(np.linalg.norm(Yn - Y0) / np.linalg.norm(Y0))
        dabs = float(np.linalg.norm(Qn - Q) / np.linalg.norm(Q))
        if record:
            hist.append({"step": k, "eig": gram(Qn)[0].tolist(), "rel_change": drel, "abs_change": dabs,
                         "selfread_imag_max": float(np.max(np.abs(sr[:, 1:]))),
                         "selfread_real": sr[:, 0].tolist()})
        Q, Y0 = Qn, Yn
        if drel < TOL:
            return Q, k, hist
    return Q, MAXSTEP, hist


def half_size_threshold():
    """自分の大きさ = 自己相関 |<q(0), q(s)>| の半値全幅 σ。s = σ/2 だけ離れた同形二波の
    グラム最小固有値 1 - |overlap| を有効階数の閾値とする（自由パラメータなし）"""
    ref = pulse(0.0, 0.0, (1, 0, 0))
    ss = np.linspace(0, np.pi, 3601)
    ov = np.array([abs(ip(ref, pulse(s, 0.0, (1, 0, 0)))[0]) for s in ss])
    sigma = 2 * ss[np.argmax(ov < 0.5)]                     # 半値全幅
    o_half = abs(ip(ref, pulse(sigma / 2, 0.0, (1, 0, 0)))[0])
    return float(sigma), float(1 - o_half)


def main():
    OUT.mkdir(exist_ok=True)
    res = {"version": "skeleton_fixedpoint_v1", "date": "2026-11-08",
           "env": {"python": platform.python_version(), "numpy": np.__version__},
           "params": {"NG": NG, "M0": M0, "PEAK_FRAC": PEAK_FRAC, "TOL": TOL, "MAXSTEP": MAXSTEP, "seed": SEED}}
    npz = {}

    # E3 閾値
    sigma, eps = half_size_threshold()
    res["E3_half_size_rule"] = {"sigma_rad": sigma, "sigma_deg": math.degrees(sigma), "threshold": eps}

    # E1 基準配置
    th = np.deg2rad([0, 25, 60, 200]); ph = np.deg2rad([0, 70, 140, 210])
    cases = {"complex": [(1, 0, 0)] * 4,
             "quaternion": [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 1)]}
    res["E1"] = {"theta_deg": [0, 25, 60, 200], "phi0_deg": [0, 70, 140, 210]}
    for name, axes in cases.items():
        Q0 = np.array([pulse(th[k], ph[k], axes[k]) for k in range(4)])
        e0, Y0 = gram(Q0)
        Qf, n, hist = iterate(Q0, record=True)
        ef, Yf = gram(Qf)
        qe0, qef = qgram_eigs(Y0), qgram_eigs(Yf)
        same, pkpos, pkh = profile(Qf); _, pk0, pkh0 = profile(Q0)
        res["E1"][name] = {"axes": axes, "real_gram_eig_initial": e0.tolist(), "real_gram_eig_fixed": ef.tolist(),
                           "quat_gram_eig_initial": qe0.tolist(), "quat_gram_eig_fixed": qef.tolist(),
                           "quat_effective_rank_fixed": int(np.sum(qef > eps)),
                           "steps_to_fixed": n, "density_profile_difference_between_waves": same,
                           "peaks_initial_deg": pk0, "peaks_fixed_deg": pkpos, "peak_heights_fixed": pkh,
                           "right_factor_h_1j": [ (Yf[0, j] / Yf[0, 0, 0]).tolist() for j in range(4)],
                           "history": hist}
        npz[f"E1_{name}_Q0"] = Q0; npz[f"E1_{name}_Qfixed"] = Qf
        npz[f"E1_{name}_Y0"] = Y0; npz[f"E1_{name}_Yfixed"] = Yf

        # E2 無名性の変換を除いた固定
        g = np.array([math.cos(.7), .6 * math.sin(.7), .8 * math.sin(.7), 0.])
        Ql = np.array([qmul(np.broadcast_to(g, q.shape), q) for q in Qf])
        Qt = np.roll(Qf, 137, axis=1)
        out = {}
        for lab, Qx in (("left_multiplication", Ql), ("translation", Qt)):
            Qn, _ = step(Qx)
            out[lab] = {"relation_change_by_transform": float(np.linalg.norm(relations(Qx) - Yf) / np.linalg.norm(Yf)),
                        "still_fixed_rel_change": float(np.linalg.norm(relations(Qn) - relations(Qx)) / np.linalg.norm(Yf)),
                        "absolute_state_moved": float(np.linalg.norm(Qx - Qf) / np.linalg.norm(Qf))}
        res["E1"][name]["E2_anonymity"] = out

    # E4 集団
    rng = np.random.default_rng(SEED); NENS = 200; ens = {}
    for name in ("complex", "quaternion"):
        eigs, ranks, steps, qranks, npk, seps, sames = [], [], [], [], [], [], []
        for _ in range(NENS):
            th = rng.uniform(0, 2 * np.pi, 4); ph = rng.uniform(0, 2 * np.pi, 4)
            if name == "complex":
                axes = [(1, 0, 0)] * 4
            else:
                v = rng.normal(size=(4, 3)); axes = [tuple(x / np.linalg.norm(x)) for x in v]
            Q0 = np.array([pulse(th[k], ph[k], axes[k]) for k in range(4)])
            Qf, n, _ = iterate(Q0)
            ef, Yf = gram(Qf); qef = qgram_eigs(Yf)
            same, pkpos, _ = profile(Qf)
            eigs.append(ef); ranks.append(int(np.sum(ef > eps))); steps.append(n)
            qranks.append(int(np.sum(qef > eps))); npk.append(len(pkpos)); sames.append(same)
            if len(pkpos) == 2: seps.append(float(min(abs(pkpos[1] - pkpos[0]), 360 - abs(pkpos[1] - pkpos[0]))))
        eigs = np.array(eigs)
        ens[name] = {"n": NENS, "real_gram_rank_counts": {str(r): ranks.count(r) for r in range(1, 5)},
                     "quat_gram_rank_counts": {str(r): qranks.count(r) for r in range(1, 5)},
                     "peak_count_counts": {str(c): npk.count(c) for c in sorted(set(npk))},
                     "max_density_difference_between_waves": float(np.max(sames)),
                     "two_peak_separation_deg_min_max": [min(seps), max(seps)] if seps else None,
                     "steps_median": float(np.median(steps)), "steps_max": int(np.max(steps)),
                     "not_converged": int(np.sum(np.array(steps) >= MAXSTEP))}
        npz[f"E4_{name}_npeaks"] = np.array(npk); npz[f"E4_{name}_qranks"] = np.array(qranks)
        npz[f"E4_{name}_eigs"] = eigs; npz[f"E4_{name}_ranks"] = np.array(ranks); npz[f"E4_{name}_steps"] = np.array(steps)
    res["E4_ensemble"] = ens

    (OUT / "skeleton_fixedpoint_v1.json").write_text(json.dumps(res, ensure_ascii=False, indent=1), encoding="utf-8")
    np.savez_compressed(OUT / "skeleton_fixedpoint_v1.npz", **npz)
    summary = {k: res[k] for k in ("E3_half_size_rule", "E4_ensemble")}
    for name in cases:
        r = res["E1"][name]
        summary[name] = {kk: r[kk] for kk in ("real_gram_eig_initial", "real_gram_eig_fixed", "quat_gram_eig_fixed",
                         "quat_effective_rank_fixed", "steps_to_fixed", "density_profile_difference_between_waves",
                         "peaks_initial_deg", "peaks_fixed_deg", "peak_heights_fixed", "right_factor_h_1j", "E2_anonymity")}
        summary[name]["selfread_imag_max_all_steps"] = max(h["selfread_imag_max"] for h in r["history"])
    print(json.dumps(summary, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
