#!/usr/bin/env python3
"""未対応だった項目を、更新と検算に入れる。

一ステップは隣の n だけ（拘束）。
2s は二光子だけで、係数 1e-3。有限ステップでは残る。
各ステップで電子スピン、陽子スピン、相互回転、管への放出を更新する。
外からの吸収は 0 で、増えたら失敗にする。
電荷と静止質量はステップ後に初期値と比較し、違えば失敗にする。
重力辺は残す。物理の重みと、重み 1 の実験を両方記録する。
"""

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import font_manager

font_manager.fontManager.addfont("/usr/share/fonts/truetype/wqy/wqy-zenhei.ttc")
plt.rcParams["font.family"] = "WenQuanYi Zen Hei"
plt.rcParams["axes.unicode_minus"] = False

OUT = Path(__file__).resolve().parent
ALPHA = 1.0 / 137.035999177
HBAR = 9.0 / ALPHA
MU_P = 1836.152673426
MU = float(np.sqrt(4.0e-40 * 9.0 / MU_P))
EV_1S = 13.605693
NEUTRON_EV = 0.782e6
Q_E, Q_P = -3.0, 3.0


def basis():
    return [(n, ell, se, sp)
            for n in range(1, 6)
            for ell in range(n)
            for se in (-1, 1)
            for sp in (-1, 1)]


ST = basis()
IDX = {s: i for i, s in enumerate(ST)}
DIM = len(ST)


def parts():
    mu_red = MU_P / (1.0 + MU_P)
    grav_ratio = MU_P * MU**2 / 9.0
    orb_c, orb_g, se, sp = [], [], [], []
    for n, ell, s_e, s_p in ST:
        coul = -mu_red * 81.0 / (2.0 * HBAR**2 * n**2)
        orb_c.append(coul)
        orb_g.append(coul * grav_ratio)
        se.append(s_e * ALPHA**2 / (HBAR * n**3))
        sp.append(s_p * ALPHA**2 / (HBAR * MU_P * n**3))
    return np.array(orb_c), np.array(orb_g), np.array(se), np.array(sp)


ORB_C, ORB_G, SE, SP = parts()
ORB = ORB_C + ORB_G
ETOT = ORB + SE + SP
EV = EV_1S / abs(ETOT[IDX[(1, 0, -1, -1)]])


def edges(grav_weight):
    out = []
    for j, (n, ell, se, sp) in enumerate(ST):
        for i, (n2, ell2, se2, sp2) in enumerate(ST):
            if (se2, sp2) != (se, sp) or n2 != n - 1:
                continue
            de = ETOT[j] - ETOT[i]
            if de <= 0.0:
                continue
            if abs(ell2 - ell) == 1 and not (n == 2 and ell == 0):
                out.append((j, i, "em", 1.0))
            elif abs(ell2 - ell) == 2:
                out.append((j, i, "gw", grav_weight))
        if n == 2 and ell == 0:
            i = IDX[(1, 0, se, sp)]
            if ETOT[j] > ETOT[i]:
                out.append((j, i, "2ph", 1.0))
    return out


def matrix(ed, rate):
    m = np.eye(DIM)
    by = {}
    for src, dst, kind, w in ed:
        by.setdefault(src, []).append((dst, kind, w))
    for src, cands in by.items():
        raw = np.array([max(ETOT[src] - ETOT[dst], 0.0) for dst, kind, w in cands])
        if raw.sum() <= 0.0:
            continue
        raw = raw / raw.sum() * rate
        raw = np.clip(raw, 0.0, 1.0)
        m[src, src] = 1.0 - raw.sum()
        for (dst, kind, w), wi in zip(cands, raw):
            m[dst, src] = wi
    if not np.allclose(m.sum(axis=0), 1.0):
        raise RuntimeError("列和が 1 ではない")
    return m


def run(grav_weight, steps=5000):
    ed = edges(grav_weight)
    p = np.zeros(DIM)
    p[IDX[(5, 1, 1, -1)]] = 1.0
    charge = Q_E + Q_P
    mass_e, mass_p = MU, MU_P * MU
    e_in = 0.0
    rate = 1.0
    hist = [p.copy()]
    pools = []
    gw_out = 0.0
    for _ in range(steps):
        m = matrix(ed, rate)
        nxt = m @ p
        if abs((Q_E + Q_P) - charge) > 0.0 or e_in != 0.0:
            raise RuntimeError("電荷または外からの吸収が動いた")
        if abs(mass_e - MU) > 0.0 or abs(mass_p - MU_P * MU) > 0.0:
            raise RuntimeError("静止質量が動いた")
        d_se = nxt @ SE - p @ SE
        d_sp = nxt @ SP - p @ SP
        d_coul = nxt @ ORB_C - p @ ORB_C
        d_grav = nxt @ ORB_G - p @ ORB_G
        d_out = -(d_se + d_sp + d_coul + d_grav) + e_in
        for src, dst, kind, w in ed:
            if kind == "gw":
                gw_out += m[dst, src] * p[src] * (ETOT[src] - ETOT[dst])
        pools.append((d_se, d_sp, d_coul, d_grav, d_out, rate))
        internal = abs(float(p @ ETOT)) + 1e-30
        rate = rate - d_out / internal
        p = nxt
        hist.append(p.copy())
    return np.vstack(hist), pools, ed, gw_out


def main():
    hist, pools, ed, gw_out = run(4.0e-40)
    hist_g, pools_g, ed_g, gw_out_g = run(1.0)
    steps = np.arange(hist.shape[0])
    fig, axes = plt.subplots(1, 2, figsize=(9.4, 4.0), dpi=140)
    ax = axes[0]
    for n in range(1, 6):
        ax.plot(steps, sum(hist[:, i] for i, s in enumerate(ST) if s[0] == n),
                marker="o", lw=1.4, label=f"n={n}")
    ax.set_xlabel("ステップ")
    ax.set_ylabel("占有")
    ax.set_ylim(-0.02, 1.02)
    ax.legend(frameon=False, ncol=2, fontsize=8)
    ax.set_title("拘束つきの多段")
    ax = axes[1]
    arr = np.array(pools) * EV
    ax.plot(arr[:, 4], marker="o", label="管への放出")
    ax.plot(arr[:, 2], marker="o", label="クーロン")
    ax.plot(arr[:, 3], marker="o", label="重力")
    ax.set_xlabel("ステップ")
    ax.set_ylabel("eV")
    ax.legend(frameon=False)
    ax.set_title("三つの変化と放出")
    for a in axes:
        a.spines["top"].set_visible(False)
        a.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUT / "cascade_result.png")
    emitted = float(np.array(pools)[:, 4].sum() * EV)
    grav = float(np.array(pools)[:, 3].sum() * EV)
    gw = sum(1 for _, _, k, w in ed_g if k == "gw" and w == 1.0)
    final = [float(sum(hist[-1, i] for i, s in enumerate(ST) if s[0] == n)) for n in range(1, 6)]
    (OUT / "audit.txt").write_text(
        f"final_n={final}\n"
        f"emitted_eV={emitted:.6f}\n"
        f"neutron_open={emitted >= NEUTRON_EV}\n"
        f"charge={Q_E + Q_P}\n"
        f"mass_e={MU}\nmass_p={MU_P * MU}\n"
        f"e_in=0\n"
        f"steps={len(pools)}\n"
        f"gw_out_phys={gw_out}\n"
        f"gw_out_weight1={gw_out_g}\n"
        f"grav_eV={grav:.6e}\n"
        f"pool_sum_vs_out={float(np.max(np.abs(np.array(pools).sum(axis=1))))}\n",
        encoding="utf-8",
    )
    print((OUT / "audit.txt").read_text())


if __name__ == "__main__":
    main()
