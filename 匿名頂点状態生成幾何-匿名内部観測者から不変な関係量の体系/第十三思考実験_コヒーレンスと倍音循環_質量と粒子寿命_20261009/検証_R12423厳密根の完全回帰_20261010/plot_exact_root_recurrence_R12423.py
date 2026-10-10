"""run_exact_root_recurrence_R12423.py の記録から図を作る（PNG と SVG）。

図（figures/）
  fig_A_steps_000-010      衝突 0〜10 の波形（原本の図 1 と同じ描き方: rho/max 対 chi/π、凡例に L と N_eff）
  fig_B_180deg_steps_052-072   交換位相 180°（k=62、完全入れ替わり）の前後 ±10 衝突
  fig_C_360deg_steps_114-134   交換位相 360°（k=124、完全回帰）の前後 ±10 衝突
  fig_D_trace_overview     全ステップの要約: 状態距離 eps_state・eps_swap（対数）、N_eff と閉形式、ピークの絶対高さ

波形の復元
  全ステップの chi 密度は、記録した座標 (alpha, beta, gamma, delta) と初期二状態の密度から
      rho_A,k = |alpha_k|^2 rho_a0 + |beta_k|^2 rho_b0,   rho_B,k = |gamma_k|^2 rho_a0 + |delta_k|^2 rho_b0
  で復元する（識別振動 m_A=1, m_B=2 により干渉項は厳密に 0）。results/rho_all_steps_R12423.npz が
  あれば、復元値との最大差を表示して照合する。

使い方
  python3 plot_exact_root_recurrence_R12423.py            # results/ を読んで figures/ に出力
"""
from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
RES = HERE / "results"
FIG = HERE / "figures"
FIG.mkdir(exist_ok=True)

plt.rcParams.update({"svg.hashsalt": "R12423",  # SVG の要素 ID を固定（再現性）
    "svg.fonttype": "none",          # 文字をパスに展開しない（SVG を小さく、テキストとして検索可能に）
    "path.simplify": True,
    "path.simplify_threshold": 1.0,
    "font.size": 8,
})

N_ROOT, M_ROOT = 124, 23
WINDOWS = {
    "fig_A_steps_000-010": (list(range(0, 11)), "R_{124,23}: collisions 0-10"),
    "fig_B_180deg_steps_052-072": (list(range(52, 73)), "R_{124,23}: around exchange phase 180° (k=62, full swap)"),
    "fig_C_360deg_steps_114-134": (list(range(114, 135)), "R_{124,23}: around exchange phase 360° (k=124, full recurrence)"),
}


def load():
    trace = list(csv.DictReader(open(RES / "trace_R12423.csv", encoding="utf-8")))
    coeff = list(csv.DictReader(open(RES / "coefficients_R12423.csv", encoding="utf-8")))
    basis = list(csv.DictReader(open(RES / "basis_density_R12423.csv", encoding="utf-8")))
    x = np.array([float(b["chi_over_pi"]) for b in basis])
    rho_a0 = np.array([float(b["rho_a0"]) for b in basis])
    rho_b0 = np.array([float(b["rho_b0"]) for b in basis])
    n = len(coeff)
    rho_a = np.empty((n, x.size))
    rho_b = np.empty((n, x.size))
    for c in coeff:
        k = int(c["collision"])
        al = complex(float(c["alpha_re"]), float(c["alpha_im"]))
        be = complex(float(c["beta_re"]), float(c["beta_im"]))
        ga = complex(float(c["gamma_re"]), float(c["gamma_im"]))
        de = complex(float(c["delta_re"]), float(c["delta_im"]))
        rho_a[k] = abs(al) ** 2 * rho_a0 + abs(be) ** 2 * rho_b0
        rho_b[k] = abs(ga) ** 2 * rho_a0 + abs(de) ** 2 * rho_b0
    npz = RES / "rho_all_steps_R12423.npz"
    if npz.exists():
        z = np.load(npz)
        print("reconstruction vs recorded rho: max |diff| =",
              max(float(np.max(np.abs(z["rho_A"] - rho_a))), float(np.max(np.abs(z["rho_B"] - rho_b)))))
    return trace, x, rho_a, rho_b


def panel_title(row):
    k = int(row["collision"])
    ph = float(row["exchange_phase_deg"])
    eps = float(row["eps_state"])
    eps_swap = float(row["eps_swap"])
    return f"k={k}   kω={ph:.2f}°   ε_state={eps:.2e}   ε_swap={eps_swap:.2e}"


def waveform_figure(name, steps, suptitle, trace, x, rho_a, rho_b):
    ncol = 3
    nrow = math.ceil(len(steps) / ncol)
    fig, axes = plt.subplots(nrow, ncol, figsize=(12, 2.3 * nrow + 0.6), sharex=True, sharey=True,
                             constrained_layout=True)
    axes = np.atleast_1d(axes).flatten()
    for ax, k in zip(axes, steps):
        row = trace[k]
        ra = rho_a[k] / np.max(rho_a[k])
        rb = rho_b[k] / np.max(rho_b[k])
        ax.plot(x, ra, lw=0.9, label=f"A L={float(row['L_A']):.3g}, N={float(row['N_eff_A']):.3g}")
        ax.plot(x, rb, lw=0.9, label=f"B L={float(row['L_B']):.3g}, N={float(row['N_eff_B']):.3g}")
        ax.set_title(panel_title(row), fontsize=7.5)
        ax.legend(fontsize=6.5, loc="upper right")
        ax.set_ylabel("rho_chi / max", fontsize=7)
    for ax in axes[len(steps):]:
        ax.set_visible(False)
    for ax in axes[max(0, len(steps) - ncol):len(steps)]:
        ax.set_xlabel("chi / pi")
    fig.suptitle(suptitle + f"   (R=cos²({M_ROOT}π/{N_ROOT})={math.cos(math.pi*M_ROOT/N_ROOT)**2:.15f}, "
                 "original update a'=normalize(r a + t b), b'=normalize(t a + r b))", fontsize=9)
    for ext in ("png", "svg"):
        fig.savefig(FIG / f"{name}.{ext}", dpi=110, metadata={"Date": None} if ext == "svg" else None)
    plt.close(fig)


def overview_title(trace):
    """全体図の見出し。要所の ε はデータから読む（固定文字列にしない）。"""
    e = {int(r["collision"]): r for r in trace}
    sw = ", ".join(f'{float(e[k]["eps_swap"]):.1e}' for k in (62, 186) if k in e)
    st = ", ".join(f'{float(e[k]["eps_state"]):.1e}' for k in (124, 248) if k in e)
    return f"R_{{124,23}}: full swap at k=62, 186 (ε_swap = {sw}) and full recurrence at k=124, 248 (ε_state = {st})"


def overview_figure(trace):
    ks = np.array([int(r["collision"]) for r in trace])
    eps = np.array([float(r["eps_state"]) for r in trace])
    eps_swap = np.array([float(r["eps_swap"]) for r in trace])
    na = np.array([float(r["N_eff_A"]) for r in trace])
    nb = np.array([float(r["N_eff_B"]) for r in trace])
    cf = np.array([float(r["closed_form_N_eff_A"]) for r in trace])
    pa = np.array([float(r["peak_raw_A"]) for r in trace])
    pb = np.array([float(r["peak_raw_B"]) for r in trace])
    fig, ax = plt.subplots(3, 1, figsize=(10, 9.5), sharex=True, constrained_layout=True)
    ax[0].semilogy(ks, np.maximum(eps, 1e-17), ".-", ms=3, lw=0.8, label="ε_state = ‖(a_k,b_k) − (a_0,b_0)‖/√2")
    ax[0].semilogy(ks, np.maximum(eps_swap, 1e-17), ".-", ms=3, lw=0.8, label="ε_swap = ‖(a_k,b_k) − (b_0,a_0)‖/√2")
    for k in (62, 124, 186, 248):
        ax[0].axvline(k, color="gray", lw=0.6, ls=":")
    ax[0].set_ylabel("state distance")
    ax[0].set_title(overview_title(trace))
    ax[0].legend(fontsize=8)
    ax[1].plot(ks, na, lw=0.9, label="N_eff A (original metrics)")
    ax[1].plot(ks, nb, lw=0.9, label="N_eff B (original metrics)")
    ax[1].plot(ks, cf, "k--", lw=0.7, label="closed form 1 + 31 sin²(kω/2)")
    ax[1].plot(ks, na + nb, color="gray", lw=0.7, label="N_eff A + B (= 33)")
    ax[1].set_ylabel("N_eff")
    ax[1].legend(fontsize=8)
    ax[2].plot(ks, pa, lw=0.9, label="A peak of sum_eta |psi|² (absolute)")
    ax[2].plot(ks, pb, lw=0.9, label="B peak of sum_eta |psi|² (absolute)")
    ax[2].set_ylabel("absolute peak height")
    ax[2].set_xlabel("collision k")
    ax[2].legend(fontsize=8)
    for ext in ("png", "svg"):
        fig.savefig(FIG / f"fig_D_trace_overview.{ext}", dpi=110, metadata={"Date": None} if ext == "svg" else None)
    plt.close(fig)


# ---------------------------------------------------------------------------
# 小さい SVG（*_compact.svg）: matplotlib を使わず同じ図を書く。
#   折れ線は Ramer–Douglas–Peucker で間引き（許容誤差 0.7 px）、座標は整数（単位 0.2 px）。
#   Google Drive / GitHub にテキストとして置くための版で、見た目は matplotlib 版と同じ内容。
# ---------------------------------------------------------------------------
A_COLOR, B_COLOR = "#1f77b4", "#ff7f0e"


def _rdp(pts, tol):
    """Ramer–Douglas–Peucker（反復版）。pts は [(x, y), ...]、tol は同じ単位の許容距離。"""
    keep = [False] * len(pts)
    keep[0] = keep[-1] = True
    stack = [(0, len(pts) - 1)]
    while stack:
        i, j = stack.pop()
        (x1, y1), (x2, y2) = pts[i], pts[j]
        dx, dy = x2 - x1, y2 - y1
        norm = math.hypot(dx, dy) or 1.0
        best, bi = -1.0, -1
        for k in range(i + 1, j):
            px, py = pts[k]
            d = abs(dy * (px - x1) - dx * (py - y1)) / norm
            if d > best:
                best, bi = d, k
        if best > tol and bi > 0:
            keep[bi] = True
            stack.append((i, bi))
            stack.append((bi, j))
    return [p for p, k in zip(pts, keep) if k]


def _path(pts):
    """整数座標の相対パス。"""
    out = []
    px, py = 0, 0
    for n, (x, y) in enumerate(pts):
        xi, yi = int(round(x)), int(round(y))
        out.append(("M%d %d" if n == 0 else "l%d %d") % ((xi, yi) if n == 0 else (xi - px, yi - py)))
        px, py = xi, yi
    return "".join(out)


def _esc(s):
    return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def compact_waveform_svg(name, steps, suptitle, trace, x, rho_a, rho_b):
    S = 5                                    # 1 px = 5 単位
    ncol, pw, ph, lm, tm, gx, gy = 3, 360, 110, 46, 30, 22, 44
    nrow = math.ceil(len(steps) / ncol)
    W = lm + ncol * (pw + gx) + 8
    H = tm + nrow * (ph + gy) + 10
    o = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W*S} {H*S}" width="{W}" height="{H}" '
         f'font-family="DejaVu Sans, Arial, sans-serif">',
         f'<rect width="{W*S}" height="{H*S}" fill="white"/>',
         f'<text x="{lm*S}" y="{14*S}" font-size="{9*S}">{_esc(suptitle)}   '
         f'(R=cos²({M_ROOT}π/{N_ROOT})={math.cos(math.pi*M_ROOT/N_ROOT)**2:.15f}, '
         f"original update a'=normalize(r a + t b), b'=normalize(t a + r b))</text>"]
    for n, k in enumerate(steps):
        r, c = divmod(n, ncol)
        x0, y0 = (lm + c * (pw + gx)) * S, (tm + r * (ph + gy)) * S
        row = trace[k]
        o.append(f'<rect x="{x0}" y="{y0}" width="{pw*S}" height="{ph*S}" fill="none" stroke="black" stroke-width="{S*0.6}"/>')
        for i, lab in enumerate(("0.0", "0.2", "0.4", "0.6", "0.8", "1.0")):
            y = y0 + ph * S - i * 0.2 * ph * S * 0.92 - 0.04 * ph * S
            o.append(f'<path d="M{x0 - 3*S} {y:.0f}h{3*S}" stroke="black" stroke-width="{S*0.6}"/>'
                     f'<text x="{x0 - 5*S}" y="{y + 2.5*S:.0f}" font-size="{6*S}" text-anchor="end">{lab}</text>')
        for xv in (-1, -0.5, 0, 0.5, 1):
            xx = x0 + (xv + 1) / 2 * pw * S
            o.append(f'<path d="M{xx:.0f} {y0 + ph*S}v{3*S}" stroke="black" stroke-width="{S*0.6}"/>')
            if r == nrow - 1 or n + ncol >= len(steps):
                o.append(f'<text x="{xx:.0f}" y="{y0 + ph*S + 11*S}" font-size="{6*S}" text-anchor="middle">{xv:g}</text>')
        o.append(f'<text x="{x0 - 34*S}" y="{y0 + ph*S/2:.0f}" font-size="{6*S}" text-anchor="middle" '
                 f'transform="rotate(-90 {x0 - 34*S} {y0 + ph*S/2:.0f})">rho_chi / max</text>')
        if r == nrow - 1 or n + ncol >= len(steps):
            o.append(f'<text x="{x0 + pw*S/2:.0f}" y="{y0 + ph*S + 20*S}" font-size="{7*S}" text-anchor="middle">chi / pi</text>')
        o.append(f'<text x="{x0 + pw*S/2:.0f}" y="{y0 - 4*S}" font-size="{6.5*S}" text-anchor="middle">{_esc(panel_title(row))}</text>')
        for rho, col, lab in ((rho_a[k], A_COLOR, f"A L={float(row['L_A']):.3g}, N={float(row['N_eff_A']):.3g}"),
                              (rho_b[k], B_COLOR, f"B L={float(row['L_B']):.3g}, N={float(row['N_eff_B']):.3g}")):
            v = rho / np.max(rho)
            pts = [(x0 + (xi + 1) / 2 * pw * S, y0 + ph * S - (0.04 + 0.92 * yi) * ph * S) for xi, yi in zip(x, v)]
            pts = _rdp(pts, 0.7 * S)
            o.append(f'<path d="{_path(pts)}" fill="none" stroke="{col}" stroke-width="{S*0.9}" stroke-linejoin="round"/>')
            dy = 9 * S if col == A_COLOR else 17 * S
            o.append(f'<path d="M{x0 + pw*S - 118*S} {y0 + dy - 2*S}h{10*S}" stroke="{col}" stroke-width="{S*0.9}"/>'
                     f'<text x="{x0 + pw*S - 105*S}" y="{y0 + dy}" font-size="{6*S}">{_esc(lab)}</text>')
    o.append("</svg>")
    (FIG / f"{name}_compact.svg").write_text("\n".join(o), encoding="utf-8")


def compact_overview_svg(trace):
    S = 5
    ks = [int(r["collision"]) for r in trace]
    series = {
        "eps": [("ε_state", A_COLOR, [max(float(r["eps_state"]), 1e-17) for r in trace]),
                ("ε_swap", B_COLOR, [max(float(r["eps_swap"]), 1e-17) for r in trace])],
        "neff": [("N_eff A", A_COLOR, [float(r["N_eff_A"]) for r in trace]),
                 ("N_eff B", B_COLOR, [float(r["N_eff_B"]) for r in trace]),
                 ("closed form 1+31sin²(kω/2)", "black", [float(r["closed_form_N_eff_A"]) for r in trace]),
                 ("A+B (=33)", "gray", [float(r["N_eff_sum"]) for r in trace])],
        "peak": [("A peak", A_COLOR, [float(r["peak_raw_A"]) for r in trace]),
                 ("B peak", B_COLOR, [float(r["peak_raw_B"]) for r in trace])],
    }
    pw, ph, lm, tm, gy = 900, 220, 70, 30, 40
    W, H = lm + pw + 20, tm + 3 * (ph + gy) + 10
    o = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W*S} {H*S}" width="{W}" height="{H}" '
         f'font-family="DejaVu Sans, Arial, sans-serif">', f'<rect width="{W*S}" height="{H*S}" fill="white"/>',
         f'<text x="{lm*S}" y="{14*S}" font-size="{9*S}">{overview_title(trace)}</text>']
    panels = [("eps", "state distance (log10)", (-17, 1)), ("neff", "N_eff", (0, 34)), ("peak", "absolute peak height", (0, 0.13))]
    for p, (key, ylab, (lo, hi)) in enumerate(panels):
        x0, y0 = lm * S, (tm + p * (ph + gy)) * S
        o.append(f'<rect x="{x0}" y="{y0}" width="{pw*S}" height="{ph*S}" fill="none" stroke="black" stroke-width="{S*0.6}"/>')
        ticks = [-16, -12, -8, -4, 0] if key == "eps" else ([0, 10, 20, 30] if key == "neff" else [0, 0.04, 0.08, 0.12])
        for t in ticks:
            y = y0 + ph * S * (1 - (t - lo) / (hi - lo))
            lab = f"1e{t}" if key == "eps" else f"{t:g}"
            o.append(f'<path d="M{x0 - 3*S} {y:.0f}h{3*S}" stroke="black" stroke-width="{S*0.6}"/>'
                     f'<text x="{x0 - 5*S}" y="{y + 2.5*S:.0f}" font-size="{6*S}" text-anchor="end">{lab}</text>')
        for kx in (0, 50, 100, 150, 200, 250):
            xx = x0 + kx / 250 * pw * S
            o.append(f'<path d="M{xx:.0f} {y0 + ph*S}v{3*S}" stroke="black" stroke-width="{S*0.6}"/>')
            if p == 2:
                o.append(f'<text x="{xx:.0f}" y="{y0 + ph*S + 11*S}" font-size="{6*S}" text-anchor="middle">{kx}</text>')
        if key == "eps":
            for kx in (62, 124, 186, 248):
                xx = x0 + kx / 250 * pw * S
                o.append(f'<path d="M{xx:.0f} {y0}v{ph*S}" stroke="gray" stroke-width="{S*0.5}" stroke-dasharray="{3*S} {3*S}"/>')
        o.append(f'<text x="{x0 - 50*S}" y="{y0 + ph*S/2:.0f}" font-size="{7*S}" text-anchor="middle" '
                 f'transform="rotate(-90 {x0 - 50*S} {y0 + ph*S/2:.0f})">{ylab}</text>')
        o.append(f'<rect x="{x0 + pw*S - 195*S}" y="{y0 + 2*S}" width="{190*S}" height="{(4 + 9*len(series[key]))*S}" fill="white" fill-opacity="0.85"/>')
        for n, (lab, col, ys) in enumerate(series[key]):
            vals = [math.log10(v) for v in ys] if key == "eps" else ys
            pts = _rdp([(x0 + k / 250 * pw * S, y0 + ph * S * (1 - (v - lo) / (hi - lo))) for k, v in zip(ks, vals)], 0.3 * S)
            dash = f' stroke-dasharray="{4*S} {3*S}"' if lab.startswith("closed") else ""
            o.append(f'<path d="{_path(pts)}" fill="none" stroke="{col}" stroke-width="{S*0.8}"{dash} stroke-linejoin="round"/>')
            o.append(f'<path d="M{x0 + pw*S - 190*S} {y0 + (9 + 9*n)*S - 2*S}h{10*S}" stroke="{col}" stroke-width="{S*0.9}"{dash}/>'
                     f'<text x="{x0 + pw*S - 177*S}" y="{y0 + (9 + 9*n)*S}" font-size="{6*S}">{_esc(lab)}</text>')
    o.append(f'<text x="{(lm + pw/2)*S:.0f}" y="{(H - 4)*S}" font-size="{7*S}" text-anchor="middle">collision k</text>')
    o.append("</svg>")
    (FIG / "fig_D_trace_overview_compact.svg").write_text("\n".join(o), encoding="utf-8")


def main():
    trace, x, rho_a, rho_b = load()
    for name, (steps, title) in WINDOWS.items():
        waveform_figure(name, steps, title, trace, x, rho_a, rho_b)
        compact_waveform_svg(name, steps, title, trace, x, rho_a, rho_b)
    overview_figure(trace)
    compact_overview_svg(trace)
    manifest = {name: {"collisions": steps, "files": [f"{name}.png", f"{name}.svg", f"{name}_compact.svg"]} for name, (steps, _) in WINDOWS.items()}
    manifest["fig_D_trace_overview"] = {"collisions": "all", "files": ["fig_D_trace_overview.png", "fig_D_trace_overview.svg", "fig_D_trace_overview_compact.svg"]}
    (FIG / "figures_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    for p in sorted(FIG.iterdir()):
        print(f"{p.name:40s} {p.stat().st_size:>9d} bytes")


if __name__ == "__main__":
    main()
