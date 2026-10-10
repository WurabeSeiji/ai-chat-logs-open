"""run_section12_R12423.py の結果（results/）から図 E〜H を作る（PNG は matplotlib、SVG は同じ図の仕様から直接書く小型版）。

  fig_E_state_tracking_R12423   §12.2 二乗閉塞 |Σψ²| とエネルギー候補 E2（閉条件、k=0..124）
  fig_F_open_vs_closed_R12423   §12.3 局在の重み w_loc=|y_P|²：閉（M=1）、プール M=2,4,16,64、全開放
  fig_G_sweep_R                 §12.4 寿命 τ_e と閉条件の最小状態距離を ω/2π に対して（印は有限位数根）
  fig_H_clock_width_R12423      §12.5 自己振幅 y_P(k) の軌跡（内部時計の読み）と N_eff（局在の読み）
"""
import csv
import html
import json
import math
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
RES, FIG = HERE / "results", HERE / "figures"
FIG.mkdir(exist_ok=True)
S = json.loads((RES / "summary_section12_R12423.json").read_text(encoding="utf-8"))
C = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#999999", "#000000"]


def tab(name):
    rows = list(csv.DictReader(open(RES / name, encoding="utf-8")))
    return {k: np.array([float(r[k]) for r in rows]) for k in rows[0] if k != "point"}, np.array([r.get("point") for r in rows])


def ser(x, y, c, l=None, w=0.9, d=None, m=False):
    return dict(x=np.asarray(x, float), y=np.asarray(y, float), c=c, l=l, w=w, d=d, m=m)


def lims(p):
    xs = np.concatenate([s["x"] for s in p["s"]])
    ys = np.concatenate([s["y"] for s in p["s"]])
    xl = p.get("xlim") or (xs.min(), xs.max())
    if p.get("ylim"):
        return xl, p["ylim"]
    if p.get("logy"):
        v = ys[ys > 0]
        return xl, (10 ** math.floor(math.log10(v.min())), 10 ** math.ceil(math.log10(v.max())))
    a, b = ys.min(), ys.max()
    return xl, (a - 0.05 * (b - a), b + 0.05 * (b - a))


def ticks(a, b, log):
    if log:
        e = list(range(math.ceil(math.log10(a) - 1e-9), math.floor(math.log10(b) + 1e-9) + 1))
        return [10.0 ** v for v in e[::max(1, math.ceil(len(e) / 6))]]
    raw = (b - a) / 5
    mag = 10 ** math.floor(math.log10(raw))
    st = next(mag * m for m in (1, 2, 5, 10) if mag * m >= raw)
    return [round(v, 10) for v in np.arange(math.ceil(a / st) * st, b + st * 1e-9, st)]


def lab(v, log):
    return "1e%d" % round(math.log10(v)) if log else "%g" % v


def png(name, sup, panels):
    fig, axs = plt.subplots(1, len(panels), figsize=(9, 3.6), constrained_layout=True)
    for ax, p in zip(axs, panels):
        for s in p["s"]:
            ax.plot(s["x"], s["y"], color=s["c"], lw=s["w"], ls="none" if s["m"] else (s["d"] or "-"), marker="o" if s["m"] else None, ms=3, label=s["l"])
        if p.get("logy"):
            ax.set_yscale("log")
        xl, yl = lims(p)
        ax.set(xlim=xl, ylim=yl, xlabel=p["x"], ylabel=p["y"], title=p["t"])
        if p.get("eq"):
            ax.set_aspect("equal")
        if any(s["l"] for s in p["s"]):
            ax.legend(fontsize=6.5, loc=p.get("loc", "best"))
    fig.suptitle(sup, fontsize=9)
    fig.savefig(FIG / f"{name}.png", dpi=130, metadata={"Software": None})
    plt.close(fig)


def svg(name, sup, panels, W=900, H=360):
    o = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}" font-family="DejaVu Sans,Arial,sans-serif" font-size="10">',
         f'<rect width="{W}" height="{H}" fill="#fff"/><text x="{W // 2}" y="16" text-anchor="middle" font-size="12">{html.escape(sup)}</text>']
    pw = W // len(panels)
    for i, p in enumerate(panels):
        (x0, x1), (y0, y1) = lims(p)
        L, T, R_, B = i * pw + 62, 46, (i + 1) * pw - 14, H - 40
        if p.get("eq"):
            sc = min((R_ - L) / (x1 - x0), (B - T) / (y1 - y0))
            R_, B = L + sc * (x1 - x0), T + sc * (y1 - y0)
        lg = p.get("logy")
        fy = (lambda v: math.log10(v)) if lg else (lambda v: v)
        X = lambda v: L + (v - x0) / (x1 - x0) * (R_ - L)
        Y = lambda v: B - (fy(min(max(v, y0), y1)) - fy(y0)) / (fy(y1) - fy(y0)) * (B - T)
        o.append(f'<rect x="{L}" y="{T}" width="{R_ - L:.0f}" height="{B - T:.0f}" fill="none" stroke="#000" stroke-width="0.8"/>')
        o.append(f'<text x="{(L + R_) / 2:.0f}" y="{T - 6}" text-anchor="middle">{html.escape(p["t"])}</text>')
        o.append(f'<text x="{(L + R_) / 2:.0f}" y="{H - 8}" text-anchor="middle">{html.escape(p["x"])}</text>')
        o.append(f'<text transform="translate({L - 46},{(T + B) / 2:.0f}) rotate(-90)" text-anchor="middle">{html.escape(p["y"])}</text>')
        for v in ticks(x0, x1, False):
            o.append(f'<path d="M{X(v):.0f} {B}v4" stroke="#000"/><text x="{X(v):.0f}" y="{B + 14}" text-anchor="middle">{lab(v, False)}</text>')
        for v in ticks(y0, y1, lg):
            o.append(f'<path d="M{L} {Y(v):.0f}h-4" stroke="#000"/><text x="{L - 6}" y="{Y(v) + 3:.0f}" text-anchor="end">{lab(v, lg)}</text>')
        for s in p["s"]:
            pts = [(round(X(a)), round(Y(b))) for a, b in zip(s["x"], s["y"]) if np.isfinite(b) and (b > 0 or not lg)]
            if s["m"]:
                d = "".join(f"M{a - 2} {b}a2 2 0 1 0 4 0a2 2 0 1 0-4 0" for a, b in pts)
                o.append(f'<path d="{d}" fill="{s["c"]}"/>')
                continue
            seg = []
            for a, b in zip(pts, pts[1:]):
                dx, dy = b[0] - a[0], b[1] - a[1]
                if seg and dy == 0 and seg[-1][1] == 0:
                    seg[-1][0] += dx
                elif dx or dy:
                    seg.append([dx, dy])
            d = f"M{pts[0][0]} {pts[0][1]}" + "".join(f"h{x}" if y == 0 else f"l{x} {y}" for x, y in seg)
            da = {":": ' stroke-dasharray="1.5 2.5"', "--": ' stroke-dasharray="5 3"'}.get(s["d"], "")
            o.append(f'<path d="{d}" fill="none" stroke="{s["c"]}" stroke-width="{s["w"]}"{da}/>')
        ent = [s for s in p["s"] if s["l"]]
        lx, ly = {"lower left": (L + 8, B - 8 - 12 * (len(ent) - 1)), "center right": (R_ - 150, (T + B) / 2 - 6 * len(ent))}.get(p.get("loc"), (R_ - 150, T + 12))
        for j, s in enumerate(ent):
            o.append(f'<path d="M{lx:.0f} {ly + 12 * j - 3:.0f}h14" stroke="{s["c"]}" stroke-width="2"/>'
                     f'<text x="{lx + 18:.0f}" y="{ly + 12 * j:.0f}" font-size="9">{html.escape(s["l"])}</text>')
    (FIG / f"{name}.svg").write_text("\n".join(o) + "\n</svg>\n", encoding="utf-8")


def make(name, sup, panels):
    png(name, sup, panels)
    svg(name, sup, panels)


def main():
    d, _ = tab("e12_2_state_tracking_R12423.csv")
    k = d["collision"]
    mx = ", ".join("%.1e" % v for v in S["E12_2"]["max_abs_sum_psi2_A_B_total"])
    make("fig_E_state_tracking_R12423", "§12.2  R = cos²(23π/124), closed two-wave exchange (pure U_R)", [
        dict(x="collision k", y="|Σψ²|  (A+B)", t=f"square closure; max A, B, A+B = {mx}", logy=True, ylim=(1e-18, 1),
             s=[ser(k, np.maximum(d["abs_sum_psi2_total"], 1e-18), C[0], w=0.8)]),
        dict(x="collision k", y="E2 = Σ n² w_n", t="energy candidate E2 (E1 = N_eff); A+B conserved", loc="center right",
             s=[ser(k, d["E2_A"], C[0], "E2 of A", 0.8), ser(k, d["E2_B"], C[1], "E2 of B", 0.8), ser(k, d["E2_A"] + d["E2_B"], C[6], "A+B", 0.8, "--")])])
    d, _ = tab("e12_3_open_vs_closed_R12423.csv")
    k = d["collision"]
    a, n, f = k <= 130, ((k <= 130) & (k % 2 == 0)) | (k > 130), d["fresh"] >= 1e-4
    pl = [ser(k[n], d[f"pool_M{M}"][n], C[i], f"pool M={M}") for i, M in enumerate((2, 4, 16, 64))]
    pl += [ser([0, 300], [1 / (M + 1)] * 2, C[i], None, 0.7, ":") for i, M in enumerate((2, 4, 16, 64))]
    make("fig_F_open_vs_closed_R12423", "§12.3  closed vs open: pool of M background channels, partner drawn at random each collision", [
        dict(x="collision k", y="w_loc = |y_P|²  (mean of 64 trials)", t="localized weight remaining in P",
             s=[ser(k[a], d["closed_M1"][a], C[5], "closed M=1", 0.6), ser(k[a], d["pool_M2"][a], C[0], "pool M=2"),
                ser(k[a], d["pool_M4"][a], C[1], "pool M=4"), ser(k[a], d["fresh"][a], C[4], "fresh (M=∞)")]),
        dict(x="collision k", y="w_loc", t="dotted: 1/(M+1)", logy=True, ylim=(1e-4, 1.2), xlim=(0, 300),
             s=pl + [ser(k[f], d["fresh"][f], C[4], "fresh: R^k")])])
    d, pt = tab("e12_4_sweep_R.csv")
    g, rt, w = pt == "grid", pt != "grid", d["omega_over_2pi"]
    e = np.maximum(d["min_eps_closed_k10_1000"], 1e-16)
    make("fig_G_sweep_R", "§12.4  sweep of R in 0.01–0.99 (not restricted to rationals); markers: finite-order roots n ≤ 8 and 124", [
        dict(x="ω/2π  (exchange phase per collision)", y="τ_e  (collisions)", t="lifetime τ_e (open pools)", loc="upper left",
             s=[ser(w[g], d["tau_e_M2"][g], C[0], "pool M=2"), ser(w[rt], d["tau_e_M2"][rt], C[0], m=True),
                ser(w[g], d["tau_e_M4"][g], C[1], "pool M=4"), ser(w[rt], d["tau_e_M4"][rt], C[1], m=True)]),
        dict(x="ω/2π", y="min ε_state  (10 ≤ k ≤ 1000)", t="closed M=1: closest return of the full state", logy=True, ylim=(1e-16, 1),
             s=[ser(w[g], e[g], C[0], "R grid", 0.8), ser(w[rt], e[rt], C[3], m=True)])])
    d, _ = tab("e12_5_clock_width_R12423.csv")
    k = d["collision"]
    cs = (("closed", "closed M=1", 0.5, C[0]), ("fresh", "fresh (M=∞)", 0.9, C[1]), ("pool_M4", "pool M=4", 0.7, C[2]))
    th = np.linspace(0, 2 * np.pi, 61)
    make("fig_H_clock_width_R12423", "§12.5  internal-clock readout and localization readout, recorded independently", [
        dict(x="Re y_P", y="Im y_P", t="y_P(k) = <b0, P_k>, k = 0..124 (reference F fixed)", eq=True, xlim=(-0.75, 1.15), ylim=(-0.7, 0.6),
             loc="lower left", s=[ser(0.5 + 0.5 * np.cos(th), 0.5 * np.sin(th), C[6], None, 0.6, ":")] +
             [ser(d[f"yP_re_{c}"], d[f"yP_im_{c}"], col, l, lw) for c, l, lw, col in cs]),
        dict(x="collision k", y="N_eff of P", t="localization readout", loc="center right",
             s=[ser(k, d[f"Neff_{c}"], col, l, 0.8) for c, l, lw, col in cs])])
    for p in sorted(FIG.iterdir()):
        print(f"{p.name:40s} {p.stat().st_size:>8d}")


if __name__ == "__main__":
    main()
