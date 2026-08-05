#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""グラフ構造の説明図: ノード（体）・辺（関係波）・隣接（体の共有）・ノード帳簿

fig_graph_N3_v1.png : N=3, M=3 の3パネル（全体／隣接／ノード帳簿）
fig_graph_N6_v1.png : N=6, M=15 の3パネル（同）
"""
import itertools

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams["font.family"] = ["Hiragino Sans", "Arial Unicode MS", "sans-serif"]
plt.rcParams["axes.unicode_minus"] = False

RED, BLUE, GRAY, GREEN = "#d62728", "#1f77b4", "#bbbbbb", "#2ca02c"


def node_pos(n, r=1.0):
    th = np.pi / 2 + 2 * np.pi * np.arange(n) / n
    return {i + 1: (r * np.cos(t), r * np.sin(t)) for i, t in enumerate(th)}


def draw_graph(ax, n, focus_edge=None, ledger_node=None, label_edges=True):
    pos = node_pos(n)
    edges = list(itertools.combinations(range(1, n + 1), 2))
    for (u, v) in edges:
        x1, y1 = pos[u]; x2, y2 = pos[v]
        if focus_edge is not None:
            if (u, v) == focus_edge:
                c, lw, z = RED, 3.2, 5
            elif u in focus_edge or v in focus_edge:
                c, lw, z = BLUE, 2.0, 4
            else:
                c, lw, z = GRAY, 1.0, 2
        elif ledger_node is not None:
            if ledger_node in (u, v):
                c, lw, z = GREEN, 2.6, 4
            else:
                c, lw, z = GRAY, 1.0, 2
        else:
            c, lw, z = "#555555", 1.6, 2
        ax.plot([x1, x2], [y1, y2], color=c, lw=lw, zorder=z,
                linestyle="--" if (focus_edge and c == GRAY) else "-")
        if label_edges:
            mx, my = 0.5 * (x1 + x2), 0.5 * (y1 + y2)
            off = 0.09 * np.hypot(mx, my) if (mx, my) != (0, 0) else 0
            ax.annotate(f"z{u}{v}", (mx, my), fontsize=8 if n > 3 else 11,
                        color=c if c != GRAY else "#999999",
                        ha="center", va="center", zorder=6,
                        bbox=dict(boxstyle="round,pad=0.12", fc="white",
                                  ec="none", alpha=0.75))
    for i, (x, y) in pos.items():
        fc = "#ffd54f" if (ledger_node == i or (focus_edge and i in focus_edge)) else "#eeeeee"
        ax.add_patch(plt.Circle((x, y), 0.13, fc=fc, ec="black", lw=1.4, zorder=8))
        ax.annotate(str(i), pos[i], ha="center", va="center",
                    fontsize=12, weight="bold", zorder=9)
    ax.set_xlim(-1.45, 1.45); ax.set_ylim(-1.45, 1.45)
    ax.set_aspect("equal"); ax.axis("off")


def make_fig(n, fname):
    m = n * (n - 1) // 2
    nb = 2 * (n - 2)
    other = m - 1 - nb
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.6))

    draw_graph(axes[0], n)
    axes[0].set_title(f"(a) 全体像  N={n}体（ノード）\n関係波（辺）M={n}·{n-1}/2={m}本が力学変数",
                      fontsize=11)

    draw_graph(axes[1], n, focus_edge=(1, 2))
    axes[1].set_title(f"(b) 隣接＝体の共有\n赤=z12  青=隣接{nb}本（体1か体2を共有）  "
                      f"灰破線=非隣接{other}本", fontsize=11)

    draw_graph(axes[2], n, ledger_node=1)
    ledger_terms = " + ".join(f"|z1{v}|²" for v in range(2, min(n + 1, 5)))
    if n > 4:
        ledger_terms += " + …"
    axes[2].set_title(f"(c) ノード帳簿（体1）\nA₁ = {ledger_terms}（緑の{n-1}本の合計）",
                      fontsize=11)

    fig.suptitle(f"N={n}体の関係波グラフ：体=ノード（変数なし・集約点）、関係波=辺（変数の実体）",
                 fontsize=12, y=1.02)
    fig.tight_layout()
    fig.savefig(fname, dpi=160, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {fname}  (M={m}, 隣接={nb}, 非隣接={other}, 検算 1+{nb}+{other}={1+nb+other})")


make_fig(3, "fig_graph_N3_v1.png")
make_fig(6, "fig_graph_N6_v1.png")
