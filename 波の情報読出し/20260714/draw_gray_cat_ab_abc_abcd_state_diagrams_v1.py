from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, Polygon, Wedge


BASE_DIR = Path(__file__).resolve().parent
OUT_DIR = BASE_DIR / "gray_cat_state_transition_figures_v1"


COLORS = {
    "ink": "#1f2933",
    "muted": "#667085",
    "grid": "#d8dee9",
    "white": "#ffffff",
    "black": "#16191f",
    "gray": "#9aa4b2",
    "blue": "#2676b8",
    "amber": "#c47b2c",
    "green": "#4b8b5b",
    "red": "#b94d4d",
    "panel": "#f7f9fc",
}


def setup_font() -> None:
    plt.rcParams["font.family"] = [
        "Hiragino Sans",
        "YuGothic",
        "Arial Unicode MS",
        "DejaVu Sans",
    ]
    plt.rcParams["axes.unicode_minus"] = False
    plt.rcParams["svg.fonttype"] = "path"


def add_label(ax, x: float, y: float, text: str, size: int = 12, weight: str = "normal", color: str = "ink", ha: str = "center") -> None:
    ax.text(x, y, text, ha=ha, va="center", fontsize=size, fontweight=weight, color=COLORS[color], linespacing=1.25)


def draw_card(ax, x: float, y: float, w: float, h: float, title: str, subtitle: str, edge: str = "grid") -> None:
    card = FancyBboxPatch(
        (x - w / 2, y - h / 2),
        w,
        h,
        boxstyle="round,pad=0.025,rounding_size=0.06",
        linewidth=1.2,
        edgecolor=COLORS[edge],
        facecolor=COLORS["panel"],
    )
    ax.add_patch(card)
    add_label(ax, x, y + h * 0.38, title, size=12, weight="bold")
    add_label(ax, x, y - h * 0.36, subtitle, size=9, color="muted")


def draw_cat(ax, x: float, y: float, r: float, state: str) -> None:
    ear_left = Polygon(
        [(x - r * 0.55, y + r * 0.55), (x - r * 0.25, y + r * 1.1), (x, y + r * 0.55)],
        closed=True,
        edgecolor=COLORS["ink"],
        linewidth=1.0,
        facecolor=COLORS["gray"] if state == "gray" else COLORS["white"],
        zorder=2,
    )
    ear_right = Polygon(
        [(x, y + r * 0.55), (x + r * 0.25, y + r * 1.1), (x + r * 0.55, y + r * 0.55)],
        closed=True,
        edgecolor=COLORS["ink"],
        linewidth=1.0,
        facecolor=COLORS["gray"] if state == "gray" else COLORS["black"] if state == "black" else COLORS["white"],
        zorder=2,
    )
    ax.add_patch(ear_left)
    ax.add_patch(ear_right)

    if state == "mix":
        ax.add_patch(Wedge((x, y), r, 90, 270, facecolor=COLORS["white"], edgecolor=COLORS["ink"], linewidth=1.3, zorder=3))
        ax.add_patch(Wedge((x, y), r, -90, 90, facecolor=COLORS["black"], edgecolor=COLORS["ink"], linewidth=1.3, zorder=3))
        ax.plot([x, x], [y - r, y + r], color=COLORS["ink"], lw=1.1, zorder=4)
    else:
        fill = COLORS[state]
        ax.add_patch(Circle((x, y), r, facecolor=fill, edgecolor=COLORS["ink"], linewidth=1.3, zorder=3))

    eye_color = COLORS["black"] if state in ("white", "gray", "mix") else COLORS["white"]
    ax.add_patch(Circle((x - r * 0.32, y + r * 0.1), r * 0.06, color=eye_color, zorder=5))
    ax.add_patch(Circle((x + r * 0.32, y + r * 0.1), r * 0.06, color=eye_color, zorder=5))
    ax.plot([x - r * 0.1, x, x + r * 0.1], [y - r * 0.18, y - r * 0.26, y - r * 0.18], color=eye_color, lw=1.0, zorder=5)


def draw_or_pair(ax, x: float, y: float, r: float) -> None:
    draw_cat(ax, x - r * 0.85, y, r * 0.72, "white")
    draw_cat(ax, x + r * 0.85, y, r * 0.72, "black")
    add_label(ax, x, y + r * 0.58, "or", size=9, weight="bold", color="muted")


def draw_observer(ax, x: float, y: float, label: str, strength: str, color: str) -> None:
    ax.add_patch(Circle((x, y), 0.16, facecolor="#ffffff", edgecolor=COLORS[color], linewidth=1.3, zorder=6))
    ax.add_patch(Circle((x, y), 0.055, facecolor=COLORS[color], edgecolor="none", zorder=7))
    add_label(ax, x, y - 0.28, label, size=10, weight="bold", color=color)
    add_label(ax, x, y - 0.45, strength, size=8, color="muted")


def draw_arrow(ax, x1: float, y1: float, x2: float, y2: float, color: str = "muted") -> None:
    arr = FancyArrowPatch(
        (x1, y1),
        (x2, y2),
        arrowstyle="-|>",
        mutation_scale=16,
        linewidth=1.4,
        color=COLORS[color],
    )
    ax.add_patch(arr)


def draw_scenario(ax, y: float, row_title: str, states: tuple[str, str, str], subtitles: tuple[str, str, str], row_color: str) -> None:
    x_positions = (2.1, 5.2, 8.3)
    for x, col_title, subtitle in zip(x_positions, ("AB", "ABC", "ABCD"), subtitles):
        draw_card(ax, x, y, 2.25, 1.85, col_title, subtitle, edge=row_color)

    for x, state in zip(x_positions, states):
        if state == "or":
            draw_or_pair(ax, x, y - 0.03, 0.48)
        else:
            draw_cat(ax, x, y - 0.03, 0.43, state)

    draw_observer(ax, x_positions[1] + 0.92, y + 0.52, "C", "弱読出し", "blue")
    draw_observer(ax, x_positions[2] + 0.92, y + 0.52, "D", "強観測", "red")
    draw_arrow(ax, 3.35, y, 3.95, y)
    draw_arrow(ax, 6.45, y, 7.05, y)
    add_label(ax, 0.6, y + 0.35, row_title, size=12, weight="bold", color=row_color, ha="left")


def draw_overview(path_stem: str, figsize: tuple[float, float] = (13.5, 8.8)) -> None:
    setup_font()
    fig, ax = plt.subplots(figsize=figsize)
    ax.set_xlim(0, 10.0)
    ax.set_ylim(0, 7.4)
    ax.axis("off")
    add_label(ax, 5.0, 7.05, "白猫・黒猫・灰色猫: AB -> ABC -> ABCD 状態遷移図", size=18, weight="bold")
    add_label(ax, 5.0, 6.65, "Cは弱く読む。Dは強く観測する。灰色準安定相と灰色固有相を分けて読む。", size=11, color="muted")

    draw_scenario(
        ax,
        5.55,
        "準安定混在",
        ("mix", "mix", "or"),
        ("白+黒を作る", "白+黒を壊さず読む", "白 or 黒へ選択"),
        "green",
    )
    draw_scenario(
        ax,
        3.45,
        "灰色固有相",
        ("gray", "gray", "gray"),
        ("灰色猫として成立", "灰色として読める", "灰色のまま残る"),
        "amber",
    )
    draw_scenario(
        ax,
        1.35,
        "Cが強すぎる場合",
        ("mix", "or", "or"),
        ("白+黒を作る", "Cで白 or 黒へ落ちる", "選択後を読む"),
        "red",
    )

    fig.savefig(OUT_DIR / f"{path_stem}.png", dpi=220, bbox_inches="tight")
    fig.savefig(OUT_DIR / f"{path_stem}.svg", bbox_inches="tight")
    plt.close(fig)


def draw_single_scenario(path_stem: str, title: str, states: tuple[str, str, str], subtitles: tuple[str, str, str], row_color: str) -> None:
    setup_font()
    fig, ax = plt.subplots(figsize=(12.0, 3.6))
    ax.set_xlim(0, 10.0)
    ax.set_ylim(0, 2.7)
    ax.axis("off")
    add_label(ax, 5.0, 2.42, title, size=17, weight="bold")
    draw_scenario(ax, 1.25, "", states, subtitles, row_color)
    fig.savefig(OUT_DIR / f"{path_stem}.png", dpi=220, bbox_inches="tight")
    fig.savefig(OUT_DIR / f"{path_stem}.svg", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    draw_overview("gray_cat_ab_abc_abcd_three_scenarios_v1")
    draw_single_scenario(
        "gray_cat_metastable_mix_to_white_or_black_v1",
        "準安定混在: 白+黒 => 白+黒 -> 白 or 黒",
        ("mix", "mix", "or"),
        ("白+黒を作る", "C弱読出しで保持", "D強観測で選択"),
        "green",
    )
    draw_single_scenario(
        "gray_cat_eigen_gray_to_gray_v1",
        "灰色固有相: 灰色 => 灰色 -> 灰色",
        ("gray", "gray", "gray"),
        ("灰色猫として成立", "C弱読出しで保持", "D強観測でも保持"),
        "amber",
    )
    draw_single_scenario(
        "gray_cat_strong_c_selects_before_d_v1",
        "Cが強すぎる場合: 白+黒 => 白 or 黒 -> 白 or 黒",
        ("mix", "or", "or"),
        ("白+黒を作る", "Cで選択が起こる", "Dは選択後を読む"),
        "red",
    )


if __name__ == "__main__":
    main()
