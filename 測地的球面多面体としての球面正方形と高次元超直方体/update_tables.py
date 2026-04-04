#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
update_tables.py  ─ マスター更新スクリプト
=============================================

このスクリプトを1回実行すると、以下をすべて自動で行う。

【1】数値計算
    内在的球面正方形・グノモン型の厳密値・級数近似・誤差を計算
【2】数表の自動差し込み
    MD ファイルの 10.2.1 / 10.2.2 節を計算結果で上書き
【3】実行スクリプトの自動コピー
    zenmi_report_spherical_polytope_complete_R2.py の内容を
    MD ファイルの 13.5 節コードブロックに自動同期
【4】PNG 図の再生成
    three_surfaces_new_R2{,_side2,_side3}.png を計算から再生成
【5】 n=4 計算実験（付録 A）の自動更新
    d_{4,m} 係数・厳密値（L=1のみ nquad）・充填率を計算し、
    MD の付録 A テーブルを自動更新

実行方法:
    python3 update_tables.py

依存ライブラリ:
    sympy, scipy, numpy, matplotlib
"""

import re
import math
from math import atan, sqrt, cos, pi, acos
from fractions import Fraction

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa

import sympy as sp
from scipy.integrate import nquad

# ══════════════════════════════════════════════════════════════
#  設定
# ══════════════════════════════════════════════════════════════

MD_PATH      = "zenmi_report_spherical_polytope_reproducible_full_R2.md"
PY_SRC_PATH  = "zenmi_report_spherical_polytope_complete_R2.py"
IMG_CONFIGS  = [
    (1.0, "three_surfaces_new_R2.png"),
    (2.0, "three_surfaces_new_R2_side2.png"),
    (3.0, "three_surfaces_new_R2_side3.png"),
]
R            = 2.0
LENGTHS      = [1.0, 2.0, 3.0]
N_INTR_TERMS = 10   # intrinsic Taylor series terms
N_PROJ_TERMS = 8    # gnomonic series terms
PLOT_GRID    = 50   # mesh resolution for 3D figures


# ══════════════════════════════════════════════════════════════
#  【計算関数】zenmi_report_spherical_polytope_complete_R2.py と同一
# ══════════════════════════════════════════════════════════════

def intrinsic_coeff_tuples(n_terms=N_INTR_TERMS):
    x = sp.symbols("x")
    f = 8 * sp.atan(1 / sp.sqrt(sp.cos(x))) - 2 * sp.pi
    ser = sp.series(f, x, 0, 2 * n_terms + 1).removeO()
    out = []
    for n in range(1, n_terms + 1):
        c = sp.simplify(ser.coeff(x, 2 * n))
        num, den = sp.fraction(c)
        out.append((int(num), int(den)))
    return out


def intrinsic_area_exact(a: float, R: float) -> float:
    x = a / R
    return (R ** 2) * (8.0 * atan(1.0 / sqrt(cos(x))) - 2.0 * pi)


def intrinsic_area_series(a: float, R: float, coeffs) -> float:
    x = a / R
    s = 0.0
    for n, (num, den) in enumerate(coeffs, start=1):
        s += (num / den) * (x ** (2 * n))
    return (R ** 2) * s


def weak_compositions(total: int, parts: int):
    if parts == 1:
        yield (total,)
        return
    for a in range(total + 1):
        for rest in weak_compositions(total - a, parts - 1):
            yield (a,) + rest


def projected_equal_edge_coeff_tuples(dim: int, n_terms=N_PROJ_TERMS):
    out = []
    for m in range(n_terms):
        total = Fraction(0, 1)
        for alpha in weak_compositions(m, dim):
            multinomial = math.factorial(m)
            for a in alpha:
                multinomial //= math.factorial(a)
            denom = 1
            for a in alpha:
                denom *= (2 * a + 1)
            total += Fraction(multinomial, denom)
        poch = Fraction(1, 1)
        for j in range(m):
            poch *= Fraction(dim + 1, 2) + j
        coeff = (Fraction((-1) ** m, 1) * poch * total
                 / Fraction(math.factorial(m) * (4 ** m), 1))
        out.append((coeff.numerator, coeff.denominator))
    return out


def projected_area_exact_equal(R: float, L: float):
    h = L / 2.0
    def integrand(y1, y2):
        return R ** 3 / (R ** 2 + y1 * y1 + y2 * y2) ** 1.5
    val, err = nquad(integrand, [[-h, h], [-h, h]])
    return val, err


def projected_area_series_equal(R: float, L: float, coeffs) -> float:
    x2 = (L * L) / (R * R)
    s = 0.0
    for m, (num, den) in enumerate(coeffs):
        s += (num / den) * (x2 ** m)
    return (L ** 2) * s


def rel_err(approx: float, exact: float) -> float:
    if exact == 0.0:
        return abs(approx - exact)
    return abs((approx - exact) / exact)


# ══════════════════════════════════════════════════════════════
#  【1】数値計算
# ══════════════════════════════════════════════════════════════

print("=" * 60)
print("【1】数値計算")
print("=" * 60)

print("  SymPy で係数を生成中...")
coeffs_intr = intrinsic_coeff_tuples(N_INTR_TERMS)
coeffs_proj = projected_equal_edge_coeff_tuples(2, N_PROJ_TERMS)

print("  厳密値・级数値・誤差を計算中...")
intrinsic_rows = []
for a in LENGTHS:
    exact  = intrinsic_area_exact(a, R)
    approx = intrinsic_area_series(a, R, coeffs_intr)
    ae     = abs(approx - exact)
    re_    = rel_err(approx, exact)
    intrinsic_rows.append((a, exact, approx, ae, re_))
    print(f"  内在的 a={a}: exact={exact:.12f}, series={approx:.12f}, err={ae:.3e}")

gnomonic_rows = []
for L in LENGTHS:
    exact, quad_err = projected_area_exact_equal(R, L)
    approx = projected_area_series_equal(R, L, coeffs_proj)
    ae     = abs(approx - exact)
    re_    = rel_err(approx, exact)
    gnomonic_rows.append((L, exact, approx, ae, re_, quad_err))
    print(f"  グノモン L={L}: exact={exact:.12f}, series={approx:.12f}, err={ae:.3e}")


# ══════════════════════════════════════════════════════════════
#  【2】MD 数表の自動差し込み
# ══════════════════════════════════════════════════════════════

print()
print("=" * 60)
print("【2】MD 数表の差し込み")
print("=" * 60)


def fmt_sci(val: float) -> str:
    """2.947e-12 → \\(2.95\\times 10^{-12}\\)"""
    s = f"{val:.2e}"
    mantissa, exp = s.split("e")
    exp_int = int(exp)
    return f"\\({mantissa}\\times 10^{{{exp_int}}}\\)"


def make_intrinsic_table(rows) -> str:
    lines = [
        "| 辺長 \\(a\\) |          厳密値 |       20 項級数 |               絶対誤差 |               相対誤差 |",
        "| ---------: | --------------: | --------------: | ---------------------: | ---------------------: |",
    ]
    for (a, exact, approx, ae, re_) in rows:
        lines.append(
            f"| {int(a):10d}"
            f" | {exact:15.12f}"
            f" | {approx:15.12f}"
            f" | {fmt_sci(ae):>22}"
            f" | {fmt_sci(re_):>22} |"
        )
    return "\n".join(lines)


def make_gnomonic_table(rows) -> str:
    lines = [
        "| 辺長 \\(L\\) |         厳密値 |       8 項級数 |                絶対誤差 |                相対誤差 |",
        "| ---------: | -------------: | -------------: | ----------------------: | ----------------------: |",
    ]
    for (L, exact, approx, ae, re_, _) in rows:
        lines.append(
            f"| {int(L):10d}"
            f" | {exact:14.12f}"
            f" | {approx:14.12f}"
            f" | {fmt_sci(ae):>23}"
            f" | {fmt_sci(re_):>23} |"
        )
    return "\n".join(lines)


new_intr_tbl = make_intrinsic_table(intrinsic_rows)
new_gnom_tbl = make_gnomonic_table(gnomonic_rows)

with open(MD_PATH, "r", encoding="utf-8") as f:
    content = f.read()

# ── 内在的テーブル置換
PAT_INTR = re.compile(
    r"(### 10\.2\.1[^\n]*\n\n)(\| 辺長 \\?\(a\\?\).*?)(?=\n\n)",
    re.DOTALL,
)
content, n1 = PAT_INTR.subn(lambda m: m.group(1) + new_intr_tbl, content)
print(f"  内在的テーブル: {n1} 箇所置換")

# ── グノモンテーブル置換
PAT_GNOM = re.compile(
    r"(### 10\.2\.2[^\n]*\n\n)(\| 辺長 \\?\(L\\?\).*?)(?=\n\n)",
    re.DOTALL,
)
content, n2 = PAT_GNOM.subn(lambda m: m.group(1) + new_gnom_tbl, content)
print(f"  グノモンテーブル: {n2} 箇所置換")


# ══════════════════════════════════════════════════════════════
#  【3】13.5 節コードブロックを py ファイルから自動同期
# ══════════════════════════════════════════════════════════════

print()
print("=" * 60)
print("【3】13.5 節スクリプト自動同期")
print("=" * 60)

with open(PY_SRC_PATH, "r", encoding="utf-8") as f:
    py_source = f.read()

# ## 13.5 節の ```python ... ``` を丸ごと置換
PAT_CODE = re.compile(
    r"(## 14\.5[^\n]*\n\n)(```python\n)(.*?)(```)",
    re.DOTALL,
)

def replace_code_block(m):
    return m.group(1) + "```python\n" + py_source + "\n```"

content, n3 = PAT_CODE.subn(replace_code_block, content)
print(f"  コードブロック: {n3} 箇所同期")

# MD を保存
with open(MD_PATH, "w", encoding="utf-8") as f:
    f.write(content)
print(f"  → {MD_PATH} を更新しました")


# ══════════════════════════════════════════════════════════════
#  【4】PNG 図の再生成
# ══════════════════════════════════════════════════════════════

print()
print("=" * 60)
print("【4】PNG 図の再生成")
print("=" * 60)

# ── 幾何ヘルパー関数 ──

def gnomonic_map(y1, y2, R):
    rho = np.sqrt(R ** 2 + y1 ** 2 + y2 ** 2)
    return R * y1 / rho, R * y2 / rho, R ** 2 / rho


def slerp(v0, v1, t):
    dot = np.clip(np.sum(v0 * v1, axis=0), -1.0, 1.0)
    omega = np.arccos(dot)
    sin_omega = np.sin(omega)
    mask = np.abs(sin_omega) < 1e-10
    c0 = np.where(mask, 1.0 - t, np.sin((1 - t) * omega) / np.where(mask, 1.0, sin_omega))
    c1 = np.where(mask, t,       np.sin(t * omega)       / np.where(mask, 1.0, sin_omega))
    return c0 * v0 + c1 * v1


def intrinsic_square_surface(a: float, R: float, N: int = PLOT_GRID):
    x = a / R
    cos_x = cos(x)
    if cos_x <= 0:
        return None, None, None
    rho = acos(sqrt(cos_x))
    phi_list = [np.pi / 4, 3 * np.pi / 4, 5 * np.pi / 4, 7 * np.pi / 4]
    verts = np.array([
        [np.sin(rho) * np.cos(phi), np.sin(rho) * np.sin(phi), np.cos(rho)]
        for phi in phi_list
    ])
    V0, V1, V2, V3 = [verts[i] for i in range(4)]
    s = np.linspace(0, 1, N)
    t = np.linspace(0, 1, N)
    S, T = np.meshgrid(s, t)
    def g(v): return v[:, None, None] * np.ones((3, N, N))
    bottom = slerp(g(V0), g(V1), S[None])
    top    = slerp(g(V3), g(V2), S[None])
    pts    = slerp(bottom, top, T[None]) * R
    return pts[0], pts[1], pts[2]


def sphere_wireframe(R: float, n_lat: int = 18, n_lon: int = 18):
    coords = []
    lat = np.linspace(0, np.pi / 2, n_lat)
    lon = np.linspace(0, 2 * np.pi, n_lon)
    for la in lat:
        coords.append((R * np.sin(la) * np.cos(lon),
                       R * np.sin(la) * np.sin(lon),
                       R * np.cos(la) * np.ones_like(lon)))
    lat2 = np.linspace(0, np.pi / 2, n_lat)
    for lo in lon:
        coords.append((R * np.sin(lat2) * np.cos(lo),
                       R * np.sin(lat2) * np.sin(lo),
                       R * np.cos(lat2)))
    return coords


def make_figure(L: float, R: float, out_path: str, N: int = PLOT_GRID):
    area_flat      = L ** 2
    area_intrinsic = intrinsic_area_exact(L, R)
    area_gnomonic, _ = projected_area_exact_equal(R, L)

    fig = plt.figure(figsize=(11, 9))
    ax  = fig.add_subplot(111, projection='3d')

    for (x, y, z) in sphere_wireframe(R):
        ax.plot(x, y, z, color='gray', linewidth=0.4, alpha=0.3)

    # 青: 平面正方形
    h_flat = L / 2
    yy = np.linspace(-h_flat, h_flat, N)
    Y1f, Y2f = np.meshgrid(yy, yy)
    ax.plot_surface(Y1f, Y2f, np.zeros_like(Y1f),
                    color='steelblue', alpha=0.45, linewidth=0.2, edgecolor='steelblue')
    corners = np.array([[-h_flat,-h_flat,0],[h_flat,-h_flat,0],
                         [h_flat,h_flat,0],[-h_flat,h_flat,0],[-h_flat,-h_flat,0]])
    ax.plot(corners[:,0], corners[:,1], corners[:,2], color='steelblue', linewidth=1.5)
    ax.text(0, 0, -0.05 * R,
            f'Blue flat square\nside={L}, area={area_flat:.6f}',
            color='steelblue', fontsize=9, ha='center', va='top')

    # 緑: グノモン型球面像
    yg = np.linspace(-h_flat, h_flat, N)
    Y1g, Y2g = np.meshgrid(yg, yg)
    Xg, Yg, Zg = gnomonic_map(Y1g, Y2g, R)
    ax.plot_surface(Xg, Yg, Zg, color='green', alpha=0.35, linewidth=0, edgecolor='none')
    for seg_y1, seg_y2 in [
        (np.linspace(-h_flat, h_flat, N), np.full(N, -h_flat)),
        (np.full(N, h_flat),              np.linspace(-h_flat, h_flat, N)),
        (np.linspace(h_flat, -h_flat, N), np.full(N, h_flat)),
        (np.full(N, -h_flat),             np.linspace(h_flat, -h_flat, N)),
    ]:
        ex, ey, ez = gnomonic_map(seg_y1, seg_y2, R)
        ax.plot(ex, ey, ez, color='darkgreen', linewidth=2.0)
    lx, ly, lz = gnomonic_map(-h_flat * 0.9, h_flat * 0.9, R)
    ax.text(float(lx), float(ly), float(lz),
            f'Green gnomonic image\narea={area_gnomonic:.6f}',
            color='darkgreen', fontsize=9, ha='center')

    # 赤: 内在的球面正方形
    Xr, Yr, Zr = intrinsic_square_surface(L, R, N=N)
    if Xr is not None:
        ax.plot_surface(Xr, Yr, Zr, color='red', alpha=0.40, linewidth=0, edgecolor='none')
        tN = N
        for edge in [(Xr[0,:], Yr[0,:], Zr[0,:]),
                     (Xr[-1,:],Yr[-1,:],Zr[-1,:]),
                     (Xr[:,0], Yr[:,0], Zr[:,0]),
                     (Xr[:,-1],Yr[:,-1],Zr[:,-1])]:
            ax.plot(edge[0], edge[1], edge[2], color='red', linewidth=2.5)
        ax.text(float(Xr[tN//2,-1])*1.05, float(Yr[tN//2,-1])*1.05, float(Zr[tN//2,-1]),
                f'Red intrinsic spherical square\narea={area_intrinsic:.6f}',
                color='red', fontsize=9, ha='left')

    lim = max(L / 2 + 0.3, R * 0.8)
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_zlim(-0.1, R + 0.1)
    ax.set_xlabel('x', labelpad=5); ax.set_ylabel('y', labelpad=5); ax.set_zlabel('z', labelpad=5)
    title = ('Blue Flat Square, Red Intrinsic Spherical Square, Green Gnomonic Image'
             if L == 1 else
             f'Comparison for Flat Square Side Length = {int(L)} (R = {int(R)})')
    ax.set_title(title, fontsize=11, pad=12)
    ax.view_init(elev=25, azim=-60)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {out_path}  (flat={area_flat:.4f}, intr={area_intrinsic:.6f}, gnom={area_gnomonic:.6f})")


for L, fname in IMG_CONFIGS:
    make_figure(L, R, fname)


# ════════════════════════════════════════════════════════════
#  「5」 n=4 計算実験（付録 A 自動更新）
# ════════════════════════════════════════════════════════════

print()
print("=" * 60)
print("「5」 n=4 計算実験（付録 A）の更新")
print("=" * 60)

_DIM4        = 4
_N4_TERMS    = 8

# ―― d_{4,m} 係数（厳密分数）
coeffs_n4 = projected_equal_edge_coeff_tuples(dim=_DIM4, n_terms=_N4_TERMS)
print(f"  d_{{4,m}} 係数: {coeffs_n4}")

# ―― L=1 のみ nquad で厳密値（収束域の深部）
print("  L=1 の nquad 実行中（小領域、数秒で完了予定）...", flush=True)
_h1 = 0.5  # L=1, h=L/2
def _integrand_n4(y4, y3, y2, y1, _h=_h1, _R=R):
    return _R**5 / (_R**2 + y1**2 + y2**2 + y3**2 + y4**2)**2.5
exact_n4_L1, qerr_n4_L1 = nquad(
    _integrand_n4,
    [[-_h1, _h1]] * 4,
    opts={'limit': 60, 'epsabs': 1e-8, 'epsrel': 1e-8}
)
print(f"  L=1: exact={exact_n4_L1:.15f}  (quad_err={qerr_n4_L1:.3e})")

# ―― L=1,2,3 の級数値
def _series_n4(L):
    x2 = L**2 / R**2
    return L**_DIM4 * sum(n / d * x2**m for m, (n, d) in enumerate(coeffs_n4))

series_n4_L1 = _series_n4(1.0)
series_n4_L2 = _series_n4(2.0)
series_n4_L3 = _series_n4(3.0)

rho_n4_L1   = exact_n4_L1 / 1.0**_DIM4
ae_n4_L1    = abs(series_n4_L1 - exact_n4_L1)
re_n4_L1    = ae_n4_L1 / exact_n4_L1

import math as _math
Lmax_n4 = 2 * R / _math.sqrt(_DIM4)
print(f"  収束上限: L < {Lmax_n4:.4f}  (=2R/√n)")
print(f"  L=1: series={series_n4_L1:.12f}, rel_err={re_n4_L1:.3e}, rho={rho_n4_L1:.8f}")
print(f"  L=2: series={series_n4_L2:.12f}（収束境界）")
print(f"  L=3: series={series_n4_L3:.12f}（発散参考値）")

# ―― 付録 A1 テーブル生成
def _frac_str(num, den):
    if den == 1:
        return f"${num}$"
    sign_str = "-" if num < 0 else ""
    return f"${sign_str}{abs(num)}/{den}$"

def make_appendix_a1_table(coeffs):
    rows = []
    rows.append("| \\(m\\) | 係数 \\(d_{4,m}\\) | 分子 | 分母 | 小数値 |")
    rows.append("| ----: | ----------------: | ---: | ---: | -----: |")
    for m, (num, den) in enumerate(coeffs):
        if den == 1:
            frac = f"\\({num}\\)"
        else:
            frac = f"\\({num}/{den}\\)"
        decimal = f"{num/den:+.8f}"
        rows.append(f"|{m:6d} | {frac} | {num} | {den} | {decimal} |")
    return "\n".join(rows)

# ―― 付録 A3 テーブル生成
def make_appendix_a3_table(exact, series, ae, re, rho):
    rows = [
        "| 項目 | 値 |",
        "| :--- | ---: |",
        f"| 厳密値 \\(V_4\\)（nquad） | {exact:.12f} |",
        f"| 8項級数値 | {series:.12f} |",
        f"| 絶対誤差 | \\({ae:.2e}\\) |",
        f"| 相対誤差 | \\({re:.2e}\\) |",
        f"| 平面超体積 \\(L^4\\) | 1.0 |",
        f"| 充填率 \\(\\rho_4\\) | {rho:.6f} |",
    ]
    return "\n".join(rows)

new_a1 = make_appendix_a1_table(coeffs_n4)
new_a3 = make_appendix_a3_table(exact_n4_L1, series_n4_L1, ae_n4_L1, re_n4_L1, rho_n4_L1)

# ―― MD の付録 A テーブルを置換
with open(MD_PATH, "r", encoding="utf-8") as f:
    content5 = f.read()

PAT_A1 = re.compile(
    r"(<!-- APPENDIX_A1_TABLE_START -->\n)(.+?)(\n<!-- APPENDIX_A1_TABLE_END -->)",
    re.DOTALL,
)
content5, na1 = PAT_A1.subn(lambda m: m.group(1) + new_a1 + m.group(3), content5)
print(f"  付録 A1 テーブル: {na1} 箇所置換")

PAT_A3 = re.compile(
    r"(<!-- APPENDIX_A3_TABLE_START -->\n)(.+?)(\n<!-- APPENDIX_A3_TABLE_END -->)",
    re.DOTALL,
)
content5, na3 = PAT_A3.subn(lambda m: m.group(1) + new_a3 + m.group(3), content5)
print(f"  付録 A3 テーブル: {na3} 箇所置換")

with open(MD_PATH, "w", encoding="utf-8") as f:
    f.write(content5)
print(f"  → {MD_PATH} (付録 A) を更新しました")

print()
print("=" * 60)
print("すべての更新が完了しました。")
print(f"  ・数表          : {MD_PATH} (10.2.1 / 10.2.2)")
print(f"  ・スクリプト同期: {MD_PATH} (14.5 節)")
print(f"  ・PNG 画像      : {', '.join(f for _,f in IMG_CONFIGS)}")
print(f"  ・付録 A        : {MD_PATH} (付録 A.1 / A.3)")
print("=" * 60)
