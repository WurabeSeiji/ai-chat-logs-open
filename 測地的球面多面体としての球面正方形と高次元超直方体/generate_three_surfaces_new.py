#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
generate_three_surfaces_new.py

zenmi_report_spherical_polytope_reproducible_full_R2.md に掲載されている
3つのPNGファイル(辺長1,2,3)を別ファイル名で再生成する。

出力ファイル名:
  three_surfaces_new_R2.png
  three_surfaces_new_R2_side2.png
  three_surfaces_new_R2_side3.png

R=2, 辺長 L=1,2,3 の場合について:
  青: 平面正方形 (z=0)
  赤: 内在的球面正方形
  緑: グノモン型球面像

実行方法:
  python3 generate_three_surfaces_new.py
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa
from math import atan, sqrt, cos, pi, acos
from scipy.integrate import nquad


# ─── 数値計算 ────────────────────────────────────────────────

def intrinsic_area_exact(a: float, R: float) -> float:
    """内在的球面正方形の厳密面積"""
    x = a / R
    return (R ** 2) * (8.0 * atan(1.0 / sqrt(cos(x))) - 2.0 * pi)


def gnomonic_area_exact(R: float, L: float) -> float:
    """グノモン型球面像の厳密面積"""
    h = L / 2.0

    def integrand(y1, y2):
        return R ** 3 / (R ** 2 + y1 * y1 + y2 * y2) ** 1.5

    val, _ = nquad(integrand, [[-h, h], [-h, h]])
    return val


# ─── 幾何変換 ────────────────────────────────────────────────

def gnomonic_map(y1, y2, R):
    """
    グノモン射影: 平面 z=R 上の点 (y1,y2) → 半径 R の球面
    F(y1,y2) = R*(y1, y2, R) / sqrt(R^2 + y1^2 + y2^2)
    """
    rho = np.sqrt(R ** 2 + y1 ** 2 + y2 ** 2)
    return R * y1 / rho, R * y2 / rho, R ** 2 / rho


def slerp(v0, v1, t):
    """球面線形補間 (shape: (3,N,...))"""
    dot = np.clip(np.sum(v0 * v1, axis=0), -1.0, 1.0)
    omega = np.arccos(dot)
    sin_omega = np.sin(omega)
    # sin_omega=0 のとき v0 をそのまま返す
    mask = np.abs(sin_omega) < 1e-10
    coeff0 = np.where(mask, 1.0 - t, np.sin((1 - t) * omega) / np.where(mask, 1.0, sin_omega))
    coeff1 = np.where(mask, t, np.sin(t * omega) / np.where(mask, 1.0, sin_omega))
    return coeff0 * v0 + coeff1 * v1


def intrinsic_square_surface(a: float, R: float, N: int = 40):
    """
    内在的球面正方形の3D座標グリッドを返す。

    測地線辺長 a の球面正方形の頂点極角: rho = arccos(sqrt(cos(a/R)))
    頂点は方位角 45°,135°,225°,315° に配置。
    内部は双線形 slerp で補完。
    """
    x = a / R
    cos_x = cos(x)
    if cos_x <= 0:
        return None, None, None

    rho = acos(sqrt(cos_x))          # 極角 (中心から頂点まで)

    # 頂点 (R 上の点): 方位角 45,135,225,315 度
    phi_list = [np.pi / 4, 3 * np.pi / 4, 5 * np.pi / 4, 7 * np.pi / 4]
    verts = np.array([
        [R * np.sin(rho) * np.cos(phi), R * np.sin(rho) * np.sin(phi), R * np.cos(rho)]
        for phi in phi_list
    ])  # shape (4,3): V0,V1,V2,V3 (反時計回り)

    # 単位ベクトル化
    V0 = verts[0] / R
    V1 = verts[1] / R
    V2 = verts[2] / R
    V3 = verts[3] / R

    s = np.linspace(0, 1, N)
    t = np.linspace(0, 1, N)
    S, T = np.meshgrid(s, t)  # (N,N)

    # slerp は axis=0 で 3 成分を一括処理するため shape を (3,N,N) に
    V0_g = V0[:, None, None] * np.ones((3, N, N))
    V1_g = V1[:, None, None] * np.ones((3, N, N))
    V2_g = V2[:, None, None] * np.ones((3, N, N))
    V3_g = V3[:, None, None] * np.ones((3, N, N))

    # bottom: slerp(V0, V1, s)  top: slerp(V3, V2, s)
    bottom = slerp(V0_g, V1_g, S[None])   # V0→V1
    top    = slerp(V3_g, V2_g, S[None])   # V3→V2

    # interior: slerp(bottom, top, t)
    pts = slerp(bottom, top, T[None]) * R  # 半径 R にスケール

    X = pts[0]
    Y = pts[1]
    Z = pts[2]
    return X, Y, Z


def sphere_wireframe(R: float, n_lat: int = 20, n_lon: int = 20):
    """球面ワイヤーフレーム用座標"""
    lat = np.linspace(0, np.pi / 2, n_lat)
    lon = np.linspace(0, 2 * np.pi, n_lon)
    coords = []
    # 緯度線
    for la in lat:
        x = R * np.sin(la) * np.cos(lon)
        y = R * np.sin(la) * np.sin(lon)
        z = R * np.cos(la) * np.ones_like(lon)
        coords.append((x, y, z))
    # 経度線
    lon2 = np.linspace(0, 2 * np.pi, n_lon)
    lat2 = np.linspace(0, np.pi / 2, n_lat)
    for lo in lon2:
        x = R * np.sin(lat2) * np.cos(lo)
        y = R * np.sin(lat2) * np.sin(lo)
        z = R * np.cos(lat2)
        coords.append((x, y, z))
    return coords


# ─── 描画 ────────────────────────────────────────────────────

def make_figure(L: float, R: float, out_path: str, N: int = 50):

    # --- 面積計算 ---
    area_flat     = L ** 2
    area_intrinsic = intrinsic_area_exact(L, R)
    area_gnomonic  = gnomonic_area_exact(R, L)

    fig = plt.figure(figsize=(11, 9))
    ax = fig.add_subplot(111, projection='3d')

    # ── 球面ワイヤーフレーム ────────────────────────────────
    for (x, y, z) in sphere_wireframe(R, n_lat=18, n_lon=18):
        ax.plot(x, y, z, color='gray', linewidth=0.4, alpha=0.3)

    # ── 青: 平面正方形 (z=0) ───────────────────────────────
    h_flat = L / 2
    yy = np.linspace(-h_flat, h_flat, N)
    Y1f, Y2f = np.meshgrid(yy, yy)
    Zf = np.zeros_like(Y1f)
    ax.plot_surface(Y1f, Y2f, Zf, color='steelblue', alpha=0.45,
                    linewidth=0.2, edgecolor='steelblue')
    # 枠線
    corners_b = np.array([[-h_flat, -h_flat, 0],
                           [ h_flat, -h_flat, 0],
                           [ h_flat,  h_flat, 0],
                           [-h_flat,  h_flat, 0],
                           [-h_flat, -h_flat, 0]])
    ax.plot(corners_b[:, 0], corners_b[:, 1], corners_b[:, 2],
            color='steelblue', linewidth=1.5)
    ax.text(0, 0, -0.05 * R,
            f'Blue flat square\nside={L}, area={area_flat:.6f}',
            color='steelblue', fontsize=9, ha='center', va='top')

    # ── 緑: グノモン型球面像 ──────────────────────────────
    yg = np.linspace(-h_flat, h_flat, N)
    Y1g, Y2g = np.meshgrid(yg, yg)
    Xg, Yg, Zg = gnomonic_map(Y1g, Y2g, R)
    ax.plot_surface(Xg, Yg, Zg, color='green', alpha=0.35,
                    linewidth=0, edgecolor='none')
    # 外縁線 (4辺)
    for seg_y1, seg_y2 in [
        (np.linspace(-h_flat, h_flat, N), np.full(N, -h_flat)),
        (np.full(N, h_flat),              np.linspace(-h_flat, h_flat, N)),
        (np.linspace(h_flat, -h_flat, N), np.full(N, h_flat)),
        (np.full(N, -h_flat),             np.linspace(h_flat, -h_flat, N)),
    ]:
        ex, ey, ez = gnomonic_map(seg_y1, seg_y2, R)
        ax.plot(ex, ey, ez, color='darkgreen', linewidth=2.0)
    # ラベル
    lx, ly, lz = gnomonic_map(-h_flat * 0.9, h_flat * 0.9, R)
    ax.text(float(lx), float(ly), float(lz),
            f'Green gnomonic image\narea={area_gnomonic:.6f}',
            color='darkgreen', fontsize=9, ha='center')

    # ── 赤: 内在的球面正方形 ─────────────────────────────
    Xr, Yr, Zr = intrinsic_square_surface(L, R, N=N)
    if Xr is not None:
        ax.plot_surface(Xr, Yr, Zr, color='red', alpha=0.40,
                        linewidth=0, edgecolor='none')
        # 外縁線 (4辺 = slerp で計算済みグリッドの境界)
        tN = N
        for edge in [
            (Xr[0, :],   Yr[0, :],   Zr[0, :]),    # bottom (t=0)
            (Xr[-1, :],  Yr[-1, :],  Zr[-1, :]),   # top    (t=1)
            (Xr[:, 0],   Yr[:, 0],   Zr[:, 0]),    # left   (s=0)
            (Xr[:, -1],  Yr[:, -1],  Zr[:, -1]),   # right  (s=1)
        ]:
            ax.plot(edge[0], edge[1], edge[2], color='red', linewidth=2.5)
        # ラベル (頂点付近)
        ax.text(float(Xr[tN // 2, -1]) * 1.05,
                float(Yr[tN // 2, -1]) * 1.05,
                float(Zr[tN // 2, -1]),
                f'Red intrinsic spherical square\narea={area_intrinsic:.6f}',
                color='red', fontsize=9, ha='left')

    # ── 軸・タイトル ──────────────────────────────────────
    lim = max(L / 2 + 0.3, R * 0.8)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_zlim(-0.1, R + 0.1)
    ax.set_xlabel('x', labelpad=5)
    ax.set_ylabel('y', labelpad=5)
    ax.set_zlabel('z', labelpad=5)

    if L == 1:
        ax.set_title(
            'Blue Flat Square, Red Intrinsic Spherical Square, Green Gnomonic Image',
            fontsize=11, pad=12)
    else:
        ax.set_title(
            f'Comparison for Flat Square Side Length = {int(L)} (R = {int(R)})',
            fontsize=11, pad=12)

    ax.view_init(elev=25, azim=-60)

    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'  Saved: {out_path}')
    print(f'    Blue  (flat)      area = {area_flat:.6f}')
    print(f'    Red   (intrinsic) area = {area_intrinsic:.6f}')
    print(f'    Green (gnomonic)  area = {area_gnomonic:.6f}')


# ─── メイン ──────────────────────────────────────────────────

if __name__ == '__main__':
    import matplotlib
    matplotlib.use('Agg')  # ヘッドレス環境用

    R = 2.0
    configs = [
        (1, 'three_surfaces_new_R2.png'),
        (2, 'three_surfaces_new_R2_side2.png'),
        (3, 'three_surfaces_new_R2_side3.png'),
    ]

    print('=== Generating 3 PNG figures ===')
    for L, fname in configs:
        print(f'\nL={L}:')
        make_figure(float(L), R, fname)

    print('\nAll done.')
