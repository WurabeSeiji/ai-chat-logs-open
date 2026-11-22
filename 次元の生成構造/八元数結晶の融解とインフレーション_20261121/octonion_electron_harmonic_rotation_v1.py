#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
結合係数で回すとき、倍音の局在の形が保たれるか。仕様は思考実験2 の §24.8〜§24.9。

  (E) 候補 2 (全部の倍音を同じ角度 φ だけ回す) を繰り返したときの形 (§24.8)
  (F) m 番目の倍音を m×φ だけ回す形を繰り返したときの形と位置 (§24.9)

φ は結合係数の和 H = Σ_j k_ij の回転角 (実部と e7 成分の atan2)。k_ij は描き戻す前の状態から先に読む。
状態は倍音 1〜6 を等しい振幅で積んだ C 型 2 個。位置 = ベース波(m=1)の中心位相。
octonion_electron_coupling_first_map_v1.py と同じフォルダで実行する。格子点数は 1024 に固定。
"""
import sys, io, contextlib
import numpy as np
sys.path.insert(0, '.')
with contextlib.redirect_stdout(io.StringIO()):
    import octonion_electron_coupling_first_map_v1 as c
m = c.m; NG = c.NG; psi = c.psi; k = c.k
np.seterr(all='ignore')                      # (F) で大きさを掛ける場合、パワーが溢れる


def step_m(S, scale=True):
    """H を先に読み、m 番目の倍音を m×(H の回転角) だけ回す。scale=True なら大きさを |H| 倍する。C 型用。"""
    out = []
    for a in S:
        H = sum(k(a, b) for b in S)
        phi = np.arctan2(H[7], H[0]); r = np.hypot(H[0], H[7]) if scale else 1.0
        Z = np.fft.fft(a[:, 0] + 1j * a[:, 7]); w = np.fft.fftfreq(NG, 1 / NG)
        z = np.fft.ifft(Z * r * np.exp(1j * np.abs(w) * phi))
        X = np.zeros((NG, 8)); X[:, 0] = z.real; X[:, 7] = z.imag; out.append(X)
    return out


def read(X, s):
    """(ベース波の中心位相[度], 各倍音の位置[度], 振幅比の外れ, 半値幅[度], 実部の山の高さ)。s は巻きの符号。"""
    Z = np.fft.fft(X[:, 0] + 1j * X[:, 7]) / NG
    co = np.array([Z[(s * j) % NG] for j in range(1, 7)])
    pos = np.degrees(-s * np.angle(co)) / np.arange(1, 7)          # m 番目の倍音の位置 = 位相/m
    d = (X ** 2).sum(1)
    re = X[:, 0] / np.sqrt(m.power(X))
    return pos[0], pos, np.abs(np.abs(co / co[0]) - 1).max(), 360 * (d > d.max() / 2).sum() / NG, re.max()


def wrap(x):
    return (x + 180) % 360 - 180


print("===== (E) 全部の倍音を同じ角度だけ回す(候補 2)を繰り返す。C 型どうし、同じ巻き =====")
for Dg in (30, 90):
    S = [c.Cst(0), c.Cst(np.radians(Dg))]
    d0 = (S[0] ** 2).sum(1); d0 = d0 / d0.mean()
    print(" 初期の中心位相の差 %d°" % Dg)
    for it in range(4):
        p0, pos, ar, hw, pk = read(S[0], -1)
        d = (S[0] ** 2).sum(1)
        print("   %d 回: 密度の形の差 %.1e 半値幅 %.1f° 振幅比の外れ %.1e | 各倍音の位置 m=1..6:"
              % (it, np.abs(d / d.mean() - d0).max(), hw, ar), np.round(pos, 2), "| 実部の山の高さ %.2f" % pk)
        S = c.step(S)

print("\n===== (F) m 番目の倍音を m×φ だけ回す。12 回繰り返す =====")
for nm, sb in (("同じ巻き", -1), ("逆巻き", +1)):
    for sc in (True, False):
        print("==", nm, " 大きさ |H| を", ("掛ける" if sc else "掛けない(回転だけ)"))
        for Dg in (0, 1, 5, 30, 90, 150, 175, 179, 180):
            S = [c.Cst(0), c.Cst(np.radians(Dg), sb)]; tr = []; al = 0; ar = 0; hw = 0
            for it in range(13):
                ra, rb = read(S[0], -1), read(S[1], sb)
                if np.isfinite(m.power(S[0])) and m.power(S[0]) > 1e-250:
                    al = max(al, np.abs(wrap((ra[1] - ra[0]) * np.arange(1, 7))).max())   # 倍音の揃いの外れ[度]
                    ar = max(ar, ra[2]); hw = max(hw, abs(ra[3] - 53.8))
                tr.append((ra[0], rb[0], m.power(S[0])))
                S = step_m(S, sc)
            print(" 初期 %3d°: 差      " % Dg + " ".join("%+7.2f" % wrap(b - a) for a, b, _ in tr))
            print("           a の位置 " + " ".join("%+7.2f" % a for a, _, _ in tr))
            print("           パワー   " + " ".join("%.3g" % p for _, _, p in tr))
            print("           形: 倍音の揃いの外れ %.1e° 振幅比の外れ %.1e 半値幅の外れ %.2f°" % (al, ar, hw))
