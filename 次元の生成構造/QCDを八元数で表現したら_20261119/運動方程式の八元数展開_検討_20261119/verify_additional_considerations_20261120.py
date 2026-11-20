#!/usr/bin/env python3
"""「運動方程式の八元数展開」追加考察 (2026-11-20) の検算。項目番号 [W1]〜[W5] は追加考察の本文に対応する。
実行: python3 verify_additional_considerations_20261120.py   (Python 3, NumPy のみ。乱数の種は固定)
"""
import itertools
import numpy as np

print("[W1] 全員と一度に回す場合と、対ごとに処理する場合のパワー (N = 4)")
rng = np.random.default_rng(3)
N = 4
z = rng.normal(size=N) + 1j * rng.normal(size=N)
th = {p: rng.uniform(0.2, 1.0) for p in itertools.combinations(range(N), 2)}
P = lambda v: float(np.sum(abs(v) ** 2))


def expm(M, terms=60):
    r = t = np.eye(len(M), dtype=complex)
    for n in range(1, terms):
        t = t @ M / n
        r = r + t
    return r


K = np.zeros((N, N))
for (i, j), t in th.items():
    K[i, j], K[j, i] = t, -t
U = expm(K.astype(complex))
print(f"   元のパワー                              : {P(z):.6f}")
print(f"   (a) 全員と一度に (一つの N×N の回転)    : {P(U @ z):.6f}")
print(f"       行ごとの |U_ij|^2 の和               : {np.round((abs(U) ** 2).sum(1), 12)}")
print(f"       対角 |U_ii|^2 (自分に残る分)         : {np.round(abs(np.diag(U)) ** 2, 4)}")
zb = z.copy()
for (i, j), t in th.items():
    c, s = np.cos(t), np.sin(t)
    zb[i] += (c * z[i] + s * z[j]) - z[i]
    zb[j] += (-s * z[i] + c * z[j]) - z[j]
print(f"   (b) 対ごとに別々に計算して足す           : {P(zb):.6f}")


def seq(order):
    v = z.copy()
    for (i, j) in order:
        t = th[(i, j)]
        c, s = np.cos(t), np.sin(t)
        a, b = v[i], v[j]
        v[i], v[j] = c * a + s * b, -s * a + c * b
    return v


o1 = list(th)
print(f"   (c) 対ごとに順番に施す                   : {P(seq(o1)):.6f}")
print(f"       順番を逆にしたときの結果の差         : {np.linalg.norm(seq(o1) - seq(o1[::-1])):.6f}")

print("\n[W2] 倍音を上に積む場合と下に積む場合の、山の数 (巻き数 12)")
psi = np.linspace(0, 360, 36000, endpoint=False)
r = np.radians(psi)
M = 12


def local_max(y):
    return np.where((y > np.roll(y, 1)) & (y >= np.roll(y, -1)))[0]


def top_peaks(y, frac=0.999):
    return int(np.sum(y[local_max(y)] > frac * y.max()))


def second(y):
    v = np.sort(y[local_max(y)])[::-1]
    return v[1] / v[0]


cases = [("巻き数 12 の波だけ", np.cos(M * r)),
         ("上に積む (12, 24, 36, 48 次)", sum(np.cos(k * M * r) for k in range(1, 5))),
         ("下に積む (12 次から 1 次まで全部)", sum(np.cos(m * r) for m in range(1, M + 1))),
         ("下に積む (12 の約数だけ: 1,2,3,4,6,12)", sum(np.cos(d * r) for d in (1, 2, 3, 4, 6, 12)))]
for name, y in cases:
    print(f"   {name:<36} 最大の山の数 = {top_peaks(y):2d}   2番目の山/最大の山 = {second(y):.3f}")

print("\n[W3] パルスの中心と幅 (中心位相 137°、1 次から M 次まで)")
psi = np.linspace(0, 360, 720000, endpoint=False)
r = np.radians(psi)


def pulse(Mx, theta):
    return sum(np.cos(m * (r - np.radians(theta))) for m in range(1, Mx + 1))


def width(y):
    i = int(np.argmax(y))
    n = len(y)
    a = b = i
    while y[a % n] > 0:
        a -= 1
    while y[b % n] > 0:
        b += 1
    return (b - a) * 360 / n


for Mx in (6, 12, 24, 48, 100):
    y = pulse(Mx, 137.0)
    print(f"   M = {Mx:3d}: 山の位置 = {psi[np.argmax(y)]:8.3f}°   山の幅 = {width(y):7.3f}°   幅 × M = {width(y) * Mx:6.1f}")

print("\n[W4] 自己相関の 180° での値 (0° を 1 とする、12 次まで)")
psi = np.linspace(0, 360, 3600, endpoint=False)
r = np.radians(psi)
for name, orders in (("奇数と偶数の両方", range(1, 13)), ("奇数だけ", range(1, 13, 2)), ("偶数だけ", range(2, 13, 2))):
    y = sum(np.cos(m * r) for m in orders)
    F = np.fft.rfft(y)
    c = np.fft.irfft(abs(F) ** 2, len(y))
    c = c / c[0]
    print(f"   {name:<12} 180° = {c[len(psi) // 2]:+.3f}   90° = {c[len(psi) // 4]:+.3f}")

print("\n[W5] 交差項の虚部 = 二つの状態が張る平行四辺形の面積")
LINES = [(1, 2, 3), (1, 4, 5), (1, 7, 6), (2, 4, 6), (2, 5, 7), (3, 4, 7), (3, 6, 5)]
f = np.zeros((8, 8, 8))
for (i, j, k) in LINES:
    for (a, b, c_), s in (((i, j, k), 1), ((j, k, i), 1), ((k, i, j), 1), ((j, i, k), -1), ((i, k, j), -1), ((k, j, i), -1)):
        f[a, b, c_] = s


def mul(x, y):
    w = np.zeros(8)
    w[0] = x[0] * y[0] - x[1:] @ y[1:]
    w[1:] = x[0] * y[1:] + y[0] * x[1:] + np.einsum("ijk,i,j->k", f[1:, 1:, 1:], x[1:], y[1:])
    return w


conj = lambda x: np.concatenate([[x[0]], -x[1:]])
rng = np.random.default_rng(5)
ok = []
for _ in range(1000):
    a, b = rng.normal(size=8), rng.normal(size=8)
    y = mul(conj(a), b)
    ang = np.arccos(a @ b / np.linalg.norm(a) / np.linalg.norm(b))
    ok.append(np.isclose(y[0], a @ b) and np.isclose(np.linalg.norm(y[1:]), np.linalg.norm(a) * np.linalg.norm(b) * np.sin(ang)))
print(f"   八元数 1000 組: Re(ā b) = |a||b|cos(角) かつ |Im(ā b)| = |a||b|sin(角) : {all(ok)}")
z1, z2 = complex(*rng.normal(size=2)), complex(*rng.normal(size=2))
print(f"   複素数: Im(z̄1 z2) = {(np.conj(z1) * z2).imag:.6f}   三角形 (原点, z1, z2) の符号付き面積の 2 倍 = {z1.real * z2.imag - z1.imag * z2.real:.6f}")
