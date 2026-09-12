#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""巡回二乗波の因数分解定理の検証（2026-09-13、木原レビュー反映）——論文4 §3–5 の一本化。

定理（巡回Fourier床 z_ij=r_d e^{i2π(i+j)/N}）:
  固定距離クラス C_d の二乗波 w_i=z_{i,i+d}^2=r_d^2 e^{i4π(2i+d)/N} は w_{i+1}=q w_i、
    q=e^{i8π/N},  L=ord(q)=N/gcd(N,4)。
  各クラス内の二乗波は正 L 角形の頂点を巡回。最小零閉塞ブロックは p=spf(L)（L の最小素因子）。
   N odd  : L=N       → p=spf(N)（3|N で3波、奇素数 N で N 波）
   N≡2(4): L=N/2     → p=spf(N/2)
   N≡4(8): L=N/4 奇  → p=spf(N/4)
   N≡0(8): L=N/4 偶  → p=2（2波＝±i 光子対）
  ゆえ 2波零閉塞が可能 ⟺ L 偶 ⟺ 8|N。
  N=8k の明示構成: L=N/4=2k, m=N/8, q^m=e^{iπ}=-1 ⇒ w_{i+m}=-w_i ⇒ z_{i+m}=±i z_i（同クラス・等振幅）。
    → グラフ探索なしで各クラスを i↔i+m で組み完全 matching が作れる。
  異クラス間: ±i 対は Δphase=π/2 が必要。位相は (2π/N)Z 上ゆえ π/2 は 4|N が必要、
    さらに位相 2π(2i+d)/N の 2i は偶なので d'-d≡N/4 (mod2) が必要。等振幅クラスは反射 r_d=r_{N/2-d}
    （数値確認）で d'-d=N/2-2d≡N/2 (mod2)。N≡4(8) は N/2 even かつ N/4 odd ⇒ 不可、8|N のみ可。
  ゆえ全空間で「厳密2波±i対が存在 ⟺ 8|N」。

読出しのみ・物理無変更・新規走行なし。r_d は既存の検証済み距離クラス振幅CSV（committed）から読む。
既存の網羅監査 even_N_summary.csv（全ペア）の candidate_pairs と独立に一致することを確認。
"""
import csv
import math
import os
from math import gcd

import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
AMP_CSV = os.path.join(BASE, '..', '高対称理論床_偶数系光子対精密解析_20260908', 'distance_class_amplitudes.csv')
SUMMARY_CSV = os.path.join(BASE, '..', '高対称理論床_偶数系光子対精密解析_20260908', 'even_N_summary.csv')
RESULTS = os.path.join(BASE, 'results')
os.makedirs(RESULTS, exist_ok=True)

log_lines = []


def log(*a):
    s = ' '.join(str(x) for x in a)
    print(s)
    log_lines.append(s)


def spf(n):
    if n < 2:
        return n
    d = 2
    while d * d <= n:
        if n % d == 0:
            return d
        d += 1
    return n


def order_root_of_unity(q, maxL):
    for L in range(1, maxL + 1):
        if abs(q**L - 1) < 1e-9:
            return L
    return None


# --- 検証済み r_d（距離クラス振幅）を読む（committed・偶数系）---
amp = {}   # amp[N][d] = r_d
with open(AMP_CSV, encoding='utf-8') as f:
    for row in csv.DictReader(f):
        amp.setdefault(int(row['N']), {})[int(row['distance_class'])] = float(row['amplitude'])

# --- 既存網羅監査の candidate_pairs（偶数系・全ペア）---
audit_pairs = {}
with open(SUMMARY_CSV, encoding='utf-8') as f:
    for row in csv.DictReader(f):
        audit_pairs[int(row['N'])] = int(row['candidate_pairs'])


def dist_class(i, j, N):
    dd = abs(i - j)
    return min(dd, N - dd)


def build_floor(N):
    """z_ij=r_d e^{i2π(i+j)/N}（辺順 triu_indices）。r_d は committed CSV。"""
    ea, eb = np.triu_indices(N, k=1)
    if N not in amp:
        return None  # 奇数系は本CSV(偶数系)に無い→位相のみ定理で扱う
    r = np.array([amp[N][dist_class(int(i), int(j), N)] for i, j in zip(ea, eb)])
    z = r * np.exp(1j * 2 * np.pi * (ea + eb) / N)
    return ea, eb, z


NS = list(range(3, 41))
log('# 巡回二乗波の因数分解定理 検証  日付 2026-09-13')
log(f'# 入力 r_d: {os.path.relpath(AMP_CSV, BASE)}（committed・検証済み）')
log('')
log('N  | gcd(N,4) L=N/gcd p=spf(L) | ord(q)=L? | 最小ブロック p (p角形和≈0) | 8|N⟺p=2 | '
    '±i対数(本検証) = 既往監査 | 反射 r_d=r_{N/2-d}')
rows = [['N', 'gcd_N_4', 'L', 'p_spf_L', 'ord_q_eq_L', 'pgon_sum_max', 'p_is_2_iff_8divN',
         'exact_pm_i_pairs', 'audit_candidate_pairs', 'match_audit', 'refl_amp_max_err',
         'n8k_construction_ok', 'n8k_matching_pairs', 'M_half']]
all_ok = True

for N in NS:
    g = gcd(N, 4)
    L = N // g
    p = spf(L)
    q = np.exp(1j * 8 * np.pi / N)
    ordq = order_root_of_unity(q, 4 * N)
    ordq_ok = (ordq == L)

    # 最小ブロック: q^{L/p} で生成する位数 p の部分群（正p角形）の和≈0
    # p<2（L=1, N=4 の退化）は多角形をなさず本機構の対象外（別機構で閉塞）。
    degenerate = (p < 2)
    if not degenerate:
        step = q**(L // p)
        pgon = np.array([step**t for t in range(p)])
        pgon_sum = float(abs(pgon.sum()))   # 正p角形の頂点和≈0
    else:
        pgon_sum = float('nan')

    p_is_2 = (p == 2)
    eight_div = (N % 8 == 0)
    p2_iff_8 = (p_is_2 == eight_div)

    # 反射振幅 r_d=r_{N/2-d}（偶数系のみ）
    refl_err = float('nan')
    exact_pairs = -1
    audit = audit_pairs.get(N, None)
    match_audit = None
    n8k_ok = None
    n8k_pairs = None
    Mhalf = None

    fl = build_floor(N)
    if fl is not None:
        ea, eb, z = fl
        M = len(z)
        Mhalf = M // 2
        # 反射対称 r_d=r_{N/2-d}
        if N % 2 == 0:
            half = N // 2
            errs = [abs(amp[N][d] - amp[N][half - d]) for d in range(1, half) if (half - d) in amp[N]]
            refl_err = max(errs) if errs else 0.0
        # 全ペア ±i 監査（eps=|z_a²+z_b²|/(|z_a|²+|z_b|²)）→ even_N_summary と一致確認
        z2 = z**2
        cnt = 0
        for a in range(M):
            for b in range(a + 1, M):
                denom = abs(z[a])**2 + abs(z[b])**2
                if abs(z2[a] + z2[b]) / denom < 1e-12:
                    cnt += 1
        exact_pairs = cnt
        match_audit = (audit is None) or (cnt == audit)

        # N=8k 明示構成: 同クラスで i→i+m（m=N/8）が z を +i 倍する（厳密）＋完全 matching
        if N % 8 == 0:
            m = N // 8
            # 辺を実頂点対（sorted）で表現。start i・距離 d の辺 = {i%N,(i+d)%N}。
            def edge_key(i, d):
                a, b = i % N, (i + d) % N
                return (min(a, b), max(a, b))
            # z を実辺（sorted 頂点対）で引く辞書
            zmap = {}
            for idx in range(M):
                i, j = int(ea[idx]), int(eb[idx])
                zmap[(min(i, j), max(i, j))] = z[idx]
            # 各クラス d の start-orbit（i→i+m, 位数 8）を4対に分割
            pairs = []
            used = set()
            rel_ok = True
            for d in range(1, N // 2 + 1):
                starts = list(range(N))
                # class N/2 は {i,i+N/2}={i+N/2,i} ゆえ start は 0..N/2-1
                if d == N // 2:
                    starts = list(range(N // 2))
                seen = set()
                for i0 in starts:
                    k0 = edge_key(i0, d)
                    if k0 in seen:
                        continue
                    # m-orbit を辿る
                    orbit = []
                    i = i0
                    while True:
                        k = edge_key(i, d)
                        if k in seen:
                            break
                        seen.add(k); orbit.append((i, k))
                        i = (i + m) % N
                    # 連続ペア (orbit[2t], orbit[2t+1]) を ±i として組む
                    for t in range(0, len(orbit) - 1, 2):
                        (ia, ka), (ib, kb) = orbit[t], orbit[t + 1]
                        za, zb = zmap[ka], zmap[kb]
                        if abs(abs(zb / za) - 1) > 1e-9 or abs(abs((zb / za).imag) - 1) > 1e-9:
                            rel_ok = False
                        if ka in used or kb in used:
                            rel_ok = False
                        used.add(ka); used.add(kb)
                        pairs.append((ka, kb))
            n8k_pairs = len(pairs)
            n8k_ok = rel_ok and (len(used) == M) and (n8k_pairs == Mhalf)

    ok = (ordq_ok and (degenerate or pgon_sum < 1e-9) and p2_iff_8
          and (match_audit is not False)
          and (math.isnan(refl_err) or refl_err < 1e-9)
          and (n8k_ok is not False))
    all_ok = all_ok and ok

    log(f"{N:>2} | gcd {g}  L={L:<3} p={p:<3} | {str(ordq_ok):>5} | pgon和 {pgon_sum:.1e} | "
        f"{str(p2_iff_8):>5} | 対数 {exact_pairs} = 監査 {audit} ({match_audit}) | "
        f"反射誤差 {refl_err:.1e}" + (f" | 8k構成 {n8k_ok} {n8k_pairs}/{Mhalf}" if N % 8 == 0 else "")
        + f" | {'OK' if ok else 'NG'}")
    rows.append([N, g, L, p, ordq_ok, f'{pgon_sum:.2e}', p2_iff_8, exact_pairs, audit,
                 match_audit, (f'{refl_err:.2e}' if not math.isnan(refl_err) else 'NA'),
                 n8k_ok, n8k_pairs, Mhalf])

with open(os.path.join(RESULTS, 'cyclic_factorization_summary.csv'), 'w', newline='', encoding='utf-8') as f:
    csv.writer(f).writerows(rows)

log('')
log('# 最小ブロックサイズ p=spf(L) の N mod 8 ごと（定理表）')
for r in ['odd', 'N≡2(4)', 'N≡4(8)', 'N≡0(8)']:
    ex = {'odd': [n for n in NS if n % 2],
          'N≡2(4)': [n for n in NS if n % 4 == 2],
          'N≡4(8)': [n for n in NS if n % 8 == 4],
          'N≡0(8)': [n for n in NS if n % 8 == 0]}[r]
    log(f"  {r:8}: " + ', '.join(f"N={n}:p={spf(n // gcd(n, 4))}" for n in ex[:6]) + (' ...' if len(ex) > 6 else ''))

log('')
log('=== 総合判定:', 'ALL OK' if all_ok else 'NG あり', '===')
with open(os.path.join(RESULTS, '実行ログ_20260913.log'), 'w', encoding='utf-8') as f:
    f.write('\n'.join(log_lines) + '\n')
import sys
sys.exit(0 if all_ok else 1)
