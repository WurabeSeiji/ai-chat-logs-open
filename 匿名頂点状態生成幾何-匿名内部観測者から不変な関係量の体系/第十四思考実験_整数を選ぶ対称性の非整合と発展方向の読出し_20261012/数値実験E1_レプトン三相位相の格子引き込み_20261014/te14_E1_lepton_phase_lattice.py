# -*- coding: utf-8 -*-
"""
第十四思考実験 数値実験 E1（2026-10-14）
課題4・課題1・課題2：空間側の拘束 C_A（四次元格子が許す回転位数）と、
発展側の拘束 C_B（三相＝(1,1,1)に垂直な面の整数位相＝アイゼンシュタイン格子 Z6）
を同一状態に同時に課したときの、引き込み・非完全ロック・準安定を調べる。

状態    : アイゼンシュタイン整数 z=a+b*omega (omega=e^{2pi i/3})
          レプトン三成分ベクトルの (1,1,1) に垂直な成分が住む A2 格子。
更新則  : z_{k+1} = round_A2( e^{2pi i/m} z_k )   （外部時計・隠れ変数なし）
          回転 = C_A の要求、丸め = C_B の要求。
m       : 対照（C_A∩C_B で両立）3, 6 ／ 試験（C_A のみ、A2 と非両立）4, 5, 8, 10, 12
読み出し（事前固定）:
  1. 軌道が落ちる周期軌道の周期 P と巻き数 W、実効回転数 rho=W/P と 1/m の差
  2. ノルム（保存量）N=a^2-ab+b^2 の変化回数と、変化の間の滞在長
  3. 丸めの同点（ほぼ等距離）の件数（数値的健全性）
出力    : te14_E1_orbits.csv, te14_E1_dwell.csv, te14_E1_summary.json
"""
import math, csv, json
SQ3 = math.sqrt(3.0)
M_CONTROL = (3, 6)
M_TEST = (4, 5, 8, 10, 12)
NMAX = 1600          # 初期点のノルム上限（半径40）
MAXIT = 2_000_000
TIE_EPS = 1e-9

def to_xy(a, b):
    return a - 0.5*b, (SQ3/2)*b

def round_a2(x, y):
    """A2 格子の最近点（候補16点の全探索、同点は辞書順、同点余裕を返す）"""
    b0 = y/(SQ3/2); a0 = x + 0.5*b0
    fa, fb = math.floor(a0), math.floor(b0)
    cands = []
    for a in (fa-1, fa, fa+1, fa+2):
        for b in (fb-1, fb, fb+1, fb+2):
            X, Y = to_xy(a, b)
            cands.append(((X-x)**2+(Y-y)**2, a, b))
    cands.sort()
    return cands[0][1], cands[0][2], cands[1][0]-cands[0][0]

def norm(a, b):
    return a*a - a*b + b*b

def run_orbit(a, b, m):
    c, s = math.cos(2*math.pi/m), math.sin(2*math.pi/m)
    seen = {}
    path = []
    ties = 0
    k = 0
    while (a, b) not in seen and k < MAXIT:
        seen[(a, b)] = k; path.append((a, b))
        x, y = to_xy(a, b)
        na, nb, gap = round_a2(c*x - s*y, s*x + c*y)
        if gap < TIE_EPS: ties += 1
        a, b = na, nb; k += 1
    if (a, b) not in seen:
        return None
    k0 = seen[(a, b)]
    cyc = path[k0:]
    P = len(cyc)
    ang = 0.0
    for i in range(P):
        x0, y0 = to_xy(*cyc[i]); x1, y1 = to_xy(*cyc[(i+1) % P])
        d = math.atan2(y1, x1) - math.atan2(y0, x0)
        ang += (d + math.pi) % (2*math.pi) - math.pi
    W = round(ang/(2*math.pi))
    ns = [norm(*p) for p in cyc]
    changes = [i for i in range(P) if ns[i] != ns[(i+1) % P]]
    dwell = []
    if changes:
        for j in range(len(changes)):
            dwell.append((changes[(j+1) % len(changes)] - changes[j]) % P or P)
    canon = min(cyc)
    return dict(transient=k0, P=P, W=W, ns_min=min(ns), ns_max=max(ns),
                n_changes=len(changes), dwell=dwell, ties=ties, cycle_id=canon)

def main():
    starts = [(a, b) for a in range(0, 41) for b in range(0, a+1)
              if 0 < norm(a, b) <= NMAX]          # 60°扇形（対称性で代表）
    rows, dwell_rows, summary = [], [], {}
    for m in M_CONTROL + M_TEST:
        cycles = {}
        locked = unlocked = 0; tie_total = 0
        for (a, b) in starts:
            r = run_orbit(a, b, m)
            if r is None: continue
            tie_total += r['ties']
            rho = r['W']/r['P']; d = rho - 1/m
            lock = abs(d) < 1e-15 and r['n_changes'] == 0
            if abs(d) < 1e-15: locked += 1
            else: unlocked += 1
            rows.append([m, a, b, norm(a, b), r['transient'], r['P'], r['W'],
                         f"{rho:.15f}", f"{d:+.6e}", r['n_changes'],
                         r['ns_min'], r['ns_max'], int(lock)])
            cid = r['cycle_id']
            if cid not in cycles:
                cycles[cid] = (r['P'], r['W'], r['n_changes'])
                for L in r['dwell']:
                    dwell_rows.append([m, cid[0], cid[1], r['P'], r['W'], L])
        summary[m] = dict(
            n_starts=len(starts), rho_eq_1_over_m=locked, rho_ne_1_over_m=unlocked,
            n_distinct_cycles=len(cycles),
            cycles_rho_ne=sorted({f"{W}/{P}" for (P, W, _) in cycles.values()
                                  if abs(W/P-1/m) > 1e-15}, key=lambda t: eval(t)),
            cycles_with_norm_change=sum(1 for v in cycles.values() if v[2] > 0),
            ties=tie_total)
        print(m, summary[m], flush=True)
    with open('te14_E1_orbits.csv', 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['m','a0','b0','N0','transient','P','W','rho','rho_minus_1_over_m',
                    'norm_changes','N_min','N_max','exact_lock'])
        w.writerows(rows)
    with open('te14_E1_dwell.csv', 'w', newline='') as f:
        w = csv.writer(f); w.writerow(['m','cycle_a','cycle_b','P','W','dwell_len']); w.writerows(dwell_rows)
    with open('te14_E1_summary.json', 'w') as f:
        json.dump({str(k): v for k, v in summary.items()}, f, ensure_ascii=False, indent=1)

if __name__ == '__main__':
    main()
