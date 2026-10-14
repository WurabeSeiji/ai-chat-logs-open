# -*- coding: utf-8 -*-
"""
第十四思考実験 数値実験 E5（2026-10-14）
電子自身の質量・エネルギーが四次元ゲージ（フルヴィッツ格子）を僅かに歪め、ゲージから外れることを確かめる。
※ 同じフォルダに E2 の te14_E2_hurwitz_124_23.py と E3 の te14_E3_ground_and_levels.py を置いて実行する。

系      : E4 と同じ（四次元フルヴィッツ格子、Z4 軸、R_{124,23} の位相回転、外部時計なし）
歪み    : 更新ごとに、状態自身のクーロン的自己エネルギー alpha/r（r=|q|、格子単位）に比例して
          丸めの基準格子を局所的に伸縮させる。  x <- u q ;  x <- x (1 + s*alpha/r) ; q <- round_H(x)
          s=+1（伸張）と s=-1（収縮）の二通りを事前に固定して両方調べる。外部の自由パラメータなし。
自己無撞着: alpha_0 = 1/x_{23/124} から始め、周期軌道の位相幾何的回転数 rho=W/P（多数ステップの平均）から
          alpha_read^{-1} = 4pi/sin^4(pi rho) を読み、その平均で alpha を更新して3回反復する。
読み出し : (1) 23/124 からの外れ（rho != 23/124 の軌道の割合）、(2) 周期内で局所回転数（124窓）が
          23/124 からずれる時間割合 f、(3) alpha_read^{-1} の分布と平均。観測値は入力しない。
乱数種  : 20261014（初期点）。半径 64 で 24 点。
"""
import math, random, json, csv
from fractions import Fraction
from te14_E2_hurwitz_124_23 import qmul, axis, round_h
from te14_E3_ground_and_levels import comps, frame
AX = 'Z4'; PHI = 2*math.pi*23/124; TARGET = Fraction(23, 124)
xread = lambda r: 4*math.pi/math.sin(math.pi*r)**4
MAXIT = 2_000_000

def rand_points(n, R, seed):
    rnd = random.Random(seed); pts = []
    for _ in range(n):
        v = [rnd.gauss(0, 1) for _ in range(4)]; s = math.sqrt(sum(x*x for x in v)); q = [R*x/s for x in v]
        pts.append(tuple(2*math.floor(x)+1 for x in q) if rnd.random() < 0.5 else tuple(2*round(x) for x in q))
    return pts

def step(p, u, alpha, s):
    q = tuple(v/2 for v in p); r = math.sqrt(sum(v*v for v in q))
    x = qmul(u, q)
    k = 1.0 + s*alpha/r if r > 0 else 1.0
    return round_h(tuple(v*k for v in x))[0]

def run(p0, alpha, s):
    I = axis(AX); sinv = frame(AX)
    u = (math.cos(PHI),) + tuple(math.sin(PHI)*v for v in I[1:])
    seen = {}; path = []; p = p0; k = 0
    while p not in seen and k < MAXIT:
        seen[p] = k; path.append(p); p = step(p, u, alpha, s); k += 1
    if p not in seen: return None
    cyc = path[seen[p]:]; P = len(cyc)
    ang = []
    for i in range(P):
        a = comps(tuple(v/2 for v in cyc[i]), sinv)[0]; b = comps(tuple(v/2 for v in cyc[(i+1) % P]), sinv)[0]
        if abs(a) < 1e-9 or abs(b) < 1e-9: return None
        ang.append((math.atan2(b.imag, b.real)-math.atan2(a.imag, a.real)+math.pi) % (2*math.pi)-math.pi)
    W = round(sum(ang)/(2*math.pi)); rho = Fraction(W, P)
    # 局所回転数（124 窓、周期的に巡回）が 23/124 からずれる時間割合
    off = 0
    if P >= 124:
        ext = ang + ang[:124]; s0 = sum(ext[:124])
        for i in range(P):
            if abs(s0/(2*math.pi) - 23) > 0.5: off += 1
            s0 += ext[i+124] - ext[i]
    else:
        off = P if rho != TARGET else 0
    rs = [math.sqrt(sum((v/2)**2 for v in c)) for c in cyc]
    return dict(P=P, rho=rho, f_off=off/P, rmin=min(rs), rmax=max(rs), transient=seen[p])

def main():
    pts = rand_points(24, 64, 20261014)
    out = {}; rows = []
    for s in (+1, -1):
        alpha = 1/xread(23/124); hist = []
        for it in range(3):
            res = [run(p, alpha, s) for p in pts]; res = [r for r in res if r]
            ai = [xread(float(r['rho'])) for r in res]
            locked = sum(1 for r in res if r['rho'] == TARGET)
            fmean = sum(r['f_off'] for r in res)/len(res)
            amean = sum(ai)/len(ai)
            hist.append(dict(iter=it, alpha_inv_in=1/alpha, n=len(res), locked=locked, f_off_mean=fmean,
                             alpha_read_inv_mean=amean, alpha_read_inv_median=sorted(ai)[len(ai)//2],
                             rhos=sorted({str(r['rho']) for r in res}, key=lambda t: Fraction(t))[:12]))
            for r, a in zip(res, ai):
                rows.append([s, it, 1/alpha, r['P'], str(r['rho']), f"{a:.9f}", f"{r['f_off']:.6f}", f"{r['rmin']:.3f}", f"{r['rmax']:.3f}", r['transient']])
            print('s', s, hist[-1], flush=True)
            alpha = 1/amean
        out[str(s)] = hist
    json.dump(out, open('te14_E5_summary.json', 'w'), ensure_ascii=False, indent=1)
    with open('te14_E5_rows.csv', 'w', newline='') as f:
        w = csv.writer(f); w.writerow(['s','iter','alpha_inv_in','P','rho','alpha_read_inv','f_off','r_min','r_max','transient']); w.writerows(rows)

if __name__ == '__main__':
    main()
