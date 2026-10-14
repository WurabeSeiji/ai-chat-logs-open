# -*- coding: utf-8 -*-
"""
第十四思考実験 数値実験 E4（2026-10-14）
四次元フルヴィッツ格子＋位相回転（E2/E3 の系、外部時計なし、位相幾何的回転数 rho_top=W/P）で、
(A) 選択：回転を 23/124 と近くの有限位数根に替え、同じ半径で「自分の比にロックする」率を比べる
(B) 引き込み：回転を 23/124+eps（eps は sqrt2 倍した無理数的なずれ）にしたとき、rho_top が 23/124 に
    引き込まれる率。比較として 23/122+eps も同様に調べる
(C) 回復：23/124 にロックした周期軌道の一点に格子1単位の擾乱を加え、再ロックまでの過渡長と、
    過渡中の局所回転数（124更新窓）から読んだ alpha_read^{-1} の平均（向き）を記録する
事前固定：軸は Z4（lambda_a^31=i の四相）。乱数種 20261014。観測値は入力しない。
(C) は遷移中の時間割合 f の「1回あたりの持続」だけを与える。f そのものには独立なスリップ頻度が要る。
※ 同じフォルダに E2 の te14_E2_hurwitz_124_23.py と E3 の te14_E3_ground_and_levels.py を置いて実行する。
"""
import math, random, csv, json
from fractions import Fraction
from te14_E2_hurwitz_124_23 import qmul, axis, round_h
from te14_E3_ground_and_levels import orbit, comps, frame
random.seed(20261014)
AX = 'Z4'
xread = lambda r: 4*math.pi/math.sin(math.pi*r)**4

def rand_point(R):
    v = [random.gauss(0, 1) for _ in range(4)]; s = math.sqrt(sum(x*x for x in v))
    q = [R*x/s for x in v]
    return tuple(2*math.floor(x)+1 for x in q) if random.random() < 0.5 else tuple(2*round(x) for x in q)

def rho_of(r):
    Ws = [w for w in r['W'] if w is not None]
    if not Ws or (len(Ws) == 2 and Ws[0] != Ws[1]): return None
    return Fraction(Ws[0], r['P'])

def part_A(rows):
    targets = ['23/124', '23/122', '23/126', '22/119', '24/131', '17/92', '18/97']
    out = {}
    for t in targets:
        fr = Fraction(t)
        for R in (32, 64, 96):
            hit = n = 0
            for _ in range(30):
                r = orbit(rand_point(R), AX, 2*math.pi*float(fr))
                if r is None: continue
                rh = rho_of(r)
                if rh is None: continue
                n += 1; hit += (rh == fr)
            out[f"{t}@R{R}"] = (hit, n); rows.append(['A', t, R, hit, n])
        print('A', t, [out[f"{t}@R{R}"] for R in (32, 64, 96)], flush=True)
    return out

def part_B(rows):
    out = {}
    for base in ('23/124', '23/122'):
        fb = Fraction(base)
        for e in (1e-3, 1e-4, 1e-5, 1e-6, -1e-6, -1e-5, -1e-4, -1e-3):
            eps = e*math.sqrt(2)
            hit = n = 0
            for R in (64, 96):
                for _ in range(15):
                    r = orbit(rand_point(R), AX, 2*math.pi*(float(fb)+eps))
                    if r is None: continue
                    rh = rho_of(r)
                    if rh is None: continue
                    n += 1; hit += (rh == fb)
            out[f"{base}{eps:+.3e}"] = (hit, n); rows.append(['B', base, f"{eps:+.3e}", hit, n])
        print('B', base, {k: v for k, v in out.items() if k.startswith(base)}, flush=True)
    return out

def local_rates(path, sinv):
    ang = []
    for i in range(len(path)-1):
        a = comps(tuple(v/2 for v in path[i]), sinv)[0]; b = comps(tuple(v/2 for v in path[i+1]), sinv)[0]
        ang.append(0.0 if abs(a) < 1e-9 or abs(b) < 1e-9 else
                   (math.atan2(b.imag, b.real)-math.atan2(a.imag, a.real)+math.pi) % (2*math.pi)-math.pi)
    return ang

def part_C(rows):
    I = axis(AX); sinv = frame(AX); phi = 2*math.pi*23/124
    u = (math.cos(phi),) + tuple(math.sin(phi)*v for v in I[1:])
    res = []; tries = 0
    while len(res) < 20 and tries < 400:
        tries += 1
        p0 = rand_point(64); r = orbit(p0, AX, phi)
        if r is None or rho_of(r) != Fraction(23, 124): continue
        p = p0
        for _ in range(r['transient']): p, _ = round_h(qmul(u, tuple(v/2 for v in p)))
        j = random.randrange(4); sgn = random.choice((-2, 2))
        pk = tuple(v+(sgn if k == j else 0) for k, v in enumerate(p))     # 格子1単位の擾乱
        seen = {}; path = []; q = pk; k = 0
        while q not in seen and k < 3_000_000:
            seen[q] = k; path.append(q); q, _ = round_h(qmul(u, tuple(v/2 for v in q))); k += 1
        k0 = seen[q]
        rr = orbit(pk, AX, phi); rh = rho_of(rr)
        tr = path[:k0+1]
        ang = local_rates(tr, sinv) if len(tr) > 1 else []
        rate = (sum(ang)/(2*math.pi*len(ang))) if ang else float('nan')
        res.append(dict(transient=k0, relock=str(rh), tr_rate=rate,
                        tr_alpha_inv=(xread(rate) if ang else float('nan'))))
        rows.append(['C', str(rh), k0, f"{rate:.9f}", f"{res[-1]['tr_alpha_inv']:.6f}"])
    print('C', res, flush=True)
    return res

def main():
    rows = []
    A = part_A(rows); B = part_B(rows); C = part_C(rows)
    with open('te14_E4_rows.csv', 'w', newline='') as f:
        w = csv.writer(f); w.writerow(['part', 'key', 'param', 'v1', 'v2']); w.writerows(rows)
    json.dump(dict(A={k: list(v) for k, v in A.items()}, B={k: list(v) for k, v in B.items()}, C=C),
              open('te14_E4_summary.json', 'w'), ensure_ascii=False, indent=1)

if __name__ == '__main__':
    main()
