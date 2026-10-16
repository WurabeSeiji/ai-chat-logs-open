# -*- coding: utf-8 -*-
"""
第十四思考実験 数値実験 E6b（2026-10-16）：E6 の系（接着した対 (x, y) を毎ステップ E₈ 格子へ丸める）の上の三つの追加試験。
E6 の te14_E6_e8_glued_selection.py の関数（nearest_pair, init_pair, winding, hurwitz_candidates, glued）をそのまま使い、
変更は「回転角」と「x と y に別々の回転を掛けること（試験 B のみ）」だけ。軸は Z4、観測値（137、α）は入力しない。

試験 A（E₈ の指数の回転の自己ロック）：回転角 2π m/30（m ∈ E₈ の指数 {1, 7, 11, 13}）と、近い角度の対照
  {7/32, 11/32, 9/28, 13/28}。x の回転数 W/P が自分の比に一致する率。mode glued（E₈）と same（y=x、接着なし＝E4 と同じ系）。
試験 B（二つの整数条件の非整合）：x に Coxeter の位相 2π·23/30、y に四分の一回転 2π/4 を掛け、毎ステップ接着条件で丸める。
  周期 P と x・y の回転数の分布を記録し、30 とも 4 とも違う整数（31、124 など）が現れるかを見る。対照は接着なし（free）。
試験 C（23/120 と 23/124 の競合）：回転角を 23/124 と 23/120 の間の点 r(t) = 23/124 + t(23/120 − 23/124)、
  t ∈ {0.25, 0.5, 0.75} に置き（t=0.5 は 23/122 に一致しないよう √2×10⁻⁶ を加える）、x の回転数がどの有理数に落ちるかを数える。
半径 64、各 20 点、周期探索の上限 3,000,000。乱数種 20261016＋試験内の番号。
"""
import math, random, csv, json, sys, time
from fractions import Fraction
from collections import Counter
from te14_E2_hurwitz_124_23 import qmul, axis, norm2x
from te14_E3_ground_and_levels import frame
from te14_E4_selection_pullin_recovery import rand_point, rho_of, AX
from te14_E6_e8_glued_selection import nearest_pair, init_pair, winding

MAXIT = 3_000_000; R0 = 64; NP = 20; SEED = 20261016

def orbit2(p0x, p0y, phix, phiy, glue):
    I = axis(AX); sinv = frame(AX)
    ux = (math.cos(phix),) + tuple(math.sin(phix)*v for v in I[1:])
    uy = (math.cos(phiy),) + tuple(math.sin(phiy)*v for v in I[1:])
    seen = {}; path = []; s = (p0x, p0y); k = 0
    while s not in seen and k < MAXIT:
        seen[s] = k; path.append(s)
        s = nearest_pair(qmul(ux, tuple(v/2 for v in s[0])), qmul(uy, tuple(v/2 for v in s[1])), glue); k += 1
    if s not in seen: return None
    k0 = seen[s]; cyc = path[k0:]; P = len(cyc)
    Wx = [winding([c[0] for c in cyc], sinv, i) for i in (0, 1)]
    Wy = [winding([c[1] for c in cyc], sinv, i) for i in (0, 1)]
    return dict(P=P, W=Wx, Wy=Wy, transient=k0)

def probe(phix, phiy, mode, n=NP):
    glue = (mode == 'glued'); res = []
    for _ in range(n):
        px, py = init_pair(R0, mode)
        r = orbit2(px, py, phix, phiy, glue)
        if r is None: res.append(None); continue
        res.append(dict(P=r['P'], transient=r['transient'], rx=rho_of(r), ry=rho_of(dict(W=r['Wy'], P=r['P']))))
    return res

def summarize(res, target=None):
    ok = [r for r in res if r is not None]
    d = dict(n=len(ok), none=len(res)-len(ok),
             periods=Counter(r['P'] for r in ok).most_common(8),
             rho_x=Counter(str(r['rx']) for r in ok).most_common(8),
             rho_y=Counter(str(r['ry']) for r in ok).most_common(8),
             transient_median=sorted(r['transient'] for r in ok)[len(ok)//2] if ok else None)
    if target is not None: d['lock_x'] = sum(1 for r in ok if r['rx'] == target)
    return d

def test_A(mode):
    out = {}
    for i, t in enumerate(['1/30', '7/30', '11/30', '13/30', '7/32', '11/32', '9/28', '13/28']):
        random.seed(SEED + i); fr = Fraction(t); t0 = time.time()
        out[t] = summarize(probe(2*math.pi*float(fr), 2*math.pi*float(fr), mode), fr)
        print('A', mode, t, out[t]['lock_x'], '/', out[t]['n'], 'none', out[t]['none'], f'{time.time()-t0:.0f}s', flush=True)
    return out

def test_B(mode):
    random.seed(SEED + 100); t0 = time.time()
    res = probe(2*math.pi*23/30, 2*math.pi/4, mode, n=30)
    s = summarize(res); s['period_factors'] = Counter(f for r in res if r for f in factor(r['P'])).most_common(10)
    print('B', mode, s, f'{time.time()-t0:.0f}s', flush=True)
    return s

def factor(n):
    f = []; p = 2
    while p*p <= n:
        while n % p == 0: f.append(p); n //= p
        p += 1
    if n > 1: f.append(n)
    return sorted(set(f))

def test_C(mode, t):
    a, b = 23/124, 23/120; r = a + t*(b-a) + (math.sqrt(2)*1e-6 if t == 0.5 else 0.0)
    random.seed(SEED + 200 + int(t*100)); t0 = time.time()
    s = summarize(probe(2*math.pi*r, 2*math.pi*r, mode))
    s['rotation'] = r
    print('C', mode, t, f'{r:.9f}', s['rho_x'], 'none', s['none'], f'{time.time()-t0:.0f}s', flush=True)
    return s

def main():
    test, mode = sys.argv[1], sys.argv[2]
    if test == 'A': out = test_A(mode)
    elif test == 'B': out = test_B(mode)
    else: out = test_C(mode, float(sys.argv[3]))
    tag = f"{test}_{mode}" + (f"_t{sys.argv[3]}" if test == 'C' else '')
    json.dump(out, open(f'te14_E6b_{tag}.json', 'w'), ensure_ascii=False, indent=1, default=str)

if __name__ == '__main__':
    main()
