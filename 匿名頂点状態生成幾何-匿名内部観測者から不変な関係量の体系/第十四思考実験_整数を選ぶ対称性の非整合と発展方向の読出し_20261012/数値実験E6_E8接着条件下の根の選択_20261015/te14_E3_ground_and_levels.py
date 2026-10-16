# -*- coding: utf-8 -*-
"""
第十四思考実験 数値実験 E3（2026-10-14）
E2 の系（四次元フルヴィッツ格子＋R_{124,23} の位相回転、外部時計なし）で、
基底状態（完全ロック）と準安定な準位の系列を分ける。
※ 同じフォルダに E2 の te14_E2_hurwitz_124_23.py を置いて実行する。

E2 からの修正：回転数を、等傾回転の二つの複素成分 z1, z2 の巻き数（整数）から求める
位相幾何的な回転数 rho_top = W/P に改める（E2 の幾何平均角 rho_geo は四次元で加法的でない）。

定義（事前固定）:
  基底状態 : rho_top = 23/124 かつ ノルム一定（周期内で保存量が変化しない）＝完全ロック
  準位     : 周期軌道の rho_top（有理数）。準位の「エネルギー」は周波数差 dE ∝ rho_top - 23/124
  読み出し : 基底状態の alpha_read^{-1} = 4pi/sin^4(pi rho_top)（K4 の対応式）、
             m_e = 2 h R_inf alpha_read^2 / c
  過渡     : 周期軌道に入るまでの過渡区間で、124 更新窓の巻き数（局所回転数）の推移
注意     : E2 の結果（近似分数 3/16, 5/27, 13/70, 18/97 の出現）を見た後に設計した。
"""
import math, random, csv, json
from fractions import Fraction
from collections import Counter, defaultdict
from te14_E2_hurwitz_124_23 import qmul, axis, round_h, norm2x, PHI_TEST
random.seed(20261014)
MAXIT = 3_000_000
TARGET = Fraction(23, 124)

def frame(axn):
    """s^{-1} を返す：s i s^{-1} = I となる単位四元数 s"""
    I = axis(axn)
    if axn == 'Z4': return (1.0, 0.0, 0.0, 0.0)
    iv = (1.0, 0.0, 0.0); Iv = I[1:]
    n = (iv[1]*Iv[2]-iv[2]*Iv[1], iv[2]*Iv[0]-iv[0]*Iv[2], iv[0]*Iv[1]-iv[1]*Iv[0])
    ln = math.sqrt(sum(x*x for x in n)); n = tuple(x/ln for x in n)
    beta = math.acos(sum(a*b for a, b in zip(iv, Iv)))
    s = (math.cos(beta/2),) + tuple(math.sin(beta/2)*x for x in n)
    return (s[0], -s[1], -s[2], -s[3])

def comps(q, sinv):
    p = qmul(sinv, q)            # s^{-1} q : 左乗 e^{I phi} は s e^{i phi} s^{-1} q
    return complex(p[0], p[1]), complex(p[2], p[3])

def orbit(p0, axn, phi=PHI_TEST):
    I = axis(axn); sinv = frame(axn)
    u = (math.cos(phi),) + tuple(math.sin(phi)*v for v in I[1:])
    seen = {}; path = []; p = p0; k = 0
    while p not in seen and k < MAXIT:
        seen[p] = k; path.append(p)
        p, _ = round_h(qmul(u, tuple(v/2 for v in p))); k += 1
    if p not in seen: return None
    k0 = seen[p]; cyc = path[k0:]; P = len(cyc)
    W = []
    for idx in (0, 1):
        tot = 0.0; ok = True
        for i in range(P):
            a = comps(tuple(v/2 for v in cyc[i]), sinv)[idx]
            b = comps(tuple(v/2 for v in cyc[(i+1) % P]), sinv)[idx]
            if abs(a) < 1e-9 or abs(b) < 1e-9: ok = False; break
            tot += (math.atan2(b.imag, b.real) - math.atan2(a.imag, a.real) + math.pi) % (2*math.pi) - math.pi
        W.append(round(tot/(2*math.pi)) if ok else None)
    ns = [norm2x(c) for c in cyc]
    nch = sum(1 for i in range(P) if ns[i] != ns[(i+1) % P])
    loc = []
    tr = path[:k0]
    if len(tr) > 248:
        ang = []
        for i in range(len(tr)-1):
            a = comps(tuple(v/2 for v in tr[i]), sinv)[0]; b = comps(tuple(v/2 for v in tr[i+1]), sinv)[0]
            ang.append(0.0 if abs(a) < 1e-9 or abs(b) < 1e-9 else
                       (math.atan2(b.imag, b.real) - math.atan2(a.imag, a.real) + math.pi) % (2*math.pi) - math.pi)
        for s0 in range(0, len(ang)-124, 124):
            loc.append(sum(ang[s0:s0+124])/(2*math.pi*124))
    return dict(P=P, W=W, nch=nch, nmin=min(ns)/4, nmax=max(ns)/4, transient=k0, loc=loc)

def main():
    RADII = [16, 24, 32, 48, 64, 96]; NPER = 40
    h = 6.62607015e-34; c = 299792458.0; Rinf = 10973731.568157; eV = 1.602176634e-19
    rows = []; levels = defaultdict(Counter); ground = []; casc = []
    for axn in ('Z4', 'Z6'):
        for R in RADII:
            for _ in range(NPER):
                v = [random.gauss(0, 1) for _ in range(4)]; sv = math.sqrt(sum(x*x for x in v))
                q = [R*x/sv for x in v]
                p0 = tuple(2*math.floor(x)+1 for x in q) if random.random() < 0.5 else tuple(2*round(x) for x in q)
                r = orbit(p0, axn)
                if r is None: continue
                Ws = [w for w in r['W'] if w is not None]
                if not Ws or (len(Ws) == 2 and Ws[0] != Ws[1]):
                    rows.append([axn, R, *p0, r['P'], 'W_inconsistent', r['W'], r['nch'], r['nmin'], r['nmax'], r['transient']]); continue
                rho = Fraction(Ws[0], r['P'])
                levels[axn][str(rho)] += 1
                is_ground = (rho == TARGET and r['nch'] == 0)
                if is_ground: ground.append((axn, R, p0))
                rows.append([axn, R, *p0, r['P'], str(rho), f"{float(rho-TARGET):+.6e}", r['nch'], r['nmin'], r['nmax'], r['transient']])
                if r['loc']: casc.append([axn, R, str(rho)] + [f"{x:.6f}" for x in r['loc'][:40]])
            print(axn, R, dict(levels[axn].most_common(8)), 'ground so far', len(ground), flush=True)
    ai_g = 4*math.pi/math.sin(math.pi*23/124)**4
    me_g = 2*h*Rinf*ai_g**2/c*c*c/eV/1e6
    out = dict(levels={k: dict(v) for k, v in levels.items()}, n_ground=len(ground),
               ground_examples=[list(g[:2])+list(g[2]) for g in ground[:10]],
               ground_alpha_read_inv=ai_g, ground_me_MeV=me_g)
    print(json.dumps(out, ensure_ascii=False))
    with open('te14_E3_orbits.csv', 'w', newline='') as f:
        w = csv.writer(f); w.writerow(['axis','R','a2','b2','c2','d2','P','rho_top','rho_minus_23_124','norm_changes','N_min','N_max','transient']); w.writerows(rows)
    with open('te14_E3_cascade.csv', 'w', newline='') as f:
        csv.writer(f).writerows(casc)
    with open('te14_E3_summary.json', 'w') as f: json.dump(out, f, ensure_ascii=False, indent=1)

if __name__ == '__main__':
    main()
