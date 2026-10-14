# -*- coding: utf-8 -*-
"""
第十四思考実験 数値実験 E2（2026-10-14）
四次元フルヴィッツ格子（整数点＋半整数点：対角線(1,1,1,1)の長さ2が整数になる四次元固有の格子）上で、
発展側の拘束 C_B：K4論文の交換作用素 U_R の反対称固有位相 lambda_a=e^{-2pi i*23/124}（R_{124,23}）による回転と、
空間側の拘束 C_A：フルヴィッツ格子への丸め（Z4 と Z6 が共存する四次元格子、位数は12以下）
を同一状態に同時に課す。外部時計なし。

状態    : 四元数 q=(a,b,c,d)、q∈Hurwitz = Z^4 ∪ (Z^4+(1/2,1/2,1/2,1/2))
更新則  : q_{k+1} = round_H( u q_k ),  u = cos(phi) + I sin(phi)   （左乗法＝等傾回転）
          phi = 2pi*23/124（試験）。対照は phi = 2pi/4, 2pi/6（格子と両立）
軸 I    : 'Z4' : I = i                （lambda_a^31 = i の90°四相の側）
          'Z6' : I = (i+j+k)/sqrt3    （三相の側：e^{I pi/3} = (1+i+j+k)/2 はフルヴィッツ単位）
初期点  : 0 < N(q) <= 9 の全フルヴィッツ点（半径3の四次元球、137セルと同じ球）
読み出し（事前固定、観測値は入力しない）:
  1. 周期軌道の周期 P、巻き数 W、実効回転数 rho = (総回転角)/(2pi P) と 23/124 の差
  2. K4論文の対応式 N(R)=4pi/(1-R)^2, R=cos^2(pi rho) による alpha_read^{-1}=4pi/sin^4(pi rho)
  3. ノルム（保存量）の変化回数と滞在長、丸めの同点数
"""
import math, csv, json, itertools
TIE = 1e-9
PHI_TEST = 2*math.pi*23/124
CASES = [('Z4', 'control_4', 2*math.pi/4), ('Z6', 'control_6', 2*math.pi/6),
         ('Z4', 'test_124_23', PHI_TEST), ('Z6', 'test_124_23', PHI_TEST)]
MAXIT = 3_000_000

def qmul(p, q):
    a1,b1,c1,d1 = p; a2,b2,c2,d2 = q
    return (a1*a2-b1*b2-c1*c2-d1*d2, a1*b2+b1*a2+c1*d2-d1*c2,
            a1*c2-b1*d2+c1*a2+d1*b2, a1*d2+b1*c2-c1*b2+d1*a2)

def axis(name):
    return (0.0, 1.0, 0.0, 0.0) if name == 'Z4' else (0.0,)+(1/math.sqrt(3),)*3

def round_h(x):
    """フルヴィッツ格子の最近点。座標は2倍した整数で返す（半整数を厳密に扱う）"""
    zi = [round(v) for v in x]
    zh = [math.floor(v)+0.5 for v in x]
    di = sum((v-w)**2 for v, w in zip(x, zi)); dh = sum((v-w)**2 for v, w in zip(x, zh))
    p = tuple(int(2*w) for w in (zi if di <= dh else zh))
    return p, abs(di-dh)

def norm2x(p):            # p は2倍座標、ノルムは sum(p^2)/4
    return sum(v*v for v in p)

def run(p0, ax, phi):
    u = (math.cos(phi),) + tuple(math.sin(phi)*v for v in ax[1:])
    seen = {}; path = []; ties = 0; p = p0; k = 0
    while p not in seen and k < MAXIT:
        seen[p] = k; path.append(p)
        x = qmul(u, tuple(v/2 for v in p))
        p, gap = round_h(x)
        if gap < TIE: ties += 1
        k += 1
    if p not in seen: return None
    cyc = path[seen[p]:]; P = len(cyc)
    ang = 0.0
    for i in range(P):
        q0 = tuple(v/2 for v in cyc[i]); q1 = tuple(v/2 for v in cyc[(i+1) % P])
        n0 = sum(v*v for v in q0); n1 = sum(v*v for v in q1)
        conj = (q0[0], -q0[1], -q0[2], -q0[3])
        r = qmul(q1, conj); s = math.sqrt(n0*n1)
        r = tuple(v/s for v in r)                       # q1 q0^{-1} を正規化
        imag = sum(r[j]*ax[j] for j in (1, 2, 3))       # 軸 I 方向成分で符号付き
        ang += math.atan2(imag, r[0])
    ns = [norm2x(c) for c in cyc]
    ch = [i for i in range(P) if ns[i] != ns[(i+1) % P]]
    dwell = [((ch[(j+1) % len(ch)]-ch[j]) % P or P) for j in range(len(ch))] if ch else []
    return dict(transient=seen[p], P=P, ang=ang, n_changes=len(ch), dwell=dwell,
                nmin=min(ns)/4, nmax=max(ns)/4, ties=ties, cid=min(cyc))

def alpha_read_inv(rho):
    return 4*math.pi/math.sin(math.pi*rho)**4

def main():
    starts = []
    for c in itertools.product(range(-3, 4), repeat=4):
        if 0 < sum(v*v for v in c) <= 9: starts.append(tuple(2*v for v in c))
    for c in itertools.product(range(-4, 4), repeat=4):
        h = tuple(2*v+1 for v in c)
        if sum(v*v for v in h) <= 36: starts.append(h)
    rows, dwell_rows, summ = [], [], {}
    for axn, label, phi in CASES:
        ax = axis(axn); key = f"{axn}_{label}"; cyc = {}; tt = 0
        target = phi/(2*math.pi)
        for p0 in starts:
            r = run(p0, ax, phi)
            if r is None: continue
            tt += r['ties']
            rho = r['ang']/(2*math.pi*r['P'])
            W = round(r['ang']/(2*math.pi))
            ai = alpha_read_inv(rho) if label.startswith('test') else float('nan')
            rows.append([key, *p0, norm2x(p0)/4, r['transient'], r['P'], W, f"{rho:.15f}",
                         f"{rho-target:+.6e}", f"{ai:.9f}", r['n_changes'], r['nmin'], r['nmax']])
            if r['cid'] not in cyc:
                cyc[r['cid']] = (r['P'], rho, r['n_changes'])
                for L in r['dwell']: dwell_rows.append([key, r['P'], L])
        rhos = [v[1] for v in cyc.values()]
        summ[key] = dict(n_starts=len(starts), n_cycles=len(cyc), ties=tt,
                         rho_exact=sum(1 for x in rhos if abs(x-target) < 1e-12),
                         rho_min=min(rhos), rho_max=max(rhos),
                         cycles_norm_change=sum(1 for v in cyc.values() if v[2] > 0),
                         periods=sorted({v[0] for v in cyc.values()})[:20])
        if label.startswith('test'):
            ai = [alpha_read_inv(x) for x in rhos]
            summ[key].update(alpha_read_inv_min=min(ai), alpha_read_inv_max=max(ai),
                             in_window_137_to_137036043=sum(1 for x in ai if 137 <= x <= 137.036042914599))
        print(key, summ[key], flush=True)
    with open('te14_E2_orbits.csv', 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['case','a2','b2','c2','d2','N0','transient','P','W','rho','rho_minus_target',
                    'alpha_read_inv','norm_changes','N_min','N_max'])
        w.writerows(rows)
    with open('te14_E2_dwell.csv', 'w', newline='') as f:
        w = csv.writer(f); w.writerow(['case','P','dwell_len']); w.writerows(dwell_rows)
    with open('te14_E2_summary.json', 'w') as f:
        json.dump(summ, f, ensure_ascii=False, indent=1)

if __name__ == '__main__':
    main()
