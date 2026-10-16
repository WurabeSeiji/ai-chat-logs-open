# -*- coding: utf-8 -*-
"""
第十四思考実験 数値実験 E6（2026-10-15）：E₈ の接着条件の下での根の選択
E4 の part A（四次元フルヴィッツ格子＋位相回転、軸 Z4、外部時計なし、回転数 rho_top=W/P）の忠実な複製に、
最小の変更として「状態を一つのフルヴィッツ整数 x から、接着条件 x − y ∈ (1+i)H を満たす対 (x, y) に替える」
ことだけを加える。接着条件は毎ステップ、回転後の丸めで課し続ける（初期値だけに課すのは不自然、著者判断 10-15）。

  状態    : (x, y) ∈ H×H、x − y ∈ (1+i)H  （この対の集合が E₈ 格子。第十四 §24.9、verify_e8.py で 240 本を確認済み）
  更新則  : (x, y) → nearest_E8( u x, u y )   u = cos(phi) + i sin(phi)、両方に同じ左乗法
  丸め    : 接着格子への最近点。x, y それぞれの近傍のフルヴィッツ点の対のうち、N(x−y) が偶数で距離の和が最小のもの
            （フルヴィッツ整数 z について z ∈ (1+i)H ⟺ N(z) が偶数）
  接着なし: 偶数条件を外すと x と y は E2 の round_h で独立に丸められる。mode 'same'（y=x）では乱数列も E4 と同じで
            x の軌道と part A の数が E4 と一致する（内部対照）。mode 'free' は y を独立に置いた対照
  読み出し: E3 の orbit と同じ巻き数 W/P を x について計算し、target と一致すればロック。y の W/P も記録
  対象    : 格子が通す mod4 の類 {17/92, 19/104, 21/116, 23/124, 25/136, 27/148, 23/120}
            窓の隣（gray 深さの図の 7 山）{23/122, 23/123, 23/125, 24/127, 22/119}
            E4 の対照 {23/126, 24/131, 18/97}
  半径    : 32, 64, 96、各 30 点。乱数種 20261015。観測値（137、α）は入力しない。
※ 同じフォルダに E2 の te14_E2_hurwitz_124_23.py、E3 の te14_E3_ground_and_levels.py、E4 の
  te14_E4_selection_pullin_recovery.py の複製を置いて実行する（E4 からは rand_point, rho_of, AX を使う）。
"""
import math, random, csv, json, sys, time
from fractions import Fraction
from te14_E2_hurwitz_124_23 import qmul, axis, norm2x, round_h
from te14_E3_ground_and_levels import comps, frame
from te14_E4_selection_pullin_recovery import rand_point, rho_of, AX

MAXIT = 3_000_000
TARGETS_MOD4 = ['17/92', '19/104', '21/116', '23/124', '25/136', '27/148', '23/120']
TARGETS_WINDOW = ['23/122', '23/123', '23/125', '24/127', '22/119']
TARGETS_E4 = ['23/126', '24/131', '18/97']
TARGETS_E4_ORIGINAL = ['23/124', '23/122', '23/126', '22/119', '24/131', '17/92', '18/97']   # E4 part A の順序
RADII = (32, 64, 96)
NPTS = 30
SEED = 20261015

def hurwitz_candidates(v):
    """v ∈ R^4 の近傍のフルヴィッツ点（2倍座標）を距離順に返す。整数点 3^4 と半整数点 3^4 の箱から取る。"""
    out = []
    base_i = [math.floor(t) for t in v]
    for d0 in (0, 1, -1):
        for d1 in (0, 1, -1):
            for d2 in (0, 1, -1):
                for d3 in (0, 1, -1):
                    zi = (base_i[0]+d0, base_i[1]+d1, base_i[2]+d2, base_i[3]+d3)
                    out.append((sum((t-w)**2 for t, w in zip(v, zi)), tuple(2*w for w in zi)))
                    zh = (base_i[0]+d0+0.5, base_i[1]+d1+0.5, base_i[2]+d2+0.5, base_i[3]+d3+0.5)
                    out.append((sum((t-w)**2 for t, w in zip(v, zh)), tuple(int(2*w) for w in zh)))
    out.sort(key=lambda t: t[0])
    return out

def glued(px, py):
    """2倍座標の対が接着条件を満たすか：N(x−y) が偶数 ⟺ 4N(x−y) = Σ(dx)^2 ≡ 0 (mod 8)"""
    return sum((a-b)**2 for a, b in zip(px, py)) % 8 == 0

def nearest_pair(vx, vy, glue):
    """(vx, vy) に最も近い (x, y)。glue=False なら E2 の round_h をそれぞれに適用（E4 と同一）。
    glue=True なら接着格子（E₈）への最近点：まず round_h の対を試し、接着条件を満たさないときだけ近傍を探す。"""
    px, _ = round_h(vx); py, _ = round_h(vy)
    if not glue or glued(px, py):
        return px, py
    cx = hurwitz_candidates(vx); cy = hurwitz_candidates(vy)
    best = None
    for dx, qx in cx:
        if best is not None and dx >= best[0]: break
        for dy, qy in cy:
            if best is not None and dx + dy >= best[0]: break
            if glued(qx, qy):
                best = (dx+dy, qx, qy); break
    return best[1], best[2]

def winding(cyc_component, sinv, idx):
    P = len(cyc_component); tot = 0.0
    for i in range(P):
        a = comps(tuple(v/2 for v in cyc_component[i]), sinv)[idx]
        b = comps(tuple(v/2 for v in cyc_component[(i+1) % P]), sinv)[idx]
        if abs(a) < 1e-9 or abs(b) < 1e-9: return None
        tot += (math.atan2(b.imag, b.real) - math.atan2(a.imag, a.real) + math.pi) % (2*math.pi) - math.pi
    return round(tot/(2*math.pi))

def orbit_pair(p0x, p0y, axn, phi, glue):
    """E3 の orbit の対版。周期 P、x と y それぞれの巻き数 W（二つの複素成分）、過渡長を返す。"""
    I = axis(axn); sinv = frame(axn)
    u = (math.cos(phi),) + tuple(math.sin(phi)*v for v in I[1:])
    seen = {}; path = []; s = (p0x, p0y); k = 0
    while s not in seen and k < MAXIT:
        seen[s] = k; path.append(s)
        vx = qmul(u, tuple(v/2 for v in s[0])); vy = qmul(u, tuple(v/2 for v in s[1]))
        s = nearest_pair(vx, vy, glue); k += 1
    if s not in seen: return None
    k0 = seen[s]; cyc = path[k0:]; P = len(cyc)
    cx = [c[0] for c in cyc]; cy = [c[1] for c in cyc]
    Wx = [winding(cx, sinv, idx) for idx in (0, 1)]
    Wy = [winding(cy, sinv, idx) for idx in (0, 1)]
    nsx = [norm2x(c) for c in cx]; nsy = [norm2x(c) for c in cy]
    return dict(P=P, W=Wx, Wy=Wy, transient=k0, nminx=min(nsx)/4, nmaxx=max(nsx)/4, nminy=min(nsy)/4, nmaxy=max(nsy)/4,
                glued_all=all(glued(a, b) for a, b in cyc))

def init_pair(R, mode):
    """mode 'same' : y = x（E4 と同じ乱数列になり、x の軌道は E4 の part A と一致する。内部対照）
       mode 'free' : y は別の rand_point(R)、接着なし
       mode 'glued': y は別の rand_point(R) を x との接着条件を満たす最近点に置く"""
    px = rand_point(R)
    if mode == 'same': return px, px
    vy = tuple(v/2 for v in rand_point(R))
    if mode == 'free': return px, tuple(int(2*v) for v in vy)
    for d, py in hurwitz_candidates(vy):
        if glued(px, py): return px, py
    raise RuntimeError('no glued partner')

def run(mode, rows, label, targets):
    glue = (mode == 'glued'); out = {}
    for t in targets:
        fr = Fraction(t)
        for R in RADII:
            hit = hity = n = none = 0; t0 = time.time()
            for _ in range(NPTS):
                p0x, p0y = init_pair(R, mode)
                r = orbit_pair(p0x, p0y, AX, 2*math.pi*float(fr), glue)
                if r is None: none += 1; continue
                rh = rho_of(r)
                if rh is None: continue
                n += 1; hit += (rh == fr)
                rhy = rho_of(dict(W=r['Wy'], P=r['P']))
                hity += (rhy == fr) if rhy is not None else 0
            out[f"{t}@R{R}"] = (hit, hity, n, none)
            rows.append([label, t, R, hit, hity, n, none, f"{time.time()-t0:.1f}"])
            print(label, t, R, (hit, hity, n), 'none', none, f"{time.time()-t0:.1f}s", flush=True)
    return out

def main():
    """引数：mode [maxit] [group]。maxit は周期探索の上限（既定 3,000,000）、group は 'all'（既定）か 'rest'
    （窓の隣＋E4 の対照）。2026-10-16 追記：接着ありでは窓の隣の根の軌道が 3,000,000 更新でも閉じないものが多く
    （23/122・半径 64 の診断で 4 本中 3 本が 150,000 更新で未閉）、1 行に数時間かかるため、残りの対象は上限を下げて
    別走行にした。閉じなかった本数は none 列に記録する。"""
    global MAXIT
    mode = sys.argv[1] if len(sys.argv) > 1 else 'glued'
    if len(sys.argv) > 2: MAXIT = int(sys.argv[2])
    group = sys.argv[3] if len(sys.argv) > 3 else 'all'
    if group.startswith('one:'):          # 2026-10-16 追記：残りの対象を 1 根ずつ並列に走らせる（乱数種は 20261015＋対象の番号）
        targets = [group[4:]]
    else:
        targets = {'all': TARGETS_MOD4 + TARGETS_WINDOW + TARGETS_E4, 'rest': TARGETS_WINDOW + TARGETS_E4,
                   'mod4': TARGETS_MOD4}[group]
    seed = SEED
    if group.startswith('one:'):
        seed = SEED + (TARGETS_MOD4 + TARGETS_WINDOW + TARGETS_E4).index(targets[0])
    suffix = mode if (group == 'all' and MAXIT == 3_000_000) else f"{mode}_{group.replace(':', '_').replace('/', '_')}_maxit{MAXIT}"
    rows = []; res = {}
    if mode == 'same':     # 内部対照：E4 part A の 7 対象、E4 と同じ順序・同じ乱数種 20261014
        random.seed(20261014); res['same'] = run('same', rows, 'same', TARGETS_E4_ORIGINAL)
    else:
        random.seed(seed); res[mode] = run(mode, rows, mode, targets)
    with open(f'te14_E6_rows_{suffix}.csv', 'w', newline='') as f:
        w = csv.writer(f); w.writerow(['mode', 'target', 'R', 'lock_x', 'lock_y', 'n', 'none', 'sec']); w.writerows(rows)
    json.dump(dict(maxit=MAXIT, group=group, seed=seed if mode != 'same' else 20261014, **{m: {k: list(v) for k, v in d.items()} for m, d in res.items()}),
              open(f'te14_E6_summary_{suffix}.json', 'w'), ensure_ascii=False, indent=1)

if __name__ == '__main__':
    main()
