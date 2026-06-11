#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 補遺33: 干渉閉合テスト
# 複合体枠: 周期 2R のトーラス、基本波 λ=2R, ν=1/(2R)  (2026-06-11 指定)
# E1: 単一セル(R=1)の奇数倍音構成とエッジ幅
# E2: 奇数入れ子(ロッキング)定理の数値確認
# E3: R=1,2,3,5 の 4D 占有場 — 面純度(per-axis odd 判定)・毛(低域)テスト・アイ開口
import numpy as np, itertools, math, sys

OUT = []
def log(s):
    print(s); OUT.append(s)

# ---------- E1: 1D 単一セル ----------
log("="*70)
log("E1: 単一セル (R=1, 周期2, 箱幅1) — 奇数純度とエッジ幅")
N = 16384
P = 2.0
x = (np.arange(N) + 0.5) / N * P          # half-pixel grid, edges at 0.5, 1.5
xc = np.minimum(x % P, P - (x % P))        # distance to 0 on torus
box = (xc < 0.5).astype(np.float64)
c = np.fft.rfft(box) / N
mag2 = np.abs(c)**2
ac = mag2.copy(); ac[0] = 0.0
odd_E  = ac[1::2].sum()
even_E = ac[2::2].sum()
log(f"  占有率 = {box.mean():.6f}  (1/2 が厳密半波条件)")
log(f"  偶数倍音エネルギー / AC 全体 = {even_E/(odd_E+even_E):.3e}  (0 なら厳密奇数純)")

def truncate_1d(c, N, nmax, odd_only=True):
    cc = np.zeros_like(c)
    cc[0] = c[0]
    for n in range(1, min(nmax, len(c)-1) + 1):
        if (not odd_only) or (n % 2 == 1):
            cc[n] = c[n]
    return np.fft.irfft(cc * N, n=N)

def edge_width_1090(f, x, edge=0.5, win=0.5):
    m = (x > edge - win) & (x < edge + win)
    xs, fs = x[m], f[m]
    # falling edge at x=edge: last crossing of 0.9 before, first crossing of 0.1 after
    try:
        i90 = np.where((fs[:-1] >= 0.9) & (fs[1:] < 0.9))[0]
        i90 = i90[np.argmin(np.abs(xs[i90] - edge))]
        x90 = xs[i90] + (0.9 - fs[i90]) * (xs[i90+1]-xs[i90]) / (fs[i90+1]-fs[i90])
        i10 = np.where((fs[:-1] >= 0.1) & (fs[1:] < 0.1))[0]
        i10 = i10[np.argmin(np.abs(xs[i10] - edge))]
        x10 = xs[i10] + (0.1 - fs[i10]) * (xs[i10+1]-xs[i10]) / (fs[i10+1]-fs[i10])
        return x10 - x90
    except Exception:
        return float('nan')

log("  打ち切り N (奇数のみ) | エッジ幅(10-90%) | 1/(2*nu_max)=1/N | オーバーシュート")
for Nh in [1, 3, 5, 9, 19, 39]:
    f = truncate_1d(c, N, Nh)
    w = edge_width_1090(f, x)
    ov = f.max() - 1.0
    log(f"    N={Nh:3d}  w={w:.4f}   pred={1.0/Nh:.4f}   overshoot={ov:+.4f}")

# ---------- E2: 奇数入れ子定理 ----------
log("="*70)
log("E2: 奇数入れ子定理 — セル梯子(周期2)が複合体梯子(周期2R)に乗る条件")
for R in [2, 3, 4, 5]:
    Pc = 2.0 * R
    Nc = 4096 * R
    xx = (np.arange(Nc) + 0.5) / Nc * Pc
    xcc = np.minimum(xx % 2.0, 2.0 - (xx % 2.0))   # 周期2でタイルしたセル波
    cellwave = (xcc < 0.5).astype(np.float64)
    cf = np.fft.rfft(cellwave) / Nc
    m2 = np.abs(cf)**2; m2[0] = 0.0
    tot = m2.sum()
    supp = np.where(m2 > tot * 1e-12)[0]
    on_mult = np.all(supp % R == 0)
    odd_frac = m2[1::2].sum() / tot
    log(f"  R={R}: 台は n=R*m のみ? {on_mult} / 奇数 n に乗る割合 = {odd_frac:.6f}"
        f"  ({'全部乗る' if odd_frac>0.999 else '乗れない' if odd_frac<1e-6 else '部分'})")

# 1D 直感: 箱幅/周期 比と奇数純度
log("  [1D 参考] 箱幅 w / 周期 P と奇数純度 (複合体の軸断面 w=(2k_max+1)/2R):")
for (w_, P_) in [(1,2), (3,4), (5,6), (9,10)]:
    Nn = 16384
    xx = (np.arange(Nn) + 0.5) / Nn * P_
    xcc = np.minimum(xx % P_, P_ - (xx % P_))
    bb = (xcc < w_/2.0).astype(np.float64)
    cf = np.fft.rfft(bb) / Nn
    m2 = np.abs(cf)**2; m2[0] = 0.0
    pur = m2[1::2].sum() / m2.sum()
    log(f"    w/P = {w_}/{P_}:  奇数純度 = {pur:.4f}")

# ---------- E3: 4D 占有場 ----------
log("="*70)
log("E3: 4D 占有場 — 面純度 / 毛テスト / アイ開口")

def centers4(R):
    K = int(math.floor(R)) + 1
    out = []
    for k in itertools.product(range(-K, K+1), repeat=4):
        if sum((abs(t)+0.5)**2 for t in k) <= R*R + 1e-9:
            out.append(k)
    return out

def occupancy_field(R, ppu):
    Pi = 2 * R
    G = Pi * ppu
    o = np.zeros((G, G, G, G), dtype=np.float64)
    cs = centers4(R)
    for k in cs:
        ax = []
        for ki in k:
            s = ((ki + R) * ppu - ppu // 2) % G
            if s + ppu <= G:
                ax.append([slice(s, s+ppu)])
            else:
                ax.append([slice(s, G), slice(0, s+ppu-G)])
        for s1 in ax[0]:
            for s2 in ax[1]:
                for s3 in ax[2]:
                    for s4 in ax[3]:
                        o[s1, s2, s3, s4] = 1.0
    return o, cs, G, Pi

def analyze(R, ppu, ncut_list, lowpass_n=3):
    o, cs, G, Pi = occupancy_field(R, ppu)
    n_cells = len(cs)
    occ_frac = o.mean()
    F = np.fft.fftn(o)
    idx = (np.fft.fftfreq(G) * G).astype(int)
    # 面純度: 全成分 n_i ∈ {0, 奇数} の AC エネルギー割合
    good1 = (idx % 2 == 1) | (idx == 0)
    g1 = good1.reshape(-1,1,1,1) & good1.reshape(1,-1,1,1) & good1.reshape(1,1,-1,1) & good1.reshape(1,1,1,-1)
    E = np.abs(F)**2
    Edc = E[0,0,0,0]
    Eac = E.sum() - Edc
    Egood = E[g1].sum() - Edc      # DC は good に含まれるので引く
    purity = Egood / Eac
    del E, g1
    # 毛テスト: |n_i|<=lowpass_n で低域化し、等体積超立方体と比較
    lp = (np.abs(idx) <= lowpass_n)
    LP = lp.reshape(-1,1,1,1) & lp.reshape(1,-1,1,1) & lp.reshape(1,1,-1,1) & lp.reshape(1,1,1,-1)
    Flp = np.where(LP, F, 0)
    olp = np.fft.ifftn(Flp).real
    B = n_cells ** 0.25            # 等体積の超立方体一辺
    xax = -R + (np.arange(G) + 0.5) / ppu
    inb = (np.abs(xax) < B/2.0).astype(np.float64)
    boxf = inb.reshape(-1,1,1,1) * inb.reshape(1,-1,1,1) * inb.reshape(1,1,-1,1) * inb.reshape(1,1,1,-1)
    Fb = np.fft.fftn(boxf)
    Fblp = np.where(LP, Fb, 0)
    blp = np.fft.ifftn(Fblp).real
    hair_dist = np.sqrt(((olp - blp)**2).sum() / (blp**2).sum())
    del Flp, Fb, Fblp, blp, boxf, LP
    # アイ開口: 帯域 |n_i|<=ncut で再構成し、占有/空セル中心値の分離
    all_slots = list(itertools.product(range(-R, R), repeat=4))
    occset = set(cs)
    def slot_index(k):
        return tuple(((ki + R) * ppu) % G for ki in k)
    occ_idx = [slot_index(k) for k in cs]
    emp_idx = [slot_index(k) for k in all_slots if k not in occset]
    eyes = []
    for ncut in ncut_list:
        m1 = (np.abs(idx) <= ncut)
        M = m1.reshape(-1,1,1,1) & m1.reshape(1,-1,1,1) & m1.reshape(1,1,-1,1) & m1.reshape(1,1,1,-1)
        ff = np.fft.ifftn(np.where(M, F, 0)).real
        vo = np.array([ff[i] for i in occ_idx])
        ve = np.array([ff[i] for i in emp_idx]) if emp_idx else np.array([0.0])
        eye = vo.min() - ve.max()
        eyes.append((ncut, vo.min(), ve.max(), eye))
        del M, ff
    return dict(R=R, ppu=ppu, G=G, n_cells=n_cells, occ_frac=occ_frac,
                purity=purity, hair_dist=hair_dist, eyes=eyes, B=B)

cases = [
    (1, 8, [1, 2, 3]),
    (2, 8, [1, 2, 3, 4, 6, 8]),
    (3, 8, [1, 2, 3, 4, 6, 9, 12, 18]),
    (5, 6, [1, 2, 3, 5, 10, 15, 20, 25, 29]),
]
expected = {1: 1, 2: 9, 3: 137, 5: 1545}
for (R, ppu, ncl) in cases:
    r = analyze(R, ppu, ncl)
    log(f"  --- R={R} (2R={2*R}, grid {r['G']}^4, ppu={ppu}) ---")
    log(f"  セル数 = {r['n_cells']} (期待値 {expected[R]}) / 占有率 = {r['occ_frac']:.4f}"
        f" (連続極限 pi^2/32 = {math.pi**2/32:.4f})")
    log(f"  面純度 (全成分 n_i∈{{0,奇}} の AC割合) = {r['purity']:.4f}")
    log(f"  毛テスト (|n|<=3 低域での等体積箱との相対L2距離, 箱一辺 B={r['B']:.3f}) = {r['hair_dist']:.4f}")
    log(f"  アイ開口 (帯域 |n_i|<=ncut, m=ncut/(2R) セル単位):")
    for (ncut, vmin, vmax, eye) in r['eyes']:
        mark = "OPEN" if eye > 0 else "closed"
        log(f"    ncut={ncut:3d} (m={ncut/(2*R):.2f})  min_occ={vmin:+.3f}  max_emp={vmax:+.3f}  eye={eye:+.3f}  {mark}")

log("="*70)
log("done")
with open('/tmp/closure_results.txt', 'w') as f:
    f.write("\n".join(OUT))
