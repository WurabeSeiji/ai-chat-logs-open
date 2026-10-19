#!/usr/bin/env python3
"""
検算 A・B・C（2026-10-19）

`octonion_closure_check_20261019.py` と同じフォルダに置いて実行する（代数の定義を読み込む）。
検討メモ `内と外の選択問題_八元数分割と読み出し_20261019.md` の §6-1、§6-4、§1.3 の
「読み」を計算結果に変えるための検算。新しい仮定は置かない。

A. 分け方の不一致と枠の回転
   二人の観測者の時空半分 H'(1) と H'(2) のずれ角 φ を連続に振り、
   三角形の枠ホロノミーを各自の分け方で読んだときの
   「時空の枠の回転」と「時空↔内部の混合」が φ にどう依存するかを測る。
   φ → 0 で 0 に落ちるか、面積比例か。

B. 対になる操作は一つか
   同じ三角形の「一周の閉じなさ」を三通りに測り、同じ量かを見る。
   (a) 枠ホロノミー ‖H − I‖
   (b) 三つの一歩の結合子 ‖(s_ef s_fg)s_ge − s_ef(s_fg s_ge)‖
   (c) 読み出しの一周 ‖u_ef + u_fg + u_ge‖（中心投影の和、平坦なら 0）

C. チャンネルが開く閾値
   内部側の虚方向が N とともに何本開くかを、複数シードで分布として見る
   （N = 6 で 4 本目、N = 9 で 7 本すべて、が乱数 1 回の偶然でないか）。

再現性：乱数シード 20261019。結果は results/checkABC.json と results/checkABC.md、
図は figures/fig3_*.svg（依存なしの手書き SVG）と fig3_*.png（matplotlib）。
"""
import json, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
exec(open(os.path.join(HERE, "octonion_closure_check_20261019.py")).read().split("results = {")[0])

SEED = 20261019
rng = np.random.default_rng(SEED)
OUT = os.path.join(HERE, "results"); FIG = os.path.join(HERE, "figures")
os.makedirs(OUT, exist_ok=True); os.makedirs(FIG, exist_ok=True)
res = {"seed": SEED}

# ---------------------------------------------------------------------------
# 共通の道具
# ---------------------------------------------------------------------------
def unit(o): return o / np.sqrt(onorm2(o))

def frame_holonomy(oe, of, og):
    """比の一歩 s_ef = o_f ō_e を三角形で一周したときの枠の写像（8x8）。"""
    sef, sfg, sge = omul(of, oconj(oe)), omul(og, oconj(of)), omul(oe, oconj(og))
    H = np.zeros((8, 8))
    for k in range(8):
        H[:, k] = omul(sge, omul(sfg, omul(sef, E[k])))
    return H

def rot_angle(M):
    return float(np.max(np.abs(np.angle(np.linalg.eigvals(M)))))

def sph_area(e, f, g):
    def ang(p, q, r):
        tq = q - (p @ q) * p; tr = r - (p @ r) * p
        return np.arccos(np.clip(tq @ tr / np.linalg.norm(tq) / np.linalg.norm(tr), -1, 1))
    return float(ang(e, f, g) + ang(f, g, e) + ang(g, e, f) - np.pi)

def im_unit(mask=None):
    v = rng.normal(size=8); v[0] = 0.0
    if mask is not None: v = v * mask
    return v / np.linalg.norm(v)

def subalgebra_from(A, B):
    """A, B（直交する単位純虚）から四元数部分代数 span(1, A, B', AB') の基底（4x8、行が基底）。"""
    Bp = B - (A @ B) * A; Bp /= np.linalg.norm(Bp)
    C = omul(A, Bp); C[0] = 0.0; C /= np.linalg.norm(C)
    return np.array([E[0], A, Bp, C])

def complement_basis(Bas):
    Q, _ = np.linalg.qr(np.hstack([Bas.T, rng.normal(size=(8, 8))]))
    return Q[:, 4:].T                                    # 4x8

def subspace_angle(Bas1, Bas2):
    """二つの 4 次元部分空間の主角のうち最大のもの（rad）。1 は共通なので 0 が一つ出る。"""
    s = np.linalg.svd(Bas1 @ Bas2.T, compute_uv=False)
    return float(np.max(np.arccos(np.clip(s, -1, 1))))

def misalignment(Bas, Cm):
    """観測者の時空半分 Bas のうち、三角形の補空間 Cm に乗る量（0..3）。"""
    return float(np.linalg.norm(Cm @ Bas.T) ** 2)

# ---------------------------------------------------------------------------
# A. 分け方の不一致と枠の回転
# ---------------------------------------------------------------------------
def studyA(n_phi=25, n_tri=40, eps=0.03):
    """
    三角形は固定の部分代数 H_T（A0, B0 から作る）の中に小さく取る。
    観測者 1 は H'(1) = H_T（整合）。観測者 2 は H_T を角 φ だけ傾けた H'(2)。
    各 φ で、観測者 2 の時空半分に現れる枠の回転と混合を、面積で割って測る。
    """
    rec = []
    for _ in range(n_tri):
        A0 = im_unit(); B0 = im_unit()
        Bas0 = subalgebra_from(A0, B0)
        # 傾ける方向：H_T の虚部の 1 本を、補空間の 1 本へ回す
        Cm0 = complement_basis(Bas0)
        a_in = Bas0[1]                      # H_T の虚部の 1 本
        a_out = Cm0[rng.integers(0, 4)]     # 補空間の 1 本
        oe, of, og = E[0], oexp(eps * A0), oexp(eps * B0)
        H = frame_holonomy(oe, of, og); area = sph_area(oe, of, og)
        for phi in np.linspace(0.0, np.pi / 2, n_phi):
            # 観測者 2 の時空半分：a_in を cos φ a_in + sin φ a_out へ傾けた四元数部分代数
            A2 = np.cos(phi) * a_in + np.sin(phi) * a_out
            A2 = A2 - A2[0] * E[0]; A2 /= np.linalg.norm(A2)
            B2 = B0 - (A2 @ B0) * A2
            if np.linalg.norm(B2) < 1e-8: continue
            B2 /= np.linalg.norm(B2)
            Bas2 = subalgebra_from(A2, B2); Cm2 = complement_basis(Bas2)
            Hqq = Bas2 @ H @ Bas2.T                       # 観測者 2 の時空ブロック
            Hqp = Bas2 @ H @ Cm2.T                        # 時空↔内部の混合
            rec.append({"phi": float(phi),
                        "subspace_angle": subspace_angle(Bas0, Bas2),
                        "misalign": misalignment(Bas2, Cm0),
                        "area": area,
                        "Hqq_dev_over_area": float(np.linalg.norm(Hqq - np.eye(4)) / area),
                        "mix_over_area": float(np.linalg.norm(Hqp) / area),
                        "angle_qq_over_area": rot_angle(Hqq) / area})
    return rec

# φ 依存の形を見るために、面積も振る
def studyA_area(n_eps=9, n_tri=25):
    rec = []
    for _ in range(n_tri):
        A0 = im_unit(); B0 = im_unit()
        Bas0 = subalgebra_from(A0, B0); Cm0 = complement_basis(Bas0)
        a_in = Bas0[1]; a_out = Cm0[rng.integers(0, 4)]
        phi = rng.uniform(0.2, np.pi / 2)
        A2 = np.cos(phi) * a_in + np.sin(phi) * a_out; A2 -= A2[0] * E[0]; A2 /= np.linalg.norm(A2)
        B2 = B0 - (A2 @ B0) * A2
        if np.linalg.norm(B2) < 1e-8: continue
        B2 /= np.linalg.norm(B2)
        Bas2 = subalgebra_from(A2, B2); Cm2 = complement_basis(Bas2)
        for eps in np.logspace(-2.5, -1, n_eps):
            oe, of, og = E[0], oexp(eps * A0), oexp(eps * B0)
            H = frame_holonomy(oe, of, og); area = sph_area(oe, of, og)
            rec.append({"phi": float(phi), "eps": float(eps), "area": area,
                        "Hqq_dev": float(np.linalg.norm(Bas2 @ H @ Bas2.T - np.eye(4))),
                        "mix": float(np.linalg.norm(Bas2 @ H @ Cm2.T))})
    return rec

recA = studyA(); recA2 = studyA_area()
phis = np.array([r["phi"] for r in recA])
hq = np.array([r["Hqq_dev_over_area"] for r in recA])
mx = np.array([r["mix_over_area"] for r in recA])
z = phis < 1e-9
res["A"] = {
    "n_samples": len(recA),
    "Hqq_dev_over_area_at_phi0_max": float(hq[z].max()),
    "mix_over_area_at_phi0_max": float(mx[z].max()),
    "Hqq_dev_over_area_at_phi_max_median": float(np.median(hq[phis > 1.5])),
    "mix_over_area_at_phi_max_median": float(np.median(mx[phis > 1.5])),
}
# 形の同定：mix/area を sin^k φ で合わせる（0 < φ < π/2 の中央値カーブ）
bins = np.linspace(0.02, np.pi / 2, 16)
med_mix, med_hq, ctr = [], [], []
for lo, hi in zip(bins[:-1], bins[1:]):
    m = (phis >= lo) & (phis < hi)
    if m.sum() > 5:
        ctr.append(0.5 * (lo + hi)); med_mix.append(np.median(mx[m])); med_hq.append(np.median(hq[m]))
ctr = np.array(ctr); med_mix = np.array(med_mix); med_hq = np.array(med_hq)
sl_mix = np.polyfit(np.log(np.sin(ctr)), np.log(med_mix), 1)
sl_hq = np.polyfit(np.log(np.sin(ctr)), np.log(med_hq), 1)
res["A"]["mix_vs_sinphi_loglog_slope"] = float(sl_mix[0])
res["A"]["Hqq_vs_sinphi_loglog_slope"] = float(sl_hq[0])
ar = np.array([r["area"] for r in recA2]); dv = np.array([r["Hqq_dev"] for r in recA2]); mi = np.array([r["mix"] for r in recA2])
res["A"]["Hqq_dev_vs_area_loglog_slope"] = float(np.polyfit(np.log(ar), np.log(dv + 1e-300), 1)[0])
res["A"]["mix_vs_area_loglog_slope"] = float(np.polyfit(np.log(ar), np.log(mi + 1e-300), 1)[0])
print("A. 分け方の不一致と枠の回転")
print(f"   φ = 0 での |Hqq−I|/面積 の最大 = {res['A']['Hqq_dev_over_area_at_phi0_max']:.2e}（整合なら 0）")
print(f"   φ = 0 での 混合/面積 の最大   = {res['A']['mix_over_area_at_phi0_max']:.2e}")
print(f"   φ ≈ π/2 での中央値：|Hqq−I|/面積 = {res['A']['Hqq_dev_over_area_at_phi_max_median']:.3f}, 混合/面積 = {res['A']['mix_over_area_at_phi_max_median']:.3f}")
print(f"   sin φ に対する傾き（log-log）：混合 {sl_mix[0]:.3f}, |Hqq−I| {sl_hq[0]:.3f}")
print(f"   面積に対する傾き（φ 固定）：|Hqq−I| {res['A']['Hqq_dev_vs_area_loglog_slope']:.3f}, 混合 {res['A']['mix_vs_area_loglog_slope']:.3f}")

# ---------------------------------------------------------------------------
# B. 対になる操作は一つか
#    三角形 (oe, of, og) について、同じ「一周の閉じなさ」を測る三つの量を比べる。
#    (a) 枠ホロノミー ‖H − I‖（比の一歩 s_ef = o_f ō_e を三つ合成した結果）
#    (b) 三つの一歩の結合子 ‖(s_ef s_fg)s_ge − s_ef(s_fg s_ge)‖（非結合性そのもの）
#    (c) 読み出しの一周 ‖u_ef + u_fg + u_ge‖（平坦なら telescoping で 0 になる量）
#    注意：三点のどれかを 1 に取ると結合子は恒等的に 0 になる（[1,y,z] = 0）。
#          ここでは三点とも 1 から離して取る。
# ---------------------------------------------------------------------------
def readout(oe, of):
    """観測者 e が相手 f を読む：関係 w = ō_e o_f のベクトル部をスカラー部で割る（中心投影）。"""
    w = omul(oconj(oe), of)
    return w[1:] / w[0]

def studyB(n=400):
    rec = []
    for _ in range(n):
        eps = 10 ** rng.uniform(-2.2, -0.9)
        A, B, C = im_unit(), im_unit(), im_unit()
        oe, of, og = oexp(eps * A), oexp(eps * B), oexp(eps * C)   # 三点とも 1 から離す
        sef, sfg, sge = omul(of, oconj(oe)), omul(og, oconj(of)), omul(oe, oconj(og))
        H = frame_holonomy(oe, of, og); area = sph_area(oe, of, og)
        holo = float(np.linalg.norm(H - np.eye(8)))
        assoc = float(np.linalg.norm(omul(omul(sef, sfg), sge) - omul(sef, omul(sfg, sge))))
        loop = float(np.linalg.norm(readout(oe, of) + readout(of, og) + readout(og, oe)))
        # 三角形自身の部分代数と、固定した観測者の分け方（span(e0..e3)）の不整合
        A2 = sef - sef[0] * E[0]; B2 = sfg - sfg[0] * E[0]
        if np.linalg.norm(A2) < 1e-12 or np.linalg.norm(B2) < 1e-12: continue
        Bas_T = subalgebra_from(A2 / np.linalg.norm(A2), B2 / np.linalg.norm(B2))
        Cm_T = complement_basis(Bas_T)
        Pq = np.zeros((8, 4)); Pq[:4, :4] = np.eye(4)
        mis = float(np.linalg.norm(Cm_T @ Pq) ** 2)
        rec.append({"eps": eps, "area": area, "holonomy": holo, "associator": assoc,
                    "readout_loop": loop, "misalign": mis,
                    "holonomy_qq": float(np.linalg.norm(H[:4, :4] - np.eye(4)))})
    return rec

recB = studyB()
def corr(a, b):
    a = np.array(a, float); b = np.array(b, float)
    if a.std() < 1e-300 or b.std() < 1e-300: return float("nan")
    return float(np.corrcoef(a, b)[0, 1])
def slope(x, y):
    x = np.array(x, float); y = np.array(y, float)
    m = (x > 0) & (y > 0)
    return float(np.polyfit(np.log(x[m]), np.log(y[m]), 1)[0]) if m.sum() > 5 else float("nan")
hol = [r["holonomy"] for r in recB]; asoc = [r["associator"] for r in recB]
loop = [r["readout_loop"] for r in recB]; mis = [r["misalign"] for r in recB]
ar = [r["area"] for r in recB]; hqq = [r["holonomy_qq"] for r in recB]
res["B"] = {
    "n_samples": len(recB),
    "corr_log_holonomy_associator": corr(np.log(hol), np.log(asoc)),
    "corr_log_holonomy_readout_loop": corr(np.log(hol), np.log(loop)),
    "corr_log_associator_readout_loop": corr(np.log(asoc), np.log(loop)),
    "ratio_associator_over_holonomy_median": float(np.median(np.array(asoc) / np.array(hol))),
    "ratio_readout_loop_over_holonomy_median": float(np.median(np.array(loop) / np.array(hol))),
    "ratio_associator_over_holonomy_iqr": [float(np.percentile(np.array(asoc) / np.array(hol), 25)),
                                           float(np.percentile(np.array(asoc) / np.array(hol), 75))],
    "holonomy_vs_area_loglog_slope": slope(ar, hol),
    "associator_vs_area_loglog_slope": slope(ar, asoc),
    "readout_loop_vs_area_loglog_slope": slope(ar, loop),
    "corr_misalign_holonomy_qq_over_area": corr(mis, np.array(hqq) / np.array(ar)),
}
print("\nB. 対になる操作（同じ一周の閉じなさを三通りに測る）")
print(f"   log-log 相関：ホロノミー vs 結合子 = {res['B']['corr_log_holonomy_associator']:+.4f}")
print(f"                 ホロノミー vs 読み出し一周 = {res['B']['corr_log_holonomy_readout_loop']:+.4f}")
print(f"                 結合子 vs 読み出し一周 = {res['B']['corr_log_associator_readout_loop']:+.4f}")
print(f"   比の中央値：結合子/ホロノミー = {res['B']['ratio_associator_over_holonomy_median']:.4f} "
      f"(四分位 {res['B']['ratio_associator_over_holonomy_iqr'][0]:.4f}–{res['B']['ratio_associator_over_holonomy_iqr'][1]:.4f})")
print(f"               読み出し一周/ホロノミー = {res['B']['ratio_readout_loop_over_holonomy_median']:.4f}")
print(f"   面積に対する傾き：ホロノミー {res['B']['holonomy_vs_area_loglog_slope']:.3f}, "
      f"結合子 {res['B']['associator_vs_area_loglog_slope']:.3f}, 読み出し一周 {res['B']['readout_loop_vs_area_loglog_slope']:.3f}")
print(f"   不整合 vs 時空枠の回転/面積 の相関 = {res['B']['corr_misalign_holonomy_qq_over_area']:+.3f}")

# ---------------------------------------------------------------------------
# C. チャンネルが開く閾値
# ---------------------------------------------------------------------------
def gen_closure_O(N, r):
    one = np.ones(N) / np.sqrt(N)
    x0 = r.normal(size=N); x0 -= one * (one @ x0)
    P = np.eye(N) - np.outer(one, one) - np.outer(x0, x0) / (x0 @ x0)
    X = np.array([P @ r.normal(size=N) for _ in range(7)])
    X *= np.sqrt((x0 @ x0) / np.sum(X * X))
    return np.column_stack([x0, X.T])

def rank_(M, tol=1e-9):
    s = np.linalg.svd(M, compute_uv=False)
    return int(np.sum(s > tol * s[0]))

def studyC(Ns=range(3, 13), trials=30):
    out = []
    for N in Ns:
        opens, shapes = [], []
        for t in range(trials):
            r = np.random.default_rng(SEED + 1000 * N + t)
            S = gen_closure_O(N, r)
            opens.append(rank_(S[:, 1:].T))
            # 形の空間の次元
            J = np.zeros((16, N * 8))
            for n in range(N):
                for k in range(8):
                    J[:8, n*8 + k] = E[k]; J[8:, n*8 + k] = 2.0 * ojordan(S[n], E[k])
            vecs = [S.reshape(-1)]
            for i in range(1, 8):
                for j in range(i + 1, 8):
                    T = np.zeros_like(S); T[:, i] = -S[:, j]; T[:, j] = S[:, i]; vecs.append(T.reshape(-1))
            shapes.append(N * 8 - rank_(J) - rank_(np.array(vecs)))
        out.append({"N": N, "open_dirs_min": int(min(opens)), "open_dirs_max": int(max(opens)),
                    "open_dirs_mean": float(np.mean(opens)),
                    "shape_dim_min": int(min(shapes)), "shape_dim_max": int(max(shapes)),
                    "shape_dim_mode": int(np.bincount(shapes).argmax()), "trials": trials})
    return out
recC = studyC()
res["C"] = recC
print("\nC. チャンネルが開く閾値（各 N で 30 シード）")
print("   N   開く虚方向(min/mean/max)   形の空間(min/mode/max)")
for r in recC:
    print(f"  {r['N']:>2}      {r['open_dirs_min']}/{r['open_dirs_mean']:.2f}/{r['open_dirs_max']}"
          f"            {r['shape_dim_min']}/{r['shape_dim_mode']}/{r['shape_dim_max']}")

# ---------------------------------------------------------------------------
# 保存
# ---------------------------------------------------------------------------
res["_recA_curve"] = {"phi": ctr.tolist(), "mix_over_area_median": med_mix.tolist(),
                      "Hqq_over_area_median": med_hq.tolist()}
with open(os.path.join(OUT, "checkABC.json"), "w", encoding="utf-8") as f:
    json.dump(res, f, ensure_ascii=False, indent=1)

# --- 依存なしの SVG
def svg_panels(path, panels, W=1080, Hh=330):
    n = len(panels); pw = W // n; m = dict(l=62, r=14, t=34, b=46)
    def tf(v, lo, hi, log):
        if log: v, lo, hi = np.log10(max(v, 1e-300)), np.log10(lo), np.log10(hi)
        return (v - lo) / (hi - lo)
    o = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{Hh}" viewBox="0 0 {W} {Hh}" font-family="Helvetica, Arial, sans-serif" font-size="11">',
         f'<rect width="{W}" height="{Hh}" fill="white"/>']
    for i, p in enumerate(panels):
        x0 = i*pw + m["l"]; x1 = (i+1)*pw - m["r"]; y0 = m["t"]; y1 = Hh - m["b"]
        xl, yl = p["xlim"], p["ylim"]
        o.append(f'<rect x="{x0}" y="{y0}" width="{x1-x0}" height="{y1-y0}" fill="none" stroke="#888"/>')
        o.append(f'<text x="{(x0+x1)/2:.0f}" y="{y0-16}" text-anchor="middle" font-size="12">{p["title"]}</text>')
        o.append(f'<text x="{(x0+x1)/2:.0f}" y="{Hh-12}" text-anchor="middle">{p["xlabel"]}</text>')
        o.append(f'<text transform="translate({i*pw+16},{(y0+y1)/2:.0f}) rotate(-90)" text-anchor="middle">{p["ylabel"]}</text>')
        xt = ([10.0**k for k in range(int(np.floor(np.log10(xl[0]))), int(np.ceil(np.log10(xl[1])))+1)]
              if p.get("xlog") else list(np.linspace(xl[0], xl[1], 5)))
        for v in xt:
            if v < xl[0]-1e-12 or v > xl[1]+1e-12: continue
            xx = x0 + tf(v, xl[0], xl[1], p.get("xlog"))*(x1-x0)
            lab = (f'1e{int(np.log10(v))}' if p.get("xlog") else f'{v:g}')
            o.append(f'<line x1="{xx:.1f}" y1="{y1}" x2="{xx:.1f}" y2="{y1+4}" stroke="#888"/><text x="{xx:.1f}" y="{y1+17}" text-anchor="middle" font-size="10">{lab}</text>')
        yt = ([10.0**k for k in range(int(np.floor(np.log10(yl[0]))), int(np.ceil(np.log10(yl[1])))+1)]
              if p.get("ylog") else list(np.linspace(yl[0], yl[1], 5)))
        for v in yt:
            if v < yl[0]-1e-12 or v > yl[1]+1e-12: continue
            yy = y1 - tf(v, yl[0], yl[1], p.get("ylog"))*(y1-y0)
            lab = (f'1e{int(np.log10(v))}' if p.get("ylog") else f'{v:g}')
            o.append(f'<line x1="{x0-4}" y1="{yy:.1f}" x2="{x0}" y2="{yy:.1f}" stroke="#888"/><text x="{x0-6}" y="{yy+3:.1f}" text-anchor="end" font-size="10">{lab}</text>')
        ly = y0 + 14
        for (xs, ys, color, kind, label) in p["series"]:
            pts = [(x0 + tf(x, xl[0], xl[1], p.get("xlog"))*(x1-x0),
                    y1 - tf(min(max(y, yl[0]), yl[1]), yl[0], yl[1], p.get("ylog"))*(y1-y0)) for x, y in zip(xs, ys)]
            if kind in ("line", "dash"):
                d = " ".join(f'{a:.1f},{b:.1f}' for a, b in pts)
                dash = ' stroke-dasharray="5 4"' if kind == "dash" else ''
                o.append(f'<polyline points="{d}" fill="none" stroke="{color}" stroke-width="1.6"{dash}/>')
            if kind in ("line", "dots"):
                r = 3 if kind == "line" else 2.2
                o.append(f'<g fill="{color}" fill-opacity="0.8">' + "".join(f'<circle cx="{a:.1f}" cy="{b:.1f}" r="{r}"/>' for a, b in pts) + '</g>')
            if label:
                o.append(f'<circle cx="{x0+10}" cy="{ly}" r="3" fill="{color}"/><text x="{x0+17}" y="{ly+4}" font-size="10">{label}</text>'); ly += 14
    o.append('</svg>')
    open(path, "w", encoding="utf-8").write("\n".join(o))

sub = rng.choice(len(recA), min(400, len(recA)), replace=False)
svg_panels(os.path.join(FIG, "fig3_misalignment_pairing_compact.svg"), [
    dict(title="A: split misalignment phi vs frame rotation (per unit area)",
         xlabel="phi  (tilt of observer 2's spacetime half)", ylabel="per unit area",
         xlim=(0, 1.6), ylim=(0, 6),
         series=[([recA[i]["phi"] for i in sub], [recA[i]["Hqq_dev_over_area"] for i in sub], "#d1247a", "dots", "|Hqq - I| / area"),
                 ([recA[i]["phi"] for i in sub], [recA[i]["mix_over_area"] for i in sub], "#2457e6", "dots", "spacetime-internal mixing / area"),
                 (ctr.tolist(), med_mix.tolist(), "#111111", "line", "median of mixing")]),
    dict(title="B: associator and readout loop vs frame holonomy", xlabel="frame holonomy |H - I|", ylabel="other measures of the same loop",
         xlog=True, ylog=True, xlim=(1e-4, 1e0), ylim=(1e-11, 1e-1),
         series=[([r["holonomy"] for r in recB], [r["associator"] for r in recB], "#d9800f", "dots", "associator of the three steps"),
                 ([r["holonomy"] for r in recB], [r["readout_loop"] for r in recB], "#2457e6", "dots", "readout loop |u_ef+u_fg+u_ge|")]),
    dict(title="C: open internal directions vs N (30 seeds each)", xlabel="N", ylabel="open imaginary directions",
         xlim=(3, 12), ylim=(0, 8),
         series=[([r["N"] for r in recC], [r["open_dirs_mean"] for r in recC], "#2457e6", "line", "mean"),
                 ([r["N"] for r in recC], [min(7, r["N"]-2) for r in recC], "#888888", "dash", "min(7, N-2)")])])

# --- matplotlib 版
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams["font.family"] = ["DejaVu Sans"]
fig, axs = plt.subplots(1, 3, figsize=(13, 3.8))
axs[0].plot(phis, hq, "o", ms=2.5, color="#d1247a", label="|Hqq-I| / area")
axs[0].plot(phis, mx, "o", ms=2.5, color="#2457e6", alpha=0.6, label="mixing / area")
axs[0].plot(ctr, med_mix, "-", color="k", lw=1.5, label="median (mixing)")
axs[0].set_xlabel("phi (misalignment of the two splits)"); axs[0].set_ylabel("per unit area")
axs[0].set_title("A: frame rotation seen by observer 2", fontsize=9); axs[0].legend(fontsize=7)
axs[1].loglog([r["holonomy"] for r in recB], [r["associator"] for r in recB], "o", ms=2.5, color="#d9800f", label="associator")
axs[1].loglog([r["holonomy"] for r in recB], [r["readout_loop"] for r in recB], "o", ms=2.5, color="#2457e6", alpha=0.6, label="readout loop")
axs[1].set_xlabel("frame holonomy |H - I|"); axs[1].set_ylabel("same loop, other measures")
axs[1].set_title(f"B: log-log corr assoc {res['B']['corr_log_holonomy_associator']:+.3f}", fontsize=9); axs[1].legend(fontsize=7)
axs[2].plot([r["N"] for r in recC], [r["open_dirs_mean"] for r in recC], "o-", color="#2457e6", label="mean open dirs")
axs[2].plot([r["N"] for r in recC], [min(7, r["N"]-2) for r in recC], "--", color="#888", label="min(7, N-2)")
axs[2].set_xlabel("N"); axs[2].set_ylabel("open imaginary directions")
axs[2].set_title("C: channel opening", fontsize=9); axs[2].legend(fontsize=7)
for a in axs: a.grid(alpha=0.3)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig3_misalignment_pairing.png"), dpi=160)
print("\nsaved:", os.path.join(OUT, "checkABC.json"), FIG)
