#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""90度整数K床の parity 二部構造による解析的一本化——検証（2026-09-12、木原レビュー反映）。

主張（論文1 D^{-1}KD=iW を使う）:
  種ラベル k_e=(i+j) mod4 について、その parity p_e=(i+j) mod2 だけで
    sin^2[(pi/2)(k_f-k_e)] = 1  (p_e != p_f),  0  (p_e == p_f)
  ゆえ W_ef = A_ef · 1[p_e != p_f]（L(K_N) を頂点番号和の偶奇で二分した二部グラフ）。
  したがって:
   (1) 反復（make_parent 不使用の §3 ラベル反復）で得た床の W は、種 parity の W と機械精度一致
       → §3 反復は床の構成には不要（parity は種のまま）。
   (2) z = D_seed · r_Perron（D_seed=diag(e^{i(pi/2)k_seed}), r=W の Perron ベクトル>0）が
       自己無撞着相対平衡（one_step 残差 ~機械精度）で、反復床とグローバル位相を除き一致。
   (3) 縮約 Perron 問題（偶数=2x2、奇数=3x3）から sigma^2・各辺振幅^2 の閉形式が任意 N で厳密に出る。
   (4) 二部 W の全スペクトル（BB^T）が閉形式（多重度つき）で任意 N。
   (5) rank K = rank W = 2(N-1)、核次元 = M-2(N-1) = (N-1)(N-4)/2（偶 N>=6 / 奇 N>=7）。
   (6) 未決の「振幅クラスの辺割当」= parity クラス、奇数残差重心も閉形式。

読出しのみ・物理無変更・新規走行なし。adjacency/one_step は物理正本（SHA照合）から抽出。
入力: ../parents_90deg_floor_analytic/ の既存 38 床（反復床）。
出力: results/parity_bipartite_summary.csv, results/実行ログ_20260912.log。
"""
import hashlib
import os
import sys

import numpy as np
import sympy as sp

BASE = os.path.dirname(os.path.abspath(__file__))
ORIG = os.path.join(BASE, '..', '..', 'N3_N40_stage123_sweep_20260905', 'run_N3_N40_stage123_v1.py')
PROG_SHA = '1abf2353fee2e4f56f05e7a6f149fd086885136beb61ab571b48a56b09691567'
FLOORS = os.path.join(BASE, '..', 'parents_90deg_floor_analytic')
RESULTS = os.path.join(BASE, 'results')
os.makedirs(RESULTS, exist_ok=True)

# --- 物理正本から adjacency/one_step を抽出（SHA照合・物理無変更） ---
src = open(ORIG, encoding='utf-8').read()
assert hashlib.sha256(src.encode()).hexdigest() == PROG_SHA, '正本SHA不一致'
g = {'__name__': 'canonical_funcs', '__file__': ORIG}
exec(compile(src.split('rows=[]; summaries=[]')[0], ORIG, 'exec'), g)
adjacency, one_step = g['adjacency'], g['one_step']

NS = list(range(3, 41))
log_lines = []


def log(*a):
    s = ' '.join(str(x) for x in a)
    print(s)
    log_lines.append(s)


def releq_residual(z, A, N):
    w = one_step(z, A, N)
    ip = np.vdot(z, w)
    return float(np.linalg.norm(w - np.exp(1j * np.angle(ip)) * z) / np.linalg.norm(z))


def perron(W):
    ev, V = np.linalg.eigh(W)
    r = V[:, -1]
    if r.sum() < 0:
        r = -r
    return float(ev[-1]), r  # (rho(W), Perron vector)


def amp2_closed_form(N, i, j):
    """辺 e={i,j} の |z_e|^2 の閉形式（parity クラス割当つき）。"""
    parity = (i + j) % 2
    if N % 2 == 0:
        # 異偶間(parity=1): 2/N^2、 同偶奇内(parity=0): 2/(N(N-2))
        return 2.0 / N**2 if parity == 1 else 2.0 / (N * (N - 2))
    else:
        D = N**2 - 2 * N - 1
        if parity == 1:            # 異偶間
            return 2.0 / ((N - 1) * (N + 1))
        elif (i % 2 == 0):         # 同偶奇内・偶側内部（i,j とも偶）
            return 2.0 * (N - 1) / ((N + 1) * D)
        else:                      # 同偶奇内・奇側内部（i,j とも奇）
            return 2.0 * (N + 1) / ((N - 1) * D)


def centroid_closed_form(N):
    D = N**2 - 2 * N - 1
    r = N % 4
    if r == 0:
        return 0.0
    if r == 2:
        return np.sqrt(2.0) / N
    if r == 1:
        return np.sqrt((N - 1) / (2.0 * (N + 1) * D))
    return np.sqrt((N + 1) / (2.0 * (N - 1) * D))  # r == 3


def spectrum_closed_form(N):
    """二部 W の非零固有値^2（BB^T の非零固有値）とその多重度リスト。偶 N>=6 / 奇 N>=7 で有効。"""
    if N % 2 == 0:
        return [(N * (N - 2), 1), (N * (N - 4) / 4.0, N - 2)]
    else:
        return [(N**2 - 2 * N - 1, 1),
                ((N - 1) * (N - 3) / 4.0, (N - 1) // 2),
                ((N - 5) * (N + 1) / 4.0, (N - 3) // 2)]


# =========================== 数値検証（全 N=3..40） ===========================
log('# 90度整数K床 parity 二部構造 解析的一本化 検証  日付 2026-09-12')
log('# 物理正本 SHA 一致:', PROG_SHA)
log('')
header = ('N  | ‖W_iter−W_parity‖  releq(z_dir)  ‖床_iter−床_dir‖(gauge) | '
         'ω=√σ² 一致 r>0 | 振幅閉形式 max誤差 | 重心 |Σz| 実測/閉形式 | '
         'rankW 2(N-1) 核 (N-1)(N-4)/2 | スペクトル閉形式')
log(header)
rows = [['N', 'M', 'W_iter_vs_parity', 'releq_z_direct', 'amp_match', 's_is_pm1',
         'floor_iter_vs_direct_gauge', 'closure_all_nonzero_modes', 'xnorm2_minus_ynorm2',
         'n_components', 'sigma2', 'sigma_match', 'r_positive', 'amp2_max_err',
         'centroid_meas', 'centroid_cf', 'centroid_err', 'rankW', 'rankW_pred', 'kernel',
         'kernel_pred', 'spectrum_match']]
all_ok = True

for N in NS:
    A = adjacency(N)
    ea, eb = np.triu_indices(N, k=1)
    M = len(ea)
    kseed = ((ea + eb) % 4).astype(int)
    p = (ea + eb) % 2                       # parity ラベル

    # parity 二部 W
    Wpar = A * (p[:, None] != p[None, :]).astype(float)

    # 反復床（既存）を読み、その W_floor を作る（gauge 不変）
    fp = os.path.join(FLOORS, f'parent_90deg_floor_N{N:05d}_analytic.npz')
    d = np.load(fp)
    Zit = d['Z0']
    phi = np.angle(Zit)
    psi = phi[None, :] - phi[:, None]
    Wfloor = A * np.sin(psi)**2
    w_iter_vs_parity = float(np.max(np.abs(Wfloor - Wpar)))

    # (2) 種ラベル直接構成 z_direct = D_seed r_Perron
    rho, r = perron(Wpar)
    r = r / np.linalg.norm(r)
    Dseed = np.exp(1j * (np.pi / 2) * kseed)
    z_dir = Dseed * r
    releq = releq_residual(z_dir, A, N)
    # 反復床と direct 床の一致。床は符号ゲージ z_e->-z_e（ラベル k->k+2, parity 不変,
    # K->SKS で Kz=iωz・W 不変）の離散自由度を持つので、グローバル位相 α と 辺ごと符号 S=±1
    # を許して整合させる（両者が同一の床＝σ_max Perron 床であることの gauge 不変判定）。
    # まず辺ごと振幅一致（Perron ベクトルの一意性）。
    amp_match = float(np.max(np.abs(np.abs(Zit) - np.abs(z_dir))))
    # グローバル位相 α を最頻比から決め、残差を ±1（符号ゲージ）に丸めて残差ノルム。
    ratio = Zit * np.conj(z_dir) / (np.abs(Zit) * np.abs(z_dir) + 1e-300)
    alpha = np.angle(np.sum(ratio))                       # 支配的グローバル位相
    s = np.real(ratio * np.exp(-1j * alpha))              # ≈ ±1（符号ゲージ）
    S = np.sign(s)
    s_is_pm1 = float(np.max(np.abs(np.abs(s) - 1.0)))     # |s|=1 か（90°格子の帰結）
    floor_iter_vs_direct = float(np.linalg.norm(Zit - np.exp(1j * alpha) * S * z_dir))

    # (3) σ² と ω 一致
    sigma2 = N * (N - 2) if N % 2 == 0 else N**2 - 2 * N - 1
    sigma_match = abs(rho - np.sqrt(sigma2)) < 1e-10
    r_positive = bool(np.all(r > 1e-12) or np.all(r < -1e-12))

    # 振幅閉形式（parity クラス割当）
    amp2_meas = np.abs(z_dir)**2
    amp2_cf = np.array([amp2_closed_form(N, int(i), int(j)) for i, j in zip(ea, eb)])
    amp2_max_err = float(np.max(np.abs(amp2_meas - amp2_cf)))

    # 重心閉形式
    centroid_meas = float(abs(z_dir.sum()))
    centroid_cf = centroid_closed_form(N)
    centroid_err = abs(centroid_meas - centroid_cf)

    # (7) 二乗閉塞は二部構造の帰結：全非零固有モードで Σz²=0（Perron に限らない）。
    #   z_e=e^{iπk_e/2}r_e より z_e²=(-1)^{p_e}r_e²、Σz²=Σ_{p=0}r²−Σ_{p=1}r²。
    #   二部 W の非零固有モード r=(x,y) は By=σx,Bᵀx=σy ⇒ σ‖x‖²=xᵀBy=(Bᵀx)ᵀy=σ‖y‖²
    #   ⇒ ‖x‖²=‖y‖²（σ≠0）ゆえ Σz²=0。全非零固有ベクトルで数値確認。
    evW_all, VW = np.linalg.eigh(Wpar)
    closure_all_max = 0.0       # 全非零モードの |Σz²|
    xy_balance_max = 0.0        # 全非零モードの |‖x‖²−‖y‖²|
    for j in range(len(evW_all)):
        if abs(evW_all[j]) > 1e-9:
            rj = VW[:, j]
            zj = Dseed * rj
            closure_all_max = max(closure_all_max, float(abs(np.sum(zj**2))))
            xnorm2 = float(np.sum(rj[p == 0]**2)); ynorm2 = float(np.sum(rj[p == 1]**2))
            xy_balance_max = max(xy_balance_max, abs(xnorm2 - ynorm2))

    # (8) 連結性：W グラフの連結成分数（N>=3 で 1 を一般証明、数値で確認）
    from scipy.sparse.csgraph import connected_components
    from scipy.sparse import csr_matrix
    ncomp = int(connected_components(csr_matrix(Wpar > 0), directed=False)[0])

    # (5) rank / 核
    evW = np.linalg.eigvalsh(Wpar)
    rankW = int(np.sum(np.abs(evW) > 1e-9))
    rankW_pred = 2 * (N - 1)
    kernel = M - rankW
    kernel_pred = (N - 1) * (N - 4) // 2
    valid_rank = (N % 2 == 0 and N >= 6) or (N % 2 == 1 and N >= 7)

    # (4) スペクトル閉形式（|固有値|² を降順、非零のみ、多重度つき）
    valid_spec = valid_rank
    spec_match = True
    lam2 = np.sort(np.abs(evW)**2)[::-1]
    lam2_nonzero = lam2[lam2 > 1e-7]
    cf = spectrum_closed_form(N)
    # 期待多重度リスト（各固有値は ±ペアで 2 回現れる）
    exp_vals = []
    for val, mult in cf:
        exp_vals += [val] * (2 * mult)
    exp_vals = np.sort(np.array(exp_vals, dtype=float))[::-1]
    if valid_spec:
        if len(exp_vals) == len(lam2_nonzero):
            spec_match = bool(np.max(np.abs(exp_vals - lam2_nonzero)) < 1e-6)
        else:
            spec_match = False

    ok = (w_iter_vs_parity < 1e-9 and releq < 1e-9 and amp_match < 1e-9
          and s_is_pm1 < 1e-9 and floor_iter_vs_direct < 1e-9
          and sigma_match and r_positive and amp2_max_err < 1e-9 and centroid_err < 1e-9
          and closure_all_max < 1e-9 and xy_balance_max < 1e-9 and ncomp == 1
          and (not valid_rank or (rankW == rankW_pred and kernel == kernel_pred))
          and (not valid_spec or spec_match))
    all_ok = all_ok and ok

    log(f"{N:>2} | {w_iter_vs_parity:.2e}  releq {releq:.1e}  振幅一致 {amp_match:.1e}  "
        f"床(gauge) {floor_iter_vs_direct:.1e} | Σz²全モード {closure_all_max:.1e}  "
        f"|‖x‖²−‖y‖²| {xy_balance_max:.1e}  連結成分 {ncomp} | "
        f"振幅CF {amp2_max_err:.1e} | |Σz| {centroid_meas:.6f}/{centroid_cf:.6f} | "
        f"rankW {rankW}={rankW_pred} 核 {kernel}={kernel_pred}"
        f"{'(N小)' if not valid_rank else ''} | spec {str(spec_match) if valid_spec else 'N小'} | {'OK' if ok else 'NG'}")
    rows.append([N, M, f'{w_iter_vs_parity:.3e}', f'{releq:.3e}', f'{amp_match:.3e}',
                 f'{s_is_pm1:.3e}', f'{floor_iter_vs_direct:.3e}', f'{closure_all_max:.3e}',
                 f'{xy_balance_max:.3e}', ncomp,
                 sigma2, sigma_match, r_positive, f'{amp2_max_err:.3e}', f'{centroid_meas:.9f}',
                 f'{centroid_cf:.9f}', f'{centroid_err:.2e}', rankW, rankW_pred, kernel, kernel_pred,
                 (spec_match if valid_spec else 'N_small')])

# CSV 保存
import csv
with open(os.path.join(RESULTS, 'parity_bipartite_summary.csv'), 'w', newline='', encoding='utf-8') as f:
    csv.writer(f).writerows(rows)

# =========================== 記号検証（任意 N を sympy で厳密） ===========================
log('')
log('# 記号検証（sympy, 任意 N）')
n = sp.symbols('n', positive=True, integer=True)

# 偶数 N: 縮約 2x2  ω x = (N-2) y, ω y = N x  → ω² = N(N-2)
Meven = sp.Matrix([[0, n - 2], [n, 0]])
ev_even = list(Meven.eigenvals().keys())
sig2_even = sp.simplify([e**2 for e in ev_even][0])
log('偶数 縮約2x2 固有値²  =', sig2_even, ' 期待 N(N-2) →',
    sp.simplify(sig2_even - n * (n - 2)) == 0)
# 正規化で振幅² (x²,y²)
x, y = sp.symbols('x y', positive=True)
w_even = sp.sqrt(n * (n - 2))
y_of_x = sp.sqrt((n - 2) / n)                       # x/y = sqrt((N-2)/N) → x = y*sqrt((N-2)/N)
# 本数: x² は N²/4 本、y² は N(N-2)/4 本。 Σ = 1
xx = sp.symbols('xx', positive=True)               # xx = x²
yy = xx * n / (n - 2)                                # y² = x² * N/(N-2)  （x²/y²=(N-2)/N）
sol = sp.solve(sp.Eq(n**2 / 4 * xx + n * (n - 2) / 4 * yy, 1), xx)[0]
x2 = sp.simplify(sol)
y2 = sp.simplify(sol * n / (n - 2))
log('偶数 振幅² 異偶間 x² =', x2, ' 期待 2/N² →', sp.simplify(x2 - 2 / n**2) == 0)
log('偶数 振幅² 同偶奇 y² =', y2, ' 期待 2/(N(N-2)) →', sp.simplify(y2 - 2 / (n * (n - 2))) == 0)

# 奇数 N: 縮約 3x3  ω x=(N-1)/2 y+(N-3)/2 z, ω y=(N-1)x, ω z=(N+1)x → ω²=N²-2N-1
Modd = sp.Matrix([[0, (n - 1) / 2, (n - 3) / 2], [n - 1, 0, 0], [n + 1, 0, 0]])
poly = sp.factor(Modd.charpoly(sp.symbols('L')).as_expr())
log('奇数 縮約3x3 特性多項式 =', poly)
sig2_odd = n**2 - 2 * n - 1
# ω² が根であることを確認: charpoly(ω)·… 直接 ω²=N²-2N-1 を検算
L = sp.symbols('L')
cp = Modd.charpoly(L).as_expr()
check_odd = sp.simplify(cp.subs(L, sp.sqrt(sig2_odd)))
log('奇数  charpoly(√(N²-2N-1)) =', sp.simplify(check_odd), ' → 0 ?', sp.simplify(check_odd) == 0)
# 奇数 振幅²
xx = sp.symbols('xx', positive=True)               # x² (異偶間)
Dodd = n**2 - 2 * n - 1
yy = (n - 1)**2 / Dodd * xx                          # y²=(N-1)²/D x²
zz = (n + 1)**2 / Dodd * xx                          # z²=(N+1)²/D x²
cnt_x = (n**2 - 1) / 4
cnt_y = (n**2 - 1) / 8
cnt_z = (n - 1) * (n - 3) / 8
sol = sp.solve(sp.Eq(cnt_x * xx + cnt_y * yy + cnt_z * zz, 1), xx)[0]
x2o = sp.simplify(sol)
y2o = sp.simplify(sol * (n - 1)**2 / Dodd)
z2o = sp.simplify(sol * (n + 1)**2 / Dodd)
log('奇数 振幅² 異偶間  x² =', x2o, ' 期待 2/((N-1)(N+1)) →',
    sp.simplify(x2o - 2 / ((n - 1) * (n + 1))) == 0)
log('奇数 振幅² 偶側内部 y² =', y2o, ' 期待 2(N-1)/((N+1)D) →',
    sp.simplify(y2o - 2 * (n - 1) / ((n + 1) * Dodd)) == 0)
log('奇数 振幅² 奇側内部 z² =', z2o, ' 期待 2(N+1)/((N-1)D) →',
    sp.simplify(z2o - 2 * (n + 1) / ((n - 1) * Dodd)) == 0)

log('')
log('=== 総合判定:', 'ALL OK' if all_ok else 'NG あり', '===')

with open(os.path.join(RESULTS, '実行ログ_20260912.log'), 'w', encoding='utf-8') as f:
    f.write('\n'.join(log_lines) + '\n')

sys.exit(0 if all_ok else 1)
