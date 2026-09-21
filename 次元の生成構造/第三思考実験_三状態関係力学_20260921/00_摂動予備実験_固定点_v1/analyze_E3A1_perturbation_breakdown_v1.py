#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E3-A1 fixed-point analysis: breakdown scan for the hierarchical two-Kepler perturbation.

No new force/readout is added.  We keep the same reciprocal linear coupling used
in E3-A0-v2 and ask only where the perturbative *hierarchical interpretation*
ceases to be accurate.

This distinguishes:
  (1) validity of the O(eps^2) eigenvalue shift,
  (2) loss of an elliptic/closed outer normal mode (lambda_out = +2),
  (3) strong mixing of inner/outer hierarchy labels.

The coupled linear system itself remains exactly diagonalizable for all eps, so
"breakdown" here does not mean mathematical divergence of the exact model.
"""
import csv, math
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
N_IN = 31
N_OUT = 127
TAU_IN = 2*math.cos(2*math.pi/N_IN)
TAU_OUT = 2*math.cos(2*math.pi/N_OUT)
DELTA = TAU_OUT - TAU_IN
MEAN = 0.5*(TAU_OUT + TAU_IN)


def exact(eps):
    rad = math.sqrt((DELTA/2)**2 + eps**2)
    lam_in = MEAN - rad
    lam_out = MEAN + rad
    shift = lam_out - TAU_OUT
    shift2 = eps**2/DELTA if DELTA != 0 else math.nan
    rel = abs(shift2/shift - 1.0) if eps > 0 and shift != 0 else 0.0
    theta = 0.5*math.atan2(2*eps, DELTA)
    leak = math.sin(theta)**2
    if abs(lam_out) < 2 - 1e-12:
        kind = 'elliptic'
    elif abs(abs(lam_out)-2) <= 1e-12:
        kind = 'parabolic_boundary'
    else:
        kind = 'hyperbolic'
    return lam_in, lam_out, shift, shift2, rel, theta, leak, kind


def eps_for_relerr(r):
    # rel error of second-order approximation is (sqrt(1+x^2)-1)/2, x=2eps/DELTA
    x = math.sqrt((1+2*r)**2 - 1)
    return 0.5*DELTA*x


def main():
    eps_ell = math.sqrt((2-TAU_IN)*(2-TAU_OUT))
    eps_mix = DELTA/2
    eps_list = [
        0.0, 5e-4, 1e-3, 2e-3, 3e-3, 4e-3, 6e-3, 8e-3,
        eps_for_relerr(0.05), 9e-3, eps_ell, 1.05*eps_ell, 1.2e-2,
        eps_for_relerr(0.10), 1.5e-2, eps_mix, 2.5e-2, 3e-2, 5e-2, 1e-1
    ]
    eps_list = sorted(set(float(x) for x in eps_list))
    rows=[]
    for eps in eps_list:
        li,lo,sh,sh2,re,th,leak,kind=exact(eps)
        rows.append({
            'eps':eps,'eps_over_delta':eps/DELTA,'lambda_in':li,'lambda_out':lo,
            'outer_kind':kind,'tau_out_shift_exact':sh,'tau_out_shift_Oeps2':sh2,
            'shift_Oeps2_rel_error':re,'mix_angle_deg':math.degrees(th),
            'inner_weight_in_outer_mode':leak
        })

    csv_path=HERE/'E3A1_perturbation_breakdown_scan_v1.csv'
    with open(csv_path,'w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0].keys()));w.writeheader();w.writerows(rows)

    thresholds={
        'tau_in':TAU_IN,'tau_out':TAU_OUT,'delta_tau':DELTA,
        'eps_1pct_shift_error':eps_for_relerr(0.01),
        'eps_5pct_shift_error':eps_for_relerr(0.05),
        'eps_10pct_shift_error':eps_for_relerr(0.10),
        'eps_outer_elliptic_boundary':eps_ell,
        'eps_delta_over_2_mix_crossover':eps_mix,
    }
    _,_,_,_,re_c,th_c,leak_c,_=exact(eps_ell)
    thresholds['shift_error_at_elliptic_boundary']=re_c
    thresholds['mix_angle_deg_at_elliptic_boundary']=math.degrees(th_c)
    thresholds['inner_weight_at_elliptic_boundary']=leak_c

    md=HERE/'E3A1_perturbation_breakdown_analysis_v1.md'
    with open(md,'w',encoding='utf-8') as f:
        f.write('# E3-A1 階層ケプラー摂動近似の崩壊点テスト v3\n\n')
        f.write('## 問い\n\n')
        f.write('E3-A0-v2 の相互線形結合は正常モードへ厳密対角化できる。そこで新しい効果を作らず、')
        f.write('「無摂動の内部/外部二体系として扱う摂動近似」がどこまで有効かだけを調べる。\n\n')
        f.write('$$H=\\begin{pmatrix}\\tau_{in}&\\varepsilon\\\\\\varepsilon&\\tau_{out}\\end{pmatrix},\\quad')
        f.write('\\lambda_\\pm=\\frac{\\tau_{in}+\\tau_{out}}2\\pm\\sqrt{(\\Delta\\tau/2)^2+\\varepsilon^2}.$$\n\n')
        f.write(f'$n_{{in}}={N_IN}$, $n_{{out}}={N_OUT}$, $\\tau_{{in}}={TAU_IN:.12f}$, ')
        f.write(f'$\\tau_{{out}}={TAU_OUT:.12f}$, $\\Delta\\tau={DELTA:.12f}$.\n\n')
        f.write('## 結果\n\n')
        f.write(f'- $O(\\varepsilon^2)$ の外部法則シフトの相対誤差が 1%: $\\varepsilon={thresholds["eps_1pct_shift_error"]:.8f}$。\n')
        f.write(f'- 同 5%: $\\varepsilon={thresholds["eps_5pct_shift_error"]:.8f}$。\n')
        f.write(f'- 同 10%: $\\varepsilon={thresholds["eps_10pct_shift_error"]:.8f}$。\n')
        f.write(f'- 外部優勢正常モードが楕円型から放物境界 $\\lambda_+=2$ に達する厳密値: $\\varepsilon_c={eps_ell:.8f}$。\n')
        f.write(f'  この点でも混合角は {math.degrees(th_c):.3f} deg、内部成分重量は {leak_c:.4f}、二次近似誤差は {re_c:.4%}。\n')
        f.write(f'- $2\\varepsilon/\\Delta\\tau=1$ の自然な混合クロスオーバー: $\\varepsilon=\\Delta\\tau/2={eps_mix:.8f}$。')
        f.write('  ここで混合角は 22.5 deg、外部モード中の内部重量は 14.64%。\n\n')
        f.write('## 解釈\n\n')
        f.write('この候補では exact system 自体は崩れない。崩れるものが三段階に分かれる。\n\n')
        f.write('1. まず $\\Delta\\tau_{out}\\simeq\\varepsilon^2/\\Delta\\tau$ という摂動展開の精度が徐々に落ちる。\n')
        f.write('2. それより明確な境界として、$\\varepsilon_c$ で外部正常モードが楕円型を離れ、閉 Kepler 軌道という解釈が失われる。\n')
        f.write('3. さらに $\\varepsilon\\sim\\Delta\\tau$ では内部/外部の固有ベクトルが強く混ざり、(ab)|c という階層ラベル自体が弱くなる。\n\n')
        f.write('特に今回のパラメータでは、強混合より先に楕円型境界へ達する。したがって最初の物理的な「崩壊点」は ')
        f.write(f'$\\varepsilon_c\\simeq {eps_ell:.5f}$ と読むのが最も非恣意的である。\n')

    # Figures are generated separately from the frozen data by plot_E3A_fixed_point_v1.py.

    print('thresholds')
    for k,v in thresholds.items(): print(k, v)
    print('csv',csv_path)
    print('report',md)

if __name__=='__main__': main()
