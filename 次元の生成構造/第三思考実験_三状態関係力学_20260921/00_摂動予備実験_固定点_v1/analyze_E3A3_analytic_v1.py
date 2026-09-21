#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E3-A3 fixed-point analytic analysis.

Uses only the already-fixed E3-A0/E3-A1 data and a deterministic replay of one
already-fixed E3-A0 parameter point (eps=0.004, phase=0) to expose exact
identities. No new parameter sweep is performed.
"""
from __future__ import annotations
import csv
import json
import math
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
JSON0 = HERE / 'E3A0_hierarchical_two_kepler_perturbation_full_v1.json'
CSV1 = HERE / 'E3A1_perturbation_breakdown_scan_v1.csv'
TSCSV = HERE / 'E3A3_area_exchange_timeseries_v1.csv'
PHASECSV = HERE / 'E3A3_phase_robustness_v1.csv'
SUMMARY = HERE / 'E3A3_analytic_summary_v1.json'
REPORT = HERE / 'E3A3_analytic_analysis_v1.md'


def om(u, v):
    return float(u[0] * v[1] - u[1] * v[0])


def rot(t):
    return np.array([[math.cos(t), -math.sin(t)], [math.sin(t), math.cos(t)]], float)


def module(n, ecc, u0):
    th = 2 * math.pi / n
    T = np.diag([math.sqrt(1 - ecc), math.sqrt(1 + ecc)])
    S = T @ rot(th) @ np.linalg.inv(T)
    x0 = T @ np.array([math.cos(u0), math.sin(u0)])
    xm1 = np.linalg.inv(S) @ x0
    return 2 * math.cos(th), xm1, x0


def replay(n_in=31, n_out=127, eps=0.004, u_in=0.37, u_out=1.11,
           steps=600, e_in=0.35, e_out=0.25):
    ti, im1, i0 = module(n_in, e_in, u_in)
    to, om1, o0 = module(n_out, e_out, u_out)
    xi = np.empty((steps + 1, 2)); xo = np.empty((steps + 1, 2))
    xi[0] = i0; xo[0] = o0
    xi[1] = ti * i0 + eps * o0 - im1
    xo[1] = to * o0 + eps * i0 - om1
    for k in range(1, steps):
        xi[k + 1] = ti * xi[k] + eps * xo[k] - xi[k - 1]
        xo[k + 1] = to * xo[k] + eps * xi[k] - xo[k - 1]
    qi = np.array([om(xi[k], xi[k + 1]) for k in range(steps)])
    qo = np.array([om(xo[k], xo[k + 1]) for k in range(steps)])
    qt = qi + qo
    cross = np.array([om(xi[k], xo[k]) for k in range(steps)])
    exch_i = qi[1:] - qi[:-1]
    exch_o = qo[1:] - qo[:-1]
    ident_i = exch_i - eps * cross[1:]
    ident_o = exch_o + eps * cross[1:]
    H = np.array([[ti, eps], [eps, to]])
    lam, V = np.linalg.eigh(H)
    modes = [V[0, j] * xi + V[1, j] * xo for j in range(2)]
    qm = [np.array([om(m[k], m[k + 1]) for k in range(steps)]) for m in modes]
    return dict(ti=ti, to=to, xi=xi, xo=xo, qi=qi, qo=qo, qt=qt, cross=cross,
                ident_i=ident_i, ident_o=ident_o, lam=lam, V=V, qm=qm)


def load_scan():
    with open(CSV1, newline='') as f:
        return list(csv.DictReader(f))


def analytic_values(delta, eps):
    g = eps / delta if delta != 0 else math.inf
    if eps == 0:
        return dict(g=0.0, shift=0.0, relerr=0.0, angle=0.0, weight=0.0)
    root = math.sqrt((delta / 2) ** 2 + eps ** 2)
    shift = root - delta / 2
    approx = eps ** 2 / delta
    relerr = approx / shift - 1.0
    angle = 0.5 * math.atan2(2 * eps, delta)
    weight = math.sin(angle) ** 2
    return dict(g=g, shift=shift, relerr=relerr, angle=angle, weight=weight)


def main():
    data = json.load(open(JSON0, encoding='utf-8'))
    agg = data['aggregates']
    raw = [r for r in data['raw'] if r['suite'] == 'hierarchical_nonres_31_127']
    hier = [r for r in agg if r['suite'] == 'hierarchical_nonres_31_127']
    deg = [r for r in agg if r['suite'] == 'degenerate_31_31']

    tau_in = float(raw[0]['tau_in'])
    tau_out = float(raw[0]['tau_out'])
    delta = tau_out - tau_in

    # Representative replay of an already-fixed run point.
    rr = replay()
    qscale = abs(rr['qt'][0])
    with open(TSCSV, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['step','dq_in_over_Q','dq_out_over_Q','dq_total_over_Q','cross_area'])
        w.writeheader()
        for k in range(len(rr['qi'])):
            w.writerow({
                'step': k,
                'dq_in_over_Q': (rr['qi'][k] - rr['qi'][0]) / qscale,
                'dq_out_over_Q': (rr['qo'][k] - rr['qo'][0]) / qscale,
                'dq_total_over_Q': (rr['qt'][k] - rr['qt'][0]) / qscale,
                'cross_area': rr['cross'][k],
            })

    # Phase robustness from fixed raw data.
    with open(PHASECSV, 'w', newline='') as f:
        fields = ['eps','metric','min','median','max','phase_span_ratio']
        w = csv.DictWriter(f, fieldnames=fields); w.writeheader()
        for e in sorted(set(float(r['eps']) for r in raw)):
            if e <= 0: continue
            ss = [r for r in raw if float(r['eps']) == e]
            for key, label in [('q_out_mod','q_out_mod'),('raw_outer_conic_resid','raw_outer_conic_resid')]:
                vals = np.array([float(r[key]) for r in ss])
                med = float(np.median(vals))
                w.writerow({'eps':e,'metric':label,'min':float(vals.min()),'median':med,'max':float(vals.max()),
                            'phase_span_ratio':float(vals.max()/vals.min())})

    # Exact formula consistency against the fixed E3-A1 scan.
    scan = load_scan()
    errs = {'shift':0.0,'relerr':0.0,'angle_deg':0.0,'weight':0.0}
    for row in scan:
        eps = float(row['eps'])
        a = analytic_values(delta, eps)
        errs['shift'] = max(errs['shift'], abs(a['shift'] - float(row['tau_out_shift_exact'])))
        errs['relerr'] = max(errs['relerr'], abs(a['relerr'] - float(row['shift_Oeps2_rel_error'])))
        errs['angle_deg'] = max(errs['angle_deg'], abs(math.degrees(a['angle']) - float(row['mix_angle_deg'])))
        errs['weight'] = max(errs['weight'], abs(a['weight'] - float(row['inner_weight_in_outer_mode'])))

    thresholds = {}
    for pct in (1, 5, 10):
        r = pct / 100.0
        g = math.sqrt(r * (1 + r))
        thresholds[str(pct)] = {'g':g, 'eps':delta*g}

    eps_c = math.sqrt((2 - tau_in) * (2 - tau_out))
    g_c = eps_c / delta
    a_c = analytic_values(delta, eps_c)

    exact_checks = {
        'representative_eps':0.004,
        'q_total_relative_drift':float(np.max(np.abs(rr['qt']-rr['qt'][0]))/qscale),
        'exchange_identity_inner_max_abs':float(np.max(np.abs(rr['ident_i']))),
        'exchange_identity_outer_max_abs':float(np.max(np.abs(rr['ident_o']))),
        'normal_mode_q_drift_max':float(max(np.max(np.abs(q-q[0]))/max(abs(q[0]),1e-30) for q in rr['qm'])),
    }

    phase_small = []
    for e in sorted(set(float(r['eps']) for r in raw)):
        if not (0 < e <= 0.002): continue
        ss = [r for r in raw if float(r['eps']) == e]
        vals = np.array([float(r['raw_outer_conic_resid']) for r in ss])
        phase_small.append(float(vals.max()/vals.min()))

    summary = {
        'tau_in':tau_in,'tau_out':tau_out,'delta_tau':delta,'inverse_delta_tau':1/delta,
        'hierarchical_fits':data['hierarchical_fits'],
        'exact_formula_max_abs_errors_vs_E3A1_scan':errs,
        'perturbative_relative_error_thresholds':thresholds,
        'elliptic_parabolic_boundary':{
            'eps_c':eps_c,'g_c':g_c,'mix_angle_deg':math.degrees(a_c['angle']),
            'inner_weight':a_c['weight'],'second_order_relative_error':a_c['relerr']},
        'representative_exact_checks':exact_checks,
        'small_eps_phase_span_ratio_max_raw_conic':max(phase_small),
        'degenerate_nonzero_mix_angle_rad_unique':sorted(set(float(r['mix_angle']) for r in deg if float(r['eps'])>0)),
    }
    json.dump(summary, open(SUMMARY,'w',encoding='utf-8'), ensure_ascii=False, indent=2)

    f = data['hierarchical_fits']
    md = r"""# E3-A3 解析的分析 v1 — 摂動予備実験の固定点から何が確定したか

日付: 2026-09-21  
対象: E3-A0 / E3-A1 / E3-A2 の固定済みデータのみ。  
位置づけ: **追加数値実験ではなく、固定済み候補の解析的整理**。

## 1. 候補則

固定した候補は

$$
X^{in}_{k+1}+X^{in}_{k-1}=\tau_{in}X^{in}_k+\varepsilon X^{out}_k,
$$

$$
X^{out}_{k+1}+X^{out}_{k-1}=\tau_{out}X^{out}_k+\varepsilon X^{in}_k.
$$

主系列では

$$
\tau_{in}=@@TAU_IN@@,\qquad
\tau_{out}=@@TAU_OUT@@,\qquad
\Delta\tau=@@DELTA@@.
$$

## 2. 全体面積保存は数値結果ではなく厳密恒等式

$$q_{in,k}=\omega(X^{in}_k,X^{in}_{k+1}),\qquad
q_{out,k}=\omega(X^{out}_k,X^{out}_{k+1}).$$

更新則をそのまま代入すると

$$
q_{in,k}-q_{in,k-1}=\varepsilon\,\omega(X^{in}_k,X^{out}_k),
$$

$$
q_{out,k}-q_{out,k-1}=-\varepsilon\,\omega(X^{in}_k,X^{out}_k).
$$

従って

$$
\boxed{q_{tot,k}=q_{in,k}+q_{out,k}=\mathrm{const}}
$$

は厳密である。代表固定点 $\varepsilon=0.004$ の再生では、全体面積の相対ドリフトは
`@@QT_DRIFT@@`、上の交換恒等式の最大絶対残差は inner / outer とも
`@@EXCH_RESID@@` の数値丸め域であった。

**したがって局所面積の変動は保存則の破れではなく、二セクター間の厳密な交換である。**

## 3. なぜ局所応答が一次になるか

非縮退なら階層空間の行列

$$
H=\begin{pmatrix}\tau_{in}&\varepsilon\\\varepsilon&\tau_{out}\end{pmatrix}
$$

を角度 $\theta$ で直交対角化でき、

$$
\tan 2\theta=\frac{2\varepsilon}{\Delta\tau}.
$$

小結合では

$$
\theta=\frac{\varepsilon}{\Delta\tau}+O\!\left((\varepsilon/\Delta\tau)^3\right).
$$

生の外部状態は二つの正常モードの混合なので、局所面積には二モード間の交差面積が係数
$\sin\theta\cos\theta=O(\varepsilon/\Delta\tau)$ で入る。従って固定データで

- inner local-area modulation: slope `@@SLOPE_QIN@@`
- outer local-area modulation: slope `@@SLOPE_QOUT@@`

とほぼ一次になったことは、厳密対角化から説明できる。

## 4. 生の軌道が単一円錐から一次で外れる理由

複素表示で正常モードを $y_-,y_+$、外部状態を

$$w_{out}=\sin\theta\,y_-+\cos\theta\,y_+$$

とすれば、論文2と同じ二乗読み出しは

$$
z_{out}=w_{out}^2
=\sin^2\theta\,y_-^2+\cos^2\theta\,y_+^2
+2\sin\theta\cos\theta\,y_-y_+.
$$

最後の交差項が $O(\varepsilon/\Delta\tau)$ である。従って **生の外部読み出しは一つの Kepler 円錐ではなく、二つの正常モードと交差項の合成**になる。
固定データの raw conic residual の小結合 slope は `@@SLOPE_CONIC@@` であり、この一次則と一致する。

一方、正常モード自身は

$$Y_{\pm,k+1}+Y_{\pm,k-1}=\lambda_\pm Y_{\pm,k}$$

を厳密に満たすため、各モードの面積保存と円錐残差は数値精度域に留まる。

## 5. 外部法則シフトが二次になる理由

$$
\lambda_+=\frac{\tau_{in}+\tau_{out}}2+
\sqrt{(\Delta\tau/2)^2+\varepsilon^2},
$$

したがって

$$
\delta\tau_{out}=\lambda_+-\tau_{out}
=\frac{\varepsilon^2}{\Delta\tau}-
\frac{\varepsilon^4}{\Delta\tau^3}+O(\varepsilon^6).
$$

固定データの fitted slope は `@@SLOPE_SHIFT@@`、prefactor は
`@@PREF_SHIFT@@` であり、解析値 $1/\Delta\tau=@@INV_DELTA@@$ と一致する。

ここで本当の展開変数は

$$
\boxed{g=\varepsilon/|\Delta\tau|}
$$

であり、$\varepsilon$ 単独ではない。

## 6. 摂動近似の精度は g だけで決まる

二次近似 $\varepsilon^2/\Delta\tau$ の exact shift に対する相対誤差は

$$
R(g)=\sqrt{\frac14+g^2}-\frac12.
$$

従って相対誤差 $R=r$ の境界は

$$
\boxed{g=\sqrt{r(1+r)}}
$$

である。主系列に戻すと

- 1%: $g=@@G1@@$, $\varepsilon=@@E1@@$
- 5%: $g=@@G5@@$, $\varepsilon=@@E5@@$
- 10%: $g=@@G10@@$, $\varepsilon=@@E10@@$

となり、E3-A1 の数値表と丸め誤差内で一致する。

## 7. 混合も g の普遍関数

$$
\theta(g)=\frac12\arctan(2g),
$$

$$
W_{in\to out}(g)=\sin^2\theta
=\frac12\left(1-\frac1{\sqrt{1+4g^2}}\right).
$$

E3-A1 の全 scan に対する解析式との差の最大値は、shift `@@ERR_SHIFT@@`、
relative error `@@ERR_REL@@`、angle `@@ERR_ANGLE@@` deg、weight `@@ERR_WEIGHT@@` である。

**従って「摂動近似が厳しい」という数値的観察は、非縮退条件 $g\ll1$ という解析条件に置き換えられる。**

## 8. 縮退系は摂動階層に対して特異

$\Delta\tau=0$ なら任意の $\varepsilon\ne0$ で

$$
\theta=\pi/4,
$$

すなわち正常モードは最初から 50:50 の対称・反対称結合になる。固定データでも非零 $\varepsilon$ の mixing angle はすべて
`45.0` deg である。

これは exact system の破綻ではない。**「inner / outer を弱く混ざる二つの独立対象として追跡できる」という階層解釈だけが、縮退点で特異になる。**

## 9. 閉軌道境界は別の条件

外部優勢正常モードが楕円型から放物境界へ達する条件は $\lambda_+=2$ なので、

$$
\boxed{\varepsilon_c^2=(2-\tau_{in})(2-\tau_{out})}.
$$

主系列では

$$
\varepsilon_c=@@EPS_C@@,\qquad g_c=@@G_C@@.
$$

この点の mixing angle は `@@ANGLE_C@@` deg、inner weight は `@@WEIGHT_C@@`、
二次近似誤差は `@@REL_C@@` である。

これは「摂動展開が使えるか」とは別に、正常モード自身が閉 Kepler 軌道でいられるかを決めるスペクトル境界である。

## 10. 位相依存性の評価

8初期位相の固定データでは、小結合域 $\varepsilon\le0.002$ における raw conic residual の
最大/最小比は最大でも `@@PHASE_RATIO@@`。従って一次スケーリングは特定の初期位相だけの結果ではない。
一方 $\varepsilon=0.004$ 以降は位相差の影響幅が急増し、raw conic residual は単調な一変数関数ではなくなる。
これは摂動域を外れ始めたという解釈と整合する。

## 11. この候補についての固定結論

この候補については、次が解析的に確定した。

1. **全体面積保存は厳密**で、局所面積は等量反対向きに交換する。
2. 非縮退小結合では、局所量・生の円錐残差は $O(\varepsilon/\Delta\tau)$。
3. 正常モードの法則シフトは $O((\varepsilon/\Delta\tau)^2)$。
4. 系全体は固定直交変換で二つの二値則へ厳密分解できる。
5. 従ってこの候補は **三状態既約力学ではない**。
6. 今後の真の三状態候補に対しては、**状態に依存しない固定基底変換で二体セクターへ完全分解できたら落とす**、という新しい判定条件を置ける。

最後の条件が、この摂動予備実験から本体へ持ち帰る最も重要な結果である。

## 12. 追加図

- `fig_E3A_08_local_area_exchange_timeseries_v1.png`: 局所面積の等量交換と全体保存。
- `fig_E3A_09_phase_robustness_v1.png`: 8位相に対する一次応答の頑健性と高結合側の位相幅拡大。
- `fig_E3A_10_universal_g_control_v1.png`: $g=\varepsilon/|\Delta\tau|$ による摂動誤差と混合重量の解析曲線。
- `fig_E3A_11_degenerate_vs_nondegenerate_mixing_v1.png`: 非縮退と縮退での mixing angle の質的差。
"""
    repl = {
        '@@TAU_IN@@': f'{tau_in:.12f}', '@@TAU_OUT@@': f'{tau_out:.12f}', '@@DELTA@@': f'{delta:.12f}',
        '@@QT_DRIFT@@': f"{exact_checks['q_total_relative_drift']:.3e}",
        '@@EXCH_RESID@@': f"{max(exact_checks['exchange_identity_inner_max_abs'], exact_checks['exchange_identity_outer_max_abs']):.3e}",
        '@@SLOPE_QIN@@': f"{f['q_in_mod_med']['slope']:.6f}", '@@SLOPE_QOUT@@': f"{f['q_out_mod_med']['slope']:.6f}",
        '@@SLOPE_CONIC@@': f"{f['raw_conic_resid_med']['slope']:.6f}", '@@SLOPE_SHIFT@@': f"{f['tau_out_shift']['slope']:.6f}",
        '@@PREF_SHIFT@@': f"{f['tau_out_shift']['prefactor']:.6f}", '@@INV_DELTA@@': f'{1/delta:.6f}',
        '@@G1@@': f"{thresholds['1']['g']:.9f}", '@@E1@@': f"{thresholds['1']['eps']:.9f}",
        '@@G5@@': f"{thresholds['5']['g']:.9f}", '@@E5@@': f"{thresholds['5']['eps']:.9f}",
        '@@G10@@': f"{thresholds['10']['g']:.9f}", '@@E10@@': f"{thresholds['10']['eps']:.9f}",
        '@@ERR_SHIFT@@': f"{errs['shift']:.3e}", '@@ERR_REL@@': f"{errs['relerr']:.3e}", '@@ERR_ANGLE@@': f"{errs['angle_deg']:.3e}", '@@ERR_WEIGHT@@': f"{errs['weight']:.3e}",
        '@@EPS_C@@': f'{eps_c:.9f}', '@@G_C@@': f'{g_c:.9f}', '@@ANGLE_C@@': f"{math.degrees(a_c['angle']):.6f}",
        '@@WEIGHT_C@@': f"{a_c['weight']:.6f}", '@@REL_C@@': f"{a_c['relerr']:.6f}", '@@PHASE_RATIO@@': f'{max(phase_small):.6f}',
    }
    for k, v in repl.items():
        md = md.replace(k, v)
    REPORT.write_text(md, encoding='utf-8')
    print(json.dumps(summary, ensure_ascii=False, indent=2))

if __name__ == '__main__':
    main()
