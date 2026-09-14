#!/usr/bin/env python3
"""SENS-2: 終端近傍 step 窓の感度検査（設計書 v1.2 §14）。

問い: 「最終 step のみ」の主推定に対し、窓（lock-40, lock-20, lock の3点）を混ぜたとき、
残余緩和方向の混入で次元推定が上振れするか。人工的な次元増加でないかを判定する。

事前宣言: 代表 N={6,10,16,20,30,40}、構成 (A,S1,K=N-1)、
比較量 d_cons(final-only, n=5N) vs d_cons(window-mixed, n=15N)。
窓データは sens_generate_window_states.py（stage1 と bit 一致検証済み）を用いる。
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import importlib.util

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('rs', HERE/'run_sens.py')
rs = importlib.util.module_from_spec(spec); spec.loader.exec_module(rs)

def main():
    rows = []
    for N in [6, 10, 16, 20, 30, 40]:
        Kmax = N - 1
        ea, eb = rs.edges(N)
        w = np.load(HERE/'window_states'/f'N{N}_window.npz')
        Zw = w['Z_window']  # (5, 3, M)
        def stars(states):
            Ps = []
            for z in states:
                for v in range(N):
                    zs = z[(ea == v) | (eb == v)]
                    H = np.sum(np.abs(zs)**2)
                    Ps.append(np.array([np.sum(zs**k)/(H**(k/2)) for k in range(1, Kmax+1)]))
            return np.array(Ps)
        P_final = stars([Zw[s, -1] for s in range(5)])
        P_mixed = stars([Zw[s, t] for s in range(5) for t in range(3)])
        (d_f, *_), _, n_f = rs.estimate(P_final, Kmax, 'A', 'S1', 1e-12)
        (d_m, *_), _, n_m = rs.estimate(P_mixed, Kmax, 'A', 'S1', 1e-12)
        rows.append({'N': N, 'n_final': n_f, 'd_final': d_f,
                     'n_window': n_m, 'd_window': d_m,
                     'window_minus_final': d_m - d_f})
        print(f'N={N:3d} d_final={d_f:.2f} (n={n_f})  d_window={d_m:.2f} (n={n_m})  '
              f'diff={d_m-d_f:+.2f}', flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(HERE/'sens2_window_results.csv', index=False)
    summary = {'max_abs_window_effect': float(df.window_minus_final.abs().max()),
               'mean_window_effect': float(df.window_minus_final.mean()),
               'conclusion_gt4_unchanged': bool((df.d_window > 4).equals(df.d_final > 4))}
    (HERE/'sens2_summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    print(json.dumps(summary, indent=2))

if __name__ == '__main__':
    main()
