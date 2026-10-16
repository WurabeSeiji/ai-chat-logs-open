# -*- coding: utf-8 -*-
"""段 3c の事後診断（判定ではない。走行 2 の後に追加）：
走行 2 の場合 B では、交差項の周波数が 248 ステップのフーリエの分解能（1/248）より小さく、周波数が η で変わるかを判別できなかった。
次数ごとの交差項 X_n(k) は単一の周波数の指数関数になるはずなので、1 ステップあたりの位相の進み arg(X_n(k+1)/X_n(k)) を
次数ごと・η ごとに直接求める。"""
import json, numpy as np
exec(open('stage3c_boost_exchange_beat.py').read().split("out = dict(cases=[])")[0])
def Bn(Cf, Cg, n): return np.trace(Cf[n].T @ np.linalg.inv(J(n)) @ Cg[n] @ J(n))
tp, to = 2*np.pi*23/124, 2*np.pi*29/156; bp, bo = qexp(tp, 'k'), qexp(to, 'k')
rows = []
for u in 'ijk':
    for eta in (0.0, 0.1, 0.2, 0.3):
        f = translate_left(f0, qexp(1j*eta, u)); g = g0
        f1 = rotate_right(f, bp); g1 = rotate_right(g, bo)
        ph = []; amp = []
        for n in range(1, N+1):
            x0, x1 = Bn(f, g, n), Bn(f1, g1, n)
            ph.append(float(np.angle(x1/x0)) if abs(x0) > 1e-14 else None); amp.append(float(abs(x0)))
        rows.append(dict(u=u, eta=eta, phase_per_step_by_degree=ph, amp_by_degree=amp))
pred = [float(np.angle(np.exp(1j*n*(tp - to)))) for n in range(1, N+1)]
spread = max(max(abs(r['phase_per_step_by_degree'][i] - rows[0]['phase_per_step_by_degree'][i]) for r in rows if r['phase_per_step_by_degree'][i] is not None) for i in range(N))
out = dict(rows=rows, predicted_n_times_dtheta=pred, max_spread_over_eta_and_axis=spread)
json.dump(out, open('stage3c_posthoc_per_degree_beat.json', 'w'), indent=1)
print('phase per step by degree (eta=0, u=i):', np.round(rows[0]['phase_per_step_by_degree'][:6], 8))
print('predicted n*(θp−θo)            :', np.round(pred[:6], 8))
print('max spread of the per-degree phase over all η and boost axes:', spread)
print('amplitude by degree, u=i, eta 0 vs 0.3:', np.round(rows[0]['amp_by_degree'][:5], 4), np.round(rows[3]['amp_by_degree'][:5], 4))
