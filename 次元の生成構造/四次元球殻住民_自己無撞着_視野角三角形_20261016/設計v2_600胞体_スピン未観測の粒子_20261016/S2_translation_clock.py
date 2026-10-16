"""S2（v2：スピン未観測の混合）並進と時計：ゲージ目盛で読んだ速さと、予算の比 <p/ω>。モードの中の予算 v²+(1/γ)²=1。位相の回転＝左の平行移動。"""
import json
import numpy as np
import m600_model as m

r = {'runs': []}
HS = [0.0, m.H_E, 2.0]
K0S = [-4, -3, -2, -1, 0, 1, 2, 3, 4]
ts = np.linspace(0, 0.6, 13)
for h in HS:
    for k0 in K0S:
        comps, md = m.packet_mixed(h, k0, 0)
        wts = sum(np.abs(c) ** 2 for c in md['coefs'])      # スピンの混合の和（総和 1）
        p = -md['k']                         # 運動量（+u 方向を正）
        v_pred = float((wts * p / md['om']).sum())
        inv_gamma = float((wts * h / md['om']).sum())
        ph = []
        for t in ts:
            cs = [md['V'] @ (c * np.exp(-1j * md['om'] * t)) for c in md['coefs']]
            ph.append(m.fibre_phase_mixed(cs))
        v_read = float(np.polyfit(ts, np.unwrap(ph), 1)[0])
        psi = comps[0]
        r['runs'].append(dict(h=h, k0=k0, v_read=v_read, v_pred=v_pred, diff=v_read - v_pred,
                              inv_gamma_clock=inv_gamma, mean_omega=float((wts * md['om']).sum()),
                              closure=float(max(abs((c ** 2).sum()) for c in comps))))
# モードの中の予算：状態から計算した右の交換の割合 = (h/ω)²
Vv, om, kk = m.joint_modes(m.H_E)
sel = om > 1e-6
Mh = np.kron(m.M, m.H_E * m.PI)
em = np.real(-1j * np.einsum('ij,ik,kj->j', Vv[:, sel].conj(), Mh, Vv[:, sel]))
r['mode_budget_err'] = float(np.abs(em / om[sel] - (m.H_E / om[sel]) ** 2).max())
# 位相の回転＝左の平行移動：ゲージの一歩（置換）で固有モードが e^{ik·36°} 倍
PgF = np.kron(np.eye(4), m.perm_left(m.G10))
vs = Vv[:, sel]
r['translation_is_phase_err'] = float(np.abs(PgF @ vs - vs * np.exp(1j * kk[sel] * m.STEP)).max())
# 判定
d = np.array([x['diff'] for x in r['runs']])
mono = True
for h in HS:
    vr = [x['v_read'] for x in r['runs'] if x['h'] == h]
    dv = np.diff(vr)
    mono &= bool(np.all(dv > 0) or np.all(dv < 0))   # 設計書：k0 に対して単調（向きは問わない）
vmax = {h: max(abs(x['v_read']) for x in r['runs'] if x['h'] == h) for h in HS}
r['max_abs_diff'] = float(np.abs(d).max())
r['monotone_in_k0'] = mono
r['max_speed_by_h'] = {str(h): v for h, v in vmax.items()}
r['PASS'] = bool(r['mode_budget_err'] < 1e-12 and r['translation_is_phase_err'] < 1e-12 and r['max_abs_diff'] <= 0.1
                 and mono and vmax[0.0] > vmax[m.H_E] > vmax[2.0])
json.dump(r, open('results/S2_translation_clock.json', 'w'), indent=1, ensure_ascii=False)
for x in r['runs']:
    print(f"h={x['h']:.3f} k0={x['k0']:+d} v_read={x['v_read']:+.3f} v_pred={x['v_pred']:+.3f} diff={x['diff']:+.3f} 1/gamma={x['inv_gamma_clock']:.3f}")
for k in ['mode_budget_err', 'translation_is_phase_err', 'max_abs_diff', 'monotone_in_k0', 'max_speed_by_h', 'PASS']:
    print(k, r[k])
