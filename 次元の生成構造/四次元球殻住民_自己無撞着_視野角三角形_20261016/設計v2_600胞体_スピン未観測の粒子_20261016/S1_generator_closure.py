"""S1 生成子と閉塞：反対称・スペクトル・ω²=λ²+h²・正の振動数の状態で Σψ²=0・保存・ゲージと可換。"""
import json
import numpy as np
from scipy.linalg import expm
import m600_model as m

r = {}
rng = np.random.default_rng(20261016)
K0 = m.K_of(0.0)
KE = m.K_of(m.H_E)
r['h_e'] = m.H_E; r['theta_e_rad'] = m.THETA_E; r['R_e=cos^2theta'] = float(np.cos(m.THETA_E) ** 2)
r['antisym_err'] = float(max(np.abs(K0 + K0.T).max(), np.abs(KE + KE.T).max()))
# スペクトル（H5 上、ゼロは H5 の外の 29×4 と定数スピノル 4）
w0 = np.linalg.eigvals(K0)
vals, cnt = np.unique(np.round(np.abs(w0.imag), 9), return_counts=True)
r['spectrum_h0'] = {str(float(v)): int(c) for v, c in zip(vals, cnt)}
r['spectrum_integer_err'] = float(np.abs(np.abs(w0.imag) - np.round(np.abs(w0.imag))).max())
wE = np.linalg.eigvals(KE)
lam = np.sort(np.abs(w0.imag)); om = np.sort(np.abs(wE.imag))
inres = lam > 1e-9   # H5 の外を除いた比較
lamr = np.sort(lam[inres]); omr = np.sort(om[om > 1e-9])
r['omega2_eq_lambda2_plus_h2_err'] = float(np.abs(np.sqrt(lamr ** 2 + m.H_E ** 2) - omr[-len(lamr):]).max())
# 正の振動数の任意の状態で Σψ² = 0（閉塞は結果）
w, Vv = np.linalg.eig(KE)
pos = w.imag > 1e-9
errs = []
for _ in range(20):
    c = rng.normal(size=pos.sum()) + 1j * rng.normal(size=pos.sum())
    psi = Vv[:, pos] @ c; psi /= np.linalg.norm(psi)
    errs.append(abs((psi ** 2).sum()))
r['closure_positive_frequency_states_max'] = float(max(errs))
# 負の振動数と混ぜると閉じない（対照）
c = rng.normal(size=len(w)) + 1j * rng.normal(size=len(w))
psi = m.PIF @ (Vv @ c); psi /= np.linalg.norm(psi)
r['closure_mixed_state_control'] = float(abs((psi ** 2).sum()))
# 時間発展で Σψ² とノルムが保存（任意の状態）
U = expm(0.37 * KE)
r['evolution_orthogonal_err'] = float(np.abs(U.T @ U - np.eye(len(U))).max())
r['evolution_sum_sq_change'] = float(abs(((U @ psi) ** 2).sum() - (psi ** 2).sum()))
# ゲージ（左の平行移動）との可換
g10 = [v for v in m.V if abs(v[0] - m.PHI / 2) < 1e-9][0]
PgF = np.kron(np.eye(4), m.perm_left(g10))
r['commutator_K_with_left_translation'] = float(np.abs(KE @ PgF @ m.PIF - PgF @ KE @ m.PIF).max())
r['PASS'] = bool(r['antisym_err'] < 1e-12 and r['spectrum_integer_err'] < 1e-12 and r['omega2_eq_lambda2_plus_h2_err'] < 1e-12
                 and r['closure_positive_frequency_states_max'] < 1e-12 and r['evolution_orthogonal_err'] < 1e-12
                 and r['evolution_sum_sq_change'] < 1e-12 and r['commutator_K_with_left_translation'] < 1e-12)
json.dump(r, open('results/S1_generator_closure.json', 'w'), indent=1, ensure_ascii=False)
for k, v in r.items():
    print(k, v)
