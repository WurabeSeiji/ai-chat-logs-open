# -*- coding: utf-8 -*-
"""第十四思考実験 E0（2026-10-14）：荷電レプトン質量の零二乗和（Koide Q=2/3）と三相位相 delta の検算
Q=2/3 は |v_perp|^2-|v_par|^2=0（v=(sqrt m_e, sqrt m_mu, sqrt m_tau)、(1,1,1)方向を虚軸）と同値。
質量: PDG 2024（MeV）。要 mpmath。"""
from mpmath import mp, mpf, sqrt, acos, findroot, pi
mp.dps = 25
me, dme = mpf('0.51099895000'), mpf('0.00000000015')
mmu, dmu = mpf('105.6583755'), mpf('0.0000023')
mtau, dtau = mpf('1776.93'), mpf('0.09')
def Q(me, mmu, mt): return (me+mmu+mt)/(sqrt(me)+sqrt(mmu)+sqrt(mt))**2
def delta(me, mmu, mt):
    M = (sqrt(me)+sqrt(mmu)+sqrt(mt))/3
    return acos((sqrt(mt)/M-1)/sqrt(2))   # sqrt m_n = M(1+sqrt2 cos(delta+2pi n/3)), tau: n=0
for mt in (mtau-dtau, mtau, mtau+dtau):
    print('m_tau', mt, 'Q-2/3', Q(me, mmu, mt)-mpf(2)/3, 'delta-2/9', delta(me, mmu, mt)-mpf(2)/9)
def tau_from_Q(me, mmu): return findroot(lambda t: Q(me, mmu, t)-mpf(2)/3, 1777)
t0 = tau_from_Q(me, mmu); d0 = delta(me, mmu, t0)-mpf(2)/9
print('Q=2/3 -> m_tau =', t0, ' delta-2/9 =', d0, ' delta/pi =', (d0+mpf(2)/9)/pi)
s_mu = delta(me, mmu+dmu, tau_from_Q(me, mmu+dmu))-mpf(2)/9-d0
s_e = delta(me+dme, mmu, tau_from_Q(me+dme, mmu))-mpf(2)/9-d0
print('sensitivity from m_mu error', s_mu, ' from m_e error', s_e, ' residual/sigma', d0/s_mu)
