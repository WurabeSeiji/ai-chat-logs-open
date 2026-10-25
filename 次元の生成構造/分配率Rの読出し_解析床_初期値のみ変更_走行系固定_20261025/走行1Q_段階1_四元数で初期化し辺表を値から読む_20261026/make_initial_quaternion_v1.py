#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""走行1Q 段階1（2026-10-26）：走行1 の床を四元数で初期化し、辺表 (ea, eb) と波 z を四元数の値だけから読み直す。

符号化（走行1 の床 npz の Z0, ea, eb を入力にする）
  関係 e＝(v, w) に対し   q_e = e^(i·πv/N) · r_e · e^(j·πw/N)
    v＝始点（左から i で回す）、w＝終点（右から j で回す）。一周 2π は 2N 目盛り、頂点は 0..N−1 の N 目盛り（半周）、
    r_e の符号は残りの N 目盛り（e^(iπ)＝−1）。
  r_e は走行1 の床から  Z0_e = e^(iφ0) · r_e · i^((v+w) mod 2)  で取る実数（符号つき）。φ0 は全体位相 1 個（スカラー）。
    この書き方が成り立つこと（虚部残差が丸めの範囲）を検査し、成り立たなければ停止する。
  成分：q = a + b i + c j + d k、 a=r cosα cosβ, b=r sinα cosβ, c=r cosα sinβ, d=r sinα sinβ（ad＝bc）。

復号（配列の添字 e も e＝v·N＋w の規則も使わない。各 q_e の 4 成分と N、φ0 だけを使う）
  α：(a,b) と (c,d) のうちノルムの大きい組の atan2 を π で割った余り → v＝round(αN/π) mod N
  β：(a,c) と (b,d) のうちノルムの大きい組の atan2 を π で割った余り → w＝round(βN/π) mod N
  r＝Re( e^(−iπv/N) · q · e^(−jπw/N) )、 z＝e^(iφ0) · r · i^((v+w) mod 2)
  ea＝v、eb＝w として npz に書く。走行系（run_sweep_from_parents_v1.py、走行1 とバイト同一）は Z0, ea, eb を読み、A は走行系が作る。

使い方: python make_initial_quaternion_v1.py RUN1_FLOOR_DIR OUT_DIR N1 N2 ..."""
import json
import math
import os
import sys

import numpy as np


def qmul(p, q):
    """四元数の積（成分 [a, b, c, d] = a + b i + c j + d k、ij＝k）。p, q は (..., 4)。"""
    a1, b1, c1, d1 = np.moveaxis(p, -1, 0)
    a2, b2, c2, d2 = np.moveaxis(q, -1, 0)
    return np.stack([a1 * a2 - b1 * b2 - c1 * c2 - d1 * d2,
                     a1 * b2 + b1 * a2 + c1 * d2 - d1 * c2,
                     a1 * c2 - b1 * d2 + c1 * a2 + d1 * b2,
                     a1 * d2 + b1 * c2 - c1 * b2 + d1 * a2], axis=-1)


def exp_i(t):
    z = np.zeros(np.shape(t) + (4,)); z[..., 0] = np.cos(t); z[..., 1] = np.sin(t); return z


def exp_j(t):
    z = np.zeros(np.shape(t) + (4,)); z[..., 0] = np.cos(t); z[..., 2] = np.sin(t); return z


def encode(Z0, ea, eb, N):
    p = (ea + eb) % 2
    m = int(np.argmax(np.abs(Z0)))
    phi0 = float(np.angle(Z0[m]) - 0.5 * math.pi * p[m])
    w = Z0 * np.exp(-1j * phi0) * (1j) ** (-p)
    r = w.real
    imag_resid = float(np.max(np.abs(w.imag)) / np.max(np.abs(Z0)))
    rq = np.zeros((len(r), 4)); rq[:, 0] = r
    Q = qmul(qmul(exp_i(math.pi * ea / N), rq), exp_j(math.pi * eb / N))
    return Q, phi0, r, imag_resid


def _angle_mod_pi(x1, y1, x2, y2):
    n1 = np.hypot(x1, y1); n2 = np.hypot(x2, y2)
    ang = np.where(n1 >= n2, np.arctan2(y1, x1), np.arctan2(y2, x2))
    return np.mod(ang, math.pi)


def decode(Q, N, phi0):
    a, b, c, d = Q[:, 0], Q[:, 1], Q[:, 2], Q[:, 3]
    alpha = _angle_mod_pi(a, b, c, d)
    beta = _angle_mod_pi(a, c, b, d)
    v = np.mod(np.rint(alpha * N / math.pi).astype(np.int64), N)
    w = np.mod(np.rint(beta * N / math.pi).astype(np.int64), N)
    core = qmul(qmul(exp_i(-math.pi * v / N), Q), exp_j(-math.pi * w / N))
    r = core[:, 0]
    z = (np.exp(1j * phi0) * r * (1j) ** ((v + w) % 2)).astype(np.complex128)
    rest = float(np.max(np.abs(core[:, 1:])))
    return v, w, r, z, alpha, beta, rest


if __name__ == '__main__':
    src, outdir = sys.argv[1], sys.argv[2]
    Ns = [int(x) for x in sys.argv[3:]]
    os.makedirs(outdir, exist_ok=True)
    manifest = {}
    print(' N     M | 符号化:虚部残差  min|r| | 復号: ea一致 eb一致 自己対数 残り成分 | max|z復号−Z0|  max|ad−bc|')
    for N in Ns:
        name = f'floor_90deg_all_N{N:05d}_analytic.npz'
        d = dict(np.load(os.path.join(src, name)))
        Z0 = np.asarray(d['Z0'], dtype=np.complex128); ea0 = d['ea'].astype(np.int64); eb0 = d['eb'].astype(np.int64)
        Q, phi0, r0, imag_resid = encode(Z0, ea0, eb0, N)
        assert imag_resid < 1e-12, f'N={N}: Z0 = e^(iφ0) r i^((v+w) mod 2) が成り立たない（虚部残差 {imag_resid:.2e}）。停止。'
        v, w, r, z, alpha, beta, rest = decode(Q, N, phi0)
        ea_ok = bool(np.array_equal(v, ea0)); eb_ok = bool(np.array_equal(w, eb0))
        n_self = int(np.sum(v == w))
        dz = float(np.max(np.abs(z - Z0)))
        adbc = float(np.max(np.abs(Q[:, 0] * Q[:, 3] - Q[:, 1] * Q[:, 2])))
        out = dict(d)
        out.update(v=z, Z0=z, theta=np.angle(z), ea=v, eb=w, Q=Q, phi0=np.array(phi0),
                   rule=np.array('RUN1Q stage1: run1 analytic floor encoded as q=e^(i pi v/N) r e^(j pi w/N); ea, eb, Z0 decoded from q values only'))
        np.savez_compressed(os.path.join(outdir, name), **out)
        chk = dict(N=N, M=int(len(Z0)), phi0=phi0, encode_imag_residual=imag_resid, min_abs_r=float(np.min(np.abs(r0))),
                   decoded_ea_equal=ea_ok, decoded_eb_equal=eb_ok, n_self_pairs_alpha_eq_beta=n_self,
                   decode_rest_components_max=rest, max_abs_z_decoded_minus_Z0=dz,
                   z_decoded_bit_identical=bool(np.array_equal(z, Z0)), max_abs_ad_minus_bc=adbc)
        manifest[name] = chk
        print(f"{N:>2} {len(Z0):>5} | {imag_resid:.1e}  {chk['min_abs_r']:.3e} | {str(ea_ok):>5} {str(eb_ok):>5} {n_self:>4} {rest:.1e} | "
              f"{dz:.1e} {str(chk['z_decoded_bit_identical']):>5}  {adbc:.1e}", flush=True)
    json.dump(manifest, open(os.path.join(outdir, 'manifest_floor_quaternion.json'), 'w'), indent=2, ensure_ascii=False)
    print('QUATERNION FLOORS DONE')
