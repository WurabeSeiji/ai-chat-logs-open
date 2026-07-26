#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""第8論文 本実験 共有観測モジュール（float64）。解釈しない（機械的観測量の算出のみ）。

観測モデル（全て機械的定義・prior整合・報告書§に明記）:
- 瞬時生成子 K(θ(Z)) の固有平面（JG 固有分解, σ降順）w_1..w_N を「方向」とする（観測のみ・状態へ非帰還）。
- q_j = sqrt( |<w_j, Z>|^2 / |Z|^2 )  … 平面 j への Z 占有振幅（j=1..8, N<8 は 0 詰め）。
- 親平面 P1（親 v の実2次元, 固定基準）: f_outside = 1 - E_P1/|Z|^2, a_outside=sqrt(max(0,f_outside))。
- 方向連続性: |<w_j(t), w_j(t-1)>|（符号・順序照合後）。回転角 = arccos。
- 流束: 平面占有エネルギー E_j の step 差の機械的 proxy（定義は関数 fluxes に明記）。
- 適応ノイズ床 η_noise: closure/norm/最小非零/精度eps/量子化補正/再射影補正/再現run差 の最大（§10.1）。
第5論文エンジン（LowRankSystem, sigma_spectrum, make_parent）を read-only import で再利用（不変更）。
"""
import sys
from pathlib import Path

import numpy as np

CODE = Path(__file__).resolve().parent
PAPER8 = CODE.parent.parent
REPO = PAPER8.parent.parent
ENGINE = REPO / "時間軸Q軸とフェルミオンの生成構造/検証_対照実験/第5論文原本_自発的分裂予備実験_v1"
sys.path.insert(0, str(ENGINE))
from run_n_scaling_lowrank_v1 import LowRankSystem, make_parent, GAMMA  # noqa: E402

TAU_FIX = 1e-8            # rank_fixed_threshold 用固定閾値（q_j/q1）
C_SPROUT = 4.0           # 萌芽係数（q_j > C_SPROUT * eta_noise）
C_ESTABLISH = 1e-3       # 成立: q_j/q1 適応閾値（Stage A で校正・以後固定）


# ---------------- 親平面（固定基準 P1） ----------------

def parent_plane(v):
    """親 v の実2次元正規直交基底 P1 (M×2)。onset_probe と同一構成。"""
    p = v.real / np.linalg.norm(v.real)
    q = v.imag - (v.imag @ p) * p
    q = q / np.linalg.norm(q)
    return np.column_stack([p, q])


def occ_real2(B, Z):
    """実基底 B (M×2) への複素 Z の占有エネルギー |B^T Zr|^2+|B^T Zi|^2。"""
    return float(np.sum((B.T @ Z.real) ** 2) + np.sum((B.T @ Z.imag) ** 2))


# ---------------- 生成子固有平面（方向） ----------------

def eigenplanes(sys_lr, Z, kmax=8):
    """K(θ(Z)) の固有平面（複素方向 w_j, σ降順, 単位ノルム）を上位 kmax 返す。

    JG の固有値 λ=±iσ。正の虚部側の固有ベクトル EV を辺空間へ持ち上げ w=W·EV, 正規化。
    戻り: sigmas(len<=N, 降順), Ws (M×k 複素, 列 j = w_j)。観測のみ。
    """
    sys_lr.set_theta(np.angle(Z))
    ev, EV = np.linalg.eig(sys_lr.J @ sys_lr.G)
    imag = ev.imag
    order = np.argsort(imag)[::-1]            # σ 降順（正の虚部が大きい順）
    sig = []
    cols = []
    for idx in order:
        if imag[idx] <= 1e-14:
            continue
        w = sys_lr.w(EV[:, idx].astype(complex))
        nw = np.linalg.norm(w)
        if nw == 0:
            continue
        w = w / nw
        sig.append(float(imag[idx]))
        cols.append(w)
        if len(cols) >= kmax:
            break
    if cols:
        Ws = np.column_stack(cols)
    else:
        Ws = np.zeros((sys_lr.m, 0), dtype=complex)
    return np.array(sig), Ws


def plane_energies(Ws, Z):
    """各固有平面 j への Z 占有エネルギー E_j = |<w_j,Z>|^2/|Z|^2（複素内積）。"""
    tot = float(np.real(np.conj(Z) @ Z))
    if Ws.shape[1] == 0 or tot == 0:
        return np.zeros(Ws.shape[1])
    proj = np.abs(Ws.conj().T @ Z) ** 2
    return proj / tot


def match_directions(Ws_prev, Ws_now, k=4):
    """w_j(t) を w(t-1) に貪欲照合（|overlap|最大・符号補正）。

    戻り: continuity[k]（|<w_j(t),matched(t-1)>|, 無ければ NaN）。
    順序は Ws_now の列順（σ降順）を direction j とみなし、直前の同順位方向との内積を取る。
    """
    kk = min(k, Ws_now.shape[1])
    cont = np.full(k, np.nan)
    if Ws_prev.shape[1] == 0:
        return cont
    for j in range(kk):
        wj = Ws_now[:, j]
        # 直前の全方向との |overlap| 最大を採用（順序交換・符号に頑健）
        ov = np.abs(Ws_prev.conj().T @ wj)
        if ov.size:
            cont[j] = float(np.max(ov))
    return cont


# ---------------- 流束（機械的 proxy, 報告書に明記） ----------------

def fluxes(E_prev, E_now):
    """平面占有エネルギー E_j の step 差に基づく機械的流束 proxy。

    net_flux_dj = E_j(now) - E_j(prev)。
    flux_parent_to_d3 = max(0, ΔE3) かつ ΔE_parent<0 のとき min(|ΔE1|,ΔE3), else 0。
    flux_d3_to_parent = max(0,-ΔE3) かつ ΔE_parent>0 のとき min(ΔE1,|ΔE3|), else 0。
    flux_d3_to_d4 = max(0, ΔE4) かつ ΔE3<0 のとき min(|ΔE3|,ΔE4), else 0。
    parent=平面1, d3=平面3, d4=平面4（0基準の index は 0,2,3）。存在しない平面は 0。
    """
    def g(E, i):
        return float(E[i]) if (E is not None and len(E) > i) else 0.0
    d1 = g(E_now, 0) - g(E_prev, 0)
    d3 = g(E_now, 2) - g(E_prev, 2)
    d4 = g(E_now, 3) - g(E_prev, 3)
    f_p_d3 = min(abs(d1), d3) if (d3 > 0 and d1 < 0) else 0.0
    f_d3_p = min(d1, abs(d3)) if (d3 < 0 and d1 > 0) else 0.0
    f_d3_d4 = min(abs(d3), d4) if (d4 > 0 and d3 < 0) else 0.0
    return dict(flux_parent_to_d3=f_p_d3, flux_d3_to_parent=f_d3_p, flux_d3_to_d4=f_d3_d4,
                net_flux_d3=d3, net_flux_d4=d4)


# ---------------- 適応ノイズ床 ----------------

def noise_floor(closure_abs, norm_error, min_nonzero_abs, prec_eps,
                quant_corr, retract_corr, repro_diff):
    """η_noise(t) = 記載量の最大（§10.1）。repro_diff は無ければ 0。"""
    return float(max(closure_abs, norm_error, min_nonzero_abs, prec_eps,
                     quant_corr, retract_corr, repro_diff))


# ---------------- 診断量 ----------------

def closure_terms(Z):
    ztz = complex(Z @ Z)                       # Z^T Z（零二乗閉鎖）
    return float(ztz.real), float(ztz.imag), float(abs(ztz))


def norm_error(Z):
    return float(abs(float(np.real(np.conj(Z) @ Z)) - 1.0))


def hermitian_antisym_norm(sys_lr, Z):
    """N が小さいとき K+K^T の Frobenius ノルム（反対称性診断）。大 N は NaN。"""
    n = sys_lr.n
    if sys_lr.m > 120:
        return float("nan")
    sys_lr.set_theta(np.angle(Z))
    ea, eb = sys_lr.ea, sys_lr.eb
    th = np.angle(Z)
    m = sys_lr.m
    A = np.zeros((m, m))
    for i in range(m):
        share = (ea == ea[i]) | (ea == eb[i]) | (eb == ea[i]) | (eb == eb[i])
        A[i, share] = 1.0
    np.fill_diagonal(A, 0.0)
    K = A * np.sin(th[None, :] - th[:, None])
    return float(np.linalg.norm(K + K.T))


def component_stats(Z):
    absZ = np.abs(Z)
    nz = absZ[absZ > 0]
    return dict(
        min_nonzero_abs_component=(float(nz.min()) if nz.size else 0.0),
        median_nonzero_abs_component=(float(np.median(nz)) if nz.size else 0.0),
        max_abs_component=float(absZ.max()),
        nonzero_real_count=int(np.sum(Z.real != 0.0)),
        nonzero_imag_count=int(np.sum(Z.imag != 0.0)),
    )
