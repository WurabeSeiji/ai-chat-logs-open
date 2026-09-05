#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""終状態の自己無撞着占有ランクの検定（読み出しのみ）。
理論的必要条件: 時計1回転写像の厳密固定点 Z[s+N]=Z[s] は exp(2πK)Z=Z を要求し、
これは Z が自分自身の生成子 K(arg Z) の**整数 σ 固有平面のみ**に占有されることと同値。
N=3 の符号反転固定点 Z[s+3]=−Z[s] は exp(2πK)Z=−Z ⟺ **半整数 σ のみ**の占有と同値。
検定: 各 N の step10000 状態について H=iK(arg Z) を密構成・固有分解し、
占有スペクトル（|V†Z|²）の (a) 有効モード数（参加比）、(b) 占有固有値の整数/半整数
からの最大偏差、(c) 上位占有モードの (w, 占有率) を機械算出する。
step0（親: 占有1モード・無理数 σ）との対比で「ランクは増えたか」に答える。
出力: check_final_selfconsistency_v1.json"""
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

BASE = Path(__file__).resolve().parent
ROOT = BASE.parents[5]
ENGINE = ROOT / "時間軸Q軸とフェルミオンの生成構造" / "検証_対照実験" / "第5論文原本_自発的分裂予備実験_v1"
sys.path.insert(0, str(ENGINE))
from run_n_scaling_lowrank_v1 import LowRankSystem

SWEEP = BASE.parent / "N3_N40_long10000_20260905"
NS = [3, 4, 5, 6, 7, 10, 40]

ledger = {}
with open(SWEEP / "SHA256SUMS.txt") as fh:
    for line in fh:
        parts = line.split()
        if len(parts) == 2:
            ledger[parts[1]] = parts[0]

out = {}
for N in NS:
    rel = f'results/hm_N{N}_den_{N}_states_10000.npz'
    h = hashlib.sha256(open(SWEEP / rel, 'rb').read()).hexdigest()
    assert ledger[rel] == h, f'INPUT GATE FAIL: {rel}'
    Z = np.asarray(np.load(SWEEP / rel)['Z'][10000], dtype=np.complex128)
    M = len(Z)
    sys_lr = LowRankSystem(N)
    sys_lr.set_theta(np.angle(Z))
    K = np.column_stack([sys_lr.kmatvec(np.eye(M)[:, j]) for j in range(M)])
    Hm = 1j * K  # エルミート（K 実反対称）
    w, V = np.linalg.eigh(Hm)
    occ = np.abs(V.conj().T @ Z) ** 2
    occ = occ / occ.sum()
    # 有効モード数（参加比）と占有モードの一覧
    reff = float(1.0 / np.sum(occ ** 2))
    idx = np.argsort(-occ)
    top = [(float(w[i]), float(occ[i])) for i in idx[:8] if occ[i] > 1e-12]
    # 占有加重の整数/半整数偏差
    dev_int = float(np.sum(occ * np.abs(w - np.round(w))))
    dev_half = float(np.sum(occ * np.abs(2 * w - np.round(2 * w))) / 2)
    n_occ = int(np.sum(occ > 1e-9))
    out[str(N)] = {'M': M, 'participation_modes': reff, 'n_modes_occ_gt_1e-9': n_occ,
                   'occ_weighted_dev_from_integer': dev_int,
                   'occ_weighted_dev_from_half_integer': dev_half,
                   'top_modes_w_occ': top}
    print(f"N={N}: 有効モード数={reff:.3f} 占有>1e-9 は {n_occ} 本 | "
          f"整数偏差={dev_int:.2e} 半整数偏差={dev_half:.2e} | 上位: "
          + ', '.join(f'(w={a:.4f}, {b:.3f})' for a, b in top[:4]), flush=True)

with open(BASE / 'check_final_selfconsistency_v1.json', 'w') as fh:
    json.dump(out, fh, indent=2)
print('ALL DONE')
