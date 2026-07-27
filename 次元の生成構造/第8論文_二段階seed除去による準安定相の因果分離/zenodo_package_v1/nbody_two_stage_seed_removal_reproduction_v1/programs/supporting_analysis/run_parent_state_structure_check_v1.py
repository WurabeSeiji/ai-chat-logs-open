#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""第8論文 草稿用補助解析 v1: 親状態（初期状態）の構造検査

冒頭の主張の追加段落「初期二方向状態の自発性」を支える実測を固定する。

検査内容:
  1. 零二乗閉鎖 Z^T Z = 0 の代数的帰結の数値確認:
     v = X + iY について ||X|| = ||Y||, X・Y = 0
  2. make_parent が生成する親状態 v の成分構造:
     全 M 成分が非零か、[X Y] の実 rank は 2 か
  3. 親状態の自己無撞着性: K(arg v) v = -i*sigma1*v の残差、
     不変平面条件 ||K X - sigma1 Y||, ||K Y + sigma1 X||
  4. 初期生成子 K(arg v) の空間分解の次元: 親平面 2 / 回転補空間 / 核

原本コードは第7論文原本 run_n_scaling_lowrank_v1.py を SHA-256 照合の上
read-only import する（Stage A2a と同一の固定原本・同一 PRNG seed）。
状態の時間発展は行わない。

出力: ../reports/parent_state_structure_check_v1.md
"""
import hashlib
import importlib.util
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.normpath(os.path.join(HERE, "..", "reports"))
os.makedirs(OUT, exist_ok=True)

ORIGINAL = (
    "/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/"
    "マイドライブ/OneDrive/GitHub/ai-chat-logs-open/時間軸Q軸とフェルミオンの生成構造/"
    "検証_対照実験/第5論文原本_自発的分裂予備実験_v1/run_n_scaling_lowrank_v1.py"
)
EXPECTED_SHA = "ba0fc19b03caf06d16e97e2cd5da499ed2ed8e95288f6dbf2062bedcaf11176d"
PARENT_PRNG_SEED = 40265722  # Stage A2a と同一
PARENT_ITERS = 1200
PARENT_TOL = 1e-12

sha = hashlib.sha256(open(ORIGINAL, "rb").read()).hexdigest()
if sha != EXPECTED_SHA:
    print("SHA-256 MISMATCH:", sha)
    sys.exit(1)

spec = importlib.util.spec_from_file_location("lowrank_original", ORIGINAL)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

lines = []
lines.append("# 草稿用補助解析 v1: 親状態（初期状態）の構造検査")
lines.append("")
lines.append("- 原本: `run_n_scaling_lowrank_v1.py`（SHA-256 `%s...` 照合 VERIFIED）" % sha[:16])
lines.append("- PRNG seed: `default_rng(%d)`、iters=%d、tol=%.0e（Stage A2a と同一）" % (PARENT_PRNG_SEED, PARENT_ITERS, PARENT_TOL))
lines.append("- 時間発展なし。親状態の構造検査のみ。")
lines.append("")

for n in (5, 40):
    rng = np.random.default_rng(PARENT_PRNG_SEED)
    sys_lr = mod.LowRankSystem(n)
    v, residual, sigma_spec = mod.make_parent(sys_lr, rng, iters=PARENT_ITERS, tol=PARENT_TOL)
    sigma_spec = np.atleast_1d(np.asarray(sigma_spec, dtype=float))
    sigma1 = float(sigma_spec[0])
    M = v.size
    X = v.real
    Y = v.imag
    nonzero = int(np.count_nonzero(np.abs(v) > 0.0))
    min_abs = float(np.min(np.abs(v)))
    rank_xy = int(np.linalg.matrix_rank(np.column_stack([X, Y])))
    # 原本の kmatvec（set_theta 済み生成子の作用）から密行列 K を構成する。
    # kmatvec は複素ベクトルに実反対称 K を成分ごとに作用させる。
    sys_lr.set_theta(np.angle(v))
    Kd = np.zeros((M, M))
    eye = np.eye(M)
    for j in range(M):
        Kd[:, j] = sys_lr.kmatvec(eye[:, j].astype(complex)).real
    antisym_err = float(np.linalg.norm(Kd + Kd.T))
    # 実測カイラリティ: iKv = mu v の mu（|mu| = sigma1）。符号は構成に依存するため実測する。
    mu = float(np.real(np.conj(v) @ (1j * (Kd @ v))))
    rx = float(np.linalg.norm(Kd @ X - mu * Y))
    ry = float(np.linalg.norm(Kd @ Y + mu * X))
    inv_plane = (rx, ry)
    sv = np.linalg.svd(Kd, compute_uv=False)
    tol_sv = sv[0] * 1e-10
    rank_K = int(np.sum(sv > tol_sv))
    ker_dim = M - rank_K
    rot_dim = rank_K - 2
    eig_dims = (2, rot_dim, ker_dim)
    lines.append("## N=%d（M=%d）" % (n, M))
    lines.append("")
    lines.append("| 量 | 値 |")
    lines.append("|:--|:--|")
    lines.append("| 非零複素成分数 / M | %d / %d |" % (nonzero, M))
    lines.append("| min |v_e| | %.6e |" % min_abs)
    lines.append("| rank[Re v, Im v] | %d |" % rank_xy)
    lines.append("| ‖Re v‖ | %.15f |" % float(np.linalg.norm(X)))
    lines.append("| ‖Im v‖ | %.15f |" % float(np.linalg.norm(Y)))
    lines.append("| Re v・Im v | %.6e |" % float(X @ Y))
    lines.append("| |v^T v|（零二乗閉鎖） | %.6e |" % abs(complex(v @ v)))
    lines.append("| |v†v − 1| | %.6e |" % abs(float(np.vdot(v, v).real) - 1.0))
    lines.append("| 固有モード残差 | %.17e |" % residual)
    lines.append("| sigma1 | %.15f |" % sigma1)
    lines.append("| 実測カイラリティ mu（iKv=mu v, |mu|=σ₁） | %.15f |" % mu)
    lines.append("| ‖K + K^T‖（反対称性） | %.6e |" % antisym_err)
    lines.append("| ‖K X − μY‖（不変平面残差） | %.6e |" % inv_plane[0])
    lines.append("| ‖K Y + μX‖（不変平面残差） | %.6e |" % inv_plane[1])
    lines.append("| 初期 K の分解次元（親平面/回転補/核） | %d / %d / %d |" % eig_dims)
    lines.append("")

lines.append("## 代数的帰結の確認")
lines.append("")
lines.append("Z = X + iY に対し Z^T Z = ‖X‖² − ‖Y‖² + 2i X・Y であるから、")
lines.append("零二乗閉鎖 Z^T Z = 0 は ‖X‖ = ‖Y‖ かつ X・Y = 0 と同値である。")
lines.append("上表はこの帰結（等ノルム・直交）と、全 M 成分非零のまま実 rank 2 に")
lines.append("集約される事実を、make_parent 出力に対して数値確認したものである。")
lines.append("")

out_path = os.path.join(OUT, "parent_state_structure_check_v1.md")
with open(out_path, "w") as fh:
    fh.write("\n".join(lines))
print("\n".join(lines))
print("written:", out_path)
