#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""全 τ のスケール系列を計算してキャッシュする v1

楕円体図に添えるスケール履歴（対数グラフ）の元データを作る。
出力: scale_series_{stem}_{side}_v1.npz
  tau        : 保存された τ
  t, R, Q    : 上位3主軸から作る楕円体の半軸
  r_rms_full : 全次元の二乗平均半径（閉鎖が保存する量）
  r_rms_3d   : 上位3主軸に射影した二乗平均半径（観測されるスケール）
  top3       : 上位3主軸が占める割合
  rank, nimag: グラム行列のランクと負固有値の本数
  spec       : 各 τ の固有値スペクトル（符号付き √λ、大きい順、自明ゼロを除く N-1 本）
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent


def compute(stem: str, side: str) -> Path:
    X = np.load(HERE / f"dump_C2_{stem}_{side}_v1.npy", mmap_mode="r")
    meta = np.load(HERE / f"dump_meta_{stem}_{side}_v1.npz")
    taus = meta["dump_taus"] if "dump_taus" in meta.files else np.arange(X.shape[0])
    M = X.shape[1]
    N = int(round((1 + np.sqrt(1 + 8 * M)) / 2))
    ia, ib = np.triu_indices(N, k=1)
    J = np.eye(N) - np.ones((N, N)) / N
    out = {k: np.zeros(len(taus)) for k in
           ("t", "R", "Q", "r_rms_full", "r_rms_3d", "top3", "rank", "nimag")}
    spec = np.zeros((len(taus), N - 1))   # 符号付き √λ（大きい順）
    for f in range(len(taus)):
        C = np.asarray(X[f])
        d = np.abs(C.reshape(M, -1).sum(axis=1))
        D2 = np.zeros((N, N)); D2[ia, ib] = d ** 2; D2[ib, ia] = d ** 2
        B = -0.5 * J @ D2 @ J; B = 0.5 * (B + B.T)
        lam, U = np.linalg.eigh(B)
        o = np.argsort(-lam); lam, U = lam[o], U[:, o]
        # 二重中心化の自明ゼロ（固有ベクトル ∝ 1）を除く。単に末尾を切ると
        # 最大の虚方向を落として自明ゼロを主軸として描いてしまう。
        j0 = int(np.argmax(np.abs(U.T @ (np.ones(N) / np.sqrt(N)))))
        keep = np.ones(N, dtype=bool); keep[j0] = False
        lam_nt = lam[keep]
        sc = max(1.0, abs(lam[0]))
        pos = np.maximum(lam, 0.0)
        V3 = (U * np.sqrt(pos))[:, :3]
        ev = np.sort(np.linalg.eigvalsh(V3.T @ V3))[::-1]
        semi = np.sqrt(np.maximum(ev, 0.0) * 3.0 / N)
        out["t"][f], out["R"][f], out["Q"][f] = semi
        out["r_rms_full"][f] = np.sqrt(pos.sum() / N)
        out["r_rms_3d"][f] = np.sqrt(pos[:3].sum() / N)
        out["top3"][f] = pos[:3].sum() / pos.sum()
        out["rank"][f] = int((lam_nt > 1e-10 * sc).sum())
        out["nimag"][f] = int((lam_nt < -1e-10 * sc).sum())
        spec[f] = np.sign(lam_nt) * np.sqrt(np.abs(lam_nt))
        if f % 500 == 0:
            print(f"  {f}/{len(taus)}", flush=True)
    p = HERE / f"scale_series_{stem}_{side}_v1.npz"
    np.savez_compressed(p, tau=taus, spec=spec, **out)
    print(f"  -> {p.name}")
    return p


if __name__ == "__main__":
    stem = sys.argv[1]
    for side in (sys.argv[2:] or ["m", "v"]):
        print(f"=== {stem} {side} ===")
        compute(stem, side)
