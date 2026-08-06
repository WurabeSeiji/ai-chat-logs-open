#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""周期表整合バッテリー: 柱実験9本を局所θエンジンで無改変追試

方法（事前記録）: 各実験スクリプトは自前で exact モジュールをロードし
ex.collision_step_exact を呼ぶ。本ハーネスはスクリプト本文を一切変更せず、
ロード後にモジュール属性 collision_step_exact を局所θ版に差し替えて main() を
実行する（engine=global は無差し替え＝公開系の再現確認）。
結果JSONは本フォルダに書かれ、公開フォルダの正本は不変。

局所θ版: P10急峻（フェルミオンマスク=偶|k|≥4・元の分類と同一）/
P12滑らか（パリティ×IRロールオフexp(−(|k|/3)⁴)）。ro は大域読出しを
そのまま返す（報告用・力学には局所θを使用）。
使い方: python3 run_consistency_battery_v1.py <engine> <script_stem>
"""
from __future__ import annotations
import importlib.util, sys, io, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent

def make_local_step(exmod, smooth):
    toy = exmod.toy
    cache = {}
    def local_step(a, b, sp):
        n, ne = sp.chi_grid_n, sp.eta_grid_n
        if "W" not in cache:
            k = np.rint(np.fft.fftfreq(n, d=1.0 / n)).astype(int)
            if smooth:
                L = np.exp(-((np.abs(k) / 3.0) ** 4))
                Wf = ((k % 2) == 0).astype(float) * (1.0 - L)
            else:
                Wf = ((np.abs(k) >= 4) & ((np.abs(k) % 2) == 0)).astype(float)
            cache["W"] = (Wf, 1.0 - Wf)
        Wf, Wb = cache["W"]
        ro = toy.theta_from_ab(a, b, sp)          # 報告用（大域読出し）
        a2 = a.reshape(n, ne); b2 = b.reshape(n, ne)
        Fa = np.fft.fft(a2, axis=0); Fb = np.fft.fft(b2, axis=0)
        f = (np.sum(np.abs(np.fft.ifft(Fa * Wf[:, None], axis=0)) ** 2, axis=1)
             + np.sum(np.abs(np.fft.ifft(Fb * Wf[:, None], axis=0)) ** 2, axis=1))
        bo = (np.sum(np.abs(np.fft.ifft(Fa * Wb[:, None], axis=0)) ** 2, axis=1)
              + np.sum(np.abs(np.fft.ifft(Fb * Wb[:, None], axis=0)) ** 2, axis=1))
        th = np.arctan2(np.sqrt(f), np.sqrt(bo + 1e-300))
        c, s_ = np.cos(th)[:, None], np.sin(th)[:, None]
        a2, b2 = c * a2 - s_ * b2, s_ * a2 + c * b2
        phi = 2.0 * (np.sin(th) ** 2)[:, None] * np.imag(np.conj(b2) * a2)
        cp, sp_ = np.cos(phi), np.sin(phi)
        a2, b2 = cp * a2 - sp_ * b2, sp_ * a2 + cp * b2
        return a2.reshape(-1), b2.reshape(-1), ro
    return local_step

def run_one(stem, engine):
    path = HERE / f"{stem}.py"
    spec = importlib.util.spec_from_file_location(f"exp_{stem}_{engine}", path)
    mod = importlib.util.module_from_spec(spec); sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    if engine != "global":
        mod.ex.collision_step_exact = make_local_step(mod.ex, engine == "smooth")
    t0 = time.time()
    mod.main()
    print(f"### {stem} [{engine}] 完了 {time.time()-t0:.0f}s")

if __name__ == "__main__":
    run_one(sys.argv[2], sys.argv[1])
