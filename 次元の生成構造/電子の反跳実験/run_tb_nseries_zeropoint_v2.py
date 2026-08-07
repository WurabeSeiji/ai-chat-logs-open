#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""N系列ゼロ点対照ラン（標準環境 Nn=16）— テストベッド N=12 選択の視覚的根拠

目的: N=4, 6, 12, 20, 50, 100 を同一条件（Nn=16・Nη=8・δ=10⁻²・T=4000）で
通し創世させ、同一形式の図で並べる。N=12 を標準テストベッドに選ぶ根拠を
(a) 時代構造 (b) 閉塞の健全性 (c) 集団時計の鋭さ の三点で視覚化する。

各 N について 2 枚:
 A) ゼロ点対照パネル（4段）
    1. f₂（空間形成史・物質 vs 真空）
    2. f_seed（物質分率）
    3. pr_n（双対占有・Nn=16 が一様＝可算ドメインなしの署名）
    4. **ゼロ閉塞 |Σxₙ²| の直接計測**——統一万能読出し関数を通さない監査量。
       宇宙の内側の観測者には読めない神の視点の量であり、読出しの健全性とは
       独立に力学の健全性を示す。総パワーで規格化して表示する。
 B) 時間–振動数ヒートマップ（物質・真空の2段）
    全関係波を合成した信号 s(τ) = Σ_{e,k,η} C2[e,k,η] の短時間フーリエ変換。
    集団時計 ω=π/72 は 1/144 cycles/step の線として現れる。

記号: τ(step) は処理ステップ（座標時間 t ではない。集団時計が毎ステップ
π/72 だけ位相を進めるので、ステップ数は固有時間の目盛に比例する）。

使い方: python3 run_tb_nseries_zeropoint_v2.py [N ...]（省略時 4 6 12 20 50 100）
"""
from __future__ import annotations
import importlib.util, json, sys, time
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams["font.family"] = ["Hiragino Sans", "sans-serif"]
plt.rcParams["axes.unicode_minus"] = False

HERE = Path(__file__).resolve().parent
UF = HERE.parent / "統一万能関数_v1"


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


ui = load("ui_ns", UF / "unified_interaction_v1.py")
G = load("ur_ns", UF / "unified_readout_v2.py")
S = load("sel_ns", UF / "selection_v1.py")

NN, NETA = 16, 8
T = 4000
DELTA = 1e-2
SEED = 2
WIN = (2000, 4000)
CLOCK_F = 1.0 / 144.0          # ω=π/72 → 1/144 cycles/step（柱1の集団時計）
NFFT, HOP = 1024, 64


def build(n, delta):
    m = n * (n - 1) // 2
    _, v, _, _, _, _, _, Z0c, wp0 = ui.abl.build_init(n, False)
    r2 = ui.gen3.make_parent(n, seed=SEED)
    Csec = np.fft.fft(r2.relation_waves, axis=1) / n
    s0 = Csec[:, 1] / np.linalg.norm(Csec[:, 1])
    C2_0 = np.zeros((m, NN, NETA), complex)
    C2_0[:, 2, 0] = Z0c
    if delta > 0:
        C2_0[:, 1, 0] = delta * s0
    p2 = C2_0[:, 2, 0].real / np.linalg.norm(C2_0[:, 2, 0].real)
    q2 = C2_0[:, 2, 0].imag - (C2_0[:, 2, 0].imag @ p2) * p2
    q2 = q2 / np.linalg.norm(q2)
    return ui.UnifiedEngine(n, C2_0, wp0), p2, q2


def run(n, delta):
    eng, p2, q2 = build(n, delta)
    H = {k: np.zeros(T) for k in ("f2", "f_seed", "pr_n", "closure", "P_tot")}
    sig = np.zeros(T, complex)
    carry = {"C_flat": None, "c_gen": None}
    for t in range(T):
        eng.step()
        C2 = eng.C2()
        pan = G.g_panel(C2, p2, q2, carry["C_flat"], carry["c_gen"])
        carry = pan["_carry"]
        H["f2"][t] = pan["f2"]
        H["f_seed"][t] = pan["f_seed"]
        H["pr_n"][t] = pan["pr_n"]
        H["P_tot"][t] = pan["P_tot"]
        # --- 監査量（統一読出しの外・内部状態の直接計測）: 零閉塞 Σxₙ²
        H["closure"][t] = abs(complex(np.sum(C2 ** 2))) / max(pan["P_tot"], 1e-300)
        # --- 全波合成信号（時間–振動数図の入力）
        sig[t] = complex(np.sum(C2))
    return H, sig


def spectrogram(sig):
    """短時間フーリエ変換（複素信号・Hann窓）。戻り: 時刻中心・振動数・振幅"""
    w = np.hanning(NFFT)
    starts = np.arange(0, len(sig) - NFFT + 1, HOP)
    Z = np.zeros((NFFT, len(starts)))
    for j, s0 in enumerate(starts):
        seg = sig[s0:s0 + NFFT] * w
        Z[:, j] = np.abs(np.fft.fftshift(np.fft.fft(seg)))
    f = np.fft.fftshift(np.fft.fftfreq(NFFT))
    tc = starts + NFFT / 2
    return tc, f, Z


def fig_panels(n, Hm, Hv):
    ts = np.arange(1, T + 1)
    fig, ax = plt.subplots(4, 1, figsize=(9, 11), sharex=True)
    ax[0].semilogy(ts, np.maximum(Hm["f2"], 1e-18), lw=0.8, label="物質 δ=10⁻²")
    ax[0].semilogy(ts, np.maximum(Hv["f2"], 1e-18), "k--", lw=0.8, label="真空 δ=0")
    ax[0].set_ylabel("f₂（空間形成史）"); ax[0].legend(fontsize=8)
    ax[0].set_title(f"N={n}（M={n*(n-1)//2}）ゼロ点対照ラン・Nn={NN}・"
                    f"統一G v2 を第0步から常時実行")
    ax[1].semilogy(ts, np.maximum(Hm["f_seed"], 1e-18), lw=0.8, label="物質")
    ax[1].semilogy(ts, np.maximum(Hv["f_seed"], 1e-18), "k--", lw=0.8, label="真空")
    ax[1].set_ylabel("f_seed（物質分率）"); ax[1].legend(fontsize=8)
    ax[2].plot(ts, Hm["pr_n"], lw=0.8, label="物質")
    ax[2].plot(ts, Hv["pr_n"], "k--", lw=0.8, label="真空")
    ax[2].axhline(NN, color="red", lw=0.8, label=f"Nn={NN}（一様＝可算ドメインなし）")
    ax[2].set_ylabel("pr_n（双対占有）"); ax[2].legend(fontsize=8)
    ax[3].semilogy(ts, np.maximum(Hm["closure"], 1e-20), lw=0.8, label="物質")
    ax[3].semilogy(ts, np.maximum(Hv["closure"], 1e-20), "k--", lw=0.8, label="真空")
    ax[3].set_ylabel("|Σxₙ²| / P_tot（監査量・直接計測）")
    ax[3].set_xlabel("τ（step）"); ax[3].legend(fontsize=8)
    fig.tight_layout()
    p = HERE / f"fig_tb_zeropoint_N{n}_Nn{NN}_v2.png"
    fig.savefig(p, dpi=130); plt.close(fig)
    return p


def fig_spectrogram(n, sm, sv):
    fig, ax = plt.subplots(2, 1, figsize=(9, 7), sharex=True, sharey=True)
    for a, sig, lab in ((ax[0], sm, "物質宇宙 δ=10⁻²"), (ax[1], sv, "真空宇宙 δ=0")):
        tc, f, Z = spectrogram(sig)
        sel = (f >= -0.005) & (f <= 0.03)
        im = a.pcolormesh(tc, f[sel], np.log10(np.maximum(Z[sel], 1e-20)),
                          shading="auto", cmap="magma")
        a.axhline(CLOCK_F, color="cyan", lw=0.9, ls="--",
                  label="1/144 cycles/step（ω=π/72・柱1集団時計）")
        a.set_ylabel(f"振動数（cycles/step）\n{lab}")
        a.legend(fontsize=8, loc="upper right")
        fig.colorbar(im, ax=a, label="log₁₀|振幅|")
    ax[1].set_xlabel("τ（step）")
    ax[0].set_title(f"N={n} 全波合成信号 s(τ)=Σ C2 の時間–振動数図（STFT・窓{NFFT}）")
    fig.tight_layout()
    p = HERE / f"fig_tb_spectrogram_N{n}_Nn{NN}_v1.png"
    fig.savefig(p, dpi=130); plt.close(fig)
    return p


def clock_sharpness(sig):
    """安定窓の合成信号スペクトルにおける集団時計線の鋭さ（Q値相当）と
    ピーク振動数。N とともに鋭くなる（~1/M）ことの定量。"""
    seg = sig[WIN[0]:WIN[1]]
    seg = seg - seg.mean()
    w = np.hanning(len(seg))
    Fq = np.fft.fftshift(np.fft.fftfreq(len(seg)))
    A = np.abs(np.fft.fftshift(np.fft.fft(seg * w)))
    band = (Fq > 0.0) & (Fq < 0.05)
    if not band.any() or A[band].max() <= 0:
        return None
    fb, Ab = Fq[band], A[band]
    i = int(np.argmax(Ab))
    half = Ab[i] / np.sqrt(2.0)
    lo = i
    while lo > 0 and Ab[lo] > half:
        lo -= 1
    hi = i
    while hi < len(Ab) - 1 and Ab[hi] > half:
        hi += 1
    width = float(fb[hi] - fb[lo])
    return {"f_peak": float(fb[i]), "fwhm": width,
            "Q": float(fb[i] / width) if width > 0 else None,
            "rel_dev_from_clock": float(abs(fb[i] - CLOCK_F) / CLOCK_F)}


def main():
    t0 = time.time()
    ns = [int(a) for a in sys.argv[1:]] or [4, 6, 12, 20, 50, 100]
    out = {"env": {"Nn": NN, "Neta": NETA, "T": T, "delta": DELTA,
                   "seed": SEED, "window": WIN, "clock_f": CLOCK_F}, "N": {}}
    for n in ns:
        t1 = time.time()
        Hm, sm = run(n, DELTA)
        Hv, sv = run(n, 0.0)
        pa = fig_panels(n, Hm, Hv)
        pb = fig_spectrogram(n, sm, sv)
        sh_m, sh_v = clock_sharpness(sm), clock_sharpness(sv)
        cw = slice(*WIN)
        rec = {"M": n * (n - 1) // 2,
               "crossing_matter": int(np.argmax(Hm["f2"] > 0.05)) + 1,
               "crossing_vacuum": int(np.argmax(Hv["f2"] > 0.05)) + 1,
               "f2_win_mean": float(np.mean(Hm["f2"][cw])),
               "f2_win_driftrel": float(abs(np.polyfit(np.arange(WIN[1]-WIN[0]),
                                                       Hm["f2"][cw], 1)[0]
                                            * (WIN[1]-WIN[0])
                                            / np.mean(Hm["f2"][cw]))),
               "fseed_final": float(Hm["f_seed"][-1]),
               "fseed_vacuum_final": float(Hv["f_seed"][-1]),
               "prn_win_mean": float(np.mean(Hm["pr_n"][cw])),
               "closure_max_matter": float(np.max(Hm["closure"])),
               "closure_max_vacuum": float(np.max(Hv["closure"])),
               "clock_matter": sh_m, "clock_vacuum": sh_v,
               "figs": [pa.name, pb.name]}
        out["N"][n] = rec
        cm = sh_m or {}
        print(f"N={n:3d} M={rec['M']:4d}: crossing 物質={rec['crossing_matter']:4d}/"
              f"真空={rec['crossing_vacuum']:4d} 閉塞max={rec['closure_max_matter']:.1e} "
              f"時計 f={cm.get('f_peak', float('nan')):.6f}"
              f"(1/144={CLOCK_F:.6f}・偏差{cm.get('rel_dev_from_clock', float('nan')):.1e}) "
              f"Q={cm.get('Q') if cm.get('Q') is None else round(cm['Q'], 1)} "
              f"[{time.time()-t1:.0f}s]")
        np.savez_compressed(HERE / f"tb_nseries_N{n}_Nn{NN}_v2.npz",
                            **{f"m_{k}": Hm[k] for k in Hm},
                            **{f"v_{k}": Hv[k] for k in Hv},
                            sig_m=sm, sig_v=sv)
        (HERE / "result_tb_nseries_zeropoint_v2.json").write_text(
            json.dumps(out, indent=1, ensure_ascii=False, default=float))
    out["runtime_sec"] = time.time() - t0
    (HERE / "result_tb_nseries_zeropoint_v2.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False, default=float))
    print(f"完了 {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
