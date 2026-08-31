#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""パス9（読出し専用）：合成波 W(t)=Σ_e z_e(t) による基底振動数・基底波長の読出し。
根拠（波長の導出方法.md の改訂方針・木原）：M 個の複素数から無ラベルで作れる正準スカラーは
Σz²（閉塞で零＝静か）と Σz（閉塞に殺されない）の 2 つ。剛体回転 z=e^{−iμt}v なら W(t)=e^{−iμt}Σv_e
——合成波は系の固有時計そのもので回る。基底振動数＝W の最強スペクトル線、基底波長＝その逆数（λ0 任意）。
内部観測量は線と線の比のみ。

走行前に固定する予測（parent_checks / pass8 の既存値から）：
 P1 親窓（step 0〜2047）：単一線、位相進み ω_early = −Δμ_parent（実測済み時計と一致）
 P2 終窓（step 37952〜40000）：単一線、ω_late = −Δμ_final（星なら −Δ·N(N−2)r̄²/2、回転方向が親と逆転）
 P3 基底振動数比 ω_late/ω_early = μ_final/μ_parent（既測の有理比 −N(N−2)/4 を無ラベルのスカラーが再現するか）
 P4 転移窓（t50 中心）：複数線（子の分裂＝うなり）
 P5 合成波の明るさ |W|/Σ|z| の終値 = |n₊−n₋|/(N−1)（pass8 の符号内訳から。ne_N15 は 7/7 → 暗い星 W≈0）
測定は (a) 位相進み ω = arg Σ_t W̄(t)W(t+1)（|W|² 重み、単一線なら機械精度）
     (b) Hann 窓 FFT の線構成（重み≥1% の線数、上位 5 線。パス8 と同じ 1% 規約）。調整パラメータなし。"""
import os, csv, json, math
import numpy as np
import sys as _sys
_sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams['font.family'] = ['Hiragino Sans', 'Arial Unicode MS', 'sans-serif']
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data"); FIG = os.path.join(ROOT, "figures")
L = 124; DELTA = 2*math.pi/L
WIN = 2048

def omega_phase(W):
    """|W|² 重み付き平均位相進み／step（単一線なら厳密）。"""
    num = np.sum(np.conj(W[:-1])*W[1:])
    return float(np.angle(num)) if abs(num) > 0 else float("nan")

def fft_lines(W):
    """Hann 窓 FFT。重み≥最大線の1% の線数と上位 5 線 (freq[rad/step], power割合)。"""
    n = len(W)
    h = np.hanning(n)
    F = np.fft.fft(W*h); P = np.abs(F)**2
    fr = np.fft.fftfreq(n)*2*math.pi
    # 隣接ビンより大きい局所ピークのみ
    pk = [(P[i], fr[i]) for i in range(n) if P[i] >= P[(i-1) % n] and P[i] >= P[(i+1) % n]]
    pk.sort(reverse=True)
    pmax = pk[0][0] if pk else 1.0
    lines = [(f, p/pmax) for p, f in pk if p >= 0.01*pmax]
    return len(lines), lines[:5]

rows = []
spectra_cache = {}
for tag in sorted(os.listdir(DATA)):
    if tag == "reference": continue
    fp = os.path.join(DATA, tag, "states_treatment.npz")
    cj = os.path.join(DATA, tag, "parent_checks.json")
    if not (os.path.exists(fp) and os.path.exists(cj)): continue
    c = json.load(open(cj)); N = c["N"]
    Z = np.load(fp)["Z"]
    W = Z.sum(axis=1)
    absZ = np.abs(Z).sum(axis=1)
    bright = np.abs(W)/np.maximum(absZ, 1e-300)
    sj = json.load(open(os.path.join(DATA, tag, "summary.json")))
    a = None
    # 窓
    We = W[:WIN]; Wl = W[-WIN:]
    om_e = omega_phase(We); om_l = omega_phase(Wl)
    nl_e, lines_e = fft_lines(We); nl_l, lines_l = fft_lines(Wl)
    # 転移窓（t50 があれば）
    ts = np.genfromtxt(os.path.join(DATA, tag, "treatment_linear124_amplitude_aware_timeseries.csv"),
                       delimiter=",", names=True)
    f = ts["H_perp"]/ts["H_total"]; i50 = np.where(f >= 0.5)[0]
    t50 = int(ts["step"][i50[0]]) if len(i50) else None
    if t50 and WIN//2 < t50 < len(W)-WIN//2:
        Wm = W[t50-WIN//2:t50+WIN//2]
        om_m = omega_phase(Wm); nl_m, lines_m = fft_lines(Wm)
    else:
        om_m = float("nan"); nl_m = None; lines_m = []
    mu_p = c["mu_new"]
    # 終状態の μ を最終状態から直接計算（H(z)z の Rayleigh 商）
    from interference_dynamics import hermitian_H
    from common import adjacency
    Zf = Z[-1]; A = adjacency(N); Hz = hermitian_H(Zf, A) @ Zf
    mu_f = float((np.vdot(Zf, Hz)/np.vdot(Zf, Zf)).real)
    pred_e = -DELTA*mu_p
    pred_l = -DELTA*mu_f if mu_f is not None else float("nan")
    ratio_meas = om_l/om_e if (np.isfinite(om_e) and abs(om_e) > 1e-12) else float("nan")
    ratio_pred = mu_f/mu_p if (mu_f is not None and abs(mu_p) > 1e-12) else float("nan")
    rows.append(dict(tag=tag, N=N, method=c["method"], t50=t50,
                     omega_early=om_e, omega_early_pred=pred_e,
                     omega_early_dev=abs(om_e-pred_e)/max(abs(pred_e), 1e-300),
                     omega_late=om_l, omega_late_pred=pred_l,
                     omega_late_dev=abs(om_l-pred_l)/max(abs(pred_l), 1e-300) if np.isfinite(pred_l) else float("nan"),
                     lambda_ratio_late_over_early=(om_e/om_l if abs(om_l) > 1e-12 else float("nan")),
                     omega_ratio_meas=ratio_meas, omega_ratio_pred=ratio_pred,
                     n_lines_early=nl_e, n_lines_mid=nl_m, n_lines_late=nl_l,
                     bright_first=float(bright[0]), bright_final=float(bright[-1]),
                     lines_mid=";".join(f"{f:+.5f}:{p:.3f}" for f, p in lines_m)))
    if tag in ("hm_N12", "mp_N13", "ne_N15", "hm_N8"):
        spectra_cache[tag] = dict(W=W, t50=t50, bright=bright)
    print(f"{tag}: ω早={om_e:+.6f}(予{pred_e:+.6f}) ω終={om_l:+.6f}(予{pred_l:+.6f}) "
          f"線数 早{nl_e}/転{nl_m}/終{nl_l} 明るさ終={bright[-1]:.3f}")

with open(os.path.join(ROOT, "results", "composite_wave_summary.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)

# 判定集計
ok_e = sum(1 for r in rows if r["omega_early_dev"] < 1e-3)
ok_l = sum(1 for r in rows if np.isfinite(r["omega_late_dev"]) and r["omega_late_dev"] < 1e-3)
multi = [r["tag"] for r in rows if r["n_lines_mid"] and r["n_lines_mid"] >= 2]
md = ["# 合成波 W(t)=Σz_e(t) による基底振動数・基底波長の読出し（パス9）", "",
      f"P1 親窓の基底振動数＝−Δμ_parent（偏差<1e-3）: {ok_e}/{len(rows)}",
      f"P2 終窓の基底振動数＝−Δμ_final（偏差<1e-3）: {ok_l}/{len(rows)}",
      f"P4 転移窓で複数線（子の分裂）: {len(multi)} 走行 → {', '.join(multi) if multi else 'なし'}", "",
      "| tag | ω早(実測) | ω早(予) | 偏差 | ω終(実測) | ω終(予) | 偏差 | λ終/λ早=ω早/ω終 | μ比(予) | 線数 早/転/終 | 明るさ終 |",
      "|---|---|---|---|---|---|---|---|---|---|---|"]
def fm(x, d=6):
    return "—" if x is None or (isinstance(x, float) and not np.isfinite(x)) else f"{x:+.{d}f}"
def fe2(x):
    return "—" if x is None or (isinstance(x, float) and not np.isfinite(x)) else f"{x:.1e}"
for r in rows:
    md.append("| " + " | ".join([r["tag"], fm(r["omega_early"]), fm(r["omega_early_pred"]), fe2(r["omega_early_dev"]),
                                 fm(r["omega_late"]), fm(r["omega_late_pred"]), fe2(r["omega_late_dev"]),
                                 fm(r["lambda_ratio_late_over_early"], 4), fm(r["omega_ratio_pred"], 4),
                                 f"{r['n_lines_early']}/{r['n_lines_mid'] if r['n_lines_mid'] is not None else '—'}/{r['n_lines_late']}",
                                 f"{r['bright_final']:.3f}"]) + " |")
with open(os.path.join(ROOT, "results", "composite_wave_reading.md"), "w") as fh:
    fh.write("\n".join(md)+"\n")

# fig11: 合成波の明るさ |W|/Σ|z| の時間発展（fig3 と同レイアウト）
col = {'mp': 'k', 'hm': '#1f77b4', 'ne': '#ff7f0e', 'rb': '#2ca02c'}
figA, axs = plt.subplots(4, 4, figsize=(17, 14)); axs = axs.ravel()
byN = {}
for tag in sorted(os.listdir(DATA)):
    if tag == "reference": continue
    fp = os.path.join(DATA, tag, "states_treatment.npz")
    if not os.path.exists(fp): continue
    Z = np.load(fp)["Z"]; W = Z.sum(axis=1)
    br = np.abs(W)/np.maximum(np.abs(Z).sum(axis=1), 1e-300)
    N = int(tag.split("N")[1]); byN.setdefault(N, []).append((tag.split("_")[0], br))
for i, N in enumerate(range(3, 17)):
    ax = axs[i]
    for m, br in byN.get(N, []):
        ax.plot(np.arange(0, len(br), 20), br[::20], color=col[m], lw=1.0, label=m)
    ax.set_ylim(0, 1.0); ax.set_title(f'N={N}', fontsize=10); ax.grid(alpha=.3); ax.legend(fontsize=6)
    if i % 4 == 0: ax.set_ylabel('|W|/Σ|z|（合成波の明るさ）')
    if i >= 10: ax.set_xlabel('step')
for j in range(14, 16): axs[j].axis('off')
figA.suptitle('合成波 W=Σz の明るさの時間発展：親=有限、星=|n₊−n₋|/(N−1)（暗い星は W≈0）', fontsize=11)
plt.tight_layout(rect=(0, 0, 1, 0.97)); plt.savefig(os.path.join(FIG, 'fig11_composite_brightness_grid.png'), dpi=130); plt.close()

# fig12: 代表スペクトル（hm_N12 早/転/終、mp_N13 転、ne_N15 終）
figB, axs = plt.subplots(1, 5, figsize=(20, 3.8))
panels = [("hm_N12", "early"), ("hm_N12", "mid"), ("hm_N12", "late"), ("mp_N13", "mid"), ("ne_N15", "late")]
for ax, (tag, win) in zip(axs, panels):
    d = spectra_cache[tag]; W = d["W"]; t50 = d["t50"]
    if win == "early": seg = W[:WIN]; ttl = f"{tag} 親窓 0–{WIN}"
    elif win == "late": seg = W[-WIN:]; ttl = f"{tag} 終窓"
    else: seg = W[t50-WIN//2:t50+WIN//2]; ttl = f"{tag} 転移窓 t50={t50}"
    F = np.fft.fftshift(np.fft.fft(seg*np.hanning(len(seg)))); fr = np.fft.fftshift(np.fft.fftfreq(len(seg))*2*math.pi)
    P = np.abs(F)**2; P /= P.max()
    ax.semilogy(fr, np.maximum(P, 1e-10), lw=0.8)
    ax.set_xlim(-0.3, 0.3); ax.set_ylim(1e-8, 2); ax.grid(alpha=.3)
    ax.set_title(ttl, fontsize=9); ax.set_xlabel('ω [rad/step]')
axs[0].set_ylabel('パワー（最大=1）')
figB.suptitle('合成波のスペクトル：親=単一線 → 転移=複数線（子の分裂） → 星=単一線（逆回転）', fontsize=10)
plt.tight_layout(rect=(0, 0, 1, 0.92)); plt.savefig(os.path.join(FIG, 'fig12_composite_spectra.png'), dpi=150); plt.close()
print("PASS9 OK")
