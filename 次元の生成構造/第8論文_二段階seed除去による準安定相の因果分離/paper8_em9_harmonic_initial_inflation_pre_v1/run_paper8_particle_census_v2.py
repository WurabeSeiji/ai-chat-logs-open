#!/usr/bin/env python3
"""粒子的グループ census v2：一覧表形式（1状態=1行×時点2行、CSV 同時保存）

v1 からの変更（木原氏指示 2026-08-04）:
    - 出力を固定幅の一覧表に変更:
      N | 状態 | 族 | 時点 | ベース波長 | 位相° | 倍音1..6 | B/F | 電荷(生)
    - 時点 = 初期（t=0 の閉包そのもの）と 緩和後（6000 tick）の2行で
      「初期在庫 → 溶解」を同じ表で読めるようにする
    - CSV（census_table_v2.csv）と JSON を保存（再現性）

列の定義:
    ベース波長: 偏差重み付き巻き数整列度 A_k = |Σ|dev|e^{ikθ}|/Σ|dev|（k≤M/2）
        の最大が 0.3 超のとき 1/k。超えなければ "—"（一様海）
    位相°: そのkのフレーム方位 arg(R_k)/k [deg]
    倍音1..6: A_1..A_6 の値（構造スペクトル。0..1）
    B/F: 行動的半周期テスト（E-M11 倍音時計 micro=段n、72tick、Re⟨Z|Z(T/2)⟩ の符号。
        F=奇/B=偶。パリティは発展で保存されるため状態につき1個）
    電荷(生): ベース波長塊の周波数/残余の周波数の生の比（緩和後のみ。
        一定にならなくてもそのまま記録——木原氏指示。塊なしは "—"）

規約: 種ラベルIFなし・read-only import・SHA-256記録。seed は v1 と同一。
"""

from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PAPER8 = HERE.parent
REPO = PAPER8.parent.parent
ABL = PAPER8 / "code" / "run_preliminary_seed_ablation_v1.py"
MPH = REPO / "次元の生成構造" / "make_parent_harmonic_unit_v1" / "make_parent_harmonic_v1.py"

spec = importlib.util.spec_from_file_location("abl_c2", ABL)
abl = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = abl
spec.loader.exec_module(abl)
spec2 = importlib.util.spec_from_file_location("mph_c2", MPH)
mph = importlib.util.module_from_spec(spec2)
sys.modules[spec2.name] = mph
spec2.loader.exec_module(mph)

T_RELAX = 6000
WIN = 500
T_HALF = 72
K_SHOW = 6
A_BASE = 0.3
FREQ_MIN = 1e-8
CFG = {5: {"H": 8, "seed": 40260801},
       40: {"H": 4, "seed": 40260802},
       300: {"H": 4, "seed": 40260803}}   # 2026-08-04 追記（既存実験と同一 seed）


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def spectrum_row(Z, m):
    """偏差重み付き巻き数スペクトル A_1..A_6、ベース波長、位相。"""
    amp = np.abs(Z); th = np.angle(Z)
    dev = amp - np.median(amp)
    w = np.abs(dev)
    k_hi = max(1, min(K_SHOW, m // 2))
    if w.sum() < 1e-8 * max(1.0, amp.sum()):
        return {"A": [0.0] * K_SHOW, "k0": None, "phase_deg": None, "s_circ": None,
                "split_rms": float(np.std(amp)), "w": w, "th": th, "dev": dev}
    A = []
    for k in range(1, K_SHOW + 1):
        A.append(float(abs(np.sum(w * np.exp(1j * k * th))) / w.sum())
                 if k <= k_hi else 0.0)
    k0 = int(np.argmax(A) + 1) if max(A) > A_BASE else None
    if k0 is not None:
        c = float(np.angle(np.sum(w * np.exp(1j * k0 * th))))
        phase = float(np.degrees(c / k0))
        s_circ = float(np.sum(dev * np.sin(2 * (k0 * th - c))))   # v1 と同一定義のC-odd量
    else:
        phase = None
        s_circ = None
    return {"A": A, "k0": k0, "phase_deg": phase, "s_circ": s_circ,
            "split_rms": float(np.std(amp)), "w": w, "th": th, "dev": dev}


def census_state(n, level, v0, wp, name, family):
    m = n * (n - 1) // 2
    row0 = spectrum_row(v0, m)                      # 初期
    sys_lr = abl.LowRankSystem(n)
    sys_lr.set_theta(np.angle(v0))
    Z = v0.copy()
    for _ in range(T_RELAX):
        for _ in range(level):
            Z, wp = abl.evolve(sys_lr, Z, wp)
    Zref = Z.copy()
    phases = np.zeros((WIN, m))
    for t in range(WIN):
        for _ in range(level):
            Z, wp = abl.evolve(sys_lr, Z, wp)
        phases[t] = np.angle(Z)
    Zh = Zref.copy(); wp2 = wp.copy()
    for _ in range(T_HALF):
        for _ in range(level):
            Zh, wp2 = abl.evolve(sys_lr, Zh, wp2)
    parity = "F" if complex(np.vdot(Zref, Zh)).real < 0 else "B"
    row1 = spectrum_row(Zref, m)                    # 緩和後
    # 電荷(生): 緩和後にベース波長塊があれば 塊周波数/残余周波数
    charge = None
    if row1["k0"] is not None:
        k0 = row1["k0"]
        w1, th1 = row1["w"], row1["th"]
        c = float(np.angle(np.sum(w1 * np.exp(1j * k0 * th1))))
        resid = np.abs(np.angle(np.exp(1j * (k0 * th1 - c))))
        member = (w1 > max(1e-3, 0.1 * w1.max())) & (resid < 0.8)
        if member.sum() >= 2 and (~member).sum() >= 1:
            fr = np.abs(np.polyfit(np.arange(WIN), np.unwrap(phases, axis=0), 1)[0])
            fS, fR = fr[member], fr[~member]
            if np.mean(fR) > FREQ_MIN:
                charge = float(np.mean(fS) / np.mean(fR))
    return [
        {"N": n, "state": name, "family": family, "time": "初期",
         "k0": row0["k0"], "phase_deg": row0["phase_deg"], "A": row0["A"],
         "split_rms": row0["split_rms"], "BF": parity, "charge_raw": None,
         "s_circ": row0["s_circ"],
         "pa": ("粒" if row0["s_circ"] > 0 else "反") if row0["s_circ"] else "—"},
        {"N": n, "state": name, "family": family, "time": "緩和後",
         "k0": row1["k0"], "phase_deg": row1["phase_deg"], "A": row1["A"],
         "split_rms": row1["split_rms"], "BF": parity, "charge_raw": charge,
         "s_circ": row1["s_circ"],
         "pa": ("粒" if row1["s_circ"] > 0 else "反") if row1["s_circ"] else "—"},
    ]


def fmt_row(r):
    wl = f"1/{r['k0']}" if r["k0"] else "—"
    ph = f"{r['phase_deg']:+7.1f}" if r["phase_deg"] is not None else "      —"
    As = " ".join(f"{a:5.2f}" for a in r["A"])
    ch = f"{r['charge_raw']:.6f}" if r["charge_raw"] is not None else "—"
    return (f"{r['N']:>3} | {r['state']:<10} | {r['family']:<7} | {r['time']:<3} | "
            f"{wl:>4} | {ph} | {As} | {r['BF']} | {r['pa']} | {ch}")


def main() -> None:
    t0 = time.time()
    ns = ([int(x) for x in sys.argv[1].split(",")] if len(sys.argv) > 1
          else list(CFG))
    print(f"粒子的グループ census v2 実行（N={ns}）")
    print(f"  import: ABL {sha256(ABL)[:16]}…  MPH {sha256(MPH)[:16]}…\n")
    rows = []
    for n in ns:
        cfg = CFG[n]
        _, v, _, _, _, _, _, Z0, wp0 = abl.build_init(n, False)
        rows += census_state(n, 1, Z0, wp0.copy(), "control", "control")
        Zh, info = mph.make_parent_harmonic(n, cfg["H"], cfg["seed"],
                                             iters=2000, restarts=10, tol=1e-12)
        for h in range(1, cfg["H"] + 1):
            lv = info["levels"][h - 1]
            fam = "N-1" if abs(lv["sigma1"] - (n - 1)) < 1e-9 else "broken"
            v0 = Zh[:, h - 1] * np.sqrt(cfg["H"])
            wp = np.random.default_rng(90000 + (h - 1)).normal(size=len(v0))
            rows += census_state(n, h, v0, wp, f"n={h}", fam)

    hdr = ("  N | 状態       | 族      | 時点 | 波長 |  位相°  | "
           "倍音1 倍音2 倍音3 倍音4 倍音5 倍音6 | B/F | 粒/反 | 電荷(生)")
    print(hdr); print("-" * len(hdr))
    cur = None
    for r in rows:
        if cur is not None and r["N"] != cur:
            print("-" * len(hdr))
        cur = r["N"]
        print(fmt_row(r))

    with open(HERE / "census_table_v2.csv", "w", newline="", encoding="utf-8") as fh:
        wtr = csv.writer(fh)
        wtr.writerow(["N", "state", "family", "time", "base_wavelength_k0",
                      "phase_deg", "A1", "A2", "A3", "A4", "A5", "A6",
                      "amp_split_rms", "BF", "particle_antiparticle_conv",
                      "s_circ", "charge_raw"])
        for r in rows:
            wtr.writerow([r["N"], r["state"], r["family"], r["time"],
                          r["k0"] if r["k0"] else "",
                          f"{r['phase_deg']:.2f}" if r["phase_deg"] is not None else "",
                          *[f"{a:.4f}" for a in r["A"]],
                          f"{r['split_rms']:.3e}", r["BF"], r["pa"],
                          f"{r['s_circ']:.3e}" if r["s_circ"] is not None else "",
                          f"{r['charge_raw']:.6f}" if r["charge_raw"] is not None else ""])
    (HERE / "census_table_v2.json").write_text(
        json.dumps({"rows": [{k: v for k, v in r.items()} for r in rows],
                    "imports": {"abl": sha256(ABL), "mph": sha256(MPH)},
                    "params": {"T_RELAX": T_RELAX, "WIN": WIN, "T_HALF": T_HALF,
                                "A_BASE": A_BASE}},
                   ensure_ascii=False, indent=2, default=float), encoding="utf-8")

    # Markdown 表（単一ソースからの再現出力）
    md = ["# 粒子的グループ census v2（N=5 / N=40）", "",
          f"生成プログラム: `run_paper8_particle_census_v2.py`"
          f"（seed 40260801/40260802・緩和 {T_RELAX} tick・決定論的に再現可能）", "",
          "- ベース波長: 偏差重み付き巻き数整列度 A_k（k≤M/2）の最大が 0.3 超のとき 1/k。"
          "「—」は一様海（構造なし）",
          "- 倍音1..6: A_1..A_6（構造スペクトル、0..1）",
          "- B/F: 行動的半周期テスト（F=奇=フェルミオン型 / B=偶=ボゾン型）",
          "- 電荷(生): 塊周波数/残余周波数の生の比（緩和後のみ。一定にならなくても生値のまま記録）",
          "- 粒/反: エンジンの時計方向（全状態共通の基準）に対する周回符号 s の規約ラベル。"
          "粒=s>0/反=s<0。絶対名は規約であり相対符号のみ物理。一様海は「—」（自己共役＝区分なし、光子と同型）",
          "",
          "| N | 状態 | 族 | 時点 | ベース波長 | 位相° | 倍音1 | 倍音2 | 倍音3 | 倍音4 | 倍音5 | 倍音6 | B/F | 粒/反 | 電荷(生) |",
          "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in rows:
        wl = f"1/{r['k0']}" if r["k0"] else "—"
        ph = f"{r['phase_deg']:+.1f}" if r["phase_deg"] is not None else "—"
        ch = f"{r['charge_raw']:.6f}" if r["charge_raw"] is not None else "—"
        md.append("| " + " | ".join(
            [str(r["N"]), r["state"], r["family"], r["time"], wl, ph]
            + [f"{a:.2f}" for a in r["A"]] + [r["BF"], r["pa"], ch]) + " |")
    (HERE / "census_table_v2.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    print(f"\nsaved: census_table_v2.csv / census_table_v2.json / census_table_v2.md "
          f"({time.time()-t0:.0f}s)")


if __name__ == "__main__":
    main()
