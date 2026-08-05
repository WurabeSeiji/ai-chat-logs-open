#!/usr/bin/env python3
"""GATE-0g: E8帯（N=238..248）の特殊値チェック——線形相の数論盲目性の外挿検定

目的（2026-08-05 木原氏指示「念のため特異点Nのチェック」）:
    E8関連特殊値の大N選択性の警戒（考察md §4.6）。GATE-0c地図（N≤144）では
    31/120/124/128/137/144 は全てnull（局所トレンド比0.93〜1.04倍）。
    本実験は地図外の 240（E8根の数）/244（指示値）/248（dim E8）と
    近傍対照 238/242/246 を同じ計器で測る。

予言（実行前固定・事後変更禁止）:
    線形相は数論盲目（GATE-0b/0c）→ 全6点の観測量（セクター末尾r・O・ρ、
    親のclass・crossing）は滑らかなトレンドの内挿・外挿上に載り、
    特殊値と対照の差は近傍散らばり以内（滑り偏差2倍規則に非該当）。
    偏差が出た場合はそれ自体が発見（線形相にE8算術が漏れている証拠）。

方法: GATE-0c と同一の計器・規約（T=4000、末尾1000、セクター標本12本、
    生成器シード2..9フォールバック、wp規約同一）。N単位で6プロセス並列。

使い方: python3 run_stage0g_e8_band_check_v1.py
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import math
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
SERIES_DIR = HERE.parent
ABL = SERIES_DIR / "第8論文_二段階seed除去による準安定相の因果分離" / "code" / "run_preliminary_seed_ablation_v1.py"
GEN3 = SERIES_DIR / "make_parent_white_managed_v1" / "make_parent_white_harmonics_n_only_v3.py"
NS = [238, 240, 242, 244, 246, 248]
SPECIAL = {240, 244, 248}
T_LONG = 4000
DELTA_O = 10
TAIL = 1000
OMEGA_ENGINE = math.pi / 72.0
WORKERS = 6

_mods = {}


def _load():
    if _mods:
        return _mods["abl"], _mods["gen3"]
    spec = importlib.util.spec_from_file_location("abl_g0g", ABL)
    abl = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = abl
    spec.loader.exec_module(abl)
    spec2 = importlib.util.spec_from_file_location("gen3_g0g", GEN3)
    gen3 = importlib.util.module_from_spec(spec2)
    sys.modules[spec2.name] = gen3
    spec2.loader.exec_module(gen3)
    _mods["abl"] = abl
    _mods["gen3"] = gen3
    return abl, gen3


def orth_plane(Z):
    p = Z.real / np.linalg.norm(Z.real)
    q = Z.imag - (Z.imag @ p) * p
    nq = np.linalg.norm(q)
    if nq < 1e-300:
        return np.column_stack([p, np.zeros_like(p)])
    return np.column_stack([p, q / nq])


def run_one(abl, n, Z0, wp, track_f):
    sys_lr = abl.LowRankSystem(n)
    sys_lr.set_theta(np.angle(Z0))
    Z = Z0 / np.linalg.norm(Z0)
    if track_f:
        p0 = Z.real / np.linalg.norm(Z.real)
        q0 = Z.imag - (Z.imag @ p0) * p0
        q0 = q0 / np.linalg.norm(q0)
    rs = np.zeros(T_LONG)
    ws = np.zeros(T_LONG)
    fs = np.zeros(T_LONG + 1) if track_f else None
    Os = []
    B_prev = orth_plane(Z)
    for t in range(1, T_LONG + 1):
        Znew, wp = abl.evolve(sys_lr, Z, wp)
        ip = np.conj(Z) @ Znew
        ph = ip / abs(ip) if abs(ip) > 0 else 1.0
        rs[t - 1] = float(np.linalg.norm(Znew - ph * Z))
        ws[t - 1] = float(np.angle(ip))
        Z = Znew
        if track_f:
            Zp = Z - p0 * (p0 @ Z) - q0 * (q0 @ Z)
            fs[t] = float(np.real(np.conj(Zp) @ Zp))
        if t % DELTA_O == 0:
            B = orth_plane(Z)
            Os.append(float(np.sum((B_prev.T @ B) ** 2) / 2.0))
            B_prev = B
    Os = np.array(Os)
    out = {"r_tail": float(np.median(rs[-TAIL:])),
           "O_tail": float(np.median(Os[-(TAIL // DELTA_O):])),
           "rho_tail": float(np.median(ws[-TAIL:]) / OMEGA_ENGINE)}
    if track_f:
        crossing = next((t for t, f in enumerate(fs) if f > 0.05), None)
        below = np.nonzero(fs < 1e-20)[0]
        out["crossing"] = crossing
        out["t_launch"] = int(below.max()) if below.size else 0
    return out


def sector_list(n):
    allowed = [k for k in range(1, n) if (2 * k) % n != 0]
    cand = [1, 2, 3, n - 1, n - 2, n // 6, n // 4, n // 3,
            (2 * n) // 5, n // 2 - 1, (3 * n) // 5, (2 * n) // 3]
    out = []
    for k in cand:
        if k in allowed and k not in out:
            out.append(k)
    return out


def process_n(n):
    abl, gen3 = _load()
    t0 = time.time()
    m = n * (n - 1) // 2
    parent = None
    gen_seed = None
    for seed in range(2, 10):
        try:
            parent = gen3.make_parent(n, seed=seed)
            gen_seed = seed
            break
        except Exception:
            continue
    if parent is None:
        return n, {"error": "make_parent failed seeds 2..9"}
    C = np.fft.fft(parent.relation_waves, axis=1) / n
    vp = parent.parent_vector / np.linalg.norm(parent.parent_vector)
    pres = run_one(abl, n, vp, np.random.default_rng(91000).normal(size=m), True)
    if pres["crossing"] is not None:
        pclass = "equilibrium-burst"
    else:
        pclass = "no-crossing(T=4000)"
    secs = {}
    for k in sector_list(n):
        Zk = C[:, k] / np.linalg.norm(C[:, k])
        wp = np.random.default_rng(92000 + k).normal(size=m)
        secs[str(k)] = run_one(abl, n, Zk, wp, False)
    med_r = float(np.median([s["r_tail"] for s in secs.values()]))
    med_O = float(np.median([s["O_tail"] for s in secs.values()]))
    med_rho = float(np.median([s["rho_tail"] for s in secs.values()]))
    return n, {"gen_seed": gen_seed, "n_sectors": len(secs),
                "parent": {**pres, "class": pclass},
                "sector_median_r": med_r, "sector_median_O": med_O,
                "sector_median_rho": med_rho,
                "sectors": secs, "runtime_sec": time.time() - t0}


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main() -> None:
    t0 = time.time()
    print(f"GATE-0g E8帯チェック N={NS}  ABL {sha256(ABL)[:16]}…  GEN3 {sha256(GEN3)[:16]}…")
    results = {}
    with ProcessPoolExecutor(max_workers=WORKERS) as ex:
        for n, res in ex.map(process_n, NS):
            results[str(n)] = res
            if "error" in res:
                print(f"  N={n}: {res['error']}", flush=True)
            else:
                tag = "★特殊" if n in SPECIAL else "対照"
                print(f"  N={n} [{tag}]: 親={res['parent']['class']} "
                      f"crossing={res['parent'].get('crossing')} "
                      f"sec r={res['sector_median_r']:.3e} O={res['sector_median_O']:.6f} "
                      f"ρ={res['sector_median_rho']:.4f} ({res['runtime_sec']:.0f}s)", flush=True)

    ok = {int(n): r for n, r in results.items() if "error" not in r}
    if len(ok) >= 4:
        ctrl = [n for n in ok if n not in SPECIAL]
        spec_ns = [n for n in ok if n in SPECIAL]
        cr = np.median([ok[n]["sector_median_r"] for n in ctrl])
        devs = {n: float(ok[n]["sector_median_r"] / cr) for n in spec_ns}
        print(f"\n対照中央値 r={cr:.3e}  特殊値の偏差: "
              + ", ".join(f"N={n}:×{d:.3f}" for n, d in devs.items()))
        nullres = all(0.5 < d < 2.0 for d in devs.values())
        print(f"判定（2倍規則）: {'null（数論盲目の外挿成立）' if nullres else '偏差あり——要精査'}")
    out = {"NS": NS, "SPECIAL": sorted(SPECIAL),
           "imports": {"abl": sha256(ABL), "generator_v3": sha256(GEN3)},
           "criteria": {"T_LONG": T_LONG, "TAIL": TAIL,
                         "prediction": "all six N on smooth trend; special/control diff within 2x"},
           "results": results, "runtime_sec": time.time() - t0}
    (HERE / "gate0g_e8_band_check_v1.json").write_text(
        json.dumps(out, ensure_ascii=False, indent=2, default=float), encoding="utf-8")
    print(f"saved ({out['runtime_sec']:.0f}s)")


if __name__ == "__main__":
    main()
