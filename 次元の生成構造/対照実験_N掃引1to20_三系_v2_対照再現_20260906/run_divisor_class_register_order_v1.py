#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""約数類定理のレジスタ位数追随検定 v1

事前登録: `事前登録_約数類定理レジスタ位数追随検定_v1.md`（本走行前に固定済み）

問い:
  前論文『波の周期表』が ne=16 で実測した約数類定理——「海中の種の性質は
  巻き数 m そのものではなく gcd(m, ne) にしか依存しない」——は、ne を
  変えたとき gcd(m, ne) に追随するか。前論文自身が §13-3 に
  「追随すれば約数類定理が一般化され、しなければ 16 固有の隠れ構造がある」
  と検証可能予言として登録している。

手続き（シリーズ内完結の再現性規約）:
  測定部は前論文の実行器 `波の周期表検討/run_pre_v10b_longcoupling_v1.py` の
  main() から**逐語コピー**した（single_winding / project_eta / shift_eta /
  band_power / comp_gram_band / 状態構成 / 窓測定）。独自の再実装はしない。
  変更したのは次の3点のみ:
    (1) eta_grid_n を引数化した（前論文は既定の 16 固定）
    (2) m_t の上限を ne/2 - 1 とした（ne=16 では 1..7 で前論文と同一）
    (3) 位数ごとの gcd 類判定と、ne=16 の保存済み結果との照合を追加した

必須対照（fail-closed）:
  ne=16 の全数値が `波の周期表検討/pre_v10b_longcoupling_result_v1.json` と
  相対 1e-12 以内で一致しなければ、以降の判定を行わず記録のみで終了する。

使い方: python3 run_divisor_class_register_order_v1.py
"""
from __future__ import annotations
import hashlib
import importlib.util
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
UIM = ROOT / "万能非弾性写像_managed_v1"
PT = ROOT / "波の周期表検討"
REF_JSON = PT / "pre_v10b_longcoupling_result_v1.json"
SRC_PRIOR = PT / "run_pre_v10b_longcoupling_v1.py"

spec = importlib.util.spec_from_file_location(
    "exact_divclass", UIM / "run_ignition_fate_exact_v3.py")
ex = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = ex
spec.loader.exec_module(ex)
v1, toy, base = ex.v1, ex.toy, ex.base

# ---- 前論文と同一の宣言値（逐語） ----
S = 8.0
J_WIN = 40
SETTLES = (2000, 4000)

# ---- 本検定の宣言値（事前登録 §3・§4・§6） ----
NE_LIST = (16, 8, 12, 32)
TOL_SAME = 1e-5      # 「一致する」の相対閾値
TOL_DIFF = 1e-3      # 「異なる」の相対閾値
TOL_REPRO = 1e-12    # ne=16 再現の相対閾値


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def run_one_ne(ne: int) -> dict:
    """前論文 main() の逐語コピー（eta_grid_n と m_t 上限のみ引数化）。"""
    params = base.Params(high_n=63, recursive_collision_count=200, eta_grid_n=ne)
    sp = base.build_source_params(params)
    n, ne_ = sp.chi_grid_n, sp.eta_grid_n
    assert ne_ == ne, f"eta_grid_n の設定が効いていない: {ne_} != {ne}"
    shape = (n, ne_)
    ms = np.arange(ne_)
    mm = np.where(ms <= ne_ // 2, ms, ms - ne_)
    eta = 2 * np.pi * np.arange(ne_) / ne_

    def single_winding(v):
        f = np.fft.fft(v.reshape(shape), axis=0, norm="ortho")
        f[0, :] = 0.0
        f[n // 2:, :] = 0.0
        return np.fft.ifft(f, axis=0, norm="ortho").reshape(v.shape)

    def project_eta(v, m_set):
        f = np.fft.fft(v.reshape(shape), axis=1, norm="ortho")
        keep = np.isin(mm, list(m_set))
        f[:, ~keep] = 0.0
        return np.fft.ifft(f, axis=1, norm="ortho").reshape(v.shape)

    def shift_eta(v, dm):
        return (v.reshape(shape) * np.exp(1j * dm * eta)[None, :]).reshape(v.shape)

    a0 = single_winding(v1.make_bundle(sp, (30, 32, 34), "A", scale=1.0)) * S
    b0 = single_winding(v1.make_bundle(sp, (30, 32, 34), "B", scale=1.0)) * S
    a0 = a0 + single_winding(v1.make_bundle(sp, (21,), "A", scale=1.0)) * (0.2 * S)
    pow0 = float(np.sum(np.abs(a0) ** 2) + np.sum(np.abs(b0) ** 2))
    a1 = project_eta(a0, {1})
    b1 = project_eta(b0, {1})
    pw = float(np.sum(np.abs(a1) ** 2) + np.sum(np.abs(b1) ** 2))
    sc = np.sqrt(pow0 / pw)
    a1 *= sc
    b1 *= sc
    sea_a = project_eta(v1.make_bundle(sp, (29, 31, 33), "A", scale=1.0) * S, {0})
    sea_b = project_eta(v1.make_bundle(sp, (29, 31, 33), "B", scale=1.0) * S, {0})
    pws = float(np.sum(np.abs(sea_a) ** 2) + np.sum(np.abs(sea_b) ** 2))
    scs = np.sqrt(0.25 * pow0 / pws)
    sea_a *= scs
    sea_b *= scs

    def band_power(a, b, m_t):
        fa = np.fft.fft(np.fft.fft(a.reshape(shape), axis=0, norm="ortho"),
                        axis=1, norm="ortho")
        fb = np.fft.fft(np.fft.fft(b.reshape(shape), axis=0, norm="ortho"),
                        axis=1, norm="ortho")
        P = (np.abs(fa) ** 2 + np.abs(fb) ** 2) / 2
        return float(P[:, mm == m_t].sum())

    def comp_gram_band(A, B, band):
        idx = np.argwhere(band)
        P = np.mean(np.abs(A) ** 2 + np.abs(B) ** 2, axis=0) / 2
        m2s, sms, szs, ws = [], [], [], []
        for (ki, mi) in idx:
            At = A[:, ki, mi]
            Bt = B[:, ki, mi]
            wA = np.angle(np.sum(At[1:] * np.conj(At[:-1])))
            wB = np.angle(np.sum(Bt[1:] * np.conj(Bt[:-1])))
            t_ = np.arange(A.shape[0])
            Ad = At * np.exp(-1j * wA * t_)
            Bd = Bt * np.exp(-1j * wB * t_)
            Gaa = np.mean(np.abs(Ad) ** 2)
            Gbb = np.mean(np.abs(Bd) ** 2)
            Gab = np.mean(Ad * np.conj(Bd))
            T = 0.5 * (Gaa + Gbb)
            if T <= 0:
                continue
            det = Gaa * Gbb - abs(Gab) ** 2
            X = Gab.real
            Y = -Gab.imag
            Z = 0.5 * (Gaa - Gbb)
            m2s.append(det / T ** 2)
            sms.append(np.sqrt(X ** 2 + Y ** 2 + Z ** 2) / T)
            szs.append(Z / T)
            ws.append(P[ki, mi])
        ws = np.array(ws) / max(sum(ws), 1e-300)
        return (float(np.sum(ws * np.array(m2s))),
                float(np.sum(ws * np.array(sms))),
                float(np.sum(ws * np.array(szs))))

    m_max = ne_ // 2 - 1                       # 事前登録 §3: Nyquist 帯は除く
    out = {"ne": ne_, "chi_grid_n": n, "J_WIN": J_WIN,
           "SETTLES": list(SETTLES), "m_range": [1, m_max], "rows": []}
    for js in SETTLES:
        print(f"  ==== ne={ne_} settle={js} ====")
        print(f"  {'m':>3} {'gcd':>3} {'質量²(補償)':>14} {'S':>10} "
              f"{'s_z':>10} {'保持率':>10}")
        for m_t in range(1, m_max + 1):
            a = shift_eta(a1, m_t - 1) + sea_a
            b = shift_eta(b1, m_t - 1) + sea_b
            p_init = band_power(a, b, m_t)
            for _ in range(js):
                a, b, _ = ex.collision_step_exact(a, b, sp)
            A = np.zeros((J_WIN, n, ne_), complex)
            B = np.zeros((J_WIN, n, ne_), complex)
            for t in range(J_WIN):
                a, b, _ = ex.collision_step_exact(a, b, sp)
                fa = np.fft.fft(a.reshape(shape), axis=0, norm="ortho")
                fb = np.fft.fft(b.reshape(shape), axis=0, norm="ortho")
                A[t] = np.fft.fft(fa, axis=1, norm="ortho")
                B[t] = np.fft.fft(fb, axis=1, norm="ortho")
            P = np.mean(np.abs(A) ** 2 + np.abs(B) ** 2, axis=0) / 2
            band = (mm[None, :] == m_t) & (P > P.max() * 1e-8)
            if band.sum() == 0:
                print(f"  {m_t:>+3}  帯空")
                continue
            m2, sm, sz = comp_gram_band(A, B, band)
            ret = float(P[:, mm == m_t].sum() / max(p_init, 1e-300))
            g = math.gcd(m_t, ne_)
            print(f"  {m_t:>+3} {g:>3} {m2:>14.9f} {sm:>10.6f} "
                  f"{sz:>+10.6f} {ret:>10.7f}")
            out["rows"].append({"settle": js, "m": m_t, "gcd": g, "mass2": m2,
                                "S": sm, "sz": sz, "retention": ret})
        a, b = sea_a.copy(), sea_b.copy()
        for _ in range(js):
            a, b, _ = ex.collision_step_exact(a, b, sp)
        A = np.zeros((J_WIN, n, ne_), complex)
        B = np.zeros((J_WIN, n, ne_), complex)
        for t in range(J_WIN):
            a, b, _ = ex.collision_step_exact(a, b, sp)
            fa = np.fft.fft(a.reshape(shape), axis=0, norm="ortho")
            fb = np.fft.fft(b.reshape(shape), axis=0, norm="ortho")
            A[t] = np.fft.fft(fa, axis=1, norm="ortho")
            B[t] = np.fft.fft(fb, axis=1, norm="ortho")
        P = np.mean(np.abs(A) ** 2 + np.abs(B) ** 2, axis=0) / 2
        band = P > P.max() * 1e-8
        m2, sm, sz = comp_gram_band(A, B, band)
        print(f"  海単独: 質量²={m2:.9f} S={sm:.6f}")
        out["rows"].append({"settle": js, "m": 0, "gcd": None, "mass2": m2,
                            "S": sm, "sz": sz, "retention": None})
    return out


def check_reproduction(rows16: list) -> dict:
    """必須対照: 保存済み ne=16 結果との照合（fail-closed）。"""
    ref = json.loads(REF_JSON.read_text())
    idx = {(r["settle"], r["m"]): r for r in ref["rows"]}
    worst, worst_key, n_cmp = 0.0, None, 0
    for r in rows16:
        key = (r["settle"], r["m"])
        if key not in idx:
            return {"pass": False, "reason": f"参照に {key} が無い"}
        q = idx[key]
        for f in ("mass2", "S", "sz", "retention"):
            if q.get(f) is None or r.get(f) is None:
                continue
            den = max(abs(q[f]), 1e-300)
            rel = abs(r[f] - q[f]) / den
            n_cmp += 1
            if rel > worst:
                worst, worst_key = rel, (key, f)
    return {"pass": bool(worst <= TOL_REPRO), "worst_rel": worst,
            "worst_at": worst_key, "n_compared": n_cmp,
            "ref_sha256": sha256(REF_JSON), "tol": TOL_REPRO}


def judge(res: dict) -> dict:
    """gcd 類と実測の一致パターンが合っているかを判定する（事前登録 §4・§5）。"""
    ne = res["ne"]
    verdicts = []
    for js in SETTLES:
        rows = [r for r in res["rows"] if r["settle"] == js and r["m"] >= 1]
        for i, ra in enumerate(rows):
            for rb in rows[i + 1:]:
                rels = []
                for f in ("mass2", "S", "retention"):
                    den = max(abs(rb[f]), 1e-300)
                    rels.append(abs(ra[f] - rb[f]) / den)
                mx = max(rels)
                same_meas = mx <= TOL_SAME
                diff_meas = mx >= TOL_DIFF
                same_pred = (ra["gcd"] == rb["gcd"])
                if same_pred:
                    ok = same_meas
                else:
                    ok = diff_meas
                verdicts.append({"settle": js, "pair": [ra["m"], rb["m"]],
                                 "gcd": [ra["gcd"], rb["gcd"]],
                                 "same_predicted": same_pred,
                                 "max_rel": mx, "ok": bool(ok)})
    return {"ne": ne, "n_pairs": len(verdicts),
            "n_fail": sum(1 for v in verdicts if not v["ok"]),
            "pass": all(v["ok"] for v in verdicts), "pairs": verdicts}


def main() -> None:
    t0 = time.time()
    out = {
        "prereg": "事前登録_約数類定理レジスタ位数追随検定_v1.md",
        "source_sha256": {
            "self": sha256(Path(__file__).resolve()),
            "prior_runner": sha256(SRC_PRIOR),
            "engine": sha256(UIM / "run_ignition_fate_exact_v3.py"),
        },
        "declared": {"S": S, "J_WIN": J_WIN, "SETTLES": list(SETTLES),
                     "NE_LIST": list(NE_LIST), "TOL_SAME": TOL_SAME,
                     "TOL_DIFF": TOL_DIFF, "TOL_REPRO": TOL_REPRO},
        "runs": {}, "judgements": {},
    }

    print("=== P1: ne=16 の再現（必須対照・fail-closed）===")
    r16 = run_one_ne(16)
    out["runs"]["16"] = r16
    rep = check_reproduction(r16["rows"])
    out["reproduction_ne16"] = rep
    print(f"\n  再現照合: 比較 {rep['n_compared']} 件・最大相対差 "
          f"{rep['worst_rel']:.3e}（許容 {TOL_REPRO:.0e}）→ "
          f"{'PASS' if rep['pass'] else 'FAIL'}")
    if not rep["pass"]:
        print(f"  最悪箇所: {rep['worst_at']}")
        print("  → 環境不一致。事前登録 §5 により以降の判定を行わない。")
        out["status"] = "reproduction_failed"
        (HERE / "result_divisor_class_register_order_v1.json").write_text(
            json.dumps(out, indent=1, ensure_ascii=False))
        return

    out["judgements"]["16"] = judge(r16)
    for ne in NE_LIST:
        if ne == 16:
            continue
        print(f"\n=== ne={ne} ===")
        r = run_one_ne(ne)
        out["runs"][str(ne)] = r
        out["judgements"][str(ne)] = judge(r)

    print("\n=== 判定（事前登録 §4・§5）===")
    print(f"{'ne':>4} {'対の数':>6} {'不一致':>6} {'判定':>6}")
    for ne in NE_LIST:
        j = out["judgements"][str(ne)]
        print(f"{ne:>4} {j['n_pairs']:>6} {j['n_fail']:>6} "
              f"{'PASS' if j['pass'] else 'FAIL':>6}")

    # P3 の要点: m=1 と m=3 の関係が位数で反転するか
    focus = {}
    for ne in NE_LIST:
        for v in out["judgements"][str(ne)]["pairs"]:
            if v["settle"] == 4000 and sorted(v["pair"]) == [1, 3]:
                focus[str(ne)] = v
    out["focus_m1_vs_m3"] = focus
    print("\n--- P3 の要点: m=1 と m=3（settle=4000）---")
    print(f"{'ne':>4} {'gcd(1,ne)':>10} {'gcd(3,ne)':>10} {'予言':>8} "
          f"{'実測最大相対差':>14} {'判定':>6}")
    for ne in NE_LIST:
        v = focus.get(str(ne))
        if v is None:
            print(f"{ne:>4} {'—':>10} {'—':>10} {'—':>8} "
                  f"{'（m=3 が範囲外）':>14} {'—':>6}")
            continue
        print(f"{ne:>4} {v['gcd'][0]:>10} {v['gcd'][1]:>10} "
              f"{'一致' if v['same_predicted'] else '不一致':>8} "
              f"{v['max_rel']:>14.3e} {'OK' if v['ok'] else 'NG':>6}")

    out["status"] = "complete"
    out["all_pass"] = all(out["judgements"][str(ne)]["pass"] for ne in NE_LIST)
    out["runtime_sec"] = time.time() - t0
    (HERE / "result_divisor_class_register_order_v1.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False))
    print(f"\n総合: {'約数類定理は ne に追随する' if out['all_pass'] else '追随しない箇所あり'}")
    print(f"完了 {out['runtime_sec']:.0f}s → result_divisor_class_register_order_v1.json")


if __name__ == "__main__":
    main()
