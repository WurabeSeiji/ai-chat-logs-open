#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""TB-2/TB-3: 標準物質宇宙のN系列——常時実行G・粒子census・固有時間の対照実験

テストベッド論文の中核実験（木原仕様）:
  各Nで標準物質宇宙（ポンプ凝縮体＋中性シード）を通し走行し、
  (1) 常時実行 G パネル（第0步から毎ステップ・一様ケイデンス・切替なし）
  (2) 粒子census——毛（電荷）×帯（k）の帳簿: 何種が・何個・どんな分布か
  (3) シードなし対照＝真空宇宙: 空間はできるが固有時間が取得できないこと
  を測る。N ∈ {4, 6, 12}（4,6=特異点対照・12=想定テストベッド。20-100は
  同一スクリプトの後続走行に登録）。

正本: HairEngine（stage3毛込みcensus・V2でコピー移植済みのものを再利用）。
宇宙 = ポンプ(k=2,毛0)=control親 ＋ シード(k=1,毛0)×δ=1e-2（V2中性レシピ）。

判定（事前固定）:
 (TB2a) 常時パネルの全時代定義性: 全ステップで全読出しが有限値（NaN/Inf なし）。
 (TB2b) 受動性: パネル併走の有無で終状態がビット単位同一（N=12で検証）。
 (TB2c) 粒子censusの決定性: 生成物質は帳簿上、予言帯（k=3×毛0 と k=1×毛0）
        のみに集中（帯外比 <1e-6）——「何が生まれるか」が完全に指定される。
 (TB2d) 固有時間の対照: シード入りでは物質時計が取得可能（生成内容の位相
        前進レート ω̂ が安定に読める: 窓間ばらつき/|ω̂| < 0.1）。
        シードなしでは物質内容が測定限界以下（クロック担体不在）で
        ω̂ が定義不能——「真空宇宙では空間はあるが固有時間は取得できない」。
使い方: python3 run_tb2_matter_universe_series_v1.py
"""
from __future__ import annotations
import importlib.util, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
MB = HERE.parent / "万能相互作用多体接続_v1"

spec3 = importlib.util.spec_from_file_location(
    "s3tb", MB / "run_stage3_sharedO_v2_and_hair_v1.py")
s3 = importlib.util.module_from_spec(spec3)
sys.modules[spec3.name] = s3
spec3.loader.exec_module(s3)
abl = s3.abl
s2 = s3.s2

spec_g3 = importlib.util.spec_from_file_location(
    "g3tb", MB / "run_genesis_v3_register_local_v1.py")
g3 = importlib.util.module_from_spec(spec_g3)
sys.modules[spec_g3.name] = g3
spec_g3.loader.exec_module(g3)

NN, NETA = 5, 8
T_RUN = 600
DELTA = 1e-2
N_LIST = (4, 6, 12)


class HairEngine(s3.VertexEngineV2):
    """正本 stage3 の HairEngine コピー移植（V2 と同一）"""

    def __init__(self, n_, C2_0, wp, Nn=NN, Neta=NETA, **kw):
        C0f = C2_0.reshape(C2_0.shape[0], -1)
        super().__init__(n_, C0f, wp, **kw)
        self.Nn, self.Neta = Nn, Neta
        ks = np.arange(Nn)
        self.odd_k = (ks % 2 == 1)
        self.even_k = (ks % 2 == 0) & (ks != 0)

    def C2(self):
        return self.C.reshape(self.m, self.Nn, self.Neta)

    def _readout(self):
        P2 = np.abs(self.C2()) ** 2
        Pk = P2.sum(axis=2)
        Av = np.zeros((self.n, self.Nn))
        np.add.at(Av, self.ia, Pk)
        np.add.at(Av, self.ib, Pk)
        Sagg = Av[self.ia] + Av[self.ib] - 2 * Pk
        comb = Pk + Sagg
        Pf = comb[:, self.odd_k].sum(axis=1)
        Pb = comb[:, self.even_k].sum(axis=1)
        th = np.arctan2(np.sqrt(np.maximum(Pf, 0)), np.sqrt(np.maximum(Pb, 0)))
        return self.scale * np.sin(th) ** 2

    def _nonlinear(self):
        R = self._readout()
        if not np.any(R > 0):
            return
        C2 = self.C2()
        W = np.fft.ifft2(C2, axes=(1, 2)) * (self.Nn * self.Neta)
        Wf = W.reshape(self.m, -1)
        rate0 = self._vertex_rate(Wf, R)
        Lmax = float(np.max(np.abs(rate0))) / max(float(np.max(np.abs(Wf))), 1e-300)
        nsub = max(1, int(np.ceil(Lmax / s2.H_MAX)))
        h = 1.0 / nsub
        for _ in range(nsub):
            k1 = self._vertex_rate(Wf, R)
            k2 = self._vertex_rate(Wf + 0.5 * h * k1, R)
            k3 = self._vertex_rate(Wf + 0.5 * h * k2, R)
            k4 = self._vertex_rate(Wf + h * k3, R)
            Wf = Wf + (h / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
        W = Wf.reshape(self.m, self.Nn, self.Neta)
        self.C = (np.fft.fft2(W, axes=(1, 2)) / (self.Nn * self.Neta)
                  ).reshape(self.m, -1)


def g_panel(eng):
    """場の万能読出しパネル（常時実行・切替なし・パラメータなし）"""
    C2 = eng.C2()
    P2 = np.abs(C2) ** 2                     # M×Nn×Nη
    P_tot = float(P2.sum())
    Pk = P2.sum(axis=(0, 2))                 # 帯別
    Pkm = P2.sum(axis=0)                     # 帯×毛の帳簿（census）
    f_odd = float(Pk[1] + Pk[3]) / max(P_tot, 1e-300)
    # 位置読出し: 奇数帯内容の双対レジスタ重心（巻きモーメント）
    W = np.fft.ifft2(C2, axes=(1, 2)) * (NN * NETA)
    mask = np.zeros(NN); mask[1] = 1.0; mask[3] = 1.0
    Codd = C2 * mask[None, :, None]
    Wo = np.fft.ifft2(Codd, axes=(1, 2)) * (NN * NETA)
    Pn = np.sum(np.abs(Wo) ** 2, axis=(0, 2))
    if Pn.sum() > 1e-300:
        z = np.sum(Pn * np.exp(2j * np.pi * np.arange(NN) / NN)) / Pn.sum()
        x_pos = float((np.angle(z) * NN / (2 * np.pi)) % NN)
        pr = float(Pn.sum() ** 2 / max(np.sum(Pn ** 2), 1e-300))
    else:
        x_pos, pr = float("nan"), float("nan")
    # 物質時計担体: 生成帯 k=3 の総振幅
    c3 = C2[:, 3, :].reshape(-1)
    return {"P_tot": P_tot, "f_odd": f_odd, "Pk": [float(v) for v in Pk],
            "census": Pkm, "x_pos": x_pos, "PR_n": pr, "c3": c3}


def run_universe(n, seeded, panel_on=True):
    _, v, _, _, _, _, _, Z0c, wp0 = abl.build_init(n, False)
    Z0c = Z0c / np.linalg.norm(Z0c)
    m = n * (n - 1) // 2
    seed_state = g3.zero_closure_state(m, np.random.default_rng(98000))
    C2_0 = np.zeros((m, NN, NETA), complex)
    C2_0[:, 2, 0] = Z0c
    if seeded:
        C2_0[:, 1, 0] = DELTA * seed_state
    eng = HairEngine(n, C2_0, wp0)
    hist = {"f": [], "x": [], "pr": [], "omega": [], "census_final": None,
            "nonfinite": 0}
    prev_c3 = None
    for t in range(1, T_RUN + 1):
        eng.step()
        if panel_on:
            p = g_panel(eng)
            vals = [p["P_tot"], p["f_odd"]] + p["Pk"]
            if not all(np.isfinite(vv) for vv in vals):
                hist["nonfinite"] += 1
            hist["f"].append(p["f_odd"])
            hist["x"].append(p["x_pos"])
            hist["pr"].append(p["PR_n"])
            if prev_c3 is not None:
                z = np.vdot(prev_c3, p["c3"])
                amp = float(np.linalg.norm(p["c3"]))
                hist["omega"].append(
                    float(np.angle(z)) if (abs(z) > 1e-30 and amp > 1e-12)
                    else float("nan"))
            prev_c3 = p["c3"].copy()
            if t == T_RUN:
                hist["census_final"] = p["census"]
    return eng, hist


def clock_metrics(hist):
    om = np.array(hist["omega"][T_RUN // 2:])
    valid = om[np.isfinite(om)]
    if len(valid) < 10:
        return None, None, 0.0
    return (float(np.mean(valid)), float(np.std(valid)),
            float(len(valid) / len(om)))


def main():
    t0 = time.time()
    out = {"params": {"Nn": NN, "Neta": NETA, "T": T_RUN, "delta": DELTA,
                      "N_list": list(N_LIST)}, "N": {}}
    for n in N_LIST:
        eng_s, h_s = run_universe(n, seeded=True)
        eng_v, h_v = run_universe(n, seeded=False)
        om_mean, om_std, om_valid = clock_metrics(h_s)
        _, _, om_valid_vac = clock_metrics(h_v)
        cen = h_s["census_final"]              # Nn×Nη
        tot = float(cen.sum())
        allowed = float(cen[1, 0] + cen[3, 0] + cen[2, 0])
        outside = (tot - allowed) / max(tot - float(cen[2, 0]), 1e-300)
        # 種の帳簿（ポンプ以外・パワー降順上位）
        book = []
        for k in range(NN):
            for mm in range(NETA):
                if (k, mm) != (2, 0) and cen[k, mm] > 1e-300:
                    book.append((float(cen[k, mm]), k, mm))
        book.sort(reverse=True)
        clock_ok = (om_mean is not None and om_valid > 0.9
                    and abs(om_std / om_mean) < 0.1) if om_mean else False
        out["N"][n] = {
            "f_final_seeded": h_s["f"][-1], "f_final_vacuum": h_v["f"][-1],
            "nonfinite_panel": h_s["nonfinite"] + h_v["nonfinite"],
            "census_top": [{"P": b[0], "k": b[1], "hair": b[2]}
                           for b in book[:6]],
            "outside_pred_ratio": float(outside),
            "clock_seeded": {"omega": om_mean, "std": om_std,
                             "valid_frac": om_valid, "ok": bool(clock_ok)},
            "clock_vacuum_valid_frac": om_valid_vac,
            "x_final": h_s["x"][-1], "PR_final": h_s["pr"][-1],
            "f_series_seeded": h_s["f"][::10],
            "f_series_vacuum": h_v["f"][::10],
        }
        print(f"N={n}: f(種)={h_s['f'][-1]:.3e} f(真空)={h_v['f'][-1]:.3e} "
              f"帳簿外比={outside:.2e} "
              f"時計(種)={'ω=%.4f±%.4f 有効%.2f' % (om_mean, om_std, om_valid) if om_mean else '取得不能'} "
              f"時計(真空)有効率={om_valid_vac:.2f} x={h_s['x'][-1]:.2f}")
        print(f"   census上位: " + "  ".join(
            f"(k={b['k']},毛={b['hair']}):{b['P']:.2e}"
            for b in out["N"][n]["census_top"][:4]))

    # TB2b 受動性（N=12・パネル有無で終状態一致）
    eng_a, _ = run_universe(12, seeded=True, panel_on=True)
    eng_b, _ = run_universe(12, seeded=True, panel_on=False)
    passive = bool(np.array_equal(eng_a.C, eng_b.C))

    tb2a = all(out["N"][n]["nonfinite_panel"] == 0 for n in N_LIST)
    tb2c = all(out["N"][n]["outside_pred_ratio"] < 1e-6 for n in N_LIST)
    tb2d = all(out["N"][n]["clock_seeded"]["ok"] for n in N_LIST) and \
        all(out["N"][n]["clock_vacuum_valid_frac"] < 0.5 for n in N_LIST) and \
        all(out["N"][n]["f_final_vacuum"] < 1e-10 for n in N_LIST)
    print(f"\n(TB2a) 常時パネル全時代定義性（NaNゼロ）: {'通過' if tb2a else '不成立'}")
    print(f"(TB2b) 受動性（パネル有無でビット同一・N=12）: {'通過' if passive else '不成立'}")
    print(f"(TB2c) census決定性（帳簿外比<1e-6 全N）: {'通過' if tb2c else '不成立'}")
    print(f"(TB2d) 固有時間の対照（種=時計取得可・真空=取得不能）: "
          f"{'通過' if tb2d else '不成立'}")
    ok = tb2a and passive and tb2c and tb2d
    verdict = ("テストベッド中核成立: 常時G・census決定性・"
               "真空宇宙は空間ありτなし" if ok else "要精査")
    print(verdict)
    for n in N_LIST:
        out["N"][n].pop("f_series_seeded", None)
        out["N"][n].pop("f_series_vacuum", None)
    out.update({"TB2a": bool(tb2a), "TB2b_passive": passive,
                "TB2c": bool(tb2c), "TB2d": bool(tb2d),
                "verdict": verdict, "runtime_sec": time.time() - t0})
    (HERE / "result_tb2_matter_universe_series_v1.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False))
    print(f"完了 {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
