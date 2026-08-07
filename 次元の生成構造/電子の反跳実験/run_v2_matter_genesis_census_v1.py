#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""V2: 物質生成の特性測定——シード依存・分布・量制御・再現性（実験環境の設計図）

問い（木原指定）: シードから生まれる物質の分布は未計測。余分な物質が大量に
ランダムに生まれる環境では実験にならない。安定した中性フェルミオンだけの
再現性ある環境（物理的真空）が作れるかをテストする。

正本: 万能相互作用多体接続_v1/run_stage3_sharedO_v2_and_hair_v1.py の
毛(η)拡張 census（HairEngine をコピー移植・和則 m*=2m_pump−m_seed・
排他比 6.1e289 実測済み）。ポンプ=control親(k=2)・種=零閉鎖状態×δ(k=1)。

掃引: 毛レシピ {中性レシピ(m_s=0,m_B=0→予言m*=0), 正本対照(m_s=2,m_B=1→m*=3)}
      × δ ∈ {1e-4, 1e-3, 1e-2, 3e-2, 1e-1}。N_GRAPH=12（分解能は本実験では
      固定し、Nでの掃引は次段に登録）。T=300・スナップショット10。

判定（事前固定）:
 (V2a) 毛制御: 中性レシピの生成内容（k=3相棒帯）の毛集中度
       P(m=0)/ΣP ≥ 0.9——「中性フェルミオンだけ」が狙って作れる。
 (V2b) 量制御: δ 掃引の中に「有限飽和」窓が存在する——f が可測に成長
       （>1e-8）しつつ後期成長率が小（|Δln f|/Δt < 1e-3）で、氾濫しない
       （f_final < 0.1）δ が少なくとも一つある。
 (V2c) 再現性: 同一構成の再走行が決定論的に一致（相対差 <1e-12）。
       種位相を変えた近傍シードでも毛集中度 ≥0.9 を維持。
 (V2d) 空間無傷: 全ノルム相対ドリフト <1e-9・ポンプ（凝縮体）のパワー比が
       制御窓の δ で過半を維持。
使い方: python3 run_v2_matter_genesis_census_v1.py
"""
from __future__ import annotations
import importlib.util, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
MB = HERE.parent / "万能相互作用多体接続_v1"

spec3 = importlib.util.spec_from_file_location(
    "s3v2", MB / "run_stage3_sharedO_v2_and_hair_v1.py")
s3 = importlib.util.module_from_spec(spec3)
sys.modules[spec3.name] = s3
spec3.loader.exec_module(s3)
abl = s3.abl
s2 = s3.s2

spec_g3 = importlib.util.spec_from_file_location(
    "g3v2", MB / "run_genesis_v3_register_local_v1.py")
g3 = importlib.util.module_from_spec(spec_g3)
sys.modules[spec_g3.name] = g3
spec_g3.loader.exec_module(g3)

N_GRAPH = 12
NN, NETA = 5, 8
T_RUN = 300
SNAP = 30
DELTAS = (1e-4, 1e-3, 1e-2, 3e-2, 1e-1)


class HairEngine(s3.VertexEngineV2):
    """正本 stage3 main() 内 HairEngine のコピー移植（Nn/Nη を引数化）"""

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


def run_case(n, m, Z0c, wp0, seed_state, m_s, m_B, delta, phase=0.0):
    C2_0 = np.zeros((m, NN, NETA), complex)
    C2_0[:, 2, m_s % NETA] = Z0c / np.linalg.norm(Z0c)
    C2_0[:, 1, m_B % NETA] = delta * seed_state * np.exp(1j * phase)
    eng = HairEngine(n, C2_0, wp0)
    norm0 = float(np.sum(np.abs(eng.C) ** 2))
    snaps = []
    for t in range(1, T_RUN + 1):
        eng.step()
        if t % SNAP == 0:
            C2 = eng.C2()
            # 生成内容 = 相棒帯 k=3（毛分解）と奇数帯全体
            P_k3 = np.sum(np.abs(C2[:, 3, :]) ** 2, axis=0)
            P_odd = float(np.sum(np.abs(C2[:, 1, :]) ** 2)
                          + np.sum(np.abs(C2[:, 3, :]) ** 2))
            P_pump = float(np.sum(np.abs(C2[:, 2, :]) ** 2))
            P_tot = float(np.sum(np.abs(C2) ** 2))
            snaps.append({"t": t, "P_k3_by_m": [float(x) for x in P_k3],
                          "f_odd": P_odd / P_tot, "pump_frac": P_pump / P_tot,
                          "norm": P_tot})
    norm_drift = abs(snaps[-1]["norm"] - norm0) / norm0
    m_star = (2 * m_s - m_B) % NETA
    Pk3 = np.array(snaps[-1]["P_k3_by_m"])
    tot3 = float(Pk3.sum())
    conc = float(Pk3[m_star] / tot3) if tot3 > 1e-300 else float("nan")
    f_series = [s["f_odd"] for s in snaps]
    # 後期成長率（ln f の傾き / step）
    if f_series[-2] > 0 and f_series[-1] > 0:
        late_growth = abs(np.log(f_series[-1] / f_series[-2])) / SNAP
    else:
        late_growth = float("nan")
    return {"m_star_pred": int(m_star), "hair_conc": conc,
            "f_final": f_series[-1], "f_series": f_series,
            "late_growth": late_growth, "pump_frac_final": snaps[-1]["pump_frac"],
            "norm_drift": norm_drift, "C_final_sum": complex(np.sum(eng.C)).real}


def main():
    t0 = time.time()
    n = N_GRAPH
    _, v, _, _, _, _, _, Z0c, wp0 = abl.build_init(n, False)
    Z0c = Z0c / np.linalg.norm(Z0c)
    m = n * (n - 1) // 2
    seed_state = g3.zero_closure_state(m, np.random.default_rng(98000))

    recipes = {"中性(ms0,mB0)": (0, 0), "正本対照(ms2,mB1)": (2, 1)}
    out = {"params": {"N": n, "M": m, "Nn": NN, "Neta": NETA, "T": T_RUN,
                      "deltas": list(DELTAS)}, "cases": {}}
    print(f"{'レシピ':>16} {'δ':>7} {'予言m*':>5} {'毛集中':>7} {'f_final':>10} "
          f"{'後期成長':>9} {'ポンプ比':>7} {'ノルム':>9}")
    for rname, (m_s, m_B) in recipes.items():
        for d in DELTAS:
            r = run_case(n, m, Z0c, wp0, seed_state, m_s, m_B, d)
            out["cases"][f"{rname}|{d:g}"] = r
            print(f"{rname:>16} {d:>7g} {r['m_star_pred']:>5} "
                  f"{r['hair_conc']:>7.4f} {r['f_final']:>10.3e} "
                  f"{r['late_growth']:>9.2e} {r['pump_frac_final']:>7.4f} "
                  f"{r['norm_drift']:>9.2e}")

    # V2a: 中性レシピの毛集中（可測に生成した δ のうち最大の集中度で判定）
    neut = [out["cases"][k] for k in out["cases"] if k.startswith("中性")]
    neut_meas = [c for c in neut if c["f_final"] > 1e-8]
    v2a = bool(neut_meas and max(c["hair_conc"] for c in neut_meas) >= 0.9)
    # V2b: 有限飽和窓
    ctrl = [c for c in neut
            if c["f_final"] > 1e-8 and c["f_final"] < 0.1
            and c["late_growth"] == c["late_growth"] and c["late_growth"] < 1e-3]
    v2b = bool(ctrl)
    # V2c: 決定論再現＋近傍シード
    m_s, m_B = recipes["中性(ms0,mB0)"]
    d_ref = 1e-2
    r1 = run_case(n, m, Z0c, wp0, seed_state, m_s, m_B, d_ref)
    r2 = run_case(n, m, Z0c, wp0, seed_state, m_s, m_B, d_ref)
    det = abs(r1["C_final_sum"] - r2["C_final_sum"]) / max(abs(r1["C_final_sum"]), 1e-300)
    rph = run_case(n, m, Z0c, wp0, seed_state, m_s, m_B, d_ref, phase=0.7)
    v2c = bool(det < 1e-12 and rph["hair_conc"] >= 0.9)
    # V2d: 空間無傷
    v2d = all(c["norm_drift"] < 1e-9 for c in neut) and \
        bool(ctrl and all(c["pump_frac_final"] > 0.5 for c in ctrl))
    print(f"\n(V2a) 毛制御（中性レシピ集中度≥0.9）: {'通過' if v2a else '不成立'}")
    print(f"(V2b) 量制御（有限飽和窓の存在）: {'通過' if v2b else '不成立'}"
          f"  制御窓δ={[f'{c!r}' for c in []] or [k.split('|')[1] for k in out['cases'] if k.startswith('中性') and out['cases'][k] in ctrl]}")
    print(f"(V2c) 再現性（決定論 {det:.2e}・近傍位相の集中 {rph['hair_conc']:.3f}）: "
          f"{'通過' if v2c else '不成立'}")
    print(f"(V2d) 空間無傷（ノルム<1e-9・制御窓でポンプ過半）: {'通過' if v2d else '不成立'}")
    ok = v2a and v2b and v2c and v2d
    verdict = ("物理的真空のレシピ成立: 中性フェルミオン狙い撃ち・量制御窓あり・"
               "決定論再現・空間無傷" if ok else "要精査")
    print(verdict)
    out.update({"V2a": bool(v2a), "V2b": bool(v2b), "V2c": bool(v2c),
                "V2d": bool(v2d), "determinism": det,
                "neighbor_conc": rph["hair_conc"],
                "verdict": verdict, "runtime_sec": time.time() - t0})
    (HERE / "result_v2_matter_genesis_census_v1.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False))
    print(f"完了 {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
