#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""TB v2（やり直し版）: 通し創世テストベッド——常時実行G・時代の語り・census・対照統計・図化

前版(run_tb2_matter_universe_series_v1)の罪: ポンプ凝縮体を初期条件に手で置き、
通し創世（潜伏→インフレーション→凝縮）を飛ばした。本版は正本 run_genesis_v1 の
構成（真空=control親スライス2＋シード・スライス2不安定枝の f₂ が空間形成史）を
毛(η)込みで忠実に踏襲し、第0步から常時パネルで全史を記録する。

宇宙（正本構成・毛拡張）: C2_0[:,2,0]=Z0c（真空・毛0）＋ C2_0[:,1,0]=δ·種（毛0）。
  物質宇宙: δ=1e-2（明瞭生成・正本探索走行と同値）
  真空宇宙: δ=0（完全無シード・内在床点火——空間はできるが物質ゼロ）

常時パネル（毎ステップ・一様・切替なし・全てパラメータフリー関係量）:
  ノルム・f₂（スライス2面外分率=空間形成史）・f_seed（物質分率）・
  一段残差 r（時計）・毛×帯 census 帳簿・奇数内容の双対レジスタ位置 x・
  物質時計 ω̂（生成内容の位相前進）

判定（事前固定）:
 (TBa) 全時代定義性: 全ステップ・全パネルで NaN/Inf ゼロ。
 (TBb) 受動性: パネル有無で終状態ビット同一（N=12物質宇宙）。
 (TBc) 時代の語り（物質宇宙）: 常時記録の中に
       潜伏→crossing（f₂>0.05）→凝縮（f₂飽和）→物質時代
       （窓別生成率 g_meta > g_latency——時代変調）が現れる。
 (TBd) 真空対照: δ=0 でも crossing は起きる（空間はできる）が、
       f_seed は床レベル（<1e-10）に留まり物質時計は取得不能。
 (TBe) census決定性: 生成内容は毛0のみ（毛外比<1e-6）。
 (TBf) 対照統計: 安定窓 [2000,4000] の全パネルの平均・std・ドリフトを
       対照データとして記録（判定でなく取得義務）。
図: fig_tb_era_N{n}（時代パネル・種vs真空）・fig_tb_census_N{n}（帳簿）・
    fig_tb_nseries（N系列比較）——全てJSONから決定論再生成可能。
使い方: python3 run_tb_genesis_universe_v2.py [N...]（省略時 4 6 12）
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
MB = HERE.parent / "万能相互作用多体接続_v1"

spec3 = importlib.util.spec_from_file_location(
    "s3tbv2", MB / "run_stage3_sharedO_v2_and_hair_v1.py")
s3 = importlib.util.module_from_spec(spec3)
sys.modules[spec3.name] = s3
spec3.loader.exec_module(s3)
abl = s3.abl
gen3 = s3.gen3
s2 = s3.s2

NN, NETA = 5, 8
T_LONG = 4000
DELTA_MATTER = 1e-2
WIN_LAT = (50, 250)
WIN_META = (2000, 4000)


class HairEngine(s3.VertexEngineV2):
    """正本 stage3 HairEngine コピー移植（V2/TB 前版と同一定義）"""

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


def window_rate(ts, lnP, lo, hi):
    m_ = (ts >= lo) & (ts < hi)
    tt = ts[m_].astype(float)
    A = np.vstack([tt, np.ones_like(tt)]).T
    coef, _, _, _ = np.linalg.lstsq(A, lnP[m_], rcond=None)
    return float(coef[0])


def build_universe(n, delta):
    """正本 run_genesis_v1 の構成の毛拡張（乱数・構成規則は正本準拠）"""
    m = n * (n - 1) // 2
    _, v, _, _, _, _, _, Z0c, wp0 = abl.build_init(n, False)
    r2 = gen3.make_parent(n, seed=2)
    Csec = np.fft.fft(r2.relation_waves, axis=1) / n
    seed_state = Csec[:, 1] / np.linalg.norm(Csec[:, 1])
    C2_0 = np.zeros((m, NN, NETA), complex)
    C2_0[:, 2, 0] = Z0c
    if delta > 0:
        C2_0[:, 1, 0] = delta * seed_state
    p2 = C2_0[:, 2, 0].real / np.linalg.norm(C2_0[:, 2, 0].real)
    q2 = C2_0[:, 2, 0].imag - (C2_0[:, 2, 0].imag @ p2) * p2
    nq = np.linalg.norm(q2)
    q2 = q2 / nq if nq > 1e-12 else np.zeros_like(p2)
    return HairEngine(n, C2_0, wp0), p2, q2


def run_history(n, delta, panel_on=True):
    eng, p2, q2 = build_universe(n, delta)
    T = T_LONG
    H = {k: np.zeros(T) for k in ("f2", "fseed", "r", "x", "omega", "norm")}
    nonfinite = 0
    crossing = None
    prev = eng.C.flatten().copy()
    prev_c3 = None
    census_final = None
    for t in range(T):
        eng.step()
        if not panel_on:
            continue
        C2 = eng.C2()
        Z2 = C2[:, 2, 0]
        d2 = float(np.real(np.conj(Z2) @ Z2))
        Zp = Z2 - p2 * (p2 @ Z2) - q2 * (q2 @ Z2)
        f2 = float(np.real(np.conj(Zp) @ Zp)) / max(d2, 1e-300)
        P2 = np.abs(C2) ** 2
        ptot = float(P2.sum())
        podd = float(P2[:, 1, :].sum() + P2[:, 3, :].sum())
        fseed = podd / max(ptot, 1e-300)
        cf = eng.C.flatten()
        ip = np.vdot(prev, cf)
        ph = ip / abs(ip) if abs(ip) > 0 else 1.0
        rres = float(np.linalg.norm(cf - ph * prev))
        prev = cf.copy()
        mask = np.zeros(NN); mask[1] = 1.0; mask[3] = 1.0
        Wo = np.fft.ifft2(C2 * mask[None, :, None], axes=(1, 2)) * (NN * NETA)
        Pn = np.sum(np.abs(Wo) ** 2, axis=(0, 2))
        if Pn.sum() > 1e-280:
            z = np.sum(Pn * np.exp(2j * np.pi * np.arange(NN) / NN)) / Pn.sum()
            xr = float((np.angle(z) * NN / (2 * np.pi)) % NN)
        else:
            xr = -1.0
        c3 = C2[:, 3, :].reshape(-1)
        if prev_c3 is not None:
            zz = np.vdot(prev_c3, c3)
            om = (float(np.angle(zz))
                  if (abs(zz) > 1e-30 and np.linalg.norm(c3) > 1e-12) else np.nan)
        else:
            om = np.nan
        prev_c3 = c3.copy()
        vals = (f2, fseed, rres, ptot)
        if not all(np.isfinite(vv) for vv in vals):
            nonfinite += 1
        H["f2"][t] = f2; H["fseed"][t] = fseed; H["r"][t] = rres
        H["x"][t] = xr; H["omega"][t] = om; H["norm"][t] = ptot
        if crossing is None and f2 > 0.05:
            crossing = t + 1
        if t == T - 1:
            census_final = P2.sum(axis=0)
    return eng, H, crossing, nonfinite, census_final


def analyze(n, Hm, Hv, cr_m, cr_v, census):
    ts = np.arange(1, T_LONG + 1)
    lnP = np.log(np.maximum(Hm["fseed"], 1e-300))
    g_lat = window_rate(ts, lnP, *WIN_LAT)
    g_meta = window_rate(ts, lnP, *WIN_META)
    om = Hm["omega"][WIN_META[0]:WIN_META[1]]
    omv = om[np.isfinite(om)]
    clock_m = {"omega": float(np.mean(omv)) if len(omv) > 10 else None,
               "std": float(np.std(omv)) if len(omv) > 10 else None,
               "valid": float(len(omv) / len(om))}
    omv2 = Hv["omega"][WIN_META[0]:WIN_META[1]]
    clock_v_valid = float(np.mean(np.isfinite(omv2)))
    tot = float(census.sum())
    hair_out = float(census[:, 1:].sum()) / max(tot, 1e-300)
    stats = {}
    for k in ("f2", "fseed", "r", "norm"):
        w = Hm[k][WIN_META[0]:WIN_META[1]]
        drift = float(np.polyfit(np.arange(len(w)), w, 1)[0]) * len(w)
        stats[k] = {"mean": float(np.mean(w)), "std": float(np.std(w)),
                    "drift_over_window": drift}
    return {"crossing_matter": cr_m, "crossing_vacuum": cr_v,
            "g_latency": g_lat, "g_metastable": g_meta,
            "clock_matter": clock_m, "clock_vacuum_valid": clock_v_valid,
            "fseed_final_matter": float(Hm["fseed"][-1]),
            "fseed_final_vacuum": float(Hv["fseed"][-1]),
            "hair_outside_ratio": hair_out,
            "control_stats_meta_window": stats}


def draw_figs(n, Hm, Hv, census, res):
    ts = np.arange(1, T_LONG + 1)
    fig, axes = plt.subplots(3, 1, figsize=(9, 9), sharex=True)
    ax = axes[0]
    ax.semilogy(ts, np.maximum(Hm["f2"], 1e-16), label="物質宇宙 δ=1e-2")
    ax.semilogy(ts, np.maximum(Hv["f2"], 1e-16), label="真空宇宙 δ=0", ls="--")
    for c, lab in ((res["crossing_matter"], "crossing(物質)"),
                   (res["crossing_vacuum"], "crossing(真空)")):
        if c:
            ax.axvline(c, color="gray", lw=0.8)
    ax.set_ylabel("f₂（空間形成史）")
    ax.legend(); ax.set_title(f"N={n} 通し創世：常時実行パネルの全史")
    ax = axes[1]
    ax.semilogy(ts, np.maximum(Hm["fseed"], 1e-16), label="物質宇宙")
    ax.semilogy(ts, np.maximum(Hv["fseed"], 1e-16), label="真空宇宙", ls="--")
    ax.set_ylabel("f_seed（物質分率）"); ax.legend()
    ax = axes[2]
    ax.plot(ts, Hm["omega"], ".", ms=1, label="物質宇宙 ω̂")
    ax.plot(ts, Hv["omega"], ".", ms=1, label="真空宇宙 ω̂（不在）")
    ax.axhline(np.pi / 72, color="red", lw=0.8, label="π/72（柱1集団時計）")
    ax.set_ylabel("物質時計 ω̂"); ax.set_xlabel("t（step）"); ax.legend()
    fig.tight_layout()
    fig.savefig(HERE / f"fig_tb_era_N{n}_v2.png", dpi=130)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(8, 4))
    im = ax.imshow(np.log10(np.maximum(census, 1e-300)), aspect="auto",
                   origin="lower", cmap="viridis")
    ax.set_xlabel("毛 η（電荷巻き）"); ax.set_ylabel("帯 k")
    ax.set_title(f"N={n} 粒子census帳簿 log₁₀P（T={T_LONG}・物質宇宙）")
    fig.colorbar(im, ax=ax)
    fig.tight_layout()
    fig.savefig(HERE / f"fig_tb_census_N{n}_v2.png", dpi=130)
    plt.close(fig)


def main():
    t0 = time.time()
    ns = [int(a) for a in sys.argv[1:]] or [4, 6, 12]
    out = {"params": {"T": T_LONG, "delta_matter": DELTA_MATTER,
                      "Nn": NN, "Neta": NETA}, "N": {}}
    for n in ns:
        t1 = time.time()
        _, Hm, cr_m, nf_m, census = run_history(n, DELTA_MATTER)
        _, Hv, cr_v, nf_v, _ = run_history(n, 0.0)
        res = analyze(n, Hm, Hv, cr_m, cr_v, census)
        res["nonfinite"] = nf_m + nf_v
        draw_figs(n, Hm, Hv, census, res)
        np.savez_compressed(HERE / f"tb_genesis_history_N{n}_v2.npz",
                            **{f"m_{k}": Hm[k] for k in Hm},
                            **{f"v_{k}": Hv[k] for k in Hv},
                            census=census)
        out["N"][n] = res
        cm = res["clock_matter"]
        print(f"N={n}: crossing 物質={cr_m} 真空={cr_v} "
              f"f_seed末 物質={res['fseed_final_matter']:.2e} "
              f"真空={res['fseed_final_vacuum']:.2e} "
              f"時計 物質=ω{cm['omega'] if cm['omega'] else float('nan'):.4f}"
              f"(有効{cm['valid']:.2f}) 真空有効={res['clock_vacuum_valid']:.2f} "
              f"毛外={res['hair_outside_ratio']:.1e} [{time.time()-t1:.0f}s]")

    # 受動性（N=最初の要素の物質宇宙）
    n0 = ns[0] if 12 not in ns else 12
    ea, _, _, _, _ = run_history(n0, DELTA_MATTER, panel_on=True)
    eb, _, _, _, _ = run_history(n0, DELTA_MATTER, panel_on=False)
    passive = bool(np.array_equal(ea.C, eb.C))

    tba = all(out["N"][n]["nonfinite"] == 0 for n in ns)
    tbc = all(out["N"][n]["crossing_matter"] is not None
              and out["N"][n]["g_metastable"] > out["N"][n]["g_latency"]
              for n in ns)
    tbd = all(out["N"][n]["crossing_vacuum"] is not None
              and out["N"][n]["fseed_final_vacuum"] < 1e-10
              and out["N"][n]["clock_vacuum_valid"] < 0.5
              and out["N"][n]["clock_matter"]["valid"] > 0.9 for n in ns)
    tbe = all(out["N"][n]["hair_outside_ratio"] < 1e-6 for n in ns)
    print(f"\n(TBa) 全時代定義性: {'通過' if tba else '不成立'}")
    print(f"(TBb) 受動性（N={n0}ビット同一）: {'通過' if passive else '不成立'}")
    print(f"(TBc) 時代の語り（crossing→時代変調 g_meta>g_lat）: "
          f"{'通過' if tbc else '不成立'}")
    print(f"(TBd) 真空対照（空間あり・物質なし・時計なし）: "
          f"{'通過' if tbd else '不成立'}")
    print(f"(TBe) census決定性（毛外<1e-6）: {'通過' if tbe else '不成立'}")
    print(f"(TBf) 対照統計: 安定窓の全パネル統計をJSONに記録")
    ok = tba and passive and tbc and tbd and tbe
    verdict = "通し創世テストベッド成立" if ok else "要精査"
    print(verdict)
    out.update({"TBa": bool(tba), "TBb_passive": passive, "TBc": bool(tbc),
                "TBd": bool(tbd), "TBe": bool(tbe),
                "verdict": verdict, "runtime_sec": time.time() - t0})
    (HERE / "result_tb_genesis_universe_v2.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False))
    print(f"完了 {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
