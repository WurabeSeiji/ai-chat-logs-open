#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""N 掃引 1→20・電子型フェルミオン（巻き m*=−3・電荷 Q=−1）——中性掃引と同一条件の対照実験

位置づけ: 中性フェルミオン掃引（run_tb_nsweep_1to20_v1.py・対照テスト済み
静止点 929554cf）と**全条件を揃えて**、シードの巻きだけを変えた対照走行。
既存の電子版データは存在しないため再現一致検定ではなく、中性版との比較対照である。

**無改変 import の方針（シリーズ規約）**
中性掃引スクリプトを本フォルダへ無改変コピー（md5 bfa5d854b637fe97c33a7148be9c7f86）
し、モジュールとして import して `run_one` / `summarize` / `fig_one` / `fig_summary`
/ `fig_matrix` をそのまま使う。走行ロジック・記録項目・判定・図は一行も書き換えない。
コピーの HERE は本フォルダなので、出力は本フォルダに出る（原本フォルダは無傷）。

**変更点は 2 つだけ**
1. 宇宙構成: `F.build_standard_universe` を、巻きを宣言できる版に差し替える
   （メモリ上の属性差し替えのみ・ファイルは改変しない）。
   和則 m* = 2·m_pump − m_seed。本実験は (m_pump, m_seed) = (0, 3) → **m* = −3**。
   ポンプ巻きは中性掃引と同一の 0 なので、差はシード巻きの 1 点のみ。
   Nη=8 ゆえ m*=−3 は毛インデックス 5（−3 mod 8）に立つ。
2. 記録の追加（力学・条件には触れない・記録のみ）:
   (E1) 巻き集中度: 相棒帯 k=3 における P(m*=−3)/ΣP（V2a と同型）。全帯版も併記。
   (E2) 電荷の可読性: 相棒帯 k=3 の優勢巻き m̂（符号付き）から Q̂ = m̂/3、
        および 3 | m̂（単独可読性・H1b と同型）。

記録は `RecordingEngine.step()` が毎步 1 回だけ追記する（`_readout`/`_nonlinear` が
内部で `C2()` を複数回呼ぶため、C2() 側での記録は步と対応しない——step 側で取る）。

条件（中性掃引と完全同一・すべて宣言値）:
  F=unified_interaction_v2.py / D=unified_dimension_v1.py /
  G=unified_readout_v3.py / S=selection_v1.py
  Nn=16・Nη=8・T=4000・δ=1e-2・親 seed=2 固定・cell=(2,0)・order=6・窓[2000,4000]
  各 N で真空対照 δ=0 も走行・N=1..20

事前登録（実行前固定）:
  (P1) 構成できない N は中性版と同一（= [1, 2, 8, 10]）。シードの巻きは
       make_parent／build_init に影響しないため。不一致なら異常。
  (P2) 空間側（f2・tau_space・align・n_eff・closure）は中性版と一致するか、
       一致しない場合はその差を記録する（巻きは毛レジスタの位置のみを変える）。
  (P3) E1 の k=3 帯集中度が ≥0.9 なら「電子型を狙って作れた」と記録する
       （中性版 V2a の閾値と同じ）。判定は記録であり合否ではない。

使い方: python3 run_nsweep_electron_v1.py [Nmin Nmax]（省略時 1 20）
"""
from __future__ import annotations

import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent

# --- 中性掃引スクリプト（無改変コピー）を import ---------------------------
_spec = importlib.util.spec_from_file_location(
    "ns_sweep", HERE / "run_tb_nsweep_1to20_v1.py")
ns = importlib.util.module_from_spec(_spec)
sys.modules["ns_sweep"] = ns
_spec.loader.exec_module(ns)

assert ns.HERE == HERE, "コピーの HERE が本フォルダでない（出力先が原本になる）"

# --- 電子レシピ（宣言値） --------------------------------------------------
M_PUMP, M_SEED = 0, 3
NETA = ns.NETA
M_STAR = (2 * M_PUMP - M_SEED)          # = −3
M_STAR_IDX = M_STAR % NETA              # = 5
assert M_STAR == -3, f"和則の帰結が −3 でない: {M_STAR}"
PARTNER_K = 3                           # 相棒帯（V2a と同じ）

_ENGINES: list = []                     # 直近に構成したエンジンを受け渡す


def _signed(idx: int) -> int:
    """毛インデックス → 符号付き巻き（Nη=8 なら 5 → −3・3 → +3）。"""
    h = NETA // 2
    return int(((idx + h) % NETA) - h)


def measure_hair(C2: np.ndarray) -> tuple:
    """(E1)(E2) の記録量。判定も選択もしない（値をそのまま返す）。"""
    P = np.abs(C2) ** 2
    tot = float(P.sum())
    Pe = P.sum(axis=(0, 1))                      # 巻きごとの全帯パワー
    conc_all = float(Pe[M_STAR_IDX] / tot) if tot > 0 else float("nan")
    if C2.shape[1] > PARTNER_K:
        P3 = P[:, PARTNER_K, :].sum(axis=0)      # 相棒帯 k=3 の巻き分布
        s3 = float(P3.sum())
        conc_k3 = float(P3[M_STAR_IDX] / s3) if s3 > 0 else float("nan")
        dom = _signed(int(np.argmax(P3))) if s3 > 0 else 0
    else:
        conc_k3, dom, s3 = float("nan"), 0, 0.0
    q_hat = dom / 3.0
    readable = 1.0 if (s3 > 0 and dom % 3 == 0) else 0.0
    return conc_all, conc_k3, float(dom), q_hat, readable, s3


class RecordingEngine(ns.F.UnifiedEngineV2):
    """力学は UnifiedEngineV2 そのまま。step 後に毛の記録を 1 回だけ追記する。"""

    def __init__(self, *a, **kw):
        super().__init__(*a, **kw)
        self._rec: list = []

    def step(self, *a, **kw):
        out = super().step(*a, **kw)
        self._rec.append(measure_hair(self.C2()))
        return out


def build_electron_universe(n, delta, Nn=5, Neta=8, seed=2):
    """F v2 の build_standard_universe と同一手順。巻きの置き場所だけを変える。"""
    m = n * (n - 1) // 2
    _, v, _, _, _, _, _, Z0c, wp0 = ns.F.abl.build_init(n, False)
    r2 = ns.F.gen3.make_parent(n, seed=seed)
    Csec = np.fft.fft(r2.relation_waves, axis=1) / n
    seed_state = Csec[:, 1] / np.linalg.norm(Csec[:, 1])
    C2_0 = np.zeros((m, Nn, Neta), complex)
    ip, isd = M_PUMP % Neta, M_SEED % Neta
    C2_0[:, 2, ip] = Z0c
    C2_0[:, 1, isd] = delta * seed_state          # δ=0 なら零（分岐不要）
    p2 = C2_0[:, 2, ip].real / np.linalg.norm(C2_0[:, 2, ip].real)
    q2 = C2_0[:, 2, ip].imag - (C2_0[:, 2, ip].imag @ p2) * p2
    with np.errstate(divide="ignore", invalid="ignore"):
        q2 = q2 / np.linalg.norm(q2)
    eng = RecordingEngine(n, C2_0, wp0)
    _ENGINES.append(eng)
    return eng, p2, q2


ns.F.build_standard_universe = build_electron_universe   # メモリ上のみ差し替え


def hair_arrays(eng) -> dict:
    a = np.array(eng._rec, dtype=float)             # (steps, 6)
    return {"conc_all": a[:, 0], "conc_k3": a[:, 1], "dom_m": a[:, 2],
            "q_hat": a[:, 3], "readable": a[:, 4], "k3_power": a[:, 5]}


def med_win(x):
    w = x[slice(*ns.WIN)]
    w = w[np.isfinite(w)]
    return float(np.median(w)) if len(w) else float("nan")


def fig_hair(n, hm, hv, rec):
    import matplotlib.pyplot as plt
    ts = np.arange(1, len(hm["conc_k3"]) + 1)
    fig, ax = plt.subplots(3, 1, figsize=(8.5, 8), sharex=True)
    ax[0].plot(ts, hm["conc_k3"], lw=0.8, label=f"k={PARTNER_K}帯 物質")
    ax[0].plot(ts, hm["conc_all"], lw=0.8, color="tab:orange", label="全帯 物質")
    ax[0].plot(ts, hv["conc_k3"], "k--", lw=0.8, label=f"k={PARTNER_K}帯 真空")
    ax[0].axhline(0.9, color="tab:green", lw=0.8, ls=":", label="0.9（V2a 閾値）")
    ax[0].set_ylabel(f"(E1) 巻き集中度 P(m={M_STAR})/ΣP")
    ax[0].legend(fontsize=7)
    ax[0].set_title(f"N={n}  電子型レシピ (m_pump,m_seed)=({M_PUMP},{M_SEED}) → "
                    f"m*={M_STAR}  集中度中央値={rec['conc_k3_med']:.4f}  "
                    f"Q̂={rec['q_hat_med']:+.4f}")
    ax[1].plot(ts, hm["dom_m"], lw=0.8, label="優勢巻き m̂（物質・k=3）")
    ax[1].axhline(M_STAR, color="tab:red", lw=0.8, ls=":", label=f"m*={M_STAR}")
    ax[1].set_ylabel("(E2) 優勢巻き m̂"); ax[1].legend(fontsize=7)
    ax[2].semilogy(ts, np.maximum(hm["k3_power"], 1e-40), lw=0.8, label="物質")
    ax[2].semilogy(ts, np.maximum(hv["k3_power"], 1e-40), "k--", lw=0.8, label="真空")
    ax[2].set_ylabel(f"k={PARTNER_K}帯 パワー"); ax[2].set_xlabel("τ（step）")
    ax[2].legend(fontsize=7)
    for a in ax:
        a.axvspan(ns.WIN[0], ns.WIN[1], color="green", alpha=0.06)
    fig.tight_layout()
    fig.savefig(HERE / f"fig_electron_hair_N{n}_v1.png", dpi=110)
    plt.close(fig)


def main():
    t0 = time.time()
    nmin = int(sys.argv[1]) if len(sys.argv) > 2 else 1
    nmax = int(sys.argv[2]) if len(sys.argv) > 2 else 20
    print(f"=== N 掃引 {nmin}→{nmax}・電子型 m*={M_STAR}（毛idx {M_STAR_IDX}）"
          f"・Nn={ns.NN}・Nη={NETA}・δ={ns.DELTA}・seed={ns.SEED}・T={ns.T} ===")
    recs, fails = [], []
    out = {"env": {"Nn": ns.NN, "Neta": NETA, "T": ns.T, "delta": ns.DELTA,
                   "seed": ns.SEED, "cell": list(ns.CELL), "order": ns.ORDER,
                   "window": list(ns.WIN),
                   "functions": ["unified_interaction_v2", "unified_dimension_v1",
                                 "unified_readout_v3", "selection_v1"],
                   "recipe": {"m_pump": M_PUMP, "m_seed": M_SEED,
                              "m_star": M_STAR, "hair_index": M_STAR_IDX,
                              "partner_band_k": PARTNER_K,
                              "sum_rule": "m* = 2*m_pump - m_seed"},
                   "base_script_md5": "bfa5d854b637fe97c33a7148be9c7f86"},
           "N": {}, "failed": {}}
    for n in range(nmin, nmax + 1):
        t1 = time.time()
        _ENGINES.clear()
        try:
            Hm, Rm, Am, Ccm, Csm = ns.run_one(n, ns.DELTA)
            eng_m = _ENGINES[-1]
            Hv, Rv, Av, Ccv, Csv = ns.run_one(n, 0.0)
            eng_v = _ENGINES[-1]
        except Exception as ex:
            msg = f"{type(ex).__name__}: {ex}"
            fails.append(n)
            out["failed"][n] = msg
            out["N"][n] = {"N": n, "M": n * (n - 1) // 2, "built": False,
                           "error": msg}
            print(f"N={n:3d} M={n*(n-1)//2:4d}: **構成不能** {msg[:70]}")
            (HERE / "result_nsweep_electron_v1.json").write_text(
                json.dumps(out, indent=1, ensure_ascii=False, default=float))
            continue
        hm, hv = hair_arrays(eng_m), hair_arrays(eng_v)
        assert len(hm["conc_k3"]) == ns.T, f"記録数 {len(hm['conc_k3'])} ≠ T"
        rec = ns.summarize(n, Hm, Rm, Am, Ccm, Csm, Hv, Av)
        rec.update({
            "conc_k3_med": med_win(hm["conc_k3"]),
            "conc_all_med": med_win(hm["conc_all"]),
            "dom_m_med": med_win(hm["dom_m"]),
            "q_hat_med": med_win(hm["q_hat"]),
            "readable_rate": med_win(hm["readable"]),
            "k3_power_med": med_win(hm["k3_power"]),
            "conc_k3_med_vacuum": med_win(hv["conc_k3"]),
            "k3_power_med_vacuum": med_win(hv["k3_power"]),
            "electron_made": bool(med_win(hm["conc_k3"]) >= 0.9),
        })
        recs.append(rec); out["N"][n] = rec
        ns.fig_one(n, Hm, Hv, Am, Ccm, Csm, rec)
        fig_hair(n, hm, hv, rec)
        np.savez_compressed(
            HERE / f"nsweep_electron_N{n}_v1.npz",
            **{f"m_{k}": Hm[k] for k in ns.KEYS},
            **{f"v_{k}": Hv[k] for k in ns.KEYS},
            m_resid=Rm, m_acq=Am, v_acq=Av,
            m_cond_closure=Ccm, m_seed_closure=Csm,
            **{f"hair_m_{k}": v for k, v in hm.items()},
            **{f"hair_v_{k}": v for k, v in hv.items()})
        print(f"N={n:3d} M={rec['M']:4d}: 空間τ={str(rec['tau_space']):>5} "
              f"物質={str(rec['matter_born']):>5} 時間τ={str(rec['tau_time']):>5} "
              f"align={rec['align_med']:.4f} n_eff={rec['n_eff_med']:.3f} "
              f"| 集中度={rec['conc_k3_med']:.4f} m̂={rec['dom_m_med']:+.1f} "
              f"Q̂={rec['q_hat_med']:+.3f} 可読={rec['readable_rate']:.2f} "
              f"[{time.time()-t1:.0f}s]")
        (HERE / "result_nsweep_electron_v1.json").write_text(
            json.dumps(out, indent=1, ensure_ascii=False, default=float))
    if recs:
        ns.fig_summary(recs, fails)
    ns.fig_matrix(recs, fails, nmin, nmax)
    out["failed_N"] = fails
    out["runtime_sec"] = time.time() - t0
    (HERE / "result_nsweep_electron_v1.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False, default=float))
    print(f"\n構成できなかった N: {fails if fails else 'なし'}")
    print(f"完了 {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
