#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""N 掃引・大 N 側（N=40, 100, 300）— 上位構造オーダの境界を探す

N=1..20 の掃引（run_tb_nsweep_1to20_v1.py）では、下端の境界（N=5 で第3次元が
立つ）しか見えず、上位構造オーダへの相転移らしき境界は現れなかった。また
align・n_eff は N に対して単調でなく（N=13, 20 で悪化）、曖昧さ 1−align は
1/M 参照線に乗らなかった。**N=12 が分子オーダかどうかは N≤20 では判定できない**
——そこで N を桁で伸ばす。

条件・関数・記録項目・図はすべて run_tb_nsweep_1to20_v1.py と同一
（同じシード seed=2 固定・Nn=16・Nη=8・δ=1e-2・T=4000・cell=(2,0)・order=6・
窓[2000,4000]）。**同じプログラム（同モジュールの関数）を再利用**して条件の
同一性を担保する。1〜20 の結果ファイルは上書きしない（別名で保存）。

構成が失敗する N（親構成の失敗など）は、その事実を結果として記録し継続する。

計算量の目安（実測 N=20/M=190 が2条件で 80 秒）:
  N=40  M=780    約 4 倍
  N=100 M=4950   約 26 倍
  N=300 M=44850  約 236 倍（数時間規模）

使い方: python3 run_tb_nsweep_large_v1.py [N ...]（省略時 40 100 300）
"""
from __future__ import annotations
import importlib.util, json, shutil, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


# 1..20 掃引と同一のプログラムを再利用する（条件の同一性を担保）
NS = load("nsweep_base", HERE / "run_tb_nsweep_1to20_v1.py")


def main():
    t0 = time.time()
    ns_list = [int(a) for a in sys.argv[1:]] or [40, 100, 300]
    print(f"=== N 掃引・大N側 {ns_list}（Nn={NS.NN}・Nη={NS.NETA}・δ={NS.DELTA}・"
          f"seed={NS.SEED} 固定・T={NS.T}）===", flush=True)
    recs, fails = [], []
    out = {"env": {"Nn": NS.NN, "Neta": NS.NETA, "T": NS.T, "delta": NS.DELTA,
                   "seed": NS.SEED, "cell": list(NS.CELL), "order": NS.ORDER,
                   "window": list(NS.WIN),
                   "functions": ["unified_interaction_v2", "unified_dimension_v1",
                                 "unified_readout_v3", "selection_v1"],
                   "base_program": "run_tb_nsweep_1to20_v1.py（同一関数を再利用）"},
           "N": {}, "failed": {}}
    res_path = HERE / "result_tb_nsweep_large_v1.json"
    for n in ns_list:
        t1 = time.time()
        print(f"N={n} 開始（M={n*(n-1)//2}）…", flush=True)
        try:
            Hm, Rm, Am, Ccm, Csm = NS.run_one(n, NS.DELTA)
            print(f"  物質宇宙 完了 [{time.time()-t1:.0f}s]", flush=True)
            Hv, Rv, Av, Ccv, Csv = NS.run_one(n, 0.0)
            print(f"  真空宇宙 完了 [{time.time()-t1:.0f}s]", flush=True)
        except Exception as ex:
            msg = f"{type(ex).__name__}: {ex}"
            fails.append(n)
            out["failed"][n] = msg
            out["N"][n] = {"N": n, "M": n * (n - 1) // 2, "built": False,
                           "error": msg}
            print(f"N={n:4d}: **構成不能** {msg[:90]}", flush=True)
            res_path.write_text(json.dumps(out, indent=1, ensure_ascii=False,
                                           default=float))
            continue
        rec = NS.summarize(n, Hm, Rm, Am, Ccm, Csm, Hv, Av)
        recs.append(rec)
        out["N"][n] = rec
        NS.fig_one(n, Hm, Hv, Am, Ccm, Csm, rec)
        np.savez_compressed(HERE / f"tb_nsweep_N{n}_v1.npz",
                            **{f"m_{k}": Hm[k] for k in NS.KEYS},
                            **{f"v_{k}": Hv[k] for k in NS.KEYS},
                            m_resid=Rm, m_acq=Am, v_acq=Av,
                            m_cond_closure=Ccm, m_seed_closure=Csm)
        print(f"N={n:4d} M={rec['M']:6d}: 空間τ={str(rec['tau_space']):>5} "
              f"物質={str(rec['matter_born']):>5} 時間τ={str(rec['tau_time']):>5} "
              f"時計定着τ={str(rec['tau_lock']):>5} "
              f"align={rec['align_med']:.4f} n_eff={rec['n_eff_med']:.3f} "
              f"凝縮体閉塞={rec['cond_closure_med']:.2e} [{time.time()-t1:.0f}s]",
              flush=True)
        res_path.write_text(json.dumps(out, indent=1, ensure_ascii=False,
                                       default=float))
    if recs:
        NS.fig_summary(recs, fails)
        shutil.move(str(HERE / "fig_nsweep_summary_v1.png"),
                    str(HERE / "fig_nsweep_summary_large_v1.png"))
    out["failed_N"] = fails
    out["runtime_sec"] = time.time() - t0
    res_path.write_text(json.dumps(out, indent=1, ensure_ascii=False, default=float))
    print(f"\n構成できなかった N: {fails if fails else 'なし'}", flush=True)
    print(f"完了 {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
