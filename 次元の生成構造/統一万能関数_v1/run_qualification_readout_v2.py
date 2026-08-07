#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""読出し v2（曖昧さ保存版）＋選択層の資格審査 — 全通過が使用条件

Q6 選択層の再現性: S∘G(v2) が G(v1) の確定値とビット単位で一致する
   （f₂・f_seed・r・census・位置 x/cover/present/pr_n・時計 ω/取得可否）。
   Nn=5（従来環境）と Nn=16（標準実験環境）の両方で検査する。
Q7 受動性: v2 パネル併走の有無で終状態がビット単位同一。
Q8 束の完全性（R7）: 位置の束が全ての巻き m=1…⌊Nn/2⌋ を重みつきで
   保持し、選択を含まないこと（両被覆の読み値が同時に取れること）。
Q9 混在度の定義健全性: 単一セル状態で per_wave_mix = 1、
   k セル等パワー状態で per_wave_mix = k（解析値との一致）。
Q10 不在の表現: 内容ゼロで content_power = 0（NaN でも真偽値でもない）
   かつ選択層 s_present が False を返す。

使い方: python3 run_qualification_readout_v2.py
"""
from __future__ import annotations
import importlib.util, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def build(ui, n, delta, Nn, Neta=8):
    """標準宇宙（任意 Nn）——unified_interaction の正本構成に一致させる。"""
    m = n * (n - 1) // 2
    _, v, _, _, _, _, _, Z0c, wp0 = ui.abl.build_init(n, False)
    r2 = ui.gen3.make_parent(n, seed=2)
    Csec = np.fft.fft(r2.relation_waves, axis=1) / n
    seed_state = Csec[:, 1] / np.linalg.norm(Csec[:, 1])
    C2_0 = np.zeros((m, Nn, Neta), complex)
    C2_0[:, 2, 0] = Z0c
    if delta > 0:
        C2_0[:, 1, 0] = delta * seed_state
    p2 = C2_0[:, 2, 0].real / np.linalg.norm(C2_0[:, 2, 0].real)
    q2 = C2_0[:, 2, 0].imag - (C2_0[:, 2, 0].imag @ p2) * p2
    nq = np.linalg.norm(q2)
    q2 = q2 / nq if nq > 1e-12 else np.zeros_like(p2)
    return ui.UnifiedEngine(n, C2_0, wp0), p2, q2


def main():
    t0 = time.time()
    ui = load("ui_q2", HERE / "unified_interaction_v1.py")
    v1 = load("ur_v1_q", HERE / "unified_readout_v1.py")
    v2 = load("ur_v2_q", HERE / "unified_readout_v2.py")
    S = load("sel_q", HERE / "selection_v1.py")
    res = {}

    # ---- Q6: S∘G(v2) ≡ G(v1)（Nn=5 と Nn=16）
    q6 = {}
    for Nn in (5, 16):
        eng, p2, q2 = build(ui, 12, 1e-2, Nn)
        c1 = {"C_flat": None, "c_gen": None}
        c2 = {"C_flat": None, "c_gen": None}
        worst = {k: 0.0 for k in ("f2", "f_seed", "r", "census", "x",
                                  "pr_n", "omega")}
        flags_ok = True
        for t in range(200):
            eng.step()
            C2 = eng.C2()
            a = v1.g_panel(C2, p2, q2, c1["C_flat"], c1["c_gen"])
            b = v2.g_panel(C2, p2, q2, c2["C_flat"], c2["c_gen"])
            c1, c2 = a["_carry"], b["_carry"]
            pos = S.s_position_maxmoment(b)
            clk = S.s_clock_acquirable(b)
            worst["f2"] = max(worst["f2"], abs(a["f2"] - b["f2"]))
            worst["f_seed"] = max(worst["f_seed"], abs(a["f_seed"] - b["f_seed"]))
            if "r" in a:
                worst["r"] = max(worst["r"], abs(a["r"] - b["r"]))
            worst["census"] = max(
                worst["census"], float(np.max(np.abs(a["census"] - b["cell_power"]))))
            # 台集合・実効本数（v1 g_wave_census との一致）
            if not np.array_equal(a["wave_count"], S.s_support_count(b)):
                flags_ok = False
            worst["census"] = max(worst["census"],
                                  float(np.max(np.abs(a["wave_pr_m"] - b["cell_pr_m"]))))
            if a["present"] != pos["present"]:
                flags_ok = False
            if a["present"]:
                worst["x"] = max(worst["x"], abs(a["x"] - pos["x"]))
                if a["cover"] != pos["cover"]:
                    flags_ok = False
                worst["pr_n"] = max(worst["pr_n"], abs(a["pr_n"] - b["pr_n"]))
            if a["acquirable"] != clk["acquirable"]:
                flags_ok = False
            if a["acquirable"]:
                worst["omega"] = max(worst["omega"], abs(a["omega"] - clk["omega"]))
        ok = flags_ok and all(v == 0.0 for v in worst.values())
        q6[f"Nn{Nn}"] = {"worst": worst, "flags_ok": flags_ok, "ok": bool(ok)}
        print(f"Q6 選択層再現（Nn={Nn}・T=200・最大差 "
              f"{max(worst.values()):.1e}・判定一致 {flags_ok}）: "
              f"{'通過' if ok else '不成立'}")
    res["Q6"] = all(q6[k]["ok"] for k in q6)
    res["Q6_detail"] = q6

    # ---- Q7: 受動性
    ea, p2, q2 = build(ui, 12, 1e-2, 16)
    eb, _, _ = build(ui, 12, 1e-2, 16)
    c = {"C_flat": None, "c_gen": None}
    for t in range(100):
        ea.step()
        p = v2.g_panel(ea.C2(), p2, q2, c["C_flat"], c["c_gen"])
        c = p["_carry"]
        eb.step()
    ok = bool(np.array_equal(ea.C, eb.C))
    res["Q7"] = ok
    print(f"Q7 受動性（v2パネル有無でビット同一・Nn=16・T=100）: "
          f"{'通過' if ok else '不成立'}")

    # ---- Q8: 束の完全性（全巻きを重みつきで保持・選択なし）
    b = v2.g_position_spectrum(ea.C2())
    Nn = ea.C2().shape[1]
    ok = (len(b["pos_m"]) == Nn // 2
          and np.array_equal(b["pos_m"], np.arange(1, Nn // 2 + 1))
          and np.all(np.isfinite(b["pos_weight"]))
          and np.allclose(b["pos_modulus"], Nn / b["pos_m"]))
    # 両被覆が同時に取れる（選択の余地が残されている）
    two = S.s_position_maxmoment(b), S.s_position_argmax(b)
    res["Q8"] = bool(ok)
    print(f"Q8 束の完全性（巻き {list(b['pos_m'][:4])}… 計{len(b['pos_m'])}本を"
          f"重みつき保持・選択子2種が適用可）: {'通過' if ok else '不成立'}")

    # ---- Q9: 混在度の定義健全性
    M, Nn9, Ne9 = 7, 6, 4
    C_single = np.zeros((M, Nn9, Ne9), complex)
    C_single[:, 1, 0] = 1.0
    m1 = v2.g_species_content(C_single)["per_wave_mix"]
    kcells = [(1, 0), (2, 0), (3, 1), (4, 2)]
    C_k = np.zeros((M, Nn9, Ne9), complex)
    for (k, e) in kcells:
        C_k[:, k, e] = 1.0
    mk = v2.g_species_content(C_k)["per_wave_mix"]
    ok = bool(np.allclose(m1, 1.0, atol=0, rtol=1e-15)
              and np.allclose(mk, len(kcells), atol=0, rtol=1e-15))
    res["Q9"] = ok
    print(f"Q9 混在度の定義健全性（単一セル→{m1[0]:.6f}＝1・"
          f"{len(kcells)}セル等パワー→{mk[0]:.6f}＝{len(kcells)}）: "
          f"{'通過' if ok else '不成立'}")

    # ---- Q10: 不在の表現
    C_empty = np.zeros((M, Nn9, Ne9), complex)
    C_empty[:, 2, 0] = 1.0        # 偶数帯のみ＝奇数帯内容は不在
    be = v2.g_position_spectrum(C_empty)
    ok = bool(be["content_power"] == 0.0 and be["pr_n"] == 0.0
              and not S.s_present(be) and np.all(np.isnan(be["pos_x"])))
    res["Q10"] = ok
    print(f"Q10 不在の表現（content_power={be['content_power']:.1e}・"
          f"重み0で表現・s_present=False）: {'通過' if ok else '不成立'}")

    allpass = all(v for k, v in res.items() if not k.endswith("detail"))
    print("読出しv2資格審査: " + ("ALL PASS — v2＋選択層は使用可"
                                  if allpass else "不成立あり"))
    res["all_pass"] = bool(allpass)
    res["runtime_sec"] = time.time() - t0
    (HERE / "qualification_readout_v2_result.json").write_text(
        json.dumps(res, indent=1, ensure_ascii=False, default=float))
    print(f"完了 {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
