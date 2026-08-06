#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""予備実験 v7: 海構成規約の確立＋読める電荷の簿記恒等式

Part A（規約）: v6でS1/S2（人工海）が Q_wind 保存を破った原因を特定する。
仮説: χ解析射影（single_winding）を全成分に課した構成のみ保存クラスに載る。
検定: 純粋+1種 + 海の4変種で Q_wind の保存を比較（J=2000）:
  sea_SW: single_winding済み中性束 + η射影{0}（規約候補）
  sea_raw: χ両符号の中性束 + η射影{0}（v4b構成・保存破れの再現）
判定 H_conv: sea_SW は CV<0.02 で保存、sea_raw は破れる → 規約確立。

Part B（簿記恒等式）: 保存構成で
  Q_wind = Σ m·P（厳密保存・v6実証）
  Q3 = Σ fold3(m)·P（読める正味電荷）
  W = Σ ((m−fold3(m))/3)·P（mod3中性複合に隠れた電荷の帳簿）
  恒等式: Q_wind = Q3 + 3W ⟹ ΔQ3 = −3ΔW（厳密のはず・数値確認）
物理的読み: 読める電荷の減少は消失ではなく、分母3時計に中性と読める
複合（m=+3 相棒=中性子的）への厳密な持ち込みである。
判定 H_ledger: |ΔQ3 + 3ΔW| / |ΔQ3| < 1e-6（全窓）。

Part C（正当な寿命比較・±1一意性の再検）: 規約準拠の海で
  S1' (+1+海SW) vs S2' (+2+海SW) の支配種残存を再測定
（v4bの比較は保存クラス外で資格が弱かった。今回は有効な計器で）。

使い方: python3 run_pre_sea_convention_ledger_v7.py
"""
from __future__ import annotations
import importlib.util, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
UIM = HERE.parent / "万能非弾性写像_managed_v1"
spec = importlib.util.spec_from_file_location("exact_v7", UIM / "run_ignition_fate_exact_v3.py")
ex = importlib.util.module_from_spec(spec); sys.modules[spec.name] = ex
spec.loader.exec_module(ex)
v1, toy, base = ex.v1, ex.toy, ex.base

S = 8.0; J_TOT = 2000; J_WIN = 40

def main() -> None:
    t0 = time.time()
    params = base.Params(high_n=63, recursive_collision_count=200)
    sp = base.build_source_params(params)
    n, ne = sp.chi_grid_n, sp.eta_grid_n
    shape = (n, ne)
    ms = np.arange(ne); mm = np.where(ms <= ne // 2, ms, ms - ne)
    fold3 = ((mm + 1) % 3) - 1
    hid = (mm - fold3) // 3
    eta = 2 * np.pi * np.arange(ne) / ne

    def single_winding(v):
        f = np.fft.fft(v.reshape(shape), axis=0, norm="ortho")
        f[0, :] = 0.0; f[n // 2:, :] = 0.0
        return np.fft.ifft(f, axis=0, norm="ortho").reshape(v.shape)

    def project_eta(v, m_set):
        f = np.fft.fft(v.reshape(shape), axis=1, norm="ortho")
        keep = np.isin(mm, list(m_set)); f[:, ~keep] = 0.0
        return np.fft.ifft(f, axis=1, norm="ortho").reshape(v.shape)

    def shift_eta(v, dm):
        return (v.reshape(shape) * np.exp(1j * dm * eta)[None, :]).reshape(v.shape)

    a0 = single_winding(v1.make_bundle(sp, (30, 32, 34), "A", scale=1.0)) * S
    b0 = single_winding(v1.make_bundle(sp, (30, 32, 34), "B", scale=1.0)) * S
    a0 = a0 + single_winding(v1.make_bundle(sp, (21,), "A", scale=1.0)) * (0.2 * S)
    pow0 = float(np.sum(np.abs(a0) ** 2) + np.sum(np.abs(b0) ** 2))
    a1 = project_eta(a0, {1}); b1 = project_eta(b0, {1})
    pw = float(np.sum(np.abs(a1) ** 2) + np.sum(np.abs(b1) ** 2))
    sc = np.sqrt(pow0 / pw); a1 *= sc; b1 *= sc

    def make_sea(sw: bool):
        sa = v1.make_bundle(sp, (29, 31, 33), "A", scale=1.0) * S
        sb = v1.make_bundle(sp, (29, 31, 33), "B", scale=1.0) * S
        if sw:
            sa = single_winding(sa); sb = single_winding(sb)
        sa = project_eta(sa, {0}); sb = project_eta(sb, {0})
        pws = float(np.sum(np.abs(sa) ** 2) + np.sum(np.abs(sb) ** 2))
        scs = np.sqrt(0.25 * pow0 / pws)
        return sa * scs, sb * scs

    sea_sw = make_sea(True)
    sea_raw = make_sea(False)
    cases = {
        "A_+1+seaSW": (a1 + sea_sw[0], b1 + sea_sw[1]),
        "B_+1+seaRAW": (a1 + sea_raw[0], b1 + sea_raw[1]),
        "C_+2+seaSW": (shift_eta(a1, 1) + sea_sw[0], shift_eta(b1, 1) + sea_sw[1]),
        "D_v3orig": (a0.copy(), b0.copy()),
    }
    out = {"J_TOT": J_TOT, "J_WIN": J_WIN, "cases": {}}
    for name, (a, b) in cases.items():
        Qw, Q3s, Ws, ret = [], [], [], []
        # 支配巻き数（初期・全k）
        fa = np.fft.fft(np.fft.fft(a.reshape(shape), axis=0, norm="ortho"), axis=1, norm="ortho")
        fb = np.fft.fft(np.fft.fft(b.reshape(shape), axis=0, norm="ortho"), axis=1, norm="ortho")
        P0m = np.sum((np.abs(fa) ** 2 + np.abs(fb) ** 2) / 2, axis=0)
        chg = [m_ for m_ in range(-6, 7) if m_ != 0]
        m_dom = max(chg, key=lambda m_: P0m[list(mm).index(m_)] if m_ in list(mm) else 0)
        idx_dom = int(np.where(mm == m_dom)[0][0])
        Pwin = np.zeros(shape)
        for j in range(J_TOT):
            a, b, _ = ex.collision_step_exact(a, b, sp)
            fa = np.fft.fft(np.fft.fft(a.reshape(shape), axis=0, norm="ortho"), axis=1, norm="ortho")
            fb = np.fft.fft(np.fft.fft(b.reshape(shape), axis=0, norm="ortho"), axis=1, norm="ortho")
            Pwin += (np.abs(fa) ** 2 + np.abs(fb) ** 2) / 2
            if (j + 1) % J_WIN == 0:
                P = Pwin / J_WIN
                Pm = np.sum(P, axis=0)
                Qw.append(float(np.sum(mm * Pm)))
                Q3s.append(float(np.sum(fold3 * Pm)))
                Ws.append(float(np.sum(hid * Pm)))
                ret.append(float(Pm[idx_dom]))
                Pwin = np.zeros(shape)
        Qw = np.array(Qw); Q3s = np.array(Q3s); Ws = np.array(Ws); ret = np.array(ret)
        cvQ = float(np.std(Qw) / max(abs(np.mean(Qw)), 1e-300))
        # 簿記恒等式: Q_wind − (Q3+3W) は恒等的に零のはず（表現の恒等式）
        ledger_res = float(np.max(np.abs(Qw - (Q3s + 3 * Ws))))
        # ΔQ3=−3ΔW の検定（Q_wind保存が成立する場合のみ意味を持つ）
        dQ3 = Q3s[-1] - Q3s[0]; dW = Ws[-1] - Ws[0]
        closure = float(abs(dQ3 + 3 * dW) / max(abs(dQ3), 1e-300))
        tau = None
        if (ret > 0).all():
            coef = np.polyfit(np.arange(len(ret), dtype=float) * J_WIN, np.log(ret), 1)
            tau = float(-1.0 / coef[0]) if coef[0] < 0 else float("inf")
        print(f"{name}: 支配m={m_dom:+d}  Q_wind CV={cvQ:.5f}  "
              f"Q3 {Q3s[0]:+.3f}→{Q3s[-1]:+.3f}  W {Ws[0]:+.3f}→{Ws[-1]:+.3f}")
        print(f"    恒等式残差max={ledger_res:.2e}  ΔQ3+3ΔW閉じ={closure:.2e}  "
              f"支配種τ={tau if tau and tau != float('inf') else '∞'}")
        out["cases"][name] = {"m_dom": m_dom, "Q_wind_cv": cvQ,
            "Q3_first": float(Q3s[0]), "Q3_last": float(Q3s[-1]),
            "W_first": float(Ws[0]), "W_last": float(Ws[-1]),
            "ledger_residual_max": ledger_res, "delta_closure": closure,
            "tau_dom": tau,
            "H_conserve": bool(cvQ < 0.02)}
    a_ = out["cases"]
    h_conv = bool(a_["A_+1+seaSW"]["H_conserve"] and not a_["B_+1+seaRAW"]["H_conserve"])
    print(f"\nH_conv（規約: 全成分χ解析射影⇒保存 / raw海⇒破れ）= {h_conv}")
    out["H_conv"] = h_conv
    if a_["A_+1+seaSW"]["H_conserve"] and a_["C_+2+seaSW"]["H_conserve"]:
        tA = a_["A_+1+seaSW"]["tau_dom"] or float("inf")
        tC = a_["C_+2+seaSW"]["tau_dom"] or float("inf")
        print(f"正当計器での寿命比較: τ(+1)={tA:.0f} vs τ(+2)={tC:.0f} → "
              f"{'+1が長寿命' if tA > tC else '+2が長寿命または同等'}")
        out["tau_pm1_vs_p2"] = [tA if tA != float('inf') else None,
                                  tC if tC != float('inf') else None]
    out["runtime_sec"] = time.time() - t0
    (HERE / "pre_sea_convention_ledger_result_v7.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False))
    print(f"完了 {out['runtime_sec']:.0f}s → pre_sea_convention_ledger_result_v7.json")

if __name__ == "__main__":
    main()
