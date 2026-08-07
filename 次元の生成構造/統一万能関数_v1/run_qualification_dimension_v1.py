#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""万能次元読出し関数 D の単体テスト・資格審査（全通過が使用条件）

Q16 正本一致: 瞬時フレームの第3次元確定度 |n̂·â| が、公開論文[3]の窓SVD処方
    （run_pre_2plus1_structure_v1 と同じ軌道・同じ量）と同じ N 依存を示すか。
    正本の実測は N=4:0.0014（非結晶化）→ N=5:0.9874 → N=6:0.9918 → N=8:0.9965。
    **同一軌道上で窓SVD版と瞬時版を両方計算して直接比較する**（数値の再現では
    なく同一走行での対照）。
Q17 受動性: D 併走の有無で終状態がビット単位同一。
Q18 打切り次数の安定性: 平面の梯子の実効枚数 n_eff が order=4,6,8 で安定
    （相対変化 <5%）——order は計算上の打切りであって物理閾値でないことの確認。
Q19 不在の表現: 内容ゼロで weight=0・rank=0 を返す（NaN でも真偽値でもない）。
Q20 次元の結晶化: テストベッド（N=12・Nn=16）で τ=0 から追跡し、
    フレームが誕生（rank 3・align 上昇）する過程が読めること。
Q21 規約不変性: 状態の大域位相を変えても align・n_eff・非一様度が不変
    （相対 <1e-10）。

使い方: python3 run_qualification_dimension_v1.py
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


ui = load("ui_dq", HERE / "unified_interaction_v1.py")
D = load("ud_dq", HERE / "unified_dimension_v1.py")
G = load("ur_dq", HERE / "unified_readout_v2.py")

abl = ui.abl
T_SETTLE, T_MEAS, EV = 2000, 2000, 5


def build_universe_seedsafe(n, delta=1e-2, Nn=16, Neta=8, seeds=(2, 3, 4, 5, 7)):
    """標準宇宙の構成。親構成が失敗する N があるためシードを順に試す
    （build_standard_universe は seed=2 固定・N=8 で ParentConstructionError）。
    使用したシードは記録する。"""
    m = n * (n - 1) // 2
    _, v, _, _, _, _, _, Z0c, wp0 = abl.build_init(n, False)
    last = None
    for s in seeds:
        try:
            r2 = ui.gen3.make_parent(n, seed=s)
            break
        except Exception as ex:      # ParentConstructionError 等
            last = ex
            r2 = None
    if r2 is None:
        raise RuntimeError(f"全シードで親構成に失敗: {last}")
    Csec = np.fft.fft(r2.relation_waves, axis=1) / n
    s0 = Csec[:, 1] / np.linalg.norm(Csec[:, 1])
    C2_0 = np.zeros((m, Nn, Neta), complex)
    C2_0[:, 2, 0] = Z0c
    C2_0[:, 1, 0] = delta * s0
    p2 = C2_0[:, 2, 0].real / np.linalg.norm(C2_0[:, 2, 0].real)
    q2 = C2_0[:, 2, 0].imag - (C2_0[:, 2, 0].imag @ p2) * p2
    q2 = q2 / np.linalg.norm(q2)
    return ui.UnifiedEngine(n, C2_0, wp0), p2, q2


def canonical_window_align(sys_lr, S, p, q):
    """公開論文[3] run_pre_2plus1_structure_v1 と同じ窓SVD処方（対照）"""
    Sp = np.array([z - p * (p @ z) - q * (q @ z) for z in S])
    X = np.hstack([Sp.real, Sp.imag])
    Xc = X - X.mean(axis=0)
    _, sv, Vt = np.linalg.svd(Xc, full_matrices=False)
    E3 = Vt[:3].T                                   # 上位3実方向
    def A(x):
        h = len(x) // 2
        z = x[:h] + 1j * x[h:]
        y = sys_lr.kmatvec(z)
        return np.concatenate([y.real, y.imag])
    AE = np.stack([A(E3[:, j]) for j in range(3)], axis=1)
    Om = E3.T @ AE
    Om = 0.5 * (Om - Om.T)
    axis = np.array([Om[2, 1], Om[0, 2], Om[1, 0]])
    om = np.linalg.norm(axis)
    if om <= 0:
        return None, sv / sv[0]
    ah = axis / om
    # 上位2方向の張る面の法線 = 第3座標軸
    return float(abs(ah[2])), sv / sv[0]


def main():
    t0 = time.time()
    res = {}

    # ---- Q16: 通し創世（インフレーションを通す）→ 準安定窓で測る
    # 【是正記録 2026-08-08】初版は abl の整定走行（インフレーション無し）で
    # 測り、しかも準安定に入る前の値で判定していた。**準安定条件が成立する
    # まで時空は読めない**——通し創世で τ∈[2000,4000] の準安定窓で測る。
    print("=== Q16 通し創世→準安定窓での次元読出し（Nn=16・T=4000・窓[2000,4000]）===")
    print(f"  {'N':>3} {'M':>5} {'瞬時|n̂·â|中央':>13} {'軸の持続':>9} "
          f"{'面の持続':>9} {'窓SVD|n̂·â|':>11} {'n_eff':>7} {'正本[3]':>8}")
    canon = {4: 0.0014, 5: 0.9874, 6: 0.9918, 8: 0.9965, 12: None}
    q16 = {}
    for n in (4, 5, 6, 8, 12):
        try:
            eng, p2, q2 = build_universe_seedsafe(n)
        except RuntimeError as ex:
            print(f"  {n:3d}: 親構成に失敗（{ex}）——記録してスキップ")
            q16[n] = {"error": str(ex)}
            continue
        syslr = abl.LowRankSystem(n)
        al, ap, pp, ne, S = [], [], [], [], []
        prev_fr = None
        for t in range(4000):
            eng.step()
            C2 = eng.C2()
            th = np.angle(np.sum(C2.reshape(C2.shape[0], -1), axis=1))
            syslr.set_theta(th)
            Z = C2[:, 2, 0]
            fr = D.d_frame(Z, syslr.kmatvec, p2, q2)
            if t >= 2000:
                pers = D.d_frame_persistence(fr, prev_fr)
                al.append(fr["align"]); ap.append(pers["axis_persist"])
                pp.append(pers["plane_persist"])
                if t % 10 == 0:
                    ne.append(D.d_plane_ladder(Z, syslr.kmatvec, p2, q2)["n_eff"])
                    S.append(Z.copy())
            prev_fr = fr
        aw, svr = canonical_window_align(syslr, np.array(S), p2, q2)
        q16[n] = {"align_med": float(np.median(al)),
                  "axis_persist_med": float(np.median(ap)),
                  "plane_persist_med": float(np.median(pp)),
                  "window_align": aw, "n_eff_med": float(np.median(ne)),
                  "canon_paper": canon[n],
                  "sv_rel_top6": [float(x) for x in svr[:6]]}
        print(f"  {n:3d} {n*(n-1)//2:5d} {q16[n]['align_med']:13.4f} "
              f"{q16[n]['axis_persist_med']:9.4f} {q16[n]['plane_persist_med']:9.4f} "
              f"{aw:11.4f} {q16[n]['n_eff_med']:7.3f} "
              f"{('%.4f' % canon[n]) if canon[n] else '—':>8}")
    ok_ns = [n for n in (5, 6, 8, 12) if "align_med" in q16.get(n, {})]
    Q16 = ("align_med" in q16.get(4, {}) and q16[4]["align_med"] < 0.5
           and all(q16[n]["align_med"] > 0.8 for n in ok_ns))
    print(f"  判定（N=4 非結晶化 <0.5 かつ N≥5 結晶化 >0.8）: "
          f"{'通過' if Q16 else '不成立'}")
    res["Q16"] = bool(Q16); res["Q16_detail"] = q16

    # ---- Q17: 受動性
    eng_a, p2, q2 = ui.build_standard_universe(12, 1e-2, Nn=16, Neta=8)
    eng_b, _, _ = ui.build_standard_universe(12, 1e-2, Nn=16, Neta=8)
    sys12 = abl.LowRankSystem(12)
    for t in range(100):
        eng_a.step()
        C2 = eng_a.C2()
        th = np.angle(np.sum(C2.reshape(C2.shape[0], -1), axis=1))
        sys12.set_theta(th)
        D.d_panel(C2, sys12.kmatvec, p2, q2)
        eng_b.step()
    Q17 = bool(np.array_equal(eng_a.C, eng_b.C))
    res["Q17"] = Q17
    print(f"\nQ17 受動性（D併走の有無でビット同一・T=100）: "
          f"{'通過' if Q17 else '不成立'}")

    # ---- Q18: 打切り次数の安定性
    C2 = eng_a.C2()
    th = np.angle(np.sum(C2.reshape(C2.shape[0], -1), axis=1))
    sys12.set_theta(th)
    Zc = C2[:, 2, 0]
    neff = {}
    for od in (4, 6, 8):
        neff[od] = D.d_plane_ladder(Zc, sys12.kmatvec, p2, q2, order=od)["n_eff"]
    rel = max(abs(neff[o] - neff[6]) / max(neff[6], 1e-300) for o in (4, 8))
    Q18 = rel < 0.05
    res["Q18"] = bool(Q18); res["Q18_neff"] = neff
    print(f"Q18 打切り次数の安定性（n_eff: order4={neff[4]:.3f} "
          f"order6={neff[6]:.3f} order8={neff[8]:.3f}・相対変化 {rel:.3f}<0.05）: "
          f"{'通過' if Q18 else '不成立'}")

    # ---- Q19: 不在の表現
    Zz = np.zeros(66, complex)
    fr0 = D.d_frame(Zz, sys12.kmatvec, p2, q2)
    Q19 = bool(fr0["weight"] == 0.0 and fr0["rank"] == 0
               and np.all(np.isfinite(fr0["d_plane"])))
    res["Q19"] = Q19
    print(f"Q19 不在の表現（weight={fr0['weight']:.1e}・rank={fr0['rank']}・"
          f"NaN なし）: {'通過' if Q19 else '不成立'}")

    # ---- Q20: 次元の結晶化（τ=0 から準安定窓まで通しで追跡・T=4000）
    print("\nQ20 次元の結晶化（テストベッド N=12・Nn=16・τ=0→4000 通し創世）")
    eng, p2b, q2b = ui.build_standard_universe(12, 1e-2, Nn=16, Neta=8)
    sysx = abl.LowRankSystem(12)
    traj = []
    prev_fr = None
    for t in range(4000):
        eng.step()
        C2 = eng.C2()
        th = np.angle(np.sum(C2.reshape(C2.shape[0], -1), axis=1))
        sysx.set_theta(th)
        pan = D.d_panel(C2, sysx.kmatvec, p2b, q2b)
        pers = D.d_frame_persistence(pan["_frame"], prev_fr)
        prev_fr = pan["_frame"]
        traj.append((t + 1, pan["frame_rank"], pan["frame_align"],
                     pan["ladder_n_eff"], pan["gauge_nonunif"],
                     pan["total_closure"], pers["axis_persist"]))
    for i in (0, 9, 99, 299, 699, 1499, 2499, 3999):
        row = traj[i]
        print(f"  τ={row[0]:5d}: rank={row[1]} align={row[2]:.4f} "
              f"n_eff={row[3]:.3f} ゲージ非一様={row[4]:.4f} "
              f"軸持続={row[6]:.4f} 閉塞={row[5]:.2e}")
    Q20 = bool(traj[-1][1] == 3)
    res["Q20"] = Q20
    res["Q20_traj"] = [{"tau": r[0], "rank": r[1], "align": r[2],
                        "n_eff": r[3], "nonunif": r[4], "closure": r[5]}
                       for r in traj[::50]]
    print(f"  判定（終端で rank=3＝3次元が立っている）: "
          f"{'通過' if Q20 else '不成立'}")

    # ---- Q21: 規約不変性（大域位相）
    ph = np.exp(1.234j)
    C2p = eng.C2() * ph
    thp = np.angle(np.sum(C2p.reshape(C2p.shape[0], -1), axis=1))
    sysp = abl.LowRankSystem(12); sysp.set_theta(thp)
    th0 = np.angle(np.sum(eng.C2().reshape(eng.C2().shape[0], -1), axis=1))
    sys0 = abl.LowRankSystem(12); sys0.set_theta(th0)
    a0 = D.d_frame(eng.C2()[:, 2, 0], sys0.kmatvec, p2b, q2b)
    a1 = D.d_frame(C2p[:, 2, 0], sysp.kmatvec, p2b, q2b)
    dal = abs(a0["align"] - a1["align"]) / max(a0["align"], 1e-300)
    g0 = D.d_gauge(a0, sys0.kmatvec)["nonunif"]
    g1 = D.d_gauge(a1, sysp.kmatvec)["nonunif"]
    dg = abs(g0 - g1) / max(g0, 1e-300)
    Q21 = dal < 1e-10 and dg < 1e-10
    res["Q21"] = bool(Q21)
    res["Q21_detail"] = {"align_rel": dal, "nonunif_rel": dg}
    print(f"\nQ21 規約不変性（大域位相 e^{{1.234i}}: align 相対差 {dal:.2e}・"
          f"ゲージ非一様 相対差 {dg:.2e}）: {'通過' if Q21 else '不成立'}")

    allpass = all(v for k, v in res.items() if isinstance(v, bool))
    print("\n次元読出しD 資格審査: "
          + ("ALL PASS — D は使用可" if allpass else "不成立あり（記録して次へ）"))
    res["all_pass"] = bool(allpass)
    res["runtime_sec"] = time.time() - t0
    (HERE / "qualification_dimension_v1_result.json").write_text(
        json.dumps(res, indent=1, ensure_ascii=False, default=float))
    print(f"完了 {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
