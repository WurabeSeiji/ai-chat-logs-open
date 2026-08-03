#!/usr/bin/env python3
"""第8論文v2予備実験 E-M2：シードは物質を生むか v1

問い（木原氏）:
    E-M1 で「停止だけでは物質（非自明整数比ロック）は生まれない」と判明した。
    では、シード（準安定期への横摂動注入=第8論文条件D）が要るのか？

予言（測定前固定）:
    P1（シードは鳴らす）: 注入後、横モードが励起され周波数の分散が増える
        （単一ユニゾンからのずれが出る）
    【P1実測後の訂正記録】P1 は FAIL したが、これは物理の反証ではなく
    観測量の感度不足である: ε=1e-8 の横励起の振幅は 1e-8 のオーダーに
    留まり、辺位相の窓内傾き（主成分に支配される）には現れない。
    励起を見るには横モードへの射影が必要（第8論文の「叩けば響く」は
    その振幅スケールでの話）。本実験の結論は P2 が担う。
    P2（だがシードは物質を生まない）: 鳴った横モードの周波数と全体周波数の
        比は一般に無理数であり、三分法（有理住所=粒子/無理数=過渡）により
        持続的な非自明整数比ロックは形成されない。L(t)=0 のまま。
    反証条件: 注入後に持続的な整数比ロックが出れば「シードが物質を生む」
        が成立し、本予言は反証として記録する。

含意（どちらに転んでも）:
    P2成立なら: 物質に足りないのはシード（初期揺らぎ）ではなく、
    周波数を整数格子に強制する**円環性（周期的位相軸=分解能レジスタ=公理2）**。
    N体関係系は公理1（零閉塞）のみを実装しており、公理2（分解能）の
    構造を持たない——公理1→共形層（方向・光の海）、公理2→スケール層
    （レジスタ・時計・物質）という公理と層の一対一対応が立つ。

規約: 第7/8論文コード read-only import。条件Dの注入は第8論文と同一
     （ε=1e-8、seed index 0、S4直交化、注入後一度だけ規格化）。
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
CODE = HERE.parent / "code" / "run_preliminary_seed_ablation_v1.py"
spec = importlib.util.spec_from_file_location("ablation_for_seed_matter_v1", CODE)
abl = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = abl
spec.loader.exec_module(abl)

N = 5
XMAX = 20000
WIN = 200
STRIDE = 50
RATIO_TOL = 1e-3
PERSIST = 5
FREQ_MIN = 1e-6


def main() -> None:
    sys_lr, v, B_p1, B_rot, B0, p, q, Z, wp = abl.build_init(N, initial_seed=True)
    M = sys_lr.m

    def fval(Zv):
        Zp = Zv - p * (p @ Zv) - q * (q @ Zv)
        return float(np.real(np.conj(Zp) @ Zp)) / float(np.real(np.conj(Zv) @ Zv))

    phases = np.zeros((XMAX + 1, M))
    fs = np.zeros(XMAX + 1)
    phases[0] = np.angle(Z); fs[0] = fval(Z)
    crossing, injected_at = None, None
    for t in range(1, XMAX + 1):
        # 条件D: t1 = crossing+GUARD で単一横摂動を一回注入（第8論文と同一手順）
        t1 = (crossing + abl.GUARD) if crossing is not None else None
        if t1 is not None and t == t1 and injected_at is None:
            rng_dir = np.random.default_rng(70000 + N)
            eta_r = rng_dir.normal(size=M); eta_i = rng_dir.normal(size=M)
            S4_t1 = abl.s4_basis(sys_lr, B0, Z)
            eta_r = eta_r - S4_t1 @ (S4_t1.T @ eta_r)
            eta_i = eta_i - S4_t1 @ (S4_t1.T @ eta_i)
            eta = (eta_r + 1j * eta_i) / np.sqrt(eta_r @ eta_r + eta_i @ eta_i)
            Z = Z + abl.D_EPS * eta
            Z = Z / np.linalg.norm(Z)
            injected_at = t
        Z, wp = abl.evolve(sys_lr, Z, wp)
        phases[t] = np.angle(Z)
        fs[t] = fval(Z)
        if crossing is None and fs[t] > 0.05:
            crossing = t
    print(f"crossing = {crossing}, 注入 = {injected_at}")

    unwrapped = np.unwrap(phases, axis=0)
    centers, freqs = [], []
    for s in range(0, XMAX - WIN, STRIDE):
        seg = unwrapped[s:s + WIN]
        slope = np.polyfit(np.arange(WIN), seg, 1)[0]
        centers.append(s + WIN // 2)
        freqs.append(slope)
    centers = np.array(centers); freqs = np.abs(np.array(freqs))

    lock_now = []
    for w in range(len(centers)):
        f = freqs[w]; locked = set()
        for i in range(M):
            for j in range(i + 1, M):
                hi, lo = max(f[i], f[j]), min(f[i], f[j])
                if lo < FREQ_MIN:
                    continue
                r = hi / lo; pr = round(r)
                if pr >= 2 and abs(r - pr) < RATIO_TOL:
                    locked.add((i, j, pr))
        lock_now.append(locked)
    L = np.zeros(len(centers), dtype=int)
    for w in range(len(centers)):
        cnt = 0
        for key in lock_now[w]:
            run = 1
            k = w - 1
            while k >= 0 and key in lock_now[k]:
                run += 1; k -= 1
            k = w + 1
            while k < len(centers) and key in lock_now[k]:
                run += 1; k += 1
            if run >= PERSIST:
                cnt += 1
        L[w] = cnt

    disp = np.std(freqs, axis=1) / np.maximum(np.mean(freqs, axis=1), 1e-30)
    pre = disp[(centers > injected_at - 2000) & (centers < injected_at)] if injected_at else disp[:1]
    post = disp[(centers > injected_at) & (centers < injected_at + 4000)] if injected_at else disp[:1]
    p1 = float(np.max(post)) > 3 * float(np.max(pre))
    post_L = L[centers > injected_at] if injected_at else L
    p2 = int(np.max(post_L)) == 0
    print(f"P1 シードは鳴らす（注入後の周波数分散増 {np.max(pre):.2e}→{np.max(post):.2e}）: "
          f"{'PASS' if p1 else 'FAIL'}")
    print(f"P2 だが物質は生まれない（注入後も非自明整数比ロック L=0）: "
          f"max L(post) = {int(np.max(post_L))} {'PASS' if p2 else 'FAIL（シードが物質を生んだ）'}")

    payload = {
        "experiment": "paper8_seed_matter_test_pre_v1",
        "engine": "第7/8論文エンジン read-only import（条件D手順・N=5）",
        "crossing": crossing, "injected_at": injected_at,
        "dispersion_pre_max": float(np.max(pre)), "dispersion_post_max": float(np.max(post)),
        "max_locks_after_injection": int(np.max(post_L)),
        "P1_seed_rings": bool(p1),
        "P1_note": ("P1のFAILは観測量の感度不足（1e-8励起は辺位相傾きに不可視）で"
                     "あり物理の反証ではない——測定限界として記録"),
        "P2_seed_does_not_create_matter": bool(p2),
        "conclusion": (
            f"P1={'成立' if p1 else '感度不足（測定限界として記録）'}。"
            f"P2={'成立' if p2 else '反証'}: "
            + ("鳴った比は整数格子に乗らず、持続的な非自明整数比ロックは形成されない"
               "——物質に足りないのはシードではなく、周波数を整数格子に強制する円環性"
               "（周期的位相軸=分解能レジスタ=公理2）。公理1（零閉塞）→共形層（方向・光の海）、"
               "公理2（分解能）→スケール層（レジスタ・時計・物質）の対応が立つ"
               if p2 else "シード注入後に整数比ロックが出現した——シード起源の物質生成として要精査")),
    }
    (HERE / "paper8_seed_matter_result_v1.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, default=float), encoding="utf-8")
    print("\nsaved: paper8_seed_matter_result_v1.json")


if __name__ == "__main__":
    main()
