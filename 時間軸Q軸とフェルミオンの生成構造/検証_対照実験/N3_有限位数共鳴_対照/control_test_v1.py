#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
対照テスト: コピーした「有限位数共鳴」コード（2体交換散乱）が公開結果を再現することの確認

N3 が依拠する 2体交換散乱の解析核（散乱行列・固有値・有限位数根 R=cos²(πm/n)）を、
コピーした phase5_eigenphase_resonance_v2.py の関数で再現し、公開論文の基準値と照合する。

公開基準値（finite_order_resonance_v1 論文より）:
  R_{124,23} = cos²(23π/124) = 0.697177927556659...
  R_{122,23} = cos²(23π/122) = 0.688363946817593...
  各共鳴 R で U^n が単位行列へ厳密復帰（matrix_return_error ≈ 0）。

【スコープ注記】本テストは根構造のみを扱う。公開バンドルに含まれる α⁻¹（N_of_R=137.03…）
同定は N3 では扱わず、本テストでも参照しない。

実行: python3 control_test_v1.py
"""

import importlib.util
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "finite_order_resonance_v1", "src",
                   "phase5_eigenphase_resonance_v2.py")
spec = importlib.util.spec_from_file_location("phase5_verified", SRC)
p5 = importlib.util.module_from_spec(spec)
sys.modules["phase5_verified"] = p5  # dataclass 解決のため登録
spec.loader.exec_module(p5)

# 公開論文の基準値（根の解析値）
PUBLISHED = [
    ("R_{124,23}", 124, 23, 0.697177927556659),
    ("R_{122,23}", 122, 23, 0.688363946817593),
]


def main():
    print("対照テスト: 2体交換散乱の解析核（phase5 の関数を import）")
    max_root_err = 0.0
    max_return_err = 0.0
    max_sym_err = 0.0
    for label, n, m, pub_R in PUBLISHED:
        R = p5.resonant_R(n, m)                 # cos²(πm/n)
        max_root_err = max(max_root_err, abs(R - pub_R))
        lam_s, lam_a = p5.exact_eigenvalues(R)  # (1, -e^{2iθ})
        max_sym_err = max(max_sym_err, abs(lam_s - 1.0))
        ret = p5.matrix_return_error(R, n)      # |U^n - I|
        max_return_err = max(max_return_err, ret)
        print(f"  {label}: resonant_R={R:.15f}  公開値={pub_R:.15f}  "
              f"|U^{n}-I|={ret:.2e}")
    print(f"根値の公開基準との最大差:   {max_root_err:.2e}（0 が予言）")
    print(f"対称固有値 λ_sym−1 の最大:  {max_sym_err:.2e}（0 が予言）")
    print(f"U^n 復帰誤差の最大:         {max_return_err:.2e}（機械精度で0）")
    ok = max_root_err < 1e-14 and max_sym_err < 1e-14 and max_return_err < 1e-9
    print("対照テスト:", "PASS（公開結果を再現）" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
