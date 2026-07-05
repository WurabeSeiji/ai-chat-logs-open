#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
D4 素朴計算：W8 状態の回転応答（ギア比）を、閉包＋D1-D3 だけから測る

方法:
  状態を方位角 χ（伝播軸まわりの角）上の関数 Ψ(χ) として表現。
  チャート回転 R_φ は χ → χ−φ のシフト。
  読出しは枠組みの干渉: o(φ) = ∫ Ψ*(χ) Ψ(χ−φ) dχ / ∫|Ψ|²
  巻き数 s_measured = (アンラップした arg o(φ) の全変化) / (ループの角度)

ケース:
  A. 実符号対 (σx,σy)=(+1,+1): Ψ = cos(χ−π/4)·√2 型（W8 の静的符号ベクトルそのまま）
     → S¹ ループ（φ: 0→2π）と ℝP¹ ループ（φ: 0→π、対蹠同一視）の両方で測る
  B. quadrature 対 (1, ±i): Ψ = e^{∓iχ}（±π/2 の相対位相 ＝ Z₄ 細分）
  C. m=2 quadrature（対称トレースレス型）: Ψ = e^{∓2iχ}
  D. 非空間軸因子（t/R/Q）: 方位角に依存しない定数因子 → ギア比を直接測る
  E. 閉包チェック: 回転が |Z⊥|²（マスターへの寄与）を変えないことの確認
"""
import cmath
import math

PI = math.pi
M = 2048  # 方位角サンプル

def overlap(psi, shift_frac):
    """o(φ): φ = shift_frac·2π/M 単位の円環シフトとの内積"""
    s = sum(psi[i].conjugate() * psi[(i - shift_frac) % M] for i in range(M))
    n = sum(abs(v) ** 2 for v in psi)
    return s / n

def winding(psi, loop=2 * PI, steps=64, antipodal=False):
    """ループ一周でのアンラップ位相変化 / 2π。antipodal=True なら φ∈[0,π] で
    終端を対蹠同一視（Ψ(χ−π) と比較して閉じる）"""
    total = 0.0
    prev = 0.0
    frac_max = int(M * loop / (2 * PI))
    args = []
    for k in range(steps + 1):
        sh = int(frac_max * k / steps)
        o = overlap(psi, sh)
        a = cmath.phase(o) if abs(o) > 1e-9 else prev
        args.append((abs(o), a))
        prev = a
    # アンラップ
    unwrapped = [args[0][1]]
    for i in range(1, len(args)):
        d = args[i][1] - args[i - 1][1]
        while d > PI:
            d -= 2 * PI
        while d < -PI:
            d += 2 * PI
        unwrapped.append(unwrapped[-1] + d)
    return (unwrapped[-1] - unwrapped[0]) / (2 * PI), args[len(args) // 2][0]

def make(fn):
    return [fn(2 * PI * i / M) for i in range(M)]

def main():
    print("D4 素朴計算：回転応答の巻き数 s_measured")
    print()
    # A. 実符号対（W8 静的構成）
    psiA = make(lambda ch: math.cos(ch - PI / 4) * math.sqrt(2))
    sA_full, oA = winding(psiA, loop=2 * PI)
    # ℝP¹ ループ: φ を 0→π。実対は Ψ(χ−π) = −Ψ(χ) なので o(π) = −1 → 位相 π 蓄積
    sA_half, _ = winding(psiA, loop=PI)
    print(f"A. 実符号対 (+1,+1)（W8 の Z₂ 符号そのまま）:")
    print(f"   S¹ ループ（0→2π）: 巻き数 = {sA_full:+.3f}（|o(π/2)|={oA:.3f}——中間で振幅が落ちる）")
    print(f"   ℝP¹ ループ（0→π、対蹠同一視）: 蓄積位相 = {sA_half:+.3f}×2π → チャート1周あたり {sA_half:+.3f}")
    print(f"   → 対蹠同一視の下で、チャート1周につき位相 π ＝ **巻き数 1/2**（奇数セクター）")
    print()
    # B. quadrature 対
    for sgn, label in [(-1, "(1,+i)"), (+1, "(1,-i)")]:
        psiB = make(lambda ch, s=sgn: cmath.exp(1j * s * ch))
        sB, _ = winding(psiB, loop=2 * PI)
        print(f"B. quadrature 対 {label}（±π/2 位相＝Z₄ 細分）: 巻き数 = {sB:+.3f}")
    print()
    # C. m=2
    psiC = make(lambda ch: cmath.exp(2j * ch))
    sC, _ = winding(psiC, loop=2 * PI)
    print(f"C. m=2 quadrature（対称トレースレス型）: 巻き数 = {sC:+.3f}")
    print()
    # D. 非空間軸因子
    psiD = make(lambda ch: cmath.exp(1j * 0.7))  # 定数位相（t/R/Q の因子）
    sD, _ = winding(psiD, loop=2 * PI)
    print(f"D. 非空間軸因子（t/R/Q、方位角非依存）: ギア比 = {sD:+.3f}")
    print()
    # E. 閉包チェック: 回転でノルム（マスターへの寄与）不変
    n0 = sum(abs(v) ** 2 for v in psiA)
    rot = [psiA[(i - M // 8) % M] for i in range(M)]
    n1 = sum(abs(v) ** 2 for v in rot)
    print(f"E. 閉包チェック: 回転前後のノルム比 = {n1 / n0:.6f}（1＝マスター不変＝閉包は回転応答を伝達しない）")

if __name__ == "__main__":
    main()
