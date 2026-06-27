#!/usr/bin/env python3
"""
干渉縞パターン生成スクリプト
振動区間 N = 波列の周期数 (有限コヒーレント波列長 N λ = N)
2つのコヒーレント光子、波長 λ=1, c=1, 光源間隔 d=W=5
L十分大きい遠方場
ドップラー効果: 2光子が同一速度 c=1 で伝播するため相対ドップラーシフトなし（共通モード）
量子効果: 検出確率密度 ∝ I(Δs)（コヒーレント状態の極限）
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from scipy.signal import find_peaks
import os

# 日本語フォント設定
rcParams['font.family'] = 'Noto Serif CJK JP'
rcParams['axes.unicode_minus'] = False

OUTPUT_DIR = "/home/workdir/artifacts"
os.makedirs(OUTPUT_DIR, exist_ok=True)

def compute_I(ds, N, d=5.0):
    """有限波列干渉の強度"""
    # 物理的に可能な光路差 |Δs| <= d
    # 重なりが存在する |Δs| < N
    support = min(N, d)
    overlap = np.maximum(0.0, N - np.abs(ds))
    I = overlap * (1.0 + np.cos(2.0 * np.pi * ds))
    # 物理的範囲外は0
    I = np.where(np.abs(ds) > d, 0.0, I)
    return I, support

def count_bright_fringes(N, d=5.0):
    """明るい干渉縞の本数を厳密カウント
    明るい最大値位置: Δs = k (k整数) かつ |k| < min(N, d) かつ overlap >0
    """
    support = min(N, d)
    # k の範囲: |k| < support
    k_min = -int(np.floor(support - 1e-9))
    k_max = int(np.floor(support - 1e-9))
    if k_max < k_min:
        return 0, []
    ks = list(range(k_min, k_max + 1))
    n = len(ks)
    return n, ks

def generate_plot(N, save_png=True, save_svg=True):
    """各 N の干渉縞パターンを図化"""
    d = 5.0
    n_fringes, ks = count_bright_fringes(N, d)
    
    support = min(N, d)
    ds = np.linspace(-support - 0.5, support + 0.5, 8000)
    I, _ = compute_I(ds, N, d)
    
    fig, ax = plt.subplots(figsize=(12, 7))
    
    ax.plot(ds, I, 'b-', linewidth=2.0, label='干渉強度 I(Δs)')
    
    # 明るい縞位置をマーク
    for k in ks:
        overlap_k = max(0.0, N - abs(k))
        I_k = overlap_k * (1 + np.cos(2 * np.pi * k))  # = 2 * overlap_k
        ax.axvline(x=k, color='red', linestyle='--', alpha=0.6, linewidth=1.2)
        ax.plot(k, I_k, 'ro', markersize=8, zorder=5)
    
    # タイトルと注記
    title = f'N={N} の干渉縞パターン\n明るい干渉縞の本数 n = {n_fringes}  (N(N-1)/2 = {N*(N-1)//2})'
    ax.set_title(title, fontsize=14, pad=20)
    
    ax.set_xlabel('光路差 Δs (λ = 1 単位)', fontsize=12)
    ax.set_ylabel('強度 I(Δs)  [任意単位]', fontsize=12)
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.set_xlim(-support - 0.5, support + 0.5)
    
    # 注記ボックス
    note_text = (
        "モデル: 有限コヒーレント波列長 = N (周期数)\n"
        "2光子コヒーレント干渉（Z = cosθ + i sinθ 表現）\n"
        "ドップラー効果: 相対シフトなし（同一 c=1 伝播）\n"
        "量子効果: 光子検出確率密度 ∝ I(Δs)\n"
        f"有効サポート = min(N, W=5) = {support}\n"
        "明るい縞位置: Δs = k (k整数, |k| < support)"
    )
    ax.text(0.02, 0.98, note_text, transform=ax.transAxes, fontsize=9,
            verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    ax.legend(loc='upper right')
    
    base = os.path.join(OUTPUT_DIR, f"interference_fringes_N{N:02d}")
    if save_png:
        plt.savefig(base + ".png", dpi=200, bbox_inches='tight')
    if save_svg:
        plt.savefig(base + ".svg", bbox_inches='tight')
    plt.close(fig)
    
    return n_fringes

# === メイン実行 ===
print("干渉縞図の生成を開始します...")

results = []
for N in range(1, 11):
    n = generate_plot(N)
    comb = N * (N - 1) // 2
    results.append((N, comb, n))
    print(f"N={N:2d} | 組み合わせ N(N-1)/2 = {comb:3d} | 実際の明るい干渉縞の本数 n = {n:2d}")

# Markdown テーブル作成
md_content = """# 有限波列干渉縞 本数検証結果 (N=1〜10)

## モデル概要
- 単一振動数、λ = 1 正規化、c = 1
- W = 5（光源間隔 d = 5）
- 各光子: 有限コヒーレント波列（振動区間 N 周期、長さ N）
- 2光子コヒーレント干渉（複素 Z = cos ϕ + i sin ϕ 表現）
- L 十分大きい遠方場
- **ドップラー効果**: 2光子が同一速度 c=1 で伝播するため**相対ドップラーシフトなし**（共通モードのみ）
- **量子効果**: 検出確率密度は古典強度 I(Δs) に比例（コヒーレント極限）

## 詳細計算式

光路差 Δs における重なり因子（波列長 N）:
$$
overlap(Δs) = \\max(0, N - |Δs|)
$$

干渉強度（時間平均）:
$$
I(Δs) = overlap(Δs) \\times (1 + \\cos(2\\pi \\Delta s / \\lambda)) \\quad (\\lambda=1)
$$

明るい干渉縞の位置（最大値）:
$$
Δs = k \\quad (k \\in \\mathbb{Z},\\ |k| < \\min(N, W=5))
$$

**干渉縞の本数 n**:
$$
n = 2 \\lfloor \\min(N, 5) \\rfloor - 1 \\quad (N \\text{ が整数の場合、厳密には } 2N-1 \\text{ または幾何制限})
$$

（k = −(m−1) から +(m−1) まで、m = floor(min(N,5))）

## 検証結果一覧表

| N | N(N-1)/2 (組み合わせ数) | 実際の明るい干渉縞の本数 n | 備考 |
|---|--------------------------|-----------------------------|------|
"""

for N, comb, n in results:
    note = ""
    if N <= 5:
        note = f"波列長 N が支配 (n = 2N-1 = {2*N-1})"
    else:
        note = f"W=5 の幾何制限が支配 (n=9)"
    md_content += f"| {N} | {comb} | {n} | {note} |\n"

md_content += """
## 結論
N=1 から 10 まで、**n = N(N-1)/2 には一致しません**。
実際の干渉縞の本数は波列長 N（または W=5 の幾何制限）で決まる**線形の量**です。
N(N-1)/2 は「ペアの総数」ですが、同じ光路差 k に複数のペアが寄与するため、観測される空間的な縞の本数には集約されます。

生成された図（PNG / SVG）は artifacts/ ディレクトリに保存されています。
各図には明るい縞位置（赤点・破線）が明示されています。
"""

with open(os.path.join(OUTPUT_DIR, "interference_fringes_summary.md"), "w", encoding="utf-8") as f:
    f.write(md_content)

print("\n=== 完了 ===")
print("生成ファイル:")
print("  - interference_fringes_N01.png ～ N10.png")
print("  - interference_fringes_N01.svg ～ N10.svg")
print("  - interference_fringes_summary.md")
print("  - このスクリプト: generate_interference_figures.py")
print(f"\nすべてのファイルは {OUTPUT_DIR} にあります。")