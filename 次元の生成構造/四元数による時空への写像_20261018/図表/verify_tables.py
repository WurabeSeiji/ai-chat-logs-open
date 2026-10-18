# -*- coding: utf-8 -*-
"""四元数による状態空間から時空への写像：論文用の表の数値を計算して確かめ、results/ に書き出す。

表 1  閉塞の条件の本数（複素数・四元数の本当の二乗・等方的な条件）
表 2  N ごとの自由度と、空間として開ける方向の数（構成した解の数値で確かめる）
表 3  最も対称な解の一覧（Σq、Σq²、二次モーメント、q⁸）と、対照としての正十六胞体
表 4  写像の各段：全投射か部分投射か（数値で確かめる）
表 5  記法
"""
import itertools, json
from pathlib import Path
import numpy as np

OUT = Path(__file__).resolve().parent / "results"
OUT.mkdir(exist_ok=True)
rng = np.random.default_rng(20261018)


def qmul(a, b):
    a0, a1, a2, a3 = a
    b0, b1, b2, b3 = b
    return np.array([a0*b0 - a1*b1 - a2*b2 - a3*b3, a0*b1 + a1*b0 + a2*b3 - a3*b2,
                     a0*b2 - a1*b3 + a2*b0 + a3*b1, a0*b3 + a1*b2 - a2*b1 + a3*b0])


def qpow(q, k):
    p = np.array([1.0, 0, 0, 0])
    for _ in range(k):
        p = qmul(p, q)
    return p


def sums(Q):
    S = Q.sum(axis=0)
    S2 = np.sum([qmul(q, q) for q in Q], axis=0)
    G = Q.T @ Q
    return S, S2, G


def sum_zero_basis(N):
    """和が 0 の超平面（N−1 次元）の正規直交基底。"""
    M = np.eye(N) - np.ones((N, N))/N
    U, s, _ = np.linalg.svd(M)
    return U[:, :N - 1]


res = {}
md = ["# 論文用の表（verify_tables.py が計算して書き出したもの）", ""]

# ---------------------------------------------------------------- 表 1
md += ["## 表 1　閉塞の条件の本数", "",
       "| 基本量 | 線形和 | 二次（長さ） | 二次（交差項） | 計 | 内容 |", "|---|---|---|---|---|---|",
       "| 複素数 z＝a＋bi | Σa＝Σb＝0（2） | Σa²＝Σb²（1） | Σab＝0（1） | 4 | Σz＝0 かつ Σz²＝0 と同じ |",
       "| 四元数・本当の二乗 Σq＝0, Σq²＝0 | Σq₀＝…＝Σq₃＝0（4） | Σq₀²＝Σq₁²＋Σq₂²＋Σq₃²（1） | Σq₀q₁＝Σq₀q₂＝Σq₀q₃＝0（3） | 8 | 1＋3 に分かれる |",
       "| 四元数・等方的 | Σq₀＝…＝Σq₃＝0（4） | Σq₀²＝Σq₁²＝Σq₂²＝Σq₃²（3） | 6 個の交差項がすべて 0（6） | 13 | 4 成分が対等 |", "",
       "等方的な 13 本を満たせば Σq₀²＝s、Σ|v|²＝3s なので、Σq² の実部は −2s ≠ 0。二つの四元数の閉塞は両立しない（表 3 の正十六胞体で確かめる）。", ""]

# ---------------------------------------------------------------- 表 2
rows2 = []
for N in range(3, 11):
    B = sum_zero_basis(N)
    k = min(3, N - 2)
    # 本当の二乗の側：q₀ を和 0 の空間の 1 方向に、空間成分を残りの方向に置く
    W = B @ np.linalg.qr(rng.normal(size=(N - 1, N - 1)))[0]
    q0 = W[:, 0]
    sp = W[:, 1:1 + k] * np.sqrt(1.0 / k)
    Q = np.zeros((N, 4)); Q[:, 0] = q0; Q[:, 1:1 + k] = sp
    S, S2, G = sums(Q)
    rank_space = np.linalg.matrix_rank(Q[:, 1:], tol=1e-9)
    ok_true = np.linalg.norm(S) < 1e-12 and np.linalg.norm(S2) < 1e-12
    # 等方の側：和 0 の空間に正規直交な 4 本が取れるか
    iso = None
    if N - 1 >= 4:
        Qi = W[:, :4]
        Si, S2i, Gi = sums(Qi)
        iso = bool(np.linalg.norm(Si) < 1e-12 and np.abs(Gi - np.eye(4)).max() < 1e-12)
    rows2.append(dict(N=N, complex_dof=N - 3, true_square_dof=4*(N - 3), iso_dof=(4*(N - 5) if N >= 5 else None),
                      space_dirs=k, measured_space_rank=int(rank_space), true_square_example_ok=bool(ok_true), iso_example_ok=iso))
res["table2"] = rows2
md += ["## 表 2　自由度と、空間として開ける方向の数", "",
       "自由度は、閉塞の条件と読めない全体の対称性（複素数では λ∈ℂ、四元数の本当の二乗の側では 3 次元の回転と実数倍、等方の側では 4 次元の回転と実数倍）を除いた形の数。空間の方向の数は min(3, N−2)：空間の 3 本は q₀ と直交し、和が 0 の N−1 次元のうち残りの N−2 次元にしか入れない。各 N で実際に解を一つ作り、閉塞が成り立つことと空間成分の階数を確かめた。", "",
       "| N | 複素数（複素） | 四元数・本当の二乗（実） | 四元数・等方（実） | 空間の方向 min(3, N−2) | 作った解の空間成分の階数 | 本当の二乗の解 | 等方の解 |",
       "|---|---|---|---|---|---|---|---|"]
for r in rows2:
    md.append(f"| {r['N']} | {r['complex_dof']} | {r['true_square_dof']} | {r['iso_dof'] if r['iso_dof'] is not None else '解なし'} | {r['space_dirs']} | {r['measured_space_rank']} | {'成立' if r['true_square_example_ok'] else '不成立'} | {'成立' if r['iso_example_ok'] else ('—' if r['iso_example_ok'] is None else '不成立')} |")
md.append("")

# ---------------------------------------------------------------- 表 3
phi = (1 + 5 ** 0.5) / 2
polys = {
    "一つの複素平面（n＝±i）": [(1, 0, 0), (-1, 0, 0)],
    "正四面体": [(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)],
    "正八面体": [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)],
    "立方体": list(itertools.product([1, -1], repeat=3)),
    "正二十面体": [(0, 1, phi), (0, -1, phi), (0, 1, -phi), (0, -1, -phi), (1, phi, 0), (-1, phi, 0), (1, -phi, 0), (-1, -phi, 0),
                (phi, 0, 1), (-phi, 0, 1), (phi, 0, -1), (-phi, 0, -1)],
}
ip = 1 / phi
dodeca = list(itertools.product([1, -1], repeat=3))
for a, b in itertools.product([1, -1], repeat=2):
    dodeca += [(0, a*ip, b*phi), (a*ip, b*phi, 0), (a*phi, 0, b*ip)]
polys["正十二面体"] = dodeca
rows3 = []
for name, P in polys.items():
    n = np.array(P, float); n /= np.linalg.norm(n, axis=1)[:, None]
    Q = np.array([s*np.r_[1, m]/np.sqrt(2) for m in n for s in (1, -1)])
    S, S2, G = sums(Q)
    diag = np.diag(G)
    off = np.abs(G - np.diag(diag)).max()
    q8 = max(np.linalg.norm(qpow(q, 8) - np.r_[1, 0, 0, 0]) for q in Q)
    q4 = max(np.linalg.norm(qpow(q, 4) + np.r_[1, 0, 0, 0]) for q in Q)
    null = max(abs(q[0]**2 - q[1:] @ q[1:]) for q in Q)
    sq_is_n = max(np.linalg.norm(qmul(q, q) - np.r_[0, q[1:]*np.sqrt(2)*np.sign(q[0])]) for q in Q)
    rows3.append(dict(name=name, N=len(Q), abs_sum_q=float(np.linalg.norm(S)), abs_sum_q2=float(np.linalg.norm(S2)),
                      gram_diag=[float(x) for x in diag], gram_offdiag_max=float(off), ratio=float(diag[0]/diag[1]) if diag[1] > 1e-12 else None,
                      max_q8_minus_1=float(q8), max_q4_plus_1=float(q4), max_null=float(null), max_q2_minus_n=float(sq_is_n)))
E = np.eye(4); C16 = np.vstack([E, -E])
S, S2, G = sums(C16)
rows3.append(dict(name="対照：正十六胞体 {±1, ±i, ±j, ±k}", N=8, abs_sum_q=float(np.linalg.norm(S)), abs_sum_q2=float(np.linalg.norm(S2)),
                  sum_q2=[float(x) for x in S2], gram_diag=[float(x) for x in np.diag(G)], gram_offdiag_max=float(np.abs(G - np.diag(np.diag(G))).max())))
res["table3"] = rows3
md += ["## 表 3　最も対称な解 q ＝ ±e^(π n/4) ＝ ±(1 ＋ n)/√2 の一覧", "",
       "n は各多面体の頂点の方向（大きさ 1）。|Σq|、|Σq²| は和の大きさ、二次モーメントは Σ qqᵀ の対角（q₀, q₁, q₂, q₃）。光円錐は max |q₀² − |v|²|、q² ＝ n は max |q² − n|。", "",
       "| 方向 n の配置 | N | |Σq| | |Σq²| | 二次モーメントの対角 | 非対角の最大 | max |q₀²−|v|²| | max |q²−n| | max |q⁴＋1| | max |q⁸−1| |",
       "|---|---|---|---|---|---|---|---|---|---|"]
for r in rows3[:-1]:
    md.append(f"| {r['name']} | {r['N']} | {r['abs_sum_q']:.1e} | {r['abs_sum_q2']:.1e} | {', '.join(f'{x:.3f}' for x in r['gram_diag'])} | {r['gram_offdiag_max']:.1e} | {r['max_null']:.1e} | {r['max_q2_minus_n']:.1e} | {r['max_q4_plus_1']:.1e} | {r['max_q8_minus_1']:.1e} |")
r = rows3[-1]
md += ["", f"対照：正十六胞体 {{±1, ±i, ±j, ±k}}（N＝8）は |Σq|＝{r['abs_sum_q']:.1e}、二次モーメントの対角 {', '.join(f'{x:.0f}' for x in r['gram_diag'])}（等方的な 13 本を満たす）だが、Σq² ＝ ({', '.join(f'{x:+.0f}' for x in r['sum_q2'])}) で 0 にならない。",
       "", "一つの複素平面の場合だけ、二次モーメントの対角が (2, 2, 0, 0) で空間の 1 方向しか開かない。正四面体から先は空間 3 方向が等方的で、比は q₀ : q₁ : q₂ : q₃ ＝ 3 : 1 : 1 : 1。", ""]

# ---------------------------------------------------------------- 表 4
# 段 1：複素数 2N 個 → 四元数 N 個（逆写像で元に戻ることを確かめる）
N = 6
Z = rng.normal(size=2*N) + 1j*rng.normal(size=2*N)
Qz = np.array([[Z[2*k].real, Z[2*k].imag, Z[2*k + 1].real, Z[2*k + 1].imag] for k in range(N)])  # q ＝ z₁ ＋ z₂ j
Zb = np.array([c for q in Qz for c in (q[0] + 1j*q[1], q[2] + 1j*q[3])])
stage1 = float(np.abs(Zb - Z).max())
# z₁ ＋ z₂ j の本当の二乗が (z₁² − |z₂|²) ＋ 2 Re(z₁)·z₂ j になることを確かめる
z1, z2 = Z[0], Z[1]
q = np.array([z1.real, z1.imag, z2.real, z2.imag])
lhs = qmul(q, q)
a = z1**2 - abs(z2)**2; b = 2*z1.real*z2
rhs = np.array([a.real, a.imag, b.real, b.imag])
cd_check = float(np.abs(lhs - rhs).max())
# 段 2：中心投影 u ＝ v/q₀。|q| を変えても、q を −q にしても像は同じ
q = rng.normal(size=4)
u = q[1:]/q[0]
scale_inv = float(np.abs((3.7*q)[1:]/(3.7*q)[0] - u).max())
sign_inv = float(np.abs((-q)[1:]/(-q)[0] - u).max())
# 3 次元の回転（両側から同じ単位を掛ける）で閉塞が保たれ、左右を別々に回すと保たれない
n = np.array(polys["正二十面体"], float); n /= np.linalg.norm(n, axis=1)[:, None]
Qi = np.array([s*np.r_[1, m]/np.sqrt(2) for m in n for s in (1, -1)])
w = rng.normal(size=4); w /= np.linalg.norm(w); wc = w*np.r_[1, -1, -1, -1]
w2 = rng.normal(size=4); w2 /= np.linalg.norm(w2); w2c = w2*np.r_[1, -1, -1, -1]
conj_rot = np.linalg.norm(sums(np.array([qmul(qmul(w, x), wc) for x in Qi]))[1])
lr_rot = np.linalg.norm(sums(np.array([qmul(qmul(w, x), w2c) for x in Qi]))[1])
res["table4"] = dict(stage1_roundtrip_max_err=stage1, cayley_dickson_square_max_err=cd_check, gnomonic_scale_invariance=scale_inv,
                     gnomonic_sign_invariance=sign_inv, sum_q2_after_conjugation=float(conj_rot), sum_q2_after_left_right=float(lr_rot))
md += ["## 表 4　写像の各段：全投射か部分投射か", "",
       "| 段 | 写像 | 読める成分 | 逃げる成分 | 判定 | 数値で確かめたこと |", "|---|---|---|---|---|---|",
       f"| 1 | {{zₙ}}（2N 個）→ {{qₙ}}（N 個）、q ＝ z₁ ＋ z₂ j | 実数 4N 個すべて | なし | 全投射 | 逆写像で元に戻る（最大誤差 {stage1:.1e}）。本当の二乗は (z₁² − |z₂|²) ＋ 2 Re(z₁)·z₂ j（誤差 {cd_check:.1e}） |",
       f"| 2 | q → 北極の接超平面への中心投影 u ＝ v/q₀ ＝ n tan θ | 比 v/q₀ の 3 成分 | 大きさ |q|、符号 ±、赤道（q₀＝0 は無限遠） | 部分投射 | 3.7 倍しても像は同じ（差 {scale_inv:.1e}）、−q でも同じ（差 {sign_inv:.1e}） |",
       f"| 3 | 全体の対称性 | — | 両側から同じ単位を掛ける 3 次元の回転と実数倍 | 読めない | 正二十面体の 24 個で、回転後の |Σq²| ＝ {conj_rot:.1e}（保たれる）。左右を別々に回すと {lr_rot:.2f}（保たれない） |", "",
       "写像の連鎖全体として部分投射である。", ""]

# ---------------------------------------------------------------- 表 5
md += ["## 表 5　記法", "", "| 記号 | 意味 |", "|---|---|",
       "| ℍ | 四元数全体 |", "| q ＝ q₀ ＋ q₁i ＋ q₂j ＋ q₃k | 四元数（qᵢ は実数）。i² ＝ j² ＝ k² ＝ ijk ＝ −1 |",
       "| s（＝ q₀）, v | スカラー部（実部）とベクトル部（純虚部、3 成分） |", "| q̄ ＝ q₀ − v | 共役。ノルム |q|² ＝ q q̄ ＝ q₀² ＋ |v|² |",
       "| q² ＝ (q₀² − |v|²) ＋ 2 q₀ v | 本当の二乗（共役をとらない） |", "| n | 大きさ 1 の純虚四元数（n² ＝ −1） |",
       "| q ＝ |q| e^(n θ) ＝ |q|(cos θ ＋ n sin θ) | 極形式。θ は北極（q ＝ 1）からの角 |", "| S³ | 大きさ 1 の四元数全体（四元円） |",
       "| u ＝ v/q₀ | 北極の接超平面 q₀ ＝ 1 の上の座標（中心投影の像） |", "| q ＝ z₁ ＋ z₂ j | 複素数 2 つによる表し方（ℍ ≅ ℂ²） |", ""]

def escape_cell_bars(line):
    """表の行で、区切りでない縦棒（絶対値の記号）を \\| にエスケープする。区切りは前後が空白か行頭・行末のもの。"""
    if not line.startswith("|") or set(line) <= set("|-"):
        return line
    out = []
    for i, ch in enumerate(line):
        if ch == "|" and not ((i == 0 or line[i-1] == " ") and (i == len(line) - 1 or line[i+1] == " ")):
            out.append("\\|")
        else:
            out.append(ch)
    return "".join(out)


md = [escape_cell_bars(l) for l in md]
(OUT / "tables.md").write_text("\n".join(md), encoding="utf-8")
(OUT / "tables.json").write_text(json.dumps(res, ensure_ascii=False, indent=1), encoding="utf-8")
# 合否の要約
checks = {
    "表2 本当の二乗の解（全 N）": all(r["true_square_example_ok"] for r in rows2),
    "表2 空間成分の階数 ＝ min(3, N−2)": all(r["measured_space_rank"] == r["space_dirs"] for r in rows2),
    "表2 等方の解（N≥5）": all(r["iso_example_ok"] for r in rows2 if r["N"] >= 5),
    "表3 |Σq|, |Σq²| < 1e-12": all(r["abs_sum_q"] < 1e-12 and r["abs_sum_q2"] < 1e-12 for r in rows3[:-1]),
    "表3 光円錐・q²＝n・q⁸＝1": all(r["max_null"] < 1e-12 and r["max_q2_minus_n"] < 1e-12 and r["max_q8_minus_1"] < 1e-12 for r in rows3[:-1]),
    "表3 正十六胞体 Σq² ＝ −4": abs(rows3[-1]["sum_q2"][0] + 4) < 1e-12,
    "表4 段 1 と段 2": stage1 < 1e-12 and cd_check < 1e-12 and scale_inv < 1e-12 and sign_inv < 1e-12,
    "表4 回転：同じ単位で保たれ、左右別で崩れる": conj_rot < 1e-12 and lr_rot > 1e-3,
}
for k, v in checks.items():
    print(("合格 " if v else "不合格 ") + k)
(OUT / "checks.txt").write_text("\n".join(("合格 " if v else "不合格 ") + k for k, v in checks.items()) + "\n", encoding="utf-8")
