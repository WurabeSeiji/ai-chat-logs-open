#!/usr/bin/env python3
"""
内と外の選択問題：補助計算（2026-10-19）

`octonion_closure_check_20261019.py` と同じフォルダに置いて実行する
（代数の定義をそのファイルから読み込む）。検討メモ
`内と外の選択問題_八元数分割と読み出し_20261019.md` の §1.2、§1.6、§5 の数値はこれで出した。

A. 安定化群の次元：G2 の中で、観測者が時空 H = span(e0..e3) を固定したとき、
   ℓ = e4 を固定したとき、両方を固定したときに残る対称性の次元。
   残った対称性が内部半分 (R, Q) にどう作用するか（R と Q を混ぜるか）。
B. 内と外の入れ替えは自己同型か。時空側として選べる四元数部分代数の空間の次元。
   同じ閉塞解を違う分け方で読んだとき、和則 Σ(q0²−|v|²) = Σ|p|² が保たれるか、
   要素ごとの m² は変わるか。
C. 二つの半分の構造的な非対称：どちらが積で閉じるか、1 はどちらにあるか、
   自己関係はどこを指すか。
"""
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
src = open(os.path.join(HERE, "octonion_closure_check_20261019.py")).read()
exec(src.split("results = {")[0])          # qmul, omul, oconj, onorm2, ojordan, oexp, E

SEED = 20261019
rng = np.random.default_rng(SEED)

def nulldim(A, tol=1e-9):
    s = np.linalg.svd(A, compute_uv=False)
    return int(np.sum(s < tol * max(1.0, s[0]))) + (A.shape[1] - len(s))

def nullbasis(A, tol=1e-9):
    u, s, vt = np.linalg.svd(A)
    return vt[int(np.sum(s > tol * s[0])):]

# ---------------------------------------------------------------------------
# A. 安定化群の次元
# ---------------------------------------------------------------------------
rows = []
for i in range(8):
    for j in range(8):
        prod = omul(E[i], E[j])
        for c in range(8):
            row = np.zeros(64)
            for m in range(8): row[c*8 + m] += prod[m]
            for m in range(8): row[m*8 + i] -= omul(E[m], E[j])[c]
            for m in range(8): row[m*8 + j] -= omul(E[i], E[m])[c]
            rows.append(row)
M = np.array(rows)                                   # 微分（自己同型の無限小）の方程式

cH = [np.eye(64)[c*8 + i] for i in range(4) for c in range(4, 8)]   # D が H を H に写す
cl = [np.eye(64)[c*8 + 4] for c in range(8)]                        # D(e4) = 0
cH0 = cH + [np.eye(64)[c*8 + i] for c in range(4) for i in range(4)]  # H 上で自明

print("A. 安定化群の次元")
print("  dim G2 =", nulldim(M), "(期待 14)")
print("  stab(H) =", nulldim(np.vstack([M, np.array(cH)])), "(期待 so(4) = 6)")
print("  stab(l) =", nulldim(np.vstack([M, np.array(cl)])), "(期待 su(3) = 8)")
print("  stab(H) ∩ stab(l) =", nulldim(np.vstack([M, np.array(cH + cl)])), "(期待 so(3) = 3)")
print("  H 上で自明な部分 =", nulldim(np.vstack([M, np.array(cH0)])), "(期待 su(2)_内部 = 3)")

ns = nullbasis(np.vstack([M, np.array(cH)]))
mixR = max(max(np.abs(v.reshape(8, 8)[4:8, 4:8][0, 1:]).max(),
               np.abs(v.reshape(8, 8)[4:8, 4:8][1:, 0]).max()) for v in ns)
print("  stab(H) の中で R(e4) と Q(e5..7) を混ぜる成分の最大 =", f"{mixR:.3f}",
      "→ 0 でない：時空だけ固定した観測者には R と Q は分離できない")
ns = nullbasis(np.vstack([M, np.array(cl)]))
mx = max(np.abs(v.reshape(8, 8)[4:8, 0:4]).max() for v in ns)
print("  su(3) が H を H の外へ写す成分の最大 =", f"{mx:.3f}", "→ su(3) は時空を保たない")

# ---------------------------------------------------------------------------
# B. 内と外の入れ替えと、分け方の自由度
# ---------------------------------------------------------------------------
def swap(o): return np.concatenate([o[4:], o[:4]])

e_auto = e_norm = e_jord = 0.0
for _ in range(2000):
    x, y = rng.normal(size=(2, 8))
    e_auto = max(e_auto, np.linalg.norm(swap(omul(x, y)) - omul(swap(x), swap(y))))
    e_jord = max(e_jord, np.linalg.norm(swap(ojordan(x, y)) - ojordan(swap(x), swap(y))))
    e_norm = max(e_norm, abs(onorm2(swap(x)) - onorm2(x)))
print("\nB. 内と外の入れ替え q ↔ p")
print(f"  積を保つか：最大誤差 {e_auto:.1f} → 自己同型ではない")
print(f"  ジョルダン積を保つか：最大誤差 {e_jord:.1f} → 保たない")
print(f"  ノルムを保つか：最大誤差 {e_norm:.1e} → 保つ")
print("  時空側に選べる四元数部分代数の空間：6 + 5 − 3 =", 6 + 5 - 3,
      "= dim G2/SO(4) =", 14 - 6)

def gen_closure(N):
    one = np.ones(N) / np.sqrt(N)
    x0 = rng.normal(size=N); x0 -= one * (one @ x0)
    P = np.eye(N) - np.outer(one, one) - np.outer(x0, x0) / (x0 @ x0)
    X = np.array([P @ rng.normal(size=N) for _ in range(7)])
    X *= np.sqrt((x0 @ x0) / np.sum(X * X))
    return np.column_stack([x0, X.T])

def rand_subalgebra():
    A = rng.normal(size=8); A[0] = 0; A /= np.linalg.norm(A)
    B = rng.normal(size=8); B[0] = 0; B -= (A @ B) * A; B /= np.linalg.norm(B)
    C = omul(A, B); C[0] = 0; C /= np.linalg.norm(C)
    return np.array([E[0], A, B, C])      # 4x8（行が基底）

S = gen_closure(10)
print(f"  閉塞の残差：|Σo| = {np.linalg.norm(S.sum(0)):.1e}, "
      f"|Σo²| = {np.linalg.norm(sum(ojordan(x, x) for x in S)):.1e}")
for t in range(4):
    Bas = rand_subalgebra()
    Q, _ = np.linalg.qr(np.hstack([Bas.T, rng.normal(size=(8, 8))])); Cm = Q[:, 4:]
    q = S @ Bas.T; p = S @ Cm
    m2 = q[:, 0]**2 - np.sum(q[:, 1:]**2, axis=1)
    n2 = np.sum(p**2, axis=1)
    print(f"  分け方 {t}: 要素ごとの (q0²−|v|²) − |p|² の最大ずれ = {np.max(np.abs(m2 - n2)):.2e};  "
          f"和則 Σ(q0²−|v|²) − Σ|p|² = {np.sum(m2) - np.sum(n2):.1e};  "
          f"Σq0 = {q[:, 0].sum():.1e}, |Σv| = {np.abs(q[:, 1:].sum(0)).max():.1e}, |Σp| = {np.abs(p.sum(0)).max():.1e}")

# ---------------------------------------------------------------------------
# C. 二つの半分の構造的な非対称
# ---------------------------------------------------------------------------
Bas = rand_subalgebra()
Q, _ = np.linalg.qr(np.hstack([Bas.T, rng.normal(size=(8, 8))])); Cm = Q[:, 4:].T
leakH = leakC = 0.0
for _ in range(500):
    x = rng.normal(size=4) @ Bas; y = rng.normal(size=4) @ Bas
    leakH = max(leakH, np.linalg.norm(omul(x, y) - Bas.T @ (Bas @ omul(x, y))))
    u = rng.normal(size=4) @ Cm; v = rng.normal(size=4) @ Cm
    leakC = max(leakC, np.linalg.norm(omul(u, v) - Cm.T @ (Cm @ omul(u, v))))
u = rng.normal(size=4) @ Cm; v = rng.normal(size=4) @ Cm; w = omul(u, v)
self_im = max(np.linalg.norm(omul(oconj(o), o)[1:]) for o in rng.normal(size=(500, 8)))
print("\nC. 二つの半分の非対称")
print(f"  H' は積で閉じる（漏れ {leakH:.1e}）")
print(f"  補空間は積で閉じない（漏れ {leakC:.1f}）")
print(f"  補空間どうしの積が H' に落ちる割合 = {np.linalg.norm(Bas.T @ (Bas @ w)) / np.linalg.norm(w):.6f}")
print(f"  1 の補空間成分 = {np.linalg.norm(Cm @ E[0]):.1e} → 1 は必ず H' 側")
print(f"  自己関係 ō o の虚部 = {self_im:.1e} → 自己関係は常に 1 を指す")
print("\n  → 内と外の選択は 8 パラメータ自由だが、二つの半分は入れ替え可能ではない。")
print("     非対称は「1 がどちらにあるか」の一点で、観測者は自分の 1 を含む側を時空として読む。")
