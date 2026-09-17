# Claude Code 指示書 — 二波99走行実験の完全やり直し v2

日付: 2026-09-17  
状態: **旧99走行系列を破棄して、新しい理論相互作用でゼロから再実験するための確定指示**

## 0. 最重要方針

旧99走行の結果・Sparse/Lanczos 実装・旧互換ゲート・旧解析結果は **再現対象にしない**。
旧系列は比較基準でも正本でもない。結果救済のための追試は行わない。

今回の実験は、理論上必要な相互作用だけを実装した新しい正本を作り、その正本を固定して99走行を最初から生成する。

禁止事項:

- 旧99走行 states を初期値・比較基準・互換ゲートに使わない
- Sparse / matrix-free / Lanczos / Krylov を使わない
- 高速化のためにアルゴリズムを置換しない
- `np.angle(z)` による振幅破棄を行わない
- 波ごとの正規化を行わない
- step ごとの再正規化を行わない
- clip / round / threshold により状態を書き換えない
- 結果を見て実験条件を変えない
- 旧解析結果を流用しない

---

## 1. 新しい理論正本相互作用

関係波を

$$
z_e=a_e+i b_e
$$

とする。

二つの関係波 $e,f$ の最小の実反対称・双線形・大域 U(1) 不変な相互作用を

$$
K_{ef}
=
\lambda A_{ef}\,\operatorname{Im}(\bar z_e z_f)
=
\lambda A_{ef}(a_e b_f-b_e a_f)
$$

と定義する。

今回の基準単位では **$\lambda=1$** とする。$\lambda$ や $M/H$ 等による追加規格化は新しい物理仮説なので導入禁止。

$A$ は従来どおり、二つの関係辺が頂点を共有するときだけ 1 となる対称隣接行列とする。

この定義から

$$
K^T=-K
$$

であり、振幅は捨てない。$z_e=0$ または $z_f=0$ ならその対の相互作用も 0 になる。

大域位相変換

$$
z_e\to e^{i\phi}z_e
$$

に対して $\bar z_e z_f$ は不変なので、相互作用も不変である。

### 重要

旧3段処理

1. $H=A\odot(\bar z z^T)$
2. $z\to e^{i\arg z}$ で振幅破棄
3. 虚部だけ取り $iK$ を作る

という手順は新正本では採用しない。

新正本では、物理的相互作用を最初から

$$
K_{ef}=A_{ef}\operatorname{Im}(\bar z_e z_f)
$$

と直接定義する。

「段3」は後処理ではなく、この実反対称双線形形式そのものに吸収する。

---

## 2. 数値実装 — Dense 全固有分解のみ

相互作用カーネルは外部フォルダから import してはならない。

この実験フォルダ内に

`interaction_kernel_theory_v1.py`

を **物理的なローカルコピーとして作成**し、99走行開始前に SHA256 を固定する。
ランナーはこの実験フォルダ内のコピーだけを使用する。

次段階の相互作用改訂実験では、このファイルをさらに別名コピーして編集する。正本 v1 自体は変更しない。

### 正本コード

```python
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import math
import numpy as np


def edges(N):
    a, b = np.triu_indices(N, k=1)
    return a.astype(np.int64), b.astype(np.int64)


def adjacency(N):
    ea, eb = edges(N)
    M = len(ea)
    A = np.zeros((M, M), dtype=np.float64)
    for e in range(M):
        share = (ea == ea[e]) | (ea == eb[e]) | (eb == ea[e]) | (eb == eb[e])
        share[e] = False
        A[e, share] = 1.0
    return A


def K_of(z, A):
    # C_ef = A_ef * conj(z_e) * z_f
    C = A * (np.conj(z)[:, None] * z[None, :])

    # 理論式 K_ef = A_ef Im(conj(z_e) z_f)
    # 有限精度でも K^T = -K を構造的に満たすよう、transpose 反対称部として構成する。
    K = ((C - C.T) / (2j)).real
    np.fill_diagonal(K, 0.0)
    return K.astype(np.float64, copy=False)


def one_step(z, A, den):
    K = K_of(z, A)

    # H=iK は Hermitian。これは物理的な追加段階ではなく、
    # exp(dtau*K) を full dense Hermitian eigendecomposition で計算するための数値表現。
    H = (1j * K).astype(np.complex128, copy=False)
    w, V = np.linalg.eigh(H)
    dtau = np.float64(2.0 * math.pi / den)
    phase = np.exp(-1j * dtau * w)
    return (V @ (phase * (V.conj().T @ z))).astype(np.complex128, copy=False)
```

### 絶対禁止

以下を追加してはならない。

```python
np.exp(1j*np.angle(z))
z / np.abs(z)
K / norm
z / norm
Lanczos
Krylov
expm_multiply
round(...)
clip(...)
```

また `K_of` を「高速だから」という理由で別形式へ書き換えてはならない。

---

## 3. 実行前の構造監査

99走行開始前に、ランダム状態と実初期値の双方で以下を確認し、JSON に保存する。

1. `A.T == A` の最大残差
2. `K.T == -K` の最大残差
3. `H.conj().T == H` の最大残差
4. 大域 U(1) 同変性
5. $z_e=0$ を含む人工状態で、その成分を介する $K$ が 0 になること
6. 一 step 後の
   - $H_0=Z^\dagger Z$
   - $Q_2=Z^T Z$
   の差
7. NaN/Inf 不在

これらは品質監査であり、結果を見て閾値や力学を変更してはならない。

---

## 4. 実験環境を固定

99走行は一つの同一環境で実行する。
開始前に以下を `environment.json` と `environment.txt` に保存する。

- OS / kernel
- CPU
- Python version
- NumPy version
- SciPy version（解析で使う場合）
- BLAS / LAPACK 情報
- `np.show_config()`
- git commit SHA
- `interaction_kernel_theory_v1.py` SHA256
- 初期値生成コード SHA256
- run code SHA256

99走行途中で Python / NumPy / BLAS / コードを変更禁止。
変更が必要になった場合は、その系列を中止して新しい実験系列として最初から開始する。

---

## 5. 初期値

二波初期条件は旧99走行で使用した数学的定義をそのまま使う。
ただし旧実験コードを外部 import しない。
初期値生成器も実験フォルダへ物理コピーし、SHA256 を固定する。

$$
U=e^{2\pi i/L},\qquad
c_a=e^{-i\pi/(2L)},\qquad
c_b=e^{+i\pi/(2L)}
$$

$$
z_{jk}(0)=c_aU^{m_a(k-j)}+c_bU^{m_b(k-j)},\qquad j<k.
$$

dtype は `complex128`、実数は `float64` とする。

---

## 6. 99走行条件

### 主系列 72 runs

$$
L\in\{8,10,12,16,20,24,32,40\}
$$

各 L について 9 組:

- 奇奇: `(1,3)`, `(1,5)`, `(3,5)`
- 偶偶: `(2,4)`, `(2,6)`, `(4,6)`
- 奇偶: `(1,2)`, `(2,3)`, `(3,6)`

主系列は `den=L`。

### gcd 系列 18 runs

$$
L\in\{12,16,20,24,32,40\}
$$

各 L について

- `(3,9)`
- `(4,8)`
- `(5,10)`

`den=L`。

### 刻み対照 9 runs

$$
L\in\{12,24,32\}
$$

各 L について

- `(1,3)`
- `(2,4)`
- `(2,3)`

`den=40`。

合計 99 runs。

---

## 7. 走行長と保存

全 run:

```text
t = 0 ... 4096
```

4097 状態を間引かず全保存する。

各 run フォルダに最低限:

- `states.npz` — 全 $Z(t)$、complex128
- `metadata.json`
- `invariants.csv`

`invariants.csv` には毎 step:

- t
- $H=Z^\dagger Z$
- Re$(Q_2)$, Im$(Q_2)$, $Q_2=Z^T Z$
- $\|K(Z)Z\|/\|Z\|$
- antisymmetry residual $\max|K+K^T|$

を記録する。

生 states が正本であり、解析結果は派生物とする。

---

## 8. 旧99走行との比較は禁止

今回の目的は旧99走行を救済することではない。

したがって実行中・実行直後に

- old vs new trajectory diff
- old fixed-point count
- old 94/5 分類
- old inflation/onset 比較
- old spectral result 比較

を行わない。

旧データはこの新系列の判断基準にしない。

---

## 9. QA

各 run で以下を機械的に確認する。

- 4097 states 存在
- NaN/Inf なし
- metadata 完備
- SHA256 記録
- H, Q2 の drift を数値として記録
- antisymmetry residual を記録

FAIL でも states を削除しない。ただし原因を記録し、その run を科学解析へ混ぜない。

結果を見て QA 閾値を変更しない。

---

## 10. 全解析は生データ完成後にゼロから行う

99走行が完了し SHA256 が確定するまで、科学的解析を開始しない。

完了後、新しい `analysis_v1/` を作り、旧解析出力を一切読み込まず、99本の新 `states.npz` だけから全解析を再生成する。

最低限、以下を新データから再計算する。

- 固定点 / 非固定点判定
- $H$, $Q_2$, $\chi=|Q_2|/H$
- 実部・虚部 Gram 構造
- 振幅分布・位相分布
- centroid
- 初期2平面からの $H_\perp/H$
- 生成子 $K$ の rank / spectrum / eigenvector structure
- 時系列・周期性・FFT 等（解析開始後に仕様を固定して実施）

旧99走行から得た「5固定点」「94離脱」等を前提に分類しない。

---

## 11. 完了報告

99走行完了時点では、物理解釈をせず以下だけ報告する。

1. 99/99 run 完了数
2. QA PASS / FAIL 数
3. 新正本 kernel SHA256
4. generator SHA256
5. run code SHA256
6. environment fingerprint
7. 各 `states.npz` SHA256
8. H/Q2 drift の最大値
9. antisymmetry residual 最大値
10. 生データ保存先

その後、木原・ChatGPT が生データと監査記録を確認してから解析段階へ進む。

---

## 12. この実験の原則

この99走行では、旧3段処理の再現を目的としない。

物理写像は最初から

$$
\boxed{K_{ef}=A_{ef}\operatorname{Im}(\bar z_e z_f)}
$$

$$
\boxed{Z_{n+1}=\exp(\Delta\tau K(Z_n))Z_n}
$$

だけである。

数値上 `H=iK` を作るのは full dense Hermitian eigendecomposition を使うための表現であり、新たな物理段階ではない。

**振幅を捨てない。実反対称性を保つ。余計な規格化を入れない。近似積分器へ置換しない。**
