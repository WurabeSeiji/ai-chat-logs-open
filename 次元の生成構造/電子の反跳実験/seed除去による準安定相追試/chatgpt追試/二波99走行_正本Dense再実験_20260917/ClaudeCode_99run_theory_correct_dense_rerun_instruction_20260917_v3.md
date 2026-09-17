# Claude Code 指示書 — 二波99走行実験の完全やり直し v3

日付: 2026-09-17  
状態: **旧99走行系列を破棄し、理論式を直接・非近似で実装してゼロから再実験するための確定指示**

## 0. 最重要原則

旧99走行の結果・Sparse/Lanczos 実装・旧互換ゲート・旧解析結果は再現対象にしない。
旧系列は比較基準でも正本でもない。結果救済のための追試は行わない。

今回の実験では、理論式を直接実装した新しい相互作用正本を作成し、その正本を固定して99走行を最初から生成する。

**最重要:** 相互作用・時間発展に対して、速度改善、低ランク化、疎近似、反復近似、収束打切り、閾値処理、正規化、丸め、補正、射影、近似指数写像を一切導入しない。

この指示書でいう「厳密」とは、指定した `float64/complex128` 環境内で、理論式どおりの dense 全行列を毎 step 構成し、全固有値・全固有ベクトルを `np.linalg.eigh` で求め、スペクトル指数写像を全成分に適用することを意味する。機械丸め誤差そのものは不可避だが、アルゴリズム上の近似を追加してはならない。

禁止事項:

- Sparse / matrix-free / Lanczos / Krylov / eigsh を使わない
- truncated SVD / randomized SVD / low-rank approximation を使わない
- `expm_multiply`、部分固有分解、反復収束型の近似指数写像を使わない
- tolerance による早期終了を使わない
- 小さい値を 0 とみなす shortcut を入れない
- `np.angle(z)` による振幅破棄を行わない
- 波ごとの正規化を行わない
- step ごとの再正規化を行わない
- clip / round / threshold により状態を書き換えない
- 対称性・反対称性をコード側で補正しない
- 結果を見て実験条件を変えない
- 旧解析結果を流用しない

---

## 1. 新しい理論正本相互作用

関係波を

$$
z_e=a_e+i b_e
$$

とする。

相互作用は理論式をそのまま

$$
\boxed{K_{ef}=A_{ef}\,\operatorname{Im}(\bar z_e z_f)}
$$

すなわち

$$
\boxed{K_{ef}=A_{ef}(a_e b_f-b_e a_f)}
$$

と定義する。

追加係数、規格化、重み、補正は導入しない。基準単位では係数は 1 とする。

`A` は従来どおり、二つの関係辺が頂点を共有するときだけ 1 となる対称隣接行列とする。

理論上、対称な `A` に対して

$$
K^T=-K
$$

が自動的に成立する。したがってコード側で

```python
(C - C.T)/(2j)
(K - K.T)/2
```

などの反対称化を行ってはならない。

反対称性は「作る」のではなく「結果として観測・監査する」。

---

## 2. 相互作用正本コード

相互作用カーネルは外部フォルダから import しない。
この実験フォルダ内に

`interaction_kernel_theory_v1.py`

を物理コピーとして作成し、99走行開始前に SHA256 を固定する。
ランナーはこのファイルだけを使用する。

正本コードは以下を基準とし、これ以外の最適化・等価変形・整理を行わない。

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
    return (A * np.imag(np.conj(z)[:, None] * z[None, :])).astype(np.float64, copy=False)


def one_step(z, A, den):
    K = K_of(z, A)

    # H=iK は exp(dtau*K) を full dense Hermitian eigendecomposition で
    # 直接評価するための数値表現であり、物理的な追加相互作用ではない。
    H = (1j * K).astype(np.complex128, copy=False)

    # 全固有値・全固有ベクトルを毎step計算する。部分固有分解は禁止。
    w, V = np.linalg.eigh(H)

    dtau = np.float64(2.0 * math.pi / den)
    phase = np.exp(-1j * dtau * w)

    return (V @ (phase * (V.conj().T @ z))).astype(np.complex128, copy=False)
```

### この関数への禁止変更

以下はすべて禁止する。

```text
np.exp(1j*np.angle(z))
z/abs(z)
K/norm
z/norm
(C-C.T)/(2j)
(K-K.T)/2
Lanczos
Krylov
eigsh
svds
randomized_svd
expm_multiply
partial eigenspectrum
Taylor truncation
series truncation
round
clip
threshold-based zeroing
if norm(K) < eps: return z
```

「数学的に同じ」「高速だから」「数値的に安定だから」という理由でも変更禁止。

---

## 3. 非近似計算の絶対条件

99走行中、相互作用と時間発展は各 step で必ず次の順序を完全実行する。

1. 現在の full state `z` から full dense `K` を直接構成する。
2. `K` の全要素を保持する。疎化しない。
3. `H = 1j*K` を full dense 行列として構成する。
4. `np.linalg.eigh(H)` により全 `M` 個の固有値と全固有ベクトルを求める。
5. 全固有値について `exp(-1j*dtau*w)` を評価する。
6. 全固有ベクトルを用いて `V @ (phase * (V.conj().T @ z))` を評価する。
7. 次状態を `complex128` のまま保存する。

以下を一切行わない。

- 収束判定
- tol 指定
- 反復回数制限
- モード打切り
- 固有値上位 k 個のみ使用
- 小固有値の削除
- 位相の丸め
- 状態の間引き
- 補間
- 平滑化
- 精度低下のための cast

99走行の途中で、計算時間が長いことを理由にアルゴリズムを変更してはならない。
時間が許容できない場合は実験を停止して報告し、勝手に近似版へ切り替えない。

---

## 4. 実行前の構造監査

99走行開始前に、ランダム状態と実初期値の双方で以下を計算し、`pre_run_structural_audit.json` に生の数値を保存する。

1. `max(abs(A - A.T))`
2. `max(abs(K + K.T))`
3. `max(abs(H - H.conj().T))`
4. 大域 U(1) 変換前後の `K` の最大差
5. `z_e=0` を含む人工状態で該当相互作用成分が 0 になること
6. 一 step 前後の `H0 = Z^dagger Z`
7. 一 step 前後の `Q2 = Z^T Z`
8. NaN/Inf 不在

**重要:** 監査値が期待と異なっても、その場で `K` を補正・対称化・正規化してはならない。
異常があれば99走行を開始せず停止し、原因を報告する。

---

## 5. 実験環境固定

99走行は一つの同一環境で実行する。
開始前に以下を保存する。

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

可能な限り計算スレッド数を 1 に固定し、実行前に値を記録する。

例:

```bash
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
```

実環境で無効な変数があってもよいが、実際に使用している BLAS/LAPACK とスレッド設定を記録する。

途中で Python / NumPy / BLAS / LAPACK / CPU / コードを変更しない。
変更が必要になった場合、その系列を中止し、新しい系列として最初から開始する。

---

## 6. 初期値

二波初期条件は数学的定義を固定する。

$$
U=e^{2\pi i/L},\quad
c_a=e^{-i\pi/(2L)},\quad
c_b=e^{+i\pi/(2L)}
$$

$$
z_{jk}(0)=c_aU^{m_a(k-j)}+c_bU^{m_b(k-j)},\quad j<k
$$

初期値生成器も実験フォルダへ物理コピーし SHA256 を固定する。
外部 import しない。

dtype は `complex128`、実数は `float64`。

初期値生成後、各 run の `z0` の SHA256 も記録する。

---

## 7. 99走行条件

### 主系列 72 runs

`L in {8,10,12,16,20,24,32,40}`

各 L について 9 組:

- 奇奇: `(1,3)`, `(1,5)`, `(3,5)`
- 偶偶: `(2,4)`, `(2,6)`, `(4,6)`
- 奇偶: `(1,2)`, `(2,3)`, `(3,6)`

主系列は `den=L`。

### gcd 系列 18 runs

`L in {12,16,20,24,32,40}`

各 L について:

- `(3,9)`
- `(4,8)`
- `(5,10)`

`den=L`。

### 刻み対照 9 runs

`L in {12,24,32}`

各 L について:

- `(1,3)`
- `(2,4)`
- `(2,3)`

`den=40`。

合計 99 runs。

---

## 8. 走行長と生データ保存

全 run:

`t = 0 ... 4096`

4097 状態を間引かず全保存する。

各 run フォルダに最低限:

- `states.npz` — 全 `Z(t)`、complex128
- `metadata.json`
- `invariants.csv`

`invariants.csv` には毎 step:

- t
- `H = Z^dagger Z`
- Re(Q2), Im(Q2), `Q2 = Z^T Z`
- `||K(Z)Z||/||Z||`
- `max(abs(K+K.T))`

を記録する。

生 `states.npz` が唯一の一次データであり、解析結果は派生物とする。

---

## 9. 解析も近似で置換しない

99走行の生データ完成と SHA256 固定が終わるまで、科学解析を開始しない。

解析開始後も、元データを縮約してから解析しない。
原則として 4097 全時刻・全成分を使用する。

禁止:

- downsampling
- smoothing
- binning
- interpolation
- randomized SVD
- truncated SVD
- partial eigenspectrum
- low-rank surrogate
- 途中データの丸め
- 図化用データを一次解析値として再利用

固有値解析は full dense spectrum を保存する。
SVDを使う場合も full SVD を使い、全 singular values を保存する。
FFTを使う場合も元時系列全点を入力し、前処理を行う場合は別解析として明記する。

rank、固定点、周期、onset 等で数値閾値が必要な場合:

1. まず生の連続値を保存する。
2. 閾値は解析仕様書で事前固定する。
3. 閾値を結果を見て変更しない。
4. 閾値依存分類と生データを混同しない。

---

## 10. QA

各 run で以下を機械的に確認する。

- 4097 states 存在
- NaN/Inf なし
- metadata 完備
- SHA256 記録
- H, Q2 drift の生値を記録
- antisymmetry residual の生値を記録

FAIL でも states を削除しない。ただし原因を記録し、その run を科学解析へ混ぜない。
結果を見て QA 閾値を変更しない。

---

## 11. 完了報告

99走行完了時点では、物理解釈をせず以下だけ報告する。

1. 99/99 run 完了数
2. QA PASS / FAIL 数
3. 新正本 kernel SHA256
4. generator SHA256
5. run code SHA256
6. environment fingerprint
7. 各 `z0` SHA256
8. 各 `states.npz` SHA256
9. H/Q2 drift 最大値
10. antisymmetry residual 最大値
11. 生データ保存先

その後、木原・ChatGPT が生データと監査記録を確認してから解析段階へ進む。

---

## 12. この実験の唯一の物理写像

この99走行で使う物理式は以下だけである。

$$
\boxed{K_{ef}=A_{ef}\operatorname{Im}(\bar z_e z_f)}
$$

$$
\boxed{Z_{n+1}=\exp(\Delta\tau K(Z_n))Z_n}
$$

この二式の間に、正規化・射影・補正・近似・低ランク化・疎化・収束打切りを挿入してはならない。

**速度より同一性を優先する。計算が重ければ停止して報告し、近似版へ変更しない。**
