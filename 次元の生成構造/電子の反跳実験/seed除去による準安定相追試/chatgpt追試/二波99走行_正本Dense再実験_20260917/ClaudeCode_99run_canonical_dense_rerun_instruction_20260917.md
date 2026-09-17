# Claude Code 実行指示書
## 二波完全関係力学 99走行 — 正本 Dense 完全再実験

作成日: 2026-09-17
目的: 2026-09-16 の99走行系列を研究結果として破棄し、数値実験を正本実装・正本数値手順に戻して最初から再実行する。

---

# 0. 最重要原則

今回の再実験は「現在の結果を救済・比較するため」ではない。
旧99走行系列は採用しない。新しい正本系列をゼロから生成する。

以下を厳守すること。

1. 高速化しない。
2. 疎行列化しない。
3. matrix-free 化しない。
4. Lanczos / Krylov / tol 収束判定を使用しない。
5. 数学的に同値であっても別実装へ書き換えない。
6. 正本関数を再記述・整形・リファクタリングしない。
7. 実験途中で測定方法を変更しない。
8. 実験途中で解析方法を変更しない。
9. 結果を見て条件・閾値・走行長を変更しない。
10. 既存の99走行データ・解析結果を新系列へ混入させない。

特に、

> 理論的同値 != 数値実験としての同一性

を今回の最優先ルールとする。

---

# 1. 正本相互作用ソース — import ではなく byte-copy する

正本:

`run_N3_N40_stage123_v1_CANONICAL_REFERENCE.py`

正本 SHA256:

```text
1abf2353fee2e4f56f05e7a6f149fd086885136beb61ab571b48a56b09691567
```

## 必須手順

### 1.1 元ファイルの SHA256 を確認

```sh
shasum -a 256 run_N3_N40_stage123_v1_CANONICAL_REFERENCE.py
```

上記 SHA256 と 1 bit でも異なれば **即停止**。実験開始禁止。

### 1.2 新しい再実験フォルダへファイルそのものをコピー

外部ディレクトリの正本を直接 import してはならない。
関数だけを copy/paste してはならない。
AST 同一・ロジック同一で済ませてはならない。

OS のファイルコピーで、正本ファイル全体をそのまま新実験フォルダへコピーする。

例:

```sh
cp -p <正本絶対パス>/run_N3_N40_stage123_v1_CANONICAL_REFERENCE.py ./run_N3_N40_stage123_v1_CANONICAL_REFERENCE.py
```

### 1.3 コピー後の SHA256 を再確認

```sh
shasum -a 256 ./run_N3_N40_stage123_v1_CANONICAL_REFERENCE.py
```

必ず

```text
1abf2353fee2e4f56f05e7a6f149fd086885136beb61ab571b48a56b09691567
```

であること。

この一致を `PRE_RUN_AUDIT.md` と `PRE_RUN_AUDIT.json` の両方へ記録する。

## 重要

今回あえて **import 元を固定するのではなく、正本の byte-copy を新実験フォルダ内に置く**。
理由は、この正本再実験終了後の次段階で、このローカル正本コピーからさらに別コピーを作り、相互作用の不要段階を一因子ずつ停止する対照実験を行うためである。

ただし **本99走行の正本コピーは一文字も変更しない**。
次段階の変更実験は必ず別ファイル・別フォルダで行う。

---

# 2. 使用する相互作用・更新写像

99走行再実験では正本 `one_step` をそのまま使用する。

正本は次の処理を行う。

```python
H=H_of(np.exp(1j*np.angle(z)),A)
H=(1j*np.imag(H)).astype(np.complex128,copy=False)
w,V=np.linalg.eigh(H)
phase=np.exp(-1j*np.float64(2.0*math.pi/den)*w)
return (V@(phase*(V.conj().T@z))).astype(np.complex128,copy=False)
```

したがって本再実験では、

- `np.exp(1j*np.angle(z))` を変更しない
- `1j*np.imag(H)` を変更しない
- `np.linalg.eigh(H)` を変更しない
- `V @ (phase * (V.conj().T @ z))` を変更しない

こと。

**SparseKernel は使用禁止。**
**Lanczos は使用禁止。**
**Krylov は使用禁止。**
**tol は導入禁止。**

---

# 3. 99走行条件

99走行の条件網そのものは 2026-09-15 の事前登録済みサーベイと同じものを使用する。
条件を追加・削除・入れ替えない。

## 主系列 72 runs

`L in {8,10,12,16,20,24,32,40}` × 9組

- 奇奇: `(1,3), (1,5), (3,5)`
- 偶偶: `(2,4), (2,6), (4,6)`
- 奇偶: `(1,2), (2,3), (3,6)`

## gcd 系列 18 runs

`L in {12,16,20,24,32,40}` ×

- `(3,9)`
- `(4,8)`
- `(5,10)`

## 刻み対照 9 runs

`L in {12,24,32}` ×

- `(1,3)`
- `(2,4)`
- `(2,3)`

`den=40`

主系列・gcd系列は `den=L`。

総計は厳密に 99 runs。

走行長:

```text
t = 0 ... 4096
```

4097 状態を全保存する。
途中停止・延長禁止。

---

# 4. 初期条件

初期条件生成規則は前99走行と同じものを使う。

```text
U = exp(2*pi*i/L)
c_a = exp(-i*pi/(2L))
c_b = exp(+i*pi/(2L))
z_jk(0) = c_a U^(m_a*(k-j)) + c_b U^(m_b*(k-j)),  j<k
```

ただし初期値生成コードも、既存正本を関数単位で再記述しない。
既存の生成器を使用する場合は、その元ファイル SHA256 とコピー後 SHA256 を記録する。

各 run の `Z[0]` SHA256 または raw bytes hash を保存し、同一条件の再実行時に bit 同一性を検査できるようにする。

---

# 5. 数値環境を固定する

実行前に以下を保存する。

- macOS version
- CPU architecture
- Python version
- NumPy version
- SciPy version
- `np.show_config()` 全出力
- BLAS/LAPACK backend
- float64 / complex128 の確認
- git commit SHA
- 正本相互作用ソース SHA256
- runner SHA256
- generator SHA256
- manifest SHA256

9/5 正本系列の記録環境と異なる項目がある場合、勝手に続行しない。
`ENVIRONMENT_DIFF.md` に差分を書いて **実行前に停止し、木原へ報告する**。

少なくとも正本記録:

```text
Python 3.9.6
NumPy 2.0.2
macOS 26.3.1 arm64
BLAS/LAPACK: macOS Accelerate
state dtype: complex128
real dtype: float64
```

を照合対象とする。

---

# 6. 測定方法も正本へ戻す

状態生成中の基本測定は、9/5正本の `plane()` / `metrics()` と同一コードを使用する。
再記述・整理・式変形禁止。

記録する基本量:

- `Hperp_frac = H_perp / H`
- `H_total`
- `global_closure = abs(z @ z) / H`

初期平面 `p,q` は run の `z0` から一度だけ構成し、走行途中で更新しない。

さらに **全状態 Z(t) を保存するため、追加解析量は後から生軌道から計算する**。
実験本体へ新しい測定ロジックを混ぜない。

---

# 7. 生データ保存

各 run について必ず保存:

```text
states.npz        # Z[0..4096], complex128, 間引き禁止
metrics.csv       # t,Hperp_frac,H_total,global_closure
metadata.json
```

`metadata.json` には最低限:

- run_id
- L, M, ma, mb, den
- dtau
- T=4096
- source SHA256
- runner SHA256
- generator SHA256
- manifest SHA256
- environment record
- Z0 hash
- states.npz SHA256

を記録する。

99 run 完了後に全ファイル SHA256 一覧を生成する。

---

# 8. 実行前ゲート

99走行を開始する前に以下を全て PASS させる。

1. 正本元 SHA256 == 指定 SHA256
2. ローカルコピー SHA256 == 正本元 SHA256
3. `diff -u` が完全空
4. runner が SparseKernel / Lanczos / Krylov を import していない
5. `one_step` がローカル byte-copy 正本の関数である
6. dtype が float64 / complex128
7. manifest が厳密に99行
8. 4096 step 指定
9. output が旧99走行フォルダと別
10. 旧99走行の states / derived data を入力として参照していない

1件でも FAIL なら実行禁止。

---

# 9. 既存99走行系列との扱い

旧系列は削除しない。
ただし新しい再実験では入力・参照・比較対象にしない。

旧系列は研究上:

```text
INVALIDATED / DO NOT USE AS SCIENTIFIC RESULT
```

として隔離する。

新系列の結果を見る前に、旧系列との比較解析を行わない。

---

# 10. 99走行完了後

99本すべてについて

- 4097 states 存在
- metrics.csv 存在
- metadata.json 存在
- SHA256 存在
- NaN/Inf なし

を確認する。

完了後は

```text
CANONICAL DENSE 99-RUN DATA COLLECTION COMPLETE
```

とだけ報告する。

**この時点では科学解析を開始しない。**
木原によるデータ完全性確認後、全解析を新生軌道だけから最初からやり直す。

---

# 11. 次段階の相互作用変更実験について — 今回はまだ実施しない

本再実験終了後、ローカルに置いた正本相互作用ソースから **別コピー** を作って変更実験を行う。
正本再実験用コピーを直接編集してはならない。

現在の正本 `one_step` では、振幅正規化は次の1行の引数内部に埋め込まれている。

```python
H=H_of(np.exp(1j*np.angle(z)),A)
```

したがって、将来「振幅正規化だけを除去」する場合、単純なコメントアウトだけでは `H` が未定義になる。
最小変更は、元行をコメントとして残し、次の1行を追加する形になる。

```python
# H=H_of(np.exp(1j*np.angle(z)),A)   # stage2: phase-only normalization OFF
H=H_of(z,A)                           # amplitude retained
```

その後の

```python
H=(1j*np.imag(H)).astype(np.complex128,copy=False)
```

は **残す**。
これは実反対称生成子を得て実直交更新にするために必要であり、これをコメントアウトすると別の相互作用になる。

また、`H_of(z,A)` で複素 H 全体を作ってから虚部を取る処理は計算上は冗長に見えるが、これを直接 `K=A*Im(conj(z)z)` に書き換えると浮動小数点演算順序が変わり、今回のような軌道敏感系では別実験になる。

よって次段階も、まずは **1行だけの最小変更** とし、直接K実装への書換えはさらに別の対照実験に分離する。

---

# 12. 禁止事項の再掲

- 正本を整形しない
- 正本を再入力しない
- 正本関数を抜粋して別実装しない
- Sparse 化しない
- Lanczos 化しない
- tol を導入しない
- 高速化しない
- 解析結果を見て条件を変えない
- 旧99走行を新99走行の代替として使用しない
- 正本再実験用コピーを次段階の編集用ファイルとして直接変更しない

以上。
