# pass1_parents.py の使い方と N を増やす方法

## 目的

`pass1_parents.py` は、干渉保存力学の初期親状態を生成・検査・保存する基準プログラムである。

現在の標準実行では `N=3` から `N=16` までを生成する。

生成法は4種類ある。

- `mp_N<N>` : `original_engine.make_parent` による親
- `hm_N<N>` : 高対称な手作り等モジュラー親
- `ne_N<N>` : 非等モジュラー親
- `rb_N<N>` : 乱数均衡親（N>=5）

このうち、現在の高対称な自己無撞着初期値として使っているのは `hm_N<N>` である。

`hm` は `state_provider.equimodular(N)` で生成し、その後 `mp` 親と同じノルムへスケールして保存する。

## 現在の N 範囲

`pass1_parents.py` の現在のループは次の通り。

```python
for N in range(3,17):
```

したがって現在の生成対象は

```text
N = 3, 4, 5, ..., 16
```

である。

`hm` の保存先はプロジェクトルートから

```text
data/hm_N3/
data/hm_N4/
...
data/hm_N16/
```

である。

各 `hm_N<N>` フォルダには、少なくとも次の親データが保存される。

```text
parent_v.npz
parent_v.csv
parent_checks.json
```

## 実行方法

`pass1_parents.py` は `program/` 内の関連モジュールを読み込むため、プロジェクト構造を保ったまま実行する。

例:

```bash
cd <project-root>
python3 program/pass1_parents.py
```

必要な主要ファイル:

```text
program/pass1_parents.py
program/original_engine.py
program/common.py
program/state_provider.py
program/interference_dynamics.py
```

## N を増やす方法

N の生成範囲は `pass1_parents.py` のこの1行で決まっている。

```python
for N in range(3,17):
```

例えば N=20 まで生成対象を増やすなら、

```python
for N in range(3,21):
```

とする。

N=40 までなら、

```python
for N in range(3,41):
```

となる。

ただし、この変更では N=3 から指定上限までを再生成する。

既存 N を再生成せず、新しい N だけを生成したい場合は、例えば N=17〜20 なら一時的に

```python
for N in range(17,21):
```

とする。

単独の N=40 だけを生成したい場合は、

```python
for N in [40]:
```

とする。

## 重要: hm の生成法を別プログラムで置き換えない

高対称 `hm` 系列を追加するときは、ランダム複素ベクトルを作る独立プログラムを使わない。

基準は `pass1_parents.py` 内の

```python
v = sp.equimodular(N)
v = v*NORM[N]/np.linalg.norm(v)
rows.append(save(
    f'hm_N{N}',
    N,
    v,
    'handmade_equimodular_'+(
        'Z3' if N==3 else
        '1factor' if N%2==0 else
        'distance_classes'
    ),
    dict()
))
```

である。

したがって N を追加するときも、同じ `state_provider.equimodular(N)` を通す。

偶数 N は 1-factor 分解、奇数 N は距離クラス、N=3 は Z3 という現在の高対称構成をそのまま継承する。

## 保存前の検査

`save()` では親状態について少なくとも以下を検査・記録する。

- 自己無撞着残差
- 全体二乗閉塞 `|Σz²|/H`
- `mu != 0`
- 局所閉塞
- 振幅範囲
- Gram/Takagi rank
- 新しい Hermitian H に対する自己無撞着量
- モノドロミー予測

受け入れ条件は現在、

```text
residual < 1e-10
global_closure < 1e-12
|mu| > 1e-6
```

である。

条件を満たさない親は `ABORT` する。

## 計算量上の注意

N を大きくすると辺数は

```text
M = N(N-1)/2
```

で増える。

さらに `pass1_parents.py` は `mp` 親の生成、ヤコビアン、固有値計算、モノドロミー解析も同時に行うため、N=3〜16の延長をそのまま大きな N まで行うと計算量とメモリ使用量が急増する。

したがって大きな N を追加する場合でも、まずは `pass1_parents.py` の生成規約を正本とし、勝手に別方式の初期値へ置き換えないこと。

特に `hm_N40` 等を作る目的で、単純なランダム複素ベクトル生成プログラムを代用してはいけない。それは `hm` 系列とは異なる初期条件である。

## 基準

高対称 hm 初期値を増やす際の基準は、

```text
pass1_parents.py
    -> state_provider.equimodular(N)
    -> mp 親と同じ norm にスケール
    -> save()
    -> data/hm_N<N>/parent_v.npz
```

である。

この系列を比較実験の標準初期条件として扱う。
