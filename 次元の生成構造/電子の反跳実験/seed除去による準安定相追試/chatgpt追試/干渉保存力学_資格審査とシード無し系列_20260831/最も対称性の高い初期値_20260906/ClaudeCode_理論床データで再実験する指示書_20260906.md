# Claude Code 指示書：高対称理論床データで失敗データを置換し再実験する

作成日: 2026-09-06

## 0. 目的

`最も対称性の高い初期値_20260906` で作成した旧 v2/v3 の「重心ゼロ親」は、N>=4 で自己無撞着固定点残差が O(1) となり、走行時に step 1 で Hperp/H が巨視的値へ跳ぶ即時緩和を起こした。

この失敗データを廃棄対象とし、ChatGPT が別途導出した **高対称・重心ゼロ・二乗ゼロ閉塞・現行力学の相対平衡**を同時に満たす「理論床」データで置換し、同じ物理正本プログラムを使って500 step実験をやり直す。

理論導出の正本:

`最も対称性の高い初期値_20260906/高対称理論床初期値の導出_N3_N40_20260906.md`

この文書を必ず先に読み、旧 `最も対称性の高い初期値の作り方.md` の v2/v3 構成を理論根拠として再利用しないこと。

---

## 1. 絶対に守る条件

### 1.1 物理プログラムは変更しない

物理正本は以下のみを使用する。

`N3_N40_stage123_sweep_20260905/run_N3_N40_stage123_v1.py`

正本 SHA256:

`1abf2353fee2e4f56f05e7a6f149fd086885136beb61ab571b48a56b09691567`

走行前に必ず SHA256 を照合すること。

変更を許可するのは、既存ラッパーと同じく以下だけ。

1. `PARENT_DIR`
2. `OUT`
3. 必要なら試走時だけ `for N in range(3,41)` の上限

`one_step`, `H_of`, `plane`, `metrics`, denominator、STEPS、dtype、保存形式、図化処理には一切手を入れない。

### 1.2 旧正常系を壊さない

以下は対照データであり、絶対に上書きしない。

`N3_N40_stage123_sweep_20260905/parents/`

`N3_N40_stage123_sweep_20260905/results/`

旧 `make_parent` 親は、同じ新走行系へ戻すと正常な指数インフレーションを再現することが ChatGPT ローカル再試験で確認済みである。

### 1.3 上書きしてよい対象

今回失敗した実験側のステージング親と結果のみ、理論床版へ更新してよい。

対象フォルダ:

`対称親v2_500step走行_20260906/`

特に:

- `parents_symmetric_staged/`
- `results/`

は、今回の理論床実験用として作り直してよい。

ただし削除前に、現在の内容が v3 失敗系列であることを manifest / SHA で確認すること。誤って旧正常系を消さないこと。

---

## 2. 理論床データの定義

旧 v2/v3 のように、

- Σz=0
- Σz^2=0
- ||z||=1

だけを先に満たして「高対称」と判断してはいけない。

主条件は現行力学に対する相対平衡:

`H(z) z = lambda z`

である。

その上で

- `||z||^2 = 1`
- `Σ z = 0`
- `Σ z^2 = 0`
- 固定点残差が機械精度級

を同時に満たすこと。

理論床系列は、導出 md に記載した巡回 Fourier 対称 Ansatz を基礎とする。

代表形:

`z_ij = r_d * exp(i * 2*pi*(i+j)/N)`

ここで `d` は巡回距離クラス。

N ごとに距離クラス部分空間へ生成子を落とし、有限次元固有値問題

`Q^(N) c = lambda c`

として解く。

`make_parent` の固定点反復は、この理論床生成には使用しない。

---

## 3. 最初に N=3..6 を再構成して検証する

いきなり N=3..40 を走らせないこと。

まず N=3,4,5,6 の4本だけ理論床親を生成し、次の全項目を検証する。

- dtype = complex128
- M = N(N-1)/2
- `abs(vdot(z,z).real - 1) < 1e-13`
- `abs(sum(z)) < 1e-12`
- `abs(sum(z*z)) < 1e-12`
- 固定点残差 `< 1e-12` を目標
- 位相・振幅が導出 md と一致
- 辺順序が正本 `np.triu_indices(N,k=1)` と一致

特に N=4,5 は、旧 v3 の誤った重心ゼロ親と混同しないこと。

旧 v3 失敗値の例:

- N=4 residual = 0.5
- N=5 residual = 0.707106781...
- N=6 residual = 1.545603...

これらのような O(1) 残差が出た場合は、その時点で停止する。走行してはいけない。

---

## 4. N=3..6 の500 step予備走行

理論床親4本を `parents_symmetric_staged/` に、正本が読む旧名

`parent_static_N{N:05d}_makeparent_20260905.npz`

で配置する。

NPZ には少なくとも正本が必要とする `Z0` を保存する。

既存ラッパー

`対称親v2_500step走行_20260906/wrapper_run_symmetric500_v1.py`

は、物理正本の SHA 照合とパス置換の仕組み自体は利用してよい。ただし旧 v2/v3 親を再生成・再コピーする `stage_symmetric()` 系の処理は使わないこと。

試走は N=3..6 のみ。

期待する挙動:

- step 0 は Hperp/H が数値床
- step 1 で O(1) へ跳ばない
- N=4,5 では数値床から数十〜数百 step をかけた指数ランプが現れる可能性がある
- N=3,6 は少なくとも500 stepで床に留まる可能性がある

この期待と異なっても結果を改変しない。まず residual・親 SHA・辺順序・正本 SHA を確認する。

---

## 5. 旧正常系との対照

必ず旧正常系と比較する。

旧親:

`N3_N40_stage123_sweep_20260905/parents/parent_static_N00003_makeparent_20260905.npz`

〜

`parent_static_N00040_makeparent_20260905.npz`

旧結果:

`N3_N40_stage123_sweep_20260905/results/timeseries_64bit_with124_N3_N40.csv`

`N3_N40_stage123_sweep_20260905/results/summary_64bit_with124_N3_N40.csv`

比較項目:

- `Hperp/H` step0
- step1
- onset `Hperp/H > 0.05`
- 初期指数成長率
- final / max
- H_total
- global_closure

理論床親では旧 make_parent 親と onset が一致する必要はない。

比較目的は「同一物理系で、初期親の枝が違うと安定性がどう変わるか」を測ることである。

---

## 6. N=3..6 が合格した後だけ N=3..40 へ拡張

N=3..6 で以下が確認できた場合のみ、N=7..40 を生成する。

合格条件:

1. 全親の固定点残差が機械精度級
2. Σz=0
3. Σz^2=0
4. ||z||=1
5. step1 即時 O(1) ジャンプなし
6. 物理正本 SHA 一致
7. 親 manifest と SHA256 台帳を保存

N=7..40 についても1本でも residual が大きければ、その N で停止し、理論床導出を再確認する。

「条件を満たさないが走らせてみる」は禁止。

---

## 7. 出力物

同じ実験フォルダ `対称親v2_500step走行_20260906/` 内に、最低限以下を残す。

### 親

`parents_symmetric_staged/`

### 親検証台帳

`theoretical_floor_parent_manifest.json`

各 N について:

- N
- M
- lambda
- fixed_point_residual
- norm2
- centroid_abs
- square_closure_abs
- sha256
- construction / branch 名

### 結果

`results/`

正本が生成する:

- `hm_N*_den_*_states_500.npz`
- `timeseries_64bit_with124_N3_N40.csv`
- `summary_64bit_with124_N3_N40.csv`
- `fig_Hperp_denominator_controls_with_124_N3_N40_stage123.png`
- `RUN_METADATA_N3_N40_stage123.json`

### 追加比較

- `compare_theoretical_floor_vs_old_makeparent.csv`
- `fig_compare_theoretical_floor_vs_old_makeparent.png`
- `再実験結果メモ_理論床_N3_N40.md`

---

## 8. 禁止事項

- `make_parent` を理論床生成器として使わない
- residual O(1) の親を「高対称」と呼ばない
- 重心ゼロだけで合格判定しない
- 物理正本を変更しない
- 結果を合わせるため振幅・位相・正規化を後調整しない
- 旧正常系を上書きしない
- N=3..6 の事前検証を飛ばして N=40 まで一括走行しない
- 異常値を「新しい発見」として先に解釈しない。まず入力資格を確認する

---

## 9. 今回の重要な教訓

旧 v3 では、N=4..6 の重心ゼロ親が

- N=4 residual = 0.5
- N=5 residual = 0.7071
- N=6 residual = 1.5456

であり、同じ物理正本に投入すると step1 で巨大な Hperp 成分を直接生成した。

一方、旧 make_parent 親を同じ新走行系へ戻すと正常な指数ランプが再現した。

したがって今回の再実験では、

**「高対称性」より先に「力学の相対平衡であること」を親の資格条件とする。**

これを破るデータは、重心ゼロ・二乗閉塞を満たしていても親として不採用とする。

