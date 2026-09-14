# DATA_INVENTORY — Phase B 終端実現多様体次元測定

**版:** v0.1 / 2026-09-14  
**設計書:** `PhaseB_終端実現多様体次元測定_設計書_v1.2.md`  
**状態:** 読み取り専用棚卸し中。数値解析・再走行は未実施。

---

## 1. 棚卸し目的

Phase B 実装前に、保存済み自己無撞着力学データについて以下を確定する。

1. 正本データフォルダ
2. N の存在範囲
3. run / seed 数
4. step 範囲
5. final step の意味
6. 辺状態データの保存形式
7. Bob-star 抽出に必要な頂点―辺対応
8. 複素値 dtype
9. N ごとの欠損
10. 重複キャンペーン / 対照再走行の区別

---

# 2. 現時点で確認済みの主要キャンペーン

## A. N3_N40_stage123_sweep_20260905

保存先:

`次元の生成構造/電子の反跳実験/seed除去による準安定相追試/chatgpt追試/干渉保存力学_資格審査とシード無し系列_20260831/N3_N40_stage123_sweep_20260905/`

確認事項:

- N=3..40 の 38 系。
- 各 N について6分母 `N-2..N+2, 124`。
- 各走行 500 step。
- 状態 npz を全保存。
- 合計 228 走行 = 38 N × 6 分母。
- 主力学 `run_N3_N40_stage123_v1.py` は N=40 正本から力学無変更で全 N へ展開。
- Δτ=2π/N 系について step0 と final=step500 の複素平面図が存在。
- static parent は N ごとに1個生成され、生成 RNG 式は `default_rng(40260722+1000N+0)`。
- このキャンペーンは **seed アンサンブルではなく N×分母 sweep** であるため、Phase B の `seed×Bob` 主標本としてそのまま用いると seed 軸が不足する。

Phase B 候補用途:

- dry-run / Bob-star 抽出検証
- final-step 定義と辺データ schema の確認
- N=6〜40 の pilot
- 分母 N の canonical branch の静的比較

注意:

- 「6分母」を seed とみなしてはならない。
- H_B1 の主推定に必要な独立終端配置アンサンブルとして十分かは未確定。

---

## B. N3_N40_long10000_20260905

保存先:

`次元の生成構造/電子の反跳実験/seed除去による準安定相追試/chatgpt追試/干渉保存力学_資格審査とシード無し系列_20260831/N3_N40_long10000_20260905/`

確認事項:

- N=3..40 の全38系。
- den=N。
- 10,000 step。
- 全状態保存。
- `hm_N*_den_*_states_10000.npz` が N ごとに存在。
- RUN_METADATA:
  - `dtype_state = complex128`
  - `dtype_real = float64`
  - `steps = 10000`
  - `N_range = [3,40]`
  - `denominators = N`
- 500-step 正本スイープの延長系。
- README では step10000 を明示的に「終状態」として検査。
- N=4 と N>=7 は時計1回転写像の固定点残差が機械精度、N=6 はほぼ凍結、N=5 は例外的に緩和継続。
- N=3..40 の final-state table が存在。

Phase B 候補用途:

- final step=10000 の定義が最も明確な canonical long-run 系。
- Bob-star 抽出・位相基準 r の pilot に適する。
- κ_terminal と generic null の比較に適する。

制約:

- 現時点で確認できる範囲では各 N 1 canonical run であり、独立 seed 数は不足。
- Bob 数 N を標本とできるが、同一 run 内 Bob は非独立なので seed-cluster bootstrap は成立しない。
- よって **主仮説の seed-ensemble dimension 推定には単独では不十分**。pilot / 構造検査用と位置づける。

---

## C. 全N6_40_SSBバッチ_20260912（175走行 SSB 系）— 棚卸し確定 2026-09-14

保存先:

`次元の生成構造/電子の反跳実験/seed除去による準安定相追試/chatgpt追試/干渉保存力学_資格審査とシード無し系列_20260831/位相ロック全N確認_相対平衡_20260912/自発的対称性の破れ_シード依存縮退真空_20260912/全N6_40_SSBバッチ_20260912/`

確認事項:

- N=6..40 の35系 × 5 seed = 175走行。床 Z0＋振幅1e-8 の異なる向き微小シード→厳密 one_step→適応ロック停止。
- 保存物は `per_N/N{N}.txt`（seed ごとの lock_step・H⊥/H・振幅比・Δ・|Δ+2π/N|）と `ssb_batch_summary.csv`（N ごとの集約8列）のみ。
- ランナー `run_ssb_batch_N6_N40_20260912.py` に np.save / savez は存在しない。

**決定的事実:** 175走行 SSB バッチには診断値のみ保存され、**終端複素状態 Z は保存されていない**。したがって、**Phase B 主推定用の既存 seed×Bob 終端状態アンサンブルは存在しない**。存在するのは決定論的な再計算手順（`run_all.sh`）だけである。

この事実は Phase B の実験条件を変えるため、解析開始前にここに固定する。

---

# 3'. npz schema 確定（2026-09-14）

`N3_N40_long10000_20260905/results/hm_N8_den_8_states_10000.npz` の読み取り検査:

| key | shape | dtype |
|---|---|---|
| `Z` | (steps+1, M) | complex128 |
| `N` | () | int64 |
| `denominator` | () | int64 |
| `steps` | () | int64 |

M=N(N−1)/2（N=8 で 28 を確認）。**辺 index ↔ 頂点 pair の対応規約は未確定**（生成プログラムの辺列挙順を pilot で確認する）。

---

# 3''. 承認済み実行順序（2026-09-14 木原指示）

$$
\text{DATA\_INVENTORY更新}
\rightarrow
\text{long10000 で読出し pilot（N=6、次元主張なし）}
\rightarrow
\text{第一段：SSBバッチ忠実再走行＋終端Z保存（別ファイル、診断全一致を合格条件）}
\rightarrow
\text{5seed 有効標本評価（n\_raw=5N vs n\_unique / n\_eff）}
$$

- seed 拡張（5→50→100）は**未承認**。拡張判断は n_eff・bootstrap 安定性・d_est(K) の分散のみで行い、結果を4へ近づける目的では行わない。
- 第一段の合格条件: 175走行全てで lock step 完全一致・H⊥/H 一致・振幅比一致・Δ 診断一致・summary CSV 全列一致・物理式/シード生成/停止条件の差分ゼロ（浮動小数点は事前宣言の機械精度許容差）。
- 元ランナーは編集せず、`run_all_save_terminal_Z.py` 等の別ファイルとする。保存内容: N, seed_id, lock_step, Z_terminal ＋ 再現照合用の既存診断量。

---

# 3. final step 定義の現時点整理

| campaign | final step | 意味 |
|---|---:|---|
| N3_N40_stage123_sweep_20260905 | 500 | sweep の保存終了時。228走行すべて同一 step 数 |
| N3_N40_long10000_20260905 | 10000 | canonical den=N の長時間終状態。凍結性検査済み |

**重要:** final step はキャンペーン間で同じ物理的緩和度を意味しない。Phase B 主標本ではキャンペーンを混ぜず、キャンペーンごとに解析する必要がある。

---

# 4. 未確認・次に探索する対象

以下は設計査読で言及されたが、現時点の棚卸しではまだ正本パス・構造を確定していない。

1. **175走行 SSB 系**
   - seed アンサンブルを持つ可能性が高く、Phase B 主標本候補として最優先。
   - N ごとの seed 数、final step、状態保存形式を確認する。

2. **2000-step 延長走行**
   - final step の定義と N 範囲を確認する。

3. **stage123 sweep 対照再走行**
   - 正本の再現性確認用。主標本には混ぜない。

4. その他 seed ensemble
   - 同一 N で複数独立初期条件を持つ保存状態を探索する。

---

# 5. 現時点の Phase B 適格性判定

| campaign | N範囲 | final | 独立seed | Bob-star状態 | 主推定適格性 |
|---|---|---:|---:|---|---|
| stage123 sweep | 3..40 | 500 | 実質1/N（分母はseedでない） | npz全保存 | pilot向け |
| long10000 | 3..40 | 10000 | 実質1/N | npz全保存 | pilot向け、主推定にはseed不足 |
| 175-run SSB | 6..40 | 適応ロック（seed毎に可変、per_N txt に記録） | 5/N | **状態未保存（診断のみ）** | 忠実再走行＋終端Z保存（第一段）後に主標本候補 |

---

# 6. 次の棚卸しゲート

Phase B コード作成前に最低限、次を確定する。

- 175-run SSB 正本パス
- N ごとの独立 seed 数
- 各 run の final step
- npz の配列 key と shape
- 辺 index ↔ 頂点 pair の対応規約
- Bob-star が状態配列から一意に抽出できること
- N=6 の1ファイルを使った **読み取りのみ**の schema dry inspection

このゲートを通るまでは intrinsic dimension の計算を開始しない。
