# 実験——最新版の `make_parent` だけを振幅込み K で自己無撞着にした版（N=5、他は最新版と同一）

**作成日:** 2026-08-28　**ベース:** 木原最新版 `N5_linear124_all3fix_seedless_parentnorm_removed_40000_20260828`。
**変更:** `program/original_engine.py` の `make_parent`（と補助関数 `_adjacency`, `_K_amplitude_aware`, `_selfconsistency_residual`）のみ。`run_amplitude_only_fix.py` は最新版と同一（diff なし）。差分は `results/diff_engine_vs_latest.patch`（make_parent の区間 1 hunk）。

make_parent：K_amp(v) = Im(v̄ᵢvⱼ) に対して K_amp(v)v = −iσ_max v（旧実装と同じカイラリティ）を、位相・振幅とも未知数として混合反復で解く。正規化なし（スケールは固有値ソルバの ‖u‖=1）。残差はスケール不変、tol=1e-10、iters=2000、restarts=20、全失敗なら abort。

## 結果

| 版 | 親残差 | H_total | baseline onset(5%) / max H⊥/H / final PR/M | treatment onset(5%) / max H⊥/H / final PR/M / 振幅範囲 |
|---|---|---|---|---|
| 最新版 無変更再実行 | 4.87e-01 | 0.8596 | 334 / 0.716 / 1.000 | 76 / 0.943 / 0.296 / 3.9e-3..0.655 |
| (a) make_parent 正規化復元 | 1.82e-09 | 1.0000 | 196 / 0.624 / 1.000 | 65 / 0.991 / 0.307 / 1.9e-3..0.695 |
| **(b) make_parent 自己無撞着（本実験）** | 2.63e-11 | 1.0000 | 7 / 0.708 / 1.000 | **None / 4e-16 / 0.486 / 3.8e-16..0.564** |

- 自己無撞着状態は 1 リスタート・100 反復で収束（残差 2.6e-11）。中身は頂点 4 を切り離した K₄（6 辺）で、4 辺は振幅 0、K_amp(v) のスペクトルは ±0.5 のみ。
- treatment：**40000 step 不動**（H⊥/H ≈ 10⁻¹⁶）。安定な固定点。潜伏も急拡大もない。
- baseline：この v は位相のみ K の固有モードではないので step 7 で離脱し等分配へ（旧力学の性質）。

## ファイル
`program/`、`data/`、`figures/`、`results/`（diff、run.log、run_progress.log）、`run_all.sh`、`SHA256SUMS.txt`

---

## 自己引用・関連過去資料（2026-09-07 追記）

本実験を、`make_parent` の内部構造について既に得られていた結果と接続するため、以下を自己引用として明記する。ここでは、**自己無撞着性から何が入力なしに自然導出されているか**を区別して参照する。

1. **`parent_state_structure_check_v1.md`**（2026-07-27、第8論文補助解析）  
   Drive ID: `1VTW453qZ3WTH7oHVISGgILQkPfHxnOYE`  
   `make_parent` 出力を時間発展なしで直接検査した正本級の構造解析。N=5 と N=40 で、全 M 成分が非零のまま実 rank 2 に集約され、
   $$
   \|\Re v\|=\|\Im v\|=1/\sqrt{2},\qquad \Re v\cdot\Im v\simeq0,
   $$
   したがって
   $$
   v^T v=\|\Re v\|^2-\|\Im v\|^2+2i\,\Re v\cdot\Im v\simeq0
   $$
   が確認された。これは **二乗ゼロ閉塞を親生成条件として課さず、実反対称生成子の円偏波固有モード構造から出現した**ことの直接確認である。

2. **`make_parent(倍音対応).md`**（2026-08-04）  
   Drive ID: `1VjrTJKafFR4vYNbOYSykuTqu8ZBxrBFu`  
   現行 `make_parent` の原理を倍音段へ拡張する際の設計書。設計原則として「**閉塞は課さず、構造から恒等的に出す**」「**初期値は選ばず、位相場 → 生成子 → 固有平面 → 位相場の不動点として解かせる**」を明記している。さらにテスト v2 では段閉塞・総閉塞・rank 2・等ノルム直交・σ_max 平面占有を確認し、自己無撞着閉包が離散的な不動点族を成すことを記録している。

3. **`ch1_static_parents_en.md`**（2026-09-05、N=3..40 再現論文 Chapter 1）  
   Drive ID: `1aYF0JPeKhu-4rxa5Fe-t2TS-R3eBUIMG`  
   `make_parent` の生成アルゴリズムを数式とコード行に対応させた再現性文書。位相差生成子
   $$
   K_{ef}=\sin(\theta_f-\theta_e)
   $$
   から小空間固有対を解き、固有ベクトルの位相を再投入して固定点を求めること、収束条件
   $$
   \|iKv-\mu v\|<\mathrm{tol}
   $$
   を明示している。90°位相差や二乗ゼロ閉塞を手で与えていないことを確認する基準文書である。

4. **`project_six_axis_dynamics_framework.md`**（2026-08-05 までの統合研究記録）  
   Drive ID: `1fFayrcpZW0kbp5AL95OOxzUrmBhxvzBN`  
   `make_parent` が落ちる自己無撞着固定点族の枝・安定性を追跡した統合記録。N=3..12 に一様族が存在すること、安定枝と `make_parent` 規約枝側の不安定性が区別されることを記録している。したがって「自己無撞着解は一種類ではなく、離散的な閉包族と安定性構造を持つ」という後続解析の参照先とする。

5. **`発見台帳_干渉保存力学_20260905起点.md`**  
   Drive ID: `1OJztVgEkyJuJlkMoMs7zwKTBkkN8wnty`  
   旧親系列の時間発展・相対平衡・位相ギャップを横断的に整理した発見台帳。増幅期に `{90°, 90°, 180°}` 型の位相ギャップが凍結した剛体回転として現れる系列が記録されている。これは、90°構造が初期値へ手で埋め込まれたものではなく、自己無撞着親とその後の力学に繰り返し現れることを確認する補助資料である。

### 現時点での因果関係の整理

過去資料と本実験を合わせると、少なくとも `make_parent` 系については、次の向きで整理するのが適切である。

$$
\text{自己無撞着固定点}
\longrightarrow
\text{実反対称 }K\text{ の円偏波固有平面}
\longrightarrow
\|\Re v\|=\|\Im v\|,\ \Re v\perp\Im v
\longrightarrow
\sum_e v_e^2=0.
$$

したがって、**二乗ゼロ閉塞はこの系列では入力拘束ではなく、自己無撞着円偏波構造からの導出結果として扱う**。

一方、`U^N=I` については、自己無撞着閉包が離散的不動点族・位相ロック・有限巡回構造を選ぶという数値証拠は蓄積しているが、**自己無撞着性のみから `U^N=I` を一般に解析導出した、とまでは本追記では主張しない**。ここは今後の導出課題として分離する。
