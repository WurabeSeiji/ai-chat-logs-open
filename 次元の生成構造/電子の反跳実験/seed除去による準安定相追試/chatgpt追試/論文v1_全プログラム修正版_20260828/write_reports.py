# -*- coding: utf-8 -*-
"""各パッケージ配下に 差異分析_修正版vs原本_20260828.md を、最上位に 全プログラム修正版_差異分析_20260828.md を生成する。"""
import os, re, json, collections, numpy as np, pandas as pd
H=os.path.dirname(os.path.abspath(__file__)); FX=os.path.join(H,"fixed"); O=os.path.join(os.path.dirname(H),"論文v1_全再現テスト_20260828","original")
cmp=json.load(open(os.path.join(H,"results","comparison.json"))); phys=open(os.path.join(H,"results","physics_comparison.md"),encoding="utf-8").read()
diff=open(os.path.join(H,"results","fix_patches.diff"),encoding="utf-8").read()
hunks=collections.Counter(); cur=None
for line in diff.splitlines():
    if line.startswith("diff -ru"): cur=line.split()[-1].replace("fixed/","")
    elif line.startswith("@@") and cur: hunks[cur]+=1
def sec(title_re):
    m=re.search(r"(^## \d+\. .*?%s.*?$)(.*?)(?=^## \d+\. |\Z)"%title_re,phys,re.S|re.M); return (m.group(1)+m.group(2)).strip() if m else ""
def csv(base,p): return pd.read_csv(os.path.join(base,p))
def js(base,p): return json.load(open(os.path.join(base,p)))
FIXES="""## 適用した修正（全パッケージ共通）

| # | 旧プログラム | 修正版 |
|---|---|---|
| FIX1 | `make_parent` で `v = v/‖v‖`（親の振幅正規化） | 除去。親の固有モード反復は位相のみの K のまま（参照 `N5_linear124_all3fix` と同じ） |
| FIX2 | 初期化 `Z = v + δg`、`Z /= ‖Z‖`（外部 seed と正規化）、乱数初期状態の正規化 | `Z = v`（seedless）、正規化なし |
| FIX3 | Cayley 変換 (I−γK)⁻¹(I+γK)、γ = tan(π/144) | 厳密な線形回転 exp(ANGLE·K)、ANGLE = 2π/144（旧と同じ名目刻み、K/σ 正規化ブランチは ANGLE/σ） |
| FIX4 | 生成子が位相のみ K_ij = sin(θ_j−θ_i)（`set_theta(np.angle(Z))`） | 振幅込み K_ij = Im(z̄_i z_j) = |z_i||z_j| sin(θ_j−θ_i)（`set_state(Z)`：c=Re z, s=Im z） |

検証（`verify_fixes_all.py`、14 エンジン全て OK）：`cayley_step` 不在、‖log U − ANGLE·K‖ ≈ 5×10⁻¹⁵、U 実直交、K(2z) = 4K(z)、参照 all3fix の `dense_K_amplitude_aware` と 3×10⁻¹⁷ で一致、低ランク実装 vs 密行列 300 step 軌道偏差 9×10⁻¹¹（`--validate`）。
"""
PK={
"K_sigma_normalization_artifact_test_N4_N5_20260826":("K/σ 正規化の有無で成長率が変わるかの対照（N=4,5）","""- 原本：raw K で N5 成長率 0.172094/step（onset 217）、normalized 0.049223（onset 757）。指数成長域が明瞭。
- 修正版：Z₀ = v（seed なし）で onset は step 1（親が修正後の力学の固定点でない）。N=4 は raw/normalized とも親平面を離れない（H⊥ = 0、onset なし）。N=5 raw の「成長率」0.000669 は指数域が無いための無意味な fit。final H⊥/H は raw 0.0736、normalized 0.3557。
- 保存則は成立（H_total = c² = 0.8596 が 1e-12 で一定、H∥+H⊥ = H_total 1e-16）。
- **§11「正規化はアーティファクトでない」の比較そのものが成立しない**（比較すべき指数成長が両ブランチに無い）。"""),
"N3_N4_complex_simplex_complete_analysis_20260826":("N=3, 4 の完全解析（等分配・120° 定理の数値側）","""- 原本：N=4 は等分配（|z|² = 1/6 ± 3e-6）、H⊥ 0.167、onset 260。N=3 は非等分配。
- 修正版：**N=3, 4 とも親平面を一切離れない**（H⊥ final 0 / 2e-16、onset なし）。|z|² は親のまま（N4: 0.083/0.167 の 2 クラス）。振幅込み K では親が（この 2 つの N で）不変平面を保つ。
- 120° の定理（§21）は数学として残るが、動力学がその状態に到達するという原本の数値主張は修正版では出ない。"""),
"N5_complex_simplex_complete_analysis_20260826":("N=5 の位相刻み検査・複素 simplex 構造（記事図 7 の根拠パッケージ）","""- 修正版で再生成できたのは `run_N5_physical_phase_step_test.py` の出力（3 csv・4 png）と `--validate`。`N5_all_steps_a_b_a2_b2_ab.csv`・4 群パターン・時間分離 milestones・インフォグラフィック等 16 件は生成プログラム未同梱のため**原本のまま（旧エンジンのデータ）**。`analyze_N5_complex_simplex_structure.py`／`plot_N5_inflation_vs_ordering.py` はこれら旧 csv を読むので、修正版でも旧結果を描いているだけ。
- 修正版の N=5（phase_by_edge より）：10 本の |z|² は等分配せず、1 本が 0.2〜0.43 を周期 ≈ 700 step で振動、他 9 本は 0.02〜0.08。スペクトルエントロピー S/lnM = 0.87（原本 1.000000）。H_total = 0.8596（= c²）。有理ロック指標は原本と同程度に「無し」。
- **記事図 7・2627/4923 step の 3+3+2+2 時間分離は、修正版では検証不能（プログラム無し）かつ、等分配が起きないので前提が崩れる。**"""),
"N6_N7_complex_simplex_complete_analysis_20260826":("N=6, 7 完全解析","""- 原本：等分配（N7 は |z|² = 1/21 が 10 桁）、H⊥ 0.79 / 0.26、onset 234 / 271。
- 修正版：onset は step 1、|z|² は 0.002〜0.27（N6）、5e-5〜0.29（N7）で等分配せず、H⊥/H_total は 0.69 / 0.84。H_total は親の c²（N6 0.62 相当の残差、N7 2.09）。rank = N−1 は保たれる。"""),
"N8_N9_complex_simplex_complete_analysis_20260826":("N=8, 9 完全解析","""- 原本：等分配（10 桁）、H⊥ 0.10 / 0.15。
- 修正版：onset step 1、|z|² は 0.006〜0.40 / 0.003〜0.27、H⊥/H_total ≈ 0.92 / 0.94（ほぼ全成分が親平面外へ）。rank N−1、|ZᵀZ| ≈ 1e-13 保存。"""),
"N10_N11_complex_simplex_complete_analysis_20260826":("N=10, 11 完全解析","""- 原本：等分配（10 桁）、H⊥ 0.10 / 0.14。
- 修正版：onset step 1、|z|² 相対ばらつき ≈ 2（等分配なし）、H⊥/H_total ≈ 0.96 / 0.99。親残差（正規化前提の残差式）14.0 / 20.7。"""),
"N12_N13_complex_simplex_complete_analysis_20260826":("N=12, 13 完全解析","""- 原本：等分配（10 桁）、H⊥ 0.098 / 0.070。
- 修正版：onset step 1、|z|² 0.0002〜0.19 / 0.0002〜0.20、H⊥/H_total ≈ 0.96 / 0.95。"""),
"N14_N15_complex_simplex_complete_analysis_20260826":("N=14, 15 完全解析","""- 原本：等分配（10 桁）、H⊥ 0.084 / 0.144。
- 修正版：onset step 1、|z|² 0.001〜0.16 / 0.0016〜0.10、H⊥/H_total ≈ 0.97 / 0.99。"""),
"N16_complex_simplex_complete_analysis_20260826":("N=16 の物理走行と位相構造","""- 原本：H_total = 1、等分配 |z|² = 1/120（15 桁）、H⊥ final 0.127。
- 修正版：親のノルム c² = 3.62（正規化除去で親のスケールが N に依存）、親残差 73.7（残差式が単位ノルム前提のため）、onset step 1、指数成長域なし（fit None）、final |z|² 0.001〜0.13、H⊥ 3.53/3.62。`analyze_and_plot_N16.py` は修正版 decompactification の csv をコピーして実行（5 図再生成）。位相構造 json・時間分離 csv・selected_steps 等 7 件は生成プログラム未同梱。"""),
"complex_simplex_decompactification_N5_N16_20260826":("記事図 1（10⁻³² → 31 桁）の源：N5/N16 の親平面外成分の対数時系列","""- 原本：N5 は H⊥ 1.1e-31 → 0.40（成長率 0.08626/step、onset 242）、N16 は 3.0e-32 → 0.127（0.23608/step、onset 94）。潜伏 → 指数増大 → 飽和。
- 修正版：H⊥(0) は丸め床（3e-32 / 6e-32）だが **step 1 で 10⁻² 台へ跳ぶ**（潜伏なし）。指数成長域が無く fit は None（報告行の None 安全化のみ追加）。N5 は H⊥ 0.43 を上限に周期 ≈ 700 step で振動、final 0.063。N16 は H⊥ 3.53/3.62。
- **記事図 1 の「31 桁のインフレーション」は修正版では存在しない。** 理由は親が修正後の力学の固定点でないこと（§決定事項 1）。"""),
"N5_gamma_continuum_test_bundle_20260825":("§12「144 は時間刻み」：刻み n_den = 144〜2304 の連続極限","""- 原本：rate_per_radian が 1.131 → 1.158 に単調収束（g(n) 則）。
- 修正版：指数成長域が無いので rate_per_radian の fit は無意味（17.97, 21.16, 0.094, …）。一方 2304 を基準にした RMSE は 0.144 → 0.095 → 0.050 → 0.018 と刻みに対し一次で収束しており、**線形回転でも軌道は刻みの連続極限を持つ**（時計の役割は保たれる）。"""),
"N5_dynamics_followup_theorems_and_stability_20260826":("記事図 2〜6・8 の源：tol 掃引・Floquet・エントロピー・モジュライ","""- tol 掃引（記事図 3）：親残差は 4 段階とも 0.4869（= σ·c·|1−c²|、正規化除去後は tol 判定が決して真にならない）。onset は全て step 1、成長率 NaN。**onset–残差則（傾き 11.616）は定義できない。**
- Floquet（記事図 4・5）：max|μ| = 1.001777（原本 1.090087）、第 2 対 1.000605、Σln|μ| = −3e-4。回転系で親はほぼ中立（固定点欠陥 2.0e-3）。**三重整合は成立しない。**
- pump depletion（記事図 2）：保存則は成立（H∥+H⊥ = H_total = 0.8596 が 1e-16）。H⊥ は最大 0.43（step 2567）、final 0.063。
- スペクトルエントロピー（記事図 6）：`analyze_followup.py` は `N5_all_steps_a_b_a2_b2_ab.csv`（未再生成・旧データ）を読むため、**修正版フォルダの entropy csv・図は旧データのまま**（無効）。修正版の値は `results/figures/cmp5_amplitudes_N5.png` で phase_by_edge から再計算：S/lnM = 0.87。
- モジュライ seed 掃引（記事図 8）：3+3+2+2 は出ない（counts 4;3;2;1、6;2;1;1 等、family moduli 0.198 / 0.05）。20 seed 版も同様。
- `time_reparameterization` は原本の Cayley 位相増分式を ANGLE·σ に置換して計算。"""),
"N3_N16_partial_zero_closure_analysis_20260826":("N=3〜16 の部分ゼロ閉包の数え上げ（最終状態の再解析、動力学なし）","""- 入力を修正版の最終状態（N=3,4,6〜16）に差し替えて再実行（N=5 だけは `N5_all_steps_a_b_a2_b2_ab.csv` が未再生成のため旧データ）。
- 修正版の最終状態は等モジュラーでないため、原本で見えた 2 辺・3 辺閉包の構造（exact cover 等）は前提を失う。出力は `out/` に生成（原本の出力先はパッケージ自身を rmtree するため変更）。"""),
"N3_N16_nontrivial_zero_closure_analysis_20260826":("非自明閉包の評価（動力学なし）","""- 同梱スクリプトは `reproduced_small_subset_classification.csv` 1 本を書くのみで、原本の 25 ファイル（MITM 表・評価表・図）は生成プログラム未同梱。SOURCE_N* は修正版の最終状態に差し替え（N=5 は旧）。"""),
"N14_N16_complete_nontrivial_zero_closure_search_20260826":("N=14〜16 の全次数閉包探索（exact k=2..4、MITM k=5,6、seeded ヒューリスティック）","""- 入力を修正版の最終状態に差し替え。exact k=2..4 の最良残差：N14 1.4e-2 / 2.0e-3 / 1.4e-4（原本 1.6e-4 / 7.2e-3 / 8.9e-6）、N15 1.3e-2 / 1.0e-3 / 1.8e-4、N16 1.8e-2 / 2.1e-4 / 3.1e-4。
- MITM k=6：N14 2.4e-5（原本の quasi-closure 9.25e-7 は**消える**）、N15 3.4e-6、N16 2.7e-6。1e-6 未満は無し。
- 焼きなましは steps/seed 未記録のため未実行（原本と同じ）。"""),
}
INDEX=[]
for pkg,(purpose,interp) in PK.items():
    files=[r for r in cmp if r["package"]==pkg]; regen=[r for r in files if not r["verdict"].startswith(("NOT_REGENERATED","MISSING"))]; notre=[r["file"] for r in files if r["verdict"].startswith("NOT_REGENERATED")]
    changed=[(f,n) for f,n in hunks.items() if f.startswith(pkg+"/")]
    L=[f"# 差異分析：修正版（FIX1-4）vs 原本 —— `{pkg}`","",f"**作成日:** 2026-08-28　**目的:** {purpose}","",
       "原本 = `論文v1_全再現テスト_20260828/original/`（旧エンジン：親正規化・seed・Cayley・位相のみ K）。修正版 = 本フォルダ（振幅込み K・線形回転・正規化なし）。ANGLE = 2π/144。","",
       "## このパッケージで変更したファイル（`results/fix_patches.diff`）","", "| ファイル | hunk 数 |","|---|---|"]+[f"| `{f.split('/',1)[1]}` | {n} |" for f,n in sorted(changed)]+["",
       "## 結果の差異（要点）","",interp,"",
       f"## ファイル別突合（再生成 {len(regen)} / 全 {len(files)}）","","| file | verdict | max|Δ| / 画素不一致率 | 備考 |","|---|---|---|---|"]
    for r in regen:
        v=r.get("max_abs_diff",r.get("pixel_mismatch_fraction","")); v=f"{v:.3g}" if isinstance(v,float) else v; note=[]
        if r.get("shape_o") and r.get("shape_o")!=r.get("shape_r"): note.append(f"shape {r['shape_o']}→{r['shape_r']}")
        if r.get("str_mismatch"): note.append(f"文字列不一致 {r['str_mismatch']}/{r['str_cells']}")
        L.append(f"| {r['file']} | {r['verdict']} | {v} | {'; '.join(note)} |")
    if notre: L+=["",f"## 未再生成（生成プログラム未同梱、原本のまま）{len(notre)} 件","","- "+"、".join(f"`{f}`" for f in notre)]
    L+=["","## 参照","","- 全体分析：`../../全プログラム修正版_差異分析_20260828.md`","- 物理量突合：`../../results/physics_comparison.md`、比較図：`../../results/figures/`","- 修正の検証：`../../results/verify_fixes_all.json`"]
    p=os.path.join(FX,pkg,"差異分析_修正版vs原本_20260828.md"); open(p,"w",encoding="utf-8").write("\n".join(L)+"\n"); INDEX.append((pkg,purpose,len(regen),len(files),len(notre)))
# ---------------------------------------------------------------- 最上位
cnt=collections.Counter(r["verdict"].split("(")[0] for r in cmp)
S=csv(FX,"N5_complex_simplex_complete_analysis_20260826/N5_phase_by_edge_5000steps.csv"); amp2=S.pivot(index="step",columns="edge_index",values="amplitude")**2; pr=amp2.div(amp2.sum(axis=1),axis=0); ent=float((-(pr*np.log(pr.where(pr>0,1))).sum(axis=1)/np.log(10)).iloc[-1])
top=f"""# 全プログラム修正版（振幅問題・回転問題を全パッケージで修正）——原本との差異分析

**作成日:** 2026-08-28　**フォルダ:** `chatgpt追試/論文v1_全プログラム修正版_20260828/`
**方針（木原指示）:** 振幅問題と回転問題を全てのプログラムで修正し、全実験を再実行・図化し直し、差異を分析してから修正方針を決める。本書はその差異分析。修正方針の決定はまだ行っていない。

対象 = 論文 v1（DOI 10.5281/zenodo.22112009）と note 記事 n07c3e4c97e3a が引用する 15 パッケージ（Zenodo 14 zip ＋ decompactification zip）。原本は `論文v1_全再現テスト_20260828/original/`。各パッケージ配下に `差異分析_修正版vs原本_20260828.md` を置いた。

{FIXES}
修正は `apply_fixes.py` で機械的に適用（193 置換・27 ファイル・185 hunk、`results/fix_patches.diff`）。修正以外の変更は、`/mnt/data` 依存のパス 4 箇所（前回と同じ）、decompactification の報告行の None 安全化（修正版で指数成長域が無いため）、閉包解析の入力 SOURCE を修正版最終状態へ差し替え、の 3 点のみ（`results/apply_fixes.log`）。

## 1. 全体の結果（旧エンジン → 修正版）

| 量 | 原本（旧エンジン） | 修正版（FIX1-4） | 出典 |
|---|---|---|---|
| 親のノルム / H_total | 1（正規化） | c は N 依存：N5 0.9272（H_total = 0.8596）、N16 1.903（3.62） | summary.json |
| 親残差（残差式は単位ノルム前提） | 1e-13 | N5 0.4869 = σc|1−c²|、N16 73.7（tol 判定は永遠に偽） | 同 |
| seedless の潜伏と onset | 潜伏あり、onset 72〜757 step（N・tol 依存） | **潜伏なし、全て step 1**（N=3,4 は親平面を離れない） | tol_sweep, global_summary |
| 指数成長 | 明瞭（N5 0.1725/step、N16 0.2361/step） | **指数域なし**（fit None / NaN） | decompactification, K_sigma |
| Floquet 最大乗数（N5） | 1.090086569（実 2 重）、三重整合 | **1.001777**（ほぼ中立、Σln|μ| = −3e-4） | floquet_spectrum |
| 厳密保存 | H_total, ZᵀZ 保存（1e-16, 1e-15） | **保存は成立**（H_total = c² が 1e-12、|ZᵀZ| ≈ 1e-13：exp(ΔsK) も実直交） | global_summary |
| 等分配・ヌル錐単体 | |z|² = 1/M（10 桁）、S/lnM = 1.000000 | **等分配しない**（|z|² 相対ばらつき ≈ 2、N5 S/lnM = {ent:.3f}） | cmp5, cmp6 |
| simplex rank | N−1 | N−1（維持） | global_summary |
| 最終クラス数（tol 1e-8） | M（全辺別クラス） | M（N=3,4 は 2） | final_classes |
| 3+3+2+2（N5） | 全 seed で成立 | **成立しない**（4;3;2;1、6;2;1;1 等） | moduli sweep |
| N14 k=6 quasi-closure 9.25e-7 | あり | **消える**（2.4e-5）。1e-6 未満の非自明閉包は N14〜16 に無し | closure search |
| 刻みの連続極限（§12） | rate が g(n) 則で収束 | rate は定義不能だが軌道の RMSE は刻みに一次で収束（時計としては健全） | gamma bundle |
| 長時間挙動（N5） | 秩序化して静止 | 1 本の |z|² が 0.2〜0.43 を周期 ≈ 700 step で振動する準周期運動 | cmp5 |

図：`results/figures/cmp1_Hperp_log_N5_N16.png`（記事図 1 相当：修正版は step 1 で 10⁻² へ跳ぶ）、`cmp2_pump_depletion_N5.png`、`cmp3_tol_sweep_onset.png`、`cmp4_floquet_spectrum.png`、`cmp5_amplitudes_N5.png`、`cmp6_N_series_final_state.png`。

## 2. 記事・論文の主張別の帰結

| 主張 | 修正版での帰結 |
|---|---|
| 記事図 1：10⁻³² から 31 桁のインフレーション | **起きない**。親は修正後の力学の固定点でないため潜伏相が無い |
| 記事図 2：実直交回転による厳密保存 | **成立**（線形回転でも実直交） |
| 記事図 3：onset–残差則 11.616 | **定義不能**（残差 0.4869 固定、onset 1） |
| 記事図 4・5：μ₁ = 1.0901、三重整合 | **不成立**（μ₁ = 1.0018） |
| 記事図 6：完全均等分配 | **不成立**（S/lnM ≈ 0.87、準周期振動） |
| 記事図 7・L154：3+3+2+2、2627/4923 step | 前提（等分配）が崩れる。生成プログラムも未同梱 |
| 記事図 8：モジュライ | 群構造が無いので対象外 |
| 論文 §11 正規化比較・§12 時計 | §11 は比較対象（指数成長）が消失、§12 は時計の一次収束のみ残る |
| 論文 §16 保存定理・§17 ヌル錐定理・§21 120° 定理 | 数学は残る。動力学が到達するという数値主張は消える |
| 論文 §23 N14 の quasi-closure | 消える |

## 3. 修正方針を決めるための論点（木原判断待ち）

1. **親の自己無撞着性**：現在の `make_parent` は位相のみ K の固有モードを（正規化なしで）返す。振幅込み K では K(v)v ≠ iσv なので親は固定点でなく、seedless 実験の設計（固定点＋線形不安定性）が成立しない。理論に忠実にするなら、親も振幅込み K で K(v)v = iσv を満たす固定点として解く必要がある（未実装）。これが最大の分岐点。
2. **親のスケール c**：正規化を外すと親のノルムが N と eig の規約で決まる（N5 0.93、N16 1.90）。振幅込み K は |z|² に比例するので、c は時間単位を変える（step 数の意味が N ごとに違う）。理論としての振幅の単位（何を 1 とするか）を決める必要がある。
3. **時計**：ANGLE = 2π/144 を採用（参照 all3fix は 2π/124）。§12 の RMSE 収束から時計の役割は健全。刻みは 1 定数で変更可。
4. **step 数**：N=5 の準周期振動は周期 ≈ 700 step で 5000 step 内に飽和しない。長時間・アンサンブルの設計が要る。
5. **未同梱プログラム**：N5 の 4 群・時間分離（記事図 7）、N16 位相構造、K/σ 集計、非自明閉包評価は修正版で作り直す必要がある（原本にプログラムが無い）。`analyze_followup.py` のエントロピーは旧 csv を読むため修正版では無効。

## 4. パッケージ一覧（各配下の差異分析 md）

| パッケージ | 目的 | 再生成 / 全 | 未再生成 |
|---|---|---|---|
""" + "\n".join(f"| `{p}` | {d} | {a}/{b} | {c} |" for p,d,a,b,c in INDEX) + f"""

機械突合の集計（336 ファイル）：{', '.join(f'{k} {v}' for k,v in sorted(cnt.items()))}。DIFF の大半は修正の帰結（軌道が別物）。

## ファイル

- `apply_fixes.py`（修正の適用）、`verify_fixes_all.py`（検証）、`run_all_fixed.sh`（`results/jobs.txt`、約 2 分）、`compare_all_fixed.py`、`compare_physics_fixed.py`、`make_comparison_figures.py`、`write_reports.py`、`run_all.sh`（全工程）
- `fixed/<pkg>/`：修正済みプログラム＋再実行データ・図＋差異分析 md
- `results/`：`fix_patches.diff`、`apply_fixes.log`、`verify_fixes_all.{{json,log}}`、`comparison.{{json,md}}`、`physics_comparison.md`、`figures/`、`logs/`
- git には 2 MB 超の CSV を含めない（再生成可、`SHA256SUMS.txt` で検証）
"""
open(os.path.join(H,"全プログラム修正版_差異分析_20260828.md"),"w",encoding="utf-8").write(top); print("reports written:",len(INDEX))
