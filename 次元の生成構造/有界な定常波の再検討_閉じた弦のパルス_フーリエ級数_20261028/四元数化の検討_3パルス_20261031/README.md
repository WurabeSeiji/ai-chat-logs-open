# 四元数化の検討（3 パルス）：index_v2c3.html の四元数化と複素数との対応表（2026-10-31）

`../複素化の検討_2パルス_20261029/index_v2c3.html`（同じ弦の上の三つの複素パルス、観測者の時空）を四元数にした。核心は q ＝ a ＋ n b ＝ A e^{nθ} で、複素数と同じ形。各パルスは平面 span(1, n_i) ≅ ℂ の中にあるので式は i → n_i の置き換えだけで、変わるのは対の交差項（スカラー部に n_i·n_j、ベクトル部に n_i×n_j）だけ。原本を `原本/` に無変更でコピーし、n を既定の k にすると原本の数値に戻ることを Python と HTML の両方で機械的に確かめた。

| ファイル | 内容 |
|---|---|
| `分析_四元数化した閉じた弦の3パルス_複素数との対応表_20261031.md` | 分析 v0.1。要点 → 約束 → 対応表（複素と四元数を横に並べる）→ 導出 §3.1〜3.7（純虚四元数の積、Sc と Vec、親のノルム、対の位置ベクトルと向きの歳差、右パネルの持ち上げ、鏡映の零閉塞と平面の共有、実のパルスに軸が見えない）→ 対照 → 軸を変えた確認値 → 実験系の仕様 → 描画記録 → 区分 → 再現 |
| `index_v2q3.html` | 実験系。index_v2c3.html の忠実コピーに、各パルスの虚部の向き n_ℓ のスライダー (i, j, k) ∈ [−1, 1]³（正規化、既定 k。i ＝ 円周・j ＝ 半径・k ＝ 垂直の名付けは読み出し写像の選択）、和を 4 成分 W・V_i・V_j・V_k で計算し 4 色で描く（黒 ＝ W 半径方向、赤 ＝ V_j 半径方向、青 ＝ V_k 垂直、緑 ＝ V_i を垂直に反対符号 −V_i。円周方向に乗せると位置と同じ向きで見えないため）、交差項を四元数の掛け算 `qmul` でスカラー部とベクトル部に分け、読める量に n_i·n_j・n_i×n_j・対の位置ベクトル vec Δ_ij（e^{n_iθ_i} e^{−n_jθ_j} の角 × 軸）・交差項のベクトル部 vec V_ij/(R_iR_j)・Σc² のベクトル部を追加。巻き数は軸が全部同じときだけ。右パネルの横軸は持ち上げ x⃗_i ＝ Θ_i n_i の相対 ±\|x⃗_io\|（n が同じなら元の Δ_io）。URL 引数に `ni1= nj1= nk1=`（2・3 も同様） |
| `原本/index_v2c3.html`、`原本/check_values_v2c3.py` | 原本の無変更コピー（SHA は原本フォルダの SHA256SUMS と同じ：491eeed8…、8c4b69fb…） |
| `check_values_v2q3.py` → `results/results_check_values_v2q3.txt` | 確認計算。原本の 9 例（軸 k）＋ 軸を変えた 10 例（同じ位置で i, j, k ／ k, k, −k ／ 45°、三角形で i, j, k、鏡映配置で k, i, i ／ k, i, −i ／ i, j, k）＋ 速度ベクトル ＋ 持ち上げ。閉じた式と掛け算の一致、R² ＝ Σ\|c_m\|² を assert。末尾で `原本/check_values_v2c3.py` を import し全項目を突き合わせる（最大差 2.1e-14、一致） |
| `shoot.sh` | ヘッドレス Chrome で描画記録 8 枚（`figs/`）と、対照用の DOM 8 本（既定・三角形・鏡映 ±45°・標準の対 × 原本／四元数版、`results/dom_*.html`） |
| `control_html.py` → `results/results_control_html.txt` | DOM の右上の読める量を項目ごとに比べる。共通の項目 32・31・33・31 個がすべて一致 |
| `results/diff_index_v2c3_v2q3.txt`、`results/diff_check_values_v2c3_v2q3.txt` | 原本との差分（変更行の全部） |
| `figs/v2q3_*.png` | 描画記録：D0（原本 v2c3_D0 と同じ配置）、tri_kkk、tri_ijk、same_ijk、mirror_k_i_i、mirror_k_i_mi、cone、cone_ijk |
| `run_all.sh` | 確認計算・差分・DOM の集計・SHA256SUMS（`--shoot` で描画記録と DOM の取り直しも） |
| `SHA256SUMS` | 一式の SHA-256 |

## 要点

    q = a + n b = A e^{nθ} = A(cos θ + n sin θ)、n² = −1。平面 span(1, n) ≅ ℂ
    n_i n_j = −n_i·n_j + n_i×n_j
    Sc[e^{n_i a} e^{−n_j b}] = cos a cos b + (n_i·n_j) sin a sin b     （n_i = n_j で cos(a − b)）
    Vec[e^{n_i a} e^{−n_j b}] = sin a cos b n_i − cos a sin b n_j − sin a sin b (n_i×n_j)
    R² = Σ_i R_i² + 2 Σ_{i<j} Σ_m c⁽ⁱ⁾c⁽ʲ⁾ Sc[e^{n_i α_i} e^{−n_j α_j}]、α_i = φ₀i − mθ_i。R に効くのはスカラー部だけ
    対の位置 g_ij = e^{n_iθ_i} e^{−n_jθ_j}：角は倍音の刻みで一定、軸は n_j のまわりを歳差 → 対では距離は決まり向きは決まらない
    同じ位置 θ = 150 で軸 i, j, k：cos Θ = 0.4311（cos·cos だけ）、|vec Δ| = 41.4°、R = 0.3958（軸 k, k, k なら 3R₁ = 0.5024）
    鏡映 φ₀ = 0, ±45：対が平面を共有（k,k,k／k,i,i）で Σc²/R² = 1、逆軸 k,i,−i で 0.4074 + (−0.2905, 0, 0)
    θ = 180・φ₀ = 0 のパルスは係数が実で軸が見えない（三角形 i,j,k で cos Θ₁₂・cos Θ₁₃ 不変、cos Θ₂₃ だけ変わる）
    対照：原本と 9 例一致（2.1e-14）、HTML の読める量 4 URL 一致、図は曲線の色だけ違う

## 再現

    sh run_all.sh            # 確認計算・差分・DOM の集計・SHA
    sh run_all.sh --shoot    # 描画記録と DOM の取り直しも（Chrome と通信が要る）

環境：macOS、Python 3.9、Google Chrome（ヘッドレス、`--headless=new --virtual-time-budget=6000`）、Plotly 2.34.0（cdnjs）。撮影ファイル名は大文字小文字を使い分けない（macOS は区別しない）。
