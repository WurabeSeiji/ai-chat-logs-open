# System B 低エネルギー標準α近傍 R 感度一覧

本表は、次の 3 つの実験結果を統合し、`R` の昇順に並べたものである。

- `log_sensitivity_alpha_low_standard_theory_v4`
- `log_sensitivity_alpha_low_standard_theory_minus1e-7_v4`
- `log_sensitivity_alpha_low_standard_theory_minus1e-6_v4`

ここでは、各点の絶対的な深さとして `fixed_step_depth_log10` を用いる。
`fixed_step_error` が小さいほど深く、`fixed_step_depth_log10` が大きいほど深い。

## 図

横軸を `R`、縦軸を強度 `fixed_step_depth_log10` として図化した。
上段は全域、下段は中央近傍の拡大である。
破線は中央値 `R = 0.69717787923100305` を示す。
系列色の縦点線は、それぞれの実験で中心として指定した `center_R` を示す。

![System B low-energy standard alpha R sensitivity](system_B_low_alpha_R_sensitivity_depth_v1.svg)

| 実験 | case | R | center_R | delta_R | fixed_step_error | depth_log10 | condition |
|---|---|---:|---:|---:|---:|---:|---|
| R-1e-6中心 | delta_m1em05 | 0.69716687923100307 | 0.69717687923100302 | -1e-05 | 3.773750e-07 | 6.423227 | phi1_s0.01_g0 |
| R-1e-7中心 | delta_m1em05 | 0.69716777923100315 | 0.69717777923100310 | -1e-05 | 3.456967e-07 | 6.461305 | phi1_s0.01_g0 |
| R-1e-6中心 | delta_m1em06 | 0.69717587923100299 | 0.69717687923100302 | -1e-06 | 6.807205e-08 | 7.167031 | phi0_s0.01_g0 |
| R-1e-7中心 | delta_m1em06 | 0.69717677923100307 | 0.69717777923100310 | -1e-06 | 3.805621e-08 | 7.419574 | phi1_s0.01_g0 |
| R-1e-6中心 | delta_m1em07 | 0.69717677923100307 | 0.69717687923100302 | -1e-07 | 3.805621e-08 | 7.419574 | phi1_s0.01_g0 |
| R-1e-6中心 | delta_m1em08 | 0.69717686923100297 | 0.69717687923100302 | -1e-08 | 3.506378e-08 | 7.455141 | phi1_s0.01_g0 |
| R-1e-6中心 | delta_m1em09 | 0.69717687823100305 | 0.69717687923100302 | -1e-09 | 3.476462e-08 | 7.458862 | phi0_s0.01_g0 |
| 厳密値中心 | delta_m1em06 | 0.69717687923100302 | 0.69717787923100305 | -1e-06 | 3.473139e-08 | 7.459278 | phi0_s0.01_g0 |
| R-1e-6中心 | center | 0.69717687923100302 | 0.69717687923100302 | 0 | 3.473139e-08 | 7.459278 | phi0_s0.01_g0 |
| R-1e-6中心 | delta_p1em08 | 0.69717688923100307 | 0.69717687923100302 | 1e-08 | 3.439902e-08 | 7.463454 | phi0_s0.01_g0 |
| R-1e-6中心 | delta_p1em07 | 0.69717697923100297 | 0.69717687923100302 | 1e-07 | 3.140861e-08 | 7.502951 | phi0_s0.01_g0 |
| R-1e-7中心 | delta_m1em07 | 0.69717767923100316 | 0.69717777923100310 | -1e-07 | 8.206714e-09 | 8.085831 | phi1_s0.01_g0 |
| R-1e-7中心 | delta_m1em08 | 0.69717776923100305 | 0.69717777923100310 | -1e-08 | 5.230914e-09 | 8.281422 | phi0_s0.01_g0 |
| R-1e-7中心 | delta_m1em09 | 0.69717777823100313 | 0.69717777923100310 | -1e-09 | 4.933425e-09 | 8.306852 | phi1_s0.01_g0 |
| 厳密値中心 | delta_m1em07 | 0.69717777923100310 | 0.69717787923100305 | -1e-07 | 4.900371e-09 | 8.309771 | phi0_s0.01_g0 |
| R-1e-7中心 | center | 0.69717777923100310 | 0.69717777923100310 | 0 | 4.900371e-09 | 8.309771 | phi0_s0.01_g0 |
| R-1e-7中心 | delta_p1em09 | 0.69717778023100307 | 0.69717777923100310 | 1e-09 | 4.867318e-09 | 8.312710 | phi0_s0.01_g0 |
| R-1e-7中心 | delta_p1em08 | 0.69717778923100315 | 0.69717777923100310 | 1e-08 | 4.569850e-09 | 8.340098 | phi1_s0.01_g0 |
| 厳密値中心 | delta_m1em08 | 0.69717786923100300 | 0.69717787923100305 | -1e-08 | 1.926420e-09 | 8.715249 | phi1_s0.01_g0 |
| 厳密値中心 | delta_m1em09 | 0.69717787823100308 | 0.69717787923100305 | -1e-09 | 1.629116e-09 | 8.788048 | phi1_s0.01_g0 |
| 厳密値中心 | delta_m1em10 | 0.69717787913100304 | 0.69717787923100305 | -1e-10 | 1.599386e-09 | 8.796047 | phi1_s0.01_g0 |
| 厳密値中心 | center | 0.69717787923100305 | 0.69717787923100305 | 0 | 1.596083e-09 | 8.796945 | phi0_s0.01_g0 |
| R-1e-7中心 | delta_p1em07 | 0.69717787923100305 | 0.69717777923100310 | 1e-07 | 1.596083e-09 | 8.796945 | phi0_s0.01_g0 |
| R-1e-6中心 | delta_p1em06 | 0.69717787923100305 | 0.69717687923100302 | 1e-06 | 1.596083e-09 | 8.796945 | phi0_s0.01_g0 |
| 厳密値中心 | delta_p1em10 | 0.69717787933100306 | 0.69717787923100305 | 1e-10 | 1.592780e-09 | 8.797844 | phi1_s0.01_g0 |
| 厳密値中心 | delta_p1em09 | 0.69717788023100302 | 0.69717787923100305 | 1e-09 | 1.563050e-09 | 8.806027 | phi0_s0.01_g0 |
| 厳密値中心 | delta_p1em08 | 0.69717788923100310 | 0.69717787923100305 | 1e-08 | 1.265767e-09 | 8.897646 | phi0_s0.01_g0 |
| 厳密値中心 | delta_p1em07 | 0.69717797923100300 | 0.69717787923100305 | 1e-07 | 1.706669e-09 | 8.767851 | phi0_s0.01_g0 |
| R-1e-7中心 | delta_p1em06 | 0.69717877923100313 | 0.69717777923100310 | 1e-06 | 2.819068e-08 | 7.549894 | phi1_s0.01_g0 |
| 厳密値中心 | delta_p1em06 | 0.69717887923100308 | 0.69717787923100305 | 1e-06 | 3.150939e-08 | 7.501560 | phi1_s0.01_g0 |
| R-1e-6中心 | delta_p1em05 | 0.69718687923100298 | 0.69717687923100302 | 1e-05 | 3.029116e-07 | 6.518684 | phi0_s0.01_g0 |
| R-1e-7中心 | delta_p1em05 | 0.69718777923100306 | 0.69717777923100310 | 1e-05 | 3.341745e-07 | 6.476027 | phi0_s0.01_g0 |
| R-1e-6中心 | delta_p1em04 | 0.69727687923100301 | 0.69717687923100302 | 0.0001 | 3.384430e-06 | 5.470514 | phi1_s0.01_g0 |

## 再現手順

次の 3 つの CSV を入力にする。

```text
波の情報読出し/20260715/minimal_system_B_gray_bugcheck_result_v1/log_sensitivity_alpha_low_standard_theory_v4/minimal_system_B_gray_direct_R_log_sensitivity_summary_v4.csv
波の情報読出し/20260715/minimal_system_B_gray_bugcheck_result_v1/log_sensitivity_alpha_low_standard_theory_minus1e-7_v4/minimal_system_B_gray_direct_R_log_sensitivity_summary_v4.csv
波の情報読出し/20260715/minimal_system_B_gray_bugcheck_result_v1/log_sensitivity_alpha_low_standard_theory_minus1e-6_v4/minimal_system_B_gray_direct_R_log_sensitivity_summary_v4.csv
```

処理手順は以下である。

1. 各 CSV を `csv.DictReader` で読む。
2. 実験ラベルを付与する。
3. `R` を数値として読み、全行を 1 つの配列へ結合する。
4. 結合した配列を `R` の昇順でソートする。
5. `fixed_step_error` と `fixed_step_depth_log10` を絶対深さの指標として表に出す。

再生成用の最小スクリプトは次である。

```python
from pathlib import Path
import csv

sources = [
    (
        "厳密値中心",
        "波の情報読出し/20260715/minimal_system_B_gray_bugcheck_result_v1/"
        "log_sensitivity_alpha_low_standard_theory_v4/"
        "minimal_system_B_gray_direct_R_log_sensitivity_summary_v4.csv",
    ),
    (
        "R-1e-7中心",
        "波の情報読出し/20260715/minimal_system_B_gray_bugcheck_result_v1/"
        "log_sensitivity_alpha_low_standard_theory_minus1e-7_v4/"
        "minimal_system_B_gray_direct_R_log_sensitivity_summary_v4.csv",
    ),
    (
        "R-1e-6中心",
        "波の情報読出し/20260715/minimal_system_B_gray_bugcheck_result_v1/"
        "log_sensitivity_alpha_low_standard_theory_minus1e-6_v4/"
        "minimal_system_B_gray_direct_R_log_sensitivity_summary_v4.csv",
    ),
]

rows = []
for label, path in sources:
    with open(Path(path), newline="") as f:
        for row in csv.DictReader(f):
            rows.append({
                "experiment": label,
                "case": row["case"],
                "R": float(row["R"]),
                "center_R": float(row["center_R"]),
                "delta_R": float(row["delta_R"]),
                "fixed_step_error": float(row["fixed_step_error"]),
                "depth_log10": float(row["fixed_step_depth_log10"]),
                "condition": row["fixed_step_condition"],
            })

rows.sort(key=lambda row: row["R"])

print("| 実験 | case | R | center_R | delta_R | fixed_step_error | depth_log10 | condition |")
print("|---|---|---:|---:|---:|---:|---:|---|")
for row in rows:
    print(
        f"| {row['experiment']} | {row['case']} | "
        f"{row['R']:.17g} | {row['center_R']:.17g} | "
        f"{row['delta_R']:.3g} | {row['fixed_step_error']:.6e} | "
        f"{row['depth_log10']:.6f} | {row['condition']} |"
    )
```
