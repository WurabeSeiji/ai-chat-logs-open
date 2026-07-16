# System B 高エネルギー標準α近傍 R 感度一覧

本表は、次の 3 つの高エネルギー側実験結果を統合し、`R` の昇順に並べたものである。

- `log_sensitivity_alpha_MZ_standard_theory_v4`
- `log_sensitivity_alpha_MZ_standard_theory_plus1e-6_v4`
- `log_sensitivity_alpha_MZ_standard_theory_plus1e-5_v4`

ここでは、各点の絶対的な深さとして `fixed_step_depth_log10` を用いる。
`fixed_step_error` が小さいほど深く、`fixed_step_depth_log10` が大きいほど深い。

## 図

横軸を `R`、縦軸を強度 `fixed_step_depth_log10` として図化した。
上段は全域、下段は中央近傍の拡大である。
破線は高エネルギー側の中央 R `R = 0.68660290255614798` を示す。
系列色の縦点線は、それぞれの実験で中心として指定した `center_R` を示す。

![System B high-energy standard alpha R sensitivity](system_B_high_alpha_R_sensitivity_depth_v1.svg)

| 実験 | case | R | center_R | delta_R | fixed_step_error | depth_log10 | condition |
|---|---|---:|---:|---:|---:|---:|---|
| MZ+1e-5中心 | delta_m1em04 | 0.68651290255614794 | 0.68661290255614793 | -0.0001 | 4.006228e-06 | 5.397264 | phi0_s0.01_g0 |
| MZ標準値中心 | delta_m1em05 | 0.68659290255614802 | 0.68660290255614798 | -1e-05 | 3.873563e-07 | 6.411889 | phi0_s0.01_g0 |
| MZ+1e-6中心 | delta_m1em05 | 0.68659390255614805 | 0.686603902556148 | -1e-05 | 3.440362e-07 | 6.463396 | phi0_s0.01_g0 |
| MZ標準値中心 | delta_m1em06 | 0.68660190255614795 | 0.68660290255614798 | -1e-06 | 3.833587e-08 | 7.416395 | phi1_s0.01_g0 |
| MZ標準値中心 | delta_m1em07 | 0.68660280255614803 | 0.68660290255614798 | -1e-07 | 1.008782e-08 | 7.996203 | phi0_s0.01_g0 |
| MZ標準値中心 | delta_m1em08 | 0.68660289255614793 | 0.68660290255614798 | -1e-08 | 1.123455e-08 | 7.949444 | phi1_s0.01_g0 |
| MZ標準値中心 | delta_m1em09 | 0.686602901556148 | 0.68660290255614798 | -1e-09 | 1.143494e-08 | 7.941766 | phi1_s0.01_g0 |
| MZ標準値中心 | center | 0.68660290255614798 | 0.68660290255614798 | 0 | 1.145817e-08 | 7.940885 | phi0_s0.01_g0 |
| MZ+1e-6中心 | delta_m1em06 | 0.68660290255614798 | 0.686603902556148 | -1e-06 | 5.007777e-08 | 7.300355 | phi0_s0.01_g0 |
| MZ+1e-5中心 | delta_m1em05 | 0.68660290255614798 | 0.68661290255614793 | -1e-05 | 3.924560e-07 | 6.406209 | phi0_s0.01_g0 |
| MZ標準値中心 | delta_p1em09 | 0.68660290355614795 | 0.68660290255614798 | 1e-09 | 1.148159e-08 | 7.939998 | phi1_s0.01_g0 |
| MZ標準値中心 | delta_p1em08 | 0.68660291255614803 | 0.68660290255614798 | 1e-08 | 1.170103e-08 | 7.931776 | phi1_s0.01_g0 |
| MZ標準値中心 | delta_p1em07 | 0.68660300255614792 | 0.68660290255614798 | 1e-07 | 1.475262e-08 | 7.831131 | phi0_s0.01_g0 |
| MZ+1e-6中心 | delta_m1em07 | 0.68660380255614806 | 0.686603902556148 | -1e-07 | 1.760582e-08 | 7.754344 | phi0_s0.01_g0 |
| MZ+1e-6中心 | delta_m1em08 | 0.68660389255614795 | 0.686603902556148 | -1e-08 | 1.805130e-08 | 7.743491 | phi1_s0.01_g0 |
| MZ+1e-6中心 | delta_m1em09 | 0.68660390155614803 | 0.686603902556148 | -1e-09 | 1.814719e-08 | 7.741191 | phi1_s0.01_g0 |
| MZ標準値中心 | delta_p1em06 | 0.686603902556148 | 0.68660290255614798 | 1e-06 | 2.908319e-08 | 7.536358 | phi0_s0.01_g0 |
| MZ+1e-6中心 | center | 0.686603902556148 | 0.686603902556148 | 0 | 1.815842e-08 | 7.740922 | phi0_s0.01_g0 |
| MZ+1e-6中心 | delta_p1em09 | 0.68660390355614798 | 0.686603902556148 | 1e-09 | 1.816977e-08 | 7.740651 | phi1_s0.01_g0 |
| MZ+1e-6中心 | delta_p1em08 | 0.68660391255614805 | 0.686603902556148 | 1e-08 | 1.827706e-08 | 7.738094 | phi0_s0.01_g0 |
| MZ+1e-6中心 | delta_p1em07 | 0.68660400255614795 | 0.686603902556148 | 1e-07 | 1.986340e-08 | 7.701946 | phi0_s0.01_g0 |
| MZ+1e-6中心 | delta_p1em06 | 0.68660490255614803 | 0.686603902556148 | 1e-06 | 2.430012e-08 | 7.614392 | phi0_s0.01_g0 |
| MZ+1e-5中心 | delta_m1em06 | 0.6866119025561479 | 0.68661290255614793 | -1e-06 | 8.691152e-08 | 7.060923 | phi0_s0.01_g0 |
| MZ+1e-5中心 | delta_m1em07 | 0.68661280255614798 | 0.68661290255614793 | -1e-07 | 5.495632e-08 | 7.259982 | phi0_s0.01_g0 |
| MZ+1e-5中心 | delta_m1em08 | 0.68661289255614788 | 0.68661290255614793 | -1e-08 | 5.255798e-08 | 7.279361 | phi1_s0.01_g0 |
| MZ標準値中心 | delta_p1em05 | 0.68661290255614793 | 0.68660290255614798 | 1e-05 | 3.171602e-07 | 6.498721 | phi1_s0.01_g0 |
| MZ+1e-5中心 | center | 0.68661290255614793 | 0.68661290255614793 | 0 | 5.195306e-08 | 7.284389 | phi1_s0.01_g0 |
| MZ+1e-5中心 | delta_p1em08 | 0.68661291255614798 | 0.68661290255614793 | 1e-08 | 5.137371e-08 | 7.289259 | phi0_s0.01_g0 |
| MZ+1e-5中心 | delta_p1em07 | 0.68661300255614788 | 0.68661290255614793 | 1e-07 | 4.731063e-08 | 7.325041 | phi0_s0.01_g0 |
| MZ+1e-6中心 | delta_p1em05 | 0.68661390255614796 | 0.686603902556148 | 1e-05 | 4.114335e-07 | 6.385700 | phi1_s0.01_g0 |
| MZ+1e-5中心 | delta_p1em06 | 0.68661390255614796 | 0.68661290255614793 | 1e-06 | 1.206355e-07 | 6.918525 | phi1_s0.01_g0 |
| MZ+1e-5中心 | delta_p1em05 | 0.68662290255614788 | 0.68661290255614793 | 1e-05 | 3.281627e-07 | 6.483911 | phi0_s0.01_g0 |
| MZ+1e-5中心 | delta_p1em04 | 0.68671290255614792 | 0.68661290255614793 | 0.0001 | 3.232288e-06 | 5.490490 | phi1_s0.01_g0 |

## 再現手順

次の 3 つの CSV を入力にする。

```text
波の情報読出し/20260715/minimal_system_B_gray_bugcheck_result_v1/log_sensitivity_alpha_MZ_standard_theory_v4/minimal_system_B_gray_direct_R_log_sensitivity_summary_v4.csv
波の情報読出し/20260715/minimal_system_B_gray_bugcheck_result_v1/log_sensitivity_alpha_MZ_standard_theory_plus1e-6_v4/minimal_system_B_gray_direct_R_log_sensitivity_summary_v4.csv
波の情報読出し/20260715/minimal_system_B_gray_bugcheck_result_v1/log_sensitivity_alpha_MZ_standard_theory_plus1e-5_v4/minimal_system_B_gray_direct_R_log_sensitivity_summary_v4.csv
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
    ('MZ標準値中心', '波の情報読出し/20260715/minimal_system_B_gray_bugcheck_result_v1/log_sensitivity_alpha_MZ_standard_theory_v4/minimal_system_B_gray_direct_R_log_sensitivity_summary_v4.csv'),
    ('MZ+1e-6中心', '波の情報読出し/20260715/minimal_system_B_gray_bugcheck_result_v1/log_sensitivity_alpha_MZ_standard_theory_plus1e-6_v4/minimal_system_B_gray_direct_R_log_sensitivity_summary_v4.csv'),
    ('MZ+1e-5中心', '波の情報読出し/20260715/minimal_system_B_gray_bugcheck_result_v1/log_sensitivity_alpha_MZ_standard_theory_plus1e-5_v4/minimal_system_B_gray_direct_R_log_sensitivity_summary_v4.csv'),
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
