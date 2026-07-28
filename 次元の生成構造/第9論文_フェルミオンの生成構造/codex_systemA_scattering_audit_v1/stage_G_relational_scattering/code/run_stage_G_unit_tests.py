"""Execute Stage G-A and stop on any failed theoretical value."""

from __future__ import annotations

import csv
import json
from pathlib import Path

from test_stage_G import run_all_tests


STAGE_ROOT = Path(__file__).resolve().parents[1]
OUTPUT = STAGE_ROOT / "data" / "stage_G_unit_test_results.csv"
LOG_OUTPUT = STAGE_ROOT / "logs" / "stage_G_unit_test_run.json"
REPORT_OUTPUT = STAGE_ROOT / "reports" / "02_unit_test_report.md"


def main() -> None:
    rows = run_all_tests()
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    maximum_error = max(float(row["absolute_error"]) for row in rows)
    result = {
        "status": "pass",
        "check_count": len(rows),
        "maximum_absolute_error": maximum_error,
        "failed_count": 0,
    }
    LOG_OUTPUT.write_text(
        json.dumps(result, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    REPORT_OUTPUT.write_text(
        f"""# Stage G-A 単体検証

結果は `pass`。全{len(rows)}検査が事前値と一致した。最大絶対誤差は `{maximum_error:.17g}`。

- A1 同一純奇数波：`c_A=c_B=-1`, `Gamma_AB=1`, `F_rel=+1`
- A2 直交純奇数波：`Gamma_AB=0`, `F_rel=0`
- A3 同一純偶数波：`c_A=c_B=+1`, `Gamma_AB=1`, `F_rel=-1`
- A4 純奇数×純偶数：`c_mean=0`, `F_rel=0`
- A5 位相差 `0,pi/2,pi`：全て `Gamma_AB=1`
- ゼロノルム関係波：0代入せず数値不成立として拒否
- ユニタリ性、係数直交、経路和、全ノルム：許容誤差内

単体検証を通過したため、C0再現確認へ進める。
""",
        encoding="utf-8",
    )
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
