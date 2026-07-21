"""一般位置ランク等号定理の厳密有理数証人検証 v1

tan半角パラメータ t_e = e (e=1..M) から単位円上の有理点
(c_e, s_e) = ((1-t^2)/(1+t^2), 2t/(1+t^2)) を取り、
K_ef = A_ef (s_f c_e - c_f s_e) を厳密有理数で構成して
分数演算のガウス消去でランクを厳密に計算する。

rank K = 2*min(N, floor(M/2)) が各Nで厳密に成立すれば、
そのNの証人配置が確定し、実解析性論法（Mityagin）により
一般位置での等号が測度零の例外を除いて成立する。
"""

from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent
RESULT_DIR = BASE_DIR / "nbody_rank_saturation_preliminary_result_v1"


def relation_pairs(body_count: int):
    return [(i, j) for i in range(body_count) for j in range(i + 1, body_count)]


def exact_rank(matrix: list) -> int:
    """分数演算ガウス消去による厳密ランク。"""
    rows = [row[:] for row in matrix]
    row_count = len(rows)
    col_count = len(rows[0]) if rows else 0
    rank = 0
    pivot_row = 0
    for col in range(col_count):
        pivot = None
        for r in range(pivot_row, row_count):
            if rows[r][col] != 0:
                pivot = r
                break
        if pivot is None:
            continue
        rows[pivot_row], rows[pivot] = rows[pivot], rows[pivot_row]
        pv = rows[pivot_row][col]
        for r in range(pivot_row + 1, row_count):
            if rows[r][col] != 0:
                factor = rows[r][col] / pv
                rows[r] = [a - factor * b for a, b in zip(rows[r], rows[pivot_row])]
        pivot_row += 1
        rank += 1
        if pivot_row == row_count:
            break
    return rank


def witness_rank(body_count: int) -> dict:
    pairs = relation_pairs(body_count)
    relation_count = len(pairs)
    cos_values = []
    sin_values = []
    for index in range(relation_count):
        t = Fraction(index + 1)
        denom = 1 + t * t
        cos_values.append((1 - t * t) / denom)
        sin_values.append(2 * t / denom)
    matrix = []
    for a in range(relation_count):
        row = []
        for b in range(relation_count):
            adjacent = a != b and bool(set(pairs[a]) & set(pairs[b]))
            if adjacent:
                value = sin_values[b] * cos_values[a] - cos_values[b] * sin_values[a]
            else:
                value = Fraction(0)
            row.append(value)
        matrix.append(row)
    rank = exact_rank(matrix)
    expected = 2 * min(body_count, relation_count // 2)
    return {
        "body_count": body_count,
        "relation_count": relation_count,
        "witness_parameters": f"t_e = e, e=1..{relation_count}",
        "exact_rank": rank,
        "expected_rank": expected,
        "witness_confirms_equality": bool(rank == expected),
    }


def main() -> None:
    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    rows = [witness_rank(n) for n in range(3, 13)]
    for row in rows:
        print(row)
    payload = {
        "description": "exact rational witness verification of generic rank equality",
        "rows": rows,
        "all_confirmed": all(row["witness_confirms_equality"] for row in rows),
    }
    with (RESULT_DIR / "nbody_exact_rational_witness_v1.json").open(
        "w", encoding="utf-8"
    ) as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)
    print("all_confirmed:", payload["all_confirmed"])


if __name__ == "__main__":
    main()
