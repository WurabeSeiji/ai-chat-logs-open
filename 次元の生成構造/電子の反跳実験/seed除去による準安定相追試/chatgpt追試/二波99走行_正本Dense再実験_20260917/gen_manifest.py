#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""99-run 事前 manifest 生成（収集仕様書 v1 §2。走行行列の唯一の定義源）。

走行前にコミットして事前登録性を固定する。走行後に本ファイル・manifest を変更しない。
QA 結果は qa_summary.json（別ファイル）に記録する。
"""
import csv
import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))

MAIN_L = [8, 10, 12, 16, 20, 24, 32, 40]
MAIN_PAIRS = [(1, 3), (1, 5), (3, 5),      # odd-odd
              (2, 4), (2, 6), (4, 6),      # even-even
              (1, 2), (2, 3), (3, 6)]      # odd-even
GCD_L = [12, 16, 20, 24, 32, 40]
GCD_PAIRS = [(3, 9), (4, 8), (5, 10)]
CTRL_L = [12, 24, 32]
CTRL_PAIRS = [(1, 3), (2, 4), (2, 3)]
CTRL_DEN = 40
T = 4096


def parity_class(ma, mb):
    pa, pb = ma % 2, mb % 2
    if pa == 1 and pb == 1:
        return 'odd-odd'
    if pa == 0 and pb == 0:
        return 'even-even'
    return 'odd-even'


def build_runs():
    runs = []
    for L in MAIN_L:
        for ma, mb in MAIN_PAIRS:
            runs.append(('main', L, ma, mb, L))
    for L in GCD_L:
        for ma, mb in GCD_PAIRS:
            runs.append(('gcd', L, ma, mb, L))
    for L in CTRL_L:
        for ma, mb in CTRL_PAIRS:
            runs.append(('den40_control', L, ma, mb, CTRL_DEN))
    rows = []
    for series, L, ma, mb, den in runs:
        rows.append(dict(
            run_id=f'L{L}_ma{ma}_mb{mb}_den{den}',
            series=series, L=L, M=L * (L - 1) // 2, ma=ma, mb=mb,
            parity_class=parity_class(ma, mb), den=den, T=T,
            gcd_L_ma=math.gcd(L, ma), gcd_L_mb=math.gcd(L, mb),
            gcd_ma_mb=math.gcd(ma, mb),
            ord_a=L // math.gcd(L, ma), ord_b=L // math.gcd(L, mb),
            ma_mod_L=ma % L, mb_mod_L=mb % L,
            degenerate_check_pass=(ma % L) != (mb % L) and (ma % L) != 0 and (mb % L) != 0))
    ids = [r['run_id'] for r in rows]
    assert len(rows) == 99, len(rows)
    assert len(set(ids)) == 99, '重複 run_id'
    assert all(r['degenerate_check_pass'] for r in rows), 'モード縮退あり'
    assert all(r['L'] % 2 == 0 for r in rows), '奇数 L 混入'
    return rows


def main():
    rows = build_runs()
    path = os.path.join(HERE, 'master_manifest.csv')
    with open(path, 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    from collections import Counter
    print(f'{len(rows)} runs -> master_manifest.csv')
    print('series:', dict(Counter(r['series'] for r in rows)))
    print('parity:', dict(Counter(r['parity_class'] for r in rows)))


if __name__ == '__main__':
    main()
