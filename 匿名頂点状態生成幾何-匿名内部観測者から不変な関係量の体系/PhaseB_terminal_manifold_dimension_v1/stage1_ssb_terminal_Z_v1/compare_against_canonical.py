#!/usr/bin/env python3
"""第一段の合格判定: 再走行出力を正本SSBバッチ出力と突合する。

合格条件（走行前に事前宣言・2026-09-14）:
  Tier-0（最強）: per_N/N{N}.txt 全35本と ssb_batch_summary.csv が正本とテキスト完全一致。
  Tier-1（数値）: Tier-0 不成立の場合、パース値で以下を全35N×5seedに要求。
    - lock_step: 完全一致
    - H⊥/H:      |diff| <= 5e-4   （正本は %.4f 印字）
    - amp比:     |diff| <= 5e-4
    - Δ:         |diff| <= 2e-6   （正本は %.6f 印字）
    - |Δ+2π/N|:  比が [0.1, 10]   （正本は %.1e 印字）
    - summary CSV: bool 4列完全一致、vac_med/vac_max |diff| <= 0.15 deg
  どちらの Tier で合格したかを必ず報告する。
"""
from pathlib import Path
import json
import re
import sys

HERE = Path(__file__).resolve().parent
CANON = Path('/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/'
             'マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/電子の反跳実験/'
             'seed除去による準安定相追試/chatgpt追試/干渉保存力学_資格審査とシード無し系列_20260831/'
             '位相ロック全N確認_相対平衡_20260912/自発的対称性の破れ_シード依存縮退真空_20260912/'
             '全N6_40_SSBバッチ_20260912')

def parse_per_n(path):
    rows = []
    for line in Path(path).read_text(encoding='utf-8').splitlines():
        m = re.match(r'\s*(\d+)\s+(\d+)\s+([\d.]+)\s+([\d.]+)\s+(-[\d.]+)\s+([\d.eE+-]+)\s*$', line)
        if m:
            rows.append((int(m.group(1)), int(m.group(2)), float(m.group(3)),
                         float(m.group(4)), float(m.group(5)), float(m.group(6))))
    return rows

def main():
    report = {'tier0_text_identical_per_N': 0, 'tier0_summary_identical': False,
              'tier1_failures': [], 'per_N_compared': 0}
    all_text_equal = True
    for N in range(6, 41):
        a = (CANON/'per_N'/f'N{N}.txt').read_text(encoding='utf-8')
        b = (HERE/'per_N'/f'N{N}.txt').read_text(encoding='utf-8')
        if a == b:
            report['tier0_text_identical_per_N'] += 1
        else:
            all_text_equal = False
            ra, rb = parse_per_n(CANON/'per_N'/f'N{N}.txt'), parse_per_n(HERE/'per_N'/f'N{N}.txt')
            if len(ra) != 5 or len(rb) != 5:
                report['tier1_failures'].append((N, 'row count', len(ra), len(rb)))
                continue
            for (s1, l1, h1, ar1, d1, e1), (s2, l2, h2, ar2, d2, e2) in zip(ra, rb):
                if l1 != l2:
                    report['tier1_failures'].append((N, s1, 'lock_step', l1, l2))
                if abs(h1-h2) > 5e-4:
                    report['tier1_failures'].append((N, s1, 'Hperp', h1, h2))
                if abs(ar1-ar2) > 5e-4:
                    report['tier1_failures'].append((N, s1, 'amp_ratio', ar1, ar2))
                if abs(d1-d2) > 2e-6:
                    report['tier1_failures'].append((N, s1, 'delta', d1, d2))
                if not (0.1 <= (e2/e1 if e1 > 0 else 1.0) <= 10.0):
                    report['tier1_failures'].append((N, s1, 'derr', e1, e2))
        report['per_N_compared'] += 1
    ca = (CANON/'ssb_batch_summary.csv').read_text()
    cb = (HERE/'ssb_batch_summary.csv').read_text()
    report['tier0_summary_identical'] = (ca == cb)
    if not report['tier0_summary_identical']:
        for la, lb in zip(ca.splitlines()[1:], cb.splitlines()[1:]):
            fa, fb = la.split(','), lb.split(',')
            if fa[:5] != fb[:5] or fa[7] != fb[7]:
                report['tier1_failures'].append(('summary-bool', la, lb))
            if abs(float(fa[5])-float(fb[5])) > 0.15 or abs(float(fa[6])-float(fb[6])) > 0.15:
                report['tier1_failures'].append(('summary-vac', la, lb))
    tier0 = all_text_equal and report['tier0_summary_identical']
    tier1 = (len(report['tier1_failures']) == 0)
    report['PASS'] = bool(tier0 or tier1)
    report['tier'] = 0 if tier0 else (1 if tier1 else None)
    (HERE/'comparison_report.json').write_text(json.dumps(report, indent=2, ensure_ascii=False, default=str),
                                               encoding='utf-8')
    print(json.dumps(report, indent=2, ensure_ascii=False, default=str))
    sys.exit(0 if report['PASS'] else 1)

if __name__ == '__main__':
    main()
