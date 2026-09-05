#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""対照再現の照合器: 本フォルダで再生成した npz/json を正本フォルダの同名ファイルと
SHA256 で全数照合し、result_paperA_claims_v1.json は凍結マニフェストの SHA とも照合する。
出力: 対照再現_照合結果_v1.json（一致/不一致/未再現/余剰の全リスト）"""
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
CANON = HERE.parent / '対照実験_N掃引1to20_三系_v2'
EXCLUDE = {'論文A_凍結マニフェスト_v1.json', '対照再現_照合結果_v1.json'}

def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for c in iter(lambda: f.read(1 << 20), b''):
            h.update(c)
    return h.hexdigest()

canon_files = {p.name: p for p in CANON.iterdir()
               if p.suffix in ('.npz', '.json') and p.name not in EXCLUDE}
here_files = {p.name: p for p in HERE.iterdir()
              if p.suffix in ('.npz', '.json') and p.name not in EXCLUDE}

match, mismatch, missing = [], [], []
for name, cp in sorted(canon_files.items()):
    if name in here_files:
        (match if sha(here_files[name]) == sha(cp) else mismatch).append(name)
    else:
        missing.append(name)
extra = sorted(set(here_files) - set(canon_files))

manifest = json.load(open(CANON / '論文A_凍結マニフェスト_v1.json'))
claims_sha_manifest = manifest['数値の記録'].get('result_paperA_claims_v1.json')
claims_here = HERE / 'result_paperA_claims_v1.json'
claims_check = (sha(claims_here) == claims_sha_manifest) if claims_here.exists() else None

out = {'n_canonical': len(canon_files), 'n_match': len(match), 'n_mismatch': len(mismatch),
       'n_missing_not_reproduced': len(missing), 'n_extra': len(extra),
       'claims_json_matches_frozen_manifest': claims_check,
       'mismatch': mismatch, 'missing': missing[:400], 'extra': extra[:100]}
with open(HERE / '対照再現_照合結果_v1.json', 'w') as f:
    json.dump(out, f, indent=2, ensure_ascii=False)
print(f"照合: 正本 {len(canon_files)} 件中 一致 {len(match)} / 不一致 {len(mismatch)} / "
      f"未再現 {len(missing)} / 余剰 {len(extra)} | claims凍結SHA一致: {claims_check}")
if mismatch:
    print('不一致（先頭20件）:', mismatch[:20])
