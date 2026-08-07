#!/bin/bash
# 重力読出し論文 Zenodo 公開スクリプト（復旧後に一括実行）
# 手順: 新規deposition → 6ファイルアップロード → メタデータ → publish
set -e
cd "$(dirname "$0")"

echo "=== 1/4 新規deposition作成 ==="
curl -sf --max-time 120 -X POST "https://zenodo.org/api/deposit/depositions" \
  -H "Authorization: Bearer $ZENODO_TOKEN" -H "Content-Type: application/json" \
  -d '{}' > zenodo_new_deposit.json
DEPID=$(python3 -c "import json;print(json.load(open('zenodo_new_deposit.json'))['id'])")
BUCKET=$(python3 -c "import json;print(json.load(open('zenodo_new_deposit.json'))['links']['bucket'])")
echo "deposit id: $DEPID"

echo "=== 2/4 ファイルアップロード（6点）==="
for f in gravity_readout_ja.md gravity_readout_en.md gravity_readout_ja.tex \
         gravity_readout_en.tex gravity_readout_ja.pdf gravity_readout_en.pdf; do
  code=$(curl -s --max-time 300 -o /dev/null -w "%{http_code}" -X PUT "$BUCKET/$f" \
    -H "Authorization: Bearer $ZENODO_TOKEN" \
    -H "Content-Type: application/octet-stream" \
    --data-binary @"$f")
  echo "  $f -> $code"
  [ "$code" = "201" ] || [ "$code" = "200" ] || { echo "UPLOAD FAILED"; exit 1; }
done

echo "=== 3/4 メタデータ設定 ==="
curl -sf --max-time 120 -X PUT "https://zenodo.org/api/deposit/depositions/$DEPID" \
  -H "Authorization: Bearer $ZENODO_TOKEN" -H "Content-Type: application/json" \
  -d @zenodo_metadata_gravity_v1.json > zenodo_metadata_response.json
echo "metadata ok"

echo "=== 4/4 公開 ==="
curl -sf --max-time 120 -X POST "https://zenodo.org/api/deposit/depositions/$DEPID/actions/publish" \
  -H "Authorization: Bearer $ZENODO_TOKEN" > zenodo_publish_response.json
python3 -c "
import json
d = json.load(open('zenodo_publish_response.json'))
print('PUBLISHED')
print('doi:', d.get('doi'))
print('conceptdoi:', d.get('conceptdoi'))
print('url:', d.get('links',{}).get('record_html'))
"
