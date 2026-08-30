#!/bin/bash
# 本文 md → PDF（pandoc + lualatex）。Drive 上での lualatex 画像不具合を避けるため /tmp/tex_compile でコンパイルして戻す。
set -e; set -o pipefail
cd "$(dirname "$0")"
NAME=複素シンプレックス基礎_N別全展開_非等モジュラー版_20260830
W=/tmp/tex_compile/$NAME; rm -rf "$W"; mkdir -p "$W"
cp "$NAME.md" "$W/"; cp -r figures "$W/"
( cd "$W" && pandoc "$NAME.md" -f markdown-implicit_figures -o "$NAME.pdf" --pdf-engine=lualatex -V documentclass=ltjsarticle -V mainfont="Hiragino Mincho ProN" -V sansfont="Hiragino Sans" -V monofont="Menlo" -V geometry:margin=20mm -V fontsize=10pt --toc --toc-depth=2 2>&1 | grep -v "Missing character" || true )
cp "$W/$NAME.pdf" ./
echo "PDF: $NAME.pdf"
