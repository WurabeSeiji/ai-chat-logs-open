#!/bin/bash
# Claude Code の memory ディレクトリを Google Drive に移動し、symlink で参照する。
#
# 目的:
#   - 複数マシンの Claude Code 間で memory を共有
#   - Claude.ai（ブラウザ/モバイル）からも同じ memory を参照可能にする
#
# 動作:
#   ~/.claude/projects/<id>/memory/  (実体)
#   → Google Drive 上に移動
#   → 元の場所には symlink を残す
#
# 既に symlink になっているプロジェクトはスキップする（再実行安全）。

set -euo pipefail

# ===== 設定 =====
GDRIVE_TARGET="/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/ClaudeCode/memory"
LOCAL_PROJECTS="$HOME/.claude/projects"

# ===== 関数 =====
print_header() {
    echo
    echo "==============================================="
    echo "  Claude Code memory → Google Drive 同期セットアップ"
    echo "==============================================="
    echo
}

print_step() {
    echo
    echo "▶ [$1] $2"
}

# ===== Pre-flight チェック =====
print_header

print_step "1/6" "Google Drive マウント確認"
if [ ! -d "$(dirname "$(dirname "$(dirname "$GDRIVE_TARGET")")")" ]; then
    echo "  ✗ Google Drive がマウントされていません"
    echo "    期待パス: /Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com"
    echo "    Google Drive for desktop を起動してから再実行してください"
    exit 1
fi
echo "  ✓ Google Drive マウント OK"

print_step "2/6" "ローカル ~/.claude/projects の存在確認"
if [ ! -d "$LOCAL_PROJECTS" ]; then
    echo "  ✗ $LOCAL_PROJECTS が存在しません"
    exit 1
fi
echo "  ✓ $LOCAL_PROJECTS 存在"

print_step "3/6" "Claude Code 関連プロセスの確認"
# Claude Code は node プロセスとして動作することが多い。VSCode 拡張から起動する場合は VSCode が親プロセス。
# ここでは警告のみ。ユーザーは事前に VSCode 終了済みの前提。
claude_procs=$(ps aux | grep -iE "claude|node.*\.claude" | grep -v grep | grep -v "$0" || true)
if [ -n "$claude_procs" ]; then
    echo "  ⚠ Claude Code 関連と思われるプロセスが残っている可能性があります:"
    echo "$claude_procs" | head -5
    echo
    echo "  続行する場合は VSCode と Claude Code セッションを全て終了してから y を入力してください。"
    read -p "  続行しますか? (y/N) " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo "  → 中止しました。"
        exit 1
    fi
else
    echo "  ✓ Claude Code 関連プロセス見当たらず"
fi

print_step "4/6" "Google Drive 側のターゲットディレクトリ作成"
mkdir -p "$GDRIVE_TARGET"
echo "  ✓ $GDRIVE_TARGET"

print_step "5/6" "各プロジェクトの memory/ を移動・symlink 化"
moved=0
skipped_symlink=0
skipped_no_memory=0
errors=0

for proj_dir in "$LOCAL_PROJECTS"/*/; do
    [ -d "$proj_dir" ] || continue
    proj_name=$(basename "$proj_dir")
    src="$proj_dir/memory"
    dst="$GDRIVE_TARGET/$proj_name"

    if [ -L "$src" ]; then
        echo "  ⊝ SKIP (既に symlink): $proj_name"
        skipped_symlink=$((skipped_symlink+1))
    elif [ -d "$src" ]; then
        if [ -e "$dst" ]; then
            echo "  ⚠ ERROR: 移動先に既に何かあります: $dst"
            echo "    手動で確認してください。スキップします。"
            errors=$((errors+1))
            continue
        fi
        echo "  → MOVE: $proj_name"
        mv "$src" "$dst"
        ln -s "$dst" "$src"
        moved=$((moved+1))
    else
        echo "  ⊝ SKIP (memory なし): $proj_name"
        skipped_no_memory=$((skipped_no_memory+1))
    fi
done

echo
echo "  集計: 移動 $moved 件 / 既 symlink $skipped_symlink 件 / memory 無 $skipped_no_memory 件 / エラー $errors 件"

print_step "6/6" "検証"
echo
echo "--- ローカル側 symlink 確認 ---"
ls -la "$LOCAL_PROJECTS"/*/memory 2>/dev/null | head -10 || echo "  (memory なし)"
echo
echo "--- Google Drive 側ディレクトリ ---"
ls "$GDRIVE_TARGET/" 2>/dev/null || echo "  (空)"
echo
echo "--- 現プロジェクトの memory ファイル数 ---"
target_proj="$GDRIVE_TARGET/-Users-kiharahanakira-Library-CloudStorage-GoogleDrive-kihara-noriaki-gmail-com--------OneDrive-GitHub-ai-chat-logs-open"
if [ -d "$target_proj" ]; then
    file_count=$(ls "$target_proj" 2>/dev/null | wc -l | tr -d ' ')
    echo "  $file_count ファイル"
else
    echo "  (移動されていません)"
fi

echo
echo "==============================================="
if [ $errors -eq 0 ]; then
    echo "  ✓ セットアップ完了"
    echo
    echo "  次のステップ:"
    echo "    1. VSCode を再起動"
    echo "    2. このプロジェクトを開いて Claude Code が正常起動するか確認"
    echo "    3. Google Drive for desktop の設定で「ミラーリング」モードを推奨"
else
    echo "  ⚠ エラーが $errors 件あります。上記ログを確認してください。"
fi
echo "==============================================="
