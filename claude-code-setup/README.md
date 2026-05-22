# Claude Code memory → Google Drive 同期セットアップ

このフォルダは、Claude Code の auto memory ディレクトリを Google Drive 上に移動し、symlink で参照するためのセットアップを提供する。

## 目的

3 経路からの memory アクセスを実現する：

```
GoogleDrive/.../OneDrive/ClaudeCode/memory/<project-id>/
  ├─ Claude Code（このマシン）         自動ロード + 書き込み（symlink 経由）
  ├─ Claude Code（別マシン）           別 namespace で自動ロード + 書き込み（読み取り共有可）
  └─ Claude.ai（ブラウザ/モバイル）    Google Drive コネクタ経由で読み取り
```

## 設計の特徴

- **書き込み競合ゼロ**: プロジェクト識別子はマシンごとの絶対パスから生成されるため、別マシンは別フォルダに書く。同時利用しても物理的に分離される
- **軽量同期**: セッション履歴 jsonl（数百MB）はローカルに残し、memory（数MB）だけクラウド同期
- **再実行安全**: 既に symlink になっているプロジェクトはスキップ。複数回実行しても破損しない
- **プライバシー境界**: 個人 Google Drive アカウント内に閉じる

## 配置

```
GoogleDrive/マイドライブ/OneDrive/
├── GitHub/                                          ← 既存（git リポジトリ群）
│   └── ai-chat-logs-open/                          ← 現プロジェクト
│       └── claude-code-setup/                       ← このフォルダ
│           ├── setup-memory-sync.sh
│           └── README.md
└── ClaudeCode/                                      ← 新規作成（スクリプト実行後）
    └── memory/                                      ← 各プロジェクトの memory が集約される
        ├── -Users-kiharahanakira-...-ai-chat-logs-open/
        ├── -Users-kiharahanakira-...-work-WF.../
        └── -Users-kiharahanakira-...-flutter-application-2/
```

## 前提条件

1. **VSCode と全 Claude Code セッションが終了している**
2. **Google Drive for desktop が起動しており、`/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com` がマウントされている**
3. **Google Drive クライアントが「ミラーリング」モードになっている**（Stream モードだと頻繁な書き込みでエラーが起きやすい）

## 実行手順

```bash
# 1. VSCode 終了

# 2. ターミナルで実行
cd /Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/claude-code-setup
./setup-memory-sync.sh

# 3. ログを確認、エラーがなければ VSCode を再起動
```

スクリプトは pre-flight チェックとして以下を行う：
- Google Drive マウント確認
- `~/.claude/projects` 存在確認
- Claude Code プロセスの残存確認（残っていれば警告 + 続行確認）

その後：
- `OneDrive/ClaudeCode/memory/` を作成
- 各プロジェクトの `memory/` を移動して symlink に置換
- 結果を表示

## 期待される出力

```
▶ [1/6] Google Drive マウント確認
  ✓ Google Drive マウント OK
▶ [2/6] ローカル ~/.claude/projects の存在確認
  ✓ ~/.claude/projects 存在
▶ [3/6] Claude Code 関連プロセスの確認
  ✓ Claude Code 関連プロセス見当たらず
▶ [4/6] Google Drive 側のターゲットディレクトリ作成
  ✓ /Users/.../OneDrive/ClaudeCode/memory
▶ [5/6] 各プロジェクトの memory/ を移動・symlink 化
  → MOVE: -Users-kiharahanakira-...-ai-chat-logs-open
  → MOVE: -Users-kiharahanakira-...-work-WF...
  → MOVE: -Users-kiharahanakira-...-flutter-application-2
  集計: 移動 3 件 / 既 symlink 0 件 / memory 無 0 件 / エラー 0 件
▶ [6/6] 検証
  (検証結果)
  ✓ セットアップ完了
```

## 動作確認

VSCode 再起動後、Claude Code に以下のように尋ねて確認する：

```
「メモリ確認、宮脇先生について何書いてある？」
```

正常に動作していれば、symlink 先（Google Drive 上）の memory が透過的に読み込まれる。

## Claude.ai 側からの参照

設定後、ブラウザ Claude.ai から以下のように指示できる：

```
「Google Drive の MEMORY.md を見て、シグマサロンの予定を教えて」
「Drive で user_contacts_morita を検索して内容を要約」
```

Claude.ai の Google Drive コネクタが有効化されていれば、同じファイル群を参照可能。

## ロールバック

万一問題があれば、symlink を解除して元に戻す：

```bash
# 各プロジェクトの symlink を実体に戻す
for proj_dir in ~/.claude/projects/*/; do
    proj_name=$(basename "$proj_dir")
    src="$proj_dir/memory"
    if [ -L "$src" ]; then
        target=$(readlink "$src")
        rm "$src"
        mv "$target" "$src"
        echo "Restored: $proj_name"
    fi
done

# Google Drive 側の空フォルダを削除
rmdir "/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/ClaudeCode/memory" 2>/dev/null
rmdir "/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/ClaudeCode" 2>/dev/null
```

## 別マシンでの利用

別マシン側でも本スクリプトを 1 回実行すれば、同じ Google Drive 上の `OneDrive/ClaudeCode/memory/` 配下に **そのマシン固有のプロジェクト識別子**で memory が集約される。

例：
```
OneDrive/ClaudeCode/memory/
├── -Users-kiharahanakira-...-ai-chat-logs-open/      ← Mac A 専用
├── -Users-otheruser-...-ai-chat-logs-open/           ← Mac B 専用
└── ...
```

別マシンの memory を参照したい時は、Claude Code に「Google Drive の `ClaudeCode/memory/-Users-otheruser-...-ai-chat-logs-open/MEMORY.md` を読んで」と指示するか、Claude.ai で同様に指定する。

## 注意事項

- **同一プロジェクトを 2 マシン同時に編集しない**：プロジェクト識別子が違えば自然分離されるが、同一識別子（同一ユーザー・同一作業パス）の場合は競合し得る
- **conflicted copy の監視**：稀に `*[Conflict]*.md` ファイルが生まれる可能性。週次で確認推奨
- **Google アカウント保護**：2 要素認証必須。memory には機微情報を含むため
