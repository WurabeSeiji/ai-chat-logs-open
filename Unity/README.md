# Unity 3D / Web (WASM) 開発環境

PC・スマートフォン両対応のインストールレス Web（WebGL / WebAssembly）3D ゲーム・シミュレーション開発用ルートディレクトリです。

- **Unity バージョン**: `6000.3.28f1` (Unity 6.3 LTS, Apple Silicon + Web Build Support)
- **テンプレート**: Universal 3D (`com.unity.template.urp-blank`)

## 構成
- `Assets/Scripts/` : C# スクリプト（`OrbitCameraController.cs`: PC マウス＋スマホタッチ両対応カメラ）
- `Assets/Editor/` : ビルド自動化スクリプト（`WebGLBuilder.cs`）
- `Assets/Scenes/` : シーンファイル（`SampleScene.unity`）
- `Packages/` : URP / Input System などのパッケージ定義
- `ProjectSettings/` : プロジェクト設定
- `build_webgl.sh` : コマンドラインからの WASM ビルドスクリプト

## WASM (WebGL) ビルド手順
```bash
bash Unity/build_webgl.sh
```
出力先：`Unity/Builds/WebGL/`（`.gitignore` により Git 管理外）

※ CLI / バッチモード実行時は `XDG_CACHE_HOME="$HOME/Library/Caches"` が必要です（`build_webgl.sh` に組込済み）。

## ローカル動作確認（WASM ビルド後）
```bash
python3 -m http.server 8080 --directory Unity/Builds/WebGL
```
ブラウザで `http://localhost:8080` を開いて確認します（WASM は `file://` 直開きでは動作しないため HTTP サーバー経由で開きます）。
