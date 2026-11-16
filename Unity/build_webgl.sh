#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
UNITY_BIN="/Applications/Unity/Hub/Editor/6000.3.28f1/Unity.app/Contents/MacOS/Unity"

export XDG_CACHE_HOME="${XDG_CACHE_HOME:-$HOME/Library/Caches}"
mkdir -p "$XDG_CACHE_HOME/bee"

echo "[1/2] Building WebGL (WASM)..."
"$UNITY_BIN" \
  -batchmode \
  -quit \
  -projectPath "$SCRIPT_DIR" \
  -buildTarget WebGL \
  -executeMethod WebGLBuilder.BuildWebGL \
  -logFile "$SCRIPT_DIR/Logs/webgl_build.log"

echo "[2/2] Build complete: $SCRIPT_DIR/Builds/WebGL"
ls -lh "$SCRIPT_DIR/Builds/WebGL"
