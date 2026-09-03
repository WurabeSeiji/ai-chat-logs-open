# Claude Code 実行指示 — N=3..33 legacy parent denominator control

## 目的
ChatGPTが作成済みの次のプログラムを、**一切編集せずそのまま実行**する。

`run_and_plot_N3_N33_legacyparent_20260903.py`

この作業でClaude Codeが行ってよいことは、入力ファイルの存在確認、実行、終了結果の報告だけである。

## 絶対禁止
- Pythonプログラムを変更しない。
- 新しい実験プログラムを書かない。
- 関数を追加・削除・書換えしない。
- 入力データを変更・再生成しない。
- パスを書換えない。
- dtype、演算順序、`np.linalg.eigh`、BLAS/LAPACK関連設定を変更しない。
- 並列化・高速化・リファクタリングをしない。
- エラー時に推測で修正しない。
- 実験条件を変更しない。
- N=34以降を実行しない。

Python exception、構文エラー、ファイル不足、exit code非0が発生した場合は、その場で停止し、エラーメッセージをそのまま報告すること。
RuntimeWarningはコードを変更せず、そのまま記録して実行を継続すること。

## 今回使用する入力

N=3..16:

`.../干渉保存力学_資格審査とシード無し系列_20260831/data/hm_N*/states_treatment.npz` の `Z[0]`

N=17..33:

`.../干渉保存力学_資格審査とシード無し系列_20260831/data/hm_N*/parent_v.npz` の `v`

**N=17..33では hm_mp_free_N3_N40_20260901 のデータを絶対に使用しない。**

出力先:

`.../干渉保存力学_資格審査とシード無し系列_20260831/ChatGPT_denominator_controls_N3_N33_legacyparent_20260903/results/`

## 実行手順

この指示書と同じフォルダに移動する。

現在の `python3` を確認する。

```bash
which python3
python3 --version
python3 -c "import numpy as np; print(np.__version__)"
```

環境を変更・インストール・更新してはいけない。

次に入力ファイルの存在だけを確認する。

```bash
python3 - <<'PY'
from pathlib import Path
root=Path('/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/電子の反跳実験/seed除去による準安定相追試/chatgpt追試/干渉保存力学_資格審査とシード無し系列_20260831')
missing=[]
for n in range(3,17):
    p=root/'data'/f'hm_N{n}'/'states_treatment.npz'
    if not p.is_file(): missing.append(str(p))
for n in range(17,34):
    p=root/'data'/f'hm_N{n}'/'parent_v.npz'
    if not p.is_file(): missing.append(str(p))
print('missing=',len(missing))
for p in missing: print(p)
PY
```

`missing=0` の場合だけ、次をそのまま実行する。

```bash
python3 run_and_plot_N3_N33_legacyparent_20260903.py
```

## 実行後

`ALL DONE` と exit code 0 を確認する。

次の出力が存在することだけ確認する。

- `results/timeseries_64bit_with124_N3_N33.csv`
- `results/summary_64bit_with124_N3_N33.csv`
- `results/fig_Hperp_denominator_controls_with_124_N3_N33.png`
- `results/RUN_METADATA_N3_N33.json`
- 状態ファイル `hm_N*_den_*_states_500.npz` が 31 N × 6 分母 = **186個**

最後に、実行したPythonとこの指示書のSHA256を表示する。

```bash
shasum -a 256 run_and_plot_N3_N33_legacyparent_20260903.py CLAUDE_CODE_RUN_INSTRUCTION_N3_N33_20260903.md
```

**それ以外の作業は行わない。**
