# Claude Code 実行指示 — N=3..40 denominator control

## 目的
ChatGPTが作成済みの次のプログラムを、**一切編集せずそのまま実行**する。

`run_and_plot_N3_N40_mixedseed_20260903.py`

この作業でClaude Codeが行ってよいことは、環境確認と上記プログラムの実行だけである。

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

エラーが発生した場合は、その場で停止し、エラーメッセージをそのまま報告すること。

## プログラムが使用する入力
N=3..16:

`.../干渉保存力学_資格審査とシード無し系列_20260831/data/hm_N*/states_treatment.npz` の `Z[0]`

N=17..40:

`.../干渉保存力学_資格審査とシード無し系列_20260831/hm_mp_free_N3_N40_20260901/data/hm_N*/parent_v.npz` の `v`

出力先:

`.../干渉保存力学_資格審査とシード無し系列_20260831/ChatGPT_denominator_controls_N3_N40_mixedseed_20260903/results/`

## 実行手順
まず、この指示書と同じフォルダに移動する。

次に、現在の実行環境を**変更せず**確認する。

```bash
{
  pwd
  which python3
  python3 -c "import sys,platform,numpy as np; print('python_executable=',sys.executable); print('python=',platform.python_version()); print('numpy=',np.__version__); print('platform=',platform.platform()); np.show_config()"
} | tee EXECUTION_ENVIRONMENT_N3_N40_20260903.txt
```

環境を変更する必要はない。この確認で使用した同じ `python3` を、そのまま実験実行にも使用する。

その後、入力ファイルの存在だけを確認する。内容を変更してはいけない。

```bash
python3 - <<'PY'
from pathlib import Path
root=Path('/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/電子の反跳実験/seed除去による準安定相追試/chatgpt追試/干渉保存力学_資格審査とシード無し系列_20260831')
missing=[]
for n in range(3,17):
    p=root/'data'/f'hm_N{n}'/'states_treatment.npz'
    if not p.is_file(): missing.append(str(p))
for n in range(17,41):
    p=root/'hm_mp_free_N3_N40_20260901'/'data'/f'hm_N{n}'/'parent_v.npz'
    if not p.is_file(): missing.append(str(p))
print('missing=',len(missing))
for p in missing: print(p)
PY
```

`missing=0` の場合だけ実行する。

```bash
python3 run_and_plot_N3_N40_mixedseed_20260903.py
```

## 実行後
プログラムの `ALL DONE` を確認する。

次のファイルが存在することだけ確認する。

- `results/timeseries_64bit_with124_N3_N40.csv`
- `results/summary_64bit_with124_N3_N40.csv`
- `results/fig_Hperp_denominator_controls_with_124_N3_N40.png`
- `results/RUN_METADATA_N3_N40.json`
- `results/hm_N3_den_1_states_500.npz` から各N・各分母の状態ファイル

最後に、実行したPythonのSHA256とこの指示書のSHA256を表示する。

```bash
shasum -a 256 run_and_plot_N3_N40_mixedseed_20260903.py CLAUDE_CODE_RUN_INSTRUCTION_N3_N40_20260903.md
```

**それ以外の作業は行わない。**
