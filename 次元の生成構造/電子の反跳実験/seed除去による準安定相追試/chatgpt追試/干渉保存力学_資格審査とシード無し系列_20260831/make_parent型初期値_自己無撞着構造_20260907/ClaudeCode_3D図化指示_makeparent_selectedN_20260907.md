# Claude Code 指示書

同じフォルダに置いた基準3D図
`N16_original_note_source_den16_step0_500_3D.html`
の表示仕様をそのまま基準にして、今回の make_parent 全N再走行データを3D図化する。

## 入力データ
`full_N3_N40_sweep/states/`

使う系列は den=N のみ。
対象Nは以下の5つに固定する。

- N=3
- N=4
- N=5
- N=17
- N=40

対応ファイル:

- `full_N3_N40_sweep/states/hm_N3_den_3_states_500.npz`
- `full_N3_N40_sweep/states/hm_N4_den_4_states_500.npz`
- `full_N3_N40_sweep/states/hm_N5_den_5_states_500.npz`
- `full_N3_N40_sweep/states/hm_N17_den_17_states_500.npz`
- `full_N3_N40_sweep/states/hm_N40_den_40_states_500.npz`

## 変更内容
基準3D図の構成・軸・軌跡表示・操作性を維持し、入力データだけ上記5系列へ変更する。
3D図化ロジック自体は勝手に変更しない。
走行系・親生成は再実行しない。
既存NPZを読むだけにする。

## 出力
5つのNそれぞれについて、基準3D図と同形式のHTMLをこのフォルダ配下へ保存する。
ファイル名には N と den を含める。

完了後、使用した図化プログラム名と出力HTML一覧を報告すること。
