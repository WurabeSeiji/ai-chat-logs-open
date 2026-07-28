# Stage G source snapshots

Stage GはStage Fの独立統合コピーだけを参照元とする。このディレクトリには、実行開始時点の次のバイト同一スナップショットを保存した。

- `parity_demodulation_STAGE_F_SNAPSHOT.py`
- `state_dependent_scattering_STAGE_F_SNAPSHOT.py`
- `system_A_stage_F_copy_STAGE_F_SNAPSHOT.py`
- `run_stage_F_31_series_check_STAGE_F_SNAPSHOT.py`
- `manifest_STAGE_F_SNAPSHOT.json`

スナップショットは変更・実行しない。新規実装は `code/` に分離する。
