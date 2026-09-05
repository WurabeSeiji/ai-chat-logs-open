Zenodo 収録 zip に関する注記（ChatGPT_denominator_controls_N40_selfcontrol_20260904）
====================================================================================

この zip は軽量版です。再現可能な大容量の状態ファイル

    results*/hm_N40_den_*_states_500.npz など（*_states_500.npz、合計約 231MB）

を Zenodo 収録から除外しています。これらは N=40 一因子実験（段1基準・段2・段3・
段2削除・段2初期化のみ・σ時計）各構成の全 step 複素状態で、プログラムから完全に
再生成できます。

除外していないもの（すべて収録）:
  - 各構成の実行プログラム（run_*.py）・図化プログラム・run_all_*.sh
  - 各 results*/ の summary CSV・timeseries CSV・RUN_METADATA・fcurve CSV・図 PNG
  - README.md・SHA256SUMS.txt（除外した npz の SHA256 も台帳には残してある）

状態 npz の再生成:
  各構成の run_*.py（例: run_N40_staticparent_imK_v1.py）または run_all_*.sh を実行すると
  results*/ に *_states_500.npz が再生成される。生成物が正本と一致することは
  SHA256SUMS.txt（除外ファイル分の SHA256 を含む）で照合できる。
