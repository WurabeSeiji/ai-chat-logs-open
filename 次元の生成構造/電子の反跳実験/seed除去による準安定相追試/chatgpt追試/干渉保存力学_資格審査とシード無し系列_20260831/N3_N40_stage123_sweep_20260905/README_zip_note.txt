Zenodo 収録 zip に関する注記（N3_N40_stage123_sweep_20260905）
================================================================

この zip は軽量版です。再現可能な大容量の状態ファイル

    results/hm_N*_den_*_states_500.npz  （228 個、合計約 460MB）

を Zenodo 収録から除外しています。これらは各走行の全 501 step の複素状態で、
プログラムから完全に再生成できます。

除外していないもの（すべて収録）:
  - プログラム 6 本（make_static_parents / run_N3_N40_stage123 / check_sweep_inputs /
    analyze_sweep_summary / plot_complex_plane / run_n_scaling_lowrank）と run_all.sh
  - 初期データ parents/ の静的親 npz 38 個＋台帳（約 632KB。これは入力なので保持）
  - results/ の timeseries CSV・summary CSV・集計 JSON・RUN_METADATA・図
  - 論文（paper_overview / paper_ch1 / paper_ch2 / paper_ch3 の md/tex/pdf 日英）
  - README.md・SHA256SUMS.txt（除外した npz の SHA256 も台帳には残してある）

状態 npz の再生成:
    ./run_all.sh
  もしくは
    python3 run_N3_N40_stage123_v1.py
  で results/ に 228 個の states_500.npz が再生成される。生成物が正本と一致することは
  SHA256SUMS.txt（除外ファイル分の SHA256 を含む）で照合できる。
