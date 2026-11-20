# 検算プログラム（2026-11-19）

「QCD を八元数で表現したら（勉強中のメモ）」の検討で使った検算。積の表はメモ (2.1) と同じ。
各プログラムは単独で動く（`python3 check_0X_....py`）。出力は `results/` に同名の .txt で保存した。乱数の種はプログラム内に固定してある。

| プログラム | 確かめたこと | 分析まとめの節 |
|---|---|---|
| check_01_memo_product_eq31_eq32_eq61.py | メモの (3.1)(3.2)(6.1) | §1 |
| check_02_first_term_component_expansion.py | 第1項の成分展開、7本の条件、ゲルマン成分との対応（SymPy 使用） | §3 |
| check_03_rotation_from_octonion_pair.py | 二つの八元数から規則に従う回転 D を作る。14個と8個 | §4 |
| check_04_triangle_loop_C_H_O.py | 交差項で三角形を一周したときのずれ（複素数・四元数・八元数） | §5 |
| check_05_sign_table_128_to_16.py | 符号の付け方128通りのうち八元数は16通り | §6 |
| check_06_remaining_rotations_14_8_3_0.py | 交差項の個数と、区別できずに残る回転の数 | §7 |
| check_07_octonion_pulse_overlap.py | 八元数パルスの重なりの距離依存と軸依存 | §8 |
| check_08_zonal_node_spacing_S1_S2_S3_S7.py | 球面の次元ごとの定在波の節の間隔 | §9 |
