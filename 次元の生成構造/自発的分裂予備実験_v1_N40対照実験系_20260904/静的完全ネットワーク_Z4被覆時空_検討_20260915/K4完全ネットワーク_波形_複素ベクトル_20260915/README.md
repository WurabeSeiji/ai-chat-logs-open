# U^4=I 基底波＋倍音1本：完全ネットワークと波形・複素ベクトル図

一回の実行で m=2（偶数倍音）と m=3（奇数倍音）の全図を生成する。

## 定義

U = exp(2πi/4) = i

S_k = (U^k, U^(m k)),  k=0,1,2,3

波形表示では振幅を1に固定し、

base = cos(theta)
harmonic = cos(m theta)
composite = cos(theta) + cos(m theta)

を描く。

## 複素ベクトル図

各離散状態 k=0,1,2,3 について
z_1 = U^k
z_m = U^(mk)
の2本を複素平面上のベクトルとして描く。

さらに関係自体を見やすくする代表図として theta=pi/4 における
z=e^(i theta), z_m=e^(i m theta)
の2ベクトルも描く。

## 完全ネットワーク

4状態 S_0...S_3 を4頂点とし、K4 の6辺を全て描く。

辺長は正規化された2波複素状態空間 C^2 の距離

d_ij = sqrt(|U^i-U^j|^2 + |U^(mi)-U^(mj)|^2)

とする。

3D座標はこの6辺距離から classical MDS で再構成するため、
見栄えだけの模式配置ではない。

## 生成ファイル（各 m=2,3）

- waves_base_harmonic_m*.png / .svg
- composite_wave_m*.png / .svg
- complex_vectors_all_states_m*.png / .svg
- complex_vectors_theta_pi4_m*.png / .svg
- K4_complete_network_m*.png / .svg
- edges_m*.csv

共通:
- all_edges_m2_m3.csv
- generate_all_K4_harmonic_figures_20260915.py
- README.md
