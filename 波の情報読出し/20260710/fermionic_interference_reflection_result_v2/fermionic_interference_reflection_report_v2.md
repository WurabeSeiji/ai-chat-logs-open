# フェルミオン的逆相核による干渉反射 実行結果 v2

## 目的

前回の完全弾性反射実験で手続き的に置いていた `q_A -> -q_A`, `q_B -> -q_B` を使わず、片側から入射する単一波束が、交換干渉チャンネルの位相だけで反射へ変換されるかを確認した。

## 判定

| 項目 | 結果 |
|---|---:|
| 外部 `q` 反転代入 | `false` |
| 片側入射初期条件 | `true` |
| `Delta_F=0` で通過 | `true` |
| `Delta_F=pi/2` で半反射 | `true` |
| `Delta_F=pi` で完全反射 | `true` |
| 局所写像位相掃引が期待式と一致 | `true` |
| 局所写像ノルム保存 | `true` |
| `Delta_F=pi` の交換干渉節 | `true` |
| 位相掃引が解析式と一致 | `true` |
| 交換経路が節形成に必要 | `true` |
| 識別振動 A 保存 | `true` |
| 識別振動 B 保存 | `true` |
| `U(pi)^2` 可逆性 | `true` |
| `U(delta)U(-delta)` 可逆性 | `true` |
| 補償付き二乗閉鎖 | `true` |
| AB/C セル差し替え | `true` |
| 最小機構判定 | `true` |

## 片側入射からの局所交換干渉反射

初期状態は左側だけに置いた単一局在波束であり、鏡像初期条件は使っていない。自由伝播で相互作用領域へ到達させ、局所窓内で偶・奇チャンネルに分解し、奇チャンネルへ内部核位相 `Delta_F` だけを与えて再合成した。

本実装は、連続ハミルトニアンを細かい時間刻みで積分する方式ではない。自由伝播で相互作用領域へ到達した波束に、局所交換干渉写像を一度作用させ、その後ふたたび自由伝播させる方式である。

| 量 | 値 |
|---|---:|
| `Delta_F=0` 反射率 | `1.8261693486616611e-19` |
| `Delta_F=0` 通過率 | `1.0000000000000002e+00` |
| `Delta_F=pi/2` 反射率 | `5.0000000000000000e-01` |
| `Delta_F=pi/2` 通過率 | `5.0000000000000022e-01` |
| `Delta_F=pi` 反射率 | `1.0000000000000004e+00` |
| `Delta_F=pi` 通過率 | `1.7939211304199106e-19` |
| 位相掃引最大誤差 | `5.5511151231257827e-16` |
| ノルム最大誤差 | `6.6613381477509392e-16` |

![single-sided scattering](fermionic_interference_single_sided_scattering_v2.png)

## 交換干渉節

直接経路と交換経路を

```math
\Psi_\Delta(1,2)=\frac{1}{\sqrt2}\left[P_A(1)P_B(2)+e^{i\Delta_F}P_A(2)P_B(1)\right]
```

として合成した。完全重なり対角 `1=2` では、`Delta_F=pi` で対角ノルムが消える。

| 量 | 値 |
|---|---:|
| `Delta_F=pi` 対角相対ノルム | `7.4987989133092880e-33` |
| 位相掃引最大誤差 | `8.8817841970012523e-16` |
| 交換あり `Delta_F=pi` | `7.4987989133092880e-33` |
| 交換なし `Delta_F=pi` | `4.9999999999999994e-01` |

![phase sweep](fermionic_interference_phase_sweep_v2.png)

## 相対位相方向の反射読出し

対照として単一自由波束を同じ相対位相座標で進めると、左半線確率はほぼゼロへ移り、通過する。一方、逆相核に対応する奇関数節を置いた波は、原点節を保ったまま左半線読出しで戻り、左向き電流へ反転する。

| 量 | 値 |
|---|---:|
| 単一自由波束の最終左確率 | `8.6915330692793126e-05` |
| 奇関数節の最終左確率 | `4.9999999999999994e-01` |
| 奇関数節の最終左電流 | `-1.1980224685843772e+00` |
| 奇関数節の最大節振幅 | `3.0038078125523204e-17` |
| 奇関数節の最大節電流 | `3.7549762132062191e-17` |

![relative dynamics](fermionic_interference_relative_dynamics_v2.png)

## 識別振動チャネル

識別振動 `eta` は反射生成チャネルではなく、保存される読出しチャネルとして扱った。これは、異なる `m_A,m_B` を交換消去チャネルへ直接混ぜると空間重なり節が壊れるためである。今回の実装では、反射はフェルミオン的逆相核の交換干渉で生成し、A/B の識別は `eta` 相関で読出す。

| 項目 | 読出し |
|---|---:|
| A の検出モード | `1` |
| B の検出モード | `2` |
| A ターゲット振幅 | `1.0000000000000000e+00` |
| B ターゲット振幅 | `1.0000000000000000e+00` |

## 可逆性と補償付き二乗閉鎖

局所交換干渉写像を二回適用し、波形が元へ戻るかを測定した。また、各サンプル係数 `x_n` に補償対 `i x_n` を付けた二乗閉鎖を、主要段階で評価した。

| 量 | 値 |
|---|---:|
| `U(pi)^2` 相対誤差 | `4.8214412843768590e-11` |
| `U(delta)U(-delta)` 最大相対誤差 | `2.2454514008125022e-16` |
| 補償付き二乗閉鎖 最大残差 | `1.2143074258005000e-17` |

## AB/C 相互作用セル差し替え

前回の AB/C 完全弾性衝突シミュレーションにおけるセル内の直接 `q` 反転を使わず、局所交換干渉写像から得た `R,T` により、

```text
q_out = q_in * (T - R)
```

として進行方向読出し量を生成した。

| 項目 | 値 |
|---|---:|
| AB/C 差し替え成立 | `true` |
| collision_cell_reached | `true` |
| post_collision_completed | `true` |
| q_A before/after | `1.0000000000000000e+00` / `-1.0000000000000004e+00` |
| q_B before/after | `-1.0000000000000000e+00` / `1.0000000000000004e+00` |
| label A initial/final | `1` / `1` |
| label B initial/final | `2` / `2` |

![ab-c replacement](fermionic_interference_ab_c_replacement_v2.png)

## 出力ファイル

| 種類 | ファイル |
|---|---|
| JSON | [fermionic_interference_reflection_result_v2.json](fermionic_interference_reflection_result_v2.json) |
| 位相掃引 CSV | [fermionic_interference_phase_sweep_v2.csv](fermionic_interference_phase_sweep_v2.csv) |
| 相対時間発展 CSV | [fermionic_interference_relative_dynamics_v2.csv](fermionic_interference_relative_dynamics_v2.csv) |
| 位相掃引図 | [fermionic_interference_phase_sweep_v2.png](fermionic_interference_phase_sweep_v2.png) |
| 相対時間発展図 | [fermionic_interference_relative_dynamics_v2.png](fermionic_interference_relative_dynamics_v2.png) |
| 片側入射散乱 CSV | [fermionic_interference_single_sided_scattering_v2.csv](fermionic_interference_single_sided_scattering_v2.csv) |
| 片側入射散乱図 | [fermionic_interference_single_sided_scattering_v2.png](fermionic_interference_single_sided_scattering_v2.png) |
| 可逆性 CSV | [fermionic_interference_reversibility_sweep_v2.csv](fermionic_interference_reversibility_sweep_v2.csv) |
| AB/C 差し替え CSV | [fermionic_interference_ab_c_replacement_timeline_v2.csv](fermionic_interference_ab_c_replacement_timeline_v2.csv) |
| AB/C 差し替え図 | [fermionic_interference_ab_c_replacement_v2.png](fermionic_interference_ab_c_replacement_v2.png) |
