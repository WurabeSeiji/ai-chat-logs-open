# 02 統合定義

## 物理散乱層

512×16全状態に対し、`a_raw=r_eff*a+t_eff*b`, `b_raw=t_eff*a+r_eff*b` を適用した。C0は `theta_eff=theta0`、反転Candidate 1は `theta_eff=theta0-kappa*rho(theta0)*(c_A+c_B)/2` であり、候補差は角度だけである。型名による条件分岐はない。

## 読出し層

A由来 `(q=+1,m_eta=1)` とB由来 `(q=-1,m_eta=2)` をeta射影し、それぞれの搬送波を除去してから半周期相関 `c_pi` と偶奇射影重みを計算した。全チャネルへの単一復調は行っていない。

## 更新系列

`existing_normalization` は正本と同じチャネル別正規化、`raw_update` はraw出力をそのまま次状態にした。raw全状態と経路ノルムを物理量、スペクトル類似度を診断量として分離した。
