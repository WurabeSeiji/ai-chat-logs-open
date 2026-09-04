# N=3..6 段階1 exact 検証

倍精度の canonical stage1 最終状態から、位相を最寄りの {1,i,-1,-i}、振幅二乗を有理数候補へ復元し、その候補を SymPy の厳密演算で元の W,J,G,K 方程式へ代入した。したがって以下の residual=0 は浮動小数点近似ではなく記号的恒等式。

|N|sigma^2|norm^2|r^2|N r^2|iKv=sigma v|Wy=v|JG y=-i sigma y|norm y=1|
|---:|---:|---:|---:|---:|:---:|:---:|:---:|:---:|
|3|2|1/5|1/15|1/5|True|True|True|False|
|4|8|2/3|1/9|4/9|True|True|True|False|
|5|14|49/57|49/570|49/114|True|True|True|False|
|6|24|6/5|2/25|12/25|True|True|True|False|

注意: これは候補の exact 検証であり、固定点枝の一意性証明ではない。枝の一意性/多重性は別の seed sweep で検査する。
