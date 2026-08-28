# N=3〜16 部分ゼロ閉包・複数閉包分解解析

## 定義

全体系の関係波を `z_e` とし、複素二乗距離 `w_e=z_e^2` を使う。部分集合 B の閉包残差を

`C(B)=|sum_{e in B} w_e| / sum_{e in B}|w_e|`

とした。C=0 がその部分集合単独の二乗ゼロ閉包である。物理更新は一切変更せず、各Nの既存step=5000生データだけを再解析した。

## 実施した探索

- 全体系のゼロ閉包残差
- 全2辺部分集合の総当たり
- 全3辺部分集合の総当たり
- C<1e-6 の3辺閉包だけを用いた全辺exact-cover探索（最大100解）
- K_N の標準round-robin matching分解
- 奇数NについてWalecki Hamilton-cycle分解


## 最重要結果

この探索では、全体系のゼロ閉包だけでなく、**全辺を複数の独立な部分ゼロ閉包へ完全分割できるNが実際に見つかった。**

- **N=4, M=6:** 3辺ゼロ閉包が8個存在し、そのうち全6辺を **3+3** の2閉包へ分割するexact coverを4通り確認。各ブロック残差は約 `4e-15`。
- **N=5, M=10:** 2辺ゼロ閉包が13個存在し、全10辺を **2+2+2+2+2** の5閉包へ分割するexact coverを12通り確認。各ブロック残差は最大でも約 `1e-8`（探索閾値 `1e-6`）。
- **N=6〜16:** 2辺・3辺の全組合せ総当たりでは `C<1e-6` の部分閉包を確認しなかった。したがって少なくとも2辺/3辺閉包という意味で、N=4,5は明瞭な特異点である。4辺以上の閉包の有無は今回の探索範囲外。

これは「Mが因数分解できるから自動的に分割される」のではなく、**実際の `z_e^2` 配置が部分和ゼロを許すNだけで複数閉包が成立する**ことを示している。

## 主要結果

- N=3, M=3: total=1.800e-15, best2=4.760e-01 (hits=0, covers=0), best3=1.800e-15 (hits=1), triple exact covers=1
- N=4, M=6: total=4.107e-15, best2=5.000e-01 (hits=0, covers=0), best3=4.020e-15 (hits=8), triple exact covers=4
- N=5, M=10: total=2.683e-15, best2=1.880e-10 (hits=13, covers=12), best3=3.333e-01 (hits=0), triple exact covers=0
- N=6, M=15: total=3.906e-15, best2=1.763e-01 (hits=0, covers=0), best3=1.182e-01 (hits=0), triple exact covers=0
- N=7, M=21: total=3.998e-15, best2=1.371e-02 (hits=0, covers=0), best3=1.934e-02 (hits=0), triple exact covers=0
- N=8, M=28: total=4.790e-15, best2=5.122e-04 (hits=0, covers=0), best3=3.291e-03 (hits=0), triple exact covers=0
- N=9, M=36: total=3.500e-16, best2=3.222e-03 (hits=0, covers=0), best3=1.260e-01 (hits=0), triple exact covers=0
- N=10, M=45: total=4.254e-15, best2=2.735e-04 (hits=0, covers=0), best3=6.252e-03 (hits=0), triple exact covers=0
- N=11, M=55: total=5.242e-15, best2=4.466e-05 (hits=0, covers=0), best3=2.539e-02 (hits=0), triple exact covers=0
- N=12, M=66: total=4.932e-16, best2=1.007e-05 (hits=0, covers=0), best3=2.102e-02 (hits=0), triple exact covers=0
- N=13, M=78: total=6.451e-16, best2=5.008e-04 (hits=0, covers=0), best3=1.002e-02 (hits=0), triple exact covers=0
- N=14, M=91: total=3.966e-15, best2=1.558e-04 (hits=0, covers=0), best3=7.225e-03 (hits=0), triple exact covers=0
- N=15, M=105: total=1.103e-15, best2=2.544e-05 (hits=0, covers=0), best3=2.757e-03 (hits=0), triple exact covers=0
- N=16, M=120: total=4.337e-16, best2=1.325e-04 (hits=0, covers=0), best3=8.069e-03 (hits=0), triple exact covers=0

## 解釈上の注意

1. Mの整数因数分解だけでは部分ゼロ閉包は保証されない。必要なのは、実際の複素ベクトル `z_e^2` の部分和がゼロになること。
2. canonical matching / Hamilton分解で残差がゼロでないことは「別の頂点置換・別の組合せ分解にも解がない」ことの証明ではない。
3. 一方、2辺・3辺探索はそのサイズについて全組合せを総当たりしているので、step=5000データに対する存在判定として直接的である。
4. triple exact cover が見つかったNでは、全体系を複数の3辺ゼロ閉包へ実際に分割できる（tol=1e-6の数値判定）。

## 次の判定軸

3辺で閉じないNについては4辺以上の部分閉包探索が次段階になる。Mが大きいNでは全組合せ総当たりが急増するため、meet-in-the-middleまたはexact-cover最適化を用いるのが適切。
