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

## 主要結果

- N=3, M=3: total=1.021e-13, best2=3.333e-01 (hits=0, covers=0), best3=1.021e-13 (hits=1), triple exact covers=1
- N=4, M=6: total=6.938e-14, best2=3.333e-01 (hits=0, covers=0), best3=1.653e-13 (hits=12), triple exact covers=6
- N=5, M=10: total=2.683e-15, best2=1.880e-10 (hits=13, covers=12), best3=3.333e-01 (hits=0), triple exact covers=0
- N=6, M=15: total=4.207e-14, best2=3.527e-01 (hits=0, covers=0), best3=2.792e-02 (hits=0), triple exact covers=0
- N=7, M=21: total=4.607e-14, best2=1.418e-01 (hits=0, covers=0), best3=4.781e-02 (hits=0), triple exact covers=0
- N=8, M=28: total=6.020e-14, best2=3.839e-02 (hits=0, covers=0), best3=2.862e-02 (hits=0), triple exact covers=0
- N=9, M=36: total=5.893e-14, best2=9.155e-02 (hits=0, covers=0), best3=1.251e-02 (hits=0), triple exact covers=0
- N=10, M=45: total=3.304e-14, best2=4.137e-02 (hits=0, covers=0), best3=1.196e-02 (hits=0), triple exact covers=0
- N=11, M=55: total=4.303e-14, best2=2.838e-02 (hits=0, covers=0), best3=3.833e-03 (hits=0), triple exact covers=0
- N=12, M=66: total=6.338e-14, best2=4.903e-02 (hits=0, covers=0), best3=7.171e-03 (hits=0), triple exact covers=0
- N=13, M=78: total=2.457e-14, best2=2.978e-02 (hits=0, covers=0), best3=1.559e-03 (hits=0), triple exact covers=0
- N=14, M=91: total=4.138e-14, best2=1.372e-02 (hits=0, covers=0), best3=2.036e-03 (hits=0), triple exact covers=0
- N=15, M=105: total=3.612e-14, best2=1.253e-02 (hits=0, covers=0), best3=1.045e-03 (hits=0), triple exact covers=0
- N=16, M=120: total=4.702e-14, best2=1.792e-02 (hits=0, covers=0), best3=2.132e-04 (hits=0), triple exact covers=0

## 解釈上の注意

1. Mの整数因数分解だけでは部分ゼロ閉包は保証されない。必要なのは、実際の複素ベクトル `z_e^2` の部分和がゼロになること。
2. canonical matching / Hamilton分解で残差がゼロでないことは「別の頂点置換・別の組合せ分解にも解がない」ことの証明ではない。
3. 一方、2辺・3辺探索はそのサイズについて全組合せを総当たりしているので、step=5000データに対する存在判定として直接的である。
4. triple exact cover が見つかったNでは、全体系を複数の3辺ゼロ閉包へ実際に分割できる（tol=1e-6の数値判定）。

## 次の判定軸

3辺で閉じないNについては4辺以上の部分閉包探索が次段階になる。Mが大きいNでは全組合せ総当たりが急増するため、meet-in-the-middleまたはexact-cover最適化を用いるのが適切。
