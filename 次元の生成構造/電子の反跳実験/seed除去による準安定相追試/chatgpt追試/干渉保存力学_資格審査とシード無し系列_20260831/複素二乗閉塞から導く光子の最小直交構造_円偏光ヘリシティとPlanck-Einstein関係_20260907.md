# 複素二乗閉塞から導く光子の最小直交構造

## ― 円偏光・ヘリシティと Planck--Einstein 関係への接続 ―

**著者:** 木原範昭\
**ORCID:** 0009-0004-6753-4020\
**DOI:** （取得後記入）\
**Concept DOI:** （取得後記入）\
**日付:** 2026-09-07

## 要旨

本稿では複素波集合の実中心投影条件 $$
\sum_{k=1}^{n}z_k^2=R^2,\qquad R\in\mathbb R
$$
を出発点として光子の最小構成を検討する。当初、局在しない光子は単一の複素数波で表現できると想定した。しかし光子に曲率半径を持たないゼロ閉塞条件
$R=0$
を課すと、単一の非零複素数では条件を満たせない。最小の非自明構成は二つの複素数であり、
$$
z_1^2+z_2^2=0
$$ から $$
z_2=\pm iz_1
$$ が得られる。したがって二成分は等振幅かつ直交し、相対位相は $\pm\pi/2$
である。

この二つの解を二次元回転に対して調べると、それぞれが回転作用素の固有状態となり
$e^{\mp i\varphi}$
の位相因子を得る。この構造は左右円偏光および光子の二つのヘリシティ
$h=\pm1$ の回転表現と一致する。さらに回転位相を時間読み出し
$\varphi=\omega t$ に写すと $e^{-i\omega t}$
が得られ、標準量子論の時間発展との対応から $E=\hbar\omega=h\nu$
に接続する。

本構成では Maxwell 方程式、電場・磁場、光子スピン、Planck
定数を出発点に置かず、複素二乗閉塞の最小非自明解から直交二成分構造を得る。

## 1. 出発点

本研究では物理状態を複素数波の集合として表し、その幾何学的閉包を $$
\sum_k z_k^2=R^2,\qquad R\in\mathbb R
$$ によって記述する。$R$
は複素ノルムではなく、複素二乗和を実数へ閉じた中心投影半径である。

$z_k=a_k+ib_k$ とすれば $$
\sum_k(a_k^2-b_k^2)+2i\sum_k a_kb_k=R^2,
$$ したがって $$
\sum_k a_kb_k=0,\qquad
\sum_k(a_k^2-b_k^2)=R^2.
$$

光子候補について $R=0$ を考えると $$
\sum_k z_k^2=0
$$ となる。

## 2. 単一複素数仮説の棄却と最小二成分

一成分なら $$
z^2=0
$$ であり、複素数体上の解は $z=0$ のみである。したがって $z\neq0$ かつ
$R=0$ の非自明な一成分状態は存在しない。

二成分では $$
z_1^2+z_2^2=0
$$ より $$
z_2^2=-z_1^2
$$ であるから $$
\boxed{z_2=\pm iz_1}.
$$ 従って $$
|z_1|=|z_2|,\qquad
\boxed{\Delta\phi=\pm\frac{\pi}{2}}.
$$

任意の非零 $z$ に対して $(z,+iz)$ と $(z,-iz)$
が解であるため、固定される本質は絶対振幅や絶対位相ではなく二成分間の直交関係である。

## 3. 曲率半径ゼロの光子

二成分全体は $$
z_1^2+z_2^2=0
$$ なので $$
R^2=0.
$$
この状態は非零の内部関係を持ちながら有限の中心投影半径を持たない。本構成では光子単独に有限の大きさ、絶対位置、絶対位相を割り当てない。定義される最小情報は
$$
\boxed{\text{二つの複素成分が直交関係にある}}
$$ ことである。

## 4. 回転から円偏光とヘリシティ

規格化して $$
\Psi_\pm=\frac1{\sqrt2}
\begin{pmatrix}1\\ \pm i\end{pmatrix}
$$ とする。二次元回転 $$
\mathcal R(\varphi)=
\begin{pmatrix}
\cos\varphi&-\sin\varphi\\
\sin\varphi&\cos\varphi
\end{pmatrix}
$$ を作用させると $$
\boxed{\mathcal R(\varphi)\Psi_+=e^{-i\varphi}\Psi_+},
$$ $$
\boxed{\mathcal R(\varphi)\Psi_-=e^{+i\varphi}\Psi_-}.
$$

したがって複素二乗ゼロ閉塞から得た二つの直交解は回転固有状態であり、無次元回転固有値は
$\pm1$
である。これは左右円偏光、および伝播軸に沿った光子の二つのヘリシティ
$h=\pm1$ の回転表現と一致する。

導出系列は $$
\boxed{
\sum z_k^2=0
\Longrightarrow
z_2=\pm iz_1
\Longrightarrow
\Delta\phi=\pm\frac{\pi}{2}
\Longrightarrow
e^{\mp i\varphi}
\Longrightarrow
h=\pm1
}
$$ となる。

## 5. 時間位相と Planck--Einstein 関係

回転位相が時間読み出しに対して一定角速度 $\omega$ で進むとする。 $$
\varphi(t)=\omega t.
$$ すると $$
\boxed{\Psi(t)=e^{-i\omega t}\Psi(0)}
$$ である。この段階まで Planck 定数は必要ない。

標準量子論のエネルギー固有状態は $$
\Psi(t)=e^{-iEt/\hbar}\Psi(0)
$$ と書かれる。両者の位相を対応させれば $$
\omega=\frac{E}{\hbar},
\qquad
E=\hbar\omega.
$$ さらに $$
\omega=2\pi\nu,\qquad \hbar=\frac{h}{2\pi}
$$ より $$
\boxed{E=h\nu}.
$$

したがって本構成で基礎的に得られるのは位相回転率 $\omega$
であり、$E=h\nu$
はその位相回転と標準的エネルギー読み出しとの対応として現れる。

空間回転の二符号はヘリシティに対応するが、両ヘリシティの正周波数光子についてエネルギーは
$E=+\hbar\omega$ である。

## 6. 既存研究との関係

Riemann--Silberstein ベクトルは電磁場を $$
\mathbf F_\pm\propto\mathbf E\pm i\mathbf B
$$ として複素化し、Maxwell
方程式と光子の二ヘリシティを簡潔に記述する。Białynicki-Birula
はこの構造を光子波動関数として体系化し、複素構造と左右ヘリシティの対応を論じている
\[1,2\]。Elbistan, Horváthy, Zhang は Riemann--Silberstein
型光子波動関数における双対性とヘリシティを扱う \[3\]。

また複素 null vector
と円偏光の関係自体には既存の記述があり、等しい直交実成分を持つ複素ベクトルが
null となる構造が論じられている \[4\]。

本稿では電場 $\mathbf E$、磁場 $\mathbf B$、Maxwell
方程式、円偏光、ヘリシティを出発条件とせず、 $$
\sum_kz_k^2=R^2
$$ という複素二乗閉塞から出発する。その $R=0$
の最小非自明解として二成分直交状態を得て、その回転表現から円偏光と二ヘリシティ構造へ接続する。

## 7. 結論

複素二乗の実中心投影条件 $$
\sum_kz_k^2=R^2
$$ に対し光子候補として $R=0$
を要求すると、単一非零複素数では解が存在しない。最小非自明解は二成分であり
$$
\boxed{z_2=\pm iz_1}
$$ となる。

したがって光子の最小構成は、絶対位置・有限の大きさ・絶対位相を必要とせず、二つの複素成分の直交関係として記述される。この二状態は
$$
\mathcal R(\varphi)\Psi_\pm=e^{\mp i\varphi}\Psi_\pm
$$ と変換し、左右円偏光および二ヘリシティ $h=\pm1$
の回転構造に一致する。

さらに $\varphi=\omega t$ とすると $$
\Psi(t)=e^{-i\omega t}\Psi(0)
$$ となり、この段階まで Planck
定数は不要である。標準量子論の時間発展との対応から $$
E=\hbar\omega=h\nu
$$ へ接続する。

従って $$
\boxed{
\text{複素二乗ゼロ閉塞}
\rightarrow
\text{最小二成分}
\rightarrow
\text{直交}
\rightarrow
\text{円偏光}
\rightarrow
\text{二ヘリシティ}
\rightarrow
\text{時間位相}
\rightarrow
E=h\nu
}
$$ という短い導出系列を得る。

## 参考文献

\[1\] I. Białynicki-Birula, "On the Wave Function of the Photon,"
*Progress in Optics*, 36, 245--294 (1996). DOI:
10.1016/S0079-6638(08)70316-0.

\[2\] I. Białynicki-Birula and Z. Białynicka-Birula, "The role of the
Riemann--Silberstein vector in classical and quantum theories of
electromagnetism," *J. Phys. A: Math. Theor.* 46, 053001 (2013). DOI:
10.1088/1751-8113/46/5/053001.

\[3\] M. Elbistan, P. A. Horváthy, P.-M. Zhang, "Duality and helicity:
the photon wave function approach," *Physics Letters A* 381, 2375--2379
(2017). DOI: 10.1016/j.physleta.2017.05.042.

\[4\] D. H. Delphenich, "Algebra of Complex Vectors and Applications in
Electromagnetic Theory and Quantum Mechanics," *Mathematics* 3(3),
781--842 (2015). DOI: 10.3390/math3030781.
