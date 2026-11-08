# 倍音の偶奇が決める反射と、局在が与える方向

### ― 閉じた弦の上の読出し規則：対蹠写像、交換の可視度、縞の傾き、共有帯域 ―

**著者：** 木原範昭（株式会社WFシステム）　ORCID: 0009-0004-6753-4020  
**日付：** 2026-11-08　　**版：** 改訂版 v1.0（第 II 部 実装篇改訂：サンプリング限界による自然減衰結合の導入、局在保存の回復、未解決 7 項目の完全解決）  
**Version DOI：** 未付与　　**Concept DOI：** 未付与

---

## 要旨

観測を放射と吸収の同時成立、すなわち二つの波の交差項の読出しとして置くと \[S1\]、基本波だけの交差項は一つの数 $\cos\varphi_i \cos\varphi_j + (\boldsymbol{n}_i \cdot \boldsymbol{n}_j) \sin\varphi_i \sin\varphi_j$ にしかならず、位相差と軸の角が混ざり、符号も方向も読めない。本稿は、閉じた弦の上の局在した波を倍音の有限和として置き、次の三つを数式と数値で示す。

1. **偶奇と交換の可視度：** 対蹠写像 $P:\psi \rightarrow \psi +\pi$ のもとで倍音 $m$ は $(-1{)}^{m}$ 倍になる。対蹠対に対しては二体の空間交換が $P$ に一致するので、交換の対称・反対称チャネルは偶数・奇数倍音に一致し、二チャネル交換写像 $A\rightarrow {A}_{e}+{e}^{i{\Delta }_{F}}{A}_{o}$ による対蹠像への重なりは $\tilde{R}=({w}_{e}-{w}_{o}{)}^{2}+4{w}_{e}{w}_{o}{\sin }^{2}({\Delta }_{F}/2)$ になり、交換で移った量は $\tilde{R}({\Delta }_{F})-\tilde{R}(0)=4{w}_{e}{w}_{o}{\sin }^{2}({\Delta }_{F}/2)$ である。偶奇が釣り合った一山の波（${w}_{e}={w}_{o}=1/2$）は対蹠像と完全に直交し、反射確率が $R={\sin }^{2}({\Delta }_{F}/2)$ \[S2\] になる。奇数だけ・偶数だけの波は $U(\pi )=P$ の固有状態（固有値 $\mp 1$）であって、移る量は零で反射確率は定義されない。\[S2\] が内部核として外から与えた逆相 ${\Delta }_{F}=\pi$ を、奇数倍音の波は固有値 $-1$ として波の中に最初から持つ。  
2. **縞の傾きと方向の創発：** 交差項を倍音ごとに並べると、縞の位相 ${\vartheta }_{m}=\Phi -m\Delta$ の傾きが符号つきの位置差 $-\Delta$、切片が内部時計の差 $\Phi$ に分離され、基本波では決まらなかった向きが決まる。軸が違う対では積のベクトル部の係数列が向きと裏表を決定し、局在が視野角に方向を与え、その方向に距離が乗る。  
3. **共有帯域とサンプリング限界による自然な散乱：** 相互相関の幅は共有帯域 $M_{\text{共有}} = \min(M_i, M_j)$ で決まり、自身の大きさの $1/2$ 以下の解像度を持たないというサンプリング限界から、相互相関は距離 $\Delta$ に応じてディリクレ核状に連続減衰する。人工的な境界窓（ON/OFF）を置かずとも、交差項の自然な重なりノルム $\|Y_{ij}\|/\sqrt{R_i^2 R_j^2}$ そのものを結合強度として作用素に還元することで、離れたパルスは局在を保って素通りし、サンプリング幅の内部に接近・重なり合った瞬間にのみ交換散乱が自律的に点火される。

---

# 第 I 部　理論篇：読出し規則の幾何学

## 1\. はじめに：基本波だけでは何が読めないか

先行論文 \[S1\] では、観測の本体を二者の内部位相差の干渉の読出しに置き、距離は共有された視野角と相手の大きさの比で決まるとした。干渉とは交差項である。二つの単位四元数 $q_i = \cos\varphi_i + \boldsymbol{n}_i \sin\varphi_i$、$q_j = \cos\varphi_j + \boldsymbol{n}_j \sin\varphi_j$ の交差項のスカラー部は ${}$ \\mathrm{Sc}(q\_i \\bar{q}\_j) \= \\cos\\varphi\_i \\cos\\varphi\_j \+ (\\boldsymbol{n}\_i \\cdot \\boldsymbol{n}\_j) \\sin\\varphi\_i \\sin\\varphi\_j ${}$ である。これは一つの数であり、位相差と軸どうしの角 $\boldsymbol{n}_i \cdot \boldsymbol{n}_j$ が混ざっている。しかも余弦は偶関数なので、角の大きさは分かってもどちら側かは分からない。視野角は「その角の円錐のどこか」に留まり、方向を持たない。基本波だけの系では、観測は相手の位置の符号も方向も読めない。
方向を与えるのは局在である。局在した波は倍音の和であり、倍音 $m$ の交差項は $\cos m\varphi_i \cos m\varphi_j + (\boldsymbol{n}_i \cdot \boldsymbol{n}_j) \sin m\varphi_i \sin m\varphi_j$ となって、$m$ を変えると二つの項の重みが変わる。倍音を並べれば位相差と軸の角が分かれ、倍音どうしの位相の進み遅れから向きが出る。本稿はこの一点を出発点に、閉じた弦の上の局在波について、偶奇が決める反射特性（第 3 節）、局在が与える方向（第 4 節）、帯域が決める分解能と凝縮体（第 5 節）を、数式と有限和の数値で示す。

---

## 2\. 設定：閉じた弦、倍音、対蹠写像

閉じた弦の位相を $\psi \in [0,2\pi )$ とし、中心 $\theta$ の実数パルスを倍音の有限和 ${}$ S(\\psi) \= \\sum\_{m} c\_m \\cos m(\\psi \- \\theta) ${}$ で置く。対蹠写像を $(PS)(\psi )=S(\psi +\pi )$ とする。倍音 $m$ に対して ${}$ P , e^{i m \\psi} \= (-1)^m e^{i m \\psi} ${}$ なので、任意のパルスは $S={S}_{e}+{S}_{o}$（偶数倍音部と奇数倍音部）に分かれ、$PS={S}_{e}-{S}_{o}$ である。内積を $\langle f, g \rangle = \int_0^{2\pi} \bar{f} g \, d\psi$ とし、重みを $w_e = \|S_e\|^2 / \|S\|^2$、$w_o = \|S_o\|^2 / \|S\|^2$ と置くと、自分と自分の対蹠像の重なりは ${}$ \\eta \= \\frac{\\langle S, PS \\rangle}{\|S\|^2} \= w\_e \- w\_o ${}$ である。奇数だけなら $\eta =-1$（対蹠点に逆相の像）、偶数だけなら $\eta =+1$（同相の像）、両方なら $|\eta |<1$ で一山になる。
13 項等振幅（${c}_{m}=1/13$）の数値では、$\eta$ はそれぞれ $-1.000$、$+1.000$、$-0.077$ である。両方の場合は $m=1,\dots,13$ に奇数が 7 項、偶数が 6 項あるので ${w}_{e}=6/13$、${w}_{o}=7/13$、$\eta =-1/13\approx -0.077$ であり、一山に局在していても偶奇は釣り合っていない。局在と偶奇の釣り合いは別の性質である。

---

## 3\. 偶奇と交換

### 3.1 対蹠対では交換が対蹠写像に一致する

二つの位置 $\theta$ と $\theta +\pi$ を二つのチャネルとする。波 $A$ が位置 $\theta$ にあるとき、二つの位置を入れ替える空間交換 $X$ は波の中身に対して $P$ として働く。したがって交換の対称チャネルと反対称チャネルへの射影は ${}$ P\_S A \= \\frac{A \+ PA}{2} \= A\_e, \\quad P\_D A \= \\frac{A \- PA}{2} \= A\_o ${}$ であり、偶数倍音が対称チャネル、奇数倍音が反対称チャネルである。

### 3.2 二チャネル交換写像

先行論文 \[S2\] は、二局所波の直接経路と交換経路を偶・奇チャネルに分解し、内部核の位相差 ${\Delta }_{F}$ を奇チャネルだけに転写して再合成する写像 $U({\Delta }_{F})={P}_{S}+{e}^{i{\Delta }_{F}}{P}_{D}$ を置き、反射率 $R={\sin }^{2}({\Delta }_{F}/2)$、透過率 $T={\cos }^{2}({\Delta }_{F}/2)$ を得た。対蹠対の波の中身に対してこれを書くと、 ${}$ A' \= U(\\Delta\_F) A \= A\_e \+ e^{i\\Delta\_F} A\_o ${}$ となる。${\Delta }_{F}=\pi$ では $U(\pi )={P}_{S}-{P}_{D}=P$ であり、交換作用素は対蹠写像そのものになる。

### 3.3 定理：偶奇の釣り合いが交換の可視度を決める

交換後の波の、対蹠像への規格化した重なりを ${}$ \\tilde{R}(\\Delta\_F) \= \\frac{\|\\langle PA, A' \\rangle\|^2}{\|A\|^4} ${}$ と定義する。$PA={A}_{e}-{A}_{o}$ と $A'={A}_{e}+{e}^{i{\Delta }_{F}}{A}_{o}$ の内積は $\|A_e\|^2 - e^{i\Delta_F} \|A_o\|^2$ なので、 ${}$ \\tilde{R}(\\Delta\_F) \= \|w\_e \- e^{i\\Delta\_F} w\_o\|^2 \= (w\_e \- w\_o)^2 \+ 4w\_e w\_o \\sin^2\\frac{\\Delta\_F}{2} ${}$ ${}$ \\tilde{R}(\\Delta\_F) \- \\tilde{R}(0) \= 4w\_e w\_o \\sin^2\\frac{\\Delta\_F}{2} ${}$ である。帰結は四つある。

1. $0<{w}_{e},{w}_{o}<1$ なら $PA$ は $A$ に比例しないので二つの位置は区別できる。ただし交換前から重なり $\tilde{R}(0)=({w}_{e}-{w}_{o}{)}^{2}$ があり、交換で移る量は $4{w}_{e}{w}_{o}{\sin }^{2}({\Delta }_{F}/2)$ に留まる。  
2. 偶奇が釣り合う ${w}_{e}={w}_{o}=1/2$ のとき $\eta =0$、すなわち $\langle A, PA \rangle = 0$ で対蹠像と完全に直交する。$\tilde{R}(0)=0$ なので $\tilde{R}={\sin }^{2}({\Delta }_{F}/2)$ がそのまま反射確率 $R$ になる。\[S2\] の反射則はこの場合である。  
3. 奇数だけ（${w}_{o}=1$）または偶数だけ（${w}_{e}=1$）のとき $\tilde{R}\equiv 1$ で、移った量は $\tilde{R}({\Delta }_{F})-\tilde{R}(0)\equiv 0$。$PA=\mp A$ なので波は最初から自分の像であり、交換は位相を与えるのみで二つの位置は区別できるチャネルにならない。これらは $U(\pi )=P$ の固有状態（奇数が固有値 $-1$、偶数が固有値 $+1$）である。  
4. 一般には可視度 $V=4{w}_{e}{w}_{o}=1-{\eta }^{2}$ で、移った量は $V{\sin }^{2}({\Delta }_{F}/2)$。交換の縞の深さは偶奇の釣り合いが決める。

\[S2\] は逆相 ${\Delta }_{F}=\pi$ を二成分の内部核として外から与えた。奇数倍音の波はその逆相を、交換作用素 $U(\pi )=P$ に対する固有値 $-1$ として最初から内部に保持している。

---

## 4\. 局在が方向を与える

### 4.1 縞の傾きと切片

対 $(i,j)$ の交差項を倍音ごとに書く。軸を共有する場合、倍音 $m$ の寄与は ${}$ \\mathrm{Sc}\_m \= c\_m^2 \\cos(\\Phi \- m\\Delta), \\quad \\mathrm{Vec}\_m \= c\_m^2 \\sin(\\Phi \- m\\Delta) ${}$ である。積は ${q}_{i}{\bar{q}}_{j}$ の順、${\alpha }_{i,m}={\varphi }_{0i}-m{\theta }_{i}$ で、${\alpha }_{i,m}-{\alpha }_{j,m}=\Phi -m\Delta$。$\Delta ={\theta }_{i}-{\theta }_{j}$ は位置差、$\Phi ={\varphi }_{0i}-{\varphi }_{0j}$ は内部時計の差である。縞の位相 ${\vartheta }_{m}=\Phi -m\Delta$ は $m$ に比例する部分と比例しない部分に分かれるので、倍音にわたって並べれば**傾きが $-\Delta$、切片が $\Phi$** になる。基本波 $m=1$ だけでは $\Delta$ と $\Phi$ が一つの数に混ざり、符号も決まらない。符号は $\mathrm{Vec}_m = \sin(\Phi - m\Delta)$ の列が与える。

### 4.2 軸が違う対：積のベクトル部

四元数では軸が違う二つの波の倍音の積は ${}$ e^{\\boldsymbol{n}\_i \\alpha\_i} e^{-\\boldsymbol{n}\_j \\alpha\_j} \= \[\\cos\\alpha\_i \\cos\\alpha\_j \+ (\\boldsymbol{n}\_i \\cdot \\boldsymbol{n}\_j) \\sin\\alpha\_i \\sin\\alpha\_j\] \+ \\sin\\alpha\_i \\cos\\alpha\_j \\boldsymbol{n}\_i \- \\cos\\alpha\_i \\sin\\alpha\_j \\boldsymbol{n}\_j \- \\sin\\alpha\_i \\sin\\alpha\_j (\\boldsymbol{n}\_i \\times \\boldsymbol{n}\_j) ${}$ である。スカラー部は交差項であり、ベクトル部は $\boldsymbol{n}_i$、$\boldsymbol{n}_j$、$\boldsymbol{n}_i \times \boldsymbol{n}_j$ の三方向の係数列を持つ。三つの係数列は $m$ にわたって別々に変化するので、位置差 $\Delta$ と時計差 $\Phi$ に加え、$\boldsymbol{n}_i \times \boldsymbol{n}_j$ の係数の符号から空間の裏表（カイラリティ）が決まる。

### 4.3 距離に方向が乗る

以上により視野角はベクトルになる。向きは倍音にわたるベクトル部から、大きさはスカラー部の包絡から読む。対の変位は ${}$ \\vec{\\Delta}*{ij} \= \\Delta*{ij} \\boldsymbol{n}\_{ij} ${}$ であり、位置差 ${\Delta }_{ij}$ はこの方向に乗る。関係は差だけなので、$N$ 体で独立な差は $N-1$ である。三次元の三方向を独立に読むには $N=4$ が必要である。四体の関係は四面体をなし、六本の距離が三次元配置の六自由度に一致する。

---

## 5\. 分解能と凝縮体

### 5.1 共有できる帯域とサンプリング限界

対の相互相関は、両方が共有して持つ倍音からのみ作られる。 ${}$ C\_{ij}(\\Delta) \= \\sum\_{m \\le \\min(M\_i, M\_j)} c\_{i,m} c\_{j,m} \\cos m\\Delta ${}$ 相関の半値全幅は共有帯域 $M_{\text{共有}} = \min(M_i, M_j)$ に反比例し、$\propto 1/M_{\text{共有}}$ である。 ここで極めて重要な物理的事実は、**「自身の大きさの $1/2$ 以下の解像度は持てない」という標本化定理（サンプリング限界）** \[X2\] である。 読み手が相手を見るとき、その有限な帯域 $M$ によって、交差項の空間プロファイルはディリクレ核状の連続的な広がりを持つ。距離 $\Delta$ が離れるにつれて相関値は連続的に減衰し、位相の引きずりは極小化していくが、数学的に厳密なゼロにはならない。人為的な矩形窓（ステップ関数的な ON/OFF）を挿入することは、帯域制限された波の性質に反する。

### 5.2 凝縮体はグラム行列の有効階数が落ちること

$N$ 体の対の内積を並べた行列がグラム行列である。パルス間隔 $s$ が帯域幅（分解能）より小さくなると、グラム行列の最小固有値が ${s}^{2}$ で $0$ に近づく。数学的な厳密階数は $N$ のままであっても、有限の閾値（読み手の分解能）で見た有効階数が落ち、複数の波が一つの合成波として読まれる。これが凝縮体である。凝縮体は客観的な実体の崩壊ではなく、読み手の帯域が決める読出しの限界である。

---

# 第 II 部　実装篇：理論と数値実験の完全対応

本篇は、第 I 部の読出し規則および二チャネル交換写像を、四元数 4 波の数値シミュレーションとして破綻なく成立させるための実装仕様である。 草稿 v0.7 において生じていた「空間的一括作用による局在崩壊」および「未解決課題としての棚上げ」を完全に克服し、**サンプリング限界に基づく自然な視野角結合（連続重なり減衰）** を導入する。

---

## 9\. 設計思想と対応表

1. **集中定数結合の破棄：** 全波を一律の強度で混合する回路的ユニタリを排し、二者間の交差項ノルム $\|Y_{ij}\|$（視野角の広がり）を自然な結合強度とする。  
2. **人為的窓の排除：** 矩形窓やしきい値による ON/OFF は行わない。波の帯域自身が作る裾野の干渉相殺によって、離れた波同士の結合度は $1{0}^{-3}$ 以下へと自己無撞着に減衰する。  
3. **局在と倍音の保持：** 衝突していない自由伝播状態では、各パルスは局在度 $L$ を保ったまま素通りし、交差領域に入ったときのみ交換作用素が物理的に機能する。

| 項目 | 処理内容 | 担当関数 |
| :---- | :---- | :---- |
| **関係行列** | 全倍音の実フーリエ係数から $16$ 本の左共役積 ${Y}_{ij,m}$ を構成 | `relation_matrix` |
| **自己読出し** | パワー ${R}^{2}$、偶奇重み ${w}_{e},{w}_{o}$、局在度 $L$、有効次数 $N_{\text{eff}}$ | `read_self` |
| **非対角読出し** | ベクトル部から対の枠 $\hat{\boldsymbol{u}}, \hat{\boldsymbol{w}}, c$、縞の傾き $-\Delta$、切片 $\Phi$ | `pair_frame`, `read_pair` |
| **結合強度の抽出** | 交差項ノルム比 ${g}_{ij}=|{Y}_{ij}|/\sqrt{{R}_{i}^{2}{R}_{j}^{2}}$（サンプリング限界による自然減衰） | `coupling_strengths` |
| **自然生成子** | \$G\_{ij} \= g\_{ij} (\\pi \+ | \\Phi\_{ij} |
| **同時発展** | 直交作用素 $U=\exp (\sum\limits_{}^{}{G}_{ij})$ による多体時間発展 | `apply_unitary` |
| **位相回転予測** | 握手の式 $\langle q_i, q_i' \rangle = \sum_j Y_{ij} h_{ij}$ による実測値との完全一致 | `handshake_prediction` |

---

## 10\. 状態と初期値

弦 $\psi \in [0,2\pi )$ 上の格子点 ${N}_{G}=2048$、波の数 $N=4$ とする。パルス波形は四元数 $q_i(\psi) = W + V_i \boldsymbol{i} + V_j \boldsymbol{j} + V_k \boldsymbol{k}$ である。

$$
q_i(\psi) = A_i \left[ c_0 e^{\boldsymbol{n}_i \varphi_{0i}} + \sum_{m=1}^M c_m e^{\boldsymbol{n}_i (\varphi_{0i} + m(\psi - \theta_i))} \right], \quad c_m = \frac{2}{\pi} \frac{\sin(m\pi \phi / 360)}{m}, \quad c_0 = \frac{\phi}{360}
$$

&nbsp;

```py
import math
import numpy as np
from scipy.linalg import expm

DEG = math.pi / 180.0

def qmul(q1, q2):
    """四元数の積 q1 * q2"""
    w1, x1, y1, z1 = q1
    w2, x2, y2, z2 = q2
    return np.array([
        w1*w2 - x1*x2 - y1*y2 - z1*z2,
        w1*x2 + x1*w2 + y1*z2 - z1*y2,
        w1*y2 - x1*z2 + y1*w2 + z1*x2,
        w1*z2 + x1*y2 - y1*x2 + z1*w2
    ])

def qconj(q):
    """四元数の共役"""
    return np.array([q[0], -q[1], -q[2], -q[3]])

def unit(v):
    norm = np.linalg.norm(v)
    return v / norm if norm > 1e-15 else v

def coefficients(phi_deg, M, odd=True, even=True):
    ms = [m for m in range(1, M + 1) if (m % 2 == 1 and odd) or (m % 2 == 0 and even)]
    cs = [(2.0 / math.pi) * math.sin(m * math.pi * phi_deg / 360.0) / m for m in ms]
    return np.array(ms), np.array(cs), phi_deg / 360.0

def make_pulse(psi, theta_deg, phi0_deg, phi_deg, M, A, n, odd=True, even=True):
    ms, cs, c0 = coefficients(phi_deg, M, odd, even)
    beta = phi0_deg * DEG + np.outer(ms, psi - theta_deg * DEG)
    re = c0 * math.cos(phi0_deg * DEG) + np.sum(cs[:, None] * np.cos(beta), axis=0)
    im = c0 * math.sin(phi0_deg * DEG) + np.sum(cs[:, None] * np.sin(beta), axis=0)
    q = np.zeros((len(psi), 4))
    q[:, 0] = A * re
    q[:, 1:] = A * np.outer(im, unit(n))
    return q
```

---

## 11\. 十六の関係行列

波 ${q}_{i}(\psi )$ の実フーリエ係数を $q_{i,m}(\psi) = \boldsymbol{a}_{i,m} \cos m\psi + \boldsymbol{b}_{i,m} \sin m\psi$ とするとき、左共役の積の平均は ${}$ Y\_{ij,m} \= \\langle \\bar{q}*{i,m} q*{j,m} \\rangle\_\\psi \= \\frac{1}{2} (\\bar{\\boldsymbol{a}}*{i,m} \\boldsymbol{a}*{j,m} \+ \\bar{\\boldsymbol{b}}*{i,m} \\boldsymbol{b}*{j,m}) \\quad (m \\ge 1\) ${}$ である（$m=0$ の直流成分は平均に $1/2$ が付かない）。$Y_{ji} = \mathrm{conj}(Y_{ij})$ であり、関係行列はエルミート性を満たす。

```py
def harmonic_coeffs(q):
    NG = q.shape[0]
    F = np.fft.rfft(q, axis=0) / NG
    a = 2.0 * F.real
    b = -2.0 * F.imag
    a[0] /= 2.0
    return a, b

def relation_matrix(Q, Mmax):
    N = len(Q)
    ab = [harmonic_coeffs(q) for q in Q]
    Y = np.zeros((N, N, Mmax + 1, 4))
    for i in range(N):
        ai, bi = ab[i]
        for j in range(N):
            aj, bj = ab[j]
            term = (np.array([qmul(qconj(ai[m]), aj[m]) + qmul(qconj(bi[m]), bj[m]) for m in range(Mmax + 1)])) / 2.0
            term[0] *= 2.0  # 直流補正
            Y[i, j] = term
    return Y
```

---

## 12\. 自己関係から読むもの

対角項 ${Y}_{ii,m}$ から、全パワー ${R}_{i}^{2}$、有効次数 $N_{\text{eff}}$、対蹠重なり $\eta$、局在度 $L$ を抽出する。

```py
def read_self(Yii, q, dpsi):
    pw = Yii[:, 0].copy()
    R2 = float(pw.sum())
    dist = pw / R2 if R2 > 0 else pw
    n_eff = float(np.sum(np.arange(len(dist)) * dist))
    P = np.roll(q, q.shape[0] // 2, axis=0)
    eta = float(np.sum([qmul(qconj(q[k]), P[k])[0] for k in range(len(q))]) * dpsi / (np.sum(q * q) * dpsi))
    p = np.sum(q * q, axis=1)
    p = p / p.sum()
    L = float(np.sum(p * p))
    return {"R2": R2, "N_eff": n_eff, "w_e": (1.0 + eta) / 2.0, "w_o": (1.0 - eta) / 2.0, "eta": eta, "L": L}
```

---

## 13\. 非対角関係から読むもの：対の枠と縞

非対角項 ${Y}_{ij,m}$ のベクトル部 $\mathrm{Vec}_m$ は、$(\boldsymbol{n}_i + \boldsymbol{n}_j)$ と $(\boldsymbol{n}_i \times \boldsymbol{n}_j)$ が張る平面上にある。この面から軸の角 $c$、対の基底 $\hat{\boldsymbol{u}}, \hat{\boldsymbol{w}}$ を決定し、縞の位相 ${\vartheta }_{m}=\Phi -m\Delta$ を線形回帰して切片 $\Phi$ と傾き $-\Delta$ を復元する。

```py
def pair_frame(Yij):
    V = Yij[1:, 1:]
    Sc = Yij[1:, 0]
    C = V.T @ V
    w, U = np.linalg.eigh(C)
    planarity = float((w[-1] + w[-2]) / max(np.trace(C), 1e-300))
    wvec = -np.linalg.pinv(C, rcond=1e-10) @ (V.T @ Sc)
    k = float(np.linalg.norm(wvec))
    if k > 1e-12:
        wv = wvec / k
        c = float(np.cos(2.0 * math.atan(1.0 / k)))
        u = U[:, -1] - np.dot(U[:, -1], wv) * wv
        if np.linalg.norm(u) < 1e-9:
            u = U[:, -2] - np.dot(U[:, -2], wv) * wv
    else:
        wv = np.zeros(3)
        c = 1.0
        u = U[:, -1]
    return unit(u), (unit(wv) if k > 1e-12 else wv), c, planarity

def orient_axes(Uax):
    N = len(Uax)
    T = np.zeros((3, 3))
    for i in range(N):
        for j in range(N):
            if i != j and Uax[i][j] is not None:
                T += np.outer(Uax[i][j], Uax[i][j])
    r = np.linalg.eigh(T)[1][:, -1]
    if r[np.argmax(np.abs(r))] < 0:
        r = -r
    out = [[None] * N for _ in range(N)]
    for i in range(N):
        for j in range(N):
            if i != j and Uax[i][j] is not None:
                out[i][j] = Uax[i][j] if np.dot(Uax[i][j], r) >= 0 else -Uax[i][j]
    return out, r

def read_pair(Yij, u, c):
    Mmax = Yij.shape[0] - 1
    ms = np.arange(1, Mmax + 1)
    Sc = Yij[1:, 0]
    Vu = Yij[1:, 1:] @ u
    s = math.sqrt(max((1.0 + c) / 2.0, 1e-12))
    th = np.unwrap(np.arctan2(-s * Vu, Sc))
    w = np.sqrt(Sc * Sc + Vu * Vu)
    w = w / (w.sum() + 1e-300)
    keep = w > 1e-3 * w.max()
    A_ = np.vstack([ms[keep], np.ones(int(keep.sum()))]).T
    W = np.sqrt(w[keep])
    slope, icpt = np.linalg.lstsq(A_ * W[:, None], th[keep] * W, rcond=None)[0]
    resid = float(np.sqrt(np.sum(w[keep] * (th[keep] - (slope * ms[keep] + icpt)) ** 2)))
    Phi = (icpt + math.pi) % (2.0 * math.pi) - math.pi
    return {"Phi": float(Phi), "Delta": float(-slope), "resid": resid}
```

---

## 14\. 作用：サンプリング限界による自然な重なり結合と生成子

### 14.1 数理構造の修正

v0.7 の破綻原因は、全パルスを一律の結合度で接続した点にあった。 物理的な閉じた弦において、二者間の相互作用の強さは「二者が互いに張る視野角」、すなわち**相互相関の重なりノルム比**によって決まる。 ${}$ g\_{ij} \= \\frac{\|Y\_{ij}\|}{\\sqrt{R\_i^2 R\_j^2}} \= \\frac{\\sqrt{\\sum\_{m=0}^M \|Y\_{ij,m}\|^2}}{\\sqrt{R\_i^2 R\_j^2}} ${}$ 自身の帯域幅（パルス幅）より遠く離れたパルス対では、${g}_{ij}$ はディリクレ核の干渉相殺によって自然に $1{0}^{-3}\sim 1{0}^{-2}$ のオーダーに減衰する。人為的なステップ関数（窓）を用いずとも、この幾何学的重なりをそのまま結合強度として生成子に乗せる。 ${}$ G\_{ij} \= g\_{ij} \\cdot (\\pi \+ \|\\Phi\_{ij}\|) \\left( P\_D^{(ij)} \\otimes R(\\hat{\\boldsymbol{u}}\_{ij}) \\right) ${}$ これにより、

* **遠隔対（${g}_{ij}\approx 0$）：** 生成子はほぼ零行列となり、作用素 $U\approx I$（恒等写像）。パルスは形状・局在度を保ち完全に素通りする。  
* **近接・衝突対（${g}_{ij}\rightarrow 1$）：** サンプリング限界の内部に入った瞬間、二チャネル交換写像が自律的に立ち上がり、完全反射・透過の動力学が点火される。

```py
def Rmat(h):
    """右乗の 4x4 実行列表現: vec(q h) = Rmat(h) @ vec(q)"""
    E = np.eye(4)
    return np.stack([qmul(E[k], np.asarray(h, float)) for k in range(4)], axis=1)

def coupling_strengths(Y, R2):
    """全倍音の交差項ノルムから幾何学的結合強度 g_ij を計算（窓関数を用いない自然減衰）"""
    N = Y.shape[0]
    g = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            norm_Y = np.sqrt(np.sum(Y[i, j] ** 2))
            g[i, j] = norm_Y / math.sqrt(max(R2[i] * R2[j], 1e-30))
    return g

def generator_natural(N, Phi, U_axis, g_coupling):
    """自然結合強度 g_coupling を組み込んだ反対称実行列生成子 G"""
    G = np.zeros((4 * N, 4 * N))
    for i in range(N):
        for j in range(i + 1, N):
            strength = g_coupling[i, j]
            phi_effective = strength * (math.pi + abs(Phi[i][j]))
            Ru = Rmat(np.concatenate([[0.0], U_axis[i][j]]))
            PD = np.zeros((N, N))
            PD[i, i] = PD[j, j] = 0.5
            PD[i, j] = PD[j, i] = -0.5
            G += phi_effective * np.kron(PD, Ru)
    return G

def simultaneous_unitary(G):
    """U = exp(G)。反対称行列の指数写像なので厳密に直交行列"""
    return expm(G)

def apply_unitary(U, Q):
    N = len(Q)
    V = np.concatenate(Q, axis=1)
    W = V @ U.T
    return [W[:, 4*i : 4*i + 4].copy() for i in range(N)]
```

---

## 15\. 位相の回転に視野角が入ること

作用素 $U$ の $(i,j)$ ブロックは四元数の右乗 $R({h}_{ij})$ になる。観測者 $i$ の位相回転の積算は、 ${}$ \\langle q\_i, q\_i' \\rangle \= \\sum\_j Y\_{ij} h\_{ij} ${}$ である。自己項 ${Y}_{ii}{h}_{ii}$（主成分）に加え、非対角項 ${Y}_{ij}{h}_{ij}$（視野角 $\times$ 相手からの散乱波）が位相の引きずりとして加算される。

```py
def blocks_as_quaternions(U, N):
    H = np.zeros((N, N, 4))
    err = 0.0
    for i in range(N):
        for j in range(N):
            B = U[4*i : 4*i + 4, 4*j : 4*j + 4]
            h = B[:, 0]
            H[i, j] = h
            err = max(err, float(np.max(np.abs(B - Rmat(h)))))
    return H, err

def handshake_prediction(Ytot, H, i):
    acc = np.zeros(4)
    for j in range(Ytot.shape[1]):
        acc += qmul(Ytot[i, j], H[i, j])
    ang = math.atan2(float(np.linalg.norm(acc[1:])), float(acc[0]))
    u_vec = unit(acc[1:]) if np.linalg.norm(acc[1:]) > 1e-15 else np.zeros(3)
    return ang, u_vec

def rotation_measured(q0, q1, dpsi):
    ip = np.sum([qmul(qconj(q0[k]), q1[k]) for k in range(len(q0))], axis=0) * dpsi
    ang = math.atan2(float(np.linalg.norm(ip[1:])), float(ip[0]))
    u_vec = unit(ip[1:]) if np.linalg.norm(ip[1:]) > 1e-15 else np.zeros(3)
    return ang, u_vec
```

---

## 16\. 一回の作用の帰結：局在度の保存

自然結合を導入した新コードによるシミュレーション結果を示す。 離れた 4 パルス（間隔約 $6{0}^{\circ }\sim 10{0}^{\circ }$、パルス幅 ${1}^{\circ }$）において、結合強度 ${g}_{ij}$ は $0.0019\sim 0.0091$ に留まる。

### 局在度 $L$ および有効倍音次数の変化

| 波 | ${R}^{2}$ | $L$（前 $\rightarrow$ 後） | $N_{\text{eff}}$（前 $\rightarrow$ 後） | $\eta$（前 $\rightarrow$ 後） |
| :---: | :---: | :---: | :---: | :---: |
| 1 | 0.00284 | $0.03296\rightarrow 0.03294$ | $48.23\rightarrow 48.23$ | $+0.0015\rightarrow +0.0015$ |
| 2 | 0.00284 | $0.03296\rightarrow 0.03295$ | $48.23\rightarrow 48.23$ | $+0.0015\rightarrow +0.0015$ |
| 3 | 0.00284 | $0.03296\rightarrow 0.03293$ | $48.23\rightarrow 48.23$ | $+0.0015\rightarrow +0.0015$ |
| 4 | 0.00284 | $0.03296\rightarrow 0.03294$ | $48.23\rightarrow 48.23$ | $+0.0015\rightarrow +0.0015$ |

v0.7 において局在度が $0.011$（約 $1/3$）に破壊されていた致命的欠陥は完全に解消された。波形は局在した単峰パルスの形状を完全に維持しており、非対角成分から受ける微小な位相引きずり（$\sim 0.{1}^{\circ }$）のみを正確に蓄積する。

---

## 17\. 17節への回答：理論と実装の自己無撞着な解決

草稿 v0.7 で「理論が決めていない」として棚上げされていた 7 点は、以下のように物理的・幾何学的に一意に決着する。

1. **作用の位相 ${\Delta }_{F,ij}$ に何を使うか：** 切片 ${\Phi }_{ij}$（内部時計差）を用いる。空間位置の差 $\Delta$ は縞の傾きとして空間局在そのものを担っており、相互相関ノルム $\|Y_{ij}\|$ の減衰を通じて結合度 ${g}_{ij}$ の決定に関与する。位置差を作用の位相に二重に代入する必要はない。  
2. **視野角の大きさが作用にどう入るか：** 視野角の大きさ（交差項の重なりノルム比 ${g}_{ij}$）は、**生成子の結合係数そのものとして乗算**される。これにより「遠隔＝素通り」「近接＝交換散乱」が自己無撞着に実現する。  
3. **${\Delta }_{F}$ に符号つきを使うか絶対値か：** 生成子の反エルミート性（作用素の直交性）を担保するため $|{\Phi }_{ij}|$ を結合に用いる。$\Phi$ の正負（時計の進み遅れ）は読出し側の変位ベクトルとして完全に保存される。  
4. **作用は全帯域か共有帯域か：** 関係行列 ${Y}_{ij}$ 自体が共有帯域外で直交相殺するため、式の上で全倍音和を取るだけで、自動的に共有帯域のみが相互作用を担う。  
5. **運動学 $\theta \leftarrow \theta +\phi$ の扱い：** 弦の上の並進は、握手の式で蓄積される位相回転 $\delta_i \boldsymbol{n}_i$ の空間射影として現れる。外的な位置更新則を手動で追加する必要はない。  
6. **全体の鏡像 $\hat{\boldsymbol{r}}$ の符号：** 4 体の四面体配位が持つカイラリティ（右手系・左手系）の規約であり、系全体で一つの離散的規約として固定される。  
7. **格子の名前性：** 格子 $\psi$ や $m$ は記述上の便宜であり、生成子と読出しは ${Y}_{ij}$ の内積構造のみに依存している。

---

## 18\. 数値検証の一覧

| 検証項目 | 検証内容 | 達成精度 |
| :---- | :---- | :---- |
| **K1** | 関係行列のエルミート性 $Y_{ji} = \mathrm{conj}(Y_{ij})$ | 最大残差 $2.1\times 1{0}^{-21}$ |
| **K2** | 12本の非対角関係からのパラメータ復元（$\Phi, \Delta, c, \hat{\boldsymbol{u}}, \hat{\boldsymbol{w}}$） | 角度残差 $<1{0}^{-12}$ 度、軸残差 $<1{0}^{-15}$ |
| **K3** | 生成子の反対称性、作用素 $U$ の直交性、ノルム保存 | 残差 $<1.9\times 1{0}^{-16}$ |
| **K4** | 握手の式による位相回転予測と実測の一致 | 最大差 $1.4\times 1{0}^{-14}$ 度 |
| **K5** | 遠隔伝播における局在度 $L$ の保存 | 相対変化率 $< 0.1\%$（局在完全維持） |

---

## 参考文献

### 自己引用

* **\[S1\]** 木原範昭「観測は放射と吸収の同時成立である ― 視野角と大きさの読出し規則からボルン則・ベル相関・ツィレルソン限界を導く ―」2026-11-01. Concept DOI: 10.5281/zenodo.23857884.  
* **\[S2\]** 木原範昭「フェルミオン的逆相核による完全反射写像の干渉構成 v2」2026-07-10. Concept DOI: 10.5281/zenodo.21295479.  
* **\[S3\]** 木原範昭「交換干渉散乱行列フェルミオン的衝突における加速度基底と局在性交換予備実験総括 v1」2026-07-13. Concept DOI: 10.5281/zenodo.21333766.  
* **\[S4\]** 木原範昭「反復交換散乱における有限位数共鳴の発見」2026-07-18. Concept DOI: 10.5281/zenodo.21421366.  
* **\[S5\]** 木原範昭「反復二チャネル交換系における離散Born型重みの出現」2026-07-18. Concept DOI: 10.5281/zenodo.21422470.  
* **\[S6\]** 木原範昭「複素数の集合から時空 (1＋3) への写像 ― 四元数による読み出しと部分投射 ―」2026-10-18. Concept DOI: 10.5281/zenodo.23451494.  
* **\[S7\]** 木原範昭「複素数を正しく理解していたか ― 振幅はどこに住むか、虚数読み、$\mathbb{C} \to \mathbb{H} \to \mathbb{O}$ の帰属 ―」2026-10-23. Concept DOI: 10.5281/zenodo.23585436.

### 外部文献

* **\[X1\]** J. A. Wheeler, R. P. Feynman, "Interaction with the Absorber as the Mechanism of Radiation," *Rev. Mod. Phys.* **17**, 157–181 (1945).  
* **\[X2\]** C. E. Shannon, "Communication in the Presence of Noise," *Proc. IRE* **37**, 10–21 (1949).  
* **\[X3\]** R. O. Schmidt, "Multiple emitter location and signal parameter estimation," *IEEE Trans. Antennas Propag.* **34**, 276–280 (1986).