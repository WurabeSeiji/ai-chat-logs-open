---
title: "奇数倍音孤立ピーク波の二コピー重ね合わせ：相対位相の波形不変性とコントラスト則"
emoji: "🌀"
type: "idea"
topics: ["数学", "フーリエ解析", "物理学", "信号処理", "観察論文"]
published: true
---

新しい観察・整理論文を公開しました。半波長区間上の奇数倍音孤立ピーク波（前稿）を、**二つのコピーに相対位相を与えて重ね合わせたとき**に何が起こるかを観察したものです。

- **Concept DOI**: https://doi.org/10.5281/zenodo.20923461
- **Version DOI (v0.1)**: https://doi.org/10.5281/zenodo.20923462
- **GitHub**: https://github.com/WurabeSeiji/ai-chat-logs-open/tree/main/平方数を基本量とした場合の検討

## 出発点：孤立ピーク波

半波長の位相区間 $\varphi\in[-\pi/2,\pi/2]$ 上で、一定振幅の奇数倍音だけを余弦で重ね合わせた波

$$
S_N(\varphi)=\sum_{m=0}^{(N-1)/2}\cos\!\big((2m+1)\varphi\big)
$$

は、中央に主ピーク $S_N(0)=(N+1)/2$ をもち両端 $\varphi=\pm\pi/2$ で零となる「孤立ピーク波」です（$N$ は最高奇数倍音次数、奇数）。両端の零は、各奇数倍音が $\cos((2m+1)(\pm\pi/2))=0$ を満たすことから、振幅や項数によらず成り立ちます。

## 二つのコピーに共通の相対位相を与える

この波の二つのコピーに、**全倍音へ共通**の相対位相 $\pm\alpha$ を与えて重ね合わせます。

$$
\psi_\alpha(\varphi)=\sum_{m=0}^{(N-1)/2}\Big[\cos\!\big((2m+1)\varphi-\alpha\big)+\cos\!\big((2m+1)\varphi+\alpha\big)\Big]
$$

ここで $\alpha$ は倍音次数 $n=2m+1$ に依存しない一定量です（倍音次数に比例してずれる量＝区間内の平行移動とは別物）。観察量は二乗振幅 $I_\alpha=|\psi_\alpha|^2$ とします。

## 主結果：波積の恒等式一本で決まる

和積の三角恒等式 $\cos(n\varphi-\alpha)+\cos(n\varphi+\alpha)=2\cos(n\varphi)\cos\alpha$ を各倍音に適用すると、$\cos\alpha$ は $n$ に依存しない共通因子なので和の外に出て、

$$
\psi_\alpha(\varphi)=2\cos\alpha\cdot S_N(\varphi),\qquad
I_\alpha(\varphi)=4\cos^2\!\alpha\cdot I_N(\varphi)
$$

が**厳密に**成り立ちます。ここから二点が直ちに従います。

### (1) 波形不変性

$\varphi$ への依存はすべて $I_N(\varphi)$ に含まれ、相対位相 $\alpha$ は前因子 $4\cos^2\!\alpha$ にしか現れません。したがってピーク値で規格化した波形は

$$
\widehat{I}_\alpha(\varphi)=\frac{I_\alpha(\varphi)}{I_\alpha(0)}=\frac{I_N(\varphi)}{I_N(0)}=\widehat{I}_N(\varphi)
$$

となり、**相対位相に依存せず、単一コピーの波形と完全に一致**します。形・両端の零・局在幅は、相対位相を変えても不変です。

### (2) コントラスト則

規格化前のピーク二乗振幅は

$$
I_\alpha(0)=(N+1)^2\cos^2\!\alpha
$$

で、同相 $\alpha=0$ の値 $(N+1)^2$ を基準にすると比は $\cos^2\!\alpha$。つまり相対位相 $2\alpha$ は、**波形を変えずに二乗振幅を $\cos^2\!\alpha$ 倍**します。$\alpha=0$ で最大、$\alpha=\pi/2$（相対位相 $\pi$）で零です。

具体例として $N=9$（$n=1,3,5,7,9$）、$\alpha=15^\circ$（相対位相 $30^\circ$）では、$I_0(0)=100$、$I_\alpha(0)=100\cos^2 15^\circ\approx 93.30$、コントラスト因子 $\cos^2 15^\circ\approx 0.9330$ です。

## 本稿の位置づけ

これらはいずれも、和積の三角恒等式を一定振幅奇数倍音和に適用しただけの初等的な事実です。新発見として主張するものではなく、物理的解釈も与えません。相対位相が「波形を変える変形」ではなく「振幅にのみ効くスカラー」として現れる、という分離を、誰でも辿れる導出として一つに整理した観察論文です。
