# E3-A3 解析的分析 v1 — 摂動予備実験の固定点から何が確定したか

日付: 2026-09-21  
対象: E3-A0 / E3-A1 / E3-A2 の固定済みデータのみ。  
位置づけ: **追加数値実験ではなく、固定済み候補の解析的整理**。

## 1. 候補則

固定した候補は

$$
X^{in}_{k+1}+X^{in}_{k-1}=\tau_{in}X^{in}_k+\varepsilon X^{out}_k,
$$

$$
X^{out}_{k+1}+X^{out}_{k-1}=\tau_{out}X^{out}_k+\varepsilon X^{in}_k.
$$

主系列では

$$
\tau_{in}=1.959059882505,\qquad
\tau_{out}=1.997552832429,\qquad
\Delta\tau=0.038492949924.
$$

## 2. 全体面積保存は数値結果ではなく厳密恒等式

$$q_{in,k}=\omega(X^{in}_k,X^{in}_{k+1}),\qquad
q_{out,k}=\omega(X^{out}_k,X^{out}_{k+1}).$$

更新則をそのまま代入すると

$$
q_{in,k}-q_{in,k-1}=\varepsilon\,\omega(X^{in}_k,X^{out}_k),
$$

$$
q_{out,k}-q_{out,k-1}=-\varepsilon\,\omega(X^{in}_k,X^{out}_k).
$$

従って

$$
\boxed{q_{tot,k}=q_{in,k}+q_{out,k}=\mathrm{const}}
$$

は厳密である。代表固定点 $\varepsilon=0.004$ の再生では、全体面積の相対ドリフトは
`9.626e-15`、上の交換恒等式の最大絶対残差は inner / outer とも
`4.450e-16` の数値丸め域であった。

**したがって局所面積の変動は保存則の破れではなく、二セクター間の厳密な交換である。**

## 3. なぜ局所応答が一次になるか

非縮退なら階層空間の行列

$$
H=\begin{pmatrix}\tau_{in}&\varepsilon\\\varepsilon&\tau_{out}\end{pmatrix}
$$

を角度 $\theta$ で直交対角化でき、

$$
\tan 2\theta=\frac{2\varepsilon}{\Delta\tau}.
$$

小結合では

$$
\theta=\frac{\varepsilon}{\Delta\tau}+O\!\left((\varepsilon/\Delta\tau)^3\right).
$$

生の外部状態は二つの正常モードの混合なので、局所面積には二モード間の交差面積が係数
$\sin\theta\cos\theta=O(\varepsilon/\Delta\tau)$ で入る。従って固定データで

- inner local-area modulation: slope `1.015987`
- outer local-area modulation: slope `0.985208`

とほぼ一次になったことは、厳密対角化から説明できる。

## 4. 生の軌道が単一円錐から一次で外れる理由

複素表示で正常モードを $y_-,y_+$、外部状態を

$$w_{out}=\sin\theta\,y_-+\cos\theta\,y_+$$

とすれば、論文2と同じ二乗読み出しは

$$
z_{out}=w_{out}^2
=\sin^2\theta\,y_-^2+\cos^2\theta\,y_+^2
+2\sin\theta\cos\theta\,y_-y_+.
$$

最後の交差項が $O(\varepsilon/\Delta\tau)$ である。従って **生の外部読み出しは一つの Kepler 円錐ではなく、二つの正常モードと交差項の合成**になる。
固定データの raw conic residual の小結合 slope は `0.990863` であり、この一次則と一致する。

一方、正常モード自身は

$$Y_{\pm,k+1}+Y_{\pm,k-1}=\lambda_\pm Y_{\pm,k}$$

を厳密に満たすため、各モードの面積保存と円錐残差は数値精度域に留まる。

## 5. 外部法則シフトが二次になる理由

$$
\lambda_+=\frac{\tau_{in}+\tau_{out}}2+
\sqrt{(\Delta\tau/2)^2+\varepsilon^2},
$$

したがって

$$
\delta\tau_{out}=\lambda_+-\tau_{out}
=\frac{\varepsilon^2}{\Delta\tau}-
\frac{\varepsilon^4}{\Delta\tau^3}+O(\varepsilon^6).
$$

固定データの fitted slope は `1.999615`、prefactor は
`25.876871` であり、解析値 $1/\Delta\tau=25.978783$ と一致する。

ここで本当の展開変数は

$$
\boxed{g=\varepsilon/|\Delta\tau|}
$$

であり、$\varepsilon$ 単独ではない。

## 6. 摂動近似の精度は g だけで決まる

二次近似 $\varepsilon^2/\Delta\tau$ の exact shift に対する相対誤差は

$$
R(g)=\sqrt{\frac14+g^2}-\frac12.
$$

従って相対誤差 $R=r$ の境界は

$$
\boxed{g=\sqrt{r(1+r)}}
$$

である。主系列に戻すと

- 1%: $g=0.100498756$, $\varepsilon=0.003868494$
- 5%: $g=0.229128785$, $\varepsilon=0.008819843$
- 10%: $g=0.331662479$, $\varepsilon=0.012766667$

となり、E3-A1 の数値表と丸め誤差内で一致する。

## 7. 混合も g の普遍関数

$$
\theta(g)=\frac12\arctan(2g),
$$

$$
W_{in\to out}(g)=\sin^2\theta
=\frac12\left(1-\frac1{\sqrt{1+4g^2}}\right).
$$

E3-A1 の全 scan に対する解析式との差の最大値は、shift `2.220e-16`、
relative error `1.336e-11`、angle `0.000e+00` deg、weight `0.000e+00` である。

**従って「摂動近似が厳しい」という数値的観察は、非縮退条件 $g\ll1$ という解析条件に置き換えられる。**

## 8. 縮退系は摂動階層に対して特異

$\Delta\tau=0$ なら任意の $\varepsilon\ne0$ で

$$
\theta=\pi/4,
$$

すなわち正常モードは最初から 50:50 の対称・反対称結合になる。固定データでも非零 $\varepsilon$ の mixing angle はすべて
`45.0` deg である。

これは exact system の破綻ではない。**「inner / outer を弱く混ざる二つの独立対象として追跡できる」という階層解釈だけが、縮退点で特異になる。**

## 9. 閉軌道境界は別の条件

外部優勢正常モードが楕円型から放物境界へ達する条件は $\lambda_+=2$ なので、

$$
\boxed{\varepsilon_c^2=(2-\tau_{in})(2-\tau_{out})}.
$$

主系列では

$$
\varepsilon_c=0.010009362,\qquad g_c=0.260031045.
$$

この点の mixing angle は `13.738616` deg、inner weight は `0.056403`、
二次近似誤差は `0.063574` である。

これは「摂動展開が使えるか」とは別に、正常モード自身が閉 Kepler 軌道でいられるかを決めるスペクトル境界である。

## 10. 位相依存性の評価

8初期位相の固定データでは、小結合域 $\varepsilon\le0.002$ における raw conic residual の
最大/最小比は最大でも `1.341374`。従って一次スケーリングは特定の初期位相だけの結果ではない。
一方 $\varepsilon=0.004$ 以降は位相差の影響幅が急増し、raw conic residual は単調な一変数関数ではなくなる。
これは摂動域を外れ始めたという解釈と整合する。

## 11. この候補についての固定結論

この候補については、次が解析的に確定した。

1. **全体面積保存は厳密**で、局所面積は等量反対向きに交換する。
2. 非縮退小結合では、局所量・生の円錐残差は $O(\varepsilon/\Delta\tau)$。
3. 正常モードの法則シフトは $O((\varepsilon/\Delta\tau)^2)$。
4. 系全体は固定直交変換で二つの二値則へ厳密分解できる。
5. 従ってこの候補は **三状態既約力学ではない**。
6. 今後の真の三状態候補に対しては、**状態に依存しない固定基底変換で二体セクターへ完全分解できたら落とす**、という新しい判定条件を置ける。

最後の条件が、この摂動予備実験から本体へ持ち帰る最も重要な結果である。

## 12. 追加図

- `fig_E3A_08_local_area_exchange_timeseries_v1.png`: 局所面積の等量交換と全体保存。
- `fig_E3A_09_phase_robustness_v1.png`: 8位相に対する一次応答の頑健性と高結合側の位相幅拡大。
- `fig_E3A_10_universal_g_control_v1.png`: $g=\varepsilon/|\Delta\tau|$ による摂動誤差と混合重量の解析曲線。
- `fig_E3A_11_degenerate_vs_nondegenerate_mixing_v1.png`: 非縮退と縮退での mixing angle の質的差。
