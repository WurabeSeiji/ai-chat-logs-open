# E3-A1 階層ケプラー摂動近似の崩壊点テスト v3

## 問い

E3-A0-v2 の相互線形結合は正常モードへ厳密対角化できる。そこで新しい効果を作らず、「無摂動の内部/外部二体系として扱う摂動近似」がどこまで有効かだけを調べる。

$$H=\begin{pmatrix}\tau_{in}&\varepsilon\\\varepsilon&\tau_{out}\end{pmatrix},\quad\lambda_\pm=\frac{\tau_{in}+\tau_{out}}2\pm\sqrt{(\Delta\tau/2)^2+\varepsilon^2}.$$

$n_{in}=31$, $n_{out}=127$, $\tau_{in}=1.959059882505$, $\tau_{out}=1.997552832429$, $\Delta\tau=0.038492949924$.

## 結果

- $O(\varepsilon^2)$ の外部法則シフトの相対誤差が 1%: $\varepsilon=0.00386849$。
- 同 5%: $\varepsilon=0.00881984$。
- 同 10%: $\varepsilon=0.01276667$。
- 外部優勢正常モードが楕円型から放物境界 $\lambda_+=2$ に達する厳密値: $\varepsilon_c=0.01000936$。
  この点でも混合角は 13.739 deg、内部成分重量は 0.0564、二次近似誤差は 6.3574%。
- $2\varepsilon/\Delta\tau=1$ の自然な混合クロスオーバー: $\varepsilon=\Delta\tau/2=0.01924647$。  ここで混合角は 22.5 deg、外部モード中の内部重量は 14.64%。

## 解釈

この候補では exact system 自体は崩れない。崩れるものが三段階に分かれる。

1. まず $\Delta\tau_{out}\simeq\varepsilon^2/\Delta\tau$ という摂動展開の精度が徐々に落ちる。
2. それより明確な境界として、$\varepsilon_c$ で外部正常モードが楕円型を離れ、閉 Kepler 軌道という解釈が失われる。
3. さらに $\varepsilon\sim\Delta\tau$ では内部/外部の固有ベクトルが強く混ざり、(ab)|c という階層ラベル自体が弱くなる。

特に今回のパラメータでは、強混合より先に楕円型境界へ達する。したがって最初の物理的な「崩壊点」は $\varepsilon_c\simeq 0.01001$ と読むのが最も非恣意的である。
