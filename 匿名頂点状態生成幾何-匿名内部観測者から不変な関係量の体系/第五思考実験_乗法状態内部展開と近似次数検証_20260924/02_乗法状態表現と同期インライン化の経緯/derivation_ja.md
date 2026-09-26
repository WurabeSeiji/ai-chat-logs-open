# 導出式: 乗法状態表現へ至る経緯

候補となった永続状態は

$$
Z=(U,P,E,H,Q)=(e^{i\chi},p,e,e^{i\phi},e^{t/\tau_0}).
$$

位相は

$$
\cos\chi=\frac{U+U^{-1}}{2},\qquad
\sin\chi=\frac{U-U^{-1}}{2i}
$$

と `arg` を使わずに内部評価できる。

時計は

$$
Q_{n+1}=Q_n\exp(\Delta t/\tau_0)
$$

と更新し、観測時だけ

$$
t=\tau_0\log Q
$$

と読む。

正の尺度 $P=p$, $E=e$ は、文字通り $e^p,e^e$ と置くのではなく、正の乗法群要素そのものとして保持する。文字通り $P=e^p,E=e^e$ とすると現行RR式の $p^{3/2},e^k$ を評価するために内部 `log` が必要になる。

`audit_five_state_full_inline.py` は、RK4中間値を現在状態だけの式にインライン化した段階を検証する。
