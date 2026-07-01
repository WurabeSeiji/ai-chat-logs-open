# Claude Code 修正指示：標準ハミルトニアンの空間微分エネルギーを保存する奇数倍音波形への変更

## 目的

既存の奇数倍音合成波形について、単純に同一振幅で奇数倍音を加算するのではなく、**標準自由粒子ハミルトニアンの空間微分エネルギー**が、基準波形 `N=1` の場合と等しく保存されるように波形定義を修正してください。

ここで対象とする 1 次元標準自由粒子ハミルトニアンは、定数係数を除けば、空間微分項

```math
\hat H = -\frac{\hbar^2}{2m}\frac{d^2}{dx^2}
```

に対応します。したがって、波形比較で保存すべきエネルギー量を次で定義します。

```math
\mathcal{E}[\psi]
=
\int_0^{2\pi}\left|\frac{d\psi}{dx}\right|^2 dx
```

基準波形を

```math
\psi_1(x)=\sin x
```

とすると、

```math
\mathcal{E}[\psi_1]
=
\int_0^{2\pi}\cos^2 x\,dx
=
\pi
```

です。以後、奇数倍音合成波形でも必ず

```math
\mathcal{E}[\psi_N]=\pi
```

となるように正規化してください。

---

## 修正後の波形定義

最高奇数倍音を

```math
N=2K-1
```

とします。ここで `K` は含める奇数倍音の項数です。

奇数倍音合成波形を次のように定義してください。

```math
\psi_N(x)
=
\sqrt{\frac{3}{K(4K^2-1)}}
\sum_{m=0}^{K-1}\sin((2m+1)x)
```

この係数

```math
a_K
=
\sqrt{\frac{3}{K(4K^2-1)}}
```

は、標準ハミルトニアンの空間微分エネルギーを `N=1` の基準波形と一致させるための正規化係数です。

---

## 閉形式

次の恒等式を使っても構いません。

```math
\sum_{m=0}^{K-1}\sin((2m+1)x)
=
\frac{\sin^2(Kx)}{\sin x}
```

したがって、閉形式では次のように実装できます。

```math
\psi_N(x)
=
\sqrt{\frac{3}{K(4K^2-1)}}
\frac{\sin^2(Kx)}{\sin x}
```

ただし、`sin(x)=0` となる点では数値的な 0 除算が発生するため、閉形式を使う場合は極限値で処理してください。安全性を優先するなら、和形式で実装してください。

---

## エネルギー保存条件の確認式

実装後、各 `K` について次を数値検証してください。

```math
\int_0^{2\pi}\left|\frac{d\psi_N}{dx}\right|^2 dx
\approx
\pi
```

解析的には、

```math
\frac{d\psi_N}{dx}
=
\sqrt{\frac{3}{K(4K^2-1)}}
\sum_{m=0}^{K-1}(2m+1)\cos((2m+1)x)
```

であり、直交性より、

```math
\int_0^{2\pi}\left|\frac{d\psi_N}{dx}\right|^2 dx
=
\pi\cdot
\frac{3}{K(4K^2-1)}
\sum_{m=0}^{K-1}(2m+1)^2
```

奇数平方和は、

```math
\sum_{m=0}^{K-1}(2m+1)^2
=
\frac{K(4K^2-1)}{3}
```

なので、

```math
\int_0^{2\pi}\left|\frac{d\psi_N}{dx}\right|^2 dx
=
\pi
```

となります。

---

## 実装上の注意

- 既存コードで奇数倍音を単純に `sum(sin((2*m+1)*x))` としている箇所を、必ず係数

```math
\sqrt{\frac{3}{K(4K^2-1)}}
```

でスケールしてください。

- `N` は最高奇数倍音、`K` は項数です。関係は次の通りです。

```math
K=\frac{N+1}{2}
```

- `N` は奇数のみを想定してください。偶数が渡された場合は、エラーにするか、最も近い奇数に丸めるのではなく、明示的に入力不正として扱ってください。

- 図示する場合は、旧定義の同一振幅合成波形ではなく、このエネルギー保存正規化済み波形を使用してください。

---

## 実装例：Python 形式の疑似コード

```python
import numpy as np

def psi_energy_preserving_odd_harmonics(x: np.ndarray, N: int) -> np.ndarray:
    """
    Energy-preserving odd-harmonic waveform.

    N: highest odd harmonic. Must be odd and >= 1.
    K: number of odd harmonic terms, K = (N + 1) // 2.
    The normalization preserves int_0^{2pi} |dpsi/dx|^2 dx = pi,
    equal to the baseline waveform sin(x).
    """
    if N < 1 or N % 2 == 0:
        raise ValueError("N must be a positive odd integer.")

    K = (N + 1) // 2
    scale = np.sqrt(3.0 / (K * (4 * K**2 - 1)))

    y = np.zeros_like(x, dtype=float)
    for m in range(K):
        n = 2 * m + 1
        y += np.sin(n * x)

    return scale * y
```

---

## 検証対象

最低限、以下の `N` で検証してください。

```text
N = 1, 3, 5, 7, 9, 13
```

期待される性質は、すべての `N` について

```math
\mathcal{E}[\psi_N]
=
\int_0^{2\pi}|\psi_N'(x)|^2dx
=
\pi
```

であることです。

---

## この修正の意味

この修正は、奇数倍音を同一振幅で無制限に加算するものではありません。標準ハミルトニアンの空間微分エネルギーが `N=1` の基準波形と同じになるように、奇数倍音合成波形を正規化するものです。

したがって、奇数倍音を増やして波形が局在化しても、定義上、標準ハミルトニアンの空間微分エネルギーは基準波形と同じ値に保存されます。
