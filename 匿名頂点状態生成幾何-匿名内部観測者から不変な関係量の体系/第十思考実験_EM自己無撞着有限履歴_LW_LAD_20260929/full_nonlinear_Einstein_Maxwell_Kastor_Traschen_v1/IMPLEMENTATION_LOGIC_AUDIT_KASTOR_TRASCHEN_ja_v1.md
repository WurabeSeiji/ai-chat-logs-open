# Kastor–Traschen 二中心 full nonlinear Einstein–Maxwell–Λ 実装
## 論文用：初期化・実行系・読出し系の主要コードと数式対応

対象プログラム：`kastor_traschen_two_center_v1.py`  
SHA-256：`893409339fc3f3509b3dbb1125d3a40bf96074ad7bdbbad7015639a6d734f19a`

---

## 0. この付録の目的

この付録の目的は、数値結果そのものではなく、**プログラムロジックが論文中で主張する物理モデルを実際に実装しているかを第三者がコードと数式を対応させて検証できるようにすること**である。

本プログラムが行っている処理は、次の3層に分けられる。

1. **初期化**：二中心 Kastor–Traschen 厳密解を定義する。
2. **実行系**：各時刻で厳密解から 3+1 幾何量と Maxwell 場を再構成し、空間微分のみを有限差分評価する。
3. **読出し系**：Hamiltonian/Gauss 制約、背景物理距離、`aa, ab, bb` 非線形項を読み出して JSON に保存する。

最重要の区別は次である。

> **このコードは Einstein–Maxwell 方程式を数値時間積分して解を生成する PDE solver ではない。Kastor–Traschen の full nonlinear exact solution を各時刻に解析式から再構成し、その場が制約式を満たすことを有限差分で独立検証する evaluator である。**

したがって、論文中では「full nonlinear Einstein–Maxwell–Λ の動的二中心厳密解を実装・評価・制約検証した」と記述できるが、「一般の charged binary を数値相対論で時間発展させた」とは記述しない。

---

# 1. 物理モデルとコードの中心変数

Kastor–Traschen 解を

\[
 ds^2=-\Omega^{-2}dt^2+a(t)^2\Omega^2\delta_{ij}dx^i dx^j,
\]

\[
 \Omega(t,\mathbf x)=1+\sum_A\frac{M_A}{a(t)r_A},
 \qquad
 a(t)=e^{Ht},
 \qquad
 \Lambda=3H^2
\]

とする。

コードでは

\[
 W\equiv a\Omega
 =a(t)+\sum_A\frac{M_A}{r_A}
\]

を基本量として使う。

この置換により 3+1 量は

\[
\gamma_{ij}=W^2\delta_{ij},
\qquad
\alpha=\frac{a}{W},
\qquad
\beta^i=0,
\]

\[
K_{ij}=-H\gamma_{ij},
\qquad
K=-3H
\]

となる。

Gauge potential を

\[
A_t=+\Omega^{-1}
\]

と取る符号規約では、電場は

\[
E_i=-\partial_i\ln W,
\qquad
E^i=-\frac{\partial_iW}{W^3}
\]

となる。

---

# 2. 初期化ロジック

## 2.1 格子・二中心・静的部分の構築

実装主要部：

```python
def run(N=49, L=12.0, H=-0.05, sep=8.0,
        m1=1.0, m2=1.0,
        times=(0.0,4.0,8.0), exclusion=1.2):
    x=np.linspace(-L,L,N); h=x[1]-x[0]
    X,Y,Z=np.meshgrid(x,x,x,indexing='ij')
    centers=[(-sep/2,0,0),(sep/2,0,0)]
    rs=[]
    for cx,cy,cz in centers:
        rs.append(np.sqrt((X-cx)**2+(Y-cy)**2+(Z-cz)**2))

    eps=1e-15
    Ustatic=m1/np.maximum(rs[0],eps)+m2/np.maximum(rs[1],eps)
    results=[]
    Lam=3*H*H
```

### 対応する数式

二中心位置を

\[
\mathbf x_a=(-d/2,0,0),
\qquad
\mathbf x_b=(+d/2,0,0)
\]

とし、

\[
r_a=|\mathbf x-\mathbf x_a|,
\qquad
r_b=|\mathbf x-\mathbf x_b|
\]

を構成する。

コード中の `Ustatic` は

\[
U_{\rm static}
=\frac{M_a}{r_a}+\frac{M_b}{r_b}
\]

に対応する。

したがって時刻依存を含む基本量は

\[
W(t,\mathbf x)=a(t)+U_{\rm static}(\mathbf x)
\]

である。

今回の既定値は

\[
M_a=M_b=1,
\qquad
Q_a=Q_b=1,
\qquad
d=8,
\qquad
H=-0.05,
\]

であり、Kastor–Traschen の extremal 条件

\[
|Q_A|=M_A
\]

を採用している。

### 実装上の意味

`Ustatic` を時刻ループの外で一度だけ作るのは、点源位置と質量パラメータが座標上固定であり、時間依存がすべて

\[
a(t)=e^{Ht}
\]

側に入る Kastor–Traschen 表現に対応する。

`eps=1e-15` は puncture でのゼロ除算回避用であり、物理モデルの新しいパラメータではない。後段で中心近傍は `exclusion` により制約評価から除外される。

---

# 3. 実行系：各時刻の exact state 再構成

## 3.1 時刻依存の生成

実装主要部：

```python
for t in times:
    a=math.exp(H*t)
    W=a+Ustatic
    alpha=a/W
```

対応式はそのまま

\[
a(t)=e^{Ht},
\]

\[
W(t,\mathbf x)
=a(t)+\frac{M_a}{r_a}+\frac{M_b}{r_b},
\]

\[
\alpha(t,\mathbf x)=\frac{a(t)}{W(t,\mathbf x)}.
\]

ここで重要なのは、時刻 `t=0,4,8` の状態が前時刻の数値状態から積分生成されていないことである。

すなわち、コードの時間処理は

\[
X(t_0)\rightarrow X(t_1)\rightarrow X(t_2)
\]

という numerical evolution ではなく、

\[
t_k\mapsto X_{\rm exact}(t_k)
\]

という **exact solution sampling** である。

---

## 3.2 空間微分の有限差分

実装：

```python
def centered_grad(f, h):
    gx = (f[2:,1:-1,1:-1]-f[:-2,1:-1,1:-1])/(2*h)
    gy = (f[1:-1,2:,1:-1]-f[1:-1,:-2,1:-1])/(2*h)
    gz = (f[1:-1,1:-1,2:]-f[1:-1,1:-1,:-2])/(2*h)
    return gx,gy,gz

def centered_lap(f,h):
    c=f[1:-1,1:-1,1:-1]
    return ((f[2:,1:-1,1:-1]+f[:-2,1:-1,1:-1]
            +f[1:-1,2:,1:-1]+f[1:-1,:-2,1:-1]
            +f[1:-1,1:-1,2:]+f[1:-1,1:-1,:-2]-6*c)/(h*h))
```

これは標準的な2次精度中心差分

\[
\partial_x f(x_i)
\approx
\frac{f_{i+1}-f_{i-1}}{2h}
\]

および7点 stencil の3次元 Laplacian

\[
\nabla^2f
\approx
\frac{f_{i+1,j,k}+f_{i-1,j,k}
+f_{i,j+1,k}+f_{i,j-1,k}
+f_{i,j,k+1}+f_{i,j,k-1}
-6f_{i,j,k}}{h^2}
\]

に対応する。

この有限差分は **exact solution を生成するためではなく、exact solution が Einstein–Maxwell 制約を満たすかを独立に検査するため**にだけ使用している。

---

# 4. 読出し系1：Einstein Hamiltonian constraint

実装主要部：

```python
gx,gy,gz=centered_grad(W,h)
lap=centered_lap(W,h)
Wi=W[1:-1,1:-1,1:-1]
grad2=gx*gx+gy*gy+gz*gz

R3=-4*lap/(Wi**3)+2*grad2/(Wi**4)
rho_em=grad2/(8*math.pi*Wi**4)

Ham=R3+6*H*H-16*math.pi*rho_em-2*Lam
```

## 4.1 3次元 Ricci scalar

空間計量は

\[
\gamma_{ij}=W^2\delta_{ij}
\]

である。

これに対する3次元 Ricci scalar は

\[
{}^{(3)}R
=-\frac{4\nabla^2W}{W^3}
+\frac{2|\nabla W|^2}{W^4}.
\]

コードの

```python
R3=-4*lap/(Wi**3)+2*grad2/(Wi**4)
```

はこの式を直接実装している。

## 4.2 電磁エネルギー密度

\[
E^i=-\frac{\partial_iW}{W^3}
\]

より

\[
E_iE^i=\frac{|\nabla W|^2}{W^4}.
\]

磁場ゼロなので

\[
\rho_{\rm EM}
=\frac{E_iE^i}{8\pi}
=\frac{|\nabla W|^2}{8\pi W^4}.
\]

コード：

```python
rho_em=grad2/(8*math.pi*Wi**4)
```

と一致する。

## 4.3 Hamiltonian constraint

Einstein–Maxwell–\(\Lambda\) の Hamiltonian constraint を

\[
{}^{(3)}R+K^2-K_{ij}K^{ij}
-16\pi\rho_{\rm EM}-2\Lambda=0
\]

とする。

\[
K_{ij}=-H\gamma_{ij}
\]

だから

\[
K=-3H,
\qquad
K^2=9H^2,
\qquad
K_{ij}K^{ij}=3H^2.
\]

したがって

\[
K^2-K_{ij}K^{ij}=6H^2.
\]

さらに

\[
\Lambda=3H^2
\]

なので

\[
6H^2-2\Lambda=0.
\]

コード

```python
Ham=R3+6*H*H-16*math.pi*rho_em-2*Lam
```

はこの制約式をそのまま評価している。

さらに式を代入すると

\[
16\pi\rho_{\rm EM}
=\frac{2|\nabla W|^2}{W^4},
\]

したがって gradient 二乗項が相殺され、

\[
\mathcal H
=-\frac{4\nabla^2W}{W^3}.
\]

点源以外では

\[
\nabla^2W=0
\]

であるため、exact solution では

\[
\mathcal H=0
\]

となる。

このことから、今回の Hamiltonian residual は実質的に **有限差分で Laplace 条件をどれだけ正確に再現できるか**を測っていることも明確である。

---

# 5. 読出し系2：Maxwell Gauss constraint

実装：

```python
Gauss=-lap/(Wi**3)
```

Maxwell Gauss constraint は、点源を除けば

\[
D_iE^i=0.
\]

\[
\sqrt{\gamma}=W^3,
\qquad
E^i=-\frac{\partial_iW}{W^3}
\]

より

\[
D_iE^i
=\frac{1}{\sqrt\gamma}
\partial_i(\sqrt\gamma E^i)
\]

\[
=\frac1{W^3}\partial_i(-\partial_iW)
=-\frac{\nabla^2W}{W^3}.
\]

したがってコード

```python
Gauss=-lap/(Wi**3)
```

は Gauss constraint を直接実装している。

Hamiltonian residual と Gauss residual が同じ \(\nabla^2W\) に支配されることも、この厳密解における Einstein sector と Maxwell sector の整合性を示す。

---

# 6. 読出し系3：Momentum constraint

コードには

```python
# exact momentum residual is identically zero for K_ij=-H gamma_ij, S_i=0
```

とあり、出力では

```python
"Momentum":{"max_abs":0.0,"rms":0.0}
```

としている。

これは numerical evaluation ではなく解析的設定である。

Momentum constraint は

\[
D_j(K^{ij}-\gamma^{ij}K)=8\pi S^i.
\]

\[
K^{ij}-\gamma^{ij}K
=(-H+3H)\gamma^{ij}
=2H\gamma^{ij}
\]

であり、metric compatibility

\[
D_j\gamma^{ij}=0
\]

と \(H=\text{const.}\) より左辺は0。

また今回の KT slice では Poynting momentum density

\[
S^i=0
\]

なので制約は解析的に満たされる。

したがって、このコードは Momentum constraint を**独立有限差分検証してはいない**。論文ではその点を明示する。

---

# 7. 読出し系4：puncture 除外と residual 統計

実装：

```python
r1=rs[0][1:-1,1:-1,1:-1]
r2=rs[1][1:-1,1:-1,1:-1]
mask=(r1>exclusion)&(r2>exclusion)

def stats(v):
    q=np.abs(v[mask])
    return {
        "max_abs":float(q.max()),
        "rms":float(np.sqrt(np.mean(q*q))),
        "median_abs":float(np.median(q))
    }

source_scale=np.abs(R3)+np.abs(16*math.pi*rho_em)+2*Lam+1e-30
normHam=np.abs(Ham)/source_scale
```

点源位置では

\[
\nabla^2\frac1r=-4\pi\delta^{(3)}(\mathbf r)
\]

となるため、vacuum/electrovac constraint

\[
\nabla^2W=0
\]

を検証する領域から puncture 近傍を除く必要がある。

そこで

\[
r_a>r_{\rm exclusion},
\qquad
r_b>r_{\rm exclusion}
\]

だけを residual 統計に用いる。

`max_abs`, `rms`, `median_abs` を同時保存するため、局所的な worst case と領域全体の平均的誤差を分離して確認できる。

---

# 8. 読出し系5：動的二中心距離

実装：

```python
bg_sep=a*sep
```

Kastor–Traschen 座標では中心の coordinate separation は `sep=d` で固定だが、背景 scale factor による物理距離の代表量は

\[
d_{\rm bg}(t)=a(t)d=e^{Ht}d
\]

である。

今回 \(H<0\) なので

\[
\frac{d}{dt}d_{\rm bg}=Hd_{\rm bg}<0,
\]

すなわち収縮枝で背景物理距離は減少する。

注意：`background_physical_separation=a*sep` は **二つの horizon 間の厳密 proper distance を積分した値ではない**。scale-factor による背景距離の読出しである。論文ではこの名称を維持し、単純に「二体間 proper distance」とは書かない。

---

# 9. 読出し系6：aa, ab, bb の非線形同時存在

実装主要部：

```python
i0=N//2
W0=float(W[i0,i0,i0])
alpha0=float(alpha[i0,i0,i0])

ua=float(m1/max(rs[0][i0,i0,i0],eps))
ub=float(m2/max(rs[1][i0,i0,i0],eps))

decomp={
    "a2":a*a,
    "linear_a":2*a*ua,
    "linear_b":2*a*ub,
    "aa":ua*ua,
    "ab":2*ua*ub,
    "bb":ub*ub,
    "sum":a*a+2*a*ua+2*a*ub+ua*ua+2*ua*ub+ub*ub,
    "W2":W0*W0
}
```

中点で

\[
u_a=\frac{M_a}{r_a},
\qquad
u_b=\frac{M_b}{r_b}
\]

とすると

\[
W=a+u_a+u_b.
\]

したがって

\[
W^2
=(a+u_a+u_b)^2
\]

\[
=a^2+2au_a+2au_b+u_a^2+2u_au_b+u_b^2.
\]

コードで

\[
aa=u_a^2,
\qquad
ab=2u_au_b,
\qquad
bb=u_b^2
\]

を別々に読み出しているが、これらを独立した相互作用として後から加えてはいない。

先に full nonlinear metric factor

\[
W^2
\]

が存在し、その展開成分として `aa`, `ab`, `bb` を**後から読出している**。

これは論文中の主張

> 「aa, ab, bb は三つの別個の力を手動合成したものではなく、一つの nonlinear field の自己項・交差項として同時に現れる」

と実装上対応している。

また

```python
"sum": ...,
"W2": W0*W0
```

を同時保存しているため、

\[
W^2
\stackrel{?}{=}
a^2+2au_a+2au_b+aa+ab+bb
\]

をデータから直接再検証できる。

---

# 10. 出力・再現性ロジック

実装：

```python
results.append({
    "t":t,
    "a":a,
    "Lambda":Lam,
    "background_physical_separation":bg_sep,
    "midpoint":{
        "W":W0,
        "alpha":alpha0,
        "decomposition":decomp
    },
    "constraints":{
        "Hamiltonian":stats(Ham),
        "Gauss":stats(Gauss),
        "Momentum":{"max_abs":0.0,"rms":0.0},
        "Hamiltonian_normalized":{
            "max":float(normHam[mask].max()),
            "rms":float(np.sqrt(np.mean(normHam[mask]**2)))
        }
    }
})
```

および

```python
if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('--N',type=int,default=49)
    ap.add_argument('--out',default='kt_two_center.json')
    args=ap.parse_args()

    d=run(N=args.N)
    with open(args.out,'w') as f:
        json.dump(d,f,indent=2)
    print(json.dumps(d,indent=2))
```

により、各解像度で同一ロジックを使って JSON を生成する。

本実験では

\[
N=33,
\quad49,
\quad65
\]

を別々に実行し、空間刻み

\[
h=0.75,
\quad0.50,
\quad0.375
\]

で residual の減少を確認した。

このため、結果は単一格子の見かけ上の一致ではなく、**格子細分化による有限差分残差の減少**として検証できる。

---

# 11. 論文主張と実装の対応表

| 論文上の主張 | 対応数式 | コード実装 | 検証可能性 |
|---|---|---|---|
| 二中心 full nonlinear Einstein–Maxwell–Λ exact geometry を使う | \(ds^2=-\Omega^{-2}dt^2+a^2\Omega^2d\mathbf x^2\) | `W=a+Ustatic`, `alpha=a/W` | コードから直接確認可能 |
| 動的状態である | \(a(t)=e^{Ht}\) | `a=math.exp(H*t)` | 各時刻 JSON から確認可能 |
| extremal charge/mass を使う | \(|Q_A|=M_A\) | 出力 metadata `masses=[m1,m2]`, `charges=[m1,m2]` | JSON から確認可能 |
| Einstein constraint を検査する | \({}^{(3)}R+K^2-K_{ij}K^{ij}-16\pi\rho-2\Lambda=0\) | `Ham=...` | N=33,49,65 で収束検証 |
| Maxwell Gauss constraint を検査する | \(D_iE^i=0\) | `Gauss=-lap/W**3` | N=33,49,65 で収束検証 |
| aa,ab,bb は一つの nonlinear field から出る | \((a+u_a+u_b)^2\) | `decomp={aa,ab,bb,...,W2}` | `sum-W2` を直接確認可能 |
| puncture の delta source を vacuum residual から除外する | \(\nabla^2(1/r)=-4\pi\delta^3\) | `mask=(r1>exclusion)&(r2>exclusion)` | コードから確認可能 |
| 時間発展 solver ではなく exact solution evaluator | \(t\mapsto X_{\rm exact}(t)\) | `for t in times: a=exp(H*t); W=...` | 前時刻依存がないことで確認可能 |

---

# 12. このプログラムから主張してよいこと

このコードと保存データから直接支持されるのは次である。

1. 二中心 Kastor–Traschen の **full nonlinear Einstein–Maxwell–Λ exact solution** を3次元格子上に実装した。
2. 動的時間依存 \(a(t)=e^{Ht}\) を含む複数時刻の exact state を評価した。
3. Einstein Hamiltonian constraint と Maxwell Gauss constraint を、exact formula の代入結果をそのままゼロと置くのではなく、**有限差分空間微分から独立に再計算した**。
4. 格子細分化に伴って residual が減少した。
5. `aa, ab, bb` が一つの nonlinear metric factor \(W^2\) に同時に含まれることをコードと出力値から再構成できる。

---

# 13. このプログラムだけからは主張してはいけないこと

以下は本コード単独では検証していない。

1. 一般の asymptotically-flat charged binary merger の3D数値相対論時間発展。
2. 一般二体軌道の生成。
3. \(\Psi_4\) による gravitational-wave extraction。
4. Maxwell Newman–Penrose scalar による electromagnetic radiation extraction。
5. horizon finder による共通 apparent horizon の形成時刻。
6. 任意の \(Q/M\) 比。Kastor–Traschen sector では各中心は extremal \(|Q|=M\) である。
7. `aa`, `ab`, `bb` の各項がそれぞれ独立の物理的 spin-2 radiation channel であること。

これらは別途保存済みの RN coupled perturbation 実験や、将来の production numerical relativity で検証すべき項目である。

---

# 14. 実装監査上の重要点

このプログラムのロジックは、結果を望ましい形へ合わせるために `aa`, `ab`, `bb` を個別に調整していない。

基本入力は

\[
(M_a,M_b,d,H,t)
\]

であり、

\[
W=a+\frac{M_a}{r_a}+\frac{M_b}{r_b}
\]

を一度構成した後、

- geometry,
- electromagnetic energy,
- constraints,
- nonlinear cross terms

をすべて同じ \(W\) から読み出している。

したがって、主要ロジックは

\[
\boxed{
(M_a,M_b,d,H,t)
\longrightarrow W
\longrightarrow
(\gamma_{ij},\alpha,E^i)
\longrightarrow
(\mathcal H,\mathcal G,aa,ab,bb)
}
\]

という一方向の構造になっている。

結果側から初期状態や相互作用係数を逆調整する処理は存在しない。

---

# 15. 再現用ファイル

同一 Google Drive フォルダに以下を保存している。

- `kastor_traschen_two_center_v1.py` — 本体プログラム
- `kt_N33.json`
- `kt_N49.json`
- `kt_N65.json` — 3解像度の生出力
- `kt_N33_stdout.txt`
- `kt_N49_stdout.txt`
- `kt_N65_stdout.txt` — 実行ログ
- `kt_convergence_summary_v1.json` — 収束サマリ
- `requirements_kt_v1.txt` — 依存関係
- `environment_probe.txt` — 実行環境記録
- `SHA256SUMS_v1.txt` — 保存物のハッシュ
- `RUN_STATUS_AND_ANALYSIS_ja_v1.md` — 実験結果分析
- 本文書 `IMPLEMENTATION_LOGIC_AUDIT_KASTOR_TRASCHEN_ja_v1.md`

---

# 16. 参考文献

D. Kastor and J. Traschen, “Cosmological Multi-Black Hole Solutions,” *Physical Review D* **47**, 5370 (1993), arXiv:hep-th/9212035.
