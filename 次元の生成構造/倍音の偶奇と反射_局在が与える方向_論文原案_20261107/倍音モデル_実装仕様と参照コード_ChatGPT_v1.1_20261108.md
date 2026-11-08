# 倍音モデル：実装仕様と再現可能な参照コード

**版：** ChatGPT v1.1、2026-11-08  
**対応論文：** 「倍音の偶奇と反射_局在と方向の読出し_改訂論文案_ChatGPT_v1.1_20261108.md」  
**対象：** 原稿 v0.7、exchange_core_v11.py、Gemini v1.0の改訂。  
**用途：** 研究用の参照実装、実装担当者への仕様、検証基準。

## 1. この文書から実行できること

この文書には、仕様に加えて実行済みPythonコード3本と、実際の結果JSONを全文収録する。設計の説明だけで終了せず、次の機能を動かせる状態とした。

1. 一般四元数波形を全Fourier係数で保存し、単峰・単一軸でなくなっても反復する。
2. 波間の共通作用、次数ごとの作用、次数間のパワー移送を同じ発展関数で比較する。
3. 偶奇、倍音パワー、コヒーレント相関、格子に依存しない局在量を計算する。
4. 既知の共通軸という限定のもとで、欠測倍音を含む位相列を読み、不可識別や別解を返す。
5. 原点共変な三次結合の小実験で、初期に零だった次数の励起を確認する。

固定生成子の検証30項目は全合格。100回の反復の最大相対ノルム二乗ずれは $2.44\times10^{-15}$、100回の逆ステップ後の相対誤差は $1.96\times10^{-15}$ だった。後半の三次結合は別の候補モデルであり、固定生成子の30項目には混ぜていない。

## 2. 状態と計算層の契約

### 2.1 権威ある状態

状態配列は実数の C[N, 2M+1, 4] とする。最後の軸は四元数の $(w,x,y,z)$、中間の軸は次の順序である。

| スロット | 基底 | 倍音次数 |
|---:|---|---:|
| 0 | $1$ | 0 |
| $2m-1$ | $\sqrt2\cos(m\psi)$ | $m$ |
| $2m$ | $\sqrt2\sin(m\psi)$ | $m$ |

ノルムは $\sum C^2=\sum_i\langle|q_i|^2\rangle_\psi$。平均測度 $d\psi/(2\pi)$ に統一する。単一の中心、位相、軸、最大ピーク位置は派生情報であり、全係数の代用品にしない。

波形復元には奇数 NG > 2M、局在量 $J=\langle|q|^4\rangle/\langle|q|^2\rangle^2$ の積分には奇数 NG > 4M を使う。純奇数の設定ではスロット0を含む偶数次数の全係数を零にする。

### 2.2 更新と読出し

~~~text
完全な係数状態 C
  → 指定した生成子による更新 advance(C, K, dt)
  → 完全な係数状態 C_next
  → 一般観測量と、適用条件付きの読出し
  → 次の更新は C_next を使用
~~~

限定読出しの失敗を、位相0や既定軸へ置換しない。参照モデルの更新は限定読出しを呼ばず、読出しが不成立でも続行する。状態依存の生成子を将来導入する場合は、限定デコーダが失敗した場合の勝手な補間ではなく、全状態に定義された作用則を別途指定する。

## 3. 関数と責務

| 関数 | 責務・注意 |
|---|---|
| to_samples / from_samples | 正規直交係数と波形の往復。指定Mを超える検出可能な尾部は無断で切り捨てず拒否する |
| modal_power | 波別・次数別のパワーを返す。全系の各次数は波の軸について合計する |
| parity_weights | DCを偶数に含め、偶奇重みと $\eta$ を返す |
| mode_relations | 全対の同次数四元数内積 $Y_{ij,m}$ |
| coherent_overlap | $|\sum_mY_{ij,m}|/\sqrt{E_iE_j}$。和の後でノルムを取る |
| localization | $J$ と旧指標 $L=J/NG$。単峰性を保証する検査ではない |
| common_axis_wave | 既知軸・位相・中心による初期条件の構成器。更新後もこの形を強制しない |
| channel_generator | 全次数共通または次数別の実反対称ブロックを埋め込む |
| slot_mixing_generator | 基底スロット間の実反対称結合。次数間移送の明示的な試験器 |
| advance | 疎行列指数の作用。Kは実験中固定する |
| left_symmetry_error | 共通の左四元数乗法との可換性を別途検査する |
| decode_common_axis | 既知共通軸・既知振幅積符号に限定した円周位相の読出し |

零パワー状態の規格化観測量はPython配列内ではNaNで返す。これを外部JSONへ保存する実用アプリでは、未定義のstatusとnullに変換する。数値0へ変換しない。

## 4. 位相読出しの仕様

### 4.1 前提

共通軸は明示入力である。四元数の横方向成分、余弦・正弦係数の回転向き、位相直線への残差を検査する。振幅積の符号 product_signs は校正情報として必須とする。同一テンプレート同士なら有効次数で+1を与えられる。未知の符号を絶対値で消してはならない。

左共役の関係から、

$$
z_m=
\frac{\operatorname{Sc}Y_{ij,m}-i\,\operatorname{Vec}Y_{ij,m}\cdot n}
{|Y_{ij,m}|}
\times\operatorname{sign}(c_{i,m}c_{j,m})
$$

を作り、$e^{i(\Phi-m\Delta)}$ に適合する。DCも校正済みで有意なら支持へ含める。支持を先に抽出し、零成分を位相展開へ入れない。

### 4.2 返り値

| status | 意味 | 呼出し側の処理 |
|---|---|---|
| identified | 前提と残差条件を満たし、支持のgcdが1 | このモデル・校正条件のもとでの候補として使用 |
| aliased | gcdが2以上 | 同値な全候補を維持。最初の候補だけを真値としない |
| insufficient_modes | 有効支持が0または1個 | 中心と内部位相を分離しない |
| missing_sign_calibration | 倍音振幅積の符号が未校正 | 読出しを確定しない |
| outside_calibrated_family | 共通既知軸・回転向きの前提外 | 一般状態として保存し、限定読出しを使わない |
| model_mismatch | 円周位相直線に適合しない | 残差を記録し、中心・時計に意味を与えない |

円周フィットはFFTで局所極大候補を列挙し、各候補を補正して最良を選ぶ有限数値法である。任意の雑音分布に対する統計的推定保証は付けていない。既定の相対支持閾値は $10^{-10}$、読出し適合許容値は $10^{-8}$。これらは物理定数ではなく、実験ごとに記録する設定値である。

多軸の未知枠を完全再構成する関数は本版にはない。旧pair_frameをSVD付きにするだけでその問題が解決したとは扱わない。一般状態では、全係数から遅延相関や複数の候補ピークを読む機能を次段として追加できる。

## 5. 三種類の作用を比較する

### 5.1 全次数共通

4N×4Nの実反対称ブロックをchannel_generatorへ一つ渡す。全ての次数に同じ作用を掛けるので、各次数の全系パワーと総強度分布が保存される。旧モデルの保存則を再現する対照実験となる。

### 5.2 次数別

辞書 {m: block} を渡す。未指定次数は生成子零。同じ次数の余弦・正弦スロットには同じブロックを使う。次数別パワーは保存するが、倍音間の干渉により全体の強度分布が変化し得る。

### 5.3 次数間

slot_mixing_generatorの辺は (source_wave, source_slot, target_wave, target_slot, rate) とする。例えばm=1とm=2の対応する余弦・正弦にrate=0.37を与える。

~~~python
import numpy as np
from reference_core import advance, modal_power, slot_mixing_generator

C = np.zeros((1, 5, 4))       # N=1, M=2
C[0, 1, 0] = 1.0            # m=1 cosine, total power 1
K = slot_mixing_generator(C.shape, [
    (0, 1, 0, 3, 0.37),      # cosine 1 -> cosine 2
    (0, 2, 0, 4, 0.37),      # sine   1 -> sine   2
])
C_next = advance(C, K, 0.8)
print(modal_power(C_next))
~~~

出力は約 [0, 0.914913147958, 0.085086852042]。表現帯域内の未占有モードを励起している。帯域の外へ無制限に次数を増やす実装ではない。固定したこの結合は基底の原点・位相規約を選ぶ能力検査であり、原点に依存しない自律的な相互作用則としては採用しない。

疎行列指数を用いるため、4N(2M+1)次元の密な指数行列全体を作る必要はない [X1]。ただし計算量・メモリ量は結合の疎密と帯域に依存し、大規模運用の性能測定は本版の検証に含まない。

## 6. 原点共変な三次結合の候補

nonlinear_probe.py は、固定生成子とは別の候補として

$$
\dot z_m=-i\kappa\sum_{a+b-c=m}z_a z_b\overline{z_c},
\qquad -M\le m\le M
$$

を直接畳込みで実装する。原点変更による $z_m$ の位相が選択則に従って変換され、ODEとして全パワーを保存する。正負両方の複素次数を持ち、reference_core.pyの非負次数・実四元数基底とは配列表現が異なる。

複素軸を四元数の $x$ 方向へ埋め込む場合は、$q(\psi)=(\operatorname{Re}z(\psi),\operatorname{Im}z(\psi),0,0)$ を作り、十分な奇数格子でfrom_samplesを通せば同じ全波形状態へ変換できる。この変換は実装接続の方法であり、複素軸が関係から創発したという意味ではない。負の複素次数は実基底の対応する同じ絶対次数に収容される。

実験条件はM=6、κ=1、初期z1=z2=1/√2、終点0.5、DOP853、rtol=1e-12、atol=1e-14。初期に零だったm=0とm=3の最終パワーは約0.0298153。101記録点のノルム二乗最大誤差は $2.28\times10^{-13}$ だった。

この数値積分器は厳密なノルム保存法でも時間対称法でもない [X2]。保存量と往復誤差を測定して報告する。局所三次位相作用の次数間移送だけでは、局在の輸送や衝突は導かれない。自律的な四波散乱の完成版として扱わない。

## 7. 既存コードからの移行

| 旧処理 | 本版での扱い |
|---|---|
| coefficientsの奇数指定にDCが残る | 初期スペクトルのm=0を偶数として管理する |
| 偶数格子Nyquistを通常倍音と同じ係数で扱う | 参照経路は奇数格子へ統一。偶数格子を将来受ける場合は別途正しい規格化を実装 |
| read_pairで全次数unwrap後にmask | 有効支持を先に選び、円周位相へ適合。支持差gcdによる別解を返す |
| pair_frameの階数1を同軸と決める | 汎用推定を保留し、校正済み共通軸デコーダを独立実装 |
| 作用後も同じパルスのパラメータを強制 | 全係数状態を保持し、適用域外を正常な読出し結果にする |
| 倍音別ノルムを先に合算するg | 倍音の四元数和を作ってからノルム |
| 全ψ共通Uだけで高次生成を期待 | 共通・次数別・次数間を別の要求として実装 |
| 毎回Uが直交だから全更新も可逆とする | 固定Kの基準系で逆を検証。状態依存系のODE・積分法は別評価 |
| 手で中心を更新して並進を補う | 波形の更新後に観測する。中心値を書き換えて現象を作らない |
| 検査用恒等式の一致を物理予測の証拠にする | 代数確認・可識別性・力学的現象・実験との比較を別の合否にする |

旧exchange_core_v11.pyは対照として残す。本版の基盤を別モジュールとして導入し、先行7/13実験のポート、位相規約、作用を照合してから移す。単独二体の写像から多体生成子の対数枝が一意に決まったとは扱わない。

## 8. 次段の受入条件

この版は固定Kの一般波形発展を実装済みとする。次の自律散乱モデルは、次を満たしてから「実装済み」とする。

1. 全一般状態に対して更新を定義し、零相関・退化した枠にも数式上の意味がある。
2. 関係からKを作る規則、保存したい量、位相・ポート・対数枝を明示する。
3. 波ラベル置換、共通四元数変換、円周原点変更について要求する共変性を検査する。
4. 状態依存ODEと数値積分法を区別し、刻み変更で誤差の振舞いを確認する。
5. 帯域端のパワーを監視し、Mを増やして観測結果が保たれる範囲を報告する。
6. 「局在保持」「接近」「衝突」「結合」「有限回帰」それぞれの観測操作と合格条件を先に定義する。

既存の山が波ラベル間で移ること、総強度分布が変わること、次数間でパワーが移ることを別々に記録する。反復が動くという基盤の完成と、全物理モデルの完成を同じ判定にしない。

## 9. 実行方法

検証環境はPython 3.12.14、NumPy 2.3.5、SciPy 1.17.0。仮想環境で必要な依存を入れ、後述の「ファイル」見出しのコードを同名の3ファイルへ保存して実行する。

~~~bash
python -m pip install numpy==2.3.5 scipy==1.17.0
python verify_reference.py
python nonlinear_probe.py --output nonlinear_results.json
~~~

Markdown全体を保存してある場合は、次のコードで付録の3本を一括抽出できる。出力先は新しいharmonic_reference_v1_1フォルダ。同名の異なるコードがあれば上書きせず停止する。

~~~bash
python - "倍音モデル_実装仕様と参照コード_ChatGPT_v1.1_20261108.md" <<'PY'
from pathlib import Path
import re, sys
text = Path(sys.argv[1]).read_text(encoding="utf-8")
blocks = re.findall(r"^### ファイル: ([a-z_]+\.py)\n\n~~~python\n(.*?)\n~~~", text, re.M | re.S)
expected = {"reference_core.py", "verify_reference.py", "nonlinear_probe.py"}
if len(blocks) != 3 or {name for name, _ in blocks} != expected:
    raise SystemExit("Expected exactly the three reference source blocks.")
destination = Path("harmonic_reference_v1_1")
destination.mkdir(exist_ok=True)
for name, code in blocks:
    target = destination / name
    data = code + "\n"
    if target.exists() and target.read_text(encoding="utf-8") != data:
        raise SystemExit(f"Existing different file: {target}")
    target.write_text(data, encoding="utf-8")
print(destination.resolve())
PY
~~~

抽出先へ移動して二つのコマンドを実行する。verify_reference.pyはreference_results.jsonを生成し、すべての検証条件を満たした場合のみ正常終了する。異なる環境で最後の数桁が変わることはある。記載結果との文字列一致ではなく、記録した公差と論理条件で比較する。

## 10. 参考資料と検証範囲

- [X1] [SciPy公式：expm_multiply](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.expm_multiply.html)
- [X2] [SciPy公式：solve_ivp](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html)
- [X3] [NumPy公式：DFTの規約](https://numpy.org/doc/stable/reference/routines.fft.html)
- [S1] [原exchange_core_v11.py](https://drive.google.com/file/d/10UbjIv4WZs-oOgc1O9eplJk8yE8C7sSP/view)
- [S2] [原check_core_v11.py](https://drive.google.com/file/d/1s22T13dYVHBGGpC65ppiH3DWyDAr_UnI/view)
- [S3] [Gemini v1.0原稿](https://drive.google.com/file/d/1y9PUDtQ_LbQbyoPwyW0fAnMiodZyer8C/view)

確認日2026-11-08。元コードの主要15関数は、以前の原稿抽出関数と構文木で一致した。元exchange_core_v11.pyのSHA-256は3ef78547e540847d8fa086f4aca28b942bbd15be1d6273748b80dd55d9ec2d18。

付録コードは実際に実行した参照ファイルをそのまま収録する。既存プロジェクト全体を改修したものではなく、移行・比較のための独立した基盤である。


## 付録A：実行済みの参照コード

次のSHA-256は各PythonブロックをUTF-8、末尾改行一つで保存した内容に対応する。

| ファイル | SHA-256 |
|---|---|
| reference_core.py | 9c3374c413501b10cad2d859b6338d2e05ec985e334dc9ee60a4786279b8b0c3 |
| verify_reference.py | e8892c13be071f7f1851265aab781f50cc1551e480ec7d36afc0d67db3672b1e |
| nonlinear_probe.py | 2983807293e0a55bf42a713409921ad67902b5ecacffbb61f1a215ecd0dc1d08 |

### ファイル: reference_core.py

~~~python
"""Finite-band quaternion waves with a fixed real-skew generator.

C.shape == (N, 2*M+1, 4); slots are [1, sqrt(2)cos(psi),
sqrt(2)sin(psi), ...]. Norms use the normalized circle measure dpsi/(2*pi).
The complete coefficient state evolves; an optional decoder never controls it.
This is an explicit reference model, not a derivation of a unique physical law.
"""
import math
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply


def qmul(a, b):
    """Hamilton product, scalar first, broadcasting over leading dimensions."""
    a, b = np.broadcast_arrays(np.asarray(a), np.asarray(b))
    s = a[..., 0]*b[..., 0] - np.sum(a[..., 1:]*b[..., 1:], axis=-1)
    v = a[..., :1]*b[..., 1:] + b[..., :1]*a[..., 1:] + np.cross(a[..., 1:], b[..., 1:])
    return np.concatenate([s[..., None], v], axis=-1)


def qconj(q):
    return np.asarray(q)*np.array([1., -1., -1., -1.])


def _state(C):
    if np.iscomplexobj(C):
        raise ValueError("Quaternion components must be real.")
    C = np.asarray(C, dtype=float)
    if C.ndim != 3 or C.shape[-1] != 4 or C.shape[0] < 1 or C.shape[1] % 2 != 1:
        raise ValueError("Expected shape (N, 2*M+1, 4).")
    if not np.all(np.isfinite(C)):
        raise ValueError("State contains nonfinite values.")
    return C


def _grid(ng, M):
    if not isinstance(ng, (int, np.integer)) or ng % 2 != 1 or ng <= 2*M:
        raise ValueError("Use an odd sample count NG > 2*M.")


def to_samples(C, ng):
    C = _state(C); N, B, _ = C.shape; M = (B-1)//2
    _grid(ng, M)
    F = np.zeros((N, ng//2+1, 4), dtype=complex)
    F[:, 0] = C[:, 0]
    F[:, 1:M+1] = (C[:, 1::2] - 1j*C[:, 2::2])/np.sqrt(2.)
    return np.fft.irfft(F*ng, n=ng, axis=1)


def from_samples(Q, M, tolerance=1e-12):
    """Reject an unresolved spectral tail; never silently truncate it."""
    if np.iscomplexobj(Q):
        raise ValueError("Quaternion components must be real.")
    Q = np.asarray(Q, dtype=float)
    if Q.ndim != 3 or Q.shape[-1] != 4 or not np.all(np.isfinite(Q)):
        raise ValueError("Expected finite samples with shape (N, NG, 4).")
    if not isinstance(M, (int, np.integer)) or M < 0:
        raise ValueError("M must be a nonnegative integer.")
    N, ng, _ = Q.shape; _grid(ng, M)
    F = np.fft.rfft(Q, axis=1)/ng
    tail = 2*np.sum(np.abs(F[:, M+1:])**2)
    energy = np.mean(np.sum(Q*Q, axis=(0, 2)))
    if tail > tolerance**2*max(energy, np.finfo(float).tiny):
        raise ValueError("Samples contain modes above M; increase the state bandwidth.")
    C = np.zeros((N, 2*M+1, 4))
    C[:, 0] = F[:, 0].real
    C[:, 1::2] = np.sqrt(2.)*F[:, 1:M+1].real
    C[:, 2::2] = -np.sqrt(2.)*F[:, 1:M+1].imag
    return C


def modal_power(C):
    C = _state(C); N, B, _ = C.shape
    p = np.sum(C*C, axis=-1)
    return np.concatenate([p[:, :1], p[:, 1:].reshape(N, (B-1)//2, 2).sum(axis=-1)], axis=1)


def parity_weights(C):
    """Return (even, odd, eta); zero-energy waves have undefined (NaN) weights."""
    p = modal_power(C); energy = p.sum(axis=1)
    even = np.divide(p[:, ::2].sum(axis=1), energy, out=np.full_like(energy, np.nan), where=energy > 0)
    odd = np.divide(p[:, 1::2].sum(axis=1), energy, out=np.full_like(energy, np.nan), where=energy > 0)
    return even, odd, even-odd


def mode_relations(C):
    C = _state(C); N, B, _ = C.shape
    slots = qmul(qconj(C[:, None]), C[None, :])
    return np.concatenate([slots[:, :, :1], slots[:, :, 1:].reshape(N, N, (B-1)//2, 2, 4).sum(axis=3)], axis=2)


def coherent_overlap(C):
    """Norm AFTER summing modes; undefined for a zero-energy member of a pair."""
    Y = mode_relations(C)
    energy = modal_power(C).sum(axis=1)
    denominator = np.sqrt(energy[:, None]*energy[None, :])
    return np.divide(np.linalg.norm(Y.sum(axis=2), axis=-1), denominator,
                     out=np.full_like(denominator, np.nan), where=denominator > 0)


def localization(C, ng=None):
    """J=<|q|^4>/<|q|^2>^2, and old discrete L=J/NG. Zero waves: NaN.

    NG>4M gives exact trapezoidal integration of this fourth-order quantity.
    """
    C = _state(C); M = (C.shape[1]-1)//2
    ng = 4*M+1 if ng is None else ng
    if ng <= 4*M:
        raise ValueError("Localization requires odd NG > 4*M.")
    density = np.sum(to_samples(C, ng)**2, axis=-1)
    energy = modal_power(C).sum(axis=1)
    J = np.divide(np.mean(density*density, axis=1), energy**2,
                  out=np.full_like(energy, np.nan), where=energy > 0)
    return J, J/ng


def common_axis_wave(amplitudes, theta, phase, axis):
    """One calibrated initial wave, angles in radians. DC is amplitudes[0].

    Pure odd parity requires amplitudes[0::2]=0, INCLUDING DC.
    This initializer is optional; subsequent states need not retain its form.
    """
    a = np.asarray(amplitudes, dtype=float); n = np.asarray(axis, dtype=float)
    if a.ndim != 1 or len(a) < 1 or n.shape != (3,) or np.linalg.norm(n) == 0:
        raise ValueError("Supply a nonempty real spectrum and a nonzero 3-axis.")
    n = n/np.linalg.norm(n); nq = np.r_[0., n]; M = len(a)-1
    alpha = phase - np.arange(M+1)*theta
    modes = a[:, None]*np.column_stack([np.cos(alpha), np.sin(alpha)[:, None]*n])
    C = np.zeros((1, 2*M+1, 4)); C[0, 0] = modes[0]
    C[0, 1::2] = modes[1:]/np.sqrt(2.)
    C[0, 2::2] = qmul(modes[1:], nq)/np.sqrt(2.)
    return _state(C)


def _check_generator(K, shape, tolerance=1e-12):
    if np.iscomplexobj(K.data if sparse.issparse(K) else K):
        raise ValueError("This engine requires a REAL skew generator.")
    K = sparse.csr_matrix(K, dtype=float); size = int(np.prod(shape))
    if K.shape != (size, size) or not np.all(np.isfinite(K.data)):
        raise ValueError("Generator shape or entries are invalid.")
    defect = K + K.T
    scale = max(1., float(np.max(np.abs(K.data), initial=0.)))
    if np.max(np.abs(defect.data), initial=0.) > tolerance*scale:
        raise ValueError("Generator must satisfy K.T = -K.")
    return K


def channel_generator(shape, blocks):
    """Embed a 4N x 4N block in each modal sector.

    A single matrix applies to ALL modes; a dict {m: matrix} is mode-dependent.
    Missing dict entries are zero. The same block acts on cos and sin of m.
    Blocks are explicit fixed inputs, not recomputed from the evolving state.
    """
    N, B, four = shape; M = (B-1)//2
    if four != 4 or B % 2 != 1:
        raise ValueError("Expected a coefficient-state shape.")
    by_mode = blocks if isinstance(blocks, dict) else {m: blocks for m in range(M+1)}
    rows, cols, data = [], [], []
    for mode, block in by_mode.items():
        if not isinstance(mode, (int, np.integer)) or not 0 <= mode <= M:
            raise ValueError("Mode index is outside the represented bandwidth.")
        block = _check_generator(block, (N, 1, 4)).tocoo()
        for slot in ([0] if mode == 0 else [2*mode-1, 2*mode]):
            rows.extend(((block.row//4*B+slot)*4+block.row%4).tolist())
            cols.extend(((block.col//4*B+slot)*4+block.col%4).tolist())
            data.extend(block.data.tolist())
    size = N*B*4
    return _check_generator(sparse.coo_matrix((data, (rows, cols)), shape=(size, size)), shape)


def slot_mixing_generator(shape, edges):
    """Edges (i,b,j,c,rate) rotate source slot (i,b) into target (j,c).

    Add K[target,source]=rate and K[source,target]=-rate for EACH quaternion
    component. b,c are real-basis slot indices, not harmonic numbers.
    m!=n mixing requires edges between different harmonic sectors.
    """
    N, B, four = shape; rows, cols, data = [], [], []
    for i, b, j, c, rate in edges:
        if four != 4 or not (0 <= i < N and 0 <= j < N and 0 <= b < B and 0 <= c < B):
            raise ValueError("Edge endpoint outside the coefficient state.")
        for a in range(4):
            source, target = (i*B+b)*4+a, (j*B+c)*4+a
            rows.extend([target, source]); cols.extend([source, target]); data.extend([rate, -rate])
    size = N*B*4
    return _check_generator(sparse.coo_matrix((data, (rows, cols)), shape=(size, size)), shape)


def advance(C, K, dt):
    """Apply exp(dt*K) without a dense exponential. Keep K fixed between steps."""
    C = _state(C); K = _check_generator(K, C.shape)
    if not np.isfinite(dt):
        raise ValueError("dt must be finite.")
    vector = expm_multiply(dt*K, C.ravel(), traceA=dt*float(K.diagonal().sum()))
    return vector.reshape(C.shape)


def left_symmetry_error(K, shape):
    """Check commutation with common left multiplication by i,j,k separately."""
    K = _check_generator(K, shape); scale = max(1., np.max(np.abs(K.data), initial=0.))
    errors = []
    for nq in np.eye(4)[1:]:
        L = np.column_stack([qmul(nq, e) for e in np.eye(4)])
        whole = sparse.kron(sparse.eye(int(np.prod(shape))//4), L, format='csr')
        defect = K@whole - whole@K
        errors.append(float(np.max(np.abs(defect.data), initial=0.)/scale))
    return errors


def decode_common_axis(C, axis, product_signs=None, i=0, j=1,
                       relative_floor=1e-10, tolerance=1e-8):
    """OPTIONAL calibrated decoder; does not participate in advance().

    Conditions: a known canonical common complex axis, positive circulation,
    and KNOWN signs of the two real harmonic envelopes' products. For the
    same template, explicitly supply +1 for all modes. Unknown signs are not
    discarded. A successful linear-phase fit is model consistency, not proof
    that arbitrary physical states belong to the initial pulse family.

    z_m=Sc(Y_m)-i*Vec(Y_m).axis fits exp(i*(Phi-m*Delta)). Active modes are
    selected BEFORE phase processing; no unwrap is used. Aliases implied by
    gcd(m-m0) are all returned, and zero/one active mode is insufficient.
    """
    C = _state(C); pair = C[[i, j]]; M = (C.shape[1]-1)//2
    n = np.asarray(axis, dtype=float)
    if n.shape != (3,) or np.linalg.norm(n) == 0:
        raise ValueError("Supply the known nonzero 3-axis.")
    n = n/np.linalg.norm(n); scale = max(np.linalg.norm(pair), np.finfo(float).tiny)
    transverse = pair[..., 1:] - (pair[..., 1:]@n)[..., None]*n
    circulation = pair[:, 2::2] - qmul(pair[:, 1::2], np.r_[0., n])
    if max(np.linalg.norm(transverse), np.linalg.norm(circulation)) > tolerance*scale:
        return {"status": "outside_calibrated_family", "solutions": []}
    Y = mode_relations(pair)[0, 1]
    z = Y[:, 0] - 1j*(Y[:, 1:]@n); magnitude = np.abs(z)
    active = magnitude > relative_floor*np.max(magnitude, initial=0.)
    modes = np.flatnonzero(active)
    if len(modes) < 2:
        return {"status": "insufficient_modes", "modes": modes.tolist(), "solutions": []}
    if product_signs is None:
        return {"status": "missing_sign_calibration", "modes": modes.tolist(), "solutions": []}
    signs = np.asarray(product_signs, dtype=float)
    if signs.shape != (M+1,) or not np.all(np.isin(signs[active], [-1., 1.])):
        raise ValueError("Each active mode needs a calibrated product sign +1 or -1.")
    phasors = z[active]*signs[active]/magnitude[active]
    weights = magnitude[active]/magnitude[active].sum()
    d = int(np.gcd.reduce(modes[1:]-modes[0])); k = (modes-modes[0])//d
    count = 1 << int(math.ceil(math.log2(max(256, 32*(int(k[-1])+1)))))
    polynomial = np.zeros(int(k[-1])+1, dtype=complex); polynomial[k] = weights*phasors
    scores = np.abs(np.fft.ifft(polynomial, n=count))
    peaks = np.flatnonzero((scores >= np.roll(scores, 1)) & (scores >= np.roll(scores, -1)))
    best_score, best_x = -1., 0.
    # Refine every sampled local peak: uneven spectra can have near-alias peaks.
    for peak in peaks:
        x = 2*np.pi*peak/count
        for _ in range(8):
            terms = weights*phasors*np.exp(1j*k*x)
            a = terms.sum(); ap = (1j*k*terms).sum(); app = (-k*k*terms).sum()
            first = 2*np.real(ap*np.conj(a)); second = 2*np.real(app*np.conj(a)+ap*np.conj(ap))
            if second >= 0:
                break
            step = first/second
            if abs(step) > 4*np.pi/count:
                break
            x -= step
            if abs(step) < 1e-13:
                break
        score = abs(np.sum(weights*phasors*np.exp(1j*k*x)))
        if score > best_score:
            best_score, best_x = score, x
    x = best_x
    wrap = lambda angle: float((angle+np.pi) % (2*np.pi)-np.pi)
    delta = wrap(x)/d
    correlation = np.sum(weights*phasors*np.exp(1j*modes*delta))
    phi = float(np.angle(correlation))
    residual = float(np.sqrt(np.sum(weights*np.angle(phasors*np.exp(-1j*(phi-modes*delta)))**2)))
    result = {"status": "model_mismatch", "modes": modes.tolist(), "gcd": d,
              "phase_rms_rad": residual, "coherence": float(abs(correlation)),
              "best_fit": {"Phi": phi, "Delta": delta}, "solutions": []}
    if residual <= tolerance:
        result["status"] = "identified" if d == 1 else "aliased"
        result["solutions"] = [{"Phi": wrap(phi+modes[0]*2*np.pi*ell/d),
                                "Delta": wrap(delta+2*np.pi*ell/d)} for ell in range(d)]
    return result
~~~

### ファイル: verify_reference.py

~~~python
"""Small deterministic verification suite for reference_core.py.

Run: python verify_reference.py
Writes reference_results.json. Uses no network and no dense large exponential.
The tests check mathematical behavior, not a proposed physical law's validity.
"""
import json
import platform
from pathlib import Path
import numpy as np
import scipy
from scipy import sparse
import reference_core as r


def relative(a, b):
    return float(np.linalg.norm(np.asarray(a)-np.asarray(b))/max(np.linalg.norm(b), np.finfo(float).tiny))


def angle_error(a, b):
    return float(abs((a-b+np.pi) % (2*np.pi)-np.pi))


def main():
    tests = {}
    def small(name, value, limit=5e-12):
        tests[name] = {"value": float(value), "upper_bound": limit, "passed": bool(value <= limit)}
    def condition(name, passed, detail):
        tests[name] = {"passed": bool(passed), "detail": detail}

    rng = np.random.default_rng(20261108)
    C = rng.normal(size=(3, 11, 4))
    Q = r.to_samples(C, 65)
    small("fft_round_trip_relative", relative(r.from_samples(Q, 5), C))
    sample_energy = np.mean(np.sum(Q*Q, axis=(0, 2)))
    small("parseval_relative", abs(sample_energy/np.sum(C*C)-1))
    direct = r.qmul(r.qconj(Q[:, None]), Q[None, :]).mean(axis=2)
    small("relations_equal_direct_inner_product", relative(r.mode_relations(C).sum(axis=2), direct))
    small("localization_grid_independent", relative(r.localization(C, 65)[0], r.localization(C, 101)[0]))
    try:
        r.from_samples(Q, 4)
        rejects_tail = False
    except ValueError:
        rejects_tail = True
    condition("no_silent_band_truncation", rejects_tail, "A nonzero m=5 tail cannot be imported as M=4.")

    odd_spectrum = np.zeros(9); odd_spectrum[1::2] = 0.5
    odd = r.common_axis_wave(odd_spectrum, .2, .3, [0, 0, 1])
    we, wo, eta = r.parity_weights(odd)
    small("pure_odd_dc_and_even_weight", max(abs(we[0]), np.linalg.norm(odd[:, 0])))
    small("pure_odd_eta_minus_one", abs(eta[0]+1))

    spectrum = np.r_[0., np.ones(8)/np.sqrt(8)]
    wave0 = r.common_axis_wave(spectrum, 0., 0., [0, 0, 1])
    identical = np.concatenate([wave0, wave0])
    separated = np.concatenate([wave0, r.common_axis_wave(spectrum, 2*np.pi/8, 0., [0, 0, 1])])
    g_same = float(r.coherent_overlap(identical)[0, 1])
    g_cancel = float(r.coherent_overlap(separated)[0, 1])
    small("coincident_overlap_is_one", abs(g_same-1))
    small("coherent_distance_cancellation", g_cancel)

    A = np.array([[0., -.2, .3], [.2, 0., -.4], [-.3, .4, 0.]])
    H = sparse.kron(A, sparse.eye(4), format='csr')
    Kcommon = r.channel_generator(C.shape, H)
    Ccommon = r.advance(C, Kcommon, .31)
    rho0 = np.sum(Q*Q, axis=(0, 2))
    rho1 = np.sum(r.to_samples(Ccommon, 65)**2, axis=(0, 2))
    small("common_generator_density_preserved", relative(rho1, rho0))
    small("common_generator_mode_power_preserved", relative(r.modal_power(Ccommon).sum(axis=0), r.modal_power(C).sum(axis=0)))
    small("common_generator_left_symmetry", max(r.left_symmetry_error(Kcommon, C.shape)))

    # Two harmonics initially share one output channel. Only m=2 is rotated
    # into the other channel, removing an inter-harmonic density cross term.
    D = np.zeros((2, 5, 4)); D[0, 1, 0] = D[0, 3, 0] = 1.
    H2 = sparse.kron(np.array([[0., -np.pi/2], [np.pi/2, 0.]]), sparse.eye(4), format='csr')
    Kmodal = r.channel_generator(D.shape, {2: H2})
    D1 = r.advance(D, Kmodal, 1.)
    density_change = relative(np.sum(r.to_samples(D1, 33)**2, axis=(0, 2)), np.sum(r.to_samples(D, 33)**2, axis=(0, 2)))
    small("mode_dependent_power_preserved", relative(r.modal_power(D1).sum(axis=0), r.modal_power(D).sum(axis=0)))
    condition("mode_dependent_density_can_change", density_change > .1, {"relative_density_change": density_change})

    F = r.common_axis_wave([0., 1., 0.], 0., 0., [0, 0, 1])
    rate, dt = .37, .8
    Kmixed = r.slot_mixing_generator(F.shape, [(0, 1, 0, 3, rate), (0, 2, 0, 4, rate)])
    F1 = r.advance(F, Kmixed, dt)
    powers = r.modal_power(F1)[0]
    expected = np.array([0., np.cos(rate*dt)**2, np.sin(rate*dt)**2])
    small("mode_1_to_2_sine_squared_transfer", np.max(np.abs(powers-expected)))
    small("mode_mixing_total_norm", abs(np.sum(F1*F1)/np.sum(F*F)-1))
    small("mode_mixing_left_symmetry", max(r.left_symmetry_error(Kmixed, F.shape)))

    repeated = F.copy(); base_norm = float(np.sum(F*F)); drifts = []
    for _ in range(100):
        repeated = r.advance(repeated, Kmixed, .03)
        drifts.append(abs(np.sum(repeated*repeated)/base_norm-1))
    for _ in range(100):
        repeated = r.advance(repeated, Kmixed, -.03)
    small("hundred_steps_max_norm_drift", max(drifts))
    small("hundred_forward_inverse_steps_return", relative(repeated, F))

    # Optional, explicitly calibrated readout. All angles returned in radians.
    phi, delta = np.deg2rad(20.), np.deg2rad(143.)
    gap = np.zeros(9); gap[[1, 4, 6]] = 1/np.sqrt(3)
    pair = np.concatenate([r.common_axis_wave(gap, delta, phi, [0, 0, 1]),
                           r.common_axis_wave(gap, 0., 0., [0, 0, 1])])
    decoded = r.decode_common_axis(pair, [0, 0, 1], np.ones(9))
    condition("gapped_coprime_phase_identified", decoded["status"] == "identified", decoded)
    if decoded["solutions"]:
        solution = decoded["solutions"][0]
        small("gapped_phase_recovery", max(angle_error(solution["Phi"], phi), angle_error(solution["Delta"], delta)))
    else:
        condition("gapped_phase_recovery", False, "No solution returned.")

    odd_pair = np.concatenate([r.common_axis_wave(odd_spectrum, np.deg2rad(30.), phi, [0, 0, 1]),
                               r.common_axis_wave(odd_spectrum, 0., 0., [0, 0, 1])])
    aliases = r.decode_common_axis(odd_pair, [0, 0, 1], np.ones(9))
    condition("odd_support_two_aliases", aliases["status"] == "aliased" and aliases["gcd"] == 2 and len(aliases["solutions"]) == 2, aliases)
    small("odd_support_true_solution_present", min(max(angle_error(s["Phi"], phi), angle_error(s["Delta"], np.deg2rad(30.))) for s in aliases["solutions"]))

    negative = gap.copy(); negative[4] *= -1
    signed_pair = np.concatenate([r.common_axis_wave(negative, delta, phi, [0, 0, 1]),
                                  r.common_axis_wave(gap, 0., 0., [0, 0, 1])])
    signs = np.ones(9); signs[4] = -1
    signed = r.decode_common_axis(signed_pair, [0, 0, 1], signs)
    wrong_sign = r.decode_common_axis(signed_pair, [0, 0, 1], np.ones(9))
    condition("known_negative_product_sign_corrected", signed["status"] == "identified", signed)
    condition("wrong_sign_model_rejected", wrong_sign["status"] == "model_mismatch", {"status": wrong_sign["status"], "phase_rms_rad": wrong_sign["phase_rms_rad"]})
    missing = r.decode_common_axis(pair, [0, 0, 1])
    condition("unknown_signs_not_guessed", missing["status"] == "missing_sign_calibration", missing)

    zero = r.decode_common_axis(np.zeros_like(pair), [0, 0, 1], np.ones(9))
    single_spectrum = np.zeros(9); single_spectrum[4] = 1.
    single_wave = r.common_axis_wave(single_spectrum, 0., 0., [0, 0, 1])
    single = r.decode_common_axis(np.concatenate([single_wave, single_wave]), [0, 0, 1], np.ones(9))
    condition("zero_modes_underdetermined", zero["status"] == "insufficient_modes", zero)
    condition("one_mode_underdetermined", single["status"] == "insufficient_modes", single)
    outside = r.decode_common_axis(C, [0, 0, 1], np.ones(6))
    condition("general_state_not_forced_into_decoder", outside["status"] == "outside_calibrated_family", outside)
    small("evolution_works_when_decoder_rejects", abs(np.sum(Ccommon*Ccommon)/np.sum(C*C)-1))

    result = {"scope": "Fixed-generator mathematical reference tests; not empirical validation of a unique interaction law.",
              "environment": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__},
              "benchmarks": {"random_seed": 20261108, "random_shape": list(C.shape),
                             "mode_transfer_rate": rate, "mode_transfer_dt": dt,
                             "mode_transfer_expected_power": expected.tolist(),
                             "mode_transfer_observed_power": powers.tolist(),
                             "coincident_g": g_same, "cancelled_g": g_cancel},
              "tests": tests, "all_passed": all(item["passed"] for item in tests.values())}
    output = Path(__file__).resolve().parent/"reference_results.json"
    output.write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False)+"\n", encoding='utf-8')
    print(json.dumps({"checks": len(tests), "all_passed": result["all_passed"],
                      "norm_drift_100_steps": max(drifts),
                      "inverse_return_relative": relative(repeated, F),
                      "mode_transfer_power": powers.tolist(),
                      "modal_density_change": density_change,
                      "coherent_cancelled_g": g_cancel}, indent=2))
    if not result["all_passed"]:
        raise RuntimeError("Reference verification failed; inspect reference_results.json.")


if __name__ == '__main__':
    main()
~~~

### ファイル: nonlinear_probe.py

~~~python
#!/usr/bin/env python3
"""Finite Galerkin cubic phase-mixing probe, not a derived physical law.

z(psi) = sum_{m=-M}^M z_m exp(i*m*psi)
dz_m/ds = -i*kappa*sum_{a+b-c=m} z_a*z_b*conj(z_c).
The convolution is evaluated exactly for the truncated coefficient arrays;
there is no FFT wraparound. The complex axis is an explicit model input,
not an axis recovered from relations. The evolution parameter is dimensionless.

The ODE conserves sum(abs(z_m)**2) and is covariant under a common shift of
the angular origin. DOP853 only approximates that flow; it is not an exactly
norm-preserving or symmetric numerical method. Measured errors are reported.
This probe does not establish spatial transport, localization, force, or
collision dynamics. No continuum/bandwidth convergence claim is made.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import scipy
from scipy.integrate import solve_ivp


def run_probe(M: int = 6, end: float = 0.5) -> dict:
    modes = np.arange(-M, M + 1)
    kappa, rtol, atol, chi = 1.0, 1e-12, 1e-14, 0.371
    z0 = np.zeros(2 * M + 1, dtype=complex)
    z0[M + 1] = z0[M + 2] = 1.0 / np.sqrt(2.0)

    def rhs(s: float, z: np.ndarray) -> np.ndarray:
        full = np.convolve(np.convolve(z, z), np.conj(z[::-1]))
        return -1j * kappa * full[2 * M:4 * M + 1]

    times = np.linspace(0.0, end, 101)
    options = dict(method="DOP853", rtol=rtol, atol=atol)
    forward = solve_ivp(rhs, (0.0, end), z0, t_eval=times, **options)
    phase = np.exp(-1j * modes * chi)
    shifted = solve_ivp(rhs, (0.0, end), phase * z0, t_eval=times, **options)
    if not forward.success or not shifted.success:
        raise RuntimeError((forward.message, shifted.message))
    backward = solve_ivp(rhs, (end, 0.0), forward.y[:, -1], **options)
    if not backward.success:
        raise RuntimeError(backward.message)

    derivative = rhs(0.0, z0)
    norm0 = float(np.vdot(z0, z0).real)
    norms = np.sum(np.abs(forward.y) ** 2, axis=0)
    final = forward.y[:, -1]
    expected0 = -1j * kappa * z0[M + 1] ** 2 * z0[M + 2].conjugate()
    expected3 = -1j * kappa * z0[M + 2] ** 2 * z0[M + 1].conjugate()
    as_pair = lambda z: [float(z.real), float(z.imag)]
    return {
        "model": "finite_Galerkin_cubic_phase_mixing",
        "M": M, "kappa": kappa, "s_final": end,
        "integrator": {"method": "DOP853", "rtol": rtol, "atol": atol,
                       "exact_norm_preservation_claimed": False,
                       "exact_numerical_reversibility_claimed": False},
        "versions": {"numpy": np.__version__, "scipy": scipy.__version__},
        "initial_norm_squared": norm0,
        "sampled_max_norm_squared_error": float(np.max(np.abs(norms - norm0))),
        "origin_shift_radians": chi,
        "rhs_origin_covariance_max_error": float(np.max(np.abs(
            rhs(0.0, phase * z0) - phase * derivative))),
        "sampled_flow_origin_covariance_max_error": float(np.max(np.abs(
            shifted.y - phase[:, None] * forward.y))),
        "forward_then_backward_max_error":
            float(np.max(np.abs(backward.y[:, -1] - z0))),
        "initial_new_mode_derivatives": {
            "m0": as_pair(derivative[M]), "m3": as_pair(derivative[M + 3]),
            "formula_max_error": float(max(abs(derivative[M] - expected0),
                                           abs(derivative[M + 3] - expected3)))},
        "initial_power_by_mode": {str(m): float(abs(z) ** 2)
                                  for m, z in zip(modes, z0)},
        "final_power_by_mode": {str(m): float(abs(z) ** 2)
                                for m, z in zip(modes, final)},
        "function_evaluations": {"forward": forward.nfev,
                                 "shifted": shifted.nfev, "backward": backward.nfev},
        "limits": ["The complex axis is supplied, not relationally derived.",
                   "Spatial transport, force and collision are not established.",
                   "No bandwidth-to-continuum convergence test is included."]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    rendered = json.dumps(run_probe(), ensure_ascii=False, indent=2)
    print(rendered)
    if args.output:
        args.output.write_text(rendered + "\n", encoding="utf-8")
~~~

## 付録B：実行結果

### reference_results.json

~~~json
{
  "scope": "Fixed-generator mathematical reference tests; not empirical validation of a unique interaction law.",
  "environment": {
    "python": "3.12.14",
    "numpy": "2.3.5",
    "scipy": "1.17.0"
  },
  "benchmarks": {
    "random_seed": 20261108,
    "random_shape": [
      3,
      11,
      4
    ],
    "mode_transfer_rate": 0.37,
    "mode_transfer_dt": 0.8,
    "mode_transfer_expected_power": [
      0.0,
      0.9149131479576741,
      0.08508685204232588
    ],
    "mode_transfer_observed_power": [
      0.0,
      0.9149131479576738,
      0.08508685204232584
    ],
    "coincident_g": 1.0000000000000002,
    "cancelled_g": 5.2441245763822065e-17
  },
  "tests": {
    "fft_round_trip_relative": {
      "value": 1.7241713484932779e-16,
      "upper_bound": 5e-12,
      "passed": true
    },
    "parseval_relative": {
      "value": 1.1102230246251565e-16,
      "upper_bound": 5e-12,
      "passed": true
    },
    "relations_equal_direct_inner_product": {
      "value": 2.6816967359087925e-16,
      "upper_bound": 5e-12,
      "passed": true
    },
    "localization_grid_independent": {
      "value": 3.1992661349303524e-16,
      "upper_bound": 5e-12,
      "passed": true
    },
    "no_silent_band_truncation": {
      "passed": true,
      "detail": "A nonzero m=5 tail cannot be imported as M=4."
    },
    "pure_odd_dc_and_even_weight": {
      "value": 0.0,
      "upper_bound": 5e-12,
      "passed": true
    },
    "pure_odd_eta_minus_one": {
      "value": 2.220446049250313e-16,
      "upper_bound": 5e-12,
      "passed": true
    },
    "coincident_overlap_is_one": {
      "value": 2.220446049250313e-16,
      "upper_bound": 5e-12,
      "passed": true
    },
    "coherent_distance_cancellation": {
      "value": 5.2441245763822065e-17,
      "upper_bound": 5e-12,
      "passed": true
    },
    "common_generator_density_preserved": {
      "value": 2.0126673792377725e-16,
      "upper_bound": 5e-12,
      "passed": true
    },
    "common_generator_mode_power_preserved": {
      "value": 2.1278496186146115e-16,
      "upper_bound": 5e-12,
      "passed": true
    },
    "common_generator_left_symmetry": {
      "value": 0.0,
      "upper_bound": 5e-12,
      "passed": true
    },
    "mode_dependent_power_preserved": {
      "value": 3.14018491736755e-16,
      "upper_bound": 5e-12,
      "passed": true
    },
    "mode_dependent_density_can_change": {
      "passed": true,
      "detail": {
        "relative_density_change": 0.6666666666666667
      }
    },
    "mode_1_to_2_sine_squared_transfer": {
      "value": 3.3306690738754696e-16,
      "upper_bound": 5e-12,
      "passed": true
    },
    "mode_mixing_total_norm": {
      "value": 2.220446049250313e-16,
      "upper_bound": 5e-12,
      "passed": true
    },
    "mode_mixing_left_symmetry": {
      "value": 0.0,
      "upper_bound": 5e-12,
      "passed": true
    },
    "hundred_steps_max_norm_drift": {
      "value": 2.4424906541753444e-15,
      "upper_bound": 5e-12,
      "passed": true
    },
    "hundred_forward_inverse_steps_return": {
      "value": 1.9646231318433117e-15,
      "upper_bound": 5e-12,
      "passed": true
    },
    "gapped_coprime_phase_identified": {
      "passed": true,
      "detail": {
        "status": "identified",
        "modes": [
          1,
          4,
          6
        ],
        "gcd": 1,
        "phase_rms_rad": 2.813125635706675e-16,
        "coherence": 1.0,
        "best_fit": {
          "Phi": 0.34906585039886573,
          "Delta": 2.4958208303518914
        },
        "solutions": [
          {
            "Phi": 0.3490658503988655,
            "Delta": 2.4958208303518914
          }
        ]
      }
    },
    "gapped_phase_recovery": {
      "value": 4.440892098500626e-16,
      "upper_bound": 5e-12,
      "passed": true
    },
    "odd_support_two_aliases": {
      "passed": true,
      "detail": {
        "status": "aliased",
        "modes": [
          1,
          3,
          5,
          7
        ],
        "gcd": 2,
        "phase_rms_rad": 2.7635588085755153e-16,
        "coherence": 1.0,
        "best_fit": {
          "Phi": 0.34906585039886545,
          "Delta": 0.5235987755982987
        },
        "solutions": [
          {
            "Phi": 0.3490658503988655,
            "Delta": 0.5235987755982987
          },
          {
            "Phi": -2.7925268031909276,
            "Delta": -2.617993877991495
          }
        ]
      }
    },
    "odd_support_true_solution_present": {
      "value": 4.440892098500626e-16,
      "upper_bound": 5e-12,
      "passed": true
    },
    "known_negative_product_sign_corrected": {
      "passed": true,
      "detail": {
        "status": "identified",
        "modes": [
          1,
          4,
          6
        ],
        "gcd": 1,
        "phase_rms_rad": 2.813125635706675e-16,
        "coherence": 1.0,
        "best_fit": {
          "Phi": 0.34906585039886573,
          "Delta": 2.4958208303518914
        },
        "solutions": [
          {
            "Phi": 0.3490658503988655,
            "Delta": 2.4958208303518914
          }
        ]
      }
    },
    "wrong_sign_model_rejected": {
      "passed": true,
      "detail": {
        "status": "model_mismatch",
        "phase_rms_rad": 0.2942516024924734
      }
    },
    "unknown_signs_not_guessed": {
      "passed": true,
      "detail": {
        "status": "missing_sign_calibration",
        "modes": [
          1,
          4,
          6
        ],
        "solutions": []
      }
    },
    "zero_modes_underdetermined": {
      "passed": true,
      "detail": {
        "status": "insufficient_modes",
        "modes": [],
        "solutions": []
      }
    },
    "one_mode_underdetermined": {
      "passed": true,
      "detail": {
        "status": "insufficient_modes",
        "modes": [
          4
        ],
        "solutions": []
      }
    },
    "general_state_not_forced_into_decoder": {
      "passed": true,
      "detail": {
        "status": "outside_calibrated_family",
        "solutions": []
      }
    },
    "evolution_works_when_decoder_rejects": {
      "value": 0.0,
      "upper_bound": 5e-12,
      "passed": true
    }
  },
  "all_passed": true
}
~~~
### nonlinear_results.json

~~~json
{
  "model": "finite_Galerkin_cubic_phase_mixing",
  "M": 6,
  "kappa": 1.0,
  "s_final": 0.5,
  "integrator": {
    "method": "DOP853",
    "rtol": 1e-12,
    "atol": 1e-14,
    "exact_norm_preservation_claimed": false,
    "exact_numerical_reversibility_claimed": false
  },
  "versions": {
    "numpy": "2.3.5",
    "scipy": "1.17.0"
  },
  "initial_norm_squared": 0.9999999999999998,
  "sampled_max_norm_squared_error": 2.278177646530821e-13,
  "origin_shift_radians": 0.371,
  "rhs_origin_covariance_max_error": 2.482534153247273e-16,
  "sampled_flow_origin_covariance_max_error": 7.28550816809042e-16,
  "forward_then_backward_max_error": 1.7019962483368879e-15,
  "initial_new_mode_derivatives": {
    "m0": [
      0.0,
      -0.3535533905932737
    ],
    "m3": [
      0.0,
      -0.3535533905932737
    ],
    "formula_max_error": 0.0
  },
  "initial_power_by_mode": {
    "-6": 0.0,
    "-5": 0.0,
    "-4": 0.0,
    "-3": 0.0,
    "-2": 0.0,
    "-1": 0.0,
    "0": 0.0,
    "1": 0.4999999999999999,
    "2": 0.4999999999999999,
    "3": 0.0,
    "4": 0.0,
    "5": 0.0,
    "6": 0.0
  },
  "final_power_by_mode": {
    "-6": 1.356834692643832e-14,
    "-5": 2.0524627744921834e-13,
    "-4": 3.0511381261660766e-11,
    "-3": 1.2910773222246196e-08,
    "-2": 3.2996012565088923e-06,
    "-1": 0.000471588590279995,
    "0": 0.029815306974815138,
    "1": 0.4697097918767958,
    "2": 0.469709792341948,
    "3": 0.02981530590533633,
    "4": 0.0004715895278059997,
    "5": 3.2991812034925007e-06,
    "6": 1.3059053477218138e-08
  },
  "function_evaluations": {
    "forward": 152,
    "shifted": 152,
    "backward": 110
  },
  "limits": [
    "The complex axis is supplied, not relationally derived.",
    "Spatial transport, force and collision are not established.",
    "No bandwidth-to-continuum convergence test is included."
  ]
}
~~~
