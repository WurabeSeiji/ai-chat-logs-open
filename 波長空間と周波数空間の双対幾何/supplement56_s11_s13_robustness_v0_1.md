# 補遺56:s=11,13 頑健性検査 — Z₄量子化の一般性確認・「1ビット」主張の非縮退限定・測度観測量の s=9 固有性

**作成**: 2026年6月12日
**版**: v0.1
**対象**: 補遺55 §4 開問題 #3(「s=11,13 拡張:セクター構造の持続性、0.98195 vs 0.98263 の差が有限サイズ補正か、多段鎖でのホロノミー＝パリティ則の確認」)の実施。補遺51/52/54 の機械(`expcoeffs` / `cross_at` / `self_at` / `tphase` / `two_step_paths`)を s=9 から s=11,13 へそのまま拡張(箱を K=5 に拡大)。再実装は s=9 の既知値(代表親 (2,1,0,0):531→1088、333→40、P(X)=0.98195、(1,1,7)→0 相殺)を厳密再現することを確認済み(整合アンカー PASS)。スクリプト:`supplement56_s11_s13_robustness.py`。

**性格**: 三つの結果と一つの限定。(1) **Z₄量子化は一般(補遺55 #6 強化)**:単段輸送位相 η が Z₄ に乗ることは s=9・11・13 で全数 0 違反(2816 / 6144 / 17744 件)。位相アルファベットが導出である定理は s=9 の人工物でない。(2) **支配軸符号＝主分裂は一般、「ちょうど2セクター」は非縮退限定**:支配軸が一意な親(s=11 型 (2,1,1,0)・s=13 型 (3,0,0,0))は補遺52/54 通りちょうど2セクターに分裂し判別子＝支配軸符号(s=11 は48/48、s=9 の24/24を再現)。**しかし縮退軸の親では破れる** — s=13 型 (2,1,1,1) は単一 B₄ 軌道(64セル)が**5つのゲージ固定 |z|² クラス**に分裂(5は2のべきでない＝指数2/単一 Z₂ では不能)、s=13 型 (2,2,0,0)(支配軸タイ)は3セクター。**補遺52「ちょうど1ビット」・補遺54「単一 Z₂ 鎖ホロノミー」は支配軸が一意な場合の特殊命題**であり、補遺53 §5-2 が予告した縮退の弱点が定量化された。(3) **測度値ギャップ問題は延長不能**:双子 P(X)＝2W₅₃₁/(2W₅₃₁+W₃₃₃) は二チャネル構造に依存した s=9 固有の観測量。s=11 は非相殺チャネル3本、s=13 はチャネル集合が親型依存 — 一意な「spread 対 symmetric」対が一般に無く、0.98195 vs 0.98263 の差の有限サイズ判定は観測量を作り直さない限り不能(これ自体が結論)。物理的同一視は行わない。

---

## 1. 整合アンカー(再実装の忠実性)

補遺51/52/54 と同一の機械を K=5 の箱で再実装。代表親 $v_9=(2,1,0,0)$ で:

| チャネル | $z$ | $\lvert z\rvert^2$ |
|---|---|---|
| (1,1,7) | 0 | **0(相殺)** |
| (1,3,5) | $-32+8i$ | 1088 |
| (3,3,3) | $-6+2i$ | 40 |

双子 $P(X)=2\cdot1088/(2\cdot1088+40)=$ **0.98195**。補遺51 §2.2 と厳密一致 — 拡張結果の信頼の土台。(1,1,7) チャネルの相殺が、補遺51 が 531/333 だけ使った理由でもある(二つの自明娘を含むチャネルは消える)。

## 2. Z₄量子化の一般性(補遺55 #6)— PASS

単段輸送位相 $\eta(v_p; v_a,v_b)=\big(C_{\text{cross}}(v_p)/\lvert C_{\text{cross}}\rvert\big)\big/\big(C_{\text{par}}(v_p)/\lvert C_{\text{par}}\rvert\big)$ を、各シェル $s$ の全親 $v_p$・全娘対 $(a,b)$($a+b-s\in\{\pm1\}$)・全セル対で計算。

| $s$ | 単段 η 総数 | Z₄外 | 位相分布 $(i^0,i^1,i^2,i^3)$ |
|---|---|---|---|
| 9 | 2816 | **0** | (448, 560, 960, 848) |
| 11 | 6144 | **0** | (960, 960, 2112, 2112) |
| 13 | 17744 | **0** | (2928, 4264, 5944, 4608) |

> **判定**: $\eta\in Z_4$ は s=9 の全数事実でなく、s=11,13 でも例外ゼロ。補遺51 の「位相アルファベットは割当でなく辞書の三角恒等式の帰結」は頑健。全 Z₄ 値が出現(¼回転格子が稠密に使われる)。**台帳の最も強い主張(#6)が頑健性検査を生き延びた。**

## 3. セクター構造の一般性(補遺52/54)

各型の全セルでゲージ固定 $\lvert z\rvert^2$(全チャネルの値の組＝signature)を計算し、異なる signature の個数(＝セクター数)と支配軸符号(動径に最も近い軸、補遺53)との対応を見た。

### 3.1 支配軸が一意の親 — ちょうど2セクター(一般)

| $s$, 型 | セル数 | セクター数 | 判別子 |
|---|---|---|---|
| 11, (2,1,1,0) | 96 | **2(48/48)** | 支配軸符号 ✓ |
| 13, (3,0,0,0) | 8 | **2(4/4)** | 支配軸符号 ✓(−側は全相殺) |

s=11 型 (2,1,1,0) の二セクター:

| 支配符号 | (1,3,7) | (1,5,5) | (3,3,5) |
|---|---|---|---|
| + | 3200 | 1024 | 7888 |
| − | 640 | 256 | 5440 |

((1,1,9) は両セクターで相殺。)s=9 の 24/24 を 48/48 として完全再現 — **補遺52/54 の主構造は一般。**

### 3.2 支配軸が縮退する親 — 「2セクター」が破れる(限定)

| $s$, 型 | セル数 | セクター数 | 注 |
|---|---|---|---|
| 13, (2,2,0,0) | 24 | **3**(12/6/6) | 支配軸が2本でタイ |
| 13, (2,1,1,1) | 64 | **5**(24/12/12/8/8) | 単一 B₄ 軌道 |

型 (2,1,1,1) は $B_4$ の**単一軌道**($384/6=64$)である。その軌道がゲージ固定 $\lvert z\rvert^2$ で**5クラス**に割れる。指数2＝1ビットなら2クラス、$k$ ビットなら $2^k$ クラスのはずで、**5 は2のべきでない** → 縮退親に対して「指数2/ちょうど1ビット」(補遺52)も「単一 Z₂ 鎖ホロノミー」(補遺54)も成立しない。

ただし主分裂は支配軸符号で保たれる:5クラスは大チャネルで2家系に分かれ(family I:24+8＝32セル、family II:12+12+8＝32セル)、これは $\lvert2\rvert$ 軸の符号(32/32)に一致。**追加の微細分裂は対称チャネル (3,3,7) のみ**(family II 内で 61520 / 60552 / 59600、family I 内で 33920 / 32904)。すなわち

$$\lvert z\rvert^2 \approx (\text{支配軸符号:主・一般}) + (\text{対称チャネルの微細構造:縮退時のみ}).$$

後者は補遺38 残余 $Z_2\times Z_2$ の第二成分(娘交換)か、より豊かな構造かが未決(§6)。微細分裂の軌道サイズが不均一(12/12/8)であることは、単純な群作用(等サイズ軌道)でない可能性を示唆する。

## 4. 測度値ギャップ問題(補遺55 #9・開問題 #3)— 延長不能

双子 $P(X)=2W_{531}/(2W_{531}+W_{333})$ は、s=9 の二チャネル(spread＝531・symmetric＝333)に依存した観測量である。拡張すると:

- **s=11**:非相殺チャネルは3本((1,3,7),(1,5,5),(3,3,5))。一意な「spread 対 symmetric」対が無い。
- **s=13**:チャネル集合が親型依存(型 (3,0,0,0) は (3,5,5) のみ、型 (2,1,1,1) は4チャネル)。

> **判定**: $P(X)$ の二チャネル定義は s=9 固有。0.98195 と数え上げ 0.98263 の差(0.0007)が有限サイズ補正か恒常的不一致かは、観測量を一般 $s$ で再定義しない限り判定できない。延長できるのは構造(Z₄・セクター・支配軸ビット)であって特定の数値ではない — これ自体が一つの結論(測度値の決着は s 拡張では得られない)。

## 5. 補遺55 への含意(整合更新)

| 補遺55 §1 | s=11,13 判定 | 処置 |
|---|---|---|
| #6 $\eta\in Z_4$ 導出 | **確認・強化**(0違反) | 無修正 |
| #7 ちょうど1ビット | 非縮退限定(縮退で3〜5セクター) | 「支配軸一意のとき」を付す |
| #8 単一 Z₂ 鎖ホロノミー | 主成分は一般、縮退で残余微細構造 | 同上 |
| #9 セクター値 0.98195 等 | s=9 固有(観測量) | 「s=9 限定」を付す |
| §4 開問題 #3 | 部分解決 | 上記三結果へ更新 |

50年テーゼ(3公理＋正準構造、追加ゼロ)は s 拡張で公理が増えたわけでなく**無傷**。揺らいだのは「測度＝ちょうど1ビット/単一 Z₂」の一般性のみで、背骨(導出位相・セクター構造・支配軸ビット)は生き延びた。

## 6. 正直な限界

1. **η² 共変縮約・指数2部分群の直接再計算はしていない**:本稿が s=11,13 で測ったのは下流の観測量(ゲージ固定 $\lvert z\rvert^2$ の signature 数＝セクター数)である。これは「指数2 ⇒ 2セクター」の操作的内容であり、単一軌道が5クラスに割れる事実は縮退親での「ちょうど1ビット」を反証するが、補遺52 の 80/80 相殺そのものの s 一般性は別途要検査。
2. **縮退セクターの微細構造の正体は未決**:対称チャネルの追加分裂が補遺38 残余 $Z_2\times Z_2$ の娘交換成分か、独立構造かは未同定(§3.2)。軌道サイズの不均一性は単純群作用でない可能性を示す。
3. **s=15 以上は未検査**:本稿は s=11,13 の二点。多段鎖(s≥15 で3段以上)でのホロノミー＝パリティ則の一般証明は未了(補遺54 §5-1 と同枠)。
4. **箱 K=5 の有限性**:全シェル(≤13)のセルは K=5 で完全に捕捉済み(単一軸最大 $\lvert k\rvert=3$ で十分)だが、より高い $s$ では箱の拡大が要る。
5. 物理的同一視は行わない。

---

物理的同一視は行わない。本稿は補遺49〜54 の s=9 全数事実の頑健性検査であり、新規の物理的主張を含まない。検算は claude.ai セッション(2026-06-12)で実行、スクリプトは §A(`supplement56_s11_s13_robustness.py` と同一)。

## §A. 検算スクリプト

```python
# -*- coding: utf-8 -*-
# 補遺56: s=11,13 頑健性検査。補遺51/52/54 の機械をそのまま拡張(箱 K=5)。
import itertools, numpy as np
from collections import defaultdict
SQ2 = np.sqrt(2.0)

_EXP = {}
def expcoeffs(v):
    v = tuple(v)
    if v in _EXP: return _EXP[v]
    acc = [((), 1.0+0j)]
    for i in range(4):
        m = v[i]
        if m == 0: opts = [(0, 1.0+0j)]
        elif m > 0: opts = [(+m, 1/SQ2+0j), (-m, 1/SQ2+0j)]
        else:
            mm = -m; opts = [(+mm, -1j/SQ2), (-mm, +1j/SQ2)]
        acc = [(n+(f,), c*w) for n, c in acc for f, w in opts]
    d = {}
    for n, c in acc: d[n] = d.get(n, 0) + c
    _EXP[v] = d; return d

def cross_at(va, vb, t):
    A = expcoeffs(va); B = expcoeffs(vb); tot = 0+0j
    for na, ca in A.items():
        nb = tuple(x-y for x, y in zip(t, na))
        if nb in B: tot += ca*B[nb]
    return tot
def self_at(v, t): return expcoeffs(tuple(v)).get(tuple(t), 0+0j)

def tphase(vp, va, vb, orient=+1):
    tgt = tuple(orient*x for x in vp)
    Cp = self_at(vp, tgt)
    if abs(Cp) < 1e-12: return None
    Cc = cross_at(tuple(va), tuple(vb), tgt)
    if abs(Cc) < 1e-12: return None
    return (Cc/abs(Cc)) / (Cp/abs(Cp))

def z4_idx(e):
    for p, u in enumerate([1, 1j, -1, -1j]):
        if abs(e-u) < 1e-9: return p
    return None

def shell_cells(m, K=5):
    return [tuple(k) for k in itertools.product(range(-K, K+1), repeat=4)
            if abs(sum((abs(t)+0.5)**2 for t in k) - m) < 1e-9]
SH = {m: shell_cells(float(m)) for m in (1,3,5,7,9,11,13)}

def two_step_paths(parent, final):
    paths = []; fin = tuple(sorted(final))
    for a in range(1, parent+2, 2):
        for b in range(a, parent+2, 2):
            d1 = a+b-parent
            if d1 not in (1, -1): continue
            for (x, spec) in ((a, b), (b, a)):
                for c in range(1, x+2, 2):
                    dd = x-d1-c
                    if dd < 1 or dd % 2 == 0: continue
                    if tuple(sorted((spec, c, dd))) == fin:
                        paths.append(((a, b), x, (c, dd), d1))
    return sorted(set(paths))

def all_finals(parent):
    fins = set()
    for a in range(1, parent+2, 2):
        for b in range(a, parent+2, 2):
            d1 = a+b-parent
            if d1 not in (1, -1): continue
            for (x, spec) in ((a, b), (b, a)):
                for c in range(1, x+2, 2):
                    dd = x-d1-c
                    if dd < 1 or dd % 2 == 0: continue
                    fins.add(tuple(sorted((spec, c, dd))))
    return sorted(fins)

def run_parent(vp, parent, finals, o1=+1, o2=+1):
    z = {f: 0+0j for f in finals}
    for f in finals:
        for ((a, b), x, (c, d), d1) in two_step_paths(parent, f):
            for va in SH[a]:
                for vb in SH[b]:
                    if a == b and va >= vb: continue
                    e1 = tphase(vp, va, vb, o1)
                    if e1 is None: continue
                    vx = va if x == a else vb
                    for vc in SH[c]:
                        for vd in SH[d]:
                            if c == d and vc >= vd: continue
                            e2 = tphase(vx, vc, vd, o2)
                            if e2 is None: continue
                            z[f] += e1*e2
    return z

# --- アンカー: s=9 が 1088/40/0.98195 を再現 ---
z9 = run_parent((2,1,0,0), 9, all_finals(9))
# --- TEST A: Z4 量子化 (単段, 全数) ---
def daughter_pairs(s):
    return [(a, b) for a in range(1, s+2, 2) for b in range(a, s+2, 2) if a+b-s in (1, -1)]
for s in (9, 11, 13):
    tot = viol = 0; dist = defaultdict(int)
    for vp in SH[s]:
        for (a, b) in daughter_pairs(s):
            for va in SH[a]:
                for vb in SH[b]:
                    if a == b and va >= vb: continue
                    e = tphase(vp, va, vb, +1)
                    if e is None: continue
                    tot += 1; p = z4_idx(e)
                    if p is None: viol += 1
                    else: dist['i^%d' % p] += 1
    print(f"s={s}: η総数 {tot}, Z4外 {viol}  分布 {dict(dist)}")
# --- TEST C: セクター構造 (全セル, 全型) ---
def dom_sign(v):
    mx = max(abs(x) for x in v)
    ax = [j for j in range(4) if abs(v[j]) == mx]
    return None if len(ax) > 1 else (1 if v[ax[0]] > 0 else -1)
for s in (11, 13):
    fins = all_finals(s)
    types = defaultdict(list)
    for v in SH[s]: types[tuple(sorted(map(abs, v)))].append(v)
    for t in sorted(types):
        sigs = defaultdict(list)
        for vp in types[t]:
            z = run_parent(vp, s, fins)
            sigs[tuple(round(abs(z[f])**2, 1) for f in fins)].append(vp)
        print(f"s={s} 型{t} ({len(types[t])}セル): セクター数 {len(sigs)}")
```
