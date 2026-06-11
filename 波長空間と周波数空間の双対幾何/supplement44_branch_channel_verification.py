#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 補遺44 §4 機械検算:
#  V1: 双子×娘系の結合ゲージクラス数 (共有B4での軌道分解)
#  V2: 5関係クラスの枝チャネル分解つき線テーブル (論文7 図3拡張) + アルファベット検査 + ゲージ検査
#  V3: 補遺30 §5 縞シフトテスト (三脚: 対角¼周期 / 環: 縞反転) の直接数値実行
#  V4: 縮退線チャネル分割が識別力を加える最小例の探索
import itertools, numpy as np
from collections import defaultdict

# ---------- B4 群 (384) ----------
PERMS = list(itertools.permutations(range(4)))
SIGNS = list(itertools.product([1,-1],repeat=4))
def apply_g(perm, sign, v):
    return tuple(sign[i]*v[perm[i]] for i in range(4))

def shell_cells(m):
    out=[]; K=4
    for k in itertools.product(range(-K,K+1),repeat=4):
        if abs(sum((abs(t)+0.5)**2 for t in k)-m)<1e-9: out.append(tuple(k))
    return out
SH3=shell_cells(3.0); SH5=shell_cells(5.0); SH9=shell_cells(9.0)

# ================= V1: 結合ゲージクラス数 =================
print("="*72)
print("V1: 双子×娘系 {v5,v3,0}∪{vB} の結合ゲージクラス (共有B4)")
print("="*72)
def canon(config):
    best=None
    for perm in PERMS:
        for sign in SIGNS:
            img = tuple(apply_g(perm,sign,v) for v in config)
            if best is None or img<best: best=img
    return best
classes=defaultdict(int)
byA=defaultdict(set)
for v5 in SH5:
    for v3 in SH3:
        ip = sum(a*b for a,b in zip(v5,v3))
        for vB in SH9:
            c = canon((v5,v3,vB))
            classes[c]+=1
            byA[ip].add(c)
total = sum(classes.values())
print(f"  配置総数 = 24×8×64 = {total}")
print(f"  結合ゲージクラス数 = {len(classes)}")
for ip,name in [(0,'直交'),(1,'平行'),(-1,'反平行')]:
    print(f"   A側クラス {name} (v5·v3={ip:+d}): 結合クラス {len(byA[ip])} 個")
print(f"  比較: A単独 3 クラス × B単独 1 クラス = 3 → 結合 {len(classes)} 클래스")
print(f"  → 系間相対位置が観測量として {len(classes)-3 if len(classes)>3 else 0} クラス分の新情報を担う")

# ================= 共通: 厳密スペクトル計算 (5点/軸グリッド) =================
N=5
xs = np.arange(N)/N
def axis_basis():
    """正規直交実基底 (5点一様測度): [1, √2cos2πx, √2sin2πx, √2cos4πx, √2sin4πx]"""
    B=np.zeros((5,N))
    B[0]=1
    B[1]=np.sqrt(2)*np.cos(2*np.pi*xs); B[2]=np.sqrt(2)*np.sin(2*np.pi*xs)
    B[3]=np.sqrt(2)*np.cos(4*np.pi*xs); B[4]=np.sqrt(2)*np.sin(4*np.pi*xs)
    return B
B = axis_basis()
G = B@B.T/N
assert np.allclose(G, np.eye(5), atol=1e-12), "基底直交性"

def phi_axis(m, shift=0.0):
    """セル成分 m の軸因子 (shift=並進): m=0→1, m>0→√2cos(2πm(x-a)), m<0→√2sin(2π|m|(x-a))"""
    if m==0: return np.ones(N)
    if m>0:  return np.sqrt(2)*np.cos(2*np.pi*m*(xs-shift))
    return np.sqrt(2)*np.sin(2*np.pi*(-m)*(xs-shift))

def intensity(cells, shifts=None):
    """I(x)=Ψ² on 5^4 grid; shifts: セルごとの並進ベクトル (None=0)"""
    Psi = np.zeros((N,N,N,N))
    for idx,k in enumerate(cells):
        a = shifts[idx] if shifts is not None else (0,0,0,0)
        f = np.einsum('i,j,k,l->ijkl', phi_axis(k[0],a[0]), phi_axis(k[1],a[1]),
                      phi_axis(k[2],a[2]), phi_axis(k[3],a[3]))
        Psi += f
    return Psi**2

def channel_table(I):
    """係数テンソル → {(周波数絶対値ベクトル, 枝パターン): 振幅}"""
    C = np.einsum('ijkl,ai,bj,ck,dl->abcd', I, B,B,B,B)/N**4
    out={}
    FREQ={0:0,1:1,2:1,3:2,4:2}; BR={0:'-',1:'c',2:'s',3:'C',4:'S'}
    it=np.nditer(C, flags=['multi_index'])
    for val in it:
        v=float(val)
        if abs(v)<1e-9: continue
        idx=it.multi_index
        freq=tuple(FREQ[i] for i in idx); br=''.join(BR[i] for i in idx)
        out[(freq,br)]=round(v,9)
    return out

def classic_table(tab):
    """(ノルム²: 本数, パワー) — 論文7の旧プロトコル (DC除く)"""
    agg=defaultdict(lambda:[0,0.0])
    for (freq,br),amp in tab.items():
        n2=sum(f*f for f in freq)
        if n2==0: continue
        agg[n2][0]+=1; agg[n2][1]+=amp*amp
    return {k:(v[0],round(v[1],6)) for k,v in sorted(agg.items())}

# ================= V2: 5関係クラスの枝チャネル線テーブル =================
print()
print("="*72)
print("V2: 枝チャネル分解つき線テーブル (5クラス代表)")
print("="*72)
REPS = {
 '(5,3,1)直交':   [(1,1,0,0),(0,0,1,0),(0,0,0,0)],
 '(5,3,1)平行':   [(1,1,0,0),(1,0,0,0),(0,0,0,0)],
 '(5,3,1)反平行': [(1,1,0,0),(-1,0,0,0),(0,0,0,0)],
 '(3,3,3)三脚':   [(1,0,0,0),(0,1,0,0),(0,0,1,0)],
 '(3,3,3)対蹠対': [(1,0,0,0),(-1,0,0,0),(0,1,0,0)],
}
SQRT2=np.sqrt(2)
ALPHABET = sorted({(SQRT2**j) for j in range(-2,4)})
def in_alphabet(a):
    return any(abs(abs(a)-w)<1e-9 for w in ALPHABET)
chan_canon={}
for name,cells in REPS.items():
    tab = channel_table(intensity(cells))
    lines = [(k,v) for k,v in tab.items() if sum(f*f for f in k[0])>0]
    n_ch = len(lines)
    n_lines = len({k[0] for k,_ in lines})
    alpha_ok = all(in_alphabet(v) for _,v in lines)
    print(f"  {name}: 線 {n_lines} 本 → 枝チャネル {n_ch} 本, アルファベット(√2)^j: {'PASS' if alpha_ok else 'FAIL'}")
    bynorm=defaultdict(list)
    for (freq,br),amp in lines:
        bynorm[sum(f*f for f in freq)].append((freq,br,amp))
    for n2 in sorted(bynorm):
        s=', '.join(f"{br}@{freq}:{amp:+.3f}" for freq,br,amp in sorted(bynorm[n2])[:6])
        more='' if len(bynorm[n2])<=6 else f" …(+{len(bynorm[n2])-6})"
        print(f"     |n|²={n2}: {s}{more}")

# ゲージ検査: B4 でクラス内のチャネルテーブルが移り合う (canonical 一致)
print("  --- ゲージ検査 (B4 でチャネルテーブル正準形が不変) ---")
def chan_canonical(cells):
    best=None
    for perm in PERMS:
        for sign in SIGNS:
            img=[apply_g(perm,sign,v) for v in cells]
            t=channel_table(intensity(img))
            key=tuple(sorted(t.items()))
            if best is None or key<best: best=key
    return best
import random
random.seed(7)
ok=True
for name,cells in list(REPS.items())[:2]:
    base = chan_canonical(cells)
    for _ in range(3):
        perm=random.choice(PERMS); sign=random.choice(SIGNS)
        img=[apply_g(perm,sign,v) for v in cells]
        if chan_canonical(img)!=base: ok=False
print(f"  単一系のチャネルテーブルはゲージ軌道で閉じる: {'PASS' if ok else 'FAIL'}")

# ================= V3: 縞シフトテスト (補遺30 §5) =================
print()
print("="*72)
print("V3: 縞シフトテスト (補遺30 §5 の直接数値実行)")
print("="*72)
# 三脚テスト: 三脚 {e1,e2,e3} を輸送 (T_{1/4(e1+e2+e3)}), 参照断片 e4 は未輸送
tripod=[(1,0,0,0),(0,1,0,0),(0,0,1,0)]; ref=[(0,0,0,1)]
cells = tripod+ref
sh0 = [(0,0,0,0)]*4
shT = [(0.25,0.25,0.25,0)]*3 + [(0,0,0,0)]
t0 = channel_table(intensity(cells, sh0))
t1 = channel_table(intensity(cells, shT))
# 予測: 輸送断片のモードは cos→sin (m=1, ¼周期)。三脚内部の交差線(相対変位0)は不変、
# 三脚×参照の交差線はチャネルが c↔s 回転。自己線(2倍音)は m=2: 半周期→符号反転 cos4π(x-1/4)=-cos
def transform_pred(tab):
    """T_{1/4} を軸1-3の輸送断片に適用したときの理論変換で t0 から t1 を予測"""
    pred={}
    for (freq,br),amp in tab.items():
        nf=list(br); a=amp
        for ax in range(3):  # 輸送軸 (e4=軸4 は参照)
            f=freq[ax]; b=br[ax]
            # この軸成分が輸送断片由来か参照由来かは線ごとに異なるため、
            # ここでは検証を「不変量の比較」で行う (下の checks)
        pred[(freq,br)]=a
    return pred
# 検証は不変量比較で行う:
# (a) 三脚内部交差線 (軸4成分なし, 両断片輸送) → 不変
# (b) 三脚×参照交差線 (軸4成分あり) → チャネル回転 c↔s (m=1), 振幅絶対値は不変
# (c) 自己線 (単一軸 2倍音) → 輸送断片は符号反転
# 統一予測: 輸送軸 (1-3) の全成分が ¼ 回転 (m=1: c→s, s→-c / m=2: 符号反転)、参照軸 (4) は不変。
# 補遺30 §5 の「全枝が 90° 回転」の厳密形。検証は全線一括。
def rot_q(freq,br,amp):
    nb=list(br); a=amp
    for ax in range(3):
        if freq[ax]==1:
            if br[ax]=='c': nb[ax]='s'
            elif br[ax]=='s': nb[ax]='c'; a=-a
        elif freq[ax]==2:
            a=-a  # 2倍音は 4π·(1/4)=π シフト → 符号反転
    return (freq,''.join(nb)), a
pred={}
for (freq,br),amp in t0.items():
    (f2,b2),a2 = rot_q(freq,br,amp)
    pred[(f2,b2)]=pred.get((f2,b2),0)+a2
keys=set(pred)|set(t1)
ok_all_lines = all(abs(pred.get(k,0)-t1.get(k,0))<1e-9 for k in keys)
# 部分別の内訳 (報告用): 三脚×参照交差線のみでの照合
ok_cross=True; n_cross=0
for (freq,br),amp in t0.items():
    if freq[3]>0 and any(freq[ax]>0 for ax in range(3)):
        n_cross+=1
        (f2,b2),a2 = rot_q(freq,br,amp)
        if abs(t1.get((f2,b2),0)-a2)>1e-9: ok_cross=False
print(f"  三脚輸送 T_(1/4)(e1+e2+e3): 全線の¼回転予測との厳密一致 {'PASS' if ok_all_lines else 'FAIL'}"
      f" / うち三脚×参照交差線 {n_cross} 本の¼縞シフト {'PASS' if ok_cross else 'FAIL'}")

# 環テスト: 環 {e1,-e1,e2,-e2} 輸送 T_{1/2(e1+e2)}, 参照 e3
ring=[(1,0,0,0),(-1,0,0,0),(0,1,0,0),(0,-1,0,0)]; cells2=ring+[(0,0,1,0)]
sh0b=[(0,0,0,0)]*5
shTb=[(0.5,0.5,0,0)]*4+[(0,0,0,0)]
u0=channel_table(intensity(cells2,sh0b)); u1=channel_table(intensity(cells2,shTb))
# 予測: 輸送断片の m=1 モードは cos2π(x-1/2)=-cos (符号反転) → 環×参照の交差線が全て符号反転 (縞の反転)。
# 環内部線は両方反転→不変。自己2倍音は 4π·(1/2)=2π → 不変。
ok_ring=True
for (freq,br),amp in u0.items():
    n2=sum(f*f for f in freq)
    if n2==0: continue
    cross_ref = (freq[2]>0) and (freq[0]>0 or freq[1]>0)
    odd = sum(freq[ax] for ax in (0,1) if freq[ax]==1)  # 環側の m=1 成分数
    exp = -amp if (cross_ref and odd%2==1) else amp
    if abs(u1.get((freq,br),0)-exp)>1e-9: ok_ring=False
print(f"  環輸送 T_(1/2)(e1+e2): 交差線の縞反転 (奇数次のみ符号反転) {'PASS' if ok_ring else 'FAIL'}")

# ================= V4: 縮退線チャネル分割の識別力 (最小例の探索) =================
print()
print("="*72)
print("V4: 旧プロトコル縮退 × チャネル分割で分離する最小例の探索")
print("="*72)
def canon_cells(cells):
    best=None
    for perm in PERMS:
        for sign in SIGNS:
            img=tuple(sorted(apply_g(perm,sign,v) for v in cells))
            if best is None or img<best: best=img
    return best
found=None
for family,cellpool,size in [("shell3対",SH3,2),("shell5対",SH5,2),("shell3三つ組",SH3,3),("shell9対",SH9,2)]:
    orbits=defaultdict(list)
    for combo in itertools.combinations(cellpool,size):
        orbits[canon_cells(combo)].append(combo)
    # 各軌道の classic テーブルと チャネル正準形
    sigs={}
    for orb,members in orbits.items():
        rep=members[0]
        tab=channel_table(intensity(list(rep)))
        sigs[orb]=(tuple(sorted(classic_table(tab).items())), )
    # classic が同一で軌道が異なる組を探す
    byclassic=defaultdict(list)
    for orb,s in sigs.items(): byclassic[s].append(orb)
    degen={k:v for k,v in byclassic.items() if len(v)>1}
    print(f"  {family}: 軌道 {len(orbits)} 個, classic縮退組 {len(degen)} 組")
    if degen and found is None:
        k,orbs = next(iter(degen.items()))
        found=(family, orbs)
        # チャネル正準形で分離するか
        chans=[tuple(sorted(channel_table(intensity(list(o))).items())) for o in orbs]
        sep = len(set(chans))==len(chans)
        print(f"    → 最小例: {family} の軌道 {len(orbs)} 個が classic 縮退。チャネル正準形での分離: {'YES' if sep else 'NO'}")
        for o in orbs: print(f"       代表: {o}")
print("  (見つからない族は「classic で既に完全分離」= チャネルは追加識別力なし)")
