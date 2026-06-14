# -*- coding: utf-8 -*-
"""
論文0.5「変位記録の仕組みの考察」検算スクリプト（区間演算版・査読1反映）
- ν 振動（奇数倍音＝論理波/矩形波）を参照に、観測量 x（離散変位）を読む。
- x 自身が ±δ(=±1/2) の幅を持つ（A2 の零点）。サンプリングも ±1/2。
- 合成は「区間の和」 [-1/2,1/2]+[-1/2,1/2]=[-1,1] ＝ ±1（分布形は問わない・追加公理ゼロ）。
- 両者は記録上 区別不能（2自由度→1自由度への射影、核は1次元）。
- ν 自身の揺らぎは省略（厳密とする）。参考: 三者なら 3×[-1/2,1/2]=[-3/2,3/2]。
出力: 図(2段) と 区間演算の検算。乱数・分布仮定は用いない。
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

for cand in ["Hiragino Sans","Hiragino Kaku Gothic ProN","YuGothic","Yu Gothic"]:
    try:
        font_manager.findfont(cand, fallback_to_default=False); plt.rcParams["font.family"]=cand; break
    except Exception: pass
plt.rcParams["axes.unicode_minus"]=False

# ---- 設定 ----
N=10
jumps=[2.3,5.7,8.4]; levels=[0,1,2,3]
delta=0.5                      # x 自身の曖昧さ ±δ（A2 の零点、δ=1/2）
phi=np.linspace(0,N,20001)     # 横軸＝巻き数（振動回数）。時間ではない。

def x_true(p):
    y=np.full_like(p,levels[0],dtype=float)
    for j,lv in zip(jumps,levels[1:]): y[p>=j]=lv
    return y

def nu_logic(p,K=15):          # 奇数倍音和 → 矩形波（論理波）。引数は巻き数 φ/2π。
    s=np.zeros_like(p)
    for k in range(1,K+1,2): s+=np.sin(2*np.pi*k*p)/k
    return (4/np.pi)*s

xt=x_true(phi); nu=nu_logic(phi)
n_samples=np.arange(0,N+1); x_rec=x_true(n_samples.astype(float))

# ---- 区間演算ヘルパ（厳密・乱数なし）----
def isum(*ints):
    lo=sum(a for a,_ in ints); hi=sum(b for _,b in ints); return (lo,hi)

Ix=(-delta,delta)       # 元値の曖昧
Is=(-0.5,0.5)           # サンプリングの曖昧
Inu=(-0.5,0.5)          # （参考）ν の曖昧

print("=== 検算1: サンプリングによる変化位置の ±1/2（区間） ===")
print("  注: ±1/2 はセル幅1の定義そのもの。中央寄り例は上界に触れない。境界例で上界への漸近を確認する。")
for j in list(jumps)+[2.999,3.001]:        # 末尾2つは境界例（指摘4）
    nf=int(np.ceil(j)); est=nf-0.5
    tag=" (境界)" if j in (2.999,3.001) else ""
    print(f"  真={j:.3f}{tag}  記録区間=({nf-1},{nf}] 幅1  推定={est:.2f}±0.5  誤差={est-j:+.3f}  真は区間内:{nf-1< j <=nf}")

print("\n=== 検算2: 合成（区間演算）===")
I2=isum(Ix,Is)
print(f"  (a) 上界(無公理): |εx+εs|<=1/2+1/2=1 → supp(ε) ⊆ {I2}  （劣加法性のみ・分布も独立も不要）")
print(f"  (b) 等号(端点到達性の下): 端点(±1/2,±1/2)に結合支持が届く ⇔ 両源が同符号方向に同時±1/2可 → supp(ε) = {I2} = ±1")
print("      独立は十分だが必要でない（完全正相関でも到達）。壊すのは完全反相関(εx=−εs→supp={0})のみ。")
print("  ※ 一様＋独立を仮定すれば三角分布・σ=1/√6 だが、A1–A3 に無い確率仮定なので本稿は用いない。")

print("\n=== 検算3: 区別不能性（情報・射影）===")
print("  記録は (εx, εs) の2自由度を ε=εx+εs の1自由度へ射影。核は1次元(εx−εs 方向)。")
print("  → 記録から (εx,εs) は復元不能＝サンプリング誤差と元値の曖昧は区別不能（分布に依らず区間で言える）。")

print("\n=== 検算4(参考): 三者すべて ±1/2 なら（区間）===")
I3=isum(Ix,Is,Inu)
print(f"  {Ix} ⊕ {Is} ⊕ {Inu} = {I3} → ±3/2（±1 ではない）。本図は ν を厳密として省略。")

# ---- 図(2段) ----
fig,(ax1,ax2)=plt.subplots(2,1,figsize=(9,6.6),sharex=True)

# 図1: 真の変位 x（±δ のグレー幅・方眼に乗らない変化）
ax1.fill_between(phi, xt-delta, xt+delta, step="post", color="0.8", label=f"x の曖昧 ±δ(=±{delta})")
ax1.step(phi,xt,where="post",color="#1f4e79",lw=2,label="x の公称値")
for n in n_samples: ax1.axvline(n,color="0.9",lw=0.8,zorder=0)
for j in jumps: ax1.plot(j,x_true(np.array([j-1e-6]))[0],"o",color="crimson",ms=5,zorder=5)
ax1.set_ylabel("変位 x"); ax1.set_yticks(levels)
ax1.set_title("図1  真の変位 x：巻き数の方眼(灰)に乗らない位置で階段変化＋ x 自身の ±δ 曖昧（グレー幅）",fontsize=9)
ax1.legend(fontsize=8,loc="upper left")

# 図2: ν論理波でサンプリングした記録（合成 ±1 の区間グレー幅）
ax2.plot(phi,0.35*nu+0.5,color="0.8",lw=0.8,label="ν 振動（奇数倍音→論理波 ON/OFF）")
xr=np.append(x_rec,x_rec[-1]); ns=np.append(n_samples,N)
ax2.fill_between(ns, xr-1.0, xr+1.0, step="post", color="0.8",
                 label="記録の合成曖昧 ±1 ＝ [−½,½]⊕[−½,½]（区間演算）")
ax2.step(ns,xr,where="post",color="#2a7a2a",lw=2,label="変位記録（公称）")
ax2.plot(n_samples,x_rec,"s",color="#2a7a2a",ms=4)
ax2.set_ylim(-1.6,4.2); ax2.set_yticks(levels)
ax2.set_xlabel("巻き数（振動回数）0〜10  ※時間ではない"); ax2.set_ylabel("変位記録")
ax2.set_title("図2  記録：サンプリング±½ と 元値±½ は区別不能 → 合成 <=±1(無公理)/=±1(端点到達)",fontsize=9)
ax2.legend(fontsize=7.5,loc="upper left")

plt.tight_layout()
out="paper0_5_fig_displacement_record.png"
plt.savefig(out,dpi=150)
print(f"\n図を保存: {out}")
