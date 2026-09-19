#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Generate the five schematic figures used in
paper_complex_square_closure_flux_condensate_readout_ja_v1_3.md.

The original figures were produced with an image-generation model, not by code.
This script is a reproducible vector-style reimplementation of the same content.

Outputs (PNG):
  1. 複素数の閉性と閉路残差.png
  2. 三体の三角閉路とベクトル関係.png
  3. 四頂点ネットワークの三角閉路比較.png
  4. 内部相殺と凝縮体の階層的成長.png
  5. 図5_思考実験の方法フローチャート.png
"""

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Rectangle, FancyBboxPatch

OUT = Path(__file__).resolve().parent

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": [
        "Noto Sans CJK JP", "Noto Sans JP", "Yu Gothic", "Hiragino Sans", "DejaVu Sans"
    ],
    "mathtext.fontset": "dejavusans",
    "axes.unicode_minus": False,
})

BLUE = "#4C78A8"
ORANGE = "#F28E2B"
RED = "#D62728"
GREEN = "#59A14F"
PURPLE = "#8F63A8"
PALE_BLUE = "#EAF2F8"
PALE_ORANGE = "#FDEBD0"
PALE_GREEN = "#E9F7EF"
PALE_PINK = "#FCE4EC"
PALE_YELLOW = "#FFF4CC"
GRAY = "#777777"


def clean_ax(ax, xlim=(-0.2, 3.2), ylim=(-0.2, 2.8), axes=True):
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_aspect("equal", adjustable="box")
    ax.axis("off")
    if axes:
        ax.annotate("", xy=(xlim[1]-0.1, 0), xytext=(0, 0), arrowprops=dict(arrowstyle="->", lw=1.2))
        ax.annotate("", xy=(0, ylim[1]-0.1), xytext=(0, 0), arrowprops=dict(arrowstyle="->", lw=1.2))
        ax.text(xlim[1]-0.05, -0.12, "Re", ha="right", va="top", fontsize=10)
        ax.text(-0.08, ylim[1]-0.02, "Im", ha="right", va="top", fontsize=10)
        ax.text(-0.08, -0.08, "0", fontsize=9)


def vec(ax, p, q, label, color=BLUE, label_offset=(0, 0), lw=2.5):
    ax.annotate("", xy=q, xytext=p, arrowprops=dict(arrowstyle="->", lw=lw, color=color))
    m=((p[0]+q[0])/2,(p[1]+q[1])/2)
    ax.text(m[0]+label_offset[0], m[1]+label_offset[1], label, color=color, fontsize=11)


def save(fig, name):
    path = OUT / name
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    print(path)


def fig1():
    fig = plt.figure(figsize=(12, 8))
    gs = fig.add_gridspec(2, 3, height_ratios=[12,1.2], hspace=0.15, wspace=0.28)
    fig.suptitle("図1　複素数の閉性から閉路残差へ", fontsize=20, y=0.98)

    ax = fig.add_subplot(gs[0,0]); clean_ax(ax)
    ax.set_title("A. 複素数の閉性", fontsize=13, pad=8)
    vec(ax,(0,0),(1.55,1.35),"z",BLUE,(0.05,0.08))
    ax.scatter([1.55],[1.35],s=20,c="black")
    ax.text(0.15,-0.45,r"$z\in\mathbb{C}\Rightarrow z^2\in\mathbb{C}$",fontsize=12)
    ax.text(0.75,-0.82,r"$Z=z^2$",fontsize=12)

    ax = fig.add_subplot(gs[0,1]); clean_ax(ax,xlim=(-0.2,3.8),ylim=(-0.2,3.0))
    ax.set_title("B. 複素関係の和", fontsize=13, pad=8)
    p0=(0,0);p1=(1.15,0.75);p2=(1.9,1.75);p3=(3.0,2.15)
    vec(ax,p0,p1,r"$Z_1$",BLUE,(0.02,0.1)); vec(ax,p1,p2,r"$Z_2$",ORANGE,(0.04,0.08)); vec(ax,p2,p3,r"$Z_3$",BLUE,(0.03,0.08))
    for p in [p1,p2,p3]: ax.scatter([p[0]],[p[1]],s=20,c="black")
    ax.text(0.2,-0.45,r"$\sum Z_n=K$",fontsize=13)
    ax.text(0.2,-0.8,r"$Z_n:=z_n^2$",fontsize=11)
    ax.text(1.9,-0.8,r"$K\in\mathbb{C}$",fontsize=11)

    ax = fig.add_subplot(gs[0,2]); clean_ax(ax,xlim=(-0.2,4.2),ylim=(-0.2,2.9))
    ax.set_title("C. 閉路と残差", fontsize=13, pad=8)
    # closed triangle
    p0=(0.1,0.1); p1=(0.9,1.65); p2=(1.65,0.3)
    vec(ax,p0,p1,r"$Z_1$",BLUE,(-0.18,0.05),2.1); vec(ax,p1,p2,r"$Z_2$",ORANGE,(0.04,0.02),2.1); vec(ax,p2,p0,r"$Z_3$",BLUE,(0,-0.18),2.1)
    ax.text(0.82,2.0,r"$K=0$",fontsize=12); ax.text(0.76,-0.45,"閉路",fontsize=11)
    # residual
    q0=(2.3,0.1); q1=(3.05,1.65); q2=(3.7,0.55); q3=(2.55,0.25)
    vec(ax,q0,q1,r"$Z_1$",BLUE,(-0.18,0.05),2.1); vec(ax,q1,q2,r"$Z_2$",ORANGE,(0.02,0.03),2.1); vec(ax,q2,q3,r"$Z_3$",BLUE,(0,-0.18),2.1)
    ax.annotate("",xy=q0,xytext=q3,arrowprops=dict(arrowstyle="->",lw=2,color=RED,linestyle="--"))
    ax.text(2.55,0.5,r"$K\neq0$",color=RED,fontsize=11)
    ax.text(2.72,-0.45,"残差",fontsize=11)

    cap=fig.add_subplot(gs[1,:]);cap.axis("off")
    cap.text(0.5,0.5,"自明な恒等式を、複素関係の和と閉路残差として読む。",ha="center",va="center",fontsize=12)
    save(fig,"複素数の閉性と閉路残差.png")


def fig2():
    fig=plt.figure(figsize=(12,8)); gs=fig.add_gridspec(2,2,height_ratios=[12,1.2],wspace=.28)
    fig.suptitle("図2　三体では一つが他の二つから決まる",fontsize=20,y=.98)

    ax=fig.add_subplot(gs[0,0]);clean_ax(ax,xlim=(-.2,3),ylim=(-.2,2.7));ax.set_title("A. 三つの複素ベクトルの三角閉路",fontsize=13)
    p0=(0,0);p1=(1.15,1.85);p2=(2.5,.4)
    vec(ax,p0,p1,r"$W_1$",BLUE,(-.15,.05));vec(ax,p1,p2,r"$W_2$",ORANGE,(.05,.05));vec(ax,p2,p0,r"$W_3$",BLUE,(0,-.2))
    ax.text(.35,-.45,r"$W_1+W_2+W_3=0$",fontsize=12)
    ax.text(.6,-.83,r"$W_3=-(W_1+W_2)$",fontsize=11)

    ax=fig.add_subplot(gs[0,1]);clean_ax(ax,xlim=(-2.2,2.2),ylim=(-1.8,2.2));ax.set_title("B. 二乗量の中心化",fontsize=13)
    pts=[(0,1.4),(1.3,-.8),(-1.35,-.65)]; cols=[BLUE,ORANGE,BLUE]; labs=[r"$W_1$",r"$W_2$",r"$W_3$"]
    for p,c,l in zip(pts,cols,labs):vec(ax,(0,0),p,l,c,(.05,.05),2.2)
    poly=Polygon(pts,closed=True,fill=False,ls="--",lw=1.4,ec=GRAY);ax.add_patch(poly)
    ax.text(-1.85,1.92,r"$Z_1+Z_2+Z_3=R^2$",fontsize=11)
    ax.text(-1.85,1.60,r"$W_k:=Z_k-R^2/3$",fontsize=11)
    ax.text(-1.85,1.28,r"$W_1+W_2+W_3=0$",fontsize=11)

    cap=fig.add_subplot(gs[1,:]);cap.axis("off");cap.text(.5,.5,"三体は、関係の閉路が最初に現れる最小の非自明系である。",ha="center",va="center",fontsize=12)
    save(fig,"三体の三角閉路とベクトル関係.png")


def fig3():
    fig=plt.figure(figsize=(14,8));gs=fig.add_gridspec(2,3,height_ratios=[12,1.2],wspace=.24)
    fig.suptitle("図3　四頂点では三角閉路だけでは完全ネットワークは導けない",fontsize=18,y=.98)
    # A
    ax=fig.add_subplot(gs[0,0]);ax.axis("off");ax.set_aspect("equal");ax.set_xlim(0,4);ax.set_ylim(0,4);ax.set_title("A. 四頂点5関係",fontsize=13)
    P={1:(.6,.7),2:(.6,3.2),3:(3.2,3.2),4:(3.2,.7)}
    ax.add_patch(Polygon([P[1],P[2],P[3]],closed=True,fc=PALE_BLUE,ec="none"));ax.add_patch(Polygon([P[1],P[3],P[4]],closed=True,fc=PALE_ORANGE,ec="none"))
    for a,b in [(1,2),(2,3),(3,4),(4,1)]: ax.plot([P[a][0],P[b][0]],[P[a][1],P[b][1]],lw=2.2,c=BLUE)
    ax.plot([P[1][0],P[3][0]],[P[1][1],P[3][1]],lw=2.2,c=ORANGE)
    for k,p in P.items():ax.scatter([p[0]],[p[1]],c="black",s=22);ax.text(p[0],p[1]+(.18 if k in [2,3] else -.3),str(k),ha="center",fontsize=11)
    ax.text(1.25,2.2,"(123)",fontsize=11);ax.text(2.1,1.45,"(134)",fontsize=11)
    ax.text(.2,.1,r"$W_{12}+W_{23}+W_{31}=0$"+"\n"+r"$W_{13}+W_{34}+W_{41}=0$",fontsize=10)
    # B
    ax=fig.add_subplot(gs[0,1]);ax.axis("off");ax.set_aspect("equal");ax.set_xlim(0,4);ax.set_ylim(0,4);ax.set_title(r"B. 完全ネットワーク $K_4$",fontsize=13)
    Q={1:(.55,1.1),2:(2,3.4),3:(3.45,1.1),4:(2,.25)}
    faces=[([Q[1],Q[2],Q[4]],PALE_BLUE,"(124)"),([Q[2],Q[3],Q[4]],PALE_ORANGE,"(234)"),([Q[1],Q[3],Q[4]],PALE_GREEN,"(134)"),([Q[1],Q[2],Q[3]],PALE_PINK,"(123)")]
    for pts,c,l in faces: ax.add_patch(Polygon(pts,closed=True,fc=c,alpha=.45,ec="none"))
    for a,b in [(1,2),(1,3),(1,4),(2,3),(2,4),(3,4)]:ax.plot([Q[a][0],Q[b][0]],[Q[a][1],Q[b][1]],c=BLUE,lw=2)
    for k,p in Q.items():ax.scatter([p[0]],[p[1]],c="black",s=22);ax.text(p[0],p[1]+(.17 if k==2 else -.25),str(k),ha="center")
    ax.text(.6,.05,"三角閉路だけからは K4 は強制されない",fontsize=10)
    # C dual
    ax=fig.add_subplot(gs[0,2]);ax.axis("off");ax.set_aspect("equal");ax.set_xlim(0,5);ax.set_ylim(0,4);ax.set_title("C. 辺が面を共有",fontsize=13)
    R={1:(.3,.6),2:(.3,3),3:(2,3),4:(2,.6)}
    for a,b in [(1,2),(2,3),(3,4),(4,1)]:ax.plot([R[a][0],R[b][0]],[R[a][1],R[b][1]],c=BLUE,lw=2)
    ax.plot([R[1][0],R[3][0]],[R[1][1],R[3][1]],c=ORANGE,lw=2)
    ax.text(.65,1.8,"対角辺: 2面",fontsize=9);ax.text(.6,.2,"外周辺: 1面",fontsize=9);ax.text(.8,-.05,"非対称",fontsize=11)
    S={1:(2.8,.8),2:(3.7,3.2),3:(4.65,.8),4:(3.7,.15)}
    for a,b in [(1,2),(1,3),(1,4),(2,3),(2,4),(3,4)]:ax.plot([S[a][0],S[b][0]],[S[a][1],S[b][1]],c=BLUE,lw=1.8)
    for p in S.values():ax.scatter([p[0]],[p[1]],c="black",s=18)
    ax.text(3.2,.0,"全辺: 2面\n対称",fontsize=10,ha="center")
    cap=fig.add_subplot(gs[1,:]);cap.axis("off");cap.text(.5,.5,"三角閉路の成立と、無名性を保つ完全対称な関係構造とは別問題である。",ha="center",va="center",fontsize=12)
    save(fig,"四頂点ネットワークの三角閉路比較.png")


def tetra(ax, center=(0,0),scale=1.0,facecolor=PALE_BLUE,label=None):
    cx,cy=center
    P=np.array([[0,1.1],[-1,-.65],[1,-.65],[0,-1.15]])*scale+np.array([cx,cy])
    faces=[(0,1,3),(0,2,3),(1,2,3),(0,1,2)]
    for f in faces:
        ax.add_patch(Polygon(P[list(f)],closed=True,fc=facecolor,alpha=.45,ec="none"))
    for a,b in [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]:ax.plot([P[a,0],P[b,0]],[P[a,1],P[b,1]],c="black",lw=1.2)
    ax.scatter(P[:,0],P[:,1],c="black",s=16)
    if label: ax.text(cx,cy-1.55*scale,label,ha="center",fontsize=10)


def fig4():
    fig=plt.figure(figsize=(14,8));gs=fig.add_gridspec(2,3,height_ratios=[12,1.2],wspace=.28)
    fig.suptitle("図4　内部相殺とフラックス：閉じた多面体を凝縮体として読む",fontsize=18,y=.98)
    ax=fig.add_subplot(gs[0,0]);ax.axis("off");ax.set_aspect("equal");ax.set_xlim(0,4);ax.set_ylim(0,4);ax.set_title("A. 隣接する二つの三角形面",fontsize=13)
    p1=(2,3.4);p2=(.4,.6);p3=(2,.6);p4=(3.6,.6)
    ax.add_patch(Polygon([p1,p2,p3],closed=True,fc=PALE_BLUE,ec="none"));ax.add_patch(Polygon([p1,p3,p4],closed=True,fc=PALE_ORANGE,ec="none"))
    for p in [p1,p2,p3,p4]:ax.scatter([p[0]],[p[1]],c="black",s=20)
    vec(ax,p2,p1,r"$Z_{12}$",BLUE,(-.2,.05),2);vec(ax,p2,p3,r"$Z_{23}$",BLUE,(0,-.18),2);vec(ax,p1,p4,r"$Z_{14}$",ORANGE,(.1,.05),2);vec(ax,p3,p4,r"$Z_{34}$",ORANGE,(0,-.18),2)
    ax.annotate("",xy=p3,xytext=p1,arrowprops=dict(arrowstyle="->",lw=2,color=BLUE));ax.annotate("",xy=p1,xytext=p3,arrowprops=dict(arrowstyle="->",lw=2,color=ORANGE))
    ax.text(2.08,2.0,r"$Z_{13}$ / $-Z_{13}$",fontsize=9);ax.text(1.9,1.75,"×",color=RED,fontsize=20)
    ax.text(.55,.15,r"$+Z_{13}+(-Z_{13})=0$",fontsize=11,bbox=dict(boxstyle="round,pad=.3",fc=PALE_PINK,ec="none"))
    ax=fig.add_subplot(gs[0,1]);ax.axis("off");ax.set_aspect("equal");ax.set_xlim(-2,2);ax.set_ylim(-2,2.2);ax.set_title("B. 閉じた多面体の面集合",fontsize=13)
    tetra(ax,(0,.35),1.2,PALE_BLUE,"凝縮体")
    ax.text(-1.5,1.65,r"$K_{f_1}$",fontsize=10);ax.text(1.0,.85,r"$K_{f_2}$",fontsize=10);ax.text(-1.55,.2,r"$K_{f_3}$",fontsize=10);ax.text(.45,-.9,r"$K_{f_4}$",fontsize=10)
    ax=fig.add_subplot(gs[0,2]);ax.axis("off");ax.set_aspect("equal");ax.set_xlim(0,5);ax.set_ylim(0,4);ax.set_title("C. 凝縮体の階層的な成長",fontsize=13)
    tetra(ax,(1,3.0),.45,PALE_BLUE,"A");tetra(ax,(1,1.95),.45,PALE_ORANGE,"B");tetra(ax,(1,.85),.45,PALE_GREEN,"C")
    for y in [3.0,1.95,.85]:ax.annotate("",xy=(2.5,1.95),xytext=(1.55,y),arrowprops=dict(arrowstyle="->",lw=1.5))
    tetra(ax,(3.65,1.95),.95,"#E8DFF5","より大きな凝縮体")
    ax.text(3.3,3.55,r"$K_f\in\mathbb{C}$",fontsize=12)
    cap=fig.add_subplot(gs[1,:]);cap.axis("off");cap.text(.5,.5,"共有された内部関係は相殺され、外部からは閉じた多面体が一つの凝縮体として現れる。",ha="center",va="center",fontsize=12)
    save(fig,"内部相殺と凝縮体の階層的成長.png")


def box(ax,xy,w,h,text,fc=PALE_BLUE,ec="black",fontsize=10):
    x,y=xy
    r=FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.03",fc=fc,ec=ec,lw=1.2)
    ax.add_patch(r);ax.text(x+w/2,y+h/2,text,ha="center",va="center",fontsize=fontsize);return r


def fig5():
    fig,ax=plt.subplots(figsize=(12,10));ax.axis("off");ax.set_xlim(0,10);ax.set_ylim(0,13)
    ax.set_title("図5　思考実験の方法：仮定を削ぎ落として残ったもの",fontsize=19,pad=18)
    steps=[("出発点： $z\in\mathbb{C}$",11.7),(r"$z^2\in\mathbb{C}$",10.6),(r"$\sum z_n^2\in\mathbb{C}$",9.5),(r"$Z_n:=z_n^2,\ \sum Z_n=K$",8.4)]
    for i,(t,y) in enumerate(steps):
        box(ax,(3.3,y),3.4,.65,t)
        if i<len(steps)-1:ax.annotate("",xy=(5,y-.05),xytext=(5,y-.45),arrowprops=dict(arrowstyle="->"))
    box(ax,(2.7,7.25),4.6,.7,"ここで何を仮定していないか？",fc="#F2F4F4")
    rejects=["完全ネットワーク","比による比較","= 判定の先取り","ゼロ閉塞の先取り","三角形の先取り"]
    xs=[.25,2.15,4.05,5.95,7.85]
    for x,t in zip(xs,rejects):
        box(ax,(x,6.05),1.7,.65,t,fc=PALE_PINK,fontsize=9)
        ax.plot([x+.15,x+1.55],[6.15,6.6],c=RED,lw=2)
        ax.plot([5, x+.85],[7.25,6.72],c=GRAY,lw=.8)
    main=[("$K=0$ → 閉塞",4.75),("$K\\neq0$ → 残差",3.75),("三項 → 閉路",2.75),("共有部分の相殺 → 面フラックス",1.75),("閉じた面集合 → 凝縮体",.75)]
    prev=(5,6.05)
    for t,y in main:
        box(ax,(3.05,y),3.9,.65,t)
        ax.annotate("",xy=(5,y+.65),xytext=prev,arrowprops=dict(arrowstyle="->"));prev=(5,y)
    box(ax,(2.0,-.3),6.0,.7,"主張：新しい方程式ではなく、最小仮定の読み方の発明",fc=PALE_YELLOW,fontsize=11)
    ax.text(5,-.8,"核心は、何を導入したかではなく、何を導入しなくても済んだかにある。",ha="center",fontsize=11)
    save(fig,"図5_思考実験の方法フローチャート.png")


if __name__ == "__main__":
    fig1(); fig2(); fig3(); fig4(); fig5()
