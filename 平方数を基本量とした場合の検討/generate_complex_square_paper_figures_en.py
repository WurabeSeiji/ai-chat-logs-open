#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""English-figure variant of generate_complex_square_paper_figures.py.

Faithful copy of the Japanese original; the only changes are
(1) label strings translated to English,
(2) output file names in English,
(3) minimal box-width / text-position adjustments where the English
    strings are longer than the Japanese ones.

Outputs (PNG):
  1. fig1_complex_closure_loop_residual_en.png
  2. fig2_three_body_triangular_loop_en.png
  3. fig3_four_vertex_triangular_loops_en.png
  4. fig4_internal_cancellation_condensate_en.png
  5. fig5_thought_experiment_flowchart_en.png
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
    fig.suptitle("Fig. 1  From the closure of complex numbers to the loop residual", fontsize=20, y=0.98)

    ax = fig.add_subplot(gs[0,0]); clean_ax(ax)
    ax.set_title("A. Closure of complex numbers", fontsize=13, pad=8)
    vec(ax,(0,0),(1.55,1.35),"z",BLUE,(0.05,0.08))
    ax.scatter([1.55],[1.35],s=20,c="black")
    ax.text(0.15,-0.45,r"$z\in\mathbb{C}\Rightarrow z^2\in\mathbb{C}$",fontsize=12)
    ax.text(0.75,-0.82,r"$Z=z^2$",fontsize=12)

    ax = fig.add_subplot(gs[0,1]); clean_ax(ax,xlim=(-0.2,3.8),ylim=(-0.2,3.0))
    ax.set_title("B. Sum of complex relations", fontsize=13, pad=8)
    p0=(0,0);p1=(1.15,0.75);p2=(1.9,1.75);p3=(3.0,2.15)
    vec(ax,p0,p1,r"$Z_1$",BLUE,(0.02,0.1)); vec(ax,p1,p2,r"$Z_2$",ORANGE,(0.04,0.08)); vec(ax,p2,p3,r"$Z_3$",BLUE,(0.03,0.08))
    for p in [p1,p2,p3]: ax.scatter([p[0]],[p[1]],s=20,c="black")
    ax.text(0.2,-0.45,r"$\sum Z_n=K$",fontsize=13)
    ax.text(0.2,-0.8,r"$Z_n:=z_n^2$",fontsize=11)
    ax.text(1.9,-0.8,r"$K\in\mathbb{C}$",fontsize=11)

    ax = fig.add_subplot(gs[0,2]); clean_ax(ax,xlim=(-0.2,4.2),ylim=(-0.2,2.9))
    ax.set_title("C. Loop and residual", fontsize=13, pad=8)
    # closed triangle
    p0=(0.1,0.1); p1=(0.9,1.65); p2=(1.65,0.3)
    vec(ax,p0,p1,r"$Z_1$",BLUE,(-0.18,0.05),2.1); vec(ax,p1,p2,r"$Z_2$",ORANGE,(0.04,0.02),2.1); vec(ax,p2,p0,r"$Z_3$",BLUE,(0,-0.18),2.1)
    ax.text(0.82,2.0,r"$K=0$",fontsize=12); ax.text(0.6,-0.45,"closed loop",fontsize=11)
    # residual
    q0=(2.3,0.1); q1=(3.05,1.65); q2=(3.7,0.55); q3=(2.55,0.25)
    vec(ax,q0,q1,r"$Z_1$",BLUE,(-0.18,0.05),2.1); vec(ax,q1,q2,r"$Z_2$",ORANGE,(0.02,0.03),2.1); vec(ax,q2,q3,r"$Z_3$",BLUE,(0,-0.18),2.1)
    ax.annotate("",xy=q0,xytext=q3,arrowprops=dict(arrowstyle="->",lw=2,color=RED,linestyle="--"))
    ax.text(2.98,0.5,r"$K\neq0$",color=RED,fontsize=11)
    ax.text(2.6,-0.45,"residual",fontsize=11)

    cap=fig.add_subplot(gs[1,:]);cap.axis("off")
    cap.text(0.5,0.5,"The trivial identity is read as a sum of complex relations and a loop residual.",ha="center",va="center",fontsize=12)
    save(fig,"fig1_complex_closure_loop_residual_en.png")


def fig2():
    fig=plt.figure(figsize=(12,8)); gs=fig.add_gridspec(2,2,height_ratios=[12,1.2],wspace=.28)
    fig.suptitle("Fig. 2  Three bodies: one is determined by the other two",fontsize=20,y=.98)

    ax=fig.add_subplot(gs[0,0]);clean_ax(ax,xlim=(-.2,3),ylim=(-.2,2.7));ax.set_title("A. Triangular loop of three complex vectors",fontsize=13)
    p0=(0,0);p1=(1.15,1.85);p2=(2.5,.4)
    vec(ax,p0,p1,r"$W_1$",BLUE,(-.15,.05));vec(ax,p1,p2,r"$W_2$",ORANGE,(.05,.05));vec(ax,p2,p0,r"$W_3$",BLUE,(0,-.2))
    ax.text(.35,-.45,r"$W_1+W_2+W_3=0$",fontsize=12)
    ax.text(.6,-.83,r"$W_3=-(W_1+W_2)$",fontsize=11)

    ax=fig.add_subplot(gs[0,1]);clean_ax(ax,xlim=(-2.2,2.2),ylim=(-1.8,2.2));ax.set_title("B. Centering the squared quantities",fontsize=13)
    pts=[(0,1.4),(1.3,-.8),(-1.35,-.65)]; cols=[BLUE,ORANGE,BLUE]; labs=[r"$W_1$",r"$W_2$",r"$W_3$"]
    for p,c,l in zip(pts,cols,labs):vec(ax,(0,0),p,l,c,(.05,.05),2.2)
    poly=Polygon(pts,closed=True,fill=False,ls="--",lw=1.4,ec=GRAY);ax.add_patch(poly)
    ax.text(-1.85,1.92,r"$Z_1+Z_2+Z_3=R^2$",fontsize=11)
    ax.text(-1.85,1.60,r"$W_k:=Z_k-R^2/3$",fontsize=11)
    ax.text(-1.85,1.28,r"$W_1+W_2+W_3=0$",fontsize=11)

    cap=fig.add_subplot(gs[1,:]);cap.axis("off");cap.text(.5,.5,"Three bodies form the smallest nontrivial system in which a loop of relations first appears.",ha="center",va="center",fontsize=12)
    save(fig,"fig2_three_body_triangular_loop_en.png")


def fig3():
    fig=plt.figure(figsize=(14,8));gs=fig.add_gridspec(2,3,height_ratios=[12,1.2],wspace=.24)
    fig.suptitle("Fig. 3  Four vertices: triangular loops alone do not force the complete network",fontsize=18,y=.98)
    # A
    ax=fig.add_subplot(gs[0,0]);ax.axis("off");ax.set_aspect("equal");ax.set_xlim(0,4);ax.set_ylim(0,4);ax.set_title("A. Four vertices, five relations",fontsize=13)
    P={1:(.6,.7),2:(.6,3.2),3:(3.2,3.2),4:(3.2,.7)}
    ax.add_patch(Polygon([P[1],P[2],P[3]],closed=True,fc=PALE_BLUE,ec="none"));ax.add_patch(Polygon([P[1],P[3],P[4]],closed=True,fc=PALE_ORANGE,ec="none"))
    for a,b in [(1,2),(2,3),(3,4),(4,1)]: ax.plot([P[a][0],P[b][0]],[P[a][1],P[b][1]],lw=2.2,c=BLUE)
    ax.plot([P[1][0],P[3][0]],[P[1][1],P[3][1]],lw=2.2,c=ORANGE)
    for k,p in P.items():ax.scatter([p[0]],[p[1]],c="black",s=22);ax.text(p[0],p[1]+(.18 if k in [2,3] else -.3),str(k),ha="center",fontsize=11)
    ax.text(1.25,2.2,"(123)",fontsize=11);ax.text(2.1,1.45,"(134)",fontsize=11)
    ax.text(.2,-.2,r"$W_{12}+W_{23}+W_{31}=0$"+"\n"+r"$W_{13}+W_{34}+W_{41}=0$",fontsize=10)
    # B
    ax=fig.add_subplot(gs[0,1]);ax.axis("off");ax.set_aspect("equal");ax.set_xlim(0,4);ax.set_ylim(0,4);ax.set_title(r"B. Complete network $K_4$",fontsize=13)
    Q={1:(.55,1.1),2:(2,3.4),3:(3.45,1.1),4:(2,.25)}
    faces=[([Q[1],Q[2],Q[4]],PALE_BLUE,"(124)"),([Q[2],Q[3],Q[4]],PALE_ORANGE,"(234)"),([Q[1],Q[3],Q[4]],PALE_GREEN,"(134)"),([Q[1],Q[2],Q[3]],PALE_PINK,"(123)")]
    for pts,c,l in faces: ax.add_patch(Polygon(pts,closed=True,fc=c,alpha=.45,ec="none"))
    for a,b in [(1,2),(1,3),(1,4),(2,3),(2,4),(3,4)]:ax.plot([Q[a][0],Q[b][0]],[Q[a][1],Q[b][1]],c=BLUE,lw=2)
    for k,p in Q.items():ax.scatter([p[0]],[p[1]],c="black",s=22);ax.text(p[0],p[1]+(.17 if k==2 else -.25),str(k),ha="center")
    ax.text(.5,-.45,r"Triangular loops alone do not force $K_4$",fontsize=9)
    # C dual
    ax=fig.add_subplot(gs[0,2]);ax.axis("off");ax.set_aspect("equal");ax.set_xlim(0,5);ax.set_ylim(0,4);ax.set_title("C. Faces shared by edges",fontsize=13)
    R={1:(.3,.6),2:(.3,3),3:(2,3),4:(2,.6)}
    for a,b in [(1,2),(2,3),(3,4),(4,1)]:ax.plot([R[a][0],R[b][0]],[R[a][1],R[b][1]],c=BLUE,lw=2)
    ax.plot([R[1][0],R[3][0]],[R[1][1],R[3][1]],c=ORANGE,lw=2)
    ax.text(.65,1.8,"diagonal: 2 faces",fontsize=9);ax.text(.5,.2,"outer edges: 1 face",fontsize=9);ax.text(.7,-.05,"asymmetric",fontsize=11)
    S={1:(2.8,.8),2:(3.7,3.2),3:(4.65,.8),4:(3.7,.15)}
    for a,b in [(1,2),(1,3),(1,4),(2,3),(2,4),(3,4)]:ax.plot([S[a][0],S[b][0]],[S[a][1],S[b][1]],c=BLUE,lw=1.8)
    for p in S.values():ax.scatter([p[0]],[p[1]],c="black",s=18)
    ax.text(3.7,-.55,"all edges: 2 faces\nsymmetric",fontsize=10,ha="center")
    cap=fig.add_subplot(gs[1,:]);cap.axis("off");cap.text(.5,.5,"Closing triangular loops and a fully symmetric relational structure preserving namelessness are separate issues.",ha="center",va="center",fontsize=12)
    save(fig,"fig3_four_vertex_triangular_loops_en.png")


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
    fig.suptitle("Fig. 4  Internal cancellation and flux: reading a closed polyhedron as a condensate",fontsize=18,y=.98)
    ax=fig.add_subplot(gs[0,0]);ax.axis("off");ax.set_aspect("equal");ax.set_xlim(0,4);ax.set_ylim(0,4);ax.set_title("A. Two adjacent triangular faces",fontsize=13)
    p1=(2,3.4);p2=(.4,.6);p3=(2,.6);p4=(3.6,.6)
    ax.add_patch(Polygon([p1,p2,p3],closed=True,fc=PALE_BLUE,ec="none"));ax.add_patch(Polygon([p1,p3,p4],closed=True,fc=PALE_ORANGE,ec="none"))
    for p in [p1,p2,p3,p4]:ax.scatter([p[0]],[p[1]],c="black",s=20)
    vec(ax,p2,p1,r"$Z_{12}$",BLUE,(-.2,.05),2);vec(ax,p2,p3,r"$Z_{23}$",BLUE,(0,-.18),2);vec(ax,p1,p4,r"$Z_{14}$",ORANGE,(.1,.05),2);vec(ax,p3,p4,r"$Z_{34}$",ORANGE,(0,-.18),2)
    ax.annotate("",xy=p3,xytext=p1,arrowprops=dict(arrowstyle="->",lw=2,color=BLUE));ax.annotate("",xy=p1,xytext=p3,arrowprops=dict(arrowstyle="->",lw=2,color=ORANGE))
    ax.text(2.12,1.35,r"$Z_{13}$ / $-Z_{13}$",fontsize=9);ax.text(1.9,1.75,"×",color=RED,fontsize=20)
    ax.text(.55,.15,r"$+Z_{13}+(-Z_{13})=0$",fontsize=11,bbox=dict(boxstyle="round,pad=.3",fc=PALE_PINK,ec="none"))
    ax=fig.add_subplot(gs[0,1]);ax.axis("off");ax.set_aspect("equal");ax.set_xlim(-2,2);ax.set_ylim(-2,2.2);ax.set_title("B. Face set of a closed polyhedron",fontsize=13)
    tetra(ax,(0,.35),1.2,PALE_BLUE,"condensate")
    ax.text(-1.5,1.65,r"$K_{f_1}$",fontsize=10);ax.text(1.0,.85,r"$K_{f_2}$",fontsize=10);ax.text(-1.55,.2,r"$K_{f_3}$",fontsize=10);ax.text(.85,-1.3,r"$K_{f_4}$",fontsize=10)
    ax=fig.add_subplot(gs[0,2]);ax.axis("off");ax.set_aspect("equal");ax.set_xlim(0,5);ax.set_ylim(0,4);ax.set_title("C. Hierarchical growth of condensates",fontsize=13)
    tetra(ax,(1,3.0),.45,PALE_BLUE);tetra(ax,(1,1.95),.45,PALE_ORANGE);tetra(ax,(1,.85),.45,PALE_GREEN)
    for yy,lab in [(3.0,"A"),(1.95,"B"),(.85,"C")]: ax.text(.25,yy,lab,ha="center",va="center",fontsize=10)
    for y in [3.0,1.95,.85]:ax.annotate("",xy=(2.5,1.95),xytext=(1.55,y),arrowprops=dict(arrowstyle="->",lw=1.5))
    tetra(ax,(3.65,1.95),.95,"#E8DFF5","larger condensate")
    ax.text(3.3,3.55,r"$K_f\in\mathbb{C}$",fontsize=12)
    cap=fig.add_subplot(gs[1,:]);cap.axis("off");cap.text(.5,.5,"Shared internal relations cancel; from outside, the closed polyhedron appears as a single condensate.",ha="center",va="center",fontsize=12)
    save(fig,"fig4_internal_cancellation_condensate_en.png")


def box(ax,xy,w,h,text,fc=PALE_BLUE,ec="black",fontsize=10):
    x,y=xy
    r=FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.03",fc=fc,ec=ec,lw=1.2)
    ax.add_patch(r);ax.text(x+w/2,y+h/2,text,ha="center",va="center",fontsize=fontsize);return r


def fig5():
    fig,ax=plt.subplots(figsize=(12,10));ax.axis("off");ax.set_xlim(0,10);ax.set_ylim(-1.3,13)
    ax.set_title("Fig. 5  Method of the thought experiment: what remains after shaving off assumptions",fontsize=16,pad=18)
    steps=[("Starting point: $z\in\mathbb{C}$",11.7),(r"$z^2\in\mathbb{C}$",10.6),(r"$\sum z_n^2\in\mathbb{C}$",9.5),(r"$Z_n:=z_n^2,\ \sum Z_n=K$",8.4)]
    for i,(t,y) in enumerate(steps):
        box(ax,(3.3,y),3.4,.65,t)
        if i<len(steps)-1:ax.annotate("",xy=(5,y-.45),xytext=(5,y-.05),arrowprops=dict(arrowstyle="->"))
    box(ax,(2.7,7.25),4.6,.7,"What is not assumed here?",fc="#F2F4F4")
    rejects=["complete\nnetwork","comparison\nby ratio","presupposed\n= test","presupposed\nzero closure","presupposed\ntriangle"]
    xs=[.25,2.15,4.05,5.95,7.85]
    for x,t in zip(xs,rejects):
        box(ax,(x,6.05),1.7,.65,t,fc=PALE_PINK,fontsize=8)
        ax.plot([x+.15,x+1.55],[6.15,6.6],c=RED,lw=2)
        ax.plot([5, x+.85],[7.25,6.72],c=GRAY,lw=.8)
    main=[("$K=0$ → closure",4.75),("$K\\neq0$ → residual",3.75),("three terms → loop",2.75),("cancellation of shared parts → face flux",1.75),("closed face set → condensate",.75)]
    prev=(5,6.05)
    for t,y in main:
        box(ax,(2.7,y),4.6,.65,t)
        ax.annotate("",xy=(5,y+.65),xytext=prev,arrowprops=dict(arrowstyle="->"));prev=(5,y)
    box(ax,(1.5,-.3),7.0,.7,"Claim: not a new equation, but a minimal-assumption way of reading",fc=PALE_YELLOW,fontsize=11)
    ax.text(5,-.8,"The core lies not in what was introduced, but in what never needed to be introduced.",ha="center",fontsize=11)
    save(fig,"fig5_thought_experiment_flowchart_en.png")


if __name__ == "__main__":
    fig1(); fig2(); fig3(); fig4(); fig5()
