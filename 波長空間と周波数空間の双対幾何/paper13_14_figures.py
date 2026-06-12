#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 論文13・14 図版 (全て実計算、英語ラベル)
import numpy as np, itertools
from collections import defaultdict
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
src=open('supplement57_convention_dependence.py',encoding='utf-8').read()
exec(src.split('# 現行(字句順序)規約')[0])

# ===== 論文13 図1: 位置=符号セクター (¼並進の4状態 + 階層格子) =====
fig,ax=plt.subplots(1,2,figsize=(12,4.2))
x=np.linspace(0,1,400)
labels=['cos (d=0)','sin (d=1)','-cos (d=2)','-sin (d=3)']
for d in range(4):
    ax[0].plot(x,np.sqrt(2)*np.cos(2*np.pi*(x-d/4)),label=labels[d],lw=1.4)
ax[0].set_xlabel('x'); ax[0].set_ylabel('wave'); ax[0].legend(fontsize=8,ncol=2)
ax[0].set_title('(a) Position = sign sector: the four quarter-translates\nof one fundamental (no new state added)')
for q in range(5): ax[0].axvline(q/4,color='gray',lw=0.4)
levels=[(0,1.0),(1,1/3),(2,1/9)]
for lv,unit in levels:
    pts=np.arange(0,1.0001,unit/4)
    ax[1].plot(pts,[lv]*len(pts),'|',ms=18-5*lv,mew=1.2)
    ax[1].text(1.02,lv,f'level {lv}: step {unit}/4 = 1/{int(4/unit)}',fontsize=8,va='center')
ax[1].set_ylim(-0.5,2.5); ax[1].set_xlim(0,1.15); ax[1].set_yticks([])
ax[1].set_xlabel('x')
ax[1].set_title('(b) Hierarchical quarter lattices $P_\\ell=(4\\cdot3^\\ell)^{-1}\\mathbb{Z}$:\nnested refinement (index 3), completion = continuum')
plt.tight_layout(); plt.savefig('paper13_fig1_sign_sector.png',dpi=150); plt.close()

# ===== 論文13 図2: 読み出し三層 11 ⊂ 13 ⊂ 21 =====
PERMS=list(itertools.permutations(range(4)))
SIGNS=list(itertools.product([1,-1],repeat=4))
G=[(p,s_) for p in PERMS for s_ in SIGNS]
def act(g,v):
    p,s_=g
    return tuple(s_[i]*v[p[i]] for i in range(4))
def canon_pair(u,v):
    best=None
    for g in G:
        a,b=act(g,u),act(g,v)
        for pr in ((a,b),(b,a)):
            if best is None or pr<best: best=pr
    return best
def unsigned(w): return max(w,tuple(-x for x in w))
def stripe_lines(u,v):
    s1=tuple(a+b for a,b in zip(u,v)); s2=tuple(a-b for a,b in zip(u,v))
    best=None
    for g in G:
        pr=tuple(sorted([unsigned(act(g,s1)),unsigned(act(g,s2))]))
        if best is None or pr<best: best=pr
    return best
def n2(w): return sum(t*t for t in w)
cells64=SH[9]
orb=set(); lines=set(); scal=set()
for i,u in enumerate(cells64):
    for j in range(i+1,len(cells64)):
        v=cells64[j]
        orb.add(canon_pair(u,v)); lines.add(stripe_lines(u,v))
        scal.add(tuple(sorted([n2(tuple(a-b for a,b in zip(u,v))),n2(tuple(a+b for a,b in zip(u,v)))])))
fig,ax=plt.subplots(figsize=(8.5,4))
layers=['scalar fringe\n(spacings only)','line-position fringe\n(unsigned line pair)','branch channels\n(complex coefficients)']
vals=[len(scal),len(lines),len(orb)]
bars=ax.barh(layers,vals,color=['tab:gray','tab:orange','tab:green'])
for b,v in zip(bars,vals): ax.text(v+0.2,b.get_y()+b.get_height()/2,str(v),va='center',fontsize=12)
ax.set_xlabel('distinguishable relation classes (all 2016 pairs of the 64 shell-9 cells)')
ax.set_title('Readout hierarchy: each layer reads strictly more of the record\n(11 ⊂ 13 ⊂ 21; branch channels separate all orbits)')
ax.set_xlim(0,24)
plt.tight_layout(); plt.savefig('paper13_fig2_readout_layers.png',dpi=150); plt.close()

# ===== 論文13 図3: 時間=導出量 (全奇数標的・時計) =====
hit=set()
K=14
for k in itertools.product(range(0,K),repeat=4):
    s4=sum((2*t+1)**2 for t in k)
    if s4%4==0 and s4//4<=199: hit.add(s4//4)
odds=list(range(1,200,2))
fig,ax=plt.subplots(1,2,figsize=(12,4))
ax[0].scatter([o for o in odds if o in hit],[1]*len([o for o in odds if o in hit]),s=14,color='tab:green',label='reachable')
miss=[o for o in odds if o not in hit]
if miss: ax[0].scatter(miss,[1]*len(miss),s=14,color='red',label='missed')
ax[0].set_yticks([]); ax[0].set_xlabel('s (norm-clock target $\\nu_t^2=s$)')
ax[0].set_title(f'(a) Targets of the derived clock: ALL odd integers\n({len(hit)}/{len(odds)} odds ≤ 199 reachable; spacing 2 — most regular possible)')
ax[0].legend(fontsize=8)
ss=np.array(sorted(hit))[:30]
ax[1].step(ss,1/np.sqrt(ss),where='post',lw=1.2)
ax[1].set_xlabel('s'); ax[1].set_ylabel('tick interval  $\\Delta t = 1/\\sqrt{s}$')
ax[1].set_title('(b) The clock is derived: tick interval is fixed by the\nconfiguration norm — no independent dial for $t$')
plt.tight_layout(); plt.savefig('paper13_fig3_derived_time.png',dpi=150); plt.close()

# ===== 論文14 図1: η の Z4 量子化 (実計算 s=9) =====
def eta_list(vp,parent,final):
    out=[]
    for ((a,b),x,(c,d),d1) in two_step_paths(parent,final):
        for ia,va in enumerate(SH[a]):
            for ib,vb in enumerate(SH[b]):
                if a==b and not ia<ib: continue
                e1=tphase(vp,va,vb)
                if e1 is None: continue
                for vx in ([va,vb] if a==b else [va if x==a else vb]):
                    for ic,vc in enumerate(SH[c]):
                        for idd,vd in enumerate(SH[d]):
                            if c==d and not ic<idd: continue
                            e2=tphase(vx,vc,vd)
                            if e2 is None: continue
                            out.append(e1*e2)
    return out
E531=eta_list((2,1,0,0),9,(1,3,5)); E333=eta_list((2,1,0,0),9,(3,3,3))
fig,ax=plt.subplots(1,2,figsize=(11,4.6))
th=np.linspace(0,2*np.pi,200)
for a_,E,nm in ((ax[0],E531,'final (5,3,1): 160 paths'),(ax[1],E333,'final (3,3,3): 36 paths')):
    a_.plot(np.cos(th),np.sin(th),color='gray',lw=0.5)
    cnt=defaultdict(int)
    for e in E: cnt[(round(e.real,6),round(e.imag,6))]+=1
    for (re,im),n in cnt.items():
        a_.scatter([re],[im],s=120,color='tab:blue')
        a_.annotate(f'×{n}',(re,im),textcoords='offset points',xytext=(10,8))
    a_.set_xlim(-1.5,1.5); a_.set_ylim(-1.5,1.5); a_.set_aspect('equal')
    a_.axhline(0,color='gray',lw=0.4); a_.axvline(0,color='gray',lw=0.4)
    a_.set_title(f'transport phase $\\eta$ — {nm}\n(every path lands exactly on $Z_4$)')
plt.tight_layout(); plt.savefig('paper14_fig1_eta_z4.png',dpi=150); plt.close()

# ===== 論文14 図2: 正準規約 (機械検証済み数値の提示) =====
fig,ax=plt.subplots(figsize=(9,4.2))
cases=['unordered\nno diagonal','unordered\nwith diagonal','ordered\nno diagonal','ordered\nwith diagonal']
w333=[40,40,160,160]
bars=ax.bar(cases,w333,color=['tab:green','tab:green','tab:red','tab:red'])
for b,v in zip(bars,w333): ax.text(b.get_x()+b.get_width()/2,v+3,str(v),ha='center')
ax.set_ylabel('$|z_{(3,3,3)}|^2$ (amplitude bookkeeping)')
ax.set_title('Canonical counting is fixed by identical-particle symmetrization alone:\nthe diagonal is irrelevant (40=40, 160=160); ordering doubles the amplitude\n(machine-verified 2×2 separation, suppl. 60 — no probability rule is selected here)')
plt.tight_layout(); plt.savefig('paper14_fig2_convention.png',dpi=150); plt.close()

# ===== 論文14 図3: ホロノミー判定 (T1/T2、検証済み数値の図解) =====
fig,ax=plt.subplots(figsize=(9.5,4.6))
ax.axis('off')
nodes={'base\n$(o_1,o_2)=(+,+)$\n(1088, 40)  class A':(0.12,0.7),
       'T2a: flip stage 1\n(320, 32)  class B':(0.55,0.92),
       'T2b: flip stage 2\n(320, 32)  class B':(0.55,0.48),
       'T2c: flip both\n(1088, 40)  class A':(0.9,0.7),
       'T1: reverse arrow\n(full conjugation)\n(1088, 40) — INVARIANT':(0.12,0.18)}
for txt,(xx,yy) in nodes.items():
    ax.text(xx,yy,txt,ha='center',va='center',fontsize=9,
            bbox=dict(boxstyle='round',fc=('lightgreen' if 'A' in txt or 'INV' in txt else 'lightcoral'),alpha=0.6))
ax.annotate('',xy=(0.42,0.88),xytext=(0.22,0.74),arrowprops=dict(arrowstyle='->'))
ax.annotate('',xy=(0.42,0.52),xytext=(0.22,0.66),arrowprops=dict(arrowstyle='->'))
ax.annotate('',xy=(0.8,0.72),xytext=(0.68,0.88),arrowprops=dict(arrowstyle='->'))
ax.annotate('',xy=(0.8,0.68),xytext=(0.68,0.52),arrowprops=dict(arrowstyle='->'))
ax.annotate('',xy=(0.12,0.3),xytext=(0.12,0.55),arrowprops=dict(arrowstyle='->',color='tab:blue'))
ax.set_title('The forced bit is NOT the arrow (T1: invariant) — it is the $Z_2$ chain holonomy of\ncanonical inversion: single-link flips toggle A↔B, double flip restores (machine, suppl. 54)')
plt.tight_layout(); plt.savefig('paper14_fig3_holonomy.png',dpi=150); plt.close()
print("figures written: 13×3, 14×3")
