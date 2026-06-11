#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Exact (computed, non-schematic) figures for papers 5-12. All text in English.
import numpy as np, itertools, math, os, sys, traceback
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from math import pi, sqrt, comb

OUT = "/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/波長空間と周波数空間の双対幾何"
plt.rcParams.update({'font.size': 10.5, 'axes.grid': True, 'grid.alpha': 0.3,
                     'figure.dpi': 110, 'savefig.dpi': 220})

def save(fig, name):
    fig.savefig(os.path.join(OUT, name + ".png"), bbox_inches='tight')
    fig.savefig(os.path.join(OUT, name + ".pdf"), bbox_inches='tight')
    plt.close(fig)
    print("saved", name)

# ---------- shared exact computations ----------
def fconv(a,b):
    n=len(a)+len(b)-1
    N=1<<max(1,(n-1)).bit_length()
    c=np.fft.irfft(np.fft.rfft(a,N)*np.fft.rfft(b,N),N)[:n]
    return np.round(c)

def N0_table(Rmax):
    """N0(R) for integer R via 1D convolution; also fractional via threshold (2R)^2."""
    lim = (2*Rmax)**2
    K = int(Rmax)+2
    ax = np.zeros(lim+1)
    for a in range(0, K+1):
        v = (2*a+1)**2
        if v <= lim: ax[v] += (1 if a == 0 else 2)
    d2 = fconv(ax, ax)[:lim+1]
    d4 = fconv(d2, d2)[:lim+1]
    cum = np.cumsum(d4)
    return cum  # cum[t] = #cells with sum(2|k|+1)^2 <= t ; N0(R)=cum[(2R)^2]

_CS={}
def count_system(system, R):
    Rmax=210
    if system not in _CS:
        if system=='torus_zero':
            vals=[((2*a+1)**2, 1 if a==0 else 2) for a in range(0, Rmax+2)]; sc=4
        elif system=='antiperiodic':
            vals=[((2*a+1)**2, 2) for a in range(0, Rmax+2)]; sc=4
        elif system=='dirichlet':
            vals=[(n**2, 1) for n in range(1, 2*Rmax+2)]; sc=4
        elif system=='torus_noshift':
            vals=[(k**2, 1 if k==0 else 2) for k in range(0, Rmax+2)]; sc=1
        L=int(sc*Rmax*Rmax)+2
        ax=np.zeros(L)
        for v,w in vals:
            if v < L: ax[v]+=w
        d2=fconv(ax,ax)[:L]; d4=fconv(d2,d2)[:L]
        _CS[system]=(np.cumsum(d4), sc)
    cum,sc=_CS[system]
    return cum[int(sc*R*R)]

CUM = N0_table(400)
def N0(R):
    return int(CUM[int(round((2*R)**2))])
def N0s(s):  # N0(sqrt(s)) for integer s
    return int(CUM[int(round(4*s))])

def g_indicator(s):
    return (pi**2/2)*s*s - N0s(s)

def partitions(n, mn=1, odd_only=False):
    if n==0:
        yield ()
        return
    start = mn
    for first in range(start, n+1):
        if odd_only and first%2==0: continue
        for rest in partitions(n-first, first, odd_only):
            yield (first,)+rest

def shell_cells4(m):
    out=[]
    K=int(math.isqrt(int(round(m))))+1
    for k in itertools.product(range(-K,K+1),repeat=4):
        if sum((abs(t)+0.5)**2 for t in k)==m: out.append(np.array(k))
    return out

def shell_count(m):
    return len(shell_cells4(float(m)))

# ============ PAPER 5 ============
def fig5_1():
    Rs = np.arange(20, 401, 10)
    gaps = [((pi**2/2)*R**4 - N0(R))/R**3 for R in Rs]
    fig, axs = plt.subplots(1, 2, figsize=(10, 3.6))
    axs[0].plot(Rs, gaps, 'o-', ms=3, lw=1, label=r'$\Delta V(R)/R^3$ (exact $N_0$)')
    axs[0].axhline(16*pi/3, color='r', ls='--', lw=1, label=r'$16\pi/3 = 16.7552$')
    axs[0].set_xlabel('$R$'); axs[0].set_ylabel(r'$\Delta V(R)/R^3$')
    axs[0].set_title('(a) Volume-gap coefficient converges to $16\\pi/3$')
    axs[0].legend()
    # (b) fitted cubic coefficients for 4 systems
    systems=[('torus_zero','Torus + zero-point $\\frac{1}{2}$'),
             ('dirichlet','Dirichlet box $[0,1]^4$'),
             ('antiperiodic','Antiperiodic torus'),
             ('torus_noshift','Torus, no shift')]
    Rf=np.arange(50,201,10); coefs=[]
    for sysk,_ in systems:
        defi=np.array([(pi**2/2)*R**4 - count_system(sysk,R) for R in Rf])
        A=np.vstack([Rf**3,Rf**2,Rf,np.ones_like(Rf)]).T
        a=np.linalg.lstsq(A,defi,rcond=None)[0][0]
        coefs.append(a)
    names=[n for _,n in systems]
    bars=axs[1].bar(range(4), coefs, color=['C0','C1','C2','C3'])
    axs[1].axhline(16*pi/3,color='r',ls='--',lw=1,label='$16\\pi/3$')
    axs[1].axhline(8*pi/3,color='m',ls=':',lw=1,label='$8\\pi/3$')
    axs[1].set_xticks(range(4)); axs[1].set_xticklabels(names, rotation=12, fontsize=8.5)
    axs[1].set_ylabel('fitted $R^3$ coefficient')
    axs[1].set_title('(b) Effective Weyl boundary term by system')
    axs[1].legend(fontsize=9)
    for i,c in enumerate(coefs):
        axs[1].text(i, c+0.4, f'{c:.2f}', ha='center', fontsize=9)
    save(fig, 'figure_paper5_1_weyl_gap_convergence')

def fig5_2():
    s=np.arange(0,41,1)
    vals=[N0s(t) if t>0 else 0 for t in s]
    fig,ax=plt.subplots(figsize=(7,3.6))
    ax.step(s, vals, where='post', lw=1.5)
    odd=[t for t in s if t%2==1 and t<=39]
    ax.plot(odd,[N0s(t) for t in odd],'ro',ms=4,label='jumps (odd $s$ only)')
    for m in [1,3,5,7,9]:
        ax.annotate(f'$[{m},{m+2})$', xy=(m+0.4, N0s(m)+4), fontsize=8)
    ax.set_xlabel(r'$s=R^2$'); ax.set_ylabel(r'$N_0(\sqrt{s})$')
    ax.set_title('Rigid tiling: each odd label owns a window of width exactly 2')
    ax.set_xlim(0,40); ax.set_ylim(0, N0s(39)*1.05); ax.legend()
    save(fig,'figure_paper5_2_rigid_tiling')

def fig5_3():
    ss=range(4,26)
    marg_all=[]; marg_odd=[]
    for s in ss:
        best_all=-1e18; best_odd=-1e18
        for p in partitions(s):
            if len(p)<2: continue
            tot=sum(g_indicator(t) for t in p)
            best_all=max(best_all,tot)
        for p in partitions(s,odd_only=True):
            if len(p)<2: continue
            if any(t%2==0 for t in p): continue
            tot=sum(g_indicator(t) for t in p)
            best_odd=max(best_odd,tot)
        marg_all.append(g_indicator(s)-best_all)
        marg_odd.append(g_indicator(s)-best_odd if best_odd>-1e17 else np.nan)
    fig,ax=plt.subplots(figsize=(7,3.8))
    ax.axhline(0,color='k',lw=0.8)
    ax.plot(list(ss),marg_all,'s-',ms=4,label='unrestricted partitions: worst margin')
    ax.plot(list(ss),marg_odd,'o-',ms=4,label='coherent (all-odd) partitions: worst margin')
    for s,m in zip(ss,marg_all):
        if m<0: ax.annotate(f's={s}',xy=(s,m),xytext=(s,m-6),ha='center',color='r',fontsize=8)
    ax.set_xlabel('$s=R^2$'); ax.set_ylabel(r'$\Delta V(s)-\max_{\rm split}\sum\Delta V(s_a)$')
    ax.set_title('Exception-free splitting theorem: coherence removes all four violations')
    ax.legend()
    save(fig,'figure_paper5_3_splitting_theorem')

# ============ PAPER 6 ============
def fig6_1():
    pts=[]
    for s in range(3,16,2):
        for p in partitions(s,odd_only=True):
            if len(p)<2 or any(t%2==0 for t in p): continue
            pts.append((1.0/s, sum(1.0/t for t in p)))
    pts=np.array(pts)
    fig,ax=plt.subplots(figsize=(5.4,4.4))
    ax.scatter(pts[:,0],pts[:,1],s=18,alpha=0.7,label='all coherent splits, odd $s\\leq 15$ (68 cases)')
    x=np.linspace(0,1.05,10)
    ax.plot(x,x,'r--',lw=1,label=r'$\sum\lambda^2$ conserved (diagonal)')
    ax.set_xlabel(r'$1/s$ (before split)'); ax.set_ylabel(r'$\sum_a 1/s_a$ (after split)')
    ax.set_title(r'$\Sigma\lambda^2$ strictly increases in every split')
    ax.legend(fontsize=9)
    save(fig,'figure_paper6_1_lambda_monotonicity')

def fig6_2():
    from collections import defaultdict
    fig,ax=plt.subplots(figsize=(8.4,4))
    shells=list(range(1,26,2))
    bottoms=np.zeros(len(shells))
    cmap=plt.get_cmap('tab20')
    shape_colors={}; ci=0
    for i,m in enumerate(shells):
        d=defaultdict(int)
        for v in shell_cells4(float(m)):
            d[tuple(sorted(np.abs(v)))]+=1
        for sh,cnt in sorted(d.items()):
            if sh not in shape_colors:
                shape_colors[sh]=cmap(ci%20); ci+=1
            ax.bar(m, cnt, bottom=bottoms[i], width=1.4, color=shape_colors[sh], edgecolor='k', linewidth=0.3)
            if cnt>=24:
                ax.text(m, bottoms[i]+cnt/2, str(tuple(int(x) for x in sh)).replace(' ',''),
                        ha='center', va='center', fontsize=6.2, rotation=90)
            bottoms[i]+=cnt
    ax.set_xlabel('shell label $m$ (odd)'); ax.set_ylabel('cells in shell')
    ax.set_title('Shell fine structure: degenerate in energy, split by $B_4$ shape orbits')
    ax.set_xticks(shells)
    save(fig,'figure_paper6_2_B4_shell_structure')

def fig6_3():
    fig,ax=plt.subplots(figsize=(5.6,5))
    q=np.linspace(-2.2,2.2,400)
    ax.plot(q,-q,'k-',lw=2,label=r'constraint $\nu\lambda=1$:  null line $u=q+p=0$')
    for c,col in [(0.5,'C0'),(1.0,'C1'),(2.0,'C2')]:
        qq=np.linspace(np.log(c)-2.0, 2.2, 300)
        ax.plot(qq, np.log(c)-qq, ls='--', lw=1, color=col,
                label=fr'$\delta_\nu\delta_\lambda={c}$ (boost orbit)' if c!=0.5 else
                      r'$\delta_\nu\delta_\lambda=\frac{1}{2}$ (zero-point minimum)')
    ax.annotate('', xy=(1.6,-1.6), xytext=(0.2,-0.2),
                arrowprops=dict(arrowstyle='->', color='r', lw=1.6))
    ax.text(1.05,-0.78,'cascade = translation\nalong $v=q-p$', color='r', fontsize=9)
    ax.annotate('', xy=(0.55,0.55), xytext=(0.15,0.15),
                arrowprops=dict(arrowstyle='->', color='b', lw=1.4))
    ax.text(0.32,0.62,'$u$: frozen by\nzero-point', color='b', fontsize=9)
    ax.set_xlabel(r'$q=\log\lambda$'); ax.set_ylabel(r'$p=\log\nu$')
    ax.set_title('Null structure of the log-conjugate plane (exact curves)')
    ax.set_aspect('equal'); ax.legend(fontsize=8, loc='lower left')
    save(fig,'figure_paper6_3_null_structure')

# ============ PAPER 7 ============
def cap(m): return shell_count(m)
def allowed_channels(s):
    """coherent (odd parts), capacity-filtered, >=2 parts"""
    caps={}
    out=[]
    for p in partitions(s,odd_only=True):
        if len(p)<2 or any(t%2==0 for t in p): continue
        ok=True
        from collections import Counter
        for m,n in Counter(p).items():
            if m not in caps: caps[m]=cap(m)
            if n>caps[m]: ok=False; break
        if ok: out.append(p)
    return out

def fig7_1():
    ss=list(range(1,26,2))
    counts=[len(allowed_channels(s)) for s in ss]
    fig,ax=plt.subplots(figsize=(7,3.6))
    cols=['#2a9d2a' if c==0 else 'C0' for c in counts]
    ax.bar(ss,counts,color=cols,edgecolor='k',linewidth=0.4)
    for s,c in zip(ss,counts):
        ax.text(s,c+0.4,str(c),ha='center',fontsize=9)
    ax.set_xlabel('$s$'); ax.set_ylabel('allowed decay channels')
    ax.set_title('Stable-species spectrum: $s=1,3,5$ have no channels; threshold at $s=7$')
    ax.set_xticks(ss)
    ax.annotate('absolutely stable', xy=(3,0.3), xytext=(2.0,6.5), color='#2a9d2a',
                arrowprops=dict(arrowstyle='->', color='#2a9d2a'))
    save(fig,'figure_paper7_1_stable_species')

def fig7_2():
    chans=[('(5,3,1)',192),('(3,3,3)',56)]
    sub=[('orthogonal',96,'C0'),('parallel',48,'C9'),('antiparallel',48,'C1'),
         ('tripod',32,'C2'),('antipodal pair',24,'C3')]
    fig,ax=plt.subplots(figsize=(6.4,3.8))
    bottom=0
    for name,w,c in sub[:3]:
        ax.bar(0,w,bottom=bottom,color=c,edgecolor='k',linewidth=0.4,label=f'(5,3,1) {name}: {w}')
        bottom+=w
    bottom=0
    for name,w,c in sub[3:]:
        ax.bar(1,w,bottom=bottom,color=c,edgecolor='k',linewidth=0.4,label=f'(3,3,3) {name}: {w}')
        bottom+=w
    ax.set_xticks([0,1]); ax.set_xticklabels(['(5,3,1)\n192 configs (77.4%)','(3,3,3)\n56 configs (22.6%)'])
    ax.set_ylabel('number of configurations')
    ax.set_title('Unique timeless branching of $s=9$ with relation-class decomposition')
    ax.legend(fontsize=8.5)
    save(fig,'figure_paper7_2_branching')

def fig7_3():
    # complete line tables (norms vs multiplicity) for the 5 relation classes
    classes={
      '(5,3,1) orthogonal':{1:1,2:2,3:4,4:3,8:2},
      '(5,3,1) parallel':{1:2,2:2,4:3,5:2,8:2},
      '(5,3,1) antiparallel':{1:1,2:2,4:1,5:2,8:2},
      '(3,3,3) tripod':{2:6,4:3},
      '(3,3,3) antipodal':{2:4,4:2},
    }
    norms=[1,2,3,4,5,8]
    fig,ax=plt.subplots(figsize=(8.2,3.9))
    w=0.15
    for i,(name,tab) in enumerate(classes.items()):
        xs=[n+(i-2)*w for n in range(len(norms))]
        ys=[tab.get(nm,0) for nm in norms]
        ax.bar(xs,ys,width=w,label=name,edgecolor='k',linewidth=0.3)
    ax.set_xticks(range(len(norms))); ax.set_xticklabels([f'$|n|^2={nm}$' for nm in norms])
    ax.set_ylabel('line multiplicity')
    ax.set_title('Power-spectrum line tables separate all five relation classes (unit-record readable)')
    ax.legend(fontsize=8.2)
    save(fig,'figure_paper7_3_line_tables')

# ============ PAPER 8 ============
def fig8_1():
    Ss=list(range(3,26))
    single=[]; bestsplit=[]
    for S in Ss:
        if S%2==1: single.append(N0s(S))
        else: single.append(np.nan)
        best=0
        for p in partitions(S,odd_only=True):
            if any(t%2==0 for t in p): continue
            from collections import Counter
            ok=all(n<=cap(m) for m,n in Counter(p).items())
            if not ok: continue
            tot=sum(N0s(t) for t in p)
            if len(p)>=2: best=max(best,tot)
        bestsplit.append(best)
    fig,ax=plt.subplots(figsize=(7,3.8))
    ax.semilogy(Ss,single,'o-',label='single fragment $N_0(\\sqrt{S})$ (odd $S$)')
    ax.semilogy(Ss,bestsplit,'s--',label='best split (capacity-filtered)')
    ax.set_xlabel('$S$'); ax.set_ylabel('total capacity')
    ax.set_title('Shared-container accounting: condensation is strictly optimal')
    ax.legend()
    save(fig,'figure_paper8_1_condensation')

def fig8_2():
    Ss=np.array([9,15,21,25,35,49,81,121,201])
    wstar=[1-(N0s(S)-S)/((pi**2/2)*S*(S-1)) for S in Ss]
    fig,ax=plt.subplots(figsize=(6.6,3.8))
    ax.loglog(Ss,wstar,'o-',label=r'$w^*(S)$ exact (two-pole formula)')
    Sc=np.linspace(8,220,200)
    ax.loglog(Sc,(32/(3*pi))/np.sqrt(Sc),'r--',lw=1,label=r'asymptote $\frac{32}{3\pi}S^{-1/2}\approx 3.40/\sqrt{S}$')
    ax.plot([9],[0.5872],'k*',ms=11,label='S=9 with exclusion rule: 0.5872')
    ax.set_xlabel('$S$'); ax.set_ylabel(r'condensation threshold $w^*$')
    ax.set_title('Jeans-type threshold: larger systems condense at smaller shared weight')
    ax.legend(fontsize=9)
    save(fig,'figure_paper8_2_jeans_threshold')

def fig8_3():
    # cascade k=3, m=14, exact N0, terminal s'=3 (stage j=13)
    k,m=3,14; S=float(k**m)
    js=np.arange(0,m)  # stop at s'=3 (j=13)
    sp=np.array([k**(m-j) for j in js],dtype=float)
    def corr(s):
        if s<=6561: return (pi**2/2)*s*s/N0s(int(s))
        return 1.0/(1.0-(32/(3*pi))/math.sqrt(s))
    a=np.array([( (8*pi**2/3)*S/s * corr(s) )**0.25 for s in sp])
    t=np.cumsum(1/np.sqrt(sp))
    fig,axs=plt.subplots(1,2,figsize=(10,3.8))
    axs[0].loglog(t,a,'o-',ms=4)
    tt=np.array([t[3],t[10]])
    axs[0].loglog(tt, a[5]*(tt/t[5])**0.5,'r--',lw=1,label='slope 1/2')
    axs[0].set_xlabel('internal time $t=\\sum 1/\\sqrt{s\'}$'); axs[0].set_ylabel('apparent scale factor $a$')
    axs[0].set_title("(a) $a\\propto t^{1/2}$ (exact $N_0$, terminal $s'=3$)")
    axs[0].legend()
    # reservoir: k=3, m=8
    k2,m2=3,8; S2=float(k2**m2); Vc=(8*pi**2/3)*S2*S2
    js2=np.arange(0,m2)  # to s'=3
    occ=np.array([ (k2**j)*N0s(int(k2**(m2-j))) for j in js2 ],dtype=float)
    gap=Vc-occ
    axs[1].stackplot(js2, occ, gap, labels=['occupied mode area $\\sum N_0$','gap (vacuum reservoir) $\\Delta V$'],
                     colors=['C0','#dddddd'], edgecolor='k', linewidths=0.4)
    axs[1].set_xlabel('cascade stage $j$'); axs[1].set_ylabel('area')
    axs[1].set_title('(b) Exact conservation: occupied + gap $= V_c$')
    axs[1].legend(loc='center right', fontsize=9)
    save(fig,'figure_paper8_3_expansion_reservoir')

# ============ PAPER 9 ============
def fig9_1():
    N=4096; P=2.0
    x=(np.arange(N)+0.5)/N*P
    xc=np.minimum(x,P-x)
    box=(xc<0.5).astype(float)
    c=np.fft.rfft(box)/N
    fig,ax=plt.subplots(figsize=(7.2,3.8))
    ax.plot(x,box,'k-',lw=1,label='cell indicator (width 1, period 2)')
    for Nh,col in [(1,'C0'),(3,'C1'),(9,'C2'),(39,'C3')]:
        cc=np.zeros_like(c); cc[0]=c[0]
        for n in range(1,Nh+1,2): cc[n]=c[n]
        f=np.fft.irfft(cc*N,n=N)
        ax.plot(x,f,lw=1,label=f'odd harmonics $\\leq {Nh}$')
    ax.set_xlabel('$x$'); ax.set_ylabel('amplitude')
    ax.set_title('A cell is exactly an odd-harmonic (half-wave) object; truncation gives Gibbs edges')
    ax.legend(fontsize=8.5); ax.set_xlim(0,2)
    save(fig,'figure_paper9_1_square_wave')

def fig9_2():
    ls=np.arange(1,21)
    shift=(9/4)/((ls+1.5)+np.sqrt(ls*(ls+3)))
    fig,axs=plt.subplots(1,2,figsize=(10,3.6))
    axs[0].plot(ls,shift,'o-',ms=4)
    axs[0].axhline(0.5,color='r',ls='--',label='censorship bound = zero-point $\\frac{1}{2}$')
    axs[0].annotate('exact saturation at $\\ell=1$', xy=(1,0.5), xytext=(4,0.45),
                    arrowprops=dict(arrowstyle='->'))
    axs[0].set_xlabel('$\\ell$'); axs[0].set_ylabel(r'anharmonic shift $\Delta\nu_\ell$')
    axs[0].set_title('(a) All curvature shifts lie at or below the resolution half-quantum')
    axs[0].legend(fontsize=9)
    # fidelity decay
    tt=np.linspace(0,8,400)
    for L,col in [(2,'C0'),(4,'C1'),(8,'C2')]:
        w=1.0/(np.arange(1,L+1)+1.5)
        dl=(9/4)/((np.arange(1,L+1)+1.5)+np.sqrt(np.arange(1,L+1)*(np.arange(1,L+1)+3)))
        F=np.abs(np.sum(w[None,:]**2*np.exp(2j*pi*dl[None,:]*tt[:,None]*(2/5)**-1*0+2j*pi*dl[None,:]*tt[:,None]),axis=1))/np.sum(w**2)
        axs[1].plot(tt,F,color=col,label=f'ladder of {L} rungs')
    axs[1].axhline(0.5,color='k',ls=':',lw=0.8)
    axs[1].set_xlabel('time (fundamental periods)'); axs[1].set_ylabel('ladder fidelity $F$')
    axs[1].set_title('(b) Without censorship: dephasing collapse in 3–6 periods')
    axs[1].legend(fontsize=9)
    save(fig,'figure_paper9_2_censorship')

def fig9_3():
    s_vals=[1,3,5,7,9,11,13,17,21,25,33,49,81,121]
    best_eye=[1.64,0.29,0.23,-0.10,0.29,0.33,0.06,0.12,0.03,-0.04,-0.09,-0.07,-0.20,-0.37]
    fig,axs=plt.subplots(1,2,figsize=(10.4,3.9))
    x=np.arange(len(s_vals))
    cols=['#2a9d2a' if s in (1,3,5) else ('#d62728' if e<0 else '#888888') for s,e in zip(s_vals,best_eye)]
    axs[0].bar(x,best_eye,color=cols,edgecolor='k',linewidth=0.4)
    axs[0].axhline(0,color='k',lw=0.8)
    axs[0].set_xticks(x); axs[0].set_xticklabels(s_vals)
    axs[0].set_xlabel('$s$'); axs[0].set_ylabel('best composite-scale eye opening')
    axs[0].set_title('(a) Three bands: certified $\\{1,3,5\\}$ / critical $7\\!-\\!21$ / forbidden $\\geq 25$')
    axs[0].annotate('certified',xy=(1,0.6),color='#2a9d2a')
    axs[0].annotate('forbidden\n(monotone worsening)',xy=(10.0,-0.32),color='#d62728',fontsize=9)
    # s=49 deep probe + control
    nc=[1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31]
    eye49=[-0.2254,-0.2012,-0.0667,-0.1164,-0.1505,-0.2255,-0.2006,-0.2223,-0.1510,-0.1078,-0.0607,-0.0686,-0.0743,-0.0731,-0.0742,-0.0698]
    axs[1].plot(nc,eye49,'o-',ms=4,label='odd (logic-wave) sector: closed at all depths')
    ncf=[1,2,3,4,5]; eyef=[-0.2254,-0.4798,-0.4187,-0.2878,0.1610]
    axs[1].plot(ncf,eyef,'s--',ms=5,color='C2',label='control: even harmonics admitted (opens)')
    axs[1].axhline(0,color='k',lw=0.8)
    axs[1].set_xlabel('band depth $n_c$'); axs[1].set_ylabel('eye opening, $s=49$')
    axs[1].set_title('(b) Closure is specific to the odd sector ($s=49$)')
    axs[1].legend(fontsize=8.6)
    save(fig,'figure_paper9_3_three_bands')

# ============ PAPER 10 ============
def great_arc(u,v,n=60):
    u=np.array(u,dtype=float); v=np.array(v,dtype=float)
    u/=np.linalg.norm(u); v/=np.linalg.norm(v)
    ang=math.acos(np.clip(u@v,-1,1))
    ts=np.linspace(0,1,n)
    w=v-(u@v)*u; w/=np.linalg.norm(w)
    return np.array([math.cos(a)*u+math.sin(a)*w for a in ang*ts])

def fig10_1():
    fig=plt.figure(figsize=(5.8,5.4))
    ax=fig.add_subplot(111,projection='3d')
    uu,vv=np.meshgrid(np.linspace(0,2*pi,40),np.linspace(0,pi,20))
    ax.plot_wireframe(np.cos(uu)*np.sin(vv),np.sin(uu)*np.sin(vv),np.cos(vv),
                      color='lightgray',linewidth=0.3,alpha=0.5)
    e1,e2,e3=np.eye(3)
    for a,b in [(e1,e2),(e2,e3),(e3,e1)]:
        arc=great_arc(a,b)
        ax.plot(arc[:,0],arc[:,1],arc[:,2],'C0-',lw=2.4)
    for v,name in [(e1,'$e_1$'),(e2,'$e_2$'),(e3,'$e_3$')]:
        ax.scatter(*v,color='k',s=24)
        ax.text(*(v*1.13),name,fontsize=12)
    # transported tangent at e1: initial points to e2, final = rotated by pi/2 (points to e3... compute)
    t0=np.array([0,1,0]); t1=np.array([0,0,1])
    ax.quiver(1,0,0,*(0.5*t0),color='g',linewidth=2)
    ax.quiver(1,0,0,*(0.5*t1),color='r',linewidth=2)
    ax.text(1.02,0.5,0.04,'initial frame vector',color='g',fontsize=9)
    ax.text(1.02,0.02,0.55,'after loop: rotated by $\\pi/2$',color='r',fontsize=9)
    ax.set_title('Tripod loop on the direction sphere:\nholonomy angle $=$ geodesic-triangle area $=\\pi/2$ (octant)')
    ax.set_box_aspect([1,1,1]); ax.set_axis_off()
    save(fig,'figure_paper10_1_tripod_holonomy')

def fig10_2():
    fig=plt.figure(figsize=(5.8,5.4))
    ax=fig.add_subplot(111,projection='3d')
    uu,vv=np.meshgrid(np.linspace(0,2*pi,40),np.linspace(0,pi,20))
    ax.plot_wireframe(np.cos(uu)*np.sin(vv),np.sin(uu)*np.sin(vv),np.cos(vv),
                      color='lightgray',linewidth=0.3,alpha=0.5)
    a=np.array([2,1,0])/sqrt(5); b=np.array([1,2,0])/sqrt(5); c=np.array([0,0,1.0])
    for u,v in [(a,b),(b,c),(c,a)]:
        arc=great_arc(u,v)
        ax.plot(arc[:,0],arc[:,1],arc[:,2],'C1-',lw=2.4)
    for v,name in [(a,'$(2,1,0)/\\sqrt{5}$'),(b,'$(1,2,0)/\\sqrt{5}$'),(c,'$e_3$')]:
        ax.scatter(*v,color='k',s=24); ax.text(*(v*1.12),name,fontsize=9)
    ax.set_title('Second exact holonomy: spherical excess $=\\arccos(4/5)$\n(irrational multiple of $\\pi$ $\\Rightarrow$ infinite holonomy group)')
    ax.set_box_aspect([1,1,1]); ax.set_axis_off()
    save(fig,'figure_paper10_2_holonomy_arccos45')

def fig10_3():
    from collections import defaultdict
    fig,ax=plt.subplots(figsize=(7.2,4))
    for m in range(1,14,2):
        d=defaultdict(int)
        for v in shell_cells4(float(m)):
            d[(int((v**2).sum()), 1 if int(np.abs(v).sum())%2 else 0, tuple(sorted(np.abs(v))))]+=1
        for (k2,odd,sh),cnt in d.items():
            ax.scatter(m,k2,s=cnt*6,color='C3' if odd else 'C0',alpha=0.8,edgecolors='k',linewidths=0.4)
            ax.annotate(str(tuple(int(x) for x in sh)).replace(' ',''),xy=(m,k2),xytext=(m+0.18,k2+0.12),fontsize=6.4)
    ax.scatter([],[],color='C3',label=r'$\varepsilon=-1$ ($\Sigma|k|$ odd)')
    ax.scatter([],[],color='C0',label=r'$\varepsilon=+1$ ($\Sigma|k|$ even)')
    ax.set_xlabel('shell label $m$'); ax.set_ylabel(r'bare norm $|k|^2$')
    ax.set_title('Excitation parity $\\varepsilon(k)=(-1)^{\\Sigma|k|}$: a gauge-invariant $Z_2$ finer than shells')
    ax.legend(fontsize=9)
    save(fig,'figure_paper10_3_excitation_parity')

# ============ PAPER 11 ============
def fig11_1():
    d=np.linspace(0,9,400)
    f=(np.sqrt(d)-1)**2/2
    fig,ax=plt.subplots(figsize=(6.8,3.9))
    ax.plot(d,f,'C0-',lw=1.8,label=r'$\Delta\nu_1(d)=(\sqrt{d}-1)^2/2$')
    ax.axhline(0.5,color='r',ls='--',lw=1,label='censorship bound $1/2$')
    ax.axvspan(0,4,color='g',alpha=0.10)
    ax.plot([0,4],[0.5,0.5],'ko',ms=6)
    di=np.arange(1,10); ax.plot(di,(np.sqrt(di)-1)**2/2,'ks',ms=4)
    ax.annotate('survival interval $[0,4]$',xy=(1.6,0.55),color='g')
    ax.annotate('saturation only\nat both endpoints',xy=(4,0.5),xytext=(5.4,0.30),
                arrowprops=dict(arrowstyle='->'))
    ax.set_xlabel('dimension $d$'); ax.set_ylabel('maximal anharmonic shift')
    ax.set_title('Both-endpoint theorem: stable composites exist only for $d\\in[0,4]$')
    ax.legend(fontsize=9); ax.set_ylim(0,2.1)
    save(fig,'figure_paper11_1_survival_interval')

def fig11_2():
    fig,axs=plt.subplots(1,3,figsize=(11.4,3.4))
    ds=np.arange(1,9)
    mods=ds%8
    cols=['#2a9d2a' if d==4 else ('#888888' if d!=8 else 'C1') for d in ds]
    axs[0].bar(ds,mods,color=cols,edgecolor='k',linewidth=0.4)
    axs[0].set_xlabel('$d$'); axs[0].set_ylabel(r'$4s$ mod $8$')
    axs[0].set_title('(a) mod-8 selector: odd integer\nlabels only at $d=4$')
    axs[0].annotate('d=4: integer odd labels',xy=(4,4.1),color='#2a9d2a',fontsize=8.5,ha='center')
    axs[0].annotate('d=8: integer even labels\n(no single states)',xy=(8,0.4),color='C1',fontsize=8,ha='right')
    dd=np.arange(0,9)
    axs[1].plot(dd,dd/4,'o-',label=r'seed value $s_{\min}(d)=d/4$')
    axs[1].axhline(1,color='r',ls='--',label=r'self-dual point $s=1$ ($\nu=\lambda=1$)')
    axs[1].plot([4],[1],'k*',ms=12)
    axs[1].set_xlabel('$d$'); axs[1].set_title('(b) Self-dual seed selects $d=4$\nfrom below')
    axs[1].legend(fontsize=8.5)
    Vd=lambda d: pi**(d/2)/math.gamma(d/2+1)
    gains=[Vd(d+1)/Vd(d) for d in range(1,8)]
    axs[2].bar(range(1,8),gains,color='C0',edgecolor='k',linewidth=0.4)
    axs[2].axhline(1,color='k',lw=0.8)
    axs[2].set_xlabel(r'step $d\to d+1$'); axs[2].set_ylabel(r'$V_{d+1}/V_d$')
    axs[2].set_title('(c) Counting drives growth:\nevery step gains capacity ($\\times\\sqrt{S}$)')
    save(fig,'figure_paper11_2_selectors')

# ============ PAPER 12 ============
def fig12_1():
    fig,axs=plt.subplots(1,3,figsize=(10.8,3.3),sharey=True)
    for ax,s in zip(axs,[1,3,5]):
        R=sqrt(s); P=2*R
        kmax={1:0,3:1,5:1}[s]
        W=2*kmax+1
        ms=np.arange(1,14)
        amp=np.abs(np.sin(pi*ms*W/P)/(pi*ms*W/P))
        odd=ms%2==1
        ax.stem(ms[odd]/P, amp[odd], basefmt=' ', linefmt='C0-', markerfmt='C0o',
                label='odd lines')
        if np.any(~odd):
            ax.stem(ms[~odd]/P, amp[~odd], basefmt=' ', linefmt='C3-', markerfmt='C3s',
                    label='even contamination')
        ax.axvline(1/(2*R),color='g',ls=':',lw=1)
        ax.set_title(f'$s={s}$:  $\\nu_0=1/(2\\sqrt{{{s}}})={1/(2*R):.3f}$\ncomb spacing $=1/\\sqrt{{{s}}}$')
        ax.set_xlabel(r'frequency $\nu$ (cell units)')
        if s==1: ax.set_ylabel('relative line amplitude')
        ax.legend(fontsize=8)
    fig.suptitle('External face of the stable species: odd combs, identity $\\nu_0 R\'=1/2$ (axial profiles, exact sinc amplitudes)',y=1.04)
    save(fig,'figure_paper12_1_species_combs')

def fig12_2():
    cells=[np.array(k) for k in itertools.product(range(-2,3),repeat=4)
           if sum((abs(t)+0.5)**2 for t in k)==9.0]
    from collections import Counter
    dd=Counter()
    for i in range(len(cells)):
        for j in range(i+1,len(cells)):
            d=cells[i]-cells[j]
            dd[int((d*d).sum())]+=1
    fig,ax=plt.subplots(figsize=(7.4,3.8))
    xs=sorted(dd)
    cols=['C0' if x%2==0 else 'C2' for x in xs]
    ax.bar(xs,[dd[x] for x in xs],color=cols,edgecolor='k',linewidth=0.4)
    for x in [1,5,9,13,17]:
        ax.plot([x],[4],'rx',ms=9,mew=2)
        ax.text(x,9,'forbidden',rotation=90,fontsize=7.2,ha='center',color='r')
    ax.set_xlabel(r'twin separation $|\Delta|^2$'); ax.set_ylabel('number of configurations')
    ax.set_title('Quantised twin separations with mod-4 selection rule '
                 r'($|\Delta|^2\equiv 1\ (\mathrm{mod}\ 4)$ forbidden)')
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color='C0',label='same-orbit pairs (even $|\\Delta|^2$)'),
                       Patch(color='C2',label=r'mixed-orbit pairs ($|\Delta|^2\equiv 3$ mod 4)')],fontsize=8.6)
    save(fig,'figure_paper12_2_separation_spectrum')

def fig12_3():
    def line_set(occ):
        lines={}
        occ=[np.array(k) for k in occ]
        for i in range(len(occ)):
            for j in range(len(occ)):
                if i==j: continue
                for sgn in (1,-1):
                    v=occ[i]+sgn*occ[j]
                    if not np.any(v): continue
                    t=tuple(v) if tuple(v)>=tuple(-v) else tuple(-v)
                    lines[t]=int(np.dot(v,v))
        return lines
    before=line_set([(2,1,0,0),(1,1,1,1)])
    after=line_set([(1,1,0,0),(0,0,1,0),(0,0,0,0),(1,1,1,1)])
    fig,axs=plt.subplots(1,2,figsize=(9.6,3.5),sharey=True)
    from collections import Counter
    for ax,L,title,dc in [(axs[0],before,'before decay: 2 lines, DC $=2$',2),
                          (axs[1],after,'after decay $9\\to(5,3,1)$: 9 lines, DC $=4$',4)]:
        cn=Counter(L.values())
        xs=sorted(cn)
        ax.bar(xs,[cn[x] for x in xs],color='C0',edgecolor='k',linewidth=0.4,width=0.5)
        ax.bar([0],[dc],color='C1',edgecolor='k',linewidth=0.4,width=0.5)
        ax.text(0,dc+0.1,'DC',ha='center',fontsize=9)
        ax.set_xlabel(r'line norm $|n|^2$'); ax.set_title(title)
    axs[0].set_ylabel('line count')
    fig.suptitle('Measurement as appearance: the record before and after the first decay '
                 r'($\Sigma\lambda^2:\ 2/9\to74/45$, strictly increasing)',y=1.04)
    save(fig,'figure_paper12_3_measurement')

ALL=[fig5_1,fig5_2,fig5_3,fig6_1,fig6_2,fig6_3,fig7_1,fig7_2,fig7_3,
     fig8_1,fig8_2,fig8_3,fig9_1,fig9_2,fig9_3,fig10_1,fig10_2,fig10_3,
     fig11_1,fig11_2,fig12_1,fig12_2,fig12_3]
fails=[]
for f in ALL:
    try:
        f()
    except Exception as e:
        fails.append((f.__name__, repr(e)))
        traceback.print_exc()
print("FAILED:", fails if fails else "none")
