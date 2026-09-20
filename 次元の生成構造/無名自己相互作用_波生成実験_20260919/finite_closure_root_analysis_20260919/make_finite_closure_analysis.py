import mpmath as mp, csv
mp.mp.dps=1200

def state_derivative(r,N):
    S=r; Q=r*r; dS=mp.mpf(1); dQ=2*r
    for _ in range(1,N):
        Snew=S*(1-Q)
        dSnew=dS*(1-Q)-S*dQ
        Qnew=Q+S*S*Q*Q
        dQnew=dQ+2*S*dS*Q*Q+2*S*S*Q*dQ
        S,Q,dS,dQ=Snew,Qnew,dSnew,dQnew
    return S,Q,dS,dQ

def rootN(N, guess):
    if N==1: return mp.mpf(1)
    x=mp.mpf(guess)
    for _ in range(100):
        S,Q,dS,dQ=state_derivative(x,N)
        dx=(Q-1)/dQ
        x-=dx
        if abs(dx) < mp.mpf('1e-1100'):
            break
    return x

roots=[mp.mpf(1)]
for N in range(2,10):
    if N==2: guess=mp.mpf('0.826')
    elif N==3: guess=mp.mpf('0.80817')
    else:
        d=roots[-2]-roots[-1]
        guess=roots[-1]-mp.mpf('14')*d**3
    roots.append(rootN(N,guess))

with open('/mnt/data/finite_closure_root_sequence.csv','w',newline='') as f:
    w=csv.writer(f)
    w.writerow(['N','r_N','x_N=r_N^2','difference_r_Nminus1_minus_r_N','log10_difference','cubic_ratio_dN_over_dPrev_cubed'])
    ds=[]
    for i,r in enumerate(roots):
        N=i+1
        if i==0:
            d=''; ld=''; cr=''
        else:
            dv=roots[i-1]-r; ds.append(dv)
            d=mp.nstr(dv,350); ld=mp.nstr(mp.log10(abs(dv)),80)
            cr=mp.nstr(ds[-1]/ds[-2]**3,80) if len(ds)>=2 else ''
        w.writerow([N,mp.nstr(r,350),mp.nstr(r*r,350),d,ld,cr])

# critical-orbit diagnostics using r_9 as accumulation approximation
r=roots[-1]
S=r; Q=r*r
with open('/mnt/data/critical_orbit_diagnostics.csv','w',newline='') as f:
    w=csv.writer(f); w.writerow(['n','S_n','Q_n','delta_n=1-Q_n','S_n_squared','delta_minus_S2','z_nplus1'])
    for n in range(1,9):
        d=1-Q
        w.writerow([n,mp.nstr(S,200),mp.nstr(Q,200),mp.nstr(d,200),mp.nstr(S*S,200),mp.nstr(d-S*S,200),mp.nstr(-S*Q,200)])
        z=-S*Q; S=S+z; Q=Q+z*z

with open('/mnt/data/finite_closure_analysis_ja.md','w') as f:
    f.write('# 有限世代閉鎖根と集積点の解析\n\n')
    f.write('更新則は元の全順序対則 $T_{ij}=(\\bar z_i z_j)z_j$。実数初期値では、全 $N^2$ 項の和を代数的にまとめると診断変数 $S_N=\\sum z_k$, $Q_N=\\sum z_k^2$ に対し\n\n')
    f.write('$$z_{N+1}=-S_NQ_N,$$\n')
    f.write('$$S_{N+1}=S_N(1-Q_N),$$\n')
    f.write('$$Q_{N+1}=Q_N+S_N^2Q_N^2.$$\n\n')
    f.write('したがって $Q_N=1$ なら $S_{N+1}=0$ となり、次の新生波 $z_{N+2}=0$、以後すべて 0 になる。これを有限世代閉鎖根 $r_N$ と呼ぶ。\n\n')
    f.write('## 根系列\n\n')
    for i,r in enumerate(roots[:8],1):
        f.write(f'- N={i}: `{mp.nstr(r,120)}`\n')
    f.write('\n根は単調減少し、差 $d_N=r_{N-1}-r_N$ は高次でほぼ三次収束する。\n\n')
    for i in range(3,9):
        dcur=roots[i-1]-roots[i]
        dprev=roots[i-2]-roots[i-1]
        f.write(f'- N={i+1}: $d_N/d_{{N-1}}^3 \\approx {mp.nstr(dcur/dprev**3,30)}$\n')
    f.write('\n極限値（N=9根で少なくとも約834桁まで近似）は\n\n')
    f.write('`'+mp.nstr(roots[-1],350)+'`\n\n')
    f.write('## 多項式次数\n\n')
    f.write('$x=r^2$, $P_N=S_N/r$ とおくと $P_1=1, Q_1=x$、\n')
    f.write('$$P_{N+1}=P_N(1-Q_N),\\qquad Q_{N+1}=Q_N+xP_N^2Q_N^2.$$\n')
    f.write('よって $\\deg Q_N=3^{N-1}$、$\\deg P_N=(3^{N-1}-1)/2$。最初の根条件は\n')
    f.write('$$N=1:\ x-1=0,$$\n')
    f.write('$$N=2:\ x^3+x-1=0,$$\n')
    f.write('$$N=3:\ x^9-2x^8+3x^7-4x^6+3x^5-2x^4+2x^3+x-1=0.$$\n\n')
    f.write('## α有限位数根との比較\n\n')
    f.write('共通点は、連続パラメータ掃引の中に有限回で残差が厳密にゼロになる離散根があり、そのため数値精度を上げるほど谷が深く狭くなる点。相違点は、α側が $U^n=I$ という周期回帰であるのに対し、今回はこちらが $Q_N=1\\Rightarrow S_{N+1}=0$ という吸収的閉鎖である点。したがって数学的機構は同一とはまだ言えないが、「有限閉鎖根がゼロ残差の特異谷を作る」というクラスは共通している。\n')
