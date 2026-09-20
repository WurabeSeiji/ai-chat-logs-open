import mpmath as mp, csv
mp.mp.dps=500

def recurrence_q_dq(r,N):
    S=r; Q=r*r; dS=mp.mpf(1); dQ=2*r
    for _ in range(1,N):
        Snew=S*(1-Q)
        dSnew=dS*(1-Q)-S*dQ
        Qnew=Q+S*S*Q*Q
        dQnew=dQ+2*S*dS*Q*Q+2*S*S*Q*dQ
        S,Q,dS,dQ=Snew,Qnew,dSnew,dQnew
    return Q,dQ

def newton_root(N, guess):
    x=mp.mpf(guess)
    for _ in range(100):
        Q,dQ=recurrence_q_dq(x,N)
        dx=(Q-1)/dQ
        x-=dx
        if abs(dx)<mp.mpf('1e-450'): break
    return x

def literal_generate(r, upto):
    zs=[mp.mpf(r)]
    for N in range(1,upto):
        total=mp.mpf('0')
        for zi in zs:
            for zj in zs:
                total += zi*zj*zj  # real case: conj(zi)=zi
        zs.append(-total)
    return zs

roots=[]
guess=mp.mpf('1')
for N in range(1,8):
    if N==1: r=mp.mpf(1)
    else:
        if N==2: guess=mp.mpf('0.826')
        else: guess=roots[-1]
        r=newton_root(N,guess)
    roots.append(r); guess=r

with open('/mnt/data/finite_closure_roots_pairwise.csv','w',newline='') as f:
    w=csv.writer(f)
    w.writerow(['N','r_N','Q_N_literal_minus_1','S_N_literal','z_Nplus1_literal','S_Nplus1_literal','z_Nplus2_literal'])
    for N,r in enumerate(roots,1):
        zs=literal_generate(r,N+2)
        QN=sum(z*z for z in zs[:N])
        SN=sum(zs[:N])
        zNp1=zs[N]
        SNp1=sum(zs[:N+1])
        zNp2=zs[N+1]
        w.writerow([N,mp.nstr(r,180),mp.nstr(QN-1,20),mp.nstr(SN,80),mp.nstr(zNp1,80),mp.nstr(SNp1,20),mp.nstr(zNp2,20)])

print(open('/mnt/data/finite_closure_roots_pairwise.csv').read())
