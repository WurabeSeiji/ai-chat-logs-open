#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
電子の候補に倍音 1〜6 を積んだ状態 2 個の、積の円周平均と結合係数。仕様は思考実験2 の §23.3〜§23.5。
octonion_pulse_map_O2_to_O_v1.py と同じフォルダで実行する: python3 octonion_electron_candidates_h6_v1.py
格子点数は 1024 に固定。
"""
import sys, numpy as np
sys.argv=[sys.argv[0],'1024']
sys.path.insert(0,'.')
import octonion_pulse_map_O2_to_O_v1 as m
NG=m.NG; psi=m.psi; mul=m.mul; conj=m.conj
M=6
def zc(o): return np.array([o[0]+1j*o[7], o[1]-1j*o[6], o[2]+1j*o[5], o[3]+1j*o[4]])
def rot(w,th):
    X=np.zeros((NG,8)); X[:,0]=np.cos(w*(psi-th)); X[:,7]=np.sin(w*(psi-th)); return X
def Cst(th,s=-1,M=M):   # (1,e7) 面, 巻数 s*m, m=1..M
    return sum(rot(s*k,th) for k in range(1,M+1))/np.sqrt(M)
def Dst(k,th,s=-1,M=M): # Q_k 面, 巻数 s*m
    V=np.zeros(8); V[k]=1
    return sum(mul(rot(s*j,th),np.broadcast_to(V,(NG,8))) for j in range(1,M+1))/np.sqrt(M)
def avg(X): return X.sum(0)/NG
def fmt(o):
    z=zc(o); return "  ".join("%.3f∠%+4.0f°"%(abs(c),np.degrees(np.angle(c))) if abs(c)>1e-9 else "   0        " for c in z)
def spec(X,tol=1e-9):
    Z=np.fft.fft(np.stack([X[:,0]+1j*X[:,7],X[:,1]-1j*X[:,6],X[:,2]+1j*X[:,5],X[:,3]+1j*X[:,4]],1),axis=0)/NG
    out=[]
    for k in range(4):
        idx=np.where(np.abs(Z[:,k])>tol)[0]; w=sorted(int(i if i<NG//2 else i-NG) for i in idx)
        out.append("%d..%d"%(w[0],w[-1]) if w else "-")
    return out
c=Cst(0); d=Dst(1,0)
for name,x in (("C型",c),("D型(Q1面)",d)):
    dn=(x**2).sum(1)
    print(name,"パワー %.4f  密度|Q|²: 最大 %.3f 最小 %.3f  半値幅 %.1f°  巻数の範囲(4面)"%(m.power(x),dn.max(),dn.min(),360*(dn>dn.max()/2).sum()/NG),spec(x))
print("密度は C型と D型で同じ:",np.allclose((c**2).sum(1),(d**2).sum(1)))
print("\n円周平均を 4 つの複素面 (1,e7) Q1 Q2 Q3 の 大きさ∠位相 で表示")
cases=[("C·C",lambda D:(Cst(0),Cst(D))),("C·C逆巻",lambda D:(Cst(0),Cst(D,+1))),
       ("D1·D1",lambda D:(Dst(1,0),Dst(1,D))),("D1·D1逆巻",lambda D:(Dst(1,0),Dst(1,D,+1))),
       ("D1·D2",lambda D:(Dst(1,0),Dst(2,D))),("C·D1",lambda D:(Cst(0),Dst(1,D)))]
for nm,fn in cases:
    print("\n=====",nm,"=====")
    for Dg in (0,10,30,60,90,120,180):
        a,b=fn(np.radians(Dg))
        print(" Δ=%3d°"%Dg)
        print("   平均 a·b      :",fmt(avg(mul(a,b))),"  巻数",spec(mul(a,b)))
        print("   平均 b·a      :",fmt(avg(mul(b,a))))
        print("   平均 a·conj(b):",fmt(avg(mul(a,conj(b)))))
        print("   平均 conj(a)·b:",fmt(avg(mul(conj(a),b))))
        print("   k_ab 右 <a,a·b>:",fmt(m.inner(a,mul(a,b))))
        print("   k_ab 左 <a,b·a>:",fmt(m.inner(a,mul(b,a))))
        P=a+b; print("   全体のパワー |a+b|² 平均 %.4f   写像後 右 %.4f 左 %.4f"%(m.power(P),m.power(mul(a,P))+m.power(mul(b,P)),m.power(mul(P,a))+m.power(mul(P,b))))
