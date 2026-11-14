# 分解能のゲージ（区画 36°）での読出し：ω を連続に動かし、最大の山の位置と、各粗い波の区画での高さを読む
import numpy as np
def S10(x): return np.sum([np.cos(n*x) for n in range(1,11)],axis=0)
phi_deg=np.linspace(-180,180,36001); phi=np.deg2rad(phi_deg); f=S10(10*phi)
def Y2(w,amp):
    r=np.deg2rad(w); return (f*(S10(phi-r)+amp*S10(phi+r)))**2
def tooth(Y,c):  # 区画の歯 c（36°の倍数）での値
    return Y[np.argmin(np.abs(phi_deg-c))]
for amp,name in [(1.0,'対称'),(np.sqrt(0.5),'非対称 0.707')]:
    print(f'== {name}')
    print('   ω     最大の山の位置   最大値     波1の最寄りの歯(値)   波2の最寄りの歯(値)   比 波1/波2')
    for w in np.arange(-72,-35+0.1,3):
        Y=Y2(w,amp); k=Y.argmax()
        c1=36*np.round(w/36); c2=36*np.round(-w/36)
        t1,t2=tooth(Y,c1),tooth(Y,c2)
        print(f'{w:6.1f}  {phi_deg[k]:8.2f}°   {Y[k]:9.1f}   {c1:6.0f}°({t1:8.1f})   {c2:6.0f}°({t2:8.1f})   {t1/t2 if t2>1e-9 else float("nan"):6.3f}')
