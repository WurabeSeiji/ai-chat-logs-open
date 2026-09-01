import numpy as np, math, csv, json, time
p='../../hm_mp_free_N3_N40_20260901/data/hm_N40/parent_v.npz'; d=np.load(p); z0=d['v'].astype(np.complex128); N=40
ea,eb=np.triu_indices(N,1); M=len(ea); A=np.zeros((M,M))
for e in range(M):
 s=(ea==ea[e])|(ea==eb[e])|(eb==ea[e])|(eb==eb[e]); s[e]=False; A[e,s]=1
pvec=z0.real.copy(); pvec/=np.linalg.norm(pvec); q=z0.imag.copy(); q-=np.dot(q,pvec)*pvec; q/=np.linalg.norm(q)
def met(z):
 h=np.vdot(z,z).real; zp=z-pvec*np.dot(pvec,z)-q*np.dot(q,z); return np.vdot(zp,zp).real/h,h,abs(z@z)/h
z=z0.copy(); rows=[]; dt=2*math.pi/N; t0=time.time(); bad=None
for t in range(151):
 hp,h,cl=met(z); rows.append((t,hp,h,cl,np.max(np.abs(z)),np.all(np.isfinite(z))))
 if not np.all(np.isfinite(z)) or abs(h-np.vdot(z0,z0).real)>1e-10:
  bad=t; break
 if t<150:
  H=A*(z.conj()[:,None]*z[None,:]); herm=np.linalg.norm(H-H.conj().T)/max(np.linalg.norm(H),1e-300)
  w,V=np.linalg.eigh(H); phase=np.exp(-1j*dt*w); zn=V@(phase*(V.conj().T@z));
  if t in (0,1,2,10,50,100): print(t,hp,h,cl,herm,np.linalg.norm(zn)/np.linalg.norm(z)-1)
  z=zn
with open('./diagnostic_0_200.csv','w',newline='') as f:
 cw=csv.writer(f); cw.writerow(['step','Hperp_frac','H_total','closure','max_abs_z','finite']); cw.writerows(rows)
print('DONE',rows[-1],'bad',bad,'sec',time.time()-t0)
