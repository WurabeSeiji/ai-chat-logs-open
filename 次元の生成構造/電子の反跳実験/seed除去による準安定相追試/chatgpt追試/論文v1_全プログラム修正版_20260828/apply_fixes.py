# -*- coding: utf-8 -*-
"""全プログラムに振幅問題・回転問題の修正を機械的に適用する（決定論的、置換件数を検証）。
 FIX1 親の振幅正規化を除去（make_parent の v = v/‖v‖）
 FIX2 初期化の外部 seed と正規化を除去（Z = v）、乱数初期状態の正規化も除去
 FIX3 Cayley 変換 → 厳密な線形回転 exp(ANGLE·K)、ANGLE = 2π/144（旧 GAMMA=tan(π/144) と同じ名目刻み）
 FIX4 位相のみ生成子 → 振幅込み生成子 K_ij = Im(z̄_i z_j)（set_state(z)：c=Re z, s=Im z）
差分は results/fix_patches.diff（original との diff -ru）に保存。"""
import os, re, sys, glob
HERE=os.path.dirname(os.path.abspath(__file__)); FX=os.path.join(HERE,"fixed"); log=[]
def sub(path, old, new, count=1, regex=False):
    s=open(path,encoding="utf-8").read()
    n = len(re.findall(old,s)) if regex else s.count(old)
    if n!=count: raise SystemExit(f"[apply_fixes] {os.path.relpath(path,FX)}: expected {count} match of {old[:60]!r}, found {n}")
    s = re.sub(old,new,s) if regex else s.replace(old,new)
    open(path,"w",encoding="utf-8").write(s); log.append(f"{os.path.relpath(path,FX)}: {old.strip()[:70]!r} -> {new.strip()[:70]!r} x{count}")
# ---------------------------------------------------------------- 1. 共通エンジン（run_n_scaling_lowrank_v1_*.py, original_RAW_K_reference.py）
ENGINES=sorted(glob.glob(os.path.join(FX,"**","run_n_scaling_lowrank_v1_*.py"),recursive=True)+glob.glob(os.path.join(FX,"**","original_RAW_K_reference.py"),recursive=True))
SET_STATE='''    def set_theta(self, theta):
        """位相のみ（|z|=1）の生成子。親の固有モード反復でのみ使用。"""
        self.set_state(np.exp(1j * np.asarray(theta, dtype=float)))

    def set_state(self, z):
        """FIX4: 振幅込み生成子 K_ij = Im(conj(z_i) z_j) = |z_i||z_j| sin(θ_j-θ_i)（隣接辺）。c=Re z, s=Im z。"""
        n = self.n
        z = np.asarray(z, dtype=complex)
        self.c = np.real(z).copy()
        self.s = np.imag(z).copy()
        CT = np.zeros((n, n))
        ST = np.zeros((n, n))
        CT[self.ea, self.eb] = self.c
        CT[self.eb, self.ea] = self.c
        ST[self.ea, self.eb] = self.s
        ST[self.eb, self.ea] = self.s
'''
OLD_SET_THETA='''    def set_theta(self, theta):
        n = self.n
        self.c = np.cos(theta)
        self.s = np.sin(theta)
        T = np.zeros((n, n))
        T[self.ea, self.eb] = theta
        T[self.eb, self.ea] = theta
        CT = np.cos(T)
        ST = np.sin(T)
        np.fill_diagonal(CT, 0.0)
        np.fill_diagonal(ST, 0.0)
'''
STEP_RAW='''    def dense_K(self):
        """現在の状態（set_state）の振幅込み生成子 K を密行列で返す（M×M、実反対称）。"""
        return np.column_stack([self.kmatvec(e) for e in np.eye(self.m)])

    def linear_rotation_step(self, z, sigma):
        """FIX3: z ← exp(ANGLE·K) z。厳密な線形回転（実直交）。K/σ 正規化なし。"""
        K = self.dense_K()
        w, V = np.linalg.eigh(1j * K)
        return V @ (np.exp(-1j * ANGLE * w) * (V.conj().T @ z))
'''
STEP_NORM='''    def dense_K(self):
        """現在の状態（set_state）の振幅込み生成子 K を密行列で返す（M×M、実反対称）。"""
        return np.column_stack([self.kmatvec(e) for e in np.eye(self.m)])

    def linear_rotation_step(self, z, sigma):
        """FIX3: z ← exp((ANGLE/σ)·K) z。厳密な線形回転、K/σ 正規化あり（比較用ブランチ）。"""
        K = self.dense_K()
        w, V = np.linalg.eigh(1j * K)
        return V @ (np.exp(-1j * (ANGLE / sigma) * w) * (V.conj().T @ z))
'''
for e in ENGINES:
    s=open(e,encoding="utf-8").read(); normalized = "gn = GAMMA / sigma" in s
    sub(e,'GAMMA = math.tan(math.pi / 144.0)\n','GAMMA = math.tan(math.pi / 144.0)  # 旧 Cayley 刻み（記録用、力学には使わない）\nANGLE = 2.0 * math.pi / 144.0    # FIX3: 線形回転の刻み角\n')
    sub(e,OLD_SET_THETA,SET_STATE)
    old_step = ('    def cayley_step(self, z, sigma):\n        """z ← (I-γK̃)^{-1}(I+γK̃) z, K̃ = K/σ。Woodbury で O(N^3)。"""\n        gn = GAMMA / sigma\n        r = z + gn * self.kmatvec(z)\n        A2 = (sigma / GAMMA) * self.J + self.G\n        rhs = self.wt(r)\n        y = np.linalg.solve(A2, rhs)\n        return r - self.w(y)\n') if normalized else \
               ('    def cayley_step(self, z, sigma):\n        """z ← (I-γK)^{-1}(I+γK) z。K/σ 正規化なし。Woodbury で O(N^3)。"""\n        gn = GAMMA\n        r = z + gn * self.kmatvec(z)\n        A2 = (1.0 / GAMMA) * self.J + self.G\n        rhs = self.wt(r)\n        y = np.linalg.solve(A2, rhs)\n        return r - self.w(y)\n')
    sub(e,old_step,STEP_NORM if normalized else STEP_RAW)
    sub(e,'            v = sys_lr.w(EV[:, idx].astype(complex))\n            v = v / np.linalg.norm(v)\n','            v = sys_lr.w(EV[:, idx].astype(complex))  # FIX1: 振幅正規化を除去\n')
    sub(e,'    g = zero_closure_kernel_seed(sys_lr, rng)\n    Z = v + delta * g\n    Z = Z / np.linalg.norm(Z)\n','    Z = v.copy()  # FIX2: 外部 seed と正規化を除去（seedless、親そのもの）\n')
    sub(e,'    Z = X + 1j * Y\n    return Z / np.linalg.norm(Z)\n','    Z = X + 1j * Y\n    return Z  # FIX1/2: 正規化を除去\n')
    # validate_against_dense: 密行列側も振幅込み K と線形回転で比較
    sub(e,'    theta0 = np.angle(Z)\n    sys_lr.set_theta(theta0)\n    Kd = A * np.sin(theta0[None, :] - theta0[:, None])\n','    sys_lr.set_state(Z)\n    Kd = A * np.imag(np.conj(Z)[:, None] * Z[None, :])  # FIX4: 振幅込み\n')
    sub(e,'''        sys_lr.set_theta(np.angle(Z_lr))
        sig = sys_lr.sigma_spectrum()[0]
        Z_lr = sys_lr.cayley_step(Z_lr, sig)

        th = np.angle(Z_d)
        Kd = A * np.sin(th[None, :] - th[:, None])
        sd = np.linalg.norm(Kd, 2)
        Kn = Kd / sd
        I = np.eye(m)
        Z_d = np.linalg.solve(I - GAMMA * Kn, (I + GAMMA * Kn) @ Z_d)
''','''        sys_lr.set_state(Z_lr)
        sig = sys_lr.sigma_spectrum()[0]
        Z_lr = sys_lr.linear_rotation_step(Z_lr, sig)

        Kd = A * np.imag(np.conj(Z_d)[:, None] * Z_d[None, :])  # FIX4
        wd, Vd = np.linalg.eigh(1j * Kd)
        ang = (ANGLE / np.linalg.norm(Kd, 2)) if %s else ANGLE
        Z_d = Vd @ (np.exp(-1j * ang * wd) * (Vd.conj().T @ Z_d))  # FIX3
''' % ("True" if normalized else "False"))
    sub(e,'        sys_lr.set_theta(np.angle(Z))\n        sig_est, wp = sys_lr.sigma_max_power(wp)\n        Z = sys_lr.cayley_step(Z, sig_est)\n','        sys_lr.set_state(Z)  # FIX4\n        sig_est, wp = sys_lr.sigma_max_power(wp)\n        Z = sys_lr.linear_rotation_step(Z, sig_est)  # FIX3\n')
    sub(e,'        sys_lr.set_theta(np.angle(Z))\n        sig_est, wp = sys_lr.sigma_max_power(wp)\n        if t % sub == 0:','        sys_lr.set_state(Z)  # FIX4\n        sig_est, wp = sys_lr.sigma_max_power(wp)\n        if t % sub == 0:')
    sub(e,'            Z = sys_lr.cayley_step(Z, sig_est)\n    t_run = time.time() - t0\n\n    sys_lr.set_theta(np.angle(Z))\n','            Z = sys_lr.linear_rotation_step(Z, sig_est)  # FIX3\n    t_run = time.time() - t0\n\n    sys_lr.set_state(Z)  # FIX4\n')
    sub(e,'    out = {"n": n, "m": m, "seed": seed, "gamma": GAMMA}','    out = {"n": n, "m": m, "seed": seed, "gamma": GAMMA, "angle": ANGLE}')
    s=open(e,encoding="utf-8").read(); assert "cayley_step" not in s and "set_theta(np.angle(Z" not in s, e
# ---------------------------------------------------------------- 2. 呼び出し側（複製エンジンを import するスクリプト）
CALLERS=[p for p in glob.glob(os.path.join(FX,"**","*.py"),recursive=True) if p not in ENGINES and not p.endswith(("run_followup_experiments.py","run_moduli_sweep_fast.py","analyze_followup.py"))]
for p in CALLERS:
    s=open(p,encoding="utf-8").read()
    for pat,rep in [(r'(\w+)\.set_theta\(np\.angle\((\w+)\)\)', r'\1.set_state(\2)  # FIX4'), (r'\.cayley_step\(', '.linear_rotation_step(')]:
        n=len(re.findall(pat,s))
        if n: sub(p,pat,rep,count=n,regex=True)
    s=open(p,encoding="utf-8").read()
    if "set_theta(theta)" in s and "run_N5_physical_phase_step_test" in p: sub(p,'        sys.set_theta(theta)\n','        sys.set_state(Z)  # FIX4\n')
    if "run_artifact_comparison" in p:
        sub(p,'    g=norm.zero_closure_kernel_seed(sys0,rng)\n    Z0=v+DELTA*g\n    Z0=Z0/np.linalg.norm(Z0)\n','    Z0=v.copy()  # FIX2: 外部 seed と正規化を除去（DELTA は不使用）\n')
        sub(p,'# Branch difference: only cayley_step K/sigma normalization.','# Branch difference: only K/sigma normalization inside the linear rotation exp((ANGLE/sigma) K) vs exp(ANGLE K).')
    if "run_N5_gamma_continuum_test" in p:
        sub(p,'    gamma = math.tan(math.pi / n_den)\n    engine.GAMMA = gamma\n    dphi = 2.0 * math.atan(gamma)\n','    gamma = math.tan(math.pi / n_den)  # 記録用\n    engine.ANGLE = 2.0 * math.pi / n_den  # FIX3: 線形回転の刻み角\n    dphi = engine.ANGLE\n')
    if "run_complex_simplex_decompactification" in p:
        sub(p,"        'N':n,'M':m,'steps':STEPS,'seed':SEED,'gamma':float(eng.GAMMA),\n        'physics':'raw K Cayley; K/sigma normalization removed',","        'N':n,'M':m,'steps':STEPS,'seed':SEED,'gamma':float(eng.GAMMA),'angle':float(eng.ANGLE),\n        'physics':'amplitude-aware K, exact linear rotation exp(ANGLE K), no normalization (FIX1-4)',")
    s=open(p,encoding="utf-8").read(); assert "cayley_step" not in s and "set_theta(np.angle(Z" not in s, p
# ---------------------------------------------------------------- 3. followup（自前エンジン）
f=os.path.join(FX,"N5_dynamics_followup_theorems_and_stability_20260826","followup_dynamics_20260826","run_followup_experiments.py")
sub(f,'GAMMA=math.tan(math.pi/144.0)\n','GAMMA=math.tan(math.pi/144.0)  # 旧 Cayley 刻み（記録用）\nANGLE=2.0*math.pi/144.0        # FIX3: 線形回転の刻み角\n')
sub(f,'''    def set_theta(self,theta):
        n=self.n; self.c=np.cos(theta); self.s=np.sin(theta)
        T=np.zeros((n,n)); T[self.ea,self.eb]=theta; T[self.eb,self.ea]=theta
        CT=np.cos(T); ST=np.sin(T); np.fill_diagonal(CT,0); np.fill_diagonal(ST,0)
''','''    def set_theta(self,theta): self.set_state(np.exp(1j*np.asarray(theta,dtype=float)))  # 位相のみ（親の固有モード反復でのみ使用）
    def set_state(self,z):
        # FIX4: 振幅込み生成子 K_ij = Im(conj(z_i) z_j)。c=Re z, s=Im z
        n=self.n; z=np.asarray(z,dtype=complex); self.c=np.real(z).copy(); self.s=np.imag(z).copy()
        CT=np.zeros((n,n)); ST=np.zeros((n,n)); CT[self.ea,self.eb]=self.c; CT[self.eb,self.ea]=self.c; ST[self.ea,self.eb]=self.s; ST[self.eb,self.ea]=self.s
''')
sub(f,'''    def cayley_step(self,z,sigma_unused=None):
        gn=GAMMA; r=z+gn*self.kmatvec(z); A2=(1.0/GAMMA)*self.J+self.G; rhs=self.wt(r); y=np.linalg.solve(A2,rhs); return r-self.w(y)
''','''    def dense_K(self): return np.column_stack([self.kmatvec(e) for e in np.eye(self.m)])
    def linear_rotation_step(self,z,sigma_unused=None):
        # FIX3: z ← exp(ANGLE·K) z（厳密な線形回転、正規化なし）
        K=self.dense_K(); w,V=np.linalg.eigh(1j*K); return V@(np.exp(-1j*ANGLE*w)*(V.conj().T@z))
''')
sub(f,'            v=sys.w(EV[:,idx].astype(complex)); v=v/np.linalg.norm(v); theta_new=np.angle(v)\n','            v=sys.w(EV[:,idx].astype(complex)); theta_new=np.angle(v)  # FIX1: 振幅正規化を除去\n')
sub(f,'                sys.set_theta(np.angle(Z)); Z=sys.cayley_step(Z)\n','                sys.set_state(Z); Z=sys.linear_rotation_step(Z)  # FIX3/4\n')
sub(f,'    sys.set_theta(np.angle(v)); Fv=sys.cayley_step(v); phase=np.vdot(v,Fv); phase=phase/abs(phase)\n','    sys.set_state(v); Fv=sys.linear_rotation_step(v); phase=np.vdot(v,Fv); phase=phase/abs(phase)  # FIX3/4\n')
sub(f,'        sys.set_theta(np.angle(z)); return sys.cayley_step(z)/phase\n','        sys.set_state(z); return sys.linear_rotation_step(z)/phase  # FIX3/4\n')
sub(f,'        for t in range(5000): sys.set_theta(np.angle(Z)); Z=sys.cayley_step(Z)\n','        for t in range(5000): sys.set_state(Z); Z=sys.linear_rotation_step(Z)  # FIX3/4\n')
s=open(f,encoding="utf-8").read(); assert "cayley_step" not in s and "set_theta(np.angle(Z" not in s and "set_theta(np.angle(z" not in s
m=os.path.join(FX,"N5_dynamics_followup_theorems_and_stability_20260826","followup_dynamics_20260826","run_moduli_sweep_fast.py")
sub(m,'        sysm.set_theta(np.angle(Z)); Z=sysm.cayley_step(Z)\n','        sysm.set_state(Z); Z=sysm.linear_rotation_step(Z)  # FIX3/4\n')
a=os.path.join(FX,"N5_dynamics_followup_theorems_and_stability_20260826","followup_dynamics_20260826","analyze_followup.py")
sub(a,'GAMMA=math.tan(math.pi/144.0)\n','GAMMA=math.tan(math.pi/144.0)\nANGLE=2.0*math.pi/144.0  # FIX3\n')
sub(a,"    if kind=='normalized': inc=np.full(len(d),2*math.atan(GAMMA))\n    else: inc=2*np.arctan(GAMMA*d.sigma_exact.to_numpy())\n","    if kind=='normalized': inc=np.full(len(d),ANGLE)  # FIX3: 線形回転では位相増分は ANGLE（正規化）/ ANGLE·σ（raw）\n    else: inc=ANGLE*d.sigma_exact.to_numpy()\n")
open(os.path.join(HERE,"results","apply_fixes.log"),"w",encoding="utf-8").write("\n".join(log)+"\n"); print(f"applied {len(log)} substitutions in {len(set(l.split(':')[0] for l in log))} files")
