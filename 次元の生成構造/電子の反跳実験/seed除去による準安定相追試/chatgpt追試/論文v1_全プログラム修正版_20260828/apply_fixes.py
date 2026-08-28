# -*- coding: utf-8 -*-
"""確定した判断（棚卸し文書・判断1〜6）を全プログラムに機械適用する。置換件数を検証し、差分を results/fix_patches.diff に保存。
 A1 親の振幅正規化を除去 / A2 外部seed δg を状態に加えない（Z=v）/ A3 初期正規化なし / A4 動力学の K を振幅込み Im(z̄ᵢzⱼ)（set_state）
 A5 zero_closure_generic の全体正規化を除去 / A6(b) σ は実際に回している K から読む（sigma_spectrum、冪反復は使わない）
 R1 Cayley → exp(ANGLE·K), ANGLE=2π/144 / R2 K/σ 枝の廃止（NORMALIZED エンジンは実行しない、γ掃引は raw）
 R3 (i) 正規化枝依存の再パラメータ化テストは raw のみ (ii) validate を修正後の力学に書換え (iii) dphi=2π/n_den
 S1 zero_closure_kernel_seed を呼ばない（定義は残す）/ S2 sigma_max_power・wp=rng.normal を呼ばない（定義は残す）/ S3 親探索の乱数は残す
 B  KMODE 環境変数で baseline（位相のみ K）/ treatment（振幅込み K）を切替え、同じプログラムを 2 回走らせる
"""
import os, re, glob
HERE=os.path.dirname(os.path.abspath(__file__)); FX=os.path.join(HERE,"fixed"); log=[]
def sub(path, old, new, count=1, regex=False):
    s=open(path,encoding="utf-8").read(); n=len(re.findall(old,s,re.M)) if regex else s.count(old)
    if n!=count: raise SystemExit(f"[apply_fixes] {os.path.relpath(path,FX)}: expected {count} of {old[:70]!r}, found {n}")
    s=re.sub(old,new,s,flags=re.M) if regex else s.replace(old,new); open(path,"w",encoding="utf-8").write(s); log.append(f"{os.path.relpath(path,FX)} | {old.strip()[:60]!r} -> {new.strip()[:60]!r} x{count}")
# ============================================================ 1. 共通エンジン（NORMALIZED_ORIGINAL は R2 廃止により実行しないので無変更）
ENGINES=sorted(p for p in glob.glob(os.path.join(FX,"**","run_n_scaling_lowrank_v1_*.py"),recursive=True)+glob.glob(os.path.join(FX,"**","original_RAW_K_reference.py"),recursive=True) if "NORMALIZED_ORIGINAL" not in p)
HDR='''GAMMA = math.tan(math.pi / 144.0)  # 旧 Cayley 刻み（記録用。力学には使わない）
ANGLE = 2.0 * math.pi / 144.0    # R1: 線形回転 exp(ANGLE·K) の刻み角（124/144 は非本質、旧名目刻みに合わせる）
KMODE = os.environ.get("KMODE", "amplitude")  # B: "amplitude" = 振幅込み K（修正後の力学）, "phase" = 位相のみ K（対照 baseline）
'''
SET_STATE='''    def set_theta(self, theta):
        """位相のみ（|z|=1）の生成子。make_parent の固有モード反復でのみ使用（S3: 親の作り方は無変更）。"""
        self.set_state(np.exp(1j * np.asarray(theta, dtype=float)), force_phase=True)

    def set_state(self, z, force_phase=False):
        """A4: 振幅込み生成子 K_ij = Im(conj(z_i) z_j)（隣接辺）。c=Re z, s=Im z で低ランク表現を組む。
        KMODE="phase"（対照 baseline）では z/|z| を使い、旧来の位相のみ K を線形回転で回す。"""
        n = self.n
        z = np.asarray(z, dtype=complex)
        if force_phase or KMODE == "phase":
            r = np.abs(z)
            z = np.where(r > 0, z / np.where(r > 0, r, 1.0), 0.0)
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
LIN='''    def dense_K(self):
        """現在の状態（set_state）の生成子 K を密行列で返す（M×M、実反対称）。"""
        return np.column_stack([self.kmatvec(e) for e in np.eye(self.m)])

    def linear_rotation_step(self, z, sigma=None):
        """R1: z ← exp(ANGLE·K) z。厳密な線形回転（実直交）。K/σ 正規化なし（R2 廃止）。"""
        K = self.dense_K()
        w, V = np.linalg.eigh(1j * K)
        return V @ (np.exp(-1j * ANGLE * w) * (V.conj().T @ z))

'''
for e in ENGINES:
    s=open(e,encoding="utf-8").read(); normalized="gn = GAMMA / sigma" in s   # abprobe（γ掃引）は K/σ 版だが R2 廃止で raw として扱う
    if "import os" not in s.split("class LowRankSystem")[0]: sub(e,"import numpy as np\n","import os\nimport numpy as np\n")
    sub(e,"GAMMA = math.tan(math.pi / 144.0)\n",HDR)
    sub(e,OLD_SET_THETA,SET_STATE)
    sub(e,"    def cayley_step(self, z, sigma):\n",LIN+"    def cayley_step(self, z, sigma):  # 旧 Cayley 更新（R1 により実行しない。定義は残す）\n")
    sub(e,"            v = sys_lr.w(EV[:, idx].astype(complex))\n            v = v / np.linalg.norm(v)\n","            v = sys_lr.w(EV[:, idx].astype(complex))  # A1: 振幅正規化を除去\n")
    sub(e,"    Z = X + 1j * Y\n    return Z / np.linalg.norm(Z)\n","    Z = X + 1j * Y\n    return Z  # A5: 全体の振幅正規化を除去（|X|=|Y| の相対比は維持）\n")
    # R3(ii) validate：密行列側も同じ力学で
    sub(e,"    theta0 = np.angle(Z)\n    sys_lr.set_theta(theta0)\n    Kd = A * np.sin(theta0[None, :] - theta0[:, None])\n",
          "    sys_lr.set_state(Z)\n    zz = Z / np.abs(Z) if KMODE == \"phase\" else Z\n    Kd = A * np.imag(np.conj(zz)[:, None] * zz[None, :])  # A4/R3(ii): 密行列側も同じ生成子\n")
    sub(e,'''    sig_lr = sys_lr.sigma_spectrum()
    sig_d = np.sort(np.linalg.eigvalsh(1j * Kd))
    sig_d = sig_d[sig_d > 1e-12][::-1]
    err_sigma = float(np.max(np.abs(sig_lr[: len(sig_d)] - sig_d) / sig_d[0]))

    # 軌道一致（両側とも厳密σ正規化）
    Z_lr = Z.copy()
    Z_d = Z.copy()
    dev = 0.0
    for _ in range(steps):
        sys_lr.set_theta(np.angle(Z_lr))
        sig = sys_lr.sigma_spectrum()[0]
        Z_lr = sys_lr.cayley_step(Z_lr, sig)

        th = np.angle(Z_d)
        Kd = A * np.sin(th[None, :] - th[:, None])
        sd = np.linalg.norm(Kd, 2)
        Kn = Kd / sd
        I = np.eye(m)
        Z_d = np.linalg.solve(I - GAMMA * Kn, (I + GAMMA * Kn) @ Z_d)
        dev = max(dev, float(np.max(np.abs(Z_lr - Z_d))))
''','''    sig_lr = np.sort(np.linalg.eigvalsh(1j * sys_lr.dense_K()))
    sig_lr = sig_lr[sig_lr > 1e-12][::-1]
    sig_d = np.sort(np.linalg.eigvalsh(1j * Kd))
    sig_d = sig_d[sig_d > 1e-12][::-1]
    err_sigma = float(np.max(np.abs(sig_lr[: len(sig_d)] - sig_d) / sig_d[0]))

    # 軌道一致：低ランク側 linear_rotation_step と密行列側 exp(ANGLE·K)（R3(ii)）
    Z_lr = Z.copy()
    Z_d = Z.copy()
    dev = 0.0
    for _ in range(steps):
        sys_lr.set_state(Z_lr)
        Z_lr = sys_lr.linear_rotation_step(Z_lr)

        zz = Z_d / np.abs(Z_d) if KMODE == "phase" else Z_d
        Kd = A * np.imag(np.conj(zz)[:, None] * zz[None, :])
        w, V = np.linalg.eigh(1j * Kd)
        Z_d = V @ (np.exp(-1j * ANGLE * w) * (V.conj().T @ Z_d))
        dev = max(dev, float(np.max(np.abs(Z_lr - Z_d))))
''')
    # onset_probe：A2/A3/S1/S2/R1/A4
    sub(e,"    g = zero_closure_kernel_seed(sys_lr, rng)\n    Z = v + delta * g\n    Z = Z / np.linalg.norm(Z)\n","    Z = v.copy()  # A2/A3/S1: 外部 seed も正規化も無し（zero_closure_kernel_seed は呼ばない）\n")
    sub(e,"    wp = rng.normal(size=sys_lr.m)\n    t0 = time.time()\n    for t in range(cap + 1):\n        h1 = abs(p @ Z) ** 2","    t0 = time.time()  # S2: 冪反復用乱数 wp は使わない\n    for t in range(cap + 1):\n        h1 = abs(p @ Z) ** 2")
    sub(e,"        sys_lr.set_theta(np.angle(Z))\n        sig_est, wp = sys_lr.sigma_max_power(wp)\n        Z = sys_lr.cayley_step(Z, sig_est)\n","        sys_lr.set_state(Z)  # A4\n        Z = sys_lr.linear_rotation_step(Z)  # R1\n")
    # relax_probe：S2/R1/A4（σ は厳密スペクトルで読む＝A6(b)）
    sub(e,"    wp = rng.normal(size=m)\n    plateau_tau = None\n","    plateau_tau = None  # S2: wp は使わない\n")
    sub(e,"        sys_lr.set_theta(np.angle(Z))\n        sig_est, wp = sys_lr.sigma_max_power(wp)\n        if t % sub == 0:\n            sig_exact = sys_lr.sigma_spectrum()\n            sigma_checks.append(abs(sig_est - sig_exact[0]) / sig_exact[0])\n",
          "        sys_lr.set_state(Z)  # A4\n        if t % sub == 0:\n            sig_exact = sys_lr.sigma_spectrum()  # A6(b): 実際の生成子の σ\n            sigma_checks.append(0.0)\n")
    sub(e,"            Z = sys_lr.cayley_step(Z, sig_est)\n    t_run = time.time() - t0\n\n    sys_lr.set_theta(np.angle(Z))\n","            Z = sys_lr.linear_rotation_step(Z)  # R1\n    t_run = time.time() - t0\n\n    sys_lr.set_state(Z)  # A4/A6(b)\n")
    sub(e,'    out = {"n": n, "m": m, "seed": seed, "gamma": GAMMA}','    out = {"n": n, "m": m, "seed": seed, "gamma": GAMMA, "angle": ANGLE, "kmode": KMODE}')
    s=open(e,encoding="utf-8").read(); assert ".cayley_step(" not in s and "set_theta(np.angle(Z" not in s and "sigma_max_power(wp)" not in s, e
# ============================================================ 2. 呼出側スクリプト
CALLERS=[p for p in glob.glob(os.path.join(FX,"**","*.py"),recursive=True) if p not in ENGINES and "NORMALIZED_ORIGINAL" not in p and not p.endswith(("run_followup_experiments.py","run_moduli_sweep_fast.py","analyze_followup.py"))]
for p in CALLERS:
    s=open(p,encoding="utf-8").read(); base=os.path.basename(p)
    for pat,rep in [(r"(\w+)\.set_theta\(np\.angle\((\w+)\)\)", r"\1.set_state(\2)  # A4"), (r"^\s*se, ?wpn? ?= ?\w+\.sigma_max_power\(wp\)\n", ""), (r"(\w+)\.cayley_step\((\w+), ?se\)", r"\1.linear_rotation_step(\2)  # R1"), (r"^\s*wp ?= ?rng\.normal\(size=\w+\)\n", "")]:
        n=len(re.findall(pat,s,re.M))
        if n: sub(p,pat,rep,count=n,regex=True); s=open(p,encoding="utf-8").read()
    if base=="run_N5_physical_phase_step_test.py": sub(p,"        sys.set_theta(theta)\n","        sys.set_state(Z)  # A4\n")
    if base=="run_artifact_comparison_N4_N5.py":
        sub(p,'norm=loadmod(HERE/"run_n_scaling_lowrank_v1_NORMALIZED_ORIGINAL.py","norm")\n','# R2 廃止: K/σ 正規化枝（NORMALIZED_ORIGINAL）は実行しない\n')
        sub(p,"        sex=float(sys.sigma_spectrum()[0])\n        rows.append((t,Ht,Hp,Ho,Ao,f,ztz.real,ztz.imag,abs(ztz),se,sex))\n        if t<STEPS:\n            Z=sys.linear_rotation_step(Z)  # R1; wp=wpn\n",
              "        sex=float(sys.sigma_spectrum()[0])  # A6(b): 実際の生成子の σ₁（sigma_est 列も同値）\n        rows.append((t,Ht,Hp,Ho,Ao,f,ztz.real,ztz.imag,abs(ztz),sex,sex))\n        if t<STEPS:\n            Z=sys.linear_rotation_step(Z)  # R1\n")
        sub(p,"    sys0=norm.LowRankSystem(N)\n    rng=np.random.default_rng(seed)\n    v,residual,sig=norm.make_parent(sys0,rng,iters=1200,tol=1e-12)\n    g=norm.zero_closure_kernel_seed(sys0,rng)\n    Z0=v+DELTA*g\n    Z0=Z0/np.linalg.norm(Z0)\n",
              "    sys0=raw.LowRankSystem(N)\n    rng=np.random.default_rng(seed)\n    v,residual,sig=raw.make_parent(sys0,rng,iters=1200,tol=1e-12)\n    Z0=v.copy()  # A2/A3/S1: 外部 seed も正規化も無し\n")
        sub(p,"    wp0=rng.normal(size=sys0.m)\n","    wp0=None  # S2\n")
        sub(p,"    sys=mod.LowRankSystem(N); Z=Z0.copy(); wp=wp0.copy(); rows=[]\n","    sys=mod.LowRankSystem(N); Z=Z0.copy(); rows=[]  # S2: wp は使わない\n")
        sub(p,'    run_branch(norm,N,Z0,p,q,wp0).to_csv(HERE/f"N{N}_normalized_K_raw_observables.csv",index=False)\n','')
        sub(p,"# Branch difference: only cayley_step K/sigma normalization.","# R2 abolished: raw branch only (linear rotation exp(ANGLE K), amplitude-aware K).")
    if base=="run_N5_gamma_continuum_test.py":
        sub(p,"    gamma = math.tan(math.pi / n_den)\n    engine.GAMMA = gamma\n    dphi = 2.0 * math.atan(gamma)\n","    gamma = math.tan(math.pi / n_den)  # 記録用\n    engine.ANGLE = 2.0 * math.pi / n_den  # R1/R3(iii): 線形回転の刻み角を掃引\n    dphi = engine.ANGLE\n")
    if base=="run_complex_simplex_decompactification.py":
        # 修正後の力学では指数成長域（1e-10<R_perp<1e-3 が 10 点以上）が無い場合があり fit が None になる。書式化だけ安全にする（物理は不変）
        sub(p,"            f\"- R_perp early log growth rate: {s['R_perp_log_growth_rate_per_step']:.6f}/step\",","            (f\"- R_perp early log growth rate: {s['R_perp_log_growth_rate_per_step']:.6f}/step\" if s['R_perp_log_growth_rate_per_step'] is not None else \"- R_perp early log growth rate: not fitted (no exponential window 1e-10<R_perp<1e-3 with >=10 points)\"),")
        sub(p,"        'N':n,'M':m,'steps':STEPS,'seed':SEED,'gamma':float(eng.GAMMA),\n        'physics':'raw K Cayley; K/sigma normalization removed',","        'N':n,'M':m,'steps':STEPS,'seed':SEED,'gamma':float(eng.GAMMA),'angle':float(eng.ANGLE),'kmode':eng.KMODE,\n        'physics':'linear rotation exp(ANGLE K); K = Im(conj z_i z_j) (KMODE=amplitude) or phase-only (KMODE=phase); no normalization, no seed',")
    s=open(p,encoding="utf-8").read(); assert ".cayley_step(" not in s and "set_theta(np.angle(Z" not in s and "sigma_max_power" not in s and "rng.normal(size" not in s, p
# ============================================================ 3. followup（自前エンジン）
f=os.path.join(FX,"N5_dynamics_followup_theorems_and_stability_20260826","followup_dynamics_20260826","run_followup_experiments.py")
sub(f,"import math, json, itertools, time\n","import math, json, itertools, time, os\n")
sub(f,"GAMMA=math.tan(math.pi/144.0)\n","GAMMA=math.tan(math.pi/144.0)  # 旧 Cayley 刻み（記録用）\nANGLE=2.0*math.pi/144.0        # R1\nKMODE=os.environ.get('KMODE','amplitude')  # B\n")
sub(f,'''    def set_theta(self,theta):
        n=self.n; self.c=np.cos(theta); self.s=np.sin(theta)
        T=np.zeros((n,n)); T[self.ea,self.eb]=theta; T[self.eb,self.ea]=theta
        CT=np.cos(T); ST=np.sin(T); np.fill_diagonal(CT,0); np.fill_diagonal(ST,0)
''','''    def set_theta(self,theta): self.set_state(np.exp(1j*np.asarray(theta,dtype=float)),force_phase=True)  # 親の固有モード反復でのみ使用（S3）
    def set_state(self,z,force_phase=False):
        # A4: 振幅込み K_ij=Im(conj(z_i)z_j)。KMODE='phase' なら z/|z|（対照 baseline）
        n=self.n; z=np.asarray(z,dtype=complex)
        if force_phase or KMODE=='phase':
            r=np.abs(z); z=np.where(r>0,z/np.where(r>0,r,1.0),0.0)
        self.c=np.real(z).copy(); self.s=np.imag(z).copy()
        CT=np.zeros((n,n)); ST=np.zeros((n,n)); CT[self.ea,self.eb]=self.c; CT[self.eb,self.ea]=self.c; ST[self.ea,self.eb]=self.s; ST[self.eb,self.ea]=self.s
''')
sub(f,"    def cayley_step(self,z,sigma_unused=None):\n","    def dense_K(self): return np.column_stack([self.kmatvec(e) for e in np.eye(self.m)])\n    def linear_rotation_step(self,z,sigma_unused=None):\n        # R1: z ← exp(ANGLE·K) z\n        K=self.dense_K(); w,V=np.linalg.eigh(1j*K); return V@(np.exp(-1j*ANGLE*w)*(V.conj().T@z))\n    def cayley_step(self,z,sigma_unused=None):  # 旧更新（実行しない）\n")
sub(f,"            v=sys.w(EV[:,idx].astype(complex)); v=v/np.linalg.norm(v); theta_new=np.angle(v)\n","            v=sys.w(EV[:,idx].astype(complex)); theta_new=np.angle(v)  # A1\n")
sub(f,"                sys.set_theta(np.angle(Z)); Z=sys.cayley_step(Z)\n","                sys.set_state(Z); Z=sys.linear_rotation_step(Z)  # A4/R1\n")
sub(f,"    sys.set_theta(np.angle(v)); Fv=sys.cayley_step(v); phase=np.vdot(v,Fv); phase=phase/abs(phase)\n","    sys.set_state(v); Fv=sys.linear_rotation_step(v); phase=np.vdot(v,Fv); phase=phase/abs(phase)  # A4/R1\n")
sub(f,"        sys.set_theta(np.angle(z)); return sys.cayley_step(z)/phase\n","        sys.set_state(z); return sys.linear_rotation_step(z)/phase  # A4/R1\n")
sub(f,"        for t in range(5000): sys.set_theta(np.angle(Z)); Z=sys.cayley_step(Z)\n","        for t in range(5000): sys.set_state(Z); Z=sys.linear_rotation_step(Z)  # A4/R1\n")
m=os.path.join(FX,"N5_dynamics_followup_theorems_and_stability_20260826","followup_dynamics_20260826","run_moduli_sweep_fast.py")
sub(m,"        sysm.set_theta(np.angle(Z)); Z=sysm.cayley_step(Z)\n","        sysm.set_state(Z); Z=sysm.linear_rotation_step(Z)  # A4/R1\n")
a=os.path.join(FX,"N5_dynamics_followup_theorems_and_stability_20260826","followup_dynamics_20260826","analyze_followup.py")
sub(a,"GAMMA=math.tan(math.pi/144.0)\n","GAMMA=math.tan(math.pi/144.0)\nANGLE=2.0*math.pi/144.0  # R1\n")
sub(a,"norm=pd.read_csv(ROOT/'N5_sigma_normalization_artifact_test/N5_normalized_K_raw_observables.csv')\nfor d,kind in [(norm,'normalized'),(raw,'raw')]:\n    if kind=='normalized': inc=np.full(len(d),2*math.atan(GAMMA))\n    else: inc=2*np.arctan(GAMMA*d.sigma_exact.to_numpy())\n",
      "norm=raw.copy()  # R2 廃止: 正規化枝は無い。raw の累積位相のみ（比較は自明に一致）\nfor d,kind in [(norm,'normalized'),(raw,'raw')]:\n    inc=ANGLE*d.sigma_exact.to_numpy()  # R3(i): 線形回転の位相増分 ANGLE·σ₁（σ₁ は実際の生成子）\n")
for p in (f,m,a):
    s=open(p,encoding="utf-8").read(); assert ".cayley_step(" not in s and "set_theta(np.angle(Z" not in s and "set_theta(np.angle(z" not in s, p
open(os.path.join(HERE,"results","apply_fixes.log"),"w",encoding="utf-8").write("\n".join(log)+"\n"); print(f"applied {len(log)} substitutions in {len(set(l.split(' | ')[0] for l in log))} files")
