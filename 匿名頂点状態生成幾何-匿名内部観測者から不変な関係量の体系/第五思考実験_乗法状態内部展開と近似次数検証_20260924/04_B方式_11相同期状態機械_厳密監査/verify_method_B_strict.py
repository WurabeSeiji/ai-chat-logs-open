#!/usr/bin/env python3
import math, json, ast, inspect, importlib.util
import numpy as np

spec=importlib.util.spec_from_file_location('bfast','/mnt/data/audit_two_paths_rk4_fast.py')
bfast=importlib.util.module_from_spec(spec); spec.loader.exec_module(bfast)
H,TAU0=bfast.H,bfast.TAU0
ref_step,init_B,micro_fast,macro_fast,coreB=bfast.ref_step,bfast.init_B,bfast.micro_fast,bfast.macro_fast,bfast.coreB

CASES=[]
for p in (24.3,60.0,120.0):
    for e in (0.20,0.45,0.55): CASES.append((p,e,1.0))
for q in (2.0,4.0,10.0): CASES.append((60.0,0.45,q))

# Static check: RHS in micro_fast must not read the newly constructed n dict.
src=inspect.getsource(micro_fast)
tree=ast.parse(src)
viol=[]
class V(ast.NodeVisitor):
    def visit_Subscript(self,node):
        # Any Load of n[...] is same-step read-after-write risk.
        if isinstance(node.ctx,ast.Load) and isinstance(node.value,ast.Name) and node.value.id=='n':
            viol.append((node.lineno,ast.unparse(node)))
        self.generic_visit(node)
V().visit(tree)

# Random one-step equivalence over physically relevant ranges.
rng=np.random.default_rng(20260924)
rand_max=np.zeros(5)
phase_fail=0
external_hold_fail=0
qcycle_fail=0
reset_fail=0
for _ in range(5000):
    P=float(rng.uniform(20,140)); E=float(rng.uniform(0.05,0.65)); qmass=float(rng.choice([1.,2.,4.,10.]))
    nu=qmass/(1+qmass)**2
    chi=float(rng.uniform(-math.pi,math.pi)); phi=float(rng.uniform(-10,10)); t=float(rng.uniform(0,2e4))
    U=complex(math.cos(chi),math.sin(chi)); Hp=complex(math.cos(phi),math.sin(phi)); Q=math.exp(t/TAU0)
    x=np.array([U,P,E,Hp,Q],complex)
    s=init_B(*x)
    core0=coreB(s).copy(); q0=s['q'].copy()
    # Walk 11 microsteps and check phase/hold invariants.
    for j in range(11):
        s2=micro_fast(s,nu)
        # q must remain one-hot and rotate exactly one slot.
        qexp=np.roll(s['q'],1)
        if not (np.array_equal(s2['q'],qexp) and np.count_nonzero(s2['q']==1.0)==1 and np.sum(s2['q'])==1.0):
            phase_fail+=1
        # External core must be held for phases 0..9.
        if j<10 and not np.array_equal(coreB(s2), coreB(s)):
            external_hold_fail+=1
        s=s2
    if not np.array_equal(s['q'],q0): qcycle_fail+=1
    # Work registers reset at commit.
    if not (np.array_equal(s['k1'],np.zeros(3)) and np.array_equal(s['k2'],np.zeros(3)) and np.array_equal(s['k3'],np.zeros(3)) and np.array_equal(s['k4'],np.zeros(3)) and np.array_equal(s['d'],np.zeros(3)) and s['dphi']==0.0):
        reset_fail+=1
    r=ref_step(x,nu); b=coreB(s)
    rand_max=np.maximum(rand_max,np.abs(b-r))

# Long-run all 12 cases.
case_rows=[]; global_max=np.zeros(5)
for p0,e0,qmass in CASES:
    nu=qmass/(1+qmass)**2
    x=np.array([1+0j,p0,e0,1+0j,1+0j],complex)
    s=init_B(*x)
    m=np.zeros(5)
    for _ in range(12000):
        x=ref_step(x,nu)
        s=macro_fast(s,nu)
        m=np.maximum(m,np.abs(coreB(s)-x))
    global_max=np.maximum(global_max,m)
    case_rows.append({'p0':p0,'e0':e0,'q':qmass,'max_abs':m.tolist()})

# Compare against saved C1 raw output.
data=np.genfromtxt('/mnt/data/experiment02_sn_2p5pn_rr_C1_N12.csv',delimiter=',',names=True)
p0=float(data['p'][0]); e0=float(data['e'][0]); nu=0.25
s=init_B(1+0j,p0,e0,1+0j,1+0j)
csv_err={'p':0.0,'e':0.0,'t':0.0,'x':0.0,'y':0.0}
for i,row in enumerate(data):
    U,P,E,Hp,Q=coreB(s)
    ce=bfast.base.trigU(U)[0]
    rad=P.real/(1+E.real*ce)
    xx=rad*Hp.real; yy=rad*Hp.imag; tt=TAU0*math.log(Q.real)
    csv_err['p']=max(csv_err['p'],abs(P.real-row['p']))
    csv_err['e']=max(csv_err['e'],abs(E.real-row['e']))
    csv_err['t']=max(csv_err['t'],abs(tt-row['t_osc']))
    csv_err['x']=max(csv_err['x'],abs(xx-row['x']))
    csv_err['y']=max(csv_err['y'],abs(yy-row['y']))
    if i+1<len(data): s=macro_fast(s,nu)

out={
 'static_new_state_reads':viol,
 'random_5000_one_macro_max_abs':rand_max.tolist(),
 'phase_fail':phase_fail,
 'external_hold_fail':external_hold_fail,
 'qcycle_fail':qcycle_fail,
 'reset_fail':reset_fail,
 'all12_global_max_abs':global_max.tolist(),
 'cases':case_rows,
 'saved_C1_csv_max_abs':csv_err,
}
with open('/mnt/data/method_B_strict_verification.json','w') as f: json.dump(out,f,indent=2)
print(json.dumps(out,indent=2))
