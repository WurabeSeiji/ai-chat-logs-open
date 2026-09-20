import csv, subprocess
from decimal import Decimal, getcontext
from pathlib import Path
getcontext().prec=80
root=Path('/mnt/data/singularity_bottom_refine_20260919')
bin=str(root/'literal_bottom_sweep_quad')
center=Decimal('0.8080759710485481360243945033000000482251')
with open(root/'quad_deep_140_best.csv','w') as outsum:
    outsum.write('stage,step_exp,r_best,final_log10_abs,last40_mean_log10_ratio,last20_mean_log10_ratio\n')
    for stage,exp in enumerate(range(29,35),start=28):
        step=Decimal(1).scaleb(-exp)
        start=center-step*5; end=center+step*5
        out=root/f'stage{stage}_quad140_1e-{exp}.csv'
        subprocess.check_call([bin,format(start,'f'),format(end,'f'),format(step,'f'),'140',str(out)])
        rows=list(csv.DictReader(open(out)))
        best=min(rows,key=lambda r:Decimal(r['final_log10_abs']))
        center=Decimal(best['r'])
        outsum.write(f"{stage},{exp},{best['r']},{best['final_log10_abs']},{best['last40_mean_log10_ratio']},{best['last20_mean_log10_ratio']}\n")
        outsum.flush()
        print(stage,exp,best['r'],best['final_log10_abs'],best['last40_mean_log10_ratio'],flush=True)
