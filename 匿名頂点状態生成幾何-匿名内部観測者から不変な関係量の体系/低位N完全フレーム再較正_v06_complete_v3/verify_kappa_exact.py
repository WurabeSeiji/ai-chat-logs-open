#!/usr/bin/env python3
import math
import pandas as pd
from pathlib import Path

BASE = Path(__file__).resolve().parent
obs = pd.read_csv(BASE/'lowN_v06_summary.csv')
rows=[]
for _,r in obs.iterrows():
    N=int(r.N); m=N-1
    if m==1:
        mean=1.0; median=1.0
    else:
        mean=(math.sqrt(math.pi)/2)*math.gamma((m+1)/2)/math.gamma((m+2)/2)
        median=math.sqrt(1-2**(-2/(m-1)))
    rows.append({
        'N':N,'m':m,
        'observed_mean_kappa':r.mean_kappa,
        'exact_mean_kappa':mean,
        'mean_error':r.mean_kappa-mean,
        'observed_median_kappa':r.median_kappa,
        'exact_median_kappa':median,
        'median_error':r.median_kappa-median,
    })
out=pd.DataFrame(rows)
out.to_csv(BASE/'kappa_exact_comparison.csv',index=False)
print(out.to_string(index=False))
