#!/usr/bin/env python3
import importlib.util
import math
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('m', HERE/'run_paper4_charged_extension_v1.py')
m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m)


def average_fixed_flux(p,e,alpha,nu,Gamma,beta,n=40000):
    chi=np.linspace(0.0,2*math.pi,n,endpoint=False)
    pg=[];pe=[];w=[]
    for z in chi:
        fl=m.radiation_fluxes(z,p,e,z,alpha,1.0,nu,Gamma,beta,True,True)
        pg.append(fl['gw_power_total']); pe.append(fl['em_power_total'])
        w.append(fl['r']**2/math.sqrt(alpha*p))
    pg=np.asarray(pg); pe=np.asarray(pe); w=np.asarray(w)
    return float(np.mean(pg*w)/np.mean(w)), float(np.mean(pe*w)/np.mean(w))


def analytic_flux(p,e,alpha,nu,Gamma,beta):
    a=p/(1-e*e)
    fg=(1+73*e*e/24+37*e**4/96)/(1-e*e)**3.5
    gw=(32/5)*nu**2*alpha**3*a**-5*fg
    em=(1/3)*Gamma*beta**2*alpha**2*(2+e*e)/(a**4*(1-e*e)**2.5)
    return gw,em

cases=[
    ('gravity',24.3,math.sqrt(1-(27/30)**2),1.0,0.25,0.0,1.0),
    ('charged_rho100_scaled',2454.3,math.sqrt(1-(27/30)**2),101.0,0.25,25.0,1.0),
    ('charged_rho100_d10',24543.0,math.sqrt(1-(27/30)**2),101.0,0.25,25.0,1.0),
]
print('case,gw_numeric,gw_analytic,gw_relerr,em_numeric,em_analytic,em_relerr')
for name,p,e,alpha,nu,Gamma,beta in cases:
    gn,en=average_fixed_flux(p,e,alpha,nu,Gamma,beta)
    ga,ea=analytic_flux(p,e,alpha,nu,Gamma,beta)
    print(name,gn,ga,abs(gn/ga-1),en,ea,(0 if ea==0 else abs(en/ea-1)),sep=',')
