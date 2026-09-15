#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from pathlib import Path
from itertools import combinations
import numpy as np
import pandas as pd
import csv
import matplotlib.pyplot as plt
import matplotlib.animation as animation

HERE=Path(__file__).resolve().parent
RAW=HERE.parent/"K4_reproducibility_raw_20260915"
OUT=HERE/"reproduced_outputs"
OUT.mkdir(parents=True,exist_ok=True)

FPS=12; HOLD=12; CYCLES=3
node_pos=np.array([[0,1.0],[1.15,0],[0,-1.0],[-1.15,0]])
edges=list(combinations(range(4),2))

def make_video(m,direction):
    direction_name="forward" if direction==1 else "reverse"
    rows=[]
    with open(RAW/"complex_states_animation_exact_i.csv", newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            if int(r["m"])==m and r["direction"]==direction_name:
                rows.append({
                    "m":int(r["m"]), "direction":r["direction"], "step":int(r["step"]), "k":int(r["k"]),
                    "theta_rad":float(r["theta_rad"]),
                    "base_re":float(r["base_re"]), "base_im":float(r["base_im"]),
                    "harmonic_re":float(r["harmonic_re"]), "harmonic_im":float(r["harmonic_im"])
                })
    q=pd.DataFrame(rows).sort_values("step").reset_index(drop=True)
    wd=pd.read_csv(RAW/f"wave_samples_mp4_1200_m{m}.csv", float_precision="round_trip")
    theta=wd["theta_rad"].to_numpy()
    y1=wd["base"].to_numpy()
    ym=wd["harmonic"].to_numpy()
    ys=wd["composite"].to_numpy()

    fig=plt.figure(figsize=(12,8))
    gs=fig.add_gridspec(2,2,height_ratios=[1,1.05],width_ratios=[1.15,.85])
    axN=fig.add_subplot(gs[0,0]); axC=fig.add_subplot(gs[0,1]); axW=fig.add_subplot(gs[1,:])

    for i,j in edges:
        temporal=((j-i)%4==1) or ((i-j)%4==1)
        axN.plot([node_pos[i,0],node_pos[j,0]],[node_pos[i,1],node_pos[j,1]],
                 linestyle="-" if temporal else "--",
                 linewidth=2.4 if temporal else 1.5)
    axN.scatter(node_pos[:,0],node_pos[:,1],s=110,facecolors="white",edgecolors="black",zorder=3)
    for k,(x,y) in enumerate(node_pos): axN.text(x+.08,y+.06,f"S{k}",fontsize=11)
    currentN,=axN.plot([],[],'o',markersize=20,markerfacecolor='none',
                       markeredgecolor='red',markeredgewidth=4,zorder=5)
    axN.set_xlim(-1.5,1.5); axN.set_ylim(-1.3,1.3); axN.set_aspect("equal"); axN.axis("off")
    # forward original used ASCII hyphen, reverse generator used em dash
    axN.set_title(f"Complete network K4 — harmonic m={m}")

    # unit circle is display scaffolding, not experimental state data
    t=np.linspace(0,2*np.pi,400)
    axC.plot(np.cos(t),np.sin(t),'--',linewidth=1)
    axC.axhline(0,linewidth=.8); axC.axvline(0,linewidth=.8)
    baseVec,=axC.plot([],[],linewidth=3,label="base")
    harmVec,=axC.plot([],[],linewidth=3,label=f"harmonic m={m}")
    baseDot,=axC.plot([],[],'ro',markersize=7)
    harmDot,=axC.plot([],[],'ro',markersize=7)
    axC.set_xlim(-1.25,1.25); axC.set_ylim(-1.25,1.25); axC.set_aspect("equal")
    axC.set_xlabel("Re"); axC.set_ylabel("Im"); axC.set_title("Complex plane")
    axC.legend(loc="upper right",fontsize=8)

    off1,offm,offs=3.0,0.0,-3.5
    axW.plot(theta,y1+off1,linewidth=2,label="base: cos(theta)")
    axW.plot(theta,ym+offm,linewidth=2,label=f"harmonic: cos({m} theta)")
    axW.plot(theta,ys+offs,linewidth=2,label="sum")
    for y in [off1,offm,offs]: axW.axhline(y,linewidth=.6)
    phaseLine=axW.axvline(0,linestyle="--",linewidth=2)
    d1,=axW.plot([],[],'ro',markersize=7); dm,=axW.plot([],[],'ro',markersize=7); ds,=axW.plot([],[],'ro',markersize=7)
    axW.set_xlim(0,2*np.pi); axW.set_ylim(-6,4.5)
    axW.set_xticks([0,np.pi/2,np.pi,3*np.pi/2,2*np.pi],["0","pi/2","pi","3pi/2","2pi"])
    axW.set_ylabel("amplitude (vertically separated)")
    axW.set_title("Base wave, harmonic, and composite waveform")
    axW.legend(loc="upper right",fontsize=8,ncol=3)

    state_text=fig.text(.5,.965,"",ha="center",fontsize=14)
    if direction==1:
        fig.suptitle("Discrete phase evolution: U^4 = I",y=.995,fontsize=16)
    else:
        fig.suptitle("Discrete phase evolution: U^4 = I (reverse)",y=.995,fontsize=16)
    fig.tight_layout(rect=[0,0,.99,.95])

    def update(frame):
        r=q.iloc[(frame//HOLD)%4]
        k=int(r["k"]); th=float(r["theta_rad"])
        z1=complex(float(r["base_re"]),float(r["base_im"]))
        zm=complex(float(r["harmonic_re"]),float(r["harmonic_im"]))
        currentN.set_data([node_pos[k,0]],[node_pos[k,1]])
        baseVec.set_data([0,z1.real],[0,z1.imag]); harmVec.set_data([0,zm.real],[0,zm.imag])
        baseDot.set_data([z1.real],[z1.imag]); harmDot.set_data([zm.real],[zm.imag])
        phaseLine.set_xdata([th,th])
        d1.set_data([th],[float(r["base_re"])+off1])
        dm.set_data([th],[float(r["harmonic_re"])+offm])
        # The discrete composite marker is the real part of z_base + z_harmonic.
        ds.set_data([th],[float(r["base_re"])+float(r["harmonic_re"])+offs])
        if direction==1:
            state_text.set_text(f"k = {k}    theta = {k} pi/2")
        else:
            state_text.set_text(f"k = {k}    theta = {k} pi/2    step = -1 mod 4")
        return currentN,baseVec,harmVec,baseDot,harmDot,phaseLine,d1,dm,ds,state_text

    ani=animation.FuncAnimation(fig,update,frames=4*HOLD*CYCLES,interval=1000/FPS,blit=False)
    writer=animation.FFMpegWriter(fps=FPS,codec="libx264",bitrate=2200,
                                  extra_args=["-pix_fmt","yuv420p","-movflags","+faststart"])
    suffix="" if direction==1 else "_reverse"
    ani.save(OUT/f"K4_harmonic_discrete_phase_m{m}{suffix}_from_raw_20260915.mp4",
             writer=writer,dpi=130)
    plt.close(fig)

if __name__ == "__main__":
    import argparse
    ap=argparse.ArgumentParser()
    ap.add_argument("--m",type=int,choices=[2,3])
    ap.add_argument("--direction",choices=["forward","reverse"])
    args=ap.parse_args()
    if args.m is not None and args.direction is not None:
        make_video(args.m,+1 if args.direction=="forward" else -1)
    else:
        for m in (2,3):
            make_video(m,+1)
            make_video(m,-1)
    print("MP4 readout complete")
