#!/usr/bin/env python3
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from itertools import combinations

OUT = Path(__file__).resolve().parent
FPS = 12
HOLD = 12
CYCLES = 3
node_pos = np.array([[0,1.0],[1.15,0],[0,-1.0],[-1.15,0]])
edges = list(combinations(range(4),2))

def state_complex(k,m):
    U=1j
    return U**k, U**(m*k)

def make_video(m,path):
    fig=plt.figure(figsize=(12,8))
    gs=fig.add_gridspec(2,2,height_ratios=[1,1.05],width_ratios=[1.15,.85])
    axN=fig.add_subplot(gs[0,0]); axC=fig.add_subplot(gs[0,1]); axW=fig.add_subplot(gs[1,:])
    for i,j in edges:
        temporal=((j-i)%4==1) or ((i-j)%4==1)
        axN.plot([node_pos[i,0],node_pos[j,0]],[node_pos[i,1],node_pos[j,1]],linestyle='-' if temporal else '--',linewidth=2.4 if temporal else 1.5)
    axN.scatter(node_pos[:,0],node_pos[:,1],s=110,facecolors='white',edgecolors='black',zorder=3)
    for k,(x,y) in enumerate(node_pos): axN.text(x+.08,y+.06,f'S{k}',fontsize=11)
    currentN,=axN.plot([],[],'o',markersize=20,markerfacecolor='none',markeredgecolor='red',markeredgewidth=4,zorder=5)
    axN.set_xlim(-1.5,1.5); axN.set_ylim(-1.3,1.3); axN.set_aspect('equal'); axN.axis('off'); axN.set_title(f'Complete network K4 - harmonic m={m}')
    t=np.linspace(0,2*np.pi,400); axC.plot(np.cos(t),np.sin(t),'--',linewidth=1); axC.axhline(0,linewidth=.8); axC.axvline(0,linewidth=.8)
    baseVec,=axC.plot([],[],linewidth=3,label='base'); harmVec,=axC.plot([],[],linewidth=3,label=f'harmonic m={m}'); baseDot,=axC.plot([],[],'ro',markersize=7); harmDot,=axC.plot([],[],'ro',markersize=7)
    axC.set_xlim(-1.25,1.25); axC.set_ylim(-1.25,1.25); axC.set_aspect('equal'); axC.set_xlabel('Re'); axC.set_ylabel('Im'); axC.set_title('Complex plane'); axC.legend(loc='upper right',fontsize=8)
    theta=np.linspace(0,2*np.pi,1200); y1=np.cos(theta); ym=np.cos(m*theta); ys=y1+ym; off1,offm,offs=3.0,0.0,-3.5
    axW.plot(theta,y1+off1,linewidth=2,label='base: cos(theta)'); axW.plot(theta,ym+offm,linewidth=2,label=f'harmonic: cos({m} theta)'); axW.plot(theta,ys+offs,linewidth=2,label='sum')
    for y in [off1,offm,offs]: axW.axhline(y,linewidth=.6)
    phaseLine=axW.axvline(0,linestyle='--',linewidth=2); d1,=axW.plot([],[],'ro',markersize=7); dm,=axW.plot([],[],'ro',markersize=7); ds,=axW.plot([],[],'ro',markersize=7)
    axW.set_xlim(0,2*np.pi); axW.set_ylim(-6,4.5); axW.set_xticks([0,np.pi/2,np.pi,3*np.pi/2,2*np.pi],['0','pi/2','pi','3pi/2','2pi']); axW.set_ylabel('amplitude (vertically separated)'); axW.set_title('Base wave, harmonic, and composite waveform'); axW.legend(loc='upper right',fontsize=8,ncol=3)
    state_text=fig.text(.5,.965,'',ha='center',fontsize=14); fig.suptitle('Discrete phase evolution: U^4 = I',y=.995,fontsize=16); fig.tight_layout(rect=[0,0,.99,.95])
    def update(frame):
        k=(frame//HOLD)%4; th=k*np.pi/2; z1,zm=state_complex(k,m)
        currentN.set_data([node_pos[k,0]],[node_pos[k,1]]); baseVec.set_data([0,z1.real],[0,z1.imag]); harmVec.set_data([0,zm.real],[0,zm.imag]); baseDot.set_data([z1.real],[z1.imag]); harmDot.set_data([zm.real],[zm.imag]); phaseLine.set_xdata([th,th]); d1.set_data([th],[np.cos(th)+off1]); dm.set_data([th],[np.cos(m*th)+offm]); ds.set_data([th],[np.cos(th)+np.cos(m*th)+offs]); state_text.set_text(f'k = {k}    theta = {k} pi/2')
        return currentN,baseVec,harmVec,baseDot,harmDot,phaseLine,d1,dm,ds,state_text
    ani=animation.FuncAnimation(fig,update,frames=4*HOLD*CYCLES,interval=1000/FPS,blit=False)
    writer=animation.FFMpegWriter(fps=FPS,codec='libx264',bitrate=2200,extra_args=['-pix_fmt','yuv420p','-movflags','+faststart'])
    ani.save(path,writer=writer,dpi=130); plt.close(fig)

if __name__ == '__main__':
    for m in (2,3):
        make_video(m,OUT/f'K4_harmonic_discrete_phase_m{m}_20260915.mp4')
