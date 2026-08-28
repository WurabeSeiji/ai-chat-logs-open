#include <bits/stdc++.h>
using namespace std;
struct E{double re,im,mag; int i,j;};
struct R{double r; vector<int> s;};
int main(int argc,char**argv){
 if(argc<4){cerr<<"N file out\n"; return 1;} int N=atoi(argv[1]); string file=argv[2], out=argv[3];
 ifstream f(file); string line; getline(f,line); vector<string> hdr; {stringstream ss(line); string x; while(getline(ss,x,',')) hdr.push_back(x);} int ci=-1,cj=-1,cr=-1,cm=-1; for(int a=0;a<(int)hdr.size();a++){if(hdr[a]=="i")ci=a;if(hdr[a]=="j")cj=a;if(hdr[a]=="z2_re")cr=a;if(hdr[a]=="z2_im")cm=a;}
 vector<E> e; while(getline(f,line)){stringstream ss(line); vector<string> v; string x; while(getline(ss,x,','))v.push_back(x); if(v.size()!=hdr.size())continue; E z; z.i=stoi(v[ci]);z.j=stoi(v[cj]);z.re=stod(v[cr]);z.im=stod(v[cm]);z.mag=hypot(z.re,z.im);e.push_back(z);} int M=e.size();
 mt19937_64 rng(20260826+N); ofstream fo(out); fo<<"N,M,k,best_residual,edges\n";
 for(int k=2;k<=N-2;k++){
   int restarts = (k<=6?5000:3000); double best=1e100; vector<int> bestS;
   vector<int> all(M); iota(all.begin(),all.end(),0);
   for(int rr=0;rr<restarts;rr++){
     shuffle(all.begin(),all.end(),rng); vector<char> in(M,0); vector<int> S; S.reserve(k); double sr=0,si=0,den=0; for(int a=0;a<k;a++){int x=all[a];in[x]=1;S.push_back(x);sr+=e[x].re;si+=e[x].im;den+=e[x].mag;}
     auto resid=[&](){return hypot(sr,si)/den;};
     double cur=resid();
     for(int iter=0;iter<100;iter++){
       double best2=cur; int bo=-1,bi=-1; double nsr=0,nsi=0,nden=0;
       for(int ai=0;ai<k;ai++){int o=S[ai]; for(int x=0;x<M;x++) if(!in[x]){ double tr=sr-e[o].re+e[x].re, ti=si-e[o].im+e[x].im, td=den-e[o].mag+e[x].mag; double r=hypot(tr,ti)/td; if(r<best2-1e-15){best2=r;bo=ai;bi=x;nsr=tr;nsi=ti;nden=td;}}}
       if(bo<0) break; int old=S[bo];in[old]=0;in[bi]=1;S[bo]=bi;sr=nsr;si=nsi;den=nden;cur=best2;
     }
     if(cur<best){best=cur;bestS=S; if(best<1e-10) break;}
     // perturb around local minimum with a few random double swaps every 20 restarts not implemented
   }
   sort(bestS.begin(),bestS.end()); fo<<N<<","<<M<<","<<k<<","<<setprecision(17)<<best<<",\""; for(size_t a=0;a<bestS.size();a++){if(a)fo<<";";fo<<e[bestS[a]].i<<"-"<<e[bestS[a]].j;}fo<<"\"\n"; cerr<<"N="<<N<<" k="<<k<<" best="<<setprecision(9)<<best<<"\n";
 }
}
