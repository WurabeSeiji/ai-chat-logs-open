#include <bits/stdc++.h>
using namespace std;struct E{double re,im,mag;int i,j;};
int main(int ac,char**av){int N=atoi(av[1]),k=atoi(av[2]);string file=av[3];long long steps=atoll(av[4]);uint64_t seed=stoull(av[5]);ifstream f(file);string line;getline(f,line);vector<string>h;{stringstream ss(line);string x;while(getline(ss,x,','))h.push_back(x);}int ci=-1,cj=-1,cr=-1,cm=-1;for(int a=0;a<h.size();a++){if(h[a]=="i")ci=a;if(h[a]=="j")cj=a;if(h[a]=="z2_re")cr=a;if(h[a]=="z2_im")cm=a;}vector<E>e;while(getline(f,line)){stringstream ss(line);vector<string>v;string x;while(getline(ss,x,','))v.push_back(x);if(v.size()!=h.size())continue;E z{stod(v[cr]),stod(v[cm]),0,stoi(v[ci]),stoi(v[cj])};z.mag=hypot(z.re,z.im);e.push_back(z);}int M=e.size();mt19937_64 rng(seed);uniform_real_distribution<double> U(0,1);vector<int>all(M);iota(all.begin(),all.end(),0);double global=1e99;vector<int>gS;
int chains=20; long long per=max(1LL,steps/chains);
for(int ch=0;ch<chains;ch++){shuffle(all.begin(),all.end(),rng);vector<char>in(M);vector<int>S;double sr=0,si=0,den=0;for(int a=0;a<k;a++){int x=all[a];in[x]=1;S.push_back(x);sr+=e[x].re;si+=e[x].im;den+=e[x].mag;}auto R=[&](){return hypot(sr,si)/den;};double cur=R(),best=cur;vector<int>bS=S;double T=0.01;
 for(long long t=0;t<per;t++){int ai=rng()%k;int x=rng()%M; if(in[x])continue;int o=S[ai];double nr=sr-e[o].re+e[x].re,ni=si-e[o].im+e[x].im,nd=den-e[o].mag+e[x].mag;double rr=hypot(nr,ni)/nd;double temp=0.02*pow(1e-5,(double)(t%200000)/200000.0); if(rr<cur || U(rng)<exp((cur-rr)/max(temp,1e-12))){in[o]=0;in[x]=1;S[ai]=x;sr=nr;si=ni;den=nd;cur=rr;if(cur<best){best=cur;bS=S;}if(cur<global){global=cur;gS=S;}}
 if((t+1)%200000==0){S=bS;fill(in.begin(),in.end(),0);sr=si=den=0;for(int y:S){in[y]=1;sr+=e[y].re;si+=e[y].im;den+=e[y].mag;}cur=R();}
 }
}
sort(gS.begin(),gS.end());cout<<setprecision(17)<<global<<"\t";for(int x:gS)cout<<e[x].i<<"-"<<e[x].j<<";";cout<<"\n";}
