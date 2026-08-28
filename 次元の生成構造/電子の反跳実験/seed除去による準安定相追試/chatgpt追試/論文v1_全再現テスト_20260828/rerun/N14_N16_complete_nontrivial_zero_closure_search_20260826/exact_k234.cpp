#include <bits/stdc++.h>
using namespace std; struct E{double re,im,mag;int i,j;};
int main(int ac,char**av){int N=atoi(av[1]); string file=av[2],out=av[3]; ifstream f(file);string line;getline(f,line);vector<string>h;{stringstream ss(line);string x;while(getline(ss,x,','))h.push_back(x);}int ci=-1,cj=-1,cr=-1,cm=-1;for(int a=0;a<h.size();a++){if(h[a]=="i")ci=a;if(h[a]=="j")cj=a;if(h[a]=="z2_re")cr=a;if(h[a]=="z2_im")cm=a;}vector<E>e;while(getline(f,line)){stringstream ss(line);vector<string>v;string x;while(getline(ss,x,','))v.push_back(x);if(v.size()!=h.size())continue;E z{stod(v[cr]),stod(v[cm]),0,stoi(v[ci]),stoi(v[cj])};z.mag=hypot(z.re,z.im);e.push_back(z);}int M=e.size();ofstream fo(out);fo<<"N,M,k,best_residual,edges\n";
for(int k=2;k<=4;k++){double best=1e99;vector<int>b;
 if(k==2){for(int a=0;a<M;a++)for(int c=a+1;c<M;c++){double r=hypot(e[a].re+e[c].re,e[a].im+e[c].im)/(e[a].mag+e[c].mag);if(r<best){best=r;b={a,c};}}}
 if(k==3){for(int a=0;a<M;a++)for(int c=a+1;c<M;c++)for(int d=c+1;d<M;d++){double r=hypot(e[a].re+e[c].re+e[d].re,e[a].im+e[c].im+e[d].im)/(e[a].mag+e[c].mag+e[d].mag);if(r<best){best=r;b={a,c,d};}}}
 if(k==4){for(int a=0;a<M;a++)for(int c=a+1;c<M;c++)for(int d=c+1;d<M;d++)for(int g=d+1;g<M;g++){double r=hypot(e[a].re+e[c].re+e[d].re+e[g].re,e[a].im+e[c].im+e[d].im+e[g].im)/(e[a].mag+e[c].mag+e[d].mag+e[g].mag);if(r<best){best=r;b={a,c,d,g};}}}
 fo<<N<<","<<M<<","<<k<<","<<setprecision(17)<<best<<",\"";for(int x:b)fo<<e[x].i<<"-"<<e[x].j<<";";fo<<"\"\n"; cerr<<N<<" k="<<k<<" "<<best<<"\n";}
}
