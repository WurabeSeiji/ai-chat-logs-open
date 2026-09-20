#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>
#include <string.h>

typedef struct { long double re, im; } C;
static inline C c_add(C a, C b){ C z={a.re+b.re,a.im+b.im}; return z; }
static inline C c_mul(C a, C b){ C z={a.re*b.re-a.im*b.im,a.re*b.im+a.im*b.re}; return z; }
static inline C c_conj(C a){ C z={a.re,-a.im}; return z; }
static inline long double c_abs(C a){ return hypotl(a.re,a.im); }

static long double parse_ld(const char *s){
    errno=0; char *e=NULL; long double v=strtold(s,&e);
    if(errno || e==s || *e!='\0'){ fprintf(stderr,"invalid long double: %s\n",s); exit(2); }
    return v;
}
static int parse_i(const char *s){ char *e=NULL; long v=strtol(s,&e,10); if(e==s || *e!='\0' || v<2){fprintf(stderr,"invalid int: %s\n",s); exit(2);} return (int)v; }

int main(int argc, char **argv){
    if(argc != 6){
        fprintf(stderr,"usage: %s r_start r_end step max_gen out_csv\n", argv[0]);
        return 2;
    }
    const long double r_start=parse_ld(argv[1]);
    const long double r_end=parse_ld(argv[2]);
    const long double step=parse_ld(argv[3]);
    const int max_gen=parse_i(argv[4]);
    const char *out_csv=argv[5];
    int num_r=(int)llroundl((r_end-r_start)/step)+1;
    if(num_r<1){ fprintf(stderr,"empty range\n"); return 2; }

    FILE *f=fopen(out_csv,"w");
    if(!f){ perror("fopen"); return 1; }
    fprintf(f,"r,generations,final_abs,final_log10_abs,min_abs,min_gen,max_abs,max_gen,sign_flips,first_flip_gen,last50_mean_log10_ratio,last50_std_log10_ratio,last20_mean_log10_ratio,last20_std_log10_ratio\n");

    for(int idx=0; idx<num_r; ++idx){
        long double r=r_start+step*(long double)idx;
        C *z=(C*)calloc((size_t)max_gen+2,sizeof(C));
        long double *absv=(long double*)calloc((size_t)max_gen+2,sizeof(long double));
        if(!z||!absv){ fprintf(stderr,"allocation failed\n"); return 3; }
        z[1].re=r; z[1].im=0.0L; absv[1]=c_abs(z[1]);
        long double min_abs=absv[1], max_abs=absv[1];
        int min_gen=1,max_gen_seen=1,sign_flips=0,first_flip=-1;
        int prev_sign=(z[1].re>0)?1:((z[1].re<0)?-1:0);

        for(int N=1; N<max_gen; ++N){
            C sum={0.0L,0.0L};
            // Literal rule: every ordered pair, self interaction included.
            for(int i=1;i<=N;++i){
                C ci=c_conj(z[i]);
                for(int j=1;j<=N;++j){
                    C t=c_mul(c_mul(ci,z[j]),z[j]);
                    sum=c_add(sum,t);
                }
            }
            z[N+1].re=-sum.re; z[N+1].im=-sum.im;
            absv[N+1]=c_abs(z[N+1]);
            if(absv[N+1]<min_abs){min_abs=absv[N+1];min_gen=N+1;}
            if(absv[N+1]>max_abs){max_abs=absv[N+1];max_gen_seen=N+1;}
            int s=(z[N+1].re>0)?1:((z[N+1].re<0)?-1:0);
            if(prev_sign!=0 && s!=0 && s!=prev_sign){ sign_flips++; if(first_flip<0) first_flip=N+1; }
            if(s!=0) prev_sign=s;
        }

        long double sums[2]={0,0}, sums2[2]={0,0};
        int counts[2]={0,0};
        int windows[2]={50,20};
        for(int w=0;w<2;++w){
            int start=max_gen-windows[w]+1; if(start<2) start=2;
            for(int n=start;n<=max_gen;++n){
                if(absv[n]>0.0L && absv[n-1]>0.0L){
                    long double lr=log10l(absv[n]/absv[n-1]);
                    sums[w]+=lr; sums2[w]+=lr*lr; counts[w]++;
                }
            }
        }
        long double mean50=NAN,std50=NAN,mean20=NAN,std20=NAN;
        if(counts[0]){ mean50=sums[0]/counts[0]; long double v=sums2[0]/counts[0]-mean50*mean50; if(v<0)v=0; std50=sqrtl(v); }
        if(counts[1]){ mean20=sums[1]/counts[1]; long double v=sums2[1]/counts[1]-mean20*mean20; if(v<0)v=0; std20=sqrtl(v); }
        long double final_abs=absv[max_gen];
        long double final_log=(final_abs>0)?log10l(final_abs):-INFINITY;
        fprintf(f,"%.18Lf,%d,%.24Le,%.24Le,%.24Le,%d,%.24Le,%d,%d,%d,%.24Le,%.24Le,%.24Le,%.24Le\n",
            r,max_gen,final_abs,final_log,min_abs,min_gen,max_abs,max_gen_seen,sign_flips,first_flip,mean50,std50,mean20,std20);
        free(z); free(absv);
    }
    fclose(f);
    return 0;
}
