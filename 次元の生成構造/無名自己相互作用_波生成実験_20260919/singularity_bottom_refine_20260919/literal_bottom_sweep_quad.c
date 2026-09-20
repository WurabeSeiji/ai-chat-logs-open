#include <stdio.h>
#include <stdlib.h>
#include <quadmath.h>
#include <math.h>
#include <errno.h>

typedef __float128 Q;
typedef struct { Q re, im; } Cq;
static inline Cq add(Cq a,Cq b){Cq z={a.re+b.re,a.im+b.im};return z;}
static inline Cq mul(Cq a,Cq b){Cq z={a.re*b.re-a.im*b.im,a.re*b.im+a.im*b.re};return z;}
static inline Cq conjq2(Cq a){Cq z={a.re,-a.im};return z;}
static inline Q absq2(Cq a){return hypotq(a.re,a.im);}
static void qstr(char *buf,size_t n,Q x){quadmath_snprintf(buf,n,"%.40Qg",x);}
static void qestr(char *buf,size_t n,Q x){quadmath_snprintf(buf,n,"%.36Qe",x);}
int main(int argc,char**argv){
 if(argc!=6){fprintf(stderr,"usage: %s start end step max_gen out.csv\n",argv[0]);return 2;}
 Q start=strtoflt128(argv[1],NULL), end=strtoflt128(argv[2],NULL), step=strtoflt128(argv[3],NULL);
 int max_gen=atoi(argv[4]); const char *out=argv[5];
 Q span=(end-start)/step; long long num=(long long)llroundq(span)+1;
 FILE*f=fopen(out,"w"); if(!f){perror("fopen");return 1;}
 fprintf(f,"r,generations,final_abs,final_log10_abs,last40_mean_log10_ratio,last20_mean_log10_ratio\n");
 for(long long idx=0;idx<num;++idx){
   Q r=start+step*(Q)idx;
   Cq *z=calloc((size_t)max_gen+2,sizeof(Cq)); Q *a=calloc((size_t)max_gen+2,sizeof(Q));
   z[1].re=r; a[1]=absq2(z[1]);
   for(int N=1;N<max_gen;++N){
     Cq sum={0,0};
     for(int i=1;i<=N;++i){ Cq ci=conjq2(z[i]);
       for(int j=1;j<=N;++j){ Cq t=mul(mul(ci,z[j]),z[j]); sum=add(sum,t); }
     }
     z[N+1].re=-sum.re; z[N+1].im=-sum.im; a[N+1]=absq2(z[N+1]);
   }
   Q s40=0,s20=0; int c40=0,c20=0;
   for(int n=max_gen-39;n<=max_gen;++n){if(n>=2&&a[n]>0&&a[n-1]>0){s40+=log10q(a[n]/a[n-1]);c40++;}}
   for(int n=max_gen-19;n<=max_gen;++n){if(n>=2&&a[n]>0&&a[n-1]>0){s20+=log10q(a[n]/a[n-1]);c20++;}}
   Q flog=a[max_gen]>0?log10q(a[max_gen]):-HUGE_VALQ;
   char br[128],ba[128],bl[128],bm40[128],bm20[128];
   qstr(br,sizeof(br),r); qestr(ba,sizeof(ba),a[max_gen]); qestr(bl,sizeof(bl),flog); qestr(bm40,sizeof(bm40),s40/c40); qestr(bm20,sizeof(bm20),s20/c20);
   fprintf(f,"%s,%d,%s,%s,%s,%s\n",br,max_gen,ba,bl,bm40,bm20);
   free(z);free(a);
 }
 fclose(f); return 0;
}
