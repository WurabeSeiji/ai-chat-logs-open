#include <stdio.h>
#include <stdlib.h>
#include <quadmath.h>
typedef __float128 Q; typedef struct{Q re,im;} C;
static inline C add(C a,C b){C z={a.re+b.re,a.im+b.im};return z;}
static inline C mul(C a,C b){C z={a.re*b.re-a.im*b.im,a.re*b.im+a.im*b.re};return z;}
static inline C cj(C a){C z={a.re,-a.im};return z;}
static inline Q ab(C a){return hypotq(a.re,a.im);}
static void qg(char*b,size_t n,Q x){quadmath_snprintf(b,n,"%.40Qg",x);} static void qe(char*b,size_t n,Q x){quadmath_snprintf(b,n,"%.36Qe",x);}
static void eval(Q r,int G,Q *fl,Q *m40){
 C*z=calloc((size_t)G+2,sizeof(C)); Q*a=calloc((size_t)G+2,sizeof(Q)); z[1].re=r;a[1]=ab(z[1]);
 for(int N=1;N<G;++N){C s={0,0};for(int i=1;i<=N;++i){C ci=cj(z[i]);for(int j=1;j<=N;++j){C t=mul(mul(ci,z[j]),z[j]);s=add(s,t);}}z[N+1].re=-s.re;z[N+1].im=-s.im;a[N+1]=ab(z[N+1]);}
 *fl=a[G]>0?log10q(a[G]):-HUGE_VALQ; Q s=0;int c=0;for(int n=G-39;n<=G;++n){if(n>=2&&a[n]>0&&a[n-1]>0){s+=log10q(a[n]/a[n-1]);c++;}}*m40=s/c;free(z);free(a);
}
int main(int argc,char**argv){if(argc!=5){fprintf(stderr,"usage center half_ulps G out.csv\n");return 2;}Q c=strtoflt128(argv[1],NULL);int h=atoi(argv[2]),G=atoi(argv[3]);FILE*f=fopen(argv[4],"w");fprintf(f,"ulp_offset,r_decimal,r_hex,final_log10_abs,last40_mean_log10_ratio\n");Q r=c;for(int k=0;k<h;k++)r=nextafterq(r,-HUGE_VALQ);for(int o=-h;o<=h;++o){Q fl,m;eval(r,G,&fl,&m);char br[128],bh[128],bf[128],bm[128];qg(br,sizeof br,r);quadmath_snprintf(bh,sizeof bh,"%Qa",r);qe(bf,sizeof bf,fl);qe(bm,sizeof bm,m);fprintf(f,"%d,%s,%s,%s,%s\n",o,br,bh,bf,bm);r=nextafterq(r,HUGE_VALQ);}fclose(f);return 0;}
