#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <time.h>

static long double parse_initial(const char *s){
    if(strcmp(s,"invsqrt2")==0) return sqrtl(0.5L);
    return strtold(s,NULL);
}

int main(int argc, char **argv){
    if(argc < 5){
        fprintf(stderr,"usage: %s INITIAL MAX_GENERATIONS CSV_PATH CONTINUE_EXACT_ZERO(0/1)\n", argv[0]);
        return 2;
    }
    long double z1 = parse_initial(argv[1]);
    int N = atoi(argv[2]);
    const char *csv_path = argv[3];
    int continue_exact_zero = atoi(argv[4]);
    if(N < 1) return 3;

    long double *z = calloc((size_t)N, sizeof(long double));
    if(!z) return 4;
    FILE *f = fopen(csv_path,"w");
    if(!f){perror("fopen"); free(z); return 5;}
    fprintf(f,"generation,new_value,new_abs,log10_abs,sign,ordered_pair_terms\n");

    z[0] = z1;
    fprintf(f,"1,%.21Le,%.21Le,%.21Le,%d,0\n", z[0], fabsl(z[0]), log10l(fabsl(z[0])), (z[0]>0)-(z[0]<0));

    int last_n = 1;
    int underflow_generation = 0;
    clock_t t0=clock();
    for(int n=1; n<N; ++n){
        long double L = 0.0L;
        /* Literal ordered-pair calculation; no factorisation of the double sum. */
        for(int i=0; i<n; ++i){
            long double zi = z[i]; /* conjugate is identical on this real-axis test */
            for(int j=0; j<n; ++j){
                long double zj = z[j];
                long double term = (zi * zj) * zj;
                L += term;
            }
        }
        z[n] = -L;
        long double a = fabsl(z[n]);
        long double lg = (a>0.0L) ? log10l(a) : -INFINITY;
        fprintf(f,"%d,%.21Le,%.21Le,%.21Le,%d,%lld\n", n+1, z[n], a, lg, (z[n]>0)-(z[n]<0), (long long)n*(long long)n);
        last_n = n+1;

        if(z[n] == 0.0L && z[n-1] != 0.0L && !continue_exact_zero){
            underflow_generation = n+1;
            break;
        }
    }
    fclose(f);
    double sec=(double)(clock()-t0)/CLOCKS_PER_SEC;
    long double last=z[last_n-1];
    printf("initial=%.21Le\n",z1);
    printf("requested_generations=%d\n",N);
    printf("completed_generations=%d\n",last_n);
    printf("last=%.21Le\n",last);
    if(last!=0.0L) printf("last_log10_abs=%.21Le\n",log10l(fabsl(last)));
    else printf("last_log10_abs=-inf\n");
    printf("underflow_generation=%d\n",underflow_generation);
    printf("LDBL_MIN=%.21Le\n",(long double)LDBL_MIN);
    printf("elapsed_cpu_sec=%.6f\n",sec);
    free(z);
    return 0;
}
