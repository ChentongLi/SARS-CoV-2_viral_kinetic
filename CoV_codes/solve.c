#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sf.h>
#include "header/init.h"
#include "header/solution.h"
#include "header/likelihood.h"
#include "header/MCMC.h"
double TT=200;

void solve(Par *p)
{
    int i=0;
    int EnT=TT/dt;
    s=malloc(sizeof(Solution)*(EnT+1));
    s[0].Ep=E0;
    s[0].Eps=p->E0s;
    s[0].v=p->v0;
    while (i<=EnT){
        s[i+1].Ep=BackEuler(dE*E0,p->beta*s[i].v+dE,s[i].Ep);
        s[i+1].Eps=BackEuler(p->beta*s[i+1].Ep,p->dEs,s[i].Eps);
        s[i+1].v=BackEuler(p->piv*s[i+1].Eps,p->dv,s[i].v);
        i++;
    }
    
}
void rwrite(FILE *fp)
{
    int i;
    int EnT=TT/dt; 
    int dii=1.0/(4*dt);
    for(i=0;i<EnT;i+=dii) fprintf(fp,"%lf,%lf,%lf,%lf\n",i*dt,s[i].Ep,s[i].Eps,s[i].v);
    free(s);
    fclose(fp); 
}

int main() 
{   
    Par *para;
    FILE *fp;
    para=malloc(sizeof(Par));
    para->dEs=0.0942;
    para->piv=0.251;
    para->dv=1.418;
    para->beta=0.0682;
    para->E0s=1.515;
    para->v0=0.096;
    //begining
    solve(para);
    fp=fopen("result/solution_original.csv","w+");
    rwrite(fp); 
    return 1;
}
