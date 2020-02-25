#define MAXMCTIMES 1e7

void SetValue (Par *out, Par *in)
{
    out->dEs=in->dEs;
    out->piv=in->piv;
    out->dv=in->dv;
    out->beta=in->beta;
    out->E0s=in->E0s;
    out->v0=in->v0;
}

void RandPar(Par *p, gsl_rng *r)
{
    // gsl_ran_lognormal(r,log(mu),10.0/mu);
    // gsl_ran_beta(r,x*n,(1-x)*n);
    double u;
    u=gsl_rng_uniform(r)-0.5;
    p->dEs*=exp(u*0.5);
    u=gsl_rng_uniform(r)-0.5;
    p->piv*=exp(u);
    u=gsl_rng_uniform(r)-0.5;
    p->dv*=exp(u*0.5);
    u=gsl_rng_uniform(r)-0.5;
    p->beta*=exp(u); 
    u=gsl_rng_uniform(r)-0.5;
    p->E0s*=exp(u);
    u=gsl_rng_uniform(r)-0.5;
    p->v0*=exp(u*0.5);

}
void printPara(FILE *fp,Par *p)
{
    fprintf(fp,"%e\t%e\t%e\t%e\t%lf\t%lf\t%lf\n",p->piv,p->beta,p->dEs,p->dv,
    p->E0s,p->v0,E0*p->beta*p->piv/p->dv/p->dEs);
}
void MCMC()
{
    const gsl_rng_type *T;
    gsl_rng *r;   
    T=gsl_rng_ranlxs0; 
    gsl_rng_default_seed = ((unsigned long)(time(NULL))); 
    r=gsl_rng_alloc(T); 

    int ncount=0;
    int id=getpid();
    double ALPHA,llhood,R;
    double pllhood=-1e50;
    char fname[30],fname2[30];

    Par *para,*tp;
    para=malloc(sizeof(Par));
    tp=malloc(sizeof(Par));
    InitPar(para);

    sprintf(fname,"result/%d_par.csv",id); 
    sprintf(fname2,"result/%d_solution.csv",id);
    FILE *fp,*fs;
    fp=fopen(fname,"w+");
    fs=fopen(fname2,"w+");
    printf("The initilization is over, then go to the calculation process, %d!\n",id);
    
    while (ncount<MAXMCTIMES){

        SetValue(tp,para);
        RandPar(para,r);
        llhood=likelihood(para);

        ALPHA=MIN(0,llhood-pllhood);
        R=log(gsl_rng_uniform(r));

        //  printf("%e, %e,%lf\n",llhood,pllhood,R);
        if (R<ALPHA){
            pllhood=llhood;
            printf("%lf\n",llhood);
            printPara(fp,para);
            PrintSolution(fs);
        }
        else{
            SetValue(para,tp);
        }
        free(s);
        ncount++;
    }
    fclose(fp);
    fclose(fs);
    gsl_rng_free(r); 
}

