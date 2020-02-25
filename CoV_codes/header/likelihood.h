double my_log(double x)
{
    double a=(1.0+x)/2.0;
    double b=sqrt(x);
    while(1){
        a=(a+b)/2.0;
        b=sqrt(a*b);
        if(fabs(a-b)<0.001) break;
    }
    return 2.0*(x-1.0)/(a+b);
}

double Priors(Par *p)
{
    double fv0=gsl_ran_lognormal_pdf(p->v0,log(2),1);
    double fdv=gsl_ran_lognormal_pdf(p->dv,log(2),1);
    double fdEs=gsl_ran_lognormal_pdf(p->dEs,log(p->dv/p->piv*p->v0),1);
    double fE0s=gsl_ran_lognormal_pdf(p->E0s,log(2),1);
    double fpiv=gsl_ran_lognormal_pdf(p->piv,log(2),1);
    double fR0=gsl_ran_lognormal_pdf(E0*p->beta*p->piv/p->dv/p->dEs,log(2),1);
    return my_log(fv0)+my_log(fdEs)+my_log(fdv)+my_log(fE0s)+my_log(fR0)+my_log(fpiv);
}

double likelihood(Par *p)
{
    int i;
    double llhood=0;
    SolveIntegral(p);
 // m=SolveIntegral(p,2);
    for(i=0;i<dl;i++)
    {
        int n=data1[i].day/dt+1;
        llhood+=data1[i].D*log(s[n].Eps)-s[n].Eps;
        //double val=data1[i].D-s[n].v;
        //llhood-=val*val;
    }
    llhood-=SumLogGamma;
    llhood+=Priors(p);
   // free(s);
    return llhood;
}
