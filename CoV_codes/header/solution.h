typedef struct Solution
{
    double Ep,Eps,v;
}Solution;
Solution *s;

double BackEuler(double A, double B, double xt)
{
    return (A*dt+xt)/(1+B*dt);
}

void  PrintSolution(FILE *fs)
{
    int i;
    int EnT=(data1[dl-1].day+1)/dt;
    int dii=1.0/(4*dt);
    for(i=0;i<EnT;i+=dii) fprintf(fs,"%lf\t", s[i].Eps);
    fprintf(fs,"\n");
}

void SolveIntegral(Par *p)
{
    int i=0;
    int EnT=(data1[dl-1].day+1)/dt;
  //  double R0=p->beta*p->piv*E0/p->dv/p->dEs;
    s=malloc(sizeof(Solution)*(EnT+1));
    s[0].Ep=E0;
    s[0].Eps=p->E0s;
    s[0].v=p->v0;
    while (i<=EnT){
        s[i+1].Ep=BackEuler(dE*E0,p->beta*s[i].v+dE,s[i].Ep);
        s[i+1].Eps=BackEuler(p->beta*s[i+1].Ep*s[i].v,p->dEs,s[i].Eps);
        s[i+1].v=BackEuler(p->piv*s[i+1].Eps,p->dv,s[i].v);
        i++;
    }
   // return s;
}
