#define MIN(x,y)  (((x)<(y)) ?  (x) : (y))
#define MAX(x,y)  (((x)>(y)) ?  (x) : (y))
double dE=1e-3;
double E0=25;
double dt=0.002;

typedef struct Par
{
    double dEs,dv;
    double piv;
    double beta;
    double E0s,v0;
}Par;

typedef struct  DataSet
{
    double day;
    double D;
   // int r;
   // double *D;
}DataSet;
DataSet *data1;
int dl;
double SumLogGamma=0;

void init()
{
	FILE *fp;
    fp=fopen("data/d19_s.csv","r");
    if (fp==NULL) {
        printf("No that files!\n");
        exit(1);
    }
    int len=32,i=0;
   // int j;
    data1=malloc(sizeof(DataSet)*len);
    while(!feof(fp)){
        fscanf(fp,"%lf,%lf\n",&data1[i].day,&data1[i].D);
       // data1[i].day+=3;
        SumLogGamma+=gsl_sf_lngamma(data1[i].D+1);
        i++;
        if (i>=len) {
            len+=32;
            data1=realloc(data1,sizeof(DataSet)*len);
        }
    }
    dl=i;
    fclose(fp);
}

void InitPar(Par *p)
{
    p->dEs=1.2;
    p->piv=1e-2;
    p->dv=1;
    p->beta=0.1; 
    p->E0s=10;
    p->v0=2.5;
}

