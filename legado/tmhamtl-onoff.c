#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <string.h>
#include <omp.h>

double espaco=500;

double nlinhas(FILE *file)
{
   int ch, prev = '\n';
   double lines = 0;
   while ( (ch = fgetc(file)) != EOF )
   {
      if ( ch == '\n' )
      {
         ++lines;
      }
      prev = ch;
   }
   rewind(file);
   if ( prev != '\n' )
   {
      ++lines;
   }
   return(lines);
}

double argfasetmof(double alpha,void * params)
{
	double temp1,temp2,temp3,temp4;
	double *pp= (double *) params;
	double d = pp[0];
	double t = pp[1];
	double g = pp[2];
	double thtm=pp[3];
	double m=pp[4];
	double w0=pp[5];
	double of=pp[6];
	temp1=2*sqrt(d*(t-of)/M_PI)*exp(-pow(alpha,2)*(d*(t-of)+pow(w0,2)/8));
	temp2=2*sqrt(d*t/M_PI)*exp(-pow(alpha,2)*(d*t+pow(w0,2)/8));
	temp3=exp(-pow(alpha*w0,2)/8)/alpha;
	temp4=((2*d*t*pow(alpha,2)-1)*gsl_sf_erfc(sqrt(d*t)*alpha)-(2*d*(t-of)*pow(alpha,2)-1)*gsl_sf_erfc(sqrt(d*(t-of))*alpha));
	return((temp1-temp2+temp3*temp4)*2*thtm*(gsl_sf_bessel_J0(sqrt(m*g)*w0*alpha)-1));
}

double argfasetm(double alpha,void * params)
{
	double temp1,temp2,temp3,temp4;
	double *pp= (double *) params;
	double d = pp[0];
	double t = pp[1];
	double g = pp[2];
	double thtm=pp[3];
	double m=pp[4];
	double w0=pp[5];
	temp1=thtm*2*exp(-pow(alpha*w0,2)/8);
	temp2=-2*sqrt(d*t/M_PI)*exp(-d*t*pow(alpha,2));
	temp3=gsl_sf_erf(sqrt(d*t)*alpha)/alpha;
	temp4=gsl_sf_erfc(sqrt(d*t)*alpha)*2*d*t*alpha;
	return(temp1*(temp2+temp3+temp4)*(gsl_sf_bessel_J0(sqrt(m*g)*w0*alpha)-1));
}

double fasetm(double t, double d, double g, double thtm,double m,double w0, double of)
{
  gsl_integration_workspace *wgsl=gsl_integration_workspace_alloc (espaco);
  double result, error;
  double p[7];
  p[0]=d;
  p[1]=t;
  p[2]=g;
  p[3]=thtm;
  p[4]=m;
  p[5]=w0;
  p[6]=of;
  gsl_function F;
  if(t<=of)
    {F.function=&argfasetm;}
    else
    {F.function=&argfasetmof;}
  F.params=p;
  gsl_integration_qag(&F,0,2e5,0,1e-4,espaco,1,wgsl,&result,&error);
  gsl_integration_workspace_free (wgsl);
  return (result);
}

double argtau(double tau,void * params)
{
	double temp1,temp2,temp3;
	double *pp2= (double *) params;
	double d=pp2[0];
	double t=pp2[1];
	double g=pp2[2];
	double m=pp2[3];
	double w0=pp2[4];
	double df=pp2[5];
	double alpha=pp2[7];
	temp1=(sqrt(df/(M_PI*tau))/d)*(exp(-pow(alpha,2)*df*tau-pow(alpha*w0,2)/8));
	temp2=gsl_sf_erf(sqrt(d*(t-tau))*alpha);
	temp3=(gsl_sf_bessel_J0(sqrt(m*g)*w0*alpha)-1);
	return(temp1*temp2*temp3);
}

double argfasear(double alpha, void * params)
{
  gsl_integration_workspace *w=gsl_integration_workspace_alloc (espaco);
  double result, error, inferior;
  double ptau[8];
  gsl_function F;
  double *pp= (double *) params;
  ptau[0]=pp[0];
  ptau[1]=pp[1];
  ptau[2]=pp[2];
  ptau[3]=pp[3];
  ptau[4]=pp[4];
  ptau[5]=pp[5];
  ptau[6]=pp[6];
  if(ptau[1]<=ptau[6])
    {inferior=0;}
    else
    {inferior=ptau[1]-ptau[6];}
  ptau[7]=alpha;
  F.function=&argtau;
  F.params=ptau;
  gsl_integration_qag(&F,inferior,ptau[1],0,1e-4,espaco,1,w,&result,&error);
  gsl_integration_workspace_free (w);
  return (result);
}

double fasear(double t,  double d, double g, double thtm,double m,double w0,double beta, double nu,double dndt,double df, double of)
{
  gsl_integration_workspace *w=gsl_integration_workspace_alloc (espaco);
  double result, error,comsemar;
  double par[7];
  par[0]=d;
  par[1]=t;
  par[2]=g;
  par[3]=m;
  par[4]=w0;
  par[5]=df;
  par[6]=of;
  if(dndt==0)
	{comsemar=0;}
  else
	{gsl_function F;
  	F.function=&argfasear;
  	F.params=par;
  	gsl_integration_qag(&F,0,2e5,0,1e-4,espaco,1,w,&result,&error);
  	gsl_integration_workspace_free (w);
	comsemar=-thtm*(2*d*dndt)/(beta*(1 + nu))*result;}
  return (comsemar);
}

double j(double g, void * params)
{
	double *pp= (double *) params;
	double t = pp[0];
	double d = pp[1];
	double thtm = pp[2];
	double v=pp[3];
	double m=pp[4];
	double w0=pp[5];
	double beta=pp[6];
	double nu=pp[7];
	double dndt=pp[8];
	double df=pp[9];
	double of=pp[10];
	return(exp(-g)*cos(v*g+fasetm(t,d,g,thtm,m,w0,of)+fasear(t,d,g,thtm,m,w0,beta,nu,dndt,df,of)));
}

double amplireal (double t, double d, double thtm,double v,double m,double w0, double beta, double nu,double dndt, double df, double of)
{
  gsl_integration_workspace *w=gsl_integration_workspace_alloc (espaco);
  double result, error;
  double p[11];
  p[0]=t;
  p[1]=d;
  p[2]=thtm;
  p[3]=v;
  p[4]=m;
  p[5]=w0;
  p[6]=beta;
  p[7]=nu;
  p[8]=dndt;
  p[9]=df;
  p[10]=of;
  gsl_function F;
  F.function=&j;
  F.params=p;
  gsl_integration_qag(&F,0,30,0,1e-4,espaco,1,w,&result,&error);
  gsl_integration_workspace_free (w);
  return (result);
}

double k(double g, void * params)
{
	double *pp= (double *) params;
	double t=pp[0];
	double d=pp[1];
	double thtm=pp[2];
	double v=pp[3];
	double m=pp[4];
	double w0=pp[5];
	double beta=pp[6];
	double nu=pp[7];
	double dndt=pp[8];
	double df=pp[9];
    double of=pp[10];
	return(exp(-g)*sin(v*g+fasetm(t,d,g,thtm,m,w0,of)+fasear(t,d,g,thtm,m,w0,beta,nu,dndt,df,of)));
}

/*Integra k, obtendo a parte complexa da amplitude*/
double amplicomplex (double t, double d, double thtm,double v,double m,double w0, double beta, double nu,double dndt,double df, double of)
{
  gsl_integration_workspace *w=gsl_integration_workspace_alloc (espaco);
  double result, error;
  double p[11];
  p[0]=t;
  p[1]=d;
  p[2]=thtm;
  p[3]=v;
  p[4]=m;
  p[5]=w0;
  p[6]=beta;
  p[7]=nu;
  p[8]=dndt;
  p[9]=df;
  p[10]=of;
  gsl_function F;
  F.function=&k;
  F.params=p;
  gsl_integration_qag(&F,0,30,0,1e-4,espaco,1,w,&result,&error);
  gsl_integration_workspace_free (w);
  return (result);
}

/*intensity corresponde ao módulo ao quadrado da amplitude complexa*/
double intensity (double t, double d, double thtm,double v,double m,double w0, double beta, double nu,double dndt, double df, double of)
{
	double temp1;
	if (t<=0.000000001)
	{
	temp1=1;
	}
	else
	{
	temp1=(pow(amplireal(t,d,thtm,v,m,w0,beta,nu,dndt,df,of),2)+pow(amplicomplex (t,d,thtm,v,m,w0,beta,nu,dndt,df,of),2))*(1+v*v);
	}
	return(temp1);
}

int main (int argc, char *argv[])
{
   double th,d,chi,m,w0,v,beta,nu,dndt,df,of;
   int i,n;
   char path[100];
   FILE *arq,*arqresult;

   d=atof(argv[2]);
   th=atof(argv[3]);
   m=atof(argv[4]);
   w0=atof(argv[5]);
   v=atof(argv[6]);
   beta=atof(argv[7]);
   nu=atof(argv[8]);
   dndt=atof(argv[9]);
   df=atof(argv[10]);
   of=atof(argv[11]);
   arq=fopen(argv[1],"r");
   sprintf(path,"result-%s-%s-%s",argv[2],argv[3],argv[1]);
   arqresult=fopen(path,"w");
   n=nlinhas(arq);
   gsl_matrix * data=gsl_matrix_alloc(n,2);
   gsl_matrix * ma=gsl_matrix_alloc(n,2);
   gsl_matrix_fscanf(arq,data);
   fclose(arq);
         #pragma omp parallel for default(shared) private(i) reduction(+:chi)
         for(i=0;i<n;i++)
         {
	  gsl_matrix_set(ma,i,0,gsl_matrix_get(data,i,0));
	  gsl_matrix_set(ma,i,1,intensity(gsl_matrix_get(ma,i,0),d,th,v,m,w0,beta,nu,dndt,df,of));
	  chi+=pow(gsl_matrix_get(data,i,1) - gsl_matrix_get(ma,i,1),2);
         }
         #pragma omp barrier
   printf("%.10lf\n",chi);
   for(i=0;i<n;i++)
   {
      fprintf(arqresult,"%lf\t%lf\n",gsl_matrix_get(ma,i,0),gsl_matrix_get(ma,i,1));
   }
   fclose(arqresult);
   gsl_matrix_free(data);
   gsl_matrix_free(ma);
   return (0);
}
