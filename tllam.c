#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sf_dawson.h>
#include <string.h>
#include <omp.h>

double ns=1.63,nf=1,E=72,coefexp=7.5,dndts=5.3,nu=0.23,qpen=0.9e-3,qpar=0.09e-3,espaco=500;

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

double argfasetl(double alpha, void * params)
{
  double qci2,hh,resto;
  double *pp= (double *) params;
  double d=pp[0];
  double t=pp[1];
  double g=pp[2];
  double m=pp[3];
  double w0=pp[4];
  double l=pp[5];
  double df=pp[6];
  qci2=pow(ns,3)*E*coefexp/(4*(1-nu));
  hh=(cosh(l*alpha)-1)/(l*alpha+sinh(l*alpha));
  resto=l*exp(-pow(w0*alpha,2)/8)*(1-exp(-d*t*pow(alpha,2)))*(gsl_sf_bessel_J0(alpha*w0*sqrt(m*g))-1)/alpha;
  return (resto*(dndts+4*(ns-nf)*(1+nu)*coefexp*hh/(l*alpha)+qci2*(qpar+3*qpen-4*(qpar*nu+qpen*(2+nu))*hh/(l*alpha))));
}

double fasetl(double t, double d,double df, double g, double thtl,double m,double w0, double l)
{
  gsl_integration_workspace *w=gsl_integration_workspace_alloc (espaco);
  double result, error;
  double p[7];
  p[0]=d;
  p[1]=t;
  p[2]=g;
  p[3]=m;
  p[4]=w0;
  p[5]=l;
  p[6]=df;
  gsl_function F;
  F.function=&argfasetl;
  F.params=p;
  gsl_integration_qag(&F,0,2e5,0,1e-5,espaco,1,w,&result,&error);
  gsl_integration_workspace_free (w);
  return (thtl*result);
}

double argfasear(double alpha, void * params)
{
  double temp0,temp1,temp2,bb1,bb2,bb3,bbexp;
  double *pp= (double *) params;
  double d=pp[0];
  double t=pp[1];
  double g=pp[2];
  double m=pp[3];
  double w0=pp[4];
  double df=pp[5];
  double cs=pp[6];
  double cf=pp[7];
  double dens=pp[8];
  double denf=pp[9];
  double k=d*cs*dens;
  double kf=df*cf*denf;
  temp0=(pow(k,2)-pow(kf,2))*d*df/(pow(k,2)*df-pow(kf,2)*d);
  temp1=exp(-pow(w0*alpha,2)/8)/((pow(k,2)-pow(kf,2))*pow(alpha,2));
  temp2=gsl_sf_bessel_J0(sqrt(m*g)*w0*alpha)-1;
  bb1=k*gsl_sf_erf(alpha*sqrt(df*t))-kf*gsl_sf_erf(alpha*sqrt(d*t));
  if(temp0>d)
  {
  bb3=(2*kf/sqrt(M_PI))*sqrt(d/(temp0-d))*exp(-pow(alpha,2)*t*d)*gsl_sf_dawson(alpha*sqrt((temp0-d)*t));
  }
  else
  {
  bbexp=-exp(-temp0*pow(alpha,2)*t);
  bb3=bbexp*kf*sqrt(d/(d-temp0))*gsl_sf_erf(alpha*sqrt((d-temp0)*t));
  }
  if(temp0>df)
  {
  bb2=(2*k/sqrt(M_PI))*sqrt(df/(temp0-df))*exp(-pow(alpha,2)*t*df)*gsl_sf_dawson(alpha*sqrt((temp0-df)*t));
  }
  else
  {
  bbexp=-exp(-temp0*pow(alpha,2)*t);
  bb2=bbexp*k*sqrt(df/(df-temp0))*gsl_sf_erf(alpha*sqrt((df-temp0)*t));
  }
  return (temp1*(bb1+bb2-bb3)*temp2);
}

double fasear(double t, double d,double df, double g, double th,double m,double w0,double cs, double cf,double dens, double denf)
{
  gsl_integration_workspace *w=gsl_integration_workspace_alloc (espaco);
  double result, error;
  double par[10];
  par[0]=d;
  par[1]=t;
  par[2]=g;
  par[3]=m;
  par[4]=w0;
  par[5]=df;
  par[6]=cs;
  par[7]=cf;
  par[8]=dens;
  par[9]=denf;
  gsl_function F;
  F.function=&argfasear;
  F.params=par;
  gsl_integration_qag(&F,0,2e5,0,1e-5,espaco,1,w,&result,&error);
  gsl_integration_workspace_free (w);
  return (th*result);
}

double j(double g, void * params)
{
	double *pp= (double *) params;
	double t = pp[0];
	double d = pp[1];
	double thsample = pp[2];
	double v=pp[3];
	double m=pp[4];
	double w0=pp[5];
	double l=pp[6];
	double thfluid=pp[7];
	double df=pp[8];
	double cs=pp[9];
	double cf=pp[10];
	double dens=pp[11];
	double denf=pp[12];
	return(exp(-g)*cos(v*g+fasetl(t,d,df,g,thsample,m,w0,l)+fasear(t,d,df,g,thfluid,m,w0,cs,cf,dens,denf)));
}

double amplireal (double t, double d,double df, double thsample,double thfluid,double v,double m,double w0, double l, double cs, double cf, double dens, double denf)
{
  gsl_integration_workspace *w=gsl_integration_workspace_alloc (espaco);
  double result, error;
  double p[13];
  p[0]=t;
  p[1]=d;
  p[2]=thsample;
  p[3]=v;
  p[4]=m;
  p[5]=w0;
  p[6]=l;
  p[7]=thfluid;
  p[8]=df;
  p[9]=cs;
  p[10]=cf;
  p[11]=dens;
  p[12]=denf;
  gsl_function F;
  F.function=&j;
  F.params=p;
  gsl_integration_qag(&F,0,30,0,1e-5,espaco,1,w,&result,&error);
  gsl_integration_workspace_free (w);
  return (result);
}

double kk(double g, void * params)
{
	double *pp= (double *) params;
	double t=pp[0];
	double d=pp[1];
	double thsample=pp[2];
	double v=pp[3];
	double m=pp[4];
	double w0=pp[5];
	double l=pp[6];
	double thfluid=pp[7];
	double df=pp[8];
	double cs=pp[9];
	double cf=pp[10];
	double dens=pp[11];
	double denf=pp[12];
	return(exp(-g)*sin(v*g+fasetl(t,d,df,g,thsample,m,w0,l)+fasear(t,d,df,g,thfluid,m,w0,cs,cf,dens,denf)));/*Ai, ó: fasetm e fasear, de novo. Puta sacanagem*/
}

/*Integra k, obtendo a parte complexa da amplitude*/
double amplicomplex (double t, double d,double df, double thsample,double thfluid,double v,double m,double w0, double l,double cs,double cf,double dens,double denf)
{
  gsl_integration_workspace *w=gsl_integration_workspace_alloc (espaco);
  double result, error;
  double p[13];
  p[0]=t;
  p[1]=d;
  p[2]=thsample;
  p[3]=v;
  p[4]=m;
  p[5]=w0;
  p[6]=l;
  p[7]=thfluid;
  p[8]=df;
  p[9]=cs;
  p[10]=cf;
  p[11]=dens;
  p[12]=denf;
  gsl_function F;
  F.function=&kk;
  F.params=p;
  gsl_integration_qag(&F,0,30,0,1e-5,espaco,1,w,&result,&error);
  gsl_integration_workspace_free (w);
  return (result);
}

/*intensity corresponde ao módulo ao quadrado da amplitude complexa*/
double intensity (double t, double d,double df, double thsample,double thfluid,double v,double m,double w0,double l,double cs,double cf,double dens,double denf)
{
	double temp1;
	if (t<=0.000000001)
	{
	temp1=1;
	}
	else
	{
	temp1=(pow(amplireal(t,d,df,thsample,thfluid,v,m,w0,l,cs,cf,dens,denf),2)+pow(amplicomplex (t,d,df,thsample,thfluid,v,m,w0,l,cs,cf,dens,denf),2))*(1+v*v);
	}
	return(temp1);
}

int main (int argc, char *argv[])
{
   double thsample,thfluid,d,chi,m,w0,v,l,df,cs,cf,dens,denf;
   int i,n;
   char path[100];
   FILE *arq,*arqresult;

   d=atof(argv[2]);
   thsample=atof(argv[3]);
   thfluid=atof(argv[4]);
   m=atof(argv[5]);
   w0=atof(argv[6]);
   v=atof(argv[7]);
   l=atof(argv[8]);
   df=atof(argv[9]);
   cs=atof(argv[10]);
   cf=atof(argv[11]);
   dens=atof(argv[12]);
   denf=atof(argv[13]);
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
        gsl_matrix_set(ma,i,1,intensity(gsl_matrix_get(ma,i,0),d,df,thsample,thfluid,v,m,w0,l,cs,cf,dens,denf));
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
