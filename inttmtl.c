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

double l=632.8e-9,dndt=-1e-6,beta=17.3e-6,nu=0.31,df=2.2e-5,espaco=500;

double thetatm(double p,double r,double beta,double nu,double l,double k)
{
    return (-(p*(1-r)*beta*(1+nu))/(l*k));
}

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

double fasetm(double t, double d, double g, double thtm,double m,double w0)
{
  gsl_integration_workspace *wgsl=gsl_integration_workspace_alloc (espaco);
  double result, error;
  double p[6];
  p[0]=d;
  p[1]=t;
  p[2]=g;
  p[3]=thtm;
  p[4]=m;
  p[5]=w0;
  gsl_function F;
  F.function=&argfasetm;
  F.params=p;
  gsl_integration_qag(&F,0,1e5,0,1e-3,espaco,1,wgsl,&result,&error);
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
	double alpha=pp2[5];
	temp1=(sqrt(df/(M_PI*tau))/d)*(exp(-pow(alpha,2)*df*tau-pow(alpha*w0,2)/8));
	temp2=gsl_sf_erf(sqrt(d*(t-tau))*alpha);
	temp3=(gsl_sf_bessel_J0(sqrt(m*g)*w0*alpha)-1);
	return(temp1*temp2*temp3);
}

double argfasear(double alpha, void * params)
{
  gsl_integration_workspace *w=gsl_integration_workspace_alloc (espaco);
  double result, error;
  double ptau[6];
  gsl_function F;
  double *pp= (double *) params;
  ptau[0]=pp[0];
  ptau[1]=pp[1];
  ptau[2]=pp[2];
  ptau[3]=pp[3];
  ptau[4]=pp[4];
  ptau[5]=alpha;
  F.function=&argtau;
  F.params=ptau;
  gsl_integration_qag(&F,0,ptau[1],0,1e-3,espaco,1,w,&result,&error);
  gsl_integration_workspace_free (w);
  return (result);
}

double fasear(double t,  double d, double g, double thtm,double m,double w0)
{
  gsl_integration_workspace *w=gsl_integration_workspace_alloc (espaco);
  double result, error;
  double par[5];
  par[0]=d;
  par[1]=t;
  par[2]=g;
  par[3]=m;
  par[4]=w0;
  gsl_function F;
  F.function=&argfasear;
  F.params=par;
  gsl_integration_qag(&F,0,1e5,0,1e-3,espaco,1,w,&result,&error);
  gsl_integration_workspace_free (w);
  return (-thtm*(2*d*dndt)/(beta*(1 + nu))*result);
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
	return(exp(-g)*cos(v*g+fasetm(t,d,g,thtm,m,w0)+fasear(t,d,g,thtm,m,w0)));
}

double amplireal (double t, double d, double thtm,double v,double m,double w0)
{
  gsl_integration_workspace *w=gsl_integration_workspace_alloc (espaco);
  double result, error;
  double p[6];
  p[0]=t;
  p[1]=d;
  p[2]=thtm;
  p[3]=v;
  p[4]=m;
  p[5]=w0;
  gsl_function F;
  F.function=&j;
  F.params=p;
  gsl_integration_qag(&F,0,20,0,1e-3,espaco,1,w,&result,&error);
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
	return(exp(-g)*sin(v*g+fasetm(t,d,g,thtm,m,w0)+fasear(t,d,g,thtm,m,w0)));/*Ai, ó: fasetm e fasear, de novo. Puta sacanagem*/
}

/*Integra k, obtendo a parte complexa da amplitude*/
double amplicomplex (double t, double d, double thtm,double v,double m,double w0)
{
  gsl_integration_workspace *w=gsl_integration_workspace_alloc (espaco);
  double result, error;
  double p[6];
  p[0]=t;
  p[1]=d;
  p[2]=thtm;
  p[3]=v;
  p[4]=m;
  p[5]=w0;
  gsl_function F;
  F.function=&k;
  F.params=p;
  gsl_integration_qag(&F,0,20,0,1e-3,espaco,1,w,&result,&error);
  gsl_integration_workspace_free (w);
  return (result);
}

/*intensityomp corresponde ao módulo ao quadrado da amplitude complexa*/
double intensityomp (double t, double d, double thtm,double v,double m,double w0)
{
	double temp[2],temp1;
	int i;
	if(t<=0.000000001)
	{
	temp1=1;
	}
	else
	{
	#pragma omp parallel for
         for(i=0;i<2;i++)
         {
		if(i==0)
		{
		temp[0]=pow(amplireal(t,d,thtm,v,m,w0),2);
		}
	  	else
		{
		temp[1]=pow(amplicomplex(t,d,thtm,v,m,w0),2);
		}
         }
         #pragma omp barrier
	temp1=(temp[0]+temp[1])*(1+v*v);
	}
	return(temp1);
}

/*intensity corresponde ao módulo ao quadrado da amplitude complexa*/
double intensity (double t, double d, double thtm,double v,double m,double w0)
{
	double temp1;
	if (t<=0.000000001)
	{
	temp1=1;
	}
	else
	{
	temp1=(pow(amplireal(t,d,thtm,v,m,w0),2)+pow(amplicomplex (t,d,thtm,v,m,w0),2))*(1+v*v);
	}
	return(temp1);
}

int main (int argc, char *argv[])
{  
   double th,d,m,w0,v,t;
   t=atof(argv[1]);
   d=atof(argv[2]);   
   th=atof(argv[3]);
   m=atof(argv[4]);
   w0=atof(argv[5]);
   v=atof(argv[6]);
   printf("%.10lf\n",intensityomp(t,d,th,v,m,w0));
   return (0);
}
