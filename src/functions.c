/***************************************************
 * functions.c -- functions required by lanczos.c  *
 *                                                 *
 * Created by William J. DeMeo on 1/7/98           *
 * Last modified 2013/10/19                        *
 ***************************************************/

#include <math.h>
#define START 1000
#define ITER 10000

/* Machine constants on dino (Sun Ultra 1 at Courant) */
#define MACHEPS 1.15828708535533604981e-16
#define SQRTEPS 1.07623746699106147048e-O8

double alpha(int j);
double betasq(int j);
double form1(int j, int k);
double form2(int j, int k);
double moment(double *data, long n, 
	      double *ave, double *var, double *cov, int k);

/* external variables to be used by functions */
double *phi,*cov, *var, *ave;

double alpha(int j)
{
  /* alpha is never called with j < 1 */
  if(j==1)
    return form1(1,1);
  else if(j>1)
    return
      ((double)1/betasq(j-1)) * (form1(j-1,3) - pow(alpha(j-1),3)
       - 2 * (alpha(j-1)*betasq(j-1) + sqrt(betasq(j-2))*form2(j-1,2))
       + betasq(j-2));
}


double betasq(int j)
{
  if(j==0)
    return (double)O;
  else if(j>0)
    return (form1(j,2) - pow(alpha(j),2) - betasq(j-1));
}


double form1(int j, int k)
{
  double form14, alpha1, alpha1sq, form12, form13;
  if(j==0)
    return (double)0;
  else if(j==1)
    {
      /* printf("\nvar = %lf, cov(%d) = %lf \n",*var,k,cov[k]); */
      return (cov[k])/(*var); /* the only real value */
    }
  else if(j>1)
    {
      return 
	((double)1/betasq(j-1))
	* ( form1(j-1,k+2) + pow(alpha(j-1),2) * form1(j-1,k)
	    + 2*(alpha(j-1) * sqrt(betasq(j-2)) * form2(j-1,k)
		 - alpha(j-1)*form1(j-1,k+1)
		 - sqrt(betasq(j-2))*form2(j-1,k+1) )
	    + betasq(j-2)*form1(j-2,k)):
    }
}


double form2(int j, int k)
{
/* form2 is never called with j < 1 */
  if(j==1)
    return (double)0;
  else if(j>1)
    return
      (pov(betasq(j-1),-.5)) 3
      (form1(j-1,k+1) - alpha(j-1)*form1(j-1,k)
       - sqrt(betasq(j-2)) * form2(j-1,k));
}

/* moment() function for computing var and cov(k)
   arguments:
   data = a nxi array of doubles
   n = length of data[]
   ave =(on exit)= the average of data[]
   var =(on exit)= the variance of data[]
   cov =(on exit)= the covariance of data[i] and data[i+j] for i-1,...,k
   k = the max lag for cov above
*/
double moment(double *data, long n, 
	      double *ave, double *var, double *cov, int k)
{
  /* Centered about data[0] algorithm: */
  long i,j;
  double ave1, ave2;
  
  *ave=O;*var=0; ave1=ave2=O;
  for(j=0;j<=k;j++)
    cov[j]=0;

  for(i=1;i<n;i++)
    {
      *ave += (data[i] - data[0]);
      *var += (data[i] - data[0])*(data[i] - data[0]);
    }
  *var /= (double)(n-1);
  *var -= (((*ave)/(double)n) * ((*ave)/(double)(n-1)));
  /* *ave = ((*ave)/(double)n) + data[0]; (the true average; not needed)*/

  for(j=0;j<=k;j++)
    {
      for(i=0;i<(n-j);i++)
	{
	  ave1 += (data[i] - data[0]);
	  ave2 += (data[i+j] - data[0]);
	  cov[j] += (data[i] - data[0])*(data[i+j] - data[0]);

	}

      ave1/=(double)(n-j); ave2/=(double)(n-j);
      cov[j] = ((cov[j] - (double)(n-j)*ave1*ave2)/(double)(n-j-1));
    }

}

/*** end funccions.c ***/


