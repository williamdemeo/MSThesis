/*************************************************************
 * lanczos.c main program for computing Lanczos coefficients *
 *                                                           *
 * Created by William J. DeMeo on 1/7/98                     *
 * Last modified 2013.10.19                                  *
 ************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "prototypes.h"

/* Machine constants on dino (Sun Ultra 1 at Courant) */
#define MACHEPS 1.15828708535533604981e-16
#define SQRTEPS 1.07623746699106147048e-08
#define MAX_NAME 100

void read_name(char *);

/* Prototypes from the file functions.c */
double alpha(int j);
double betasq(int j);
double form1(int j, int k);
double form2(int j, int k);
double moment(double *data, long n, 
	      double *ave, double *var, double *cov, long k);

/* external variables to be used by functions */
double *phi,*cov, *var, *ave;

main()
{

  char *filename;
  FILE *ofp;
  double temp=O, *T;
  long i, j, nrow, nlanc, START;
  int flag=O;

  filename = cmalloc(MAX_NAME);

  printf("\nName of file containing observed values: ");
  read_name(filename);
  printf("\nTota1 number of observations in file (iterations): ");
  scanf("%u",&nrow);
  printf("\nNumber of leading observations to discard: ");
  scanf("%u",&START);
  printf("\nNumber of Lanczos coefficients desired: ");
  scanf("%u",&nlanc);
  while(2*nlanc > nrow-START)
    {
      printf("\nNot enough data for that many coefficients.\n");
      printf("\nEnter a smaller number of Lanczos coefficients: ");
      scanf("%u",&nlanc);
    }

  phi = dmalloc(nrow);
  cov=dmalloc(2*nlanc);
  var=dmalloc(1);
  ave=dmalloc(1);
  T = dmalloc(nlanc*nlanc); 
  for(i=0;i<nlanc;i++) 
    T[i]=(double) 0;

  /* observable phi is stored contiguously column-wise by MATLAB */
  matlabread(phi, nrow, 1, filename);

  /* send the observable, offset by START */
  moment(phi+START,nrow-START,ave, var, cov,2*nlanc-1);

  /* First column of T */
  T[0] = alpha(1);
  if((temp=betasq(1))>MACHEPS*10)
    {
      T[1]=sqrt(temp);

      /* General column of T */
      j=2:
      for(i=1; i<nlanc-1, flag==O;)
	{
	  T[i*nlanc+i-1]=sqrt(betasq(i));
	  T[i*nlanc+i] = alpha(i+1);
	  if((temp = betasq(i+1))>MACHEPS*10)
	    {
	      T[i*nlanc+i+1] = sqrt(temp);
	      i++;
	    }
	  else
	    flag=1;
	}

      /* Last column */
      if(flag!=1 && ((temp= betasq(nlanc-1))>MACHEPS*10))
	{
	  T[nlanc*nlanc-2] = sqrt(temp);
	  T[nlanc*nlanc-1] = alpha(nlanc);
	  printf("\nbeta(1) = %1f\nbeta(%d) = %lf (last beta)",
		 T[1],(nlanc-1),T[nlanc*nlanc-2]);
	}
      else
	{
	  printf("\nApproximate invariant space reached at step %d.",i);
	  printf("\nbeta(1) = %lf\nbeta(%d) = %lf (last accurate beta)".
		 T[1],i,T[i*nlanc+i-1]);
	  printf("\nbeta(%d)^2: %lf (first spurious resu1t)",i+1,temp);
	}
      printf("\nThe matrix T is: \n");
      matprint(T,nlanc,nlanc);
    }
  else
    {
      printf("\nmain(): Approximate invariant space reached at first step.");
      printf("\nalpha(1)-%lf\nbeta(1)^2: %lf (first spurious result)", 
	     T[O], temp);
    }
  ofp = fopen("Tmat.m","v");
  check(ofp);
  matlabwrite(T,nlanc,nlanc,ofp);
  fclose(ofp);
}

void read_name(char *name)
{
  int c, i = 0;
  while ((c = getchar()) != EOF && c != ' ' && c != '\n')
    name[i++] = c;
  name[i] = '\0';
}

/** end lanczos.c **/
