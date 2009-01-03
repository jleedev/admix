/**************************************************************************

  program name: cal_const.c

  function to compute the constant part of the probability function for
  testing the H-W equilibrium

  constant = log N! - log (2N)! + sum(i) log n(i)! 

  status: OK

  date: 12/3/89

**************************************************************************/


#include "hwe.h"

double cal_const ( no_allele, n, total )

int no_allele;
int n[MAX_ALLELE];
int total;

{
  double constant;
  register int i;
  double log_factorial();

  constant = log_factorial ( total ) - log_factorial ( 2*total );

  for ( i = 0; i < no_allele; ++i ) 
    constant += log_factorial ( n[i] );

  return ( constant );

}
