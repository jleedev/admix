/* log_factorial.c */

/************************************************************************

 function to calculate log ( k! )

 status: OK

 date: 8/16/89

************************************************************************/

#include  "hwe.h"

double log_factorial ( k )

int k;

{

  register double result; 

  if ( k == 0 )
    result = 0.0;
  else
    result = log( (double)k ) + log_factorial ( k - 1 );

  return ( result );

}
