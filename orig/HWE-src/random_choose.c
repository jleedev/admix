/* random_choose.c */

/***********************************************************************

  function to randomly choose two integer numbers, k1 and k2, between 0 
  and k - 1.  ( 0 <= k1 < k2 < k )

  status: OK

  date: 9/21/89

************************************************************************/

#include "hwe.h"

void random_choose ( k1, k2, k )
int *k1, *k2, k;

{
  register int temp, i, not_find;
  double drand48();
  int work[MAX_ALLELE];

  for ( i = 0; i < k; ++i )
    work[i] = i;

   *k1 = drand48() * k;
   
   --k;

  for ( i = *k1; i < k; ++i )
    work[i] = i + 1;

  not_find = 1;
  
  while ( not_find ) {
    i = drand48() * k;
    *k2 = work[i];
    not_find = 0;
  }
    
  if ( *k1 > *k2 ) {
    temp = *k1;
    *k1 = *k2;
    *k2 = temp;
  }

}
