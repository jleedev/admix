/************************************************************************

  program name: cal_n.c

  funtion to calculate n(i)

  status:  OK

  date: 12/6/89

*************************************************************************/

#include "hwe.h"

void cal_n ( no_allele, a, n )

int no_allele;
int a[LENGTH];
int n[MAX_ALLELE];

{
  register int i, j, l;

  for ( i = 0; i < no_allele; ++i ) {
    l = LL(i, i);
    n[i] = a[l];
    
    for ( j = 0; j < no_allele; ++j ) {
	 l = L(i, j);
	 n[i] += a[l];
    }

  }

}
