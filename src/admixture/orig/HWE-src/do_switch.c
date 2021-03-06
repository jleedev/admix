/*********************************************************************

  program name: do_switch.c

  function to make switch according to given switchability and switch
  type.

  status:  in progress

  date: 1/14/90

**********************************************************************/

#include "hwe.h"

void do_switch ( a, index, type )

int a[LENGTH];
Index index;
int type;

{
  register int k11, k22, k12, k21;

  k11 = L(index.i1, index.j1);
  k12 = L(index.i1, index.j2);
  k21 = L(index.i2, index.j1);
  k22 = L(index.i2, index.j2);


  if ( type == 0 ) {  /* D-switch */
    --a[k11];
    --a[k22];
    ++a[k12];
    ++a[k21];
  } else {     /* R-switch */
    ++a[k11];
    ++a[k22];
    --a[k12];
    --a[k21];
  }
}
