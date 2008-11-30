/************************************************************************

  program name: read_data.c

  funtion to read data

  status: in progress

  date: 12/5/89

*************************************************************************/

#include "hwe.h"

read_data ( a, no_allele, total, sample, infile )

int a[LENGTH];
int *no_allele, *total;
struct randomization *sample;
FILE **infile;

{
  register int i, j, l, err = 1;

  *total = 0;

  if( fscanf(*infile, "%d", no_allele) != 1) {
    fprintf(stderr, "Please supply number of alleles\n");
    return ( err );
  }

  if ( *no_allele < 3 ) {
    fprintf(stderr, "***Error! Number of alleles less than 3. \n");
    return ( err );
  }

  for ( i = 0; i < *no_allele; ++i ) {
    for ( j = 0; j <= i; ++j ) {
	 l = LL(i, j);
	 fscanf (*infile, "%d ", &a[l]);
	 *total += a[l];
    }
  }

  if( fscanf(*infile, "%d %d %d \n", &sample->step, 
                   &sample->group, &sample->size) != 3 ) {
    fprintf( stderr, " Please supply parameters.\n" );
    return ( err );
  }

  if ( sample->step < 1 || sample->group <= 1 ) {
    fprintf( stderr, "***Error in parameter specification.\n" );
    return ( err );
  }

  return ( 0 );

}
