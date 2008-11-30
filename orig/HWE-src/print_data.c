/************************************************************************

  program name: print_data.c

  functio to print out the data

  status: 

  date: 12/7/89

*************************************************************************/

#include "hwe.h"

void print_data ( a, no_allele, sample, outfile )

int a[LENGTH];
int no_allele;
struct randomization sample;
FILE **outfile;


{

  register i, j, k, l;
  char line[120];

  line[0] = '-';

  k = 1;

  fprintf (*outfile, "Observed genotype frequencies: \n\n");
  
  for ( i = 0; i < no_allele; ++i ) {

    for ( j = k; j < k + 5; ++j ) 
	 line[j] = '-';

    line[j] = STR_END;
    k = j;

    fprintf (*outfile, "%s\n", line);

    fprintf (*outfile, "|");

    for ( j = 0; j <= i; ++j ) {
	 l = LL(i, j);
	 fprintf(*outfile, "%4d|", a[l]);
    }
    
    fprintf (*outfile, "\n");
  }

  fprintf (*outfile, "%s\n\n", line);
  fprintf (*outfile, "Total number of alleles: %2d\n\n", no_allele);

  fprintf(*outfile, "Number of initial steps: %d\n", sample.step);
  fprintf(*outfile, "Number of chunks: %d\n", sample.group);
  fprintf(*outfile, "Size of each chunk: %d\n\n", sample.size); 

}
