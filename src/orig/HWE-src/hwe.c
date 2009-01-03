/***********************************************************************

  program name: hwe.c
  
  main program


  Status:

  Date: 1/14/90

************************************************************************/

#include "hwe.h"
#include "func.h"

main(argc,argv)  /* correct execution of this program:
                    hwe infile outfile               */

int argc;
char *argv[];

{

  int a[LENGTH], n[MAX_ALLELE];
  double ln_p_observed, ln_p_simulated, p_mean, p_square;
  double constant, p_simulated, total_step;
  int no_allele, total, counter, actual_switch;
  Index index;
  struct randomization sample;
  struct outcome result;
  FILE *infile, *outfile;
  long t1, time();

  register int i, j;

  if ( check_file ( argc, argv, &infile, &outfile ) )
    exit ( 1 );

  time(&t1);

  if ( read_data ( a, &no_allele, &total, &sample, &infile ) )
    exit (2);

  print_data ( a, no_allele, sample, &outfile );

/*  cal_n ( no_allele, a, n );

  constant = cal_const ( no_allele, n, total );

  ln_p_observed = ln_p_value ( a, no_allele, constant );  */

  ln_p_observed = 0.0;

  ln_p_simulated = ln_p_observed;  

  p_mean = p_square = (double) 0.0;

  result.p_value = result.se = (double) 0.0;  /* initialization */

  result.swch_count[0] = result.swch_count[1] = result.swch_count[2] = 0;

/*  cal_exp( n, total, no_allele, a); */

  for ( i = 0; i < sample.step; ++i ) {  /* de-memorization for given steps */

    select_index ( &index, no_allele );

/*  for ( i = 0; i < no_allele; ++i ) {

    for ( j = 0; j <= i; ++j ) {
         l = LL(i, j);
         printf("%4d", a[l]);
	  }
    
    printf ("\n");
  }      */

    ln_p_simulated = cal_prob(a, index, ln_p_simulated, &actual_switch);

    ++result.swch_count[actual_switch];
  }

  for ( i = 0; i < sample.group; ++i ) {  

    counter = 0;

    for ( j = 0; j < sample.size; ++j ) {

    select_index ( &index, no_allele );

    ln_p_simulated = cal_prob(a, index, ln_p_simulated, &actual_switch);

    if ( ln_p_simulated <= ln_p_observed )
	 ++counter;

    ++result.swch_count[actual_switch];

  }
    p_simulated = (double) counter  / sample.size;
    p_mean += p_simulated;
    p_square += p_simulated * p_simulated;

  }
  
  p_mean /= sample.group;
  result.p_value = p_mean;
  result.se = p_square / ((double) sample.group)/(sample.group - 1.0)
    - p_mean / ( sample.group - 1.0 ) * p_mean;
  result.se = sqrt ( result.se );

  total_step = sample.step + sample.group * sample.size;

  fprintf(outfile, "Randomization test P-value: %7.4g  (%7.4g) \n",
		result.p_value, result.se);
  fprintf(outfile, "Percentage of partial switches: %6.2f \n",
		result.swch_count[1] / total_step * 100);
  fprintf(outfile, "Percentage of full switches: %6.2f \n",
		result.swch_count[2] / total_step * 100);
  fprintf(outfile, "Percentage of all switches: %6.2f \n",
		(result.swch_count[1] + result.swch_count[2]) / total_step * 100 );

  stamp_time ( t1, &outfile );

}
