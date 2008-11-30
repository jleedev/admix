****************************************************************************
Makefile
****************************************************************************

NAME   = hwe

CFLAGS = -O
CC      = cc
LIB = -lm

OBJS = $(NAME).o cal_const.o cal_n.o cal_prob.o\
		 check_file.o do_switch.o\
		 ln_p_value.o log_factorial.o print_data.o\
		random_choose.o  read_data.o select_index.o\
		 stamp_time.o test_switch.o 

$(OBJS): hwe.h func.h
$(NAME): $(OBJS)
	$(CC) $(CFLAGS) -o $(NAME) $(OBJS) $(LIB)

lixc:
	lint $(NAME).c cal_const.c cal_n.c cal_prob.c\
	check_file.c do_switch.c\ 
	ln_p_value.c log_factorial.c print_data.c\
	random_choose.c  read_data.c select_index.c stamp_time.c\
	test_switch.c 

print: *.h *.o
	pr $? | lpr -Pstat
		 touch print

clear:
	rm -f *.o


/***************************************************************************

  func.h 

  header file containing the function names.

  status: OK

  date: 1/15/90

***************************************************************************/

char *ctime();

double rans();
double log_factorial();
double ln_p_value();
double cal_prob();
double cal_const();


void print_data();
void get_interval();
void select_index();
void cal_n();
void stamp_time();
/***************************************************************************

  program name: hwe.h

  header file for hwe.c

  status: 

  date: 1/14/90

***************************************************************************/

#include  <stdio.h>
#include  <math.h>
#include  <strings.h>

#define  MAX_ALLELE    20
#define  LENGTH        MAX_ALLELE * ( MAX_ALLELE + 1 ) / 2
#define  STR_END       '\0'

#define  MIN(x, y)     ((x) < (y)) ? (x) : (y)
#define  RATIO(u, v)   ( (double) (u) ) / ( 1.0 + (double) (v) ) 
#define  TRANS(x)     (MIN(1.0, x))/2.0  /* transition probability */

#define  LL(a, b)      a * ( a + 1 ) / 2  + b
#define  L(a, b)       ( a < b ) ? b*(b + 1)/2 + a : a*(a+1)/2 + b

#define  EXPECT(a,b,c) ((double) a) / ((double) c) * ((double) b) / 2.0


typedef struct _Index
{
  int i1;
  int i2;
  int j1;
  int j2;
  int type;
  double cst;
} Index;

struct outcome 
{
  double p_value;  /* mean p-value */
  double se;       /* standard error of the p-value */
  int swch_count[3];  /* switch counts for partial and full switch */
};

struct randomization 
{
  int group; /* total number of chunks */
  int size;  /* size of a chunk */
  int step;  /* number of steps to de-memerization */
};

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







/***************************************************************************

  program name: cal_prob.c

  function to calculate 

  status: 

  date: 12/7/89

***************************************************************************/

#include "hwe.h"

double cal_prob ( a, index, ln_p_old, actual_switch )

int a[LENGTH];
Index index;
double ln_p_old;
int *actual_switch;

{

  double p1_ratio, p2_ratio;
  register double ln_p_new;
  double rand_num;
  int switch_ind, type;
  double drand48();
  void test_switch(), do_switch();

  *actual_switch = 0;

/* determine the switchability and direction of switch for given face */

  test_switch( a, index, &switch_ind, &type, &p1_ratio, &p2_ratio);
/*  printf("%d %d %d swhind=%d\n",index.i1, index.i2,index.i3,switch_ind);*/

  switch (switch_ind)
    {
    case 0:           /* non-switchable */

      ln_p_new = ln_p_old;  /* retain the pattern, probability unchanged */
      break;

    case 1:           /* partially-switchable */

	 if ( type == 1 )
	   p1_ratio = p2_ratio;
      rand_num = drand48();      

      if ( rand_num < TRANS( p1_ratio ) ) {  /* switch w/ transition P TRANS */
           do_switch ( a, index, type );
           ln_p_new = ln_p_old + log (p1_ratio);  /* ln P_after-switch */
           *actual_switch = 1;
         } else                   /* remain the same w/ P = 1 - TRANS */
         ln_p_new = ln_p_old;           /* probability unchanged */
       break;     

    default:          /* fully switchable */
         rand_num = drand48(); 

         if ( rand_num <= TRANS(p1_ratio)) {
           do_switch ( a, index, 0 ); /* D-switch */
           ln_p_new = ln_p_old + log (p1_ratio);  /* ln P_after-switch */
           *actual_switch = 2;
         } else if ( rand_num <= TRANS(p1_ratio) + TRANS(p2_ratio) ) {
           do_switch ( a, index, 1 ); /* R-switch */
           ln_p_new = ln_p_old + log (p2_ratio);
           *actual_switch = 2;
         } else
           ln_p_new = ln_p_old;
         break;    
    }

  return (ln_p_new);
}

/* check_file.c */

/************************************************************************

  function to check execution command and file

  status:

  date: 11/14/89

*************************************************************************/

#include "hwe.h"

check_file ( argc, argv, infile, outfile )

int argc;
char *argv[];
FILE **infile, **outfile;
{

  int exit_value = 0;

  /* file manipulation */

  if ( argc != 3 ) {
    fprintf (stderr, "Bad commond.\nCorrect usage: hwe infile outfile.\n");
    exit_value = 1;
  }

  if ( ( *infile = fopen (argv[1], "r")) == (FILE *) NULL ) {
    fprintf (stderr, "Can't read %s\n", argv[1]);
    exit_value = 2;
  }

  if( ( *outfile = fopen (argv[2], "w")) == (FILE *) NULL ) {
    fprintf (stderr, "Can't write %s\n", argv[2]);
    exit_value = 3;
  }

  return (exit_value);

}
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
/*************************************************************************

  program name: ln_p_value.c

  function to compute the log p value for given genetype frequencies

  status: 

  date: 12/6/89

*************************************************************************/

#include "hwe.h"

double ln_p_value ( a, no_allele, constant )

int a[LENGTH];
int no_allele;
double constant;

{
  register int i, j, l, temp;
  register double ln_prob;
  double log_factorial();
  
  ln_prob = constant;
  temp = 0;

  for ( i = 0; i < no_allele; ++i ) {
    for ( j = 0; j < i; ++j ) {
	 l = LL(i, j);
	 temp += a[l];
      ln_prob -= log_factorial ( a[l] );
    }

	 l = LL(i, i);
	 ln_prob -= log_factorial ( a[l] );
  }

  ln_prob += temp * log ( 2.0 );

  return ( ln_prob );

}

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



/***********************************************************************

  program name: select_index.c

  function to randomly choose three integers i1, i2, i3 with
  0 <= i1 < i2 <= no_allele, 0 <= i3 <= no_allele

  status: OK

  date: 1/14/90

************************************************************************/

#include "hwe.h"

void select_index ( index, no_allele )

Index *index;
int no_allele;

{

  void random_choose();

  int i1, i2, j1, j2;
  int k = 0;
  int l = 0;

/* generate row indices */

  random_choose ( &i1, &i2, no_allele );

  index->i1 = i1;
  index->i2 = i2;

/* generate column indices */

  random_choose ( &j1, &j2, no_allele );

  index->j1 = j1;
  index->j2 = j2;

/* calculate Delta = d(i1,j1) + d(i1,j2) + d(i2,j1) + d(i2,j2) */

  if ( i1 == j1 )
    ++k;

  if ( i1 == j2 )
    ++k;

  if ( i2 == j1 )
    ++k;

  if ( i2 == j2 )
    ++k;

  index->type = k;
  
  if ( ( i1 == j1 ) || ( i2 == j2 ) )
    ++l;

  index->cst = ( l == 1 ) ? pow(2.0, (double) k) : pow(2.0, - (double) k);
}
/***********************************************************************

  program name: stamp_time.c

  function to stamp the date and time the program is executed.


  Status:

  Date: 12/10/89 

************************************************************************/

#include "hwe.h"

void stamp_time ( t1, outfile)

long t1;
FILE **outfile;

{
  char *ctime();
  long t2, now;
  long time();
  
  time(&t2);
  t2 -= t1;
  time(&now);

  fprintf (*outfile, "\nTotal elapsed time: %d''\n", t2);
  fprintf (*outfile, "Date and time: %s\n", ctime(&now));

}
/***************************************************************************

  program name: test_switch.c

  function to determine the switchibility 
  
  switch_ind = 0 if non-switchable, 1 if partially-switchable, 2 if switchable.

  And it returns the switch type if switchable.

  switch_type = 0 if D-switchable and 1 if R-switchable;
  
  But switch_dir = 0 if switch_ind = 2.

  In addition, it returns the probability ratio if appropriete.

  status:

  date: 1/14/90

****************************************************************************/

#include "hwe.h"

void test_switch ( a, index, switch_ind, switch_type, p1_rt, p2_rt)

int a[LENGTH];
Index index;
int *switch_ind, *switch_type; /* switchability and type of switch */
double *p1_rt, *p2_rt; /* probability ratio */

{
  register int k11, k22, k12, k21;

  *switch_ind = 0;

  k11 = L(index.i1, index.j1);
  k22 = L(index.i2, index.j2);
  k12 = L(index.i1, index.j2);
  k21 = L(index.i2, index.j1);

  if ( index.type <= 1 ) { /* type = 0, 1 */
    if ( a[k11] > 0 && a[k22] > 0 ) {
	 *switch_ind = 1;
	 *switch_type = 0; /* D-switchable */
	 *p1_rt = RATIO(a[k11], a[k12]) *  RATIO(a[k22], a[k21]) * index.cst;
    }
    if ( a[k12] > 0 && a[k21] > 0 ) {
	 *switch_ind += 1;
	 *switch_type = 1; /* R-switchable */
      *p2_rt = RATIO(a[k12], a[k11]) *  RATIO(a[k21], a[k22]) / index.cst;
    }

  } else {                  /* type = 2 */
    if ( a[k11] > 0 && a[k22] > 0 ) {
      *switch_ind = 1;
      *switch_type = 0; /* D-switchable */
	 *p1_rt = RATIO(a[k11],a[k12] + 1.0)*RATIO(a[k22],a[k12]) * index.cst;
    }
    if ( a[k12] > 1 ) {
      *switch_ind += 1;
      *switch_type = 1; /* R-switchable */
      *p2_rt = RATIO(a[k12],a[k11]) * RATIO(a[k12] - 1,a[k22]) / index.cst;
    }
    
  }

}

