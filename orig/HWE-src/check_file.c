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
