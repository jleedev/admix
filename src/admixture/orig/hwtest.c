/*
 * hwtest.c
 *
 * Written by Thomas Dyer
 * Copyright (c) 2003, 2005
 *
 *
 * This program tests a marker locus for Hardy-Weinberg equilibrium. If
 * the marker has three or more alleles, an input file is created for the
 * program HWE (Guo and Thompson, 1992) which is then run to perform the
 * exact test of Hardy-Weinberg proportion for multiple alleles. For a
 * marker with fewer than three alleles, a chi-square test is performed.
 * There is a limit of 40 alleles per marker (MXALL). If necessary, this
 * limit can be increased in the source code and the program recompiled.
 *
 *
 * Usage:  hwtest [-au] [-p pop_id] [-o outfile] allfrq genfrq marker
 *
 *         allfrq       allele frequencies input file
 *         genfrq       genotype frequencies input file
 *         marker       marker to be tested
 *
 *         options:
 *           -a           include affecteds only
 *           -u           include unaffecteds only
 *           -p pop_id    include population pop_id only
 *           -o outfile   output file (required for program HWE)
 *
 *   Note that if the marker being tested has more than 2 alleles, in
 *   which case program HWE will be run, an output file must be specified.
 *
 *
 * File formats:
 *
 *   There are two blank-delimited input files - one for allele counts,
 *   the other for genotype counts. The first three fields on each line
 *   of these files are: marker name, affection status, and population ID.
 *   These fields are followed by an allele (or genotype) name, the number
 *   of occurrences of that allele (or genotype) among individuals within
 *   the specified population and having the specified affection status,
 *   and the corresponding allele (genotype) frequency. If the affection
 *   status field contains a hyphen (-), the count and frequency are for
 *   unaffected and affected indviduals combined. If the population ID
 *   field contains a hyphen, the count and frequency are for all
 *   populations combined. For example, the following line from a genotype
 *   counts file shows that there are 17 affected individuals in all
 *   populations combined who have the genotype 290/312 at marker D19S571,
 *   yielding a genotype frequency of 0.18478.
 *
 *       D19S571 a - 290 312 17 0.18478
 *
 *
 *   The HWE input file is created in free format. The first line of the
 *   file contains the number of alleles, n. The next n lines contain the
 *   lower triangular entries of the observed genotype matrix. The first of
 *   these lines has a single entry, the observed number of occurrences of
 *   the genotype a1/a1, where a1 is the first of the n alleles. The next
 *   line has two entries, the number of occurrences of the genotypes a2/a1
 *   and a2/a2. And so forth. The last line of the input file contains the
 *   number of dememorization steps, the number of batches, and the size of
 *   each batch.
 *
 *
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define MXALL	40	/* max # alleles per marker	*/

#define ALLBLK  256	/* # array elements malloc'd at a time		*/

#define DMSTEPS	2000	/* # dememorization steps	*/
#define NBATCH	1000	/* number of batches		*/
#define BSIZE	10000	/* batch size			*/


struct Marker {
    char *name;		/* marker name			*/
    int nall;		/* number of alleles		*/
    char **alleles;	/* allele names			*/
} ;

int nmrk;		/* number of markers		*/
struct Marker *mrk;

int **acnt;		/* marker allele counts		*/
int ***gcnt;		/* marker genotype counts	*/

int aff_only = 0, unaff_only = 0;
char *pop_id = 0;


void show_usage (char *);
void read_allfreq_file (char *);
void read_genfreq_file (char *);
int get_ndx (char *, char **, int);
int add_name (char *, char ***, int *);
void *allocMem (size_t);


main (int argc, char **argv)
{
    int i, j, line, rc;
    int mrkn, nall;
    double freq;
    double n, n_A, n_a, n_AA, n_Aa, n_aa;
    double chi, df = 1, pchis_();
    char cmd[1000];
    char *recp, rec[10000];
    char *ofile = 0;
    FILE *fpo;

    int errflg = 0;
    extern char *optarg;
    extern int optind, optopt;

    /* gather command line arguments */
    while ((i = getopt(argc, argv, ":aup:o:")) != -1) {
        switch (i) {
        case 'a':
            aff_only = 1;
            break;
        case 'u':
            unaff_only = 1;
            break;
        case 'p':
            pop_id = optarg;
            break;
        case 'o':
            ofile = optarg;
            break;
        case ':':
            fprintf(stderr, "option -%c requires an operand\n", optopt);
            errflg++;
            break;
        case '?':
            fprintf(stderr, "unrecognized option: -%c\n", optopt);
            errflg++;
        }
    }

    if (argc - optind != 3 || errflg) {
        show_usage(argv[0]);
        exit(1);
    }

    read_allfreq_file(argv[optind]);

    mrkn = -1;
    for (i = 0; i < nmrk; i++) {
        if (!strcmp(argv[optind+2], mrk[i].name)) {
            mrkn = i;
            break;
        }
    }
    if (mrkn == -1) {
        fprintf(stderr, "marker %s not found in allele frequency file\n",
                argv[optind+2]);
        exit(1);
    }

    if (!mrk[mrkn].nall) {
        fprintf(stderr,
    "no marker data available for specified population and affection status\n");
        exit(1);
    }

    read_genfreq_file(argv[optind+1]);

    /* perform chi-square test of HWE */
    if (mrk[mrkn].nall <= 2)
    {
        if (ofile) {
            fpo = fopen(ofile, "w");
            if (!fpo) {
                fprintf(stderr, "cannot open output file %s\n", ofile);
                exit(1);
            }
        }
        else
            fpo = stdout;

        n_AA = gcnt[mrkn][0][0];
        if (!n_AA)
            n_AA = 0.00000001;
        n_Aa = gcnt[mrkn][0][1] + gcnt[mrkn][1][0];
        if (!n_Aa)
            n_Aa = 0.00000001;
        n_aa = gcnt[mrkn][1][1];
        if (!n_aa)
            n_aa = 0.00000001;

        n = n_AA + n_Aa + n_aa;
        n_A = 2*n_AA + n_Aa;
        n_a = 2*n_aa + n_Aa;

        chi = -2*(n*log(n) + n_A*log(n_A) + n_a*log(n_a) + n_Aa*log(2.)
                  - 2*n*log(2*n) - n_AA*log(n_AA) - n_Aa*log(n_Aa)
                  - n_aa*log(n_aa));
        if (chi <= 0) {
            fprintf(fpo, "%s: chi = 0  p = 1\n", mrk[mrkn].name);
        }
        else {
            fprintf(fpo, "%s: chi = %g  p = %g\n", mrk[mrkn].name, chi,
                    1 - pchis_(&chi, &df));
        }
    }

    /* set up input files and run program HWE */
    else
    {
        if (!ofile) {
            fprintf(stderr, "an output file is required (-o option)\n");
            exit(1);
        }

        fpo = fopen("hwe.in", "w");
        if (!fpo) {
            fprintf(stderr, "cannot open hwe.in\n");
            exit(1);
        }

        nall = 0;
        for (i = 0; i < mrk[mrkn].nall; i++)
            if (acnt[mrkn][i]) nall++;
        fprintf(fpo, "%d\n", nall);

        for (i = 0; i < mrk[mrkn].nall; i++) {
            if (acnt[mrkn][i]) {
                for (j = 0; j < i; j++) {
                    if (acnt[mrkn][j]) {
                        fprintf(fpo, "%d ", gcnt[mrkn][i][j] + gcnt[mrkn][j][i]);
                    }
                }
                fprintf(fpo, "%d\n", gcnt[mrkn][i][i]);
            }
        }

        fprintf(fpo, "%d %d %d\n", DMSTEPS, NBATCH, BSIZE);
        fclose(fpo);

        sprintf(cmd, "hwe hwe.in %s\n", ofile);
        rc = system(cmd);
        if (rc) {
            fprintf(stderr, "%s: failed to run program HWE\n", argv[0]);
            exit(1);
        }
    }
}

void
show_usage (char *prog)
{
    printf("usage: %s [-au] [-p pop_id] allfrq genfrq marker\n\n", prog);
    printf("   allfrq       allele frequencies input file\n");
    printf("   genfrq       genotype frequencies input file\n");
    printf("   marker       marker to be tested\n");
    printf("\n   options:\n");
    printf("     -a           include affecteds only\n");
    printf("     -u           include unaffecteds only\n");
    printf("     -p pop_id    include population pop_id only\n");
    printf("     -o outfile   output file (required for program HWE)\n");
    exit(1);
}

void read_allfreq_file (char *frqfile)
{
    int i, j, k, line;
    int imrk, iall;
    char currmrk[1000];
    char *recp, rec[10000];
    FILE *fp;

    fp = fopen(frqfile, "r");
    if (!fp) {
        fprintf(stderr, "cannot open %s\n", frqfile);
        exit(1);
    }

    nmrk = 0;
    *currmrk = 0;
    line = 0;
    while (fgets(rec, sizeof(rec), fp)) {
        line++;
        if (!(recp = strtok(rec, " \t\n"))) {
            fprintf(stderr, "%s: missing marker name, line %d\n", frqfile,
                    line);
            exit(1);
        }
        if (strcmp(recp, currmrk)) {
            nmrk++;
            strcpy(currmrk, recp);
        }
    }

    mrk = (struct Marker *) allocMem(nmrk*sizeof(struct Marker));

    acnt = (int **) allocMem(nmrk*sizeof(int *));
    for (i = 0; i < nmrk; i++) {
        acnt[i] = (int *) allocMem(MXALL*sizeof(int));
        for (j = 0; j < MXALL; j++)
            acnt[i][j] = 0;
    }

    gcnt = (int ***) allocMem(nmrk*sizeof(int **));
    for (i = 0; i < nmrk; i++) {
        gcnt[i] = (int **) allocMem(MXALL*sizeof(int *));
        for (j = 0; j < MXALL; j++) {
            gcnt[i][j] = (int *) allocMem(MXALL*sizeof(int));
            for (k = 0; k < MXALL; k++)
                gcnt[i][j][k] = 0;
        }
    }

    rewind(fp);
    imrk = -1;
    *currmrk = 0;
    line = 0;
    while (fgets(rec, sizeof(rec), fp)) {
        line++;

        if (!(recp = strtok(rec, " \t\n"))) {
            fprintf(stderr, "%s: missing marker name, line %d\n", frqfile,
                    line);
            exit(1);
        }
        if (strcmp(recp, currmrk)) {
            imrk++;
            mrk[imrk].name =
                (char *) allocMem(strlen(rec));	/* rec contains \n */
            sscanf(rec, "%s", mrk[imrk].name);
            mrk[imrk].nall = 0;
            strcpy(currmrk, recp);
        }

        if (!(recp = strtok(NULL, " \t\n"))) {
            fprintf(stderr, "%s: missing affection status, line %d\n", frqfile,
                    line);
            exit(1);
        }

        if (unaff_only) {
            if (strcmp(recp, "U") && strcmp(recp, "1"))
                continue;
        }
        else if (aff_only) {
            if (strcmp(recp, "A") && strcmp(recp, "2"))
                continue;
        }
        else if (strcmp(recp, "-"))
            continue;

        if (!(recp = strtok(NULL, " \t\n"))) {
            fprintf(stderr, "%s: missing population identifier, line %d\n",
                    frqfile, line);
            exit(1);
        }

        if (pop_id) {
            if (strcmp(recp, pop_id))
                continue;
        }
        else if (strcmp(recp, "-"))
            continue;

        if (!(recp = strtok(NULL, " \t\n"))) {
            fprintf(stderr, "%s: missing marker allele, line %d\n", frqfile,
                    line);
            exit(1);
        }
        if ((iall = get_ndx(recp, mrk[imrk].alleles, mrk[imrk].nall)) == -1)
            iall = add_name(recp, &mrk[imrk].alleles, &mrk[imrk].nall);

        if (!(recp = strtok(NULL, " \t\n"))) {
            fprintf(stderr, "%s: missing allele count, line %d\n", frqfile,
                    line);
            exit(1);
        }
        acnt[imrk][iall] = atoi(recp);
    }

    fclose(fp);
}

void read_genfreq_file (char *frqfile)
{
    int i, j, line;
    int imrk, a1, a2;
    double gfreq;
    char *recp, rec[10000];
    FILE *fp;

    fp = fopen(frqfile, "r");
    if (!fp) {
        fprintf(stderr, "cannot open %s\n", frqfile);
        exit(1);
    }

    line = 0;
    while (fgets(rec, sizeof(rec), fp)) {
        line++;

        if (!(recp = strtok(rec, " \t\n"))) {
            fprintf(stderr, "%s: missing marker name, line %d\n", frqfile,
                    line);
            exit(1);
        }
        imrk = -1;
        for (i = 0; i < nmrk; i++) {
            if (!strcmp(mrk[i].name, recp)) {
                imrk = i;
                break;
            }
        }
        if (imrk == -1) {
            fprintf(stderr,
            "%s: marker %s not found in allele frequencies file, line %d\n",
                    frqfile, recp, line);
            exit(1);
        }

        if (!(recp = strtok(NULL, " \t\n"))) {
            fprintf(stderr, "%s: missing affection status, line %d\n", frqfile,
                    line);
            exit(1);
        }

        if (unaff_only) {
            if (strcmp(recp, "U") && strcmp(recp, "1"))
                continue;
        }
        else if (aff_only) {
            if (strcmp(recp, "A") && strcmp(recp, "2"))
                continue;
        }
        else if (strcmp(recp, "-"))
            continue;

        if (!(recp = strtok(NULL, " \t\n"))) {
            fprintf(stderr, "%s: missing population identifier, line %d\n",
                    frqfile, line);
            exit(1);
        }

        if (pop_id) {
            if (strcmp(recp, pop_id))
                continue;
        }
        else if (strcmp(recp, "-"))
            continue;

        if (!(recp = strtok(NULL, " \t\n"))) {
            fprintf(stderr, "%s: missing marker allele, line %d\n", frqfile,
                    line);
            exit(1);
        }
        a1 = get_ndx(recp, mrk[imrk].alleles, mrk[imrk].nall);
        if (a1 == -1) {
            fprintf(stderr,
    "%s: marker %s allele %s not found in allele frequencies file, line %d\n",
                    frqfile, mrk[imrk].name, recp, line);
            exit(1);
        }

        if (!(recp = strtok(NULL, " \t\n"))) {
            fprintf(stderr, "%s: missing marker allele, line %d\n", frqfile,
                    line);
            exit(1);
        }
        a2 = get_ndx(recp, mrk[imrk].alleles, mrk[imrk].nall);
        if (a2 == -1) {
            fprintf(stderr,
    "%s: marker %s allele %s not found in allele frequencies file, line %d\n",
                    frqfile, mrk[imrk].name, recp, line);
            exit(1);
        }

        if (!(recp = strtok(NULL, " \t\n"))) {
            fprintf(stderr, "%s: missing genotype count, line %d\n", frqfile,
                    line);
            exit(1);
        }
        gcnt[imrk][a1][a2] = atoi(recp);
    }

    fclose(fp);
}

int
get_ndx (char *str, char **array, int nelem)
{
    int i;

    for (i = 0; i < nelem; i++) {
        if (!strcmp(str, array[i]))
            return i;
    }

    return -1;
}

int
add_name (char *str, char ***array, int *nelem)
{
    int nblk;
    char **p;

    if (!*nelem) {
        *array = (char **) malloc(ALLBLK*sizeof(char *));
        if (!*array) {
            fprintf(stderr, "not enough memory\n");
            exit(1);
        }
    }
    else if (!(*nelem%ALLBLK)) {
        nblk = (*nelem)/ALLBLK + 1;
        p = (char **) realloc(*array, nblk*ALLBLK*sizeof(char *));
        if (!p) {
            fprintf(stderr, "not enough memory\n");
            exit(1);
        }
        *array = p;
    }

    (*array)[*nelem] = (char *) allocMem(strlen(str)+1);
    strcpy((*array)[*nelem], str);
    (*nelem)++;

    return *nelem - 1;
}

void
*allocMem (size_t nbytes)
{
    void *ptr;
    ptr = (void *) malloc(nbytes);
    if (!ptr) {
        fprintf(stderr, "not enough memory\n");
        exit(1);
    }
    return ptr;
}
