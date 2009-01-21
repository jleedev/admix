/*
 * ldtest.c
 *
 * Written by Thomas Dyer
 * Copyright (c) 2003, 2004, 2005
 *
 *
 * This program tests for allelic linkage disequilibrium between marker
 * loci. A test of 2-locus LD is performed for each ordered combination
 * of two markers among the set of markers in the allle frequencies file,
 * i.e. the first marker in the frequencies file is compared with the
 * second marker, the second marker is compared with the third marker,
 * the third marker with the fourth, and so on. If three or more markers
 * are specified, a test of 3-locus LD is also performed for each ordered
 * combination of three markers.
 *
 * Note that if the allele frequencies file includes markers located on
 * more than one chromosome, LD tests will be made across chromosomal
 * boundaries. Such tests should of course be ignored. The set of markers
 * tested does not have to include all markers in the freuqencies file,
 * nor do the combinations of markers tested have to follow frequencies
 * file order. A list of markers to be tested can be supplied and those
 * markers will be compared in the order in which they occur in the list.
 *
 *
 * Usage:  ldtest [-m mlist] [-M missval] pedfile allfreq outfile
 *
 *         pedfile      pedigree file (marker genotypes)
 *         allfreq      allele frequencies file
 *         outfile      output file
 *
 *         options:
 *           -m mlist     read markers to be tested from file mlist
 *           -M missval   missing allele value in quotes, e.g. "0"
 *
 *   The marker list file consists of marker names, one name per line.
 *
 *   The missing allele value is used in genotypes to denote an untyped
 *   allele. The default missing value is an asterisk (*).
 *
 *
 * File formats:
 *
 *   The pedigree file is blank- or tab-delimited and consists of one
 *   line per individual. Each line contains the following fields:
 *   family ID, individual ID, affection status (coded U/A or 1/2 for
 *   unaffected/affected), population ID, sex (coded M/F or 1/2 for
 *   male/female), and marker genotypes, one for each marker, in the
 *   order in which the markers appear in the allele frequencies file.
 *   Each marker genotype consists of a pair of alleles separated by
 *   a blank or tab. Individual, family, and population IDs can be
 *   arbitrary character strings but cannot contain embedded blanks or
 *   tabs. In the following example line from a pedigree file, there
 *   are two marker genotypes for a male Caucasian individual who is
 *   unaffected.
 *
 *       Blow Joe 1 Caucasian M A a 123 131
 *
 *
 *   The second input file contains allele frequencies. The first three
 *   fields on each line of this file are: marker name, affection status,
 *   and population ID. These fields are followed by an allele name, the
 *   number of occurrences of the allele among individuals within the
 *   specified population and having the specified affection status, and
 *   the corresponding allele frequency. If the affection status field
 *   contains a hyphen (-), the count and frequency are for unaffected
 *   and affected individuals combined. If the population ID field
 *   contains a hyphen, the count and frequency are for all populations
 *   combined. For example, the following line shows that, for marker
 *   D19S571, there are 35 copies of allele 290 among the unaffected
 *   individuals in all populations combined, which corresponds to an
 *   allele frequency of 0.07623.
 *
 *       D19S571 u - 290 35 0.07623
 *
 *
 *   The output file reports, for each combination of alleles from each
 *   set of two (or three) markers, the frequency of that combination of
 *   alleles, the allelic disequilibrium, the chi-square statistic for
 *   testing whether the disequilibrium differs significantly from zero,
 *   and 1 minus the p-value associated with that statistic.
 *
 *
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define PHASE_UNKNOWN

#define MXPOP	3	/* max # populations		*/
#define MXALL	40	/* max # alleles per marker	*/

#define MISSVAL	"*"	/* missing value		*/

#define ALLBLK  256	/* # array elements malloc'd at a time		*/


struct Marker {
    char *name;		/* marker name			*/
    int nall;		/* number of alleles		*/
    char **alleles;	/* allele names			*/
    int *all_sort;	/* sorted order of allele names	*/
} ;

struct Marker *mrk;

int nmrk;		/* number of markers		*/
double **afreq;		/* marker allele frequencies	*/
double **afrequ;	/* allele freqs - unaffecteds	*/
double **afreqa;	/* allele freqs - affecteds	*/
double ***afreqp;	/* allele freqs by population	*/
double ***afreqpu;	/* freqs by pop - unaffecteds	*/
double ***afreqpa;	/* freqs by pop - affecteds	*/

int nind;		/* number of individuals	*/
int *pop;		/* population			*/
int *aff;		/* affection status		*/
int **all1, **all2;	/* alleles at each marker	*/


/*
 *  Only one copy of each identifier is stored, i.e. each entry
 *  in the following arrays is unique.
 */

int npop = 0;		/* number of populations	*/
char **pops;		/* population identifiers	*/

int nfamid = 0;		/* number of family IDs		*/
char **famids;		/* family IDs			*/

int nid = 0;		/* number of individual IDs	*/
char **ids;		/* individual IDs		*/

int ntst, *mtst;
double **pi, **tau, ****d2;
double **piu, **tauu, ****d2u;
double **pia, **taua, ****d2a;
double ***pip, ***taup, *****d2p;
double ***pipu, ***taupu, *****d2pu;
double ***pipa, ***taupa, *****d2pa;
double ***d3;


void show_usage (char *);
void read_allfreq_file (char *);
void read_pedigree_file (char *, char *);
void do_2locus_tests (FILE *, int, int);
void do_3locus_tests (FILE *, int, int);
void setup_storage (double *****, double ***, double ***, double **);
int get_ndx (char *, char **, int);
int add_name (char *, char ***, int *);
void sort_names(char **, int, int *);
void *allocMem (size_t);

#define min(a,b)	((a) <= (b) ? (a) : (b))


main (int argc, char **argv)
{
    int i, j, k, l;
    int errflg = 0;
    char *mfile = 0;
    char missval[10] = "";
    char *recp, rec[10000];
    FILE *fp;
    extern char *optarg;
    extern int optind, optopt;

    /* gather command line arguments */
    while ((i = getopt(argc, argv, ":m:M:")) != -1) {
        switch (i) {
        case 'm':
            mfile = optarg;
            break;
        case 'M':
            strncpy(missval, optarg, sizeof(missval)-1);
            missval[sizeof(missval)-1] = 0;
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

    if (optind > argc - 2 || errflg) {
        show_usage(argv[0]);
    }

    if (!strlen(missval)) {
        strncpy(missval, MISSVAL, sizeof(missval)-1);
        missval[sizeof(missval)-1] = 0;
    }

    read_allfreq_file(argv[optind+1]);


    /* set up list of markers to be tested */
    mtst = (int *) allocMem(nmrk*sizeof(int));
    if (!mfile) {
        for (i = 0; i < nmrk; i++) {
            mtst[i] = i;
        }
        ntst = nmrk;
    }
    else {
        fp = fopen(mfile, "r");
        if (!fp) {
            fprintf(stderr, "cannot open marker list %s\n", mfile);
            exit(1);
        }

        ntst = 0;
        while (fgets(rec, sizeof(rec), fp)) {
            recp = strtok(rec, "\n");
            mtst[ntst] = -1;
            for (j = 0; j < nmrk; j++) {
                if (!strcmp(recp, mrk[j].name)) {
                    mtst[ntst] = j;
                    break;
                }
            }
            if (mtst[ntst] == -1) {
                fprintf(stderr, "marker %s not found in frequencies file\n",
                        recp);
                exit(1);
            }
            ntst++;
        }

        if (ntst == 0) {
            fprintf(stderr, "marker list is empty\n");
            exit(1);
        }

        fclose(fp);
    }

    if (ntst < 2) {
        fprintf(stderr, "at least two markers must be specified\n");
        exit(1);
    }

    read_pedigree_file(argv[optind], missval);

    /* allocate and initialize storage */
    setup_storage(&d2, &pi, &tau, afreq);
    setup_storage(&d2u, &piu, &tauu, afrequ);
    setup_storage(&d2a, &pia, &taua, afreqa);

    d2p = (double *****) allocMem(npop*sizeof(double ****));
    pip = (double ***) allocMem(npop*sizeof(double **));
    taup = (double ***) allocMem(npop*sizeof(double **));
    for (i = 0; i < npop; i++) {
        setup_storage(&d2p[i], &pip[i], &taup[i], afreqp[i]);
    }

    d2pu = (double *****) allocMem(npop*sizeof(double ****));
    pipu = (double ***) allocMem(npop*sizeof(double **));
    taupu = (double ***) allocMem(npop*sizeof(double **));
    for (i = 0; i < npop; i++) {
        setup_storage(&d2pu[i], &pipu[i], &taupu[i], afreqpu[i]);
    }

    d2pa = (double *****) allocMem(npop*sizeof(double ****));
    pipa = (double ***) allocMem(npop*sizeof(double **));
    taupa = (double ***) allocMem(npop*sizeof(double **));
    for (i = 0; i < npop; i++) {
        setup_storage(&d2pa[i], &pipa[i], &taupa[i], afreqpa[i]);
    }

    d3 = (double ***) allocMem(nmrk*sizeof(double **));
    for (i = 0; i < MXALL; i++) {
        d3[i] = (double **) allocMem(MXALL*sizeof(double *));
        for (j = 0; j < MXALL; j++) {
            d3[i][j] = (double *) allocMem(MXALL*sizeof(double));
        }
    }

    fp = fopen(argv[optind+2], "w");
    if (!fp) {
        fprintf(stderr, "cannot open output file %s\n", argv[optind+2]);
        exit(1);
    }

/*
 *  Perform LD test for all pairs of adjacent markers specified in
 *  the marker list file.
 */

    fprintf(fp, "2-LOCUS LD TESTS\n");
    fprintf(fp, "================\n\n");
    fprintf(fp, "TOTAL SAMPLE");
    do_2locus_tests(fp, -1, 0);
    fprintf(fp, "\nUNAFFECTED");
    do_2locus_tests(fp, -1, 1);
    fprintf(fp, "\nAFFECTED");
    do_2locus_tests(fp, -1, 2);

    for (i = 0; i < npop; i++) {
        fprintf(fp, "\nPOPULATION: %s\n", pops[i]);
        fprintf(fp, "\n   ALL");
        do_2locus_tests(fp, i, 0);
        fprintf(fp, "\n   UNAFFECTED", pops[i]);
        do_2locus_tests(fp, i, 1);
        fprintf(fp, "\n   AFFECTED", pops[i]);
        do_2locus_tests(fp, i, 2);
    }

    if (ntst == 2) {
        fclose(fp);
        exit(0);
    }

/*
 *  Perform LD test for the set of markers in a 3-marker "sliding window"
 *  moving from the beginning to the end of the list of markers.
 */
/*
    fprintf(fp, "\n\n3-LOCUS LD TESTS\n");
    fprintf(fp, "================\n\n");
    fprintf(fp, "TOTAL SAMPLE");
    do_3locus_tests(fp, -1, 0);
    fprintf(fp, "\nUNAFFECTED");
    do_3locus_tests(fp, -1, 1);
    fprintf(fp, "\nAFFECTED");
    do_3locus_tests(fp, -1, 2);

    for (i = 0; i < npop; i++) {
        fprintf(fp, "\nPOPULATION: %s\n", pops[i]);
        fprintf(fp, "\n   ALL");
        do_3locus_tests(fp, i, 0);
        fprintf(fp, "\n   UNAFFECTED", pops[i]);
        do_3locus_tests(fp, i, 1);
        fprintf(fp, "\n   AFFECTED", pops[i]);
        do_3locus_tests(fp, i, 2);
    }
*/
    fclose(fp);
}

void
show_usage (char *prog)
{
    printf("usage: %s [-m mlist] [-M missval] pedfile allfreq outfile\n\n", prog);
    printf("   pedfile      pedigree file\n");
    printf("   allfreq      allele frequencies file\n");
    printf("   outfile      output file\n");
    printf("\n   options:\n");
    printf("     -m mlist     read sets of markers from file mlist\n");
    printf("     -M missval   missing allele value in quotes\n");
    exit(1);
}

void read_allfreq_file (char *frqfile)
{
    int i, line;
    int imrk, iaff, ipop, iall;
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
            fprintf(stderr, "%s: missing marker name, line %d\n", frqfile, line);
            exit(1);
        }
        if (strcmp(recp, currmrk)) {
            nmrk++;
            strcpy(currmrk, recp);
        }
    }

    mrk = (struct Marker *) allocMem(nmrk*sizeof(struct Marker));

    afreq = (double **) allocMem(nmrk*sizeof(double *));
    for (imrk = 0; imrk < nmrk; imrk++)
        afreq[imrk] = (double *) allocMem(MXALL*sizeof(double));

    afrequ = (double **) allocMem(nmrk*sizeof(double *));
    for (imrk = 0; imrk < nmrk; imrk++)
        afrequ[imrk] = (double *) allocMem(MXALL*sizeof(double));

    afreqa = (double **) allocMem(nmrk*sizeof(double *));
    for (imrk = 0; imrk < nmrk; imrk++)
        afreqa[imrk] = (double *) allocMem(MXALL*sizeof(double));

    afreqp = (double ***) allocMem(MXPOP*sizeof(double **));
    for (ipop = 0; ipop < MXPOP; ipop++) {
        afreqp[ipop] = (double **) allocMem(nmrk*sizeof(double *));
        for (imrk = 0; imrk < nmrk; imrk++)
            afreqp[ipop][imrk] = (double *) allocMem(MXALL*sizeof(double));
    }

    afreqpu = (double ***) allocMem(MXPOP*sizeof(double **));
    for (ipop = 0; ipop < MXPOP; ipop++) {
        afreqpu[ipop] = (double **) allocMem(nmrk*sizeof(double *));
        for (imrk = 0; imrk < nmrk; imrk++)
            afreqpu[ipop][imrk] = (double *) allocMem(MXALL*sizeof(double));
    }

    afreqpa = (double ***) allocMem(MXPOP*sizeof(double **));
    for (ipop = 0; ipop < MXPOP; ipop++) {
        afreqpa[ipop] = (double **) allocMem(nmrk*sizeof(double *));
        for (imrk = 0; imrk < nmrk; imrk++)
            afreqpa[ipop][imrk] = (double *) allocMem(MXALL*sizeof(double));
    }

    rewind(fp);
    imrk = -1;
    *currmrk = 0;
    line = 0;
    while (fgets(rec, sizeof(rec), fp)) {
        line++;

        if (!(recp = strtok(rec, " \t\n"))) {
            fprintf(stderr, "%s: missing marker name, line %d\n", frqfile, line);
            exit(1);
        }
        if (strcmp(recp, currmrk)) {
            imrk++;
            mrk[imrk].name = (char *) allocMem(strlen(recp)+1);
            sscanf(rec, "%s", mrk[imrk].name);
            mrk[imrk].nall = 0;
            strcpy(currmrk, recp);
        }

        if (!(recp = strtok(NULL, " \t\n"))) {
            fprintf(stderr, "%s: missing affection status, line %d\n", frqfile,
                    line);
            exit(1);
        }
        if (!strcmp(recp, "-"))
            iaff = 0;
        else if (!strcmp(recp, "U"))
            iaff = 1;
        else if (!strcmp(recp, "A"))
            iaff = 2;

        if (!(recp = strtok(NULL, " \t\n"))) {
            fprintf(stderr, "%s: missing population identifier, line %d\n",
                    frqfile, line);
            exit(1);
        }
        if (!strcmp(recp, "-")) {
            ipop = -1;
        }
        else if ((ipop = get_ndx(recp, pops, npop)) == -1) {
            if (npop == MXPOP) {
                fprintf(stderr, "too many populations, MXPOP = %d\n", MXPOP);
                exit(1);
            }
            ipop = add_name(recp, &pops, &npop);
        }

        if (!(recp = strtok(NULL, " \t\n"))) {
            fprintf(stderr, "%s: missing marker allele, line %d\n", frqfile,
                    line);
            exit(1);
        }
        if ((iall = get_ndx(recp, mrk[imrk].alleles, mrk[imrk].nall)) == -1)
            iall = add_name(recp, &mrk[imrk].alleles, &mrk[imrk].nall);

        if (!(recp = strtok(NULL, " \t\n"))) {
            fprintf(stderr, "%s: missing allele count, line %d\n", frqfile, line);
            exit(1);
        }

        if (!(recp = strtok(NULL, " \t\n"))) {
            fprintf(stderr, "%s: missing allele frequency, line %d\n", frqfile,
                    line);
            exit(1);
        }
        if (ipop == -1) {
            if (iaff == 0)
                afreq[imrk][iall] = atof(recp);
            else if (iaff == 1)
                afrequ[imrk][iall] = atof(recp);
            else if (iaff == 2)
                afreqa[imrk][iall] = atof(recp);
        }
        else {
            if (iaff == 0)
                afreqp[ipop][imrk][iall] = atof(recp);
            else if (iaff == 1)
                afreqpu[ipop][imrk][iall] = atof(recp);
            else if (iaff == 2)
                afreqpa[ipop][imrk][iall] = atof(recp);
        }

    }

    for (i = 0; i < nmrk; i++) {
        mrk[i].all_sort = (int *) allocMem(mrk[i].nall*sizeof(int));
        sort_names(mrk[i].alleles, mrk[i].nall, mrk[i].all_sort);
    }

    fclose(fp);
}

void read_pedigree_file (char *pedfile, char *missval)
{
    int i, j, k, l, line;
    int famid, id, sex;
    char *recp, rec[10000];
    FILE *fp;

    fp = fopen(pedfile, "r");
    if (!fp) {
        fprintf(stderr, "cannot open pedigree file %s\n", pedfile);
        exit(1);
    }

    nind = 0;
    while (fgets(rec, sizeof(rec), fp)) {
        nind++;
    }

    pop = (int *) allocMem(nind*sizeof(int));
    aff = (int *) allocMem(nind*sizeof(int));

    all1 = (int **) allocMem(nind*sizeof(int *));
    all2 = (int **) allocMem(nind*sizeof(int *));
    for (i = 0; i < nind; i++) {
        all1[i] = (int *) allocMem(nmrk*sizeof(int));
        all2[i] = (int *) allocMem(nmrk*sizeof(int));
    }

    rewind(fp);
    i = 0;
    line = 0;
    while (fgets(rec, sizeof(rec), fp)) {
        line++;

        if (!(recp = strtok(rec, " \t\n"))) {
            fprintf(stderr, "%s: missing family ID, line %d\n", pedfile, line);
            exit(1);
        }
        if ((famid = get_ndx(recp, famids, nfamid)) == -1)
            famid = add_name(recp, &famids, &nfamid);

        if (!(recp = strtok(NULL, " \t\n"))) {
            fprintf(stderr, "%s: missing ID, line %d\n", pedfile, line);
            exit(1);
        }
        if ((id = get_ndx(recp, ids, nid)) == -1)
            id = add_name(recp, &ids, &nid);

        if (!(recp = strtok(NULL, " \t\n"))) {
            fprintf(stderr, "%s: missing affection status, line %d\n", pedfile,
                    line);
            exit(1);
        }
        if (!strcmp(recp, "U") || !strcmp(recp, "u") || !strcmp(recp, "1"))
            aff[i] = 1;
        else if (!strcmp(recp, "A") || !strcmp(recp, "a") || !strcmp(recp, "2"))
            aff[i] = 2;
        else {
            fprintf(stderr,
    "%s: invalid affection status [%s], line %d: must be coded U/A or 1/2\n",
                    pedfile, recp, line);
            exit(1);
        }

        if (!(recp = strtok(NULL, " \t\n"))) {
            fprintf(stderr, "%s: missing population identifier, line %d\n",
                    pedfile, line);
            exit(1);
        }
        if ((pop[i] = get_ndx(recp, pops, npop)) == -1) {
            fprintf(stderr,
    "%s: population identifier %s not found in allele frequencies file, line %d\n",
                    pedfile, recp, line);
            exit(1);
        }

        if (!(recp = strtok(NULL, " \t\n"))) {
            fprintf(stderr, "%s: missing sex code, line %d\n", pedfile, line);
            exit(1);
        }
        if (!strcmp(recp, "M") || !strcmp(recp, "m") || !strcmp(recp, "1"))
            sex = 1;
        else if (!strcmp(recp, "F") || !strcmp(recp, "f") || !strcmp(recp, "2"))
            sex = 2;
        else {
            fprintf(stderr,
            "%s: invalid sex code [%s], line %d: must be coded M/F or 1/2\n",
                    pedfile, recp, line);
            exit(1);
        }

        /* read in the marker alleles */
        for (j = 0; j < nmrk; j++) {
            if (!(recp = strtok(NULL, " \t\n")))
            {
                fprintf(stderr, "%s: missing allele, marker %s, line %d\n",
                        pedfile, mrk[j].name, line);
                exit(1);
            }
            if (!strcmp(recp, missval))
            {
                all1[i][j] = -1;
            }
            else if ((all1[i][j] =
                 get_ndx(recp, mrk[j].alleles, mrk[j].nall)) == -1)
            {
                fprintf(stderr, "%s: unknown allele %s, marker %s, line %d\n",
                        pedfile, recp, mrk[j].name, line);
                exit(1);
            }

            if (!(recp = strtok(NULL, " \t\n")))
            {
                fprintf(stderr, "%s: missing allele, marker %s, line %d\n",
                        pedfile, mrk[j].name, line);
                exit(1);
            }
            if (!strcmp(recp, missval))
            {
                all2[i][j] = -1;
            }
            else if ((all2[i][j] =
                 get_ndx(recp, mrk[j].alleles, mrk[j].nall)) == -1)
            {
                fprintf(stderr, "%s: unknown allele %s, marker %s, line %d\n",
                        pedfile, recp, mrk[j].name, line);
                exit(1);
            }

            if (mrk[j].nall > MXALL) {
                fprintf(stderr, "marker %s has too many alleles, MXALL = %d\n",
                        mrk[j].name, MXALL);
                exit(1);
            }
        }

        i++;
    }

    fclose(fp);
}

void do_2locus_tests (FILE *fp, int ipop, int iaff)
{
    int i, j, ii, jj, n, nfrq;
    int m1, m2, mm1, mm2;
    int acnt[2][MXALL];
    int hcnt2[MXALL][MXALL];
    double afrq[2][MXALL];
    double hfrq2[MXALL][MXALL];
    double ****td2, **tafreq, **tpi;
    double freq, sum, asum;
    double dmax, dprime;
    double chi, var, df = 1, pchis_();
    char buf[10], indent[10];

    if (ipop == -1) {
        if (iaff == 1) {
            tafreq = afrequ;
            tpi = piu;
            td2 = d2u;
        }
        else if (iaff == 2) {
            tafreq = afreqa;
            tpi = pia;
            td2 = d2a;
        }
        else {
            tafreq = afreq;
            tpi = pi;
            td2 = d2;
        }
        strcpy(indent, "   ");
    }
    else {
        if (iaff == 1) {
            tafreq = afreqpu[ipop];
            tpi = pipu[ipop];
            td2 = d2pu[ipop];
        }
        else if (iaff == 2) {
            tafreq = afreqpa[ipop];
            tpi = pipa[ipop];
            td2 = d2pa[ipop];
        }
        else {
            tafreq = afreqp[ipop];
            tpi = pip[ipop];
            td2 = d2p[ipop];
        }
        strcpy(indent, "      ");
    }

    n = 0;
    for (i = 0; i < nind; i++) {
        if (ipop != -1 && pop[i] != ipop)
            continue;
        if (iaff != 0 && aff[i] != iaff)
            continue;
        n++;
    }
    fprintf(fp, "  (N = %d)\n", n);

    if (!n) return;	/* empty subgroup, so nothing to do	*/

    /* test pairs of adjacent markers */
    for (mm1 = 0; mm1 < ntst - 1; mm1++) {
        for (mm2 = mm1 + 1; mm2 < mm1 + 2; mm2++) {
            m1 = mtst[mm1];
            m2 = mtst[mm2];

            for (i = 0; i < MXALL; i++) {
                acnt[0][i] = acnt[1][i] = 0;
                afrq[0][i] = afrq[1][i] = 0;
                for (j = 0; j < MXALL; j++) {
                    hcnt2[i][j] = 0;
                    hfrq2[i][j] = 0;
                }
            }

            for (i = 0; i < nind; i++) {
                if (ipop != -1 && pop[i] != ipop)
                    continue;
                if (iaff != 0 && aff[i] != iaff)
                    continue;
                if (all1[i][m1] != -1) {
                    if (all1[i][m2] != -1)
                        acnt[0][all1[i][m1]]++;
                        acnt[1][all1[i][m2]]++;
                        hcnt2[all1[i][m1]][all1[i][m2]]++;
#ifdef PHASE_UNKNOWN
                    if (all2[i][m2] != -1)
                        hcnt2[all1[i][m1]][all2[i][m2]]++;
#endif
                }
                if (all2[i][m1] != -1) {
#ifdef PHASE_UNKNOWN
                    if (all1[i][m2] != -1)
                        hcnt2[all2[i][m1]][all1[i][m2]]++;
#endif
                    if (all2[i][m2] != -1)
                        acnt[0][all2[i][m1]]++;
                        acnt[1][all2[i][m2]]++;
                        hcnt2[all2[i][m1]][all2[i][m2]]++;
                }
            }

            sum = 0;
            for (i = 0; i < mrk[m1].nall; i++) {
                sum += acnt[0][i];
            }
            asum = 0;
            for (i = 0; i < mrk[m1].nall; i++) {
                afrq[0][i] = acnt[0][i]/sum;
            }

            sum = 0;
            for (i = 0; i < mrk[m2].nall; i++) {
                sum += acnt[1][i];
            }
            asum = 0;
            for (i = 0; i < mrk[m2].nall; i++) {
                afrq[1][i] = acnt[1][i]/sum;
            }

            nfrq = 0;
            sum = 0;
            for (i = 0; i < mrk[m1].nall; i++) {
                for (j = 0; j < mrk[m2].nall; j++) {
                    sum += hcnt2[i][j];
                    nfrq++;
                }
            }

            dprime = 0;
            asum = 0;
            for (i = 0; i < mrk[m1].nall; i++) {
                if (!afrq[0][i]) continue;
                for (j = 0; j < mrk[m2].nall; j++) {
                    if (!afrq[1][j]) continue;
                    hfrq2[i][j] = hcnt2[i][j]/sum;
                    td2[m1][m2][i][j] = hfrq2[i][j] - afrq[0][i]*afrq[1][j];
                    if (td2[m1][m2][i][j] < 0)
                        dmax = min(afrq[0][i]*afrq[1][j],
                                   (1-afrq[0][i])*(1-afrq[1][j]));
                    else
                        dmax = min(afrq[0][i]*(1-afrq[1][j]),
                                   (1-afrq[0][i])*afrq[1][j]);
                    if (dmax)
                        dprime +=
                            afrq[0][i]*afrq[1][j]*fabs(td2[m1][m2][i][j]/dmax);
                }
            }

            fprintf(fp, "\n%sMARKERS: %s %s\n", indent, mrk[m1].name,
                    mrk[m2].name);
            fprintf(fp, "%sD' = %g\n", indent, dprime);
            fprintf(fp, "%sALL1  ALL2    H2FREQ      D2       CHI2     1-PVAL\n",
                    indent);

            for (ii = 0; ii < mrk[m1].nall; ii++) {
                i = mrk[m1].all_sort[ii];
                for (jj = 0; jj < mrk[m2].nall; jj++) {
                    j = mrk[m2].all_sort[jj];
                    if (tpi[m1][i] && tpi[m2][j]) {
                        chi = n * pow(td2[m1][m2][i][j],2.)/(tpi[m1][i]*tpi[m2][j]);
                        if (chi == 0) {
                            fprintf(fp, "%s%-5s %-5s  %8.6f %9.6f %9.6f  %8.6f\n",
                                    indent, mrk[m1].alleles[i], mrk[m2].alleles[j],
                                    hfrq2[i][j], td2[m1][m2][i][j], chi, 1.);
                        }
                        else {
                            fprintf(fp, "%s%-5s %-5s  %8.6f %9.6f %9.6f  %8.6f\n",
                                    indent, mrk[m1].alleles[i], mrk[m2].alleles[j],
                                    hfrq2[i][j], td2[m1][m2][i][j], chi,
                                    1 - pchis_(&chi, &df));
                        }
                    } else {
                        fprintf(fp, "%s%-5s %-5s  %8.6f %9.6f   ******    ******\n",
                                indent, mrk[m1].alleles[i], mrk[m2].alleles[j],
                                hfrq2[i][j], td2[m1][m2][i][j]);
                    }
                }
            }
        }
    }
}

void do_3locus_tests (FILE *fp, int ipop, int iaff)
{
    int i, j, k, ii, jj, kk, n, nfrq;
    int m, m1, m2, m3, mm1, mm2;
    int hcnt3[MXALL][MXALL][MXALL];
    double hfrq3[MXALL][MXALL][MXALL];
    double ****td2, **tafreq, **tpi, **ttau;
    double freq, sum, asum;
    double chi, var, df = 1, pchis_();
    char buf[10], indent[10];

    if (ipop == -1) {
        if (iaff == 1) {
            tafreq = afrequ;
            tpi = piu;
            ttau = tauu;
            td2 = d2u;
        }
        else if (iaff == 2) {
            tafreq = afreqa;
            tpi = pia;
            ttau = taua;
            td2 = d2a;
        }
        else {
            tafreq = afreq;
            tpi = pi;
            ttau = tau;
            td2 = d2;
        }
        strcpy(indent, "   ");
    }
    else {
        if (iaff == 1) {
            tafreq = afreqpu[ipop];
            tpi = pipu[ipop];
            ttau = taupu[ipop];
            td2 = d2pu[ipop];
        }
        else if (iaff == 2) {
            tafreq = afreqpa[ipop];
            tpi = pipa[ipop];
            ttau = taupa[ipop];
            td2 = d2pa[ipop];
        }
        else {
            tafreq = afreqp[ipop];
            tpi = pip[ipop];
            ttau = taup[ipop];
            td2 = d2p[ipop];
        }
        strcpy(indent, "      ");
    }

    n = 0;
    for (i = 0; i < nind; i++) {
        if (ipop != -1 && pop[i] != ipop)
            continue;
        if (iaff != 0 && aff[i] != iaff)
            continue;
        n++;
    }
    fprintf(fp, "  (N = %d)\n", n);

    if (!n) return;	/* empty subgroup, so nothing to do	*/

    for (m = 0; m < ntst - 2; m++) {
        m1 = mtst[m];
        m2 = mtst[m+1];
        m3 = mtst[m+2];
        fprintf(fp, "\n%sMARKERS: %s %s %s\n", indent, mrk[m1].name, mrk[m2].name,
                mrk[m3].name);
        fprintf(fp, "%sALL1  ALL2  ALL3    H3FREQ      D3       CHI2     1-PVAL\n",
                indent);

        for (i = 0; i < MXALL; i++) {
            for (j = 0; j < MXALL; j++) {
                for (k = 0; k < MXALL; k++) {
                    hcnt3[i][j][k] = 0;
                    hfrq3[i][j][k] = 0;
                }
            }
        }

        for (i = 0; i < nind; i++) {
            if (ipop != -1 && pop[i] != ipop)
                continue;
            if (iaff != 0 && aff[i] != iaff)
                continue;
            if (all1[i][m1] != -1) {
                if (all1[i][m2] != -1) {
                    if (all1[i][m3] != -1)
                        hcnt3[all1[i][m1]][all1[i][m2]][all1[i][m3]]++;
#ifdef PHASE_UNKNOWN
                    if (all2[i][m3] != -1)
                        hcnt3[all1[i][m1]][all1[i][m2]][all2[i][m3]]++;
#endif
                }
#ifdef PHASE_UNKNOWN
                if (all2[i][m2] != -1) {
                    if (all1[i][m3] != -1)
                        hcnt3[all1[i][m1]][all2[i][m2]][all1[i][m3]]++;
                    if (all2[i][m3] != -1)
                        hcnt3[all1[i][m1]][all2[i][m2]][all2[i][m3]]++;
                }
#endif
            }
            if (all2[i][m1] != -1) {
#ifdef PHASE_UNKNOWN
                if (all1[i][m2] != -1) {
                    if (all1[i][m3] != -1)
                        hcnt3[all2[i][m1]][all1[i][m2]][all1[i][m3]]++;
                    if (all2[i][m3] != -1)
                        hcnt3[all2[i][m1]][all1[i][m2]][all2[i][m3]]++;
                }
#endif
                if (all2[i][m2] != -1) {
#ifdef PHASE_UNKNOWN
                    if (all1[i][m3] != -1)
                        hcnt3[all2[i][m1]][all2[i][m2]][all1[i][m3]]++;
#endif
                    if (all2[i][m3] != -1)
                        hcnt3[all2[i][m1]][all2[i][m2]][all2[i][m3]]++;
                }
            }
        }

        nfrq = 0;
        sum = 0;
        for (i = 0; i < mrk[m1].nall; i++) {
            for (j = 0; j < mrk[m2].nall; j++) {
                for (k = 0; k < mrk[m3].nall; k++) {
                    if (hcnt3[i][j][k]) {
                        sum += hcnt3[i][j][k];
                        nfrq++;
                    }
                }
            }
        }

        asum = 0;
        for (ii = 0; nfrq && ii < mrk[m1].nall; ii++) {
            i = mrk[m1].all_sort[ii];
            for (jj = 0; nfrq && jj < mrk[m2].nall; jj++) {
                j = mrk[m2].all_sort[jj];
                for (kk = 0; nfrq && kk < mrk[m3].nall; kk++) {
                    k = mrk[m3].all_sort[kk];
                    if (hcnt3[i][j][k]) {
                        nfrq--;
                        if (nfrq) {
                            sprintf(buf, "%8.6f", hcnt3[i][j][k]/sum);
                            sscanf(buf, "%lf", &hfrq3[i][j][k]);
                            asum += hfrq3[i][j][k];
                        }
                        else
                            hfrq3[i][j][k] = 1 - asum;

                        d3[i][j][k] = hfrq3[i][j][k]
                                - tafreq[m1][i] * td2[m2][m3][j][k]
                                - tafreq[m2][j] * td2[m1][m3][i][k]
                                - tafreq[m3][k] * td2[m1][m2][i][j]
                                - tafreq[m1][i] * tafreq[m2][j] * tafreq[m3][k];
                        var = (tpi[m1][i] * tpi[m2][j] * tpi[m3][k]
                               + 6 * td2[m1][m2][i][j]
                                   * td2[m2][m3][j][k] * td2[m1][m3][i][k]
                               + tpi[m1][i] *
                                 (ttau[m2][j] * ttau[m3][k] * td2[m2][m3][j][k]
                                  - pow(td2[m2][m3][j][k],2.))
                               + tpi[m2][j] *
                                 (ttau[m1][i] * ttau[m3][k] * td2[m1][m3][i][k]
                                  - pow(td2[m1][m3][i][k],2.))
                               + tpi[m3][k] *
                                 (ttau[m1][i] * ttau[m2][j] * td2[m1][m2][i][j]
                                  - pow(td2[m1][m2][i][j],2.))
                               + d3[i][j][k] *
                                 (ttau[m1][i] * ttau[m2][j] * ttau[m3][k]
                                  - 2 * ttau[m1][i] * td2[m2][m3][j][k]
                                  - 2 * ttau[m2][j] * td2[m1][m3][i][k]
                                  - 2 * ttau[m3][k] * td2[m1][m2][i][j]
                                  - d3[i][j][k]))/n;
                        chi = pow(d3[i][j][k],2.)/var;
if (var < 0) {
printf("ta1i=%g ta2j=%g ta3k=%g\n", tafreq[m1][i], tafreq[m2][j], tafreq[m3][k]);
printf("tp1i=%g tp2j=%g tp3k=%g\n", tpi[m1][i], tpi[m2][j], tpi[m3][k]);
printf("tt1i=%g tt2j=%g tt3k=%g\n", ttau[m1][i], ttau[m2][j], ttau[m3][k]);
printf("td23jk=%g td13ik=%g td12ij=%g\n", td2[m2][m3][j][k], td2[m1][m3][i][k], td2[m1][m2][i][j]);
printf("n=%d d3=%g hf=%g var=%g chi=%g\n", n, d3[i][j][k], hfrq3[i][j][k], var, chi);
fflush(stdout);
}
                        if (chi == 0) {
                            fprintf(fp, "%s%-5s %-5s %-5s  %8.6f %9.6f %9.6f %.8g\n",
                                    indent, mrk[m1].alleles[i], mrk[m2].alleles[j],
                                    mrk[m3].alleles[k], hfrq3[i][j][k],
                                    d3[i][j][k], chi, 1.);
                        }
                        else {
                            fprintf(fp, "%s%-5s %-5s %-5s  %8.6f %9.6f %9.6f %.8g\n",
                                    indent, mrk[m1].alleles[i], mrk[m2].alleles[j],
                                    mrk[m3].alleles[k], hfrq3[i][j][k],
                                    d3[i][j][k], chi, 1 - pchis_(&chi, &df));
                        }
                    }
                }
            }
        }
    }
}

void setup_storage(double *****d2_ptr, double ***pi_ptr, double ***tau_ptr,
                   double **freq)
{
    int i, j, k;
    double **pi, **tau, ****d2;

    d2 = (double ****) allocMem(nmrk*sizeof(double ***));
    pi = (double **) allocMem(nmrk*sizeof(double *));
    tau = (double **) allocMem(nmrk*sizeof(double *));
    for (i = 0; i < nmrk; i++) {
        d2[i] = (double ***) allocMem(nmrk*sizeof(double **));
        for (j = 0; j < nmrk; j++) {
            d2[i][j] = (double **) allocMem(MXALL*sizeof(double *));
            for (k = 0; k < MXALL; k++) {
                d2[i][j][k] = (double *) allocMem(MXALL*sizeof(double));
            }
        }
        pi[i] = (double *) allocMem(MXALL*sizeof(double));
        tau[i] = (double *) allocMem(MXALL*sizeof(double));
        for (j = 0; j < mrk[i].nall; j++) {
            pi[i][j] = freq[i][j] * (1 - freq[i][j]);
            tau[i][j] = 1 - 2*freq[i][j];
        }
    }

    *d2_ptr = d2;
    *pi_ptr = pi;
    *tau_ptr = tau;
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
    if (!*nelem) {
        *array = (char **) malloc(ALLBLK*sizeof(char *));
        if (!*array) {
            fprintf(stderr, "not enough memory\n");
            exit(1);
        }
    }
    else if (!(*nelem)%ALLBLK) {
        *array = (char **) realloc(*array, ALLBLK*sizeof(char *));
        if (!*array) {
            fprintf(stderr, "not enough memory\n");
            exit(1);
        }
    }

    (*array)[*nelem] = (char *) allocMem(strlen(str)+1);
    strcpy((*array)[*nelem], str);
    (*nelem)++;

    return *nelem - 1;
}

void
sort_names (char **array, int nelem, int *order)
{
    int i, b, t, tmp, done;

    for (i = 0; i < nelem; i++)
        order[i] = i;

    b = nelem - 1;
    done = 0;
    while (!done) {
        t = -1;
        for (i = 0; i < b; i++) {
            if (strcmp(array[order[i]], array[order[i+1]]) > 0) {
                tmp = order[i+1];
                order[i+1] = order[i];
                order[i] = tmp;
                t = i;
            }
        }
        if (t == -1)
            done = 1;
        else
            b = t;
    }
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
