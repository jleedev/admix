/*
 * hetfst.c
 *
 * Written by Thomas Dyer
 * Copyright (c) 2003, 2004, 2005
 *
 *
 * This program calculates marker heterozygosity and Wright's locus-
 * specific F statistic, F_st. Marker heterozygosity and its variance
 * are computed for each population and for the total sample. F_st is
 * derived both analytically and by jackknifing, and is reported for
 * affected individuals, unaffected individuals, and the total sample.
 * Jackknifing also yields an estimate of the variance in F_st.
 *
 *
 * Usage:  hetfst [-M missval] pedfile genfrq hetout fstout
 *
 *         pedfile      pedigree file (marker genotypes)
 *         genfrq       genotype frequencies file
 *         hetout       heterozygosity output file
 *         fstout       F statistics output file
 *
 *         options:
 *           -M missval   missing allele value in quotes, e.g. "0"
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
 *   order specified in the locus file. Each marker genotype consists
 *   of a pair of alleles separated by a blank or tab. Individual,
 *   family, and population IDs can be arbitrary character strings but
 *   cannot contain embedded blanks or tabs. In the following example
 *   line from a pedigree file, there are two marker genotypes for a
 *   male Caucasian individual who is unaffected.
 *
 *       Blow Joe 1 Caucasian M A a 123 131
 *
 *
 *   The second input file contains genotype frequencies. The first
 *   three fields on each line of this file are: marker name, affection
 *   status, and population ID. These fields are followed by the names
 *   of the pair of alleles which make up the genotype, the number of
 *   occurrences of the genotype among individuals within the specified
 *   population and having the specified affection status, and the
 *   corresponding genotype frequency. If the affection status field
 *   contains a hyphen (-), the count and frequency are for unaffected
 *   and affected individuals combined. If the population ID field
 *   contains a hyphen, the count and frequency are for all populations
 *   combined. For example, the following line shows that there are 17
 *   affected individuals in all populations combined who have the
 *   genotype 290/312 at marker D19S571, yielding a genotype frequency
 *   of 0.18478.
 *
 *       D19S571 a - 290 312 17 0.18478
 *
 *
 *   The first output file reports the heterozygosity and its variance
 *   for each marker, first for the total sample and then for each of
 *   the populations. The second output file gives the analytic and
 *   jackknifed estimates of F_st and its variance for each marker,
 *   first for all individuals in the sample and then for unaffecteds
 *   and affecteds.
 *
 *
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define MXPOP	3	/* max # populations		*/
#define MXALL	40	/* max # alleles per marker	*/

#define MISSVAL	"*"	/* missing value		*/

#define ALLBLK  256	/* # array elements malloc'd at a time		*/


struct Marker {
    char *name;		/* marker name			*/
    int nall;		/* number of alleles		*/
    char **alleles;	/* allele names			*/
} ;

struct Marker *mrk;

int nmrk;		/* number of markers		*/
double *het;		/* marker heterozygosity	*/
double **phet;		/* heterozygosity by population	*/
double *var;		/* variance of heterozygosity	*/

int nindt;		/* number of individuals		*/
int nindp[MXPOP];	/* # individuals by population		*/
int ninda[2];		/* # individuals by affection status	*/

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


void show_usage (char *);
void read_genfreq_file (char *);
void read_pedigree_file (char *, char *);
int calc_f (int, int, int, double *, double *, double *);
int get_ndx (char *, char **, int);
int add_name (char *, char ***, int *);
void *allocMem (size_t);


main (int argc, char **argv)
{
    int i, j, err;
    double fis, ujfis, *jfis, jfisvar;
    double fit, ujfit, *jfit, jfitvar;
    double fst, ujfst, *jfst, jfstvar;
    FILE *fpf, *fph;

    int errflg = 0;
    char missval[10] = "";
    extern char *optarg;
    extern int optind, optopt;

    /* gather command line arguments */
    while ((i = getopt(argc, argv, ":M:")) != -1) {
        switch (i) {
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

    if (argc - optind != 4 || errflg) {
        show_usage(argv[0]);
        exit(1);
    }

    if (!strlen(missval)) {
        strncpy(missval, MISSVAL, sizeof(missval)-1);
        missval[sizeof(missval)-1] = 0;
    }

    read_genfreq_file(argv[optind+1]);

    read_pedigree_file(argv[optind], missval);

    jfis = (double *) allocMem(nindt*sizeof(double));
    jfit = (double *) allocMem(nindt*sizeof(double));
    jfst = (double *) allocMem(nindt*sizeof(double));

    /* compute standard and jackknifed F statistics for total sample */
    fpf = fopen(argv[optind+3], "w");
    fprintf(fpf, "TOTAL SAMPLE (N = %d)\n", nindt);
    fprintf(fpf,
    "MARKER   FIS     FIT     FST     J_FIS   J_SE    J_FIT   J_SE    J_FST   J_SE\n");

    for (j = 0; j < nmrk; j++) {
        calc_f(j, -1, 0, &ujfis, &ujfit, &ujfst);	/* standard F stats */
        fis = 0;
        fit = 0;
        fst = 0;
        for (i = 0; i < nindt; i++) {
            /* jackknifed F stats */
            err = calc_f(j, i, 0, &jfis[i], &jfit[i], &jfst[i]);
            if (err) break;
            fis += jfis[i];
            fit += jfit[i];
            fst += jfst[i];
        }
        if (err) {
            fprintf(fpf,
                    "%-8s %7.4f %7.4f %7.4f Sample too small for jackknifing.\n",
                    mrk[j].name, ujfis, ujfit, ujfst);
        }
        else {
            jfisvar = 0;
            jfitvar = 0;
            jfstvar = 0;
            for (i = 0; i < nindt; i++) {
                jfisvar += pow(jfis[i] - fis/nindt, 2.);
                jfitvar += pow(jfit[i] - fit/nindt, 2.);
                jfstvar += pow(jfst[i] - fst/nindt, 2.);
            }
            fprintf(fpf,
                    "%-8s %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f\n",
                    mrk[j].name, ujfis, ujfit, ujfst,
                    nindt*ujfis - (nindt-1)*fis/nindt, sqrt((nindt-1)*jfisvar/nindt),
                    nindt*ujfit - (nindt-1)*fit/nindt, sqrt((nindt-1)*jfitvar/nindt),
                    nindt*ujfst - (nindt-1)*fst/nindt, sqrt((nindt-1)*jfstvar/nindt));
        }
    }

    /* compute standard and jackknifed F statistics for unaffecteds */
    fprintf(fpf, "\nUNAFFECTED (N = %d)\n", ninda[0]);
    fprintf(fpf,
    "MARKER   FIS     FIT     FST     J_FIS   J_SE    J_FIT   J_SE    J_FST   J_SE\n");

    for (j = 0; j < nmrk; j++) {
        calc_f(j, -1, 1, &ujfis, &ujfit, &ujfst);	/* standard F stats */
        fis = 0;
        fit = 0;
        fst = 0;
        for (i = 0; i < nindt; i++) {
            /* jackknifed F stats */
            if (aff[i] == 1) {
                err = calc_f(j, i, 1, &jfis[i], &jfit[i], &jfst[i]);
                if (err) break;
                fis += jfis[i];
                fit += jfit[i];
                fst += jfst[i];
            }
        }
        if (err) {
            fprintf(fpf, "%-8s %7.4f %7.4f %7.4f Sample too small for jackknifing.\n",
                    mrk[j].name, ujfis, ujfit, ujfst);
        }
        else {
            jfisvar = 0;
            jfitvar = 0;
            jfstvar = 0;
            for (i = 0; i < nindt; i++) {
                if (aff[i] == 1) {
                    jfisvar += pow(jfis[i] - fis/ninda[0], 2.);
                    jfitvar += pow(jfit[i] - fit/ninda[0], 2.);
                    jfstvar += pow(jfst[i] - fst/ninda[0], 2.);
                }
            }
            fprintf(fpf,
                    "%-8s %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f\n",
                    mrk[j].name, ujfis, ujfit, ujfst,
                    ninda[0]*ujfis - (ninda[0]-1)*fis/ninda[0],
                    sqrt((ninda[0]-1)*jfisvar/ninda[0]),
                    ninda[0]*ujfit - (ninda[0]-1)*fit/ninda[0],
                    sqrt((ninda[0]-1)*jfitvar/ninda[0]),
                    ninda[0]*ujfst - (ninda[0]-1)*fst/ninda[0],
                    sqrt((ninda[0]-1)*jfstvar/ninda[0]));
        }
    }

    /* compute standard and jackknifed F statistics for affecteds */
    fprintf(fpf, "\nAFFECTED (N = %d)\n", ninda[1]);
    fprintf(fpf,
    "MARKER   FIS     FIT     FST     J_FIS   J_SE    J_FIT   J_SE    J_FST   J_SE\n");

    for (j = 0; j < nmrk; j++) {
        calc_f(j, -1, 2, &ujfis, &ujfit, &ujfst);	/* standard F stats */
        fis = 0;
        fit = 0;
        fst = 0;
        for (i = 0; i < nindt; i++) {
            /* jackknifed F stats */
            if (aff[i] == 2) {
                err = calc_f(j, i, 2, &jfis[i], &jfit[i], &jfst[i]);
                if (err) break;
                fis += jfis[i];
                fit += jfit[i];
                fst += jfst[i];
            }
        }
        if (err) {
            fprintf(fpf,
                    "%-8s %7.4f %7.4f %7.4f Sample too small for jackknifing.\n",
                    mrk[j].name, ujfis, ujfit, ujfst);
        }
        else {
            jfisvar = 0;
            jfitvar = 0;
            jfstvar = 0;
            for (i = 0; i < nindt; i++) {
                if (aff[i] == 2) {
                    jfisvar += pow(jfis[i] - fis/ninda[1], 2.);
                    jfitvar += pow(jfit[i] - fit/ninda[1], 2.);
                    jfstvar += pow(jfst[i] - fst/ninda[1], 2.);
                }
            }
            fprintf(fpf,
                    "%-8s %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f\n",
                    mrk[j].name, ujfis, ujfit, ujfst,
                    ninda[1]*ujfis - (ninda[1]-1)*fis/ninda[1],
                    sqrt((ninda[1]-1)*jfisvar/ninda[1]),
                    ninda[1]*ujfit - (ninda[1]-1)*fit/ninda[1],
                    sqrt((ninda[1]-1)*jfitvar/ninda[1]),
                    ninda[1]*ujfst - (ninda[1]-1)*fst/ninda[1],
                    sqrt((ninda[1]-1)*jfstvar/ninda[1]));
        }
    }

    fclose(fpf);

    fph = fopen(argv[optind+2], "w");

    /* compute heterozygosity for total sample */
    fprintf(fph, "TOTAL SAMPLE\n");
    fprintf(fph, "MARKER   HETERO   VAR(HET) S.E.\n");
    for (i = 0; i < nmrk; i++) {
        var[i] = het[i]*(1 - het[i])/nindt;
        fprintf(fph, "%-8s %8.6f %8.6f %8.6f\n", mrk[i].name, het[i], var[i],
                sqrt(var[i]));
    }

    /* compute heterozygosity by population */
    for (i = 0; i < npop; i++) {
        fprintf(fph, "\nPOPULATION %d\n", i+1);
        fprintf(fph, "MARKER   HETERO   VAR(HET) S.E.\n");
        for (j = 0; j < nmrk; j++) {
            var[j] = phet[j][i]*(1 - phet[j][i])/nindp[i];
            fprintf(fph, "%-8s %8.6f %8.6f %8.6f\n", mrk[j].name, phet[j][i],
                    var[j], sqrt(var[j]));
        }
    }

    fclose(fph);
}

void
show_usage (char *prog)
{
    printf("usage: %s [-M missval] pedfile genfrq hetout fstout\n\n", prog);
    printf("   pedfile      pedigree file\n");
    printf("   genfrq       genotype frequencies file\n");
    printf("   hetout       heterozygosity output file\n");
    printf("   fstout       F statistics output file\n");
    printf("\n   options:\n");
    printf("     -M missval   missing allele value in quotes\n");
    exit(1);
}

void read_genfreq_file (char *frqfile)
{
    int i, j, line;
    int imrk, ipop, a1, a2;
    double gfreq;
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
    het = (double *) allocMem(nmrk*sizeof(double));
    var = (double *) allocMem(nmrk*sizeof(double));
    phet = (double **) allocMem(nmrk*sizeof(double));
    for (i = 0; i < nmrk; i++)
        phet[i] = (double *) allocMem(MXPOP*sizeof(double));
 
    for (i = 0; i < nmrk; i++)
        het[i] = 0;

    for (i = 0; i < MXPOP; i++) {
        nindp[i] = 0;
        for (j = 0; j < nmrk; j++)
            phet[j][i] = 0;
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
        if (strcmp(recp, "-"))
            continue;

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
        if ((a1 = get_ndx(recp, mrk[imrk].alleles, mrk[imrk].nall)) == -1)
            a1 = add_name(recp, &mrk[imrk].alleles, &mrk[imrk].nall);

        if (!(recp = strtok(NULL, " \t\n"))) {
            fprintf(stderr, "%s: missing marker allele, line %d\n", frqfile,
                    line);
            exit(1);
        }
        if ((a2 = get_ndx(recp, mrk[imrk].alleles, mrk[imrk].nall)) == -1)
            a2 = add_name(recp, &mrk[imrk].alleles, &mrk[imrk].nall);

        if (a1 == a2)
            continue;

        if (!(recp = strtok(NULL, " \t\n"))) {
            fprintf(stderr, "%s: missing genotype count, line %d\n", frqfile,
                    line);
            exit(1);
        }

        if (!(recp = strtok(NULL, " \t\n"))) {
            fprintf(stderr, "%s: missing genotype frequency, line %d\n",
                    frqfile, line);
            exit(1);
        }
        gfreq = atof(recp);

        if (ipop == -1)
            het[imrk] += gfreq;
        else {
            phet[imrk][ipop] += gfreq;
            nindp[ipop]++;
        }
    }

    fclose(fp);
}

void read_pedigree_file (char *pedfile, char *missval)
{
    int i, j, line;
    int famid, id, sex;
    char *recp, rec[1024];
    FILE *fp;

    fp = fopen(pedfile, "r");
    if (!fp) {
        fprintf(stderr, "cannot open pedigree file %s\n", pedfile);
        exit(1);
    }

    nindt = 0;
    while (fgets(rec, sizeof(rec), fp)) {
        nindt++;
    }

    ninda[0] = ninda[1] = 0;

    pop = (int *) allocMem(nindt*sizeof(int));
    aff = (int *) allocMem(nindt*sizeof(int));

    all1 = (int **) allocMem(nindt*sizeof(int *));
    all2 = (int **) allocMem(nindt*sizeof(int *));
    for (i = 0; i < nindt; i++) {
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
        ninda[aff[i]-1]++;

        if (!(recp = strtok(NULL, " \t\n"))) {
            fprintf(stderr, "%s: missing population identifier, line %d\n",
                    pedfile,
                   line);
            exit(1);
        }
        if ((pop[i] = get_ndx(recp, pops, npop)) == -1) {
            fprintf(stderr,
        "%s: population identifier %s not found in genotype frequencies file, line %d\n",
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
            if (!(recp = strtok(NULL, " \t\n"))) {
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
                fprintf(stderr,
    "%s: marker %s allele %s not found in genotype frequencies file, line %d\n",
                        pedfile, mrk[j].name, recp, line);
                exit(1);
            }

            if (!(recp = strtok(NULL, " \t\n"))) {
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
                fprintf(stderr,
    "%s: marker %s allele %s not found in genotype frequencies file, line %d\n",
                        pedfile, mrk[j].name, recp, line);
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

/*
 *
 *  Compute locus-specific F statistics, weighted by (1/#populations).
 *
 *     calc_f(m, n, iaff, *fis, *fit, *fst)
 *
 *        m     locus
 *        n     population (-1 = include all)
 *        iaff  0 = include all, 1 = unaffect only, 2 = affect only
 *        fis   F_is
 *        fit   F_it
 *        fst   F_st
 *
 *     Returns  0 if F statistics computed OK, non-zero otherwise
 *
 */

int calc_f (int m, int n, int iaff, double *fis, double *fit, double *fst)
{
    int i, j, ipop;
    int nidp[MXPOP];
    int acnt[MXPOP][MXALL];
    int gcnt[MXPOP][MXALL];

    double hetp0, hetps, hetpt;
    double sum, hets[MXPOP];
    double x, xbar;

    for (i = 0; i < npop; i++) {
        nidp[i] = 0;
        for (j = 0; j < MXALL; j++) {
            acnt[i][j] = 0;
            gcnt[i][j] = 0;
        }
    }

    /* count #individuals by population
       count #alleles and #homozygotes at locus m by population */

    for (i = 0; i < nindt; i++) {

        /* drop individual n */
        if (i == n)
            continue;

        /* include/exclude by affection status */
        if (iaff && aff[i] != iaff)
            continue;

        ipop = pop[i];
        nidp[ipop]++;
        if (all1[i][m] != -1)
            acnt[ipop][all1[i][m]]++;
        if (all2[i][m] != -1)
            acnt[ipop][all2[i][m]]++;
        if (all1[i][m] != -1 && all1[i][m] == all2[i][m])
            gcnt[ipop][all1[i][m]]++;
    }

    /* compute hetS_i for all alleles at locus m */
    for (i = 0; i < npop; i++) {
        if (nidp[i] == 0) {
            return 1;
        }
        sum = 0;
        for (j = 0; j < MXALL; j++) {
            sum += pow((double)acnt[i][j]/(2*nidp[i]), 2.);
        }
        hets[i] = 1 - sum;
    }

    /* compute het0' */
    sum = 0;
    for (i = 0; i < npop; i++) {
        for (j = 0; j < MXALL; j++) {
            sum += gcnt[i][j]/nidp[i];
        }
    }
    hetp0 = 1 - sum/npop;

    /* compute hetS' */
    sum = 0;
    for (i = 0; i < npop; i++) {
        sum += (2*nidp[i]*hets[i])/(2*nidp[i] - 1);
    }
    hetps = sum/npop;

    /* compute hetT' */
    hetpt = 0;
    for (j = 0; j < MXALL; j++) {
        sum = 0;
        xbar = 0;
        for (i = 0; i < npop; i++) {
            x = (double)acnt[i][j]/(2*nidp[i]);
            xbar += x;
            sum += x*(1 - x)/(2*nidp[i] - 1);
        }
        sum = sum/(npop*npop);
        xbar = xbar/npop;
        hetpt += xbar*(1 - xbar) + sum;
    }

    *fis = 1 - hetp0/hetps;
    *fit = 1 - hetp0/hetpt;
    *fst = 1 - hetps/hetpt;

    return 0;
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
