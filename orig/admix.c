/*
 * admix.c
 *
 * Written by Thomas Dyer
 * Copyright (c) 2003, 2004, 2005
 *
 *
 * This program reads the marker genotype data for a set of individuals
 * and computes maximum likelihood estimates of population admixture.
 * These estimates are derived using the allele frequencies observed in
 * the parent populations.
 *
 *
 * Usage:  admix [-q] [-M missval] [-m mlist] [-g gridint] locfile
 *               pedfile outfile
 *
 *         locfile      locus file (marker info)
 *         pedfile      pedigree file (marker genotypes)
 *         outfile      output file
 *
 *         options:
 *           -M missval   missing allele value in quotes, e.g. "0"
 *           -m mlist     read sets of markers from file mlist
 *           -g gridint   starting interval for grid search
 *           -q           suppress warning messages
 *
 *   The missing allele value is used in genotypes to denote an untyped
 *   allele. The default missing value is an asterisk (*).
 *
 *   Maximum likelihood estimates are found by searching over a grid
 *   of admixture proportions. The maximization procedure begins with
 *   a relatively coarse grid, determines the grid point for which the
 *   likelihood is maximized, and then conducts a search in the vicinity
 *   of this point with a finer grid (the interval is divided by 10).
 *   This process is repeated until the grid interval is less than or
 *   equal to 0.001. The default starting grid interval is 0.01.
 *
 * File formats:
 *
 *   The locus file contains a set of lines for each marker. The first
 *   line in each set contains the marker name. The lines which follow
 *   contain per-population marker allele frequency information, one
 *   line per allele.  Each allele info line is blank- or tab-delimited
 *   and contains the allele name followed by one or more population
 *   allele frequencies.  The order of the allele frequencies is fixed,
 *   i.e. the i-th allele frequency on every allele info line is for the
 *   same population (the i-th population) and every allele info line
 *   must contain the same number of allele frequencies. Marker names
 *   and allele names can be arbitrary character strings but cannot
 *   contain embedded blanks or tabs.
 *
 *   The pedigree file is blank- or tab-delimited and consists of one
 *   line per individual. Each line contains the following fields:
 *   family ID, individual ID, affection status (coded U/A or 1/2 for
 *   unaffected/affected), population ID, sex (coded M/F or 1/2 for
 *   male/female), and marker genotypes, one for each marker, in the
 *   order specified in the locus file. Each marker genotype consists
 *   of a pair of alleles separated by a blank or tab. Individual,
 *   family, and population IDs can be arbitrary character strings but
 *   cannot contain embedded blanks or tabs.
 *
 *   In the example input files below, there are two marker genotypes
 *   for one individual.
 *
 *     Locus file:
 *       marker1
 *       D23S999
 *
 *     Pedigree file:
 *       Blow Joe 1 Caucasian M A a 123 131
 *
 *
 *   The output file reports the maximum likelihood estimate of
 *   population admixture, and its standard error, for each individual.
 *
 *
 */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define MIDLEN	20	/* max length of IDs		*/
#define MXALL	40	/* max # alleles per marker	*/
#define MXPOP	3	/* number of populations	*/

#define MISSVAL	"*"	/* missing value		*/

#define ALLBLK	256	/* # array elements malloc'd at a time	*/

#define TINY	0.0000001


struct Marker {
    char *name;			/* marker name			*/
    int nall;			/* number of alleles		*/
    char **alleles;		/* allele names			*/
    double afreq[MXPOP][MXALL];	/* by-population allele freqs	*/
    int all[2];			/* individual marker genotypes	*/
} ;

int npop;		/* number of populations	*/
int nmrk;		/* number of markers		*/
struct Marker *mrk;


void show_usage (char *);
void read_locus_file (char *, int);
int get_ndx (char *, char **, int);
int add_name (char *, char ***, int *);
void *allocMem (size_t);


main(int argc, char **argv)
{
    FILE *fpi, *fpm, *fpo;
    char *recp, rec[10000];

    int i, ii, j, k, l, line;
    int aff, sex, iall, done;
    int ntst, *mtst;
    char famid[MIDLEN+1], id[MIDLEN+1];

    double *delta, *m, *info, *work;
    int nose, *ipvt, noinfo;
    double *mlo, *mhi, *maxlm, *se;
    double sum, prob, loglike, maxlike;
    float gridint, startint = 0.01;

    int iter, npop1, one = 1;
    double denom, det[2];
    double diff, sum_m, sum_v;

    int errflg = 0, nowarn = 0;
    char *pedfile, *mfile = 0;
    char missval[10] = "";
    extern char *optarg;
    extern int optind, optopt;

    /* gather command line arguments */
    while ((i = getopt(argc, argv, ":qM:m:g:")) != -1) {
        switch (i) {
        case 'M':
            strncpy(missval, optarg, sizeof(missval)-1);
            missval[sizeof(missval)-1] = 0;
            break;
        case 'm':
            mfile = optarg;
            break;
        case 'q':
            nowarn = 1;
            break;
        case 'g':
            if (sscanf(optarg, "%f", &startint) != 1) {
                fprintf(stderr,
                        "option -%c requires a floating point operand\n", optopt);
                errflg++;
            }
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

    if (!strlen(missval)) {
        strncpy(missval, MISSVAL, sizeof(missval)-1);
        missval[sizeof(missval)-1] = 0;
    }

    read_locus_file(argv[optind], nowarn);

    /* set up list of markers to be tested */
    mtst = (int *) allocMem(nmrk*sizeof(int));
    if (!mfile) {
        for (i = 0; i < nmrk; i++) {
            mtst[i] = i;
        }
        ntst = nmrk;
    }
    else {
        fpm = fopen(mfile, "r");
        if (!fpm) {
            fprintf(stderr, "cannot open marker list %s\n", mfile);
            exit(1);
        }

        ntst = 0;
        while (fgets(rec, sizeof(rec), fpm)) {
            recp = strtok(rec, "\n");
            if (!recp) continue;	/* skip blank line */

            mtst[ntst] = -1;
            for (j = 0; j < nmrk; j++) {
                if (!strcmp(recp, mrk[j].name)) {
                    mtst[ntst] = j;
                    break;
                }
            }
            if (mtst[ntst] == -1) {
                fprintf(stderr, "marker %s not found in locus file\n",
                        argv[optind]);
                exit(1);
            }
            for (j = 0; j < ntst; j++) {
                if (mtst[ntst] == mtst[j]) {
                    fprintf(stderr,
    "marker %s occurs more than once in marker list, but was used only once\n",
                            recp);
                    break;
                }
            }
            if (j == ntst) {
                 ntst++;
            }
        }
        if (ntst == 0) {
            fprintf(stderr, "marker list is empty\n");
            exit(1);
        }
        fclose(fpm);
    }

    /* count number of individuals and allocate storage */
    pedfile = argv[optind+1];
    fpi = fopen(pedfile, "r");
    if (!fpi) {
        fprintf(stderr, "cannot open pedigree file %s\n", pedfile);
        exit(1);
    }

    npop1 = npop - 1;

    delta = (double *) allocMem(npop1*sizeof(double));
    m = (double *) allocMem(npop1*sizeof(double));
    info = (double *) allocMem(npop1*npop1*sizeof(double));
    ipvt = (int *) allocMem(npop1*sizeof(int));
    work = (double *) allocMem(npop1*sizeof(double));

    mlo = (double *) allocMem(npop1*sizeof(double));
    mhi = (double *) allocMem(npop1*sizeof(double));
    maxlm = (double *) allocMem(npop1*sizeof(double));
    se = (double *) allocMem(npop1*sizeof(double));

    fpo = fopen(argv[optind+2], "w");
    if (!fpo) {
        fprintf(stderr, "cannot open %s\n", argv[optind+2]);
        exit(1);
    }

    fprintf(fpo, "FAMID    ID       MLE_1 SE_1  ");
    for (i = 2; i <= npop; i++)
        fprintf(fpo, "MLE_%d SE_%d  ", i, i);
    fprintf(fpo, "\n");


/*
 *  For each individual in the pedigree file, read in the genotype data
 *  and compute ML estimate of that individual's admixture based on the
 *  population-specific allele frequencies in the locus file.
 */

    line = 0;
    while (fgets(rec, sizeof(rec), fpi)) {
        line++;
        if (!(recp = strtok(rec, " \t\n"))) {
            fprintf(stderr, "%s: missing family ID, line %d\n", pedfile, line);
            exit(1);
        }
        strncpy(famid, recp, MIDLEN);
        famid[MIDLEN] = 0;

        if (!(recp = strtok(NULL, " \t\n"))) {
            fprintf(stderr, "%s: missing ID, line %d\n", pedfile, line);
            exit(1);
        }
        strncpy(id, recp, MIDLEN);
        id[MIDLEN] = 0;

        if (!(recp = strtok(NULL, " \t\n")))
        {
            fprintf(stderr, "%s: missing affection status, line %d\n", pedfile,
                    line);
            exit(1);
        }
        if (!strcmp(recp, "U") || !strcmp(recp, "u") || !strcmp(recp, "1"))
            aff = 1;
        else if (!strcmp(recp, "A") || !strcmp(recp, "a") || !strcmp(recp, "2"))
            aff = 2;
        else
        {
            fprintf(stderr,
    "%s: invalid affection status [%s], line %d: must be coded U/A or 1/2\n",
                    pedfile, recp, line);
            exit(1);
        }

        if (!(recp = strtok(NULL, " \t\n")))
        {
            fprintf(stderr, "%s: missing population identifier, line %d\n",
                    pedfile, line);
            exit(1);
        }

        if (!(recp = strtok(NULL, " \t\n")))
        {
            fprintf(stderr, "%s: missing sex code, line %d\n", pedfile, line);
            exit(1);
        }
        if (!strcmp(recp, "M") || !strcmp(recp, "m") || !strcmp(recp, "1"))
            sex = 1;
        else if (!strcmp(recp, "F") || !strcmp(recp, "f") || !strcmp(recp, "2"))
            sex = 2;
        else
        {
            fprintf(stderr,
            "%s: invalid sex code [%s], line %d: must be coded M/F or 1/2\n",
                    pedfile, recp, line);
            exit(1);
        }

        noinfo = 1;
        for (i = 0; i < nmrk; i++) {
            for (j = 0; j < 2; j++) {
                if (!(recp = strtok(NULL, " \t\n"))) {
                    fprintf(stderr, "%s: missing allele, marker %d, line %d\n",
                            pedfile, i+1, line);
                    exit(1);
                }
                if (!strcmp(recp, missval)) {
                    iall = -1;
                }
                else {
                    iall = get_ndx(recp, mrk[i].alleles, mrk[i].nall);
                    if (iall == -1 && !nowarn) {
                        fprintf(stderr,
                "Warning: unknown allele %s for marker %s on line %d of %s\n",
                                recp, mrk[i].name, line, pedfile);
                        fprintf(stderr,
                "    Alleles not found in locus file are treated as missing.\n");
                    }
                }
                if (iall != -1) {
                    noinfo = 0;
                }
                mrk[i].all[j] = iall;
            }
        }

        if (noinfo) {
            fprintf(fpo, "%-8s %-8s NOINFO\n", famid, id);
            continue;
        }

        /* conduct grid search to find max likelihood */
        for (i = 0; i < npop1; i++) {
            m[i] = 0;
            mlo[i] = 0;
            mhi[i] = 1;
            maxlm[i] = 0;
        }

        maxlike = -1.e300;
        gridint = startint;
        done = 0;
        while (!done) {
            loglike = 0;
            for (ii = 0; ii < ntst; ii++) {
                i = mtst[ii];
                for (j = 0; j < 2; j++) {
                    if (mrk[i].all[j] == -1) continue;
                        prob = mrk[i].afreq[npop1][mrk[i].all[j]];
                    for (k = 0; k < npop1; k++) {
                        prob += m[k]*(mrk[i].afreq[k][mrk[i].all[j]] -
                                      mrk[i].afreq[npop1][mrk[i].all[j]]);
                    }
                    loglike += log(prob);
                }
            }

            if (loglike > maxlike) {
                maxlike = loglike;
                for (i = 0; i < npop1; i++)
                    maxlm[i] = m[i];
            }

            for (i = npop1 - 1; i >= 0; i--) {
                sum = 0;
                for (j = 0; j < npop1; j++) {
                    if (j != i)
                        sum += m[j];
                }
                m[i] += gridint;
                if (m[i] > mhi[i] + TINY || m[i] > 1 + TINY - sum)
                    m[i] = mlo[i];
                else
                    break;
            }
            if (i == -1) {
                if (gridint < 0.001 + TINY)
                    done = 1;
                else {
                   for (i = 0; i < npop1; i++) {
                       mlo[i] = maxlm[i] > gridint + TINY ?
                                maxlm[i] - gridint : 0;
                       mhi[i] = maxlm[i] < 1 - TINY - gridint ?
                                maxlm[i] + gridint : 1;
                       m[i] = mlo[i];
                   }
                   gridint *= 0.1;
                   maxlike = -1.e300;
                }
            }
        }

        for (i = 0; i < npop1; i++) {
            m[i] = maxlm[i];
            for (j = 0; j < npop1; j++) {
                info[i*npop1+j] = 0;
            }
        }

        /* compute information matrix */
        for (ii = 0; ii < ntst; ii++) {
            i = mtst[ii];
            for (j = 0; j < mrk[i].nall; j++) {
                denom = mrk[i].afreq[npop1][j];
                for (k = 0; k < npop1; k++) {
                    delta[k] = mrk[i].afreq[k][j] - mrk[i].afreq[npop1][j];
                    denom += m[k]*delta[k];
                }
                if (denom > TINY) {
                    for (k = 0; k < npop1; k++) {
                        for (l = 0; l < npop1; l++) {
                            info[k*npop1+l] += 2*delta[k]*delta[l]/denom;
                        }
                    }
                }
            }
        }

        fprintf(fpo, "%-8s %-8s", famid, id);
        sum_m = 0;
        for (i = 0; i < npop1; i++) {
            sum_m += m[i];
        }
        if (sum_m > 1) sum_m = 1;

        /* invert information matrix to get variance-covariance matrix */
        dgefa_(info, &npop1, &npop1, ipvt, &nose);
        if (!nose) {
            dgedi_(info, &npop1, &npop1, ipvt, det, work, &one);
            sum_v = 0;
            for (i = 0; i < npop1; i++) {
                se[i] = sqrt(info[i*npop]);
                sum_v += info[i*npop];
                for (j = 0; j < i; j++) {
                    sum_v -= 2*info[i*npop1+j];
                }
            }
        }

        for (i = 0; i < npop1; i++) {
            fprintf(fpo, " %5.3f", m[i]);
            if (!nose)
                fprintf(fpo, " %5.3f", se[i]);
            else
                fprintf(fpo, " ******");
        }
        fprintf(fpo, " %5.3f", 1 - sum_m);
        if (!nose)
            fprintf(fpo, " %5.3f", sqrt(sum_v));
        else
            fprintf(fpo, " ******", 1 - sum_m);
        fprintf(fpo, "\n");
        fflush(fpo);
    }

    fclose(fpo);
    fclose(fpi);
}

void
show_usage (char *prog)
{
    printf("usage: %s [-q] [-M missval] locfile pedfile outfile\n\n", prog);
    printf("   locfile      locus file\n");
    printf("   pedfile      pedigree file\n");
    printf("   outfile      output file\n");
    printf("\n   options:\n");
    printf("     -M missval   missing allele value in quotes\n");
    printf("     -m mlist     read sets of markers from file mlist\n");
    printf("     -g gridint   starting interval for grid search\n");
    printf("     -q           suppress warning messages\n");
    exit(1);
}

void read_locus_file (char *locfile, int nowarn)
{
    char *recp, rec[10000];
    char buf[10000];
    int i, j, k, line;
    int nfld, tnpop, iall;
    double freq[MXPOP];
    FILE *fp;

    fp = fopen(locfile, "r");
    if (!fp) {
        fprintf(stderr, "cannot open locus file %s\n", locfile);
        exit(1);
    }

    nmrk = 0;
    npop = 0;
    line = 0;
    while (fgets(rec, sizeof(rec), fp)) {
        line++;
        nfld = sscanf(rec, "%s %s", buf, buf);
        if (nfld == 1)
            nmrk++;
        else if (nfld == 2) {
            recp = strtok(rec, " \t\n");
            tnpop = 0;
            while (recp = strtok(NULL, " \t\n"))
                tnpop++;
            if (!npop)
                npop = tnpop;
            if (tnpop != npop) {
                fprintf(stderr,
                        "%s: inconsistent number of populations, line %d\n",
                        locfile, line);
                exit(1);
            }
        }
        else {
            fprintf(stderr, "%s: blank lines not allowed, line %d\n", locfile,
                    line);
            exit(1);
        }
    }

    if (npop > MXPOP) {
        fprintf(stderr, "%s: too many populations, MXPOP = %d\n", locfile,
                MXPOP);
        exit(1);
    }

    mrk = (struct Marker *) allocMem(nmrk*sizeof(struct Marker));

    rewind(fp);
    fgets(rec, sizeof(rec), fp);
    mrk[0].name = (char *) allocMem(strlen(rec));	/* rec contains \n */
    sscanf(rec, "%s", mrk[0].name);

    line = 1;
    for (i = 0; i < nmrk; i++) {
        mrk[i].nall = 0;
        for (j = 0; j < npop; j++) freq[j] = 0;

        while (fgets(rec, sizeof(rec), fp)) {
            line++;
            if (sscanf(rec, "%s %s", buf, buf) == 2) {
                recp = strtok(rec, " \t\n");
                iall = add_name(recp, &mrk[i].alleles, &mrk[i].nall);

                for (j = 0; j < npop; j++) {
                    if (!(recp = strtok(NULL, " \t\n")))
                    {
                        fprintf(stderr,
                                "%s: missing allele frequency, line %d\n",
                                locfile, line);
                        exit(1);
                    }
                    if (sscanf(recp, "%lf", &mrk[i].afreq[j][iall]) != 1 ||
                        mrk[i].afreq[j][iall] < 0 || mrk[i].afreq[j][iall] > 1)
                    {
                        fprintf(stderr,
                                "%s: invalid allele frequency, line %d\n",
                                locfile, line);
                        exit(1);
                    }
                    freq[j] += mrk[i].afreq[j][iall];
                }
            }
            else {
                mrk[i+1].name = (char *) allocMem(strlen(rec));
                sscanf(rec, "%s", mrk[i+1].name);
                break;
            }
        }

        for (j = 0; j < npop; j++) {
            if (!nowarn && (freq[j] < 1 - TINY || freq[j] > 1 + TINY)) {
                fprintf(stderr,
        "Warning: allele frequencies sum to %f for marker %s, population %d\n",
                       freq[j], mrk[i].name, j+1);
                fprintf(stderr, "   Frequencies being adjusted to sum to 1.\n");
            }
            for (k = 0; k < mrk[i].nall; k++) {
                mrk[i].afreq[j][k] /= freq[j];
            }
        }
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
