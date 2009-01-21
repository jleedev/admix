/*
 * inform.c
 *
 * Written by Thomas Dyer
 * Copyright (c) 2005
 *
 *
 * This program computes the informativeness of a set of markers for
 * estimating individual admixture from a set of ancestral populations.
 * The expected information with respect to a particular set of
 * admixture proportions is given by the determinant of a matirx with
 * elements equal to the negative expected second partial derivative
 * of the log likelihood function (the information matrix), which is
 * derived using the marker allele frequencies observed in the parent
 * populations.
 *
 *
 * Usage:  inform [-aq] [-p plist] [-m mlist] [-g gridint] [-o outfile]
 *                locfile [m1 ...]
 *
 *         locfile      population-specific allele frequencies
 *         m1 ...       admixture proportions (optional)
 *
 *         options:
 *           -p plist     read sets of admix proportions from file plist
 *           -m mlist     read sets of markers from file mlist
 *           -a           use all markers simultaneously
 *           -g gridint   grid interval
 *           -o outfile   output file
 *           -q           suppress warning messages
 *
 *   If admixture proportions are not specified on the command line,
 *   informativeness is calculated at each point in an equally-spaced
 *   grid of proportions.
 *
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
 */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define MXALL	40	/* max # alleles per marker	*/
#define MXPOP	3	/* number of populations	*/

#define ALLBLK	256	/* # array elements malloc'd at a time	*/

#define TINY	0.0000001

struct Marker {
    char *name;			/* marker name			*/
    int nall;			/* number of alleles		*/
    char **alleles;		/* allele names			*/
    double afreq[MXPOP][MXALL];	/* by-population allele freqs	*/
    int *all[2];		/* individual marker genotypes	*/
} ;

int npop;		/* number of populations	*/
int nmrk;		/* number of markers		*/
struct Marker *mrk;

float gridint = .01;
double *delta, *info, *work;
int *ipvt;

void compute_inform (int, int, double **, int, int *, FILE *);
void show_usage (char *);
void read_locus_file (char *, int);
int add_name (char *, char ***, int *);
void *allocMem (size_t);


main(int argc, char **argv)
{
    int i, j, line;
    int npop1;
    int grid, nprp, ntst, *nmtst, **mtst;
    int errflg = 0, useall = 0, nowarn = 0;
    double **m, sum;
    char *mfile = 0, *pfile = 0, *ofile = 0;
    char *recp, rec[10000];
    FILE *pfp, *mfp, *ofp;
    extern char *optarg;
    extern int optind, optopt;

    /* gather command line arguments */
    while ((i = getopt(argc, argv, ":aqm:p:g:o:")) != -1) {
        switch (i) {
        case 'a':
            useall = 1;
            break;
        case 'q':
            nowarn = 1;
            break;
        case 'm':
            mfile = optarg;
            break;
        case 'p':
            pfile = optarg;
            break;
        case 'o':
            ofile = optarg;
            break;
        case 'g':
            if (sscanf(optarg, "%f", &gridint) != 1) {
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

    if (optind == argc || errflg) {
        show_usage(argv[0]);
    }

    read_locus_file(argv[optind], nowarn);
    npop1 = npop - 1;

    nprp = 0;
    if (pfile || (argc - optind == npop)) {
        grid = 0;
        if (argc - optind == npop)
            nprp++;

        if (pfile) {
            pfp = fopen(pfile, "r");
            if (!pfp) {
                fprintf(stderr, "cannot open admix proportions file %s\n", pfile);
                exit(1);
            }
            while (fgets(rec, sizeof(rec), pfp)) {
                nprp++;
            }
            rewind(pfp);
        }
    }
    else if (argc - optind == 1) {
        grid = 1;
    }
    else {
        show_usage(argv[0]);
    }

    /* allocate and initialize storage */
    if (grid) {
        m = (double **) allocMem(sizeof(double *));
        m[0] = (double *) allocMem(npop1*sizeof(double));
        for (i = 0; i < npop1; i++) {
            m[0][i] = 0;
        }
    }
    else {
        m = (double **) allocMem(nprp*sizeof(double *));
        for (i = 0; i < nprp; i++) {
            m[i] = (double *) allocMem(npop1*sizeof(double));
        }
        nprp = 0;

        /* admix proportions specified on command line? */
        if (argc - optind == npop) {
            sum = 0;
            for (j = 0; j < npop1; j++) {
                m[nprp][j] = atof(argv[optind+1+j]);
                sum += m[nprp][j];
            }
            if (sum > 1) {
                fprintf(stderr,
                        "command line admix proportions sum to more than 1\n");
                exit(1);
            }
            nprp++;
        }

        /* if necessary, read admix proportions from input file */
        if (pfile) {
            line = 0;
            while (fgets(rec, sizeof(rec), pfp)) {
                line++;
                recp = strtok(rec, " \t\n");
                if (!recp) continue;	/* skip blank line */

                if (sscanf(recp, "%lf", &m[nprp][0]) != 1) {
                    fprintf(stderr,
                        "invalid or missing admix proportion on line %d of %s\n",
                            line, pfile);
                    exit(1);
                }

                sum = m[nprp][0];
                j = 1;
                while (recp = strtok(NULL, " \t\n")) {
                    if (sscanf(recp, "%lf", &m[nprp][j]) != 1) {
                        fprintf(stderr,
                                "invalid admix proportion on line %d of %s\n",
                                line, pfile);
                        exit(1);
                    }
                    sum += m[nprp][j];
                    j++;
                }

                if (j != npop1) {
                    fprintf(stderr,
                            "not enough admix proportions on line %d of %s\n",
                            line, pfile);
                    exit(1);
                }

                if (sum > 1) {
                    fprintf(stderr,
                        "admix proportions sum to more than 1 on line %d of %s\n",
                            line, pfile);
                    exit(1);
                }
                nprp++;
            }

            fclose(pfp);
        }
    }

    /* set up storage for a list of the sets of markers to be tested */
    if (!mfile && !useall) {
        ntst = nmrk;
        mtst = (int **) allocMem(ntst*sizeof(int *));
        nmtst = (int *) allocMem(ntst*sizeof(int));
        for (i = 0; i < ntst; i++) {
            mtst[i] = (int *) allocMem(sizeof(int));
            mtst[i][0] = i;
            nmtst[i] = 1;
        }
    }
    else {
        ntst = 0;
        if (useall) ntst++;

        if (mfile) {
            mfp = fopen(mfile, "r");
            if (!mfp) {
                fprintf(stderr, "cannot open marker list %s\n", mfile);
                exit(1);
            }
            while (fgets(rec, sizeof(rec), mfp)) {
                ntst++;
            }
            rewind(mfp);
        }

        mtst = (int **) allocMem(ntst*sizeof(int *));
        nmtst = (int *) allocMem(ntst*sizeof(int));
        for (i = 0; i < ntst; i++) {
            mtst[i] = (int *) allocMem(nmrk*sizeof(int));
        }

        if (useall) {
            nmtst[0] = nmrk;
            for (i = 0; i < nmrk; i++) {
                mtst[0][i] = i;
            }
        }
    }

    /* if necessary, read marker sets from input file */
    if (mfile) {
        ntst = 0;
        if (useall) ntst++;

        while (fgets(rec, sizeof(rec), mfp)) {
            recp = strtok(rec, " \t\n");
            if (!recp) continue;	/* skip blank line */

            mtst[ntst][0] = -1;
            for (j = 0; j < nmrk; j++) {
                if (!strcmp(recp, mrk[j].name)) {
                    mtst[ntst][0] = j;
                    break;
                }
            }
            if (mtst[ntst][0] == -1) {
                fprintf(stderr,
                        "marker %s not found in frequencies file\n", argv[optind]);
                exit(1);
            }
            nmtst[ntst] = 1;

            while (recp = strtok(NULL, " \t\n")) {
                mtst[ntst][nmtst[ntst]] = -1;
                for (j = 0; j < nmrk; j++) {
                    if (!strcmp(recp, mrk[j].name)) {
                        mtst[ntst][nmtst[ntst]] = j;
                        break;
                    }
                }
                if (mtst[ntst][nmtst[ntst]] == -1) {
                    fprintf(stderr,
                            "marker %s not found in frequencies file\n", recp);
                    exit(1);
                }
                for (j = 0; j < nmtst[ntst]; j++) {
                    if (mtst[ntst][nmtst[ntst]] == mtst[ntst][j]) {
                        fprintf(stderr,
    "marker %s occurs more than once in marker set %d, but was used only once\n",
                                recp, ntst + 1);
                        break;
                    }
                }
                if (j == nmtst[ntst]) {
                    nmtst[ntst]++;
                }
            }
            ntst++;
        }

        fclose(mfp);
    }

    delta = (double *) allocMem(npop1*sizeof(double));
    info = (double *) allocMem(npop1*npop1*sizeof(double));
    ipvt = (int *) allocMem(npop1*sizeof(int));
    work = (double *) allocMem(npop1*sizeof(double));

    if (ofile) {
        ofp = fopen(ofile, "w");
        if (!ofp) {
            fprintf(stderr, "cannot open output file %s\n", ofile);
            exit(1);
        }
    }
    else
        ofp = stdout;

/*
 *  For each set of markers to be tested, compute the informativeness of that
 *  set of markers for admixture, either over a grid of admixture proportions
 *  or at each of a specified set of proportions.
 */
    compute_inform (grid, nprp, m, nmtst[0], mtst[0], ofp);
    for (i = 1; i < ntst; i++) {
        fprintf(ofp, "\n");
        compute_inform (grid, nprp, m, nmtst[i], mtst[i], ofp);
    }

    fclose(ofp);
}

void
compute_inform (int grid, int nprp, double *m[npop-1], int ntst, int *mtst, FILE *fp)
{
    int i, ii, j, k, l;
    int npop1, iprp;
    int nose, job = 11;
    int done, ok;
    double denom, sum, det[2];

    fprintf(fp, "Marker Set:");
    if (ntst == nmrk)
        fprintf(fp, " all markers\n");
    else {
        for (i = 0; i < ntst; i++)
            fprintf(fp, " %s", mrk[mtst[i]].name);
        fprintf(fp, "\n");
    }
    for (i = 1; i <= npop; i++)
        fprintf(fp, "  M%d  ", i);
    fprintf(fp, " INFORM\n");
    npop1 = npop - 1;

    if (grid) nprp = 1;
    for (iprp = 0; iprp < nprp; iprp++) {
        done = 0;
        while (!done) {
            for (i = 0; i < npop1; i++) {
                for (j = 0; j < npop1; j++) {
                    info[i*npop1+j] = 0;
                }
            }

            ok = 1;
            for (ii = 0; ii < ntst; ii++) {
                i = mtst[ii];
                for (j = 0; j < mrk[i].nall; j++) {
                    sum = 0;
                    for (k = 0; k < npop; k++)
                        sum += mrk[i].afreq[k][j];
                    if (!sum) continue;

                    denom = mrk[i].afreq[npop1][j];
                    for (k = 0; k < npop1; k++) {
                        delta[k] = mrk[i].afreq[k][j] - mrk[i].afreq[npop1][j];
                        denom += m[iprp][k]*delta[k];
                    }
                    if (denom > TINY) {
                        for (k = 0; k < npop1; k++) {
                            for (l = 0; l < npop1; l++) {
                                info[k*npop1+l] += 2*delta[k]*delta[l]/denom;
                            }
                        }
                    }
                    else {
                        ok = 0;
                        break;
                    }
                }
                if (!ok) break;
            }

            if (ok) {
                dgefa_(info, &npop1, &npop1, ipvt, &nose);
                if (!nose) {
                    dgedi_(info, &npop1, &npop1, ipvt, det, work, &job);
                }
            }

            sum = 0;
            for (i = 0; i < npop1; i++) {
                fprintf(fp, "%5.3f ", m[iprp][i]);
                sum += m[iprp][i];
            }
            if (sum > 1) sum = 1;
            fprintf(fp, "%5.3f ", 1 - sum);
            if (!ok)
                fprintf(fp, "********\n");
            else if (!nose)
                fprintf(fp, "%g\n", det[0]*pow(10.,det[1]));
            else
                fprintf(fp, "information matrix can't be factored\n");

            if (grid) {
                for (i = npop1 - 1; i >= 0; i--) {
                    sum = 0;
                    for (j = 0; j < npop1; j++) {
                        if (j != i)
                            sum += m[iprp][j];
                    }
                    m[iprp][i] += gridint;
                    if (m[iprp][i] > 1 + TINY - sum)
                        m[iprp][i] = 0;
                    else
                        break;
                }
                if (i == -1)
                    done = 1;
            }
            else
                done = 1;
        }
    }
}

void
show_usage (char *prog)
{
    printf(
"usage: %s [-aq] [-p plist] [-m mlist] [-g gridint] [-o outfile] locfile [m1 ...]\n\n",
           prog);
    printf("   locfile      population-specific allele frequencies\n");
    printf("   m1 m2 ...    admixture proportions (optional)\n");
    printf("\n   options:\n");
    printf("     -p plist     read sets of admix proportions from file plist\n");
    printf("     -m mlist     read sets of markers from file mlist\n");
    printf("     -a           use all markers simultaneously\n");
    printf("     -g gridint   grid interval\n");
    printf("     -o outfile   output file\n");
    printf("     -q           suppress warning messages\n");
    exit(1);
}

void
read_locus_file (char *locfile, int nowarn)
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
            fprintf(stderr,
                    "%s: blank lines not allowed, line %d\n", locfile, line);
            exit(1);
        }
    }

    if (npop > MXPOP) {
        fprintf(stderr, "%s: too many populations, MXPOP = %d\n", locfile, MXPOP);
        exit(1);
    }

    mrk = (struct Marker *) allocMem(nmrk*sizeof(struct Marker));

    rewind(fp);
    fgets(rec, sizeof(rec), fp);
    mrk[0].name = (char *) allocMem(strlen(rec)); /* rec contains \n */
    sscanf(rec, "%s", mrk[0].name);
    line = 1;
    for (i = 0; i < nmrk; i++) {
        mrk[i].nall = 0;
        for (j = 0; j < npop; j++) freq[j] = 0;

        while (fgets(rec, sizeof(rec), fp) &&
                   sscanf(rec, "%s %s", buf, buf) == 2)
        {
            line++;
            recp = strtok(rec, " \t\n");
            iall = add_name(recp, &mrk[i].alleles, &mrk[i].nall);

            for (j = 0; j < npop; j++) {
                if (!(recp = strtok(NULL, " \t\n")))
                {
                    fprintf(stderr, "%s: missing allele frequency, line %d\n",
                            locfile, line);
                    exit(1);
                }
                if (sscanf(recp, "%lf", &mrk[i].afreq[j][iall]) != 1 ||
                         mrk[i].afreq[j][iall] < 0 || mrk[i].afreq[j][iall] > 1)
                {
                    fprintf(stderr, "%s: invalid allele frequency, line %d\n",
                            locfile, line);
                    exit(1);
                }
                freq[j] += mrk[i].afreq[j][iall];
            }
        }

        for (j = 0; j < npop; j++) {
            if (!nowarn && (freq[j] < 1 - TINY || freq[j] > 1 + TINY)) {
                fprintf(stderr,
        "Warning: allele frequencies sum to %f for marker %d, population %d\n",
                        freq[j], i+1, j+1);
                fprintf(stderr, "   Frequencies being adjusted to sum to 1.\n");
            }
            for (k = 0; k < mrk[i].nall; k++) {
                mrk[i].afreq[j][k] /= freq[j];
            }
        }

        if (i < nmrk - 1) {
            mrk[i+1].name = (char *) allocMem(strlen(rec)); /* rec contains \n */
            sscanf(rec, "%s", mrk[i+1].name);
        }
    }

    fclose(fp);
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
