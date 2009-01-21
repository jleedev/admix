/*
 * count.c
 *
 * Written by Thomas Dyer
 * Copyright (c) 2003, 2004, 2005
 *
 *
 * This program reads the genotype data for a set of individuals and
 * generates a tally of the marker alleles and genotypes. Allele and
 * genotype counts are broken down both by population and affection
 * status. The counts are stored in files, one for allele counts and
 * one for genotype counts, which are used by various other programs.
 * There is a limit of 2 populations (MXPOP) and 40 alleles per marker
 * (MXALL). If necessary, these limits can be increased in the source
 * code and the program recompiled.
 *
 *
 * Usage:  count [-M missval] locfile pedfile afrqout gfrqout
 *
 *         locfile      locus file (marker info)
 *         pedfile      pedigree file (marker genotypes)
 *         afrqout      allele frequencies output file
 *         gfrqout      genotype frequencies output file
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
 *   The locus file is simply a list of marker names, one name per
 *   line, in the order that the corresponding marker genotypes occur
 *   in the pedigree file. If a locus file containing allele frequencies
 *   (as described in the documentation for program admix) is used,
 *   the lines containing allele frequencies are simply ignored.
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
 *   There are two blank-delimited output files - one for allele counts,
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
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define MXPOP	4	/* max # populations		*/
#define MXALL	40	/* max # alleles per marker	*/

#define MISSVAL	"*"	/* missing value		*/

#define ALLBLK	256	/* # array elements malloc'd at a time		*/


struct Pop {
    int acntu[MXALL];		/* allele counts for unaffecteds	*/
    int acnta[MXALL];		/* allele counts for affecteds		*/
    int gcntu[MXALL][MXALL];	/* genotype counts for unaffecteds	*/
    int gcnta[MXALL][MXALL];	/* genotype counts for affecteds	*/
} ;

struct Marker {
    char *name;			/* marker name				*/
    int nall;			/* number of alleles			*/
    char **alleles;		/* allele names				*/
    int *all_sort;		/* sorted order of allele names		*/
    struct Pop pop[MXPOP];	/* by-population allele/genotype counts	*/
} ;

int nmrk;		/* number of markers		*/
struct Marker *mrk;


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
void read_locus_file (char *);
void read_pedigree_file (char *, char *);
void do_allele_freqs (struct Marker *, int, int, FILE *);
void do_genotype_freqs (struct Marker *, int, int, FILE *);
int get_ndx (char *, char **, int);
int add_name (char *, char ***, int *);
void sort_names(char **, int, int *);
void *allocMem (size_t);


main (int argc, char **argv)
{
    int i, j;
    char *recp, rec[10000];
    char buf[10];
    FILE *fpa, *fpg;

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

    read_locus_file(argv[optind]);

    read_pedigree_file(argv[optind+1], missval);

    fpa = fopen(argv[optind+2], "w");
    if (!fpa) {
        fprintf(stderr, "cannot open %s\n", argv[optind+2]);
        exit(1);
    }

    fpg = fopen(argv[optind+3], "w");
    if (!fpg) {
        fprintf(stderr, "cannot open %s\n", argv[optind+3]);
        exit(1);
    }


/*
 *  For each locus, sum up allele counts, then divide to get frequencies.
 *  Force the frequencies into a fixed-width format, sum the first N - 1
 *  formatted frequencies where N = #alleles, and subtract this sum from
 *  1 to get the frequency of the Nth allele. This makes the frequencies
 *  sum to 1.
 */

    for (i = 0; i < nmrk; i++) {

        /* unaffected */
	do_allele_freqs(&mrk[i], 1, -1, fpa);

        /* unaffected, by population */
        for (j = 0; j < npop; j++) {
	    do_allele_freqs(&mrk[i], 1, j, fpa);
        }

        /* affected */
	do_allele_freqs(&mrk[i], 2, -1, fpa);

        /* affected, by population */
        for (j = 0; j < npop; j++) {
	    do_allele_freqs(&mrk[i], 2, j, fpa);
        }

        /* by population */
        for (j = 0; j < npop; j++) {
	    do_allele_freqs(&mrk[i], 0, j, fpa);
        }

        /* everybody */
	do_allele_freqs(&mrk[i], 0, -1, fpa);
    }


/*
 *  Genotype frequencies.
 */

    for (i = 0; i < nmrk; i++) {

        /* unaffected */
	do_genotype_freqs(&mrk[i], 1, -1, fpg);

        /* unaffected, by population */
        for (j = 0; j < npop; j++) {
	    do_genotype_freqs(&mrk[i], 1, j, fpg);
        }

        /* affected */
	do_genotype_freqs(&mrk[i], 2, -1, fpg);

        /* affected, by population */
        for (j = 0; j < npop; j++) {
	    do_genotype_freqs(&mrk[i], 2, j, fpg);
        }

        /* by population */
        for (j = 0; j < npop; j++) {
	    do_genotype_freqs(&mrk[i], 0, j, fpg);
        }

        /* everybody */
	do_genotype_freqs(&mrk[i], 0, -1, fpg);
    }
}

void
show_usage (char *prog)
{
    printf("usage: %s [-M missval] locfile pedfile afrqout gfrqout\n\n", prog);
    printf("   locfile      locus file\n");
    printf("   pedfile      pedigree file\n");
    printf("   afrqout      allele frequencies output file\n");
    printf("   gfrqout      genotype frequencies output file\n");
    printf("\n   options:\n");
    printf("     -M missval   missing allele value in quotes\n");
    exit(1);
}

void read_locus_file (char *locfile)
{
    char *recp, rec[10000];
    char buf[10000];
    int i;
    FILE *fp;

    fp = fopen(locfile, "r");
    if (!fp) {
        fprintf(stderr, "cannot open locus file %s\n", locfile);
        exit(1);
    }

    nmrk = 0;
    while (fgets(rec, sizeof(rec), fp)) {
        if (sscanf(rec, "%s %s", buf, buf) == 1)
            nmrk++;
    }

    mrk = (struct Marker *) allocMem(nmrk*sizeof(struct Marker));

    rewind(fp);
    fgets(rec, sizeof(rec), fp);
    for (i = 0; i < nmrk; i++) {
        mrk[i].name = (char *) allocMem(strlen(rec));	/* rec contains \n */
        sscanf(rec, "%s", mrk[i].name);
        mrk[i].nall = 0;
        while (fgets(rec, sizeof(rec), fp) &&
               sscanf(rec, "%s %s", buf, buf) == 2) ;
    }

    fclose(fp);
}

void read_pedigree_file (char *pedfile, char *missval)
{
    char *recp, rec[10000];
    int line = 0;
    int famid, id, sex, aff, pop;
    int i, j, k, l;
    int *all[2];
    FILE *fp;

    fp = fopen(pedfile, "r");
    if (!fp) {
        fprintf(stderr, "cannot open pedigree file %s\n", pedfile);
        exit(1);
    }

    all[0] = (int *) allocMem(nmrk*sizeof(int));
    all[1] = (int *) allocMem(nmrk*sizeof(int));

    while (fgets(rec, sizeof(rec), fp)) {
        line++;
/*
    This way of reading the data file (using strtok) assumes that the
    fields are ALWAYS separated by blanks or tabs.
*/
        if (!(recp = strtok(rec, " \t\n")))
        {
            fprintf(stderr, "%s: missing family ID, line %d\n", pedfile, line);
            exit(1);
        }
        if ((famid = get_ndx(recp, famids, nfamid)) == -1)
            famid = add_name(recp, &famids, &nfamid);

        if (!(recp = strtok(NULL, " \t\n")))
        {
            fprintf(stderr, "%s: missing ID, line %d\n", pedfile, line);
            exit(1);
        }
        if ((id = get_ndx(recp, ids, nid)) == -1)
            id = add_name(recp, &ids, &nid);

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

        if ((pop = get_ndx(recp, pops, npop)) == -1) {
            if (npop == MXPOP) {
                fprintf(stderr, "too many populations, MXPOP = %d\n", MXPOP);
                exit(1);
            }

            pop = add_name(recp, &pops, &npop);
            for (j = 0; j < nmrk; j++) {
                for (k = 0; k < MXALL; k++) {
                    mrk[j].pop[pop].acntu[k] = 0;
                    mrk[j].pop[pop].acnta[k] = 0;
                    for (l = 0; l < MXALL; l++) {
                        mrk[j].pop[pop].gcnta[k][l] = 0;
                        mrk[j].pop[pop].gcnta[l][k] = 0;
                        mrk[j].pop[pop].gcntu[k][l] = 0;
                        mrk[j].pop[pop].gcntu[l][k] = 0;
                    }
                }
            }
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

        /* read in the marker alleles */
        for (i = 0; i < nmrk; i++) {
            if (!(recp = strtok(NULL, " \t\n")))
            {
                fprintf(stderr, "%s: missing allele, marker %s, line %d\n",
                        pedfile, mrk[i].name, line);
                exit(1);
            }
            if (!strcmp(recp, missval))
            {
                all[0][i] = -1;
            }
            else if ((all[0][i] =
                 get_ndx(recp, mrk[i].alleles, mrk[i].nall)) == -1)
            {
                all[0][i] = add_name(recp, &mrk[i].alleles, &mrk[i].nall);
            }

            if (!(recp = strtok(NULL, " \t\n")))
            {
                fprintf(stderr, "%s: missing allele, marker %s, line %d\n",
                        pedfile, mrk[i].name, line);
                exit(1);
            }
            if (!strcmp(recp, missval))
            {
                all[1][i] = -1;
            }
            else if ((all[1][i] =
                 get_ndx(recp, mrk[i].alleles, mrk[i].nall)) == -1)
            {
                all[1][i] = add_name(recp, &mrk[i].alleles, &mrk[i].nall);
            }

            if (mrk[i].nall > MXALL) {
                fprintf(stderr, "marker %s has too many alleles, MXALL = %d\n",
                        mrk[i].name, MXALL);
                exit(1);
            }
        }

        /* count marker alleles and genotypes */
        for (i = 0; i < nmrk; i++) {
            if (all[0][i] != -1) {
                if (aff == 2)
                    mrk[i].pop[pop].acnta[all[0][i]]++;
                else
                    mrk[i].pop[pop].acntu[all[0][i]]++;
            }
            if (all[1][i] != -1) {
                if (aff == 2)
                    mrk[i].pop[pop].acnta[all[1][i]]++;
                else
                    mrk[i].pop[pop].acntu[all[1][i]]++;
            }
            if (all[0][i] != -1 && all[1][i] != -1) {
                if (aff == 2)
                    mrk[i].pop[pop].gcnta[all[0][i]][all[1][i]]++;
                else
                    mrk[i].pop[pop].gcntu[all[0][i]][all[1][i]]++;
            }
        }
    }

    fclose(fp);

    for (i = 0; i < nmrk; i++) {
        mrk[i].all_sort = (int *) allocMem(mrk[i].nall*sizeof(int));
        sort_names(mrk[i].alleles, mrk[i].nall, mrk[i].all_sort);
    }
}

void do_allele_freqs (struct Marker *m, int do_aff, int do_pop, FILE *fp)
{
    int j, k, kk, tsum, fnd, nfrq;
    double sum, num, asum;
    char affect, *popid, buf[10];

    switch (do_aff) {
    case 0:
	affect = '-';
        break;
    case 1:
	affect = 'U';
        break;
    case 2:
	affect = 'A';
        break;
    }

    if (do_pop == -1)
        popid = "-";
    else
        popid = pops[do_pop];

    nfrq = 0;
    sum = 0;
    for (k = 0; k < m->nall; k++) {
        fnd = 0;
        for (j = 0; j < npop; j++) {
            if (do_pop >= 0 && j != do_pop) continue;
            if (do_aff != 2 && m->pop[j].acntu[k]) {
                fnd = 1;
                sum += m->pop[j].acntu[k];
            }
            if (do_aff != 1 && m->pop[j].acnta[k]) {
                fnd = 1;
                sum += m->pop[j].acnta[k];
            }
        }
        if (fnd) nfrq++;
    }

    asum = 0;
    for (k = 0; nfrq && k < m->nall; k++) {
        kk = m->all_sort[k];
        tsum = 0;
        for (j = 0; j < npop; j++) {
            if (do_pop >= 0 && j != do_pop) continue;
            if (do_aff != 2)
                tsum += m->pop[j].acntu[kk];
            if (do_aff != 1)
                tsum += m->pop[j].acnta[kk];
        }
        if (tsum) {
            nfrq--;
            if (nfrq) {
                sprintf(buf, "%7.5f", tsum/sum);
                sscanf(buf, "%lf", &num);
                asum += num;
                fprintf(fp, "%s %c %s %s %d %7.5f\n", m->name, affect, popid,
                        m->alleles[kk], tsum, num);
            }
            else
                fprintf(fp, "%s %c %s %s %d %7.5f\n", m->name, affect, popid,
                        m->alleles[kk], tsum, 1. - asum);
        }
    }
}

void do_genotype_freqs (struct Marker *m, int do_aff, int do_pop, FILE *fp)
{
    int j, k, kk, l, ll, tsum, fnd, nfrq;
    double sum, num, asum;
    char affect, *popid, buf[10];

    switch (do_aff) {
    case 0:
	affect = '-';
        break;
    case 1:
	affect = 'U';
        break;
    case 2:
	affect = 'A';
        break;
    }

    if (do_pop == -1)
        popid = "-";
    else
        popid = pops[do_pop];

    nfrq = 0;
    sum = 0;
    for (k = 0; k < m->nall; k++) {
        for (l = 0; l < k; l++) {
            fnd = 0;
            for (j = 0; j < npop; j++) {
                if (do_pop >= 0 && j != do_pop) continue;
                if (do_aff != 2 &&
                    (m->pop[j].gcntu[k][l] || m->pop[j].gcntu[l][k]))
                {
                    fnd = 1;
                    sum += m->pop[j].gcntu[k][l] + m->pop[j].gcntu[l][k];
                }
                if (do_aff != 1 &&
                    (m->pop[j].gcnta[k][l] || m->pop[j].gcnta[l][k]))
                {
                    fnd = 1;
                    sum += m->pop[j].gcnta[k][l] + m->pop[j].gcnta[l][k];
                }
            }
            if (fnd) nfrq++;
        }

        fnd = 0;
        for (j = 0; j < npop; j++) {
            if (do_pop >= 0 && j != do_pop) continue;
            if (do_aff != 2 && m->pop[j].gcntu[k][k]) {
                fnd = 1;
                sum += m->pop[j].gcntu[k][k];
            }
            if (do_aff != 1 && m->pop[j].gcnta[k][k]) {
                fnd = 1;
                sum += m->pop[j].gcnta[k][k];
            }
        }
        if (fnd) nfrq++;
    }

    asum = 0;
    for (k = 0; nfrq && k < m->nall; k++) {
        kk = m->all_sort[k];
        for (l = 0; nfrq && l < k; l++) {
            ll = m->all_sort[l];
            tsum = 0;
            for (j = 0; j < npop; j++) {
                if (do_pop >= 0 && j != do_pop) continue;
                if (do_aff != 2)
                    tsum += m->pop[j].gcntu[kk][ll] + m->pop[j].gcntu[ll][kk];
                if (do_aff != 1)
                    tsum += m->pop[j].gcnta[kk][ll] + m->pop[j].gcnta[ll][kk];
            }
            if (tsum) {
                nfrq--;
                if (nfrq) {
                    sprintf(buf, "%7.5f", tsum/sum);
                    sscanf(buf, "%lf", &num);
                    asum += num;
                    fprintf(fp, "%s %c %s %s %s %d %7.5f\n", m->name,
                            affect, popid, m->alleles[ll], m->alleles[kk],
                            tsum, num);
                }
                else
                    fprintf(fp, "%s %c %s %s %s %d %7.5f\n", m->name,
                            affect, popid, m->alleles[ll], m->alleles[kk],
                            tsum, 1. - asum);
            }
        }

        tsum = 0;
        for (j = 0; j < npop; j++) {
            if (do_pop >= 0 && j != do_pop) continue;
            if (do_aff != 2)
                tsum += m->pop[j].gcntu[kk][kk];
            if (do_aff != 1)
                tsum += m->pop[j].gcnta[kk][kk];
        }
        if (tsum) {
            nfrq--;
            if (nfrq) {
                sprintf(buf, "%7.5f", tsum/sum);
                sscanf(buf, "%lf", &num);
                asum += num;
                fprintf(fp, "%s %c %s %s %s %d %7.5f\n", m->name, affect,
                        popid, m->alleles[kk], m->alleles[kk], tsum, num);
            }
            else
                fprintf(fp, "%s %c %s %s %s %d %7.5f\n", m->name, affect,
                        popid, m->alleles[kk], m->alleles[kk], tsum,
                        1. - asum);
        }
    }
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
