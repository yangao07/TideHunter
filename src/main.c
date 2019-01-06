#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <time.h>
#include <math.h>
#include "utils.h"
#include "mini_tandem.h"

const char PROG[20] = "miniTandem";
const char VERSION[20] = "1.0.0";
//const char DATE[20] = "2018-12-01";
const char CONTACT[30] = "yangaoucla@gmail.com";


const struct option mini_tandem_opt [] = {
    { "kmer-length", 1, NULL, 'k' },
    { "window-size", 1, NULL, 'w' },
    { "step-size", 1, NULL, 's' },
    { "minimal-m", 1, NULL, 'm' },
    { "HPC-kmer", 0, NULL, 'H' },
    { "max-diverg", 1, NULL, 'e' },
    { "min-period", 1, NULL, 'p' },
    { "max-period", 1, NULL, 'P' },

    { "rep-range", 1, NULL, 'r' },

    { "splint-seq", 1, NULL, 'S' },
    { "detail-out", 1, NULL, 'd' },

    { "thread", 1, NULL, 't' },
    { 0, 0, 0, 0}
};

static int usage(void)
{
    err_printf("\n");
	err_printf("%s: Tandem repeat detection and consensus calling from noisy concatemeric long-read\n\n", PROG);

    time_t t; time(&t);
    err_printf("Version: %s Build date: %s", VERSION, ctime(&t));
    err_printf("Contact: %s\n\n", CONTACT);

    err_printf("Usage:   %s [options] in.fa/fq > cons_out.fastq\n\n", PROG);

	err_printf("Options: \n");
    err_printf("         -t --thread      [INT]    number of threads to use. [%d]\n", THREAD_N);
    err_printf("         -k --kmer-length [INT]    k-mer length (no larger than 16). [%d]\n", KMER_SIZE); // TODO largest kmer len
    err_printf("         -w --window-size [INT]    window size. [%d]\n", KMER_WSIZE);
    err_printf("         -s --step-size   [INT]    step size. [%d]\n", KMER_SSIZE);
    err_printf("         -m --minimal-m   [INT]    number of minimal k-mer to keep in each window. [%d]\n", KMER_MINM);
    err_printf("         -e --max-diverg  [INT]    maximum allowed divergence rate between two consecutive repeats. [%.2f]\n", MAX_DIV);
    err_printf("         -H --HPC-kmer             use homopolymer-compressed k-mer. [False]\n");
    err_printf("         -p --min-period  [INT]    minimum period size of tandem repeat. (>=%d) [%d]\n", 2, MIN_PERIOD);
    err_printf("         -P --max-period  [INT]    maximum period size of tandem repeat. (<=%d) [%d]\n", MAX_PERIOD, MAX_PERIOD);

//  err_printf("         -r --rep-range   [INT]    maximum range to find tandem repeat. [%d]\n", REP_RANGE); 
//  err_printf("                                   (-1: no limit, tandem repeat can span the whole sequence)\n");

    err_printf("         -S --splint-seq  [STR]    splint sequence in FASTA/FASTQ format. [NULL]\n");
    err_printf("         -d --detail-out  [STR]    detailed information of each consensus. [NULL]\n");
    err_printf("                                   (start, end, score, etc.)\n");

	err_printf("\n");
	return 1;
}

// TODO add score para 
int main(int argc, char *argv[])
{
    mini_tandem_para *mtp = mini_tandem_init_para();
    int c;
    while ((c = getopt_long(argc, argv, "k:w:m:Hs:r:e:p:P:S:d:t:",mini_tandem_opt, NULL)) >= 0) {
        switch(c)
        {
            case 'k': mtp->k = atoi(optarg); break;
            case 'w': mtp->w = atoi(optarg); break;
            case 'm': mtp->m = atoi(optarg); break;
            case 's': mtp->s = atoi(optarg); break;
            case 'H': mtp->hpc = 1; break;
            case 'e': mtp->max_div = atof(optarg); break;
            case 'p': mtp->min_p = atoi(optarg);
                      if (mtp->min_p < MIN_PERIOD) {
                          err_printf("Error: -p --min-period(%d) needs to be >= %d.\n", mtp->min_p, MIN_PERIOD); 
                          goto End;
                      }
                      break;
            case 'P': mtp->max_p = atoi(optarg); 
                      if (mtp->max_p > MAX_PERIOD) {
                          err_printf("Error: -P --max-period(%d) needs to be <= %d.\n", mtp->max_p, MAX_PERIOD); 
                          goto End;
                      }
                      break;
          //case 'r': mtp->max_range = atoi(optarg); break;
            case 'd': mtp->detail_fp = xopen(optarg, "w"); break;
            case 'S': mtp->splint_fn = strdup(optarg); break;
            case 't': mtp->n_thread = atoi(optarg); break;
            default:
                      err_printf("Error: unknown option: -%c %s.\n", c, optarg);
                      goto End;
        }
    }
	if (argc < 2) return usage();

    mtp->div_exp=exp(mtp->k * mtp->max_div);
    mini_tandem(argv[optind], mtp);
End:
    mini_tandem_free_para(mtp);
    double sys_t, usr_t; usr_sys_cputime(&usr_t, &sys_t); err_func_printf(__func__, "User: %.3f sec; Sys: %.3f sec.\n", usr_t, sys_t);
    return 0;
}
