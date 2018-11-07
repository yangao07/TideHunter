#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include "utils.h"
#include "mini_tandem.h"

const char PROG[20] = "miniTandem";
const char VERSION[20] = "1.0.0";
const char DATE[20] = "2018-12-01";
const char CONTACT[30] = "yangaoucla@gmail.com";


const struct option mini_tandem_opt [] = {
    { "kmer-length", 1, NULL, 'k' },
    { "window-size", 1, NULL, 'w' },
    { "minimal-m", 1, NULL, 'm' },
    { "rep-range", 1, NULL, 'r' },
    { "bucket-sig", 1, NULL, 'b' },
    { "bucket-thd", 1, NULL, 'u' },

    { "detail-out", 1, NULL, 'd' },

    { "thread", 1, NULL, 't' },
    { 0, 0, 0, 0}
};

static int usage(void)
{
    err_printf("\n");
	err_printf("Program: %s : Consensus calling from noisy concatemeric long-read.\n", PROG);
    err_printf("Version: %s. Date: %s\n", VERSION, DATE);
    err_printf("Contact: %s\n", CONTACT);
    err_printf("Usage:   %s [options] in.fa/fq > cons_out.fastq\n\n", PROG);

	err_printf("Options: \n");
    err_printf("         -t --thread      [INT]    number of threads to use. [%d]\n", THREAD_N);
    err_printf("         -k --kmer-length [INT]    k-mer length (no larger than 16). [%d]\n", KMER_SIZE); // TODO largest kmer len
    err_printf("         -w --window-size [INT]    window size. [%d]\n", KMER_WSIZE);
    err_printf("         -m --minimal-m   [INT]    number of minimal k-mer to keep in each window. [%d]\n", KMER_MINM);
    err_printf("         -b --bucket-sig  [DOU]    variance of each hit-bucket. [%.2f]\n", HIT_BKT_SIG);
    err_printf("         -u --bucket-thd  [INT]    minimum bucket size to keep. [%d]\n", HIT_BKT_SIZE);
    err_printf("         -r --rep-range   [INT]    maximum range to find tandem repeat. [%d]\n", REP_RANGE); 
    err_printf("                                   (-1 means no limit, tandem repeat can span the whole sequence)\n");

    err_printf("         -d --detail-out  [STR]    detailed information of each consensus. [NULL]\n");
    err_printf("                                   (start, end, score, etc.)\n");

	err_printf("\n");
	return 1;
}

// TODO 
// input splint/primer sequence
int main(int argc, char *argv[])
{
    mini_tandem_para *mtp = mini_tandem_init_para();
    int c;
    while ((c = getopt_long(argc, argv, "k:w:m:b:u:r:d:t:",mini_tandem_opt, NULL)) >= 0) {
        switch(c)
        {
            case 'k': mtp->k = atoi(optarg); break;
            case 'w': mtp->w = atoi(optarg); break;
            case 'm': mtp->m = atoi(optarg); break;
            case 'r': mtp->max_range = atoi(optarg); break;
            case 'b': mtp->sigma = atof(optarg); break;
            case 'u': mtp->bucket_T = atoi(optarg); break;
            case 'd': mtp->detail_fp = xopen(optarg, "w"); break;
            case 't': mtp->n_thread = atoi(optarg); break;
            default:
                      err_printf("Error: unknown option: -%c %s.\n", c, optarg);
        }
    }
	if (argc < 2) return usage();

    mini_tandem(argv[optind], mtp);
    mini_tandem_free_para(mtp);
    return 0;
}
