#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include "utils.h"
#include "mini_tandem.h"

const char PROG[20] = "NanoCirc";
const char VERSION[20] = "1.0.0";
const char DATE[20] = "2017-12-15";
const char CONTACT[30] = "yangaoucla@gmail.com";


const struct option mini_tandem_opt [] = {
    { "kmer-length", 1, NULL, 'k' },

    { "thread", 1, NULL, 't' },
    { 0, 0, 0, 0}
};

static int usage(void)
{
    err_printf("\n");
	err_printf("Program: %s (Nanopore Circular RNA consensus caller)\n", PROG);
    err_printf("Version: %s, Date: %s\n", VERSION, DATE);
    err_printf("Contact: %s\n", CONTACT);
    err_printf("Usage:   %s [options] in.fa/fq > cons_out.fa\n\n", PROG);
	err_printf("Options: \n");
    err_printf("         -k --kmer-length [INT]    k-mer length (no larger than XXX). [%d]\n", KMER_SIZE); // TODO largest kmer len

    err_printf("         -t --thread      [INT]    number of threads used.. [1]\n");
	err_printf("\n");
	return 1;
}

mini_tandem_para *mini_tandem_init_para(void) {
    mini_tandem_para *ncp = (mini_tandem_para*)_err_malloc(sizeof(mini_tandem_para));
    ncp->k = KMER_SIZE;
    ncp->n_thread = 1;
    return ncp;
}


/*
int main(int argc, char *argv[])
{
    mini_tandem_para *ncp = mini_tandem_init_para();
    int c;
    while ((c = getopt_long(argc, argv, "k:t:",mini_tandem_opt, NULL)) >= 0) {
        switch(c)
        {
            case 'k': ncp->k = atoi(optarg); break;
            case 't': ncp->n_thread = atoi(optarg); break;
            default:
                      err_printf("Error: unknown option: -%c %s.\n", c, optarg);
        }
    }
	if (argc < 2) return usage();

    mini_tandem_core(argv[1], ncp);

    free(ncp);
    return 0;
}*/
