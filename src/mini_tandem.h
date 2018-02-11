#ifndef _MINI_TANDEM_H_
#define _MINI_TANDEM_H_

typedef struct {
    int k;
    int n_thread;
} mini_tandem_para;

#define CHUNK_READ_N 10000
#define KMER_SIZE 12

int mini_tandem_core(const char *read_fn, mini_tandem_para *ncp);

#endif
