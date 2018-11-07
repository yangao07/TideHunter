#ifndef _MINI_TANDEM_H_
#define _MINI_TANDEM_H_

typedef struct {
    int k, w, m; // k-mer length, window size, selected minimum m hash values
                 // keep all k-mer when w == m
    double sigma; int bucket_T; // hit bucket range: P * sigma; filter out bucket with size < bucket_T
    int max_range; // max range to find tandem repeat, -1 for no limit
    FILE *detail_fp; //char detail_out[1024];
    int n_thread;
} mini_tandem_para;

#define THREAD_N 4
#define CHUNK_READ_N 10000

#define KMER_SIZE 8
#define KMER_WSIZE 1
#define KMER_MINM 1 // 1: minimizer
#define HIT_BKT_SIG 0.005
#define HIT_BKT_SIZE 2
#define REP_RANGE -1


int mini_tandem(const char *read_fn, mini_tandem_para *mtp);

mini_tandem_para *mini_tandem_init_para(void);

void mini_tandem_free_para(mini_tandem_para *mtp);

#endif
