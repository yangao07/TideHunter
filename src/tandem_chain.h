#ifndef TANDEM_CHAIN_H
#define TANDEM_CHAIN_H

#include "tide_hunter.h"
#include "tandem_hit.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int i:32, j:32;
} cell_t;

typedef struct {
    int32_t x,y,z;
} triple_t;

typedef struct {
    int i:32, j:32;
    int score;
} dp_score_t;

typedef struct {
    int from_i, from_j; // i : index of end postion, j : index of hit of end
    int start, end, mem_l; int score;
    int8_t is_tracked;
} dp_t;

typedef struct {
    cell_t *cell;
    int len, score;
    int est_ch_i, est_period, est_start;
    int max_period, min_period;
} chain_t;

typedef struct {
    int start, end;
} rep_reg_t;

//void radix_sort_hash(hash_t *beg, hash_t *end);

int hash_partition(char *seq, int seq_len, tandem_seq_t *tandem_seq, mini_tandem_para *mtp);
int mini_tandem_core(kseq_t *read_seq, tandem_seq_t *tandem_seq, mini_tandem_para *mtp);
int tandem_chain(int seq_len, hash_t *hash_hit, int hash_hit_n, mini_tandem_para *mtp, dp_t ***_dp, int *tot_N,chain_t **post_chain, int *post_ch_m);

#ifdef __cplusplus
}
#endif

#endif
