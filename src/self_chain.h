#ifndef SELF_CHAIN_H
#define SELF_CHAIN_H

#ifdef __cplusplus
extern "C" {
#endif

// TODO use int32_t 
typedef int64_t hash_t;
//typedef struct {int64_t x:32,y:32;} hash_t; // TODO
#define _8mask 0xff
#define _16mask 0xffff
#define _32mask 0xffffffff

typedef struct {
    int min_p, max_p;
    int start_i, end_i;
    int tot_n;
} hit_bucket_t;

typedef struct {
    int i:32, j:32;
} cell_t;

typedef struct {
    int32_t x,y,z;
} triple_t;

typedef struct {
    int i:32, j:32;
    double score;
} dp_score_t;

typedef struct {
    int from_i, from_j; // i : index of end postion, j : index of hit of end
    int start, end;
    double period, score;
    int8_t is_tracked;
} self_dp_t;

typedef struct {
    cell_t *chain;
    int len;
} chain_t;

typedef struct {
    int start, end;
} rep_reg_t;

//void radix_sort_hash(hash_t *beg, hash_t *end);

int hash_partition(char *seq, int seq_len, mini_tandem_para *mtp);

#ifdef __cplusplus
}
#endif

#endif
