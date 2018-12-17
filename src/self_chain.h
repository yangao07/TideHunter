#ifndef SELF_CHAIN_H
#define SELF_CHAIN_H

#ifdef __cplusplus
extern "C" {
#endif

typedef int64_t hash_t;
#define _8mask 0xff
#define _16mask 0xffff
#define _32mask 0xffffffff

#define MIN_SIM_PERIOD_RAT 0.98
#define MAX_SIM_PERIOD_RAT 1.02

typedef struct {
    int min_p, max_p;
    int start_i, end_i;
    int tot_n;
} hit_bucket_t;

typedef struct {
    int i:32, j:32;
} cell_t;

// TODO total_hits_of_kmer: used as weight of kmer
// hash_hit1: period:32 | end:32
#define _set_hash_hit1(start, end) ((((end)-(start)) << 32) | ((end) & _32mask))
#define _get_hash_hit1_end(hash_hit, i) ((hash_hit)[i] & _32mask)
#define _get_hash_hit1_period(hash_hit, i) ((hash_hit)[i] >> 32)
#define _get_hash_hit1_start(hash_hit, i) (_get_hash_hit1_end((hash_hit), i) - _get_hash_hit1_period((hash_hit), i))

//TODO max allowed period size: pow(2,16)
//TODO deal with K >= pow(2,16) : use 128_t end:64 | period:32 | K:32
// hash_hit2: end:32 | period:16 | K:16
// K: MEM hit length
#define _set_hash_hit2(end, period, k) (((end) << 32) | ((period) << 16) | k)
#define _get_hash_hit2_end(hash_hit, i) ((hash_hit)[i] >> 32)
#define _get_hash_hit2_period(hash_hit, i) (((hash_hit)[i] & _32mask) >> 16)
#define _get_hash_hit2_start(hash_hit, i) (_get_hash_hit2_end((hash_hit), i) - _get_hash_hit2_period((hash_hit), i))
#define _get_hash_hit2_meml(hash_hit, i) ((hash_hit)[i] & _16mask)


#define _get_hash_hit_end(hash_hit, i) ((hash_hit)[i] >> 32)
#define _get_hash_hit_period(hash_hit, i) (((hash_hit)[i] & _32mask) >> 16)
#define _get_hash_hit_seed_id(hash_hit, i) ((hash_hit)[i] & _16mask)
#define _get_hash_hit_start(hash_hit, i) (_get_hash_hit_end((hash_hit), i) - _get_hash_hit_period((hash_hit), i))

#define _set_hash_hit(start, end, hi) (((end) << 32) | (((end)-(start)) << 16) | hi)

typedef struct {
    int32_t x,y,z;
} triple_t;

typedef struct {
    int i:32, j:32;
    int score;
} dp_score_t;

typedef struct {
    int from_i, from_j; // i : index of end postion, j : index of hit of end
    int start, end, mem_l;
    double period; int score;
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
