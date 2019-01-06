#ifndef SELF_CHAIN_H
#define SELF_CHAIN_H

#include "mini_tandem.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef int64_t hash_t;
#define _8mask 0xff
#define _16mask 0xffff
#define _32mask 0xffffffff

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
//TODO deal with M >= pow(2,16) : use 128_t end:64 | period:32 | M:32
// mem_hash_hit: end:32 | period:16 | M:16
// M: MEM hit length
#define _set_mem_hash_hit(end, period, m) (((end) << 32) | ((period) << 16) | m)
#define _get_mem_hash_hit_end(hash_hit, i) ((hash_hit)[i] >> 32)
#define _get_mem_hash_hit_period(hash_hit, i) (((hash_hit)[i] & _32mask) >> 16)
#define _get_mem_hash_hit_start(hash_hit, i) (_get_mem_hash_hit_end((hash_hit), i) - _get_mem_hash_hit_period((hash_hit), i))
#define _get_mem_hash_hit_meml(hash_hit, i) ((hash_hit)[i] & _16mask)


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
    int start, end, mem_l; int score;
    int8_t is_tracked;
} dp_t;

typedef struct {
    cell_t *cell;
    int len;
} chain_t;

typedef struct {
    int start, end;
} rep_reg_t;

//void radix_sort_hash(hash_t *beg, hash_t *end);

int hash_partition(char *seq, int seq_len, tandem_seq_t *tandem_seq, int8_t **hit_array, int *array_m, mini_tandem_para *mtp);

#ifdef __cplusplus
}
#endif

#endif
