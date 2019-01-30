#ifndef TANDEM_HIT_H
#define TANDEM_HIT_H
#include <stdint.h>
#include "tide_hunter.h"

typedef int64_t hash_t;
#define _8mask 0xff
#define _16mask 0xffff
#define _32mask 0xffffffff

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

#ifdef __cplusplus
extern "C" {
#endif

int collect_tandem_repeat_hit(uint8_t *bseq, int seq_len, mini_tandem_para *mtp, hash_t **hit_h);

#ifdef __cplusplus
}
#endif

#endif