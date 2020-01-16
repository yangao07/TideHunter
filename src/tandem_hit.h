#ifndef TANDEM_HIT_H
#define TANDEM_HIT_H
#include <stdint.h>
#include "tide_hunter.h"

typedef uint64_t hash_t;
#define _8mask 0xff
#define _16mask 0xffff
#define _32mask 0xffffffff

// TODO total_hits_of_kmer: used as weight of kmer

// end:32 | period:32
#define _set_hash_hit(end, period) (((end) << 32) | (period) )
#define _get_hash_hit_end(hash_hit, i) ((hash_hit)[i] >> 32)
#define _get_hash_hit_period(hash_hit, i) ((uint32_t)(hash_hit)[i])
#define _get_hash_hit_start(hash_hit, i) (_get_hash_hit_end((hash_hit), i) - _get_hash_hit_period((hash_hit), i))

#ifdef __cplusplus
extern "C" {
#endif

int collect_tandem_repeat_hit(uint8_t *bseq, int seq_len, mini_tandem_para *mtp, hash_t **hit_h);

#ifdef __cplusplus
}
#endif

#endif
