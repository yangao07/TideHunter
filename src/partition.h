#ifndef PARTITION_H
#define PARTITION_H
#include "tandem_chain.h"

#ifdef __cplusplus
extern "C" {
#endif

int *get_partition_pos_with_narrow_global_alignment(uint8_t *bseq, int seq_len, dp_t **dp, chain_t ch, mini_tandem_para *mtp, int *par_n);

#ifdef __cplusplus
}
#endif

#endif
