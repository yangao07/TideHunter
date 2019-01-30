#ifndef GEN_CONS_H
#define GEN_CONS_H

#include "tide_hunter.h"

#ifdef __cplusplus
extern "C" {
#endif

void seqs_msa(int seq_len, uint8_t *bseq, int par_n, int *par_pos, tandem_seq_t *tseq, mini_tandem_para *mtp);

#ifdef __cplusplus
}
#endif

#endif
