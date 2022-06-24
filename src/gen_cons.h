#ifndef GEN_CONS_H
#define GEN_CONS_H

#include "tidehunter.h"

#ifdef __cplusplus
extern "C" {
#endif

void single_copy_full_len_seq(int seq_len, char *seq, tandem_seq_t *tseq, mini_tandem_para *mtp);
void seqs_msa(int seq_len, uint8_t *bseq, int par_n, int *par_pos, tandem_seq_t *tseq, mini_tandem_para *mtp, abpoa_t *ab, abpoa_para_t *abpt);

#ifdef __cplusplus
}
#endif

#endif
