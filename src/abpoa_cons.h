#ifndef _ABPOA_ALIGN_H
#define _ABPOA_ALIGN_H

#include <stdint.h>
#include "tidehunter.h"

#ifdef __cplusplus
extern "C" {
#endif

abpoa_para_t *mt_abpoa_init_para(mini_tandem_para *mtp);
int abpoa_gen_cons(abpoa_t *ab, abpoa_para_t *abpt, uint8_t *bseqs, int seq_len, int *pos, int pos_n, uint8_t *cons_bseq, uint8_t *cons_qual, int min_cov);

#ifdef __cplusplus
}
#endif

#endif
