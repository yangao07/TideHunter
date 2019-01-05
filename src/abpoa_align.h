#ifndef _ABPOA_ALIGN_H
#define _ABPOA_ALIGN_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

int abpoa_msa(uint8_t *bseqs, int seq_len, int *pos, int pos_n, char *cons_seq);

#ifdef __cplusplus
}
#endif

#endif
