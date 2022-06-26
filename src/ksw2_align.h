#ifndef _KSW2_ALIGN_H
#define _KSW2_ALIGN_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

int ksw2_backtrack_left_end(int n_cigar, uint32_t *cigar, int qlen, int tlen, int q_left_ext);

int ksw2_global(const uint8_t *query, int qlen, const uint8_t *target, int tlen);
int ksw2_global_with_cigar(const uint8_t *query, int qlen, const uint8_t *target, int tlen, int *n_cigar, uint32_t **cigar);

void ksw2_left_ext(const uint8_t *query, int qlen, const uint8_t *target, int tlen, int *max_q, int *max_t);
void ksw2_right_ext(const uint8_t *query, int qlen, const uint8_t *target, int tlen, int *max_q, int *max_t);
int ksw2_left_extend(const uint8_t *query, int qlen, const uint8_t *target, int tlen, int *mqe_tlen);
int ksw2_right_extend(const uint8_t *query, int qlen, const uint8_t *target, int tlen, int *mqe_tlen);

#ifdef __cplusplus
}
#endif

#endif
