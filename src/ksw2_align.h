#ifndef _KSW2_ALIGN_H
#define _KSW2_ALIGN_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

int ksw2_global(const uint8_t *query, int qlen, const uint8_t *target, int tlen);

int ksw2_left_ext(const uint8_t *query, int qlen, const uint8_t *target, int tlen);
int ksw2_ext(const uint8_t *query, int qlen, const uint8_t *target, int tlen);

#ifdef __cplusplus
}
#endif

#endif
