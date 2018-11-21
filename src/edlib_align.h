#ifndef _EDLIB_ALIGN_H
#define _EDLIB_ALIGN_H

#ifdef __cplusplus
extern "C" {
#endif

int edlib_align_HW(char *query, int qlen, char *target, int tlen, int *start, int *end);

#ifdef __cplusplus
}
#endif

#endif
