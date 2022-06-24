#ifndef _EDLIB_ALIGN_H
#define _EDLIB_ALIGN_H

#ifdef __cplusplus
extern "C" {
#endif

int edlib_align_HW(char *query, int qlen, char *target, int tlen, int *start, int *end, int k);

// int edlib_align_NW(char *query, int qlen, char *target, int tlen);

#ifdef __cplusplus
}
#endif

#endif
