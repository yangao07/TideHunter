#ifndef POA_H
#define POA_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    // score matrix
    int m; int8_t *mat;
    int8_t match, mismatch, gap_open, gap_ext;
    int bw; // band width
    int zdrop, end_bonus; // from minimap2
    // alignment mode
    int align_mode; // 0: global, 1: local, 2: extend
    // bits number
    int8_t score_n, id_n;
} poa_para_t;

poa_para_t *poa_para_init(void);
void poa_para_free(poa_para_t *ppt);

// int poa_main(const char *seq_fn, poa_para_t *ppt) { TODO
int poa_main(int seq_n, char (*seq)[100], poa_para_t *ppt);

#ifdef __cplusplus
}
#endif

#endif
