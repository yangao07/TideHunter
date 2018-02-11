#ifndef POA_ALIGN_H
#define POA_ALIGN_H

#include "poa.h"
#include "poa_graph.h"

#define POA_MATCH 5
#define POA_MISMATCH 4
#define POA_GAP_OPEN 0
#define POA_GAP_EXT 8

#define POA_GLOBAL_FLAG 0
#define POA_LOCAL_FLAG 1
#define POA_EXTEND_FLAG 2

#define POA_CIGAR_STR "MIDXSH"
#define POA_CMATCH     0
#define POA_CINS       1
#define POA_CDEL       2
#define POA_CDIFF      3
#define POA_CSOFT_CLIP 4
#define POA_CHARD_CLIP 5


#define INF_16_MIN -0x4000
#define INF_16_MAX 0x7fff
#define INF_32_MIN -0x40000000
#define INF_32_MAX 0x7fffffff
#define INF_64_MIN -0x4000000000000000
#define INF_64_MAX 0x7fffffffffffffff

#ifdef __cplusplus
extern "C" {
#endif

void gen_simple_mat(int m, int8_t *mat, int8_t match, int8_t mismatch);
//poa_graph_cigar_t *poa_init_graph_cigar(int n);
//void poa_free_graph_cigar(poa_graph_cigar_t *poa_cigar, int n);


/* Adaptive banded global partial order graph alignment */
int poa_ada_global_align_sequence_with_graph(poa_graph_t *graph, uint8_t *query, int qlen, poa_para_t *ppt, int *n_cigar, poa_cigar_t **graph_cigar);

/* Banded global partial order graph alignment */
int poa_global_align_sequence_with_graph(poa_graph_t *graph, uint8_t *query, int qlen, poa_para_t *ppt, int *n_cigar, poa_cigar_t **graph_cigar);

int poa_local_align_sequence_with_graph(poa_graph_t *graph, uint8_t *query, int qlen, poa_para_t *ppt, int *n_cigar, poa_cigar_t **graph_cigar);

int poa_extend_align_sequence_with_graph(poa_graph_t *graph, uint8_t *query, int qlen, poa_para_t *ppt, int *n_cigar, poa_cigar_t **graph_cigar);


#ifdef __cplusplus
}
#endif

#endif
