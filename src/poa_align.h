#ifndef POA_ALIGN_H
#define POA_ALIGN_H

#include "poa.h"
#include "poa_graph.h"

// SPOA
//#define POA_MATCH  5
//#define POA_MISMATCH  4
//#define POA_GAP_OPEN  0
//#define POA_GAP_EXT  8

#define POA_MATCH  1
#define POA_MISMATCH  3
#define POA_GAP_OPEN  5
#define POA_GAP_EXT  2

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

#define HID_OFF_SET 0
#define HOP_OFF_SET   29
#define EID_OFF_SET 31
#define EOP_OFF_SET   60
#define FOP_OFF_SET   62

#define INF_16_MIN -0x2000
#define INF_16_MAX 0x7fff
#define INF_32_MIN -0x20000000
#define INF_32_MAX 0x7fffffff
#define INF_64_MIN -0x2000000000000000
#define INF_64_MAX 0x7fffffffffffffff

#define _dp_cell(i, j, col_n, dp_matrix) dp_matrix[i * col_n + j]

#define GET_DP_BEGIN(graph, w, i) MAX_OF_TWO(0, MIN_OF_TWO(poa_graph_index_to_min_rank(graph, i), qlen - poa_graph_index_to_max_remain(graph, i))-w)
#define GET_DP_END(graph, w, i) MIN_OF_TWO(qlen, MAX_OF_TWO(poa_graph_index_to_max_rank(graph, i), qlen - poa_graph_index_to_min_remain(graph, i))+w);

#ifdef __cplusplus
extern "C" {
#endif

void gen_simple_mat(int m, int8_t *mat, int8_t match, int8_t mismatch);
//poa_graph_cigar_t *poa_init_graph_cigar(int n);
//void poa_free_graph_cigar(poa_graph_cigar_t *poa_cigar, int n);

/* Adaptive banded global partial order graph alignment */
int poa_ada_global_align_sequence_with_graph(poa_graph_t *graph, uint8_t *query, int qlen, poa_para_t *ppt, int *n_cigar, poa_cigar_t **graph_cigar);

/* Adaptive banded global partial order graph alignment */
int poa_ada_forefront_global_align_sequence_with_graph(poa_graph_t *graph, uint8_t *query, int qlen, poa_para_t *ppt, int *n_cigar, poa_cigar_t **graph_cigar);


/* Banded global partial order graph alignment */
int poa_global_align_sequence_with_graph(poa_graph_t *graph, uint8_t *query, int qlen, poa_para_t *ppt, int *n_cigar, poa_cigar_t **graph_cigar);

static inline poa_cigar_t *poa_push_cigar(int *n_cigar, int *m_cigar, poa_cigar_t *cigar, int op, int len, int32_t node_id, int32_t query_id) {
    poa_cigar_t l = len;
    if (*n_cigar == 0 || (op != POA_CINS && op != POA_CSOFT_CLIP && op != POA_CHARD_CLIP) || op != (cigar[(*n_cigar)-1] & 0xf)) {
        if (*n_cigar == *m_cigar) {
            *m_cigar = *m_cigar? (*m_cigar)<<1 : 4;
            cigar = (poa_cigar_t*)_err_realloc(cigar, (*m_cigar) * sizeof(poa_cigar_t));
        }
        poa_cigar_t n_id = node_id, q_id = query_id;
        if (op == POA_CMATCH || op == POA_CDIFF) 
            cigar[(*n_cigar)++] = n_id << 34 | q_id << 4 | op;
        else if (op == POA_CINS || op == POA_CSOFT_CLIP || op == POA_CHARD_CLIP) 
            cigar[(*n_cigar)++] = q_id << 34 | l << 4 | op;
        else if (op == POA_CDEL)
            cigar[(*n_cigar)++] = n_id << 34 | l << 4 | op;
        else
            err_fatal(__func__, "Unknown cigar operation: %s\n", op);
    } else cigar[(*n_cigar)-1] += l << 4;

    return cigar;
}

static inline poa_cigar_t *poa_reverse_cigar(int n_cigar, poa_cigar_t *cigar) {
    int i; poa_cigar_t tmp;
    for (i = 0; i < n_cigar >> 1; ++i) {
        tmp = cigar[i];
        cigar[i] = cigar[n_cigar-1-i];
        cigar[n_cigar-1-i] = tmp;
    }
    return cigar;
}

static inline void poa_print_cigar(int n_cigar, poa_cigar_t *cigar) {
    int i, op, len, node_id, query_id;
    int n[6] = {0, 0, 0, 0, 0, 0};
    for (i = 0; i < n_cigar; ++i) {
        op = cigar[i] & 0xf; node_id = (int)(cigar[i] >> 34); 
        len = query_id = (int)(cigar[i] >> 4) & 0x3fffffff;
        if (op == POA_CMATCH || op == POA_CDIFF) {
            printf("1%c:%d,%d\t", POA_CIGAR_STR[op], node_id, query_id);
            n[op] += 1;
        } else if (op == POA_CDEL) {
            printf("%d%c:%d\t", len, POA_CIGAR_STR[op], node_id);
            n[op] += len;
        } else if (op == POA_CINS || op == POA_CSOFT_CLIP || op == POA_CHARD_CLIP) { 
            printf("%d%c:%d\t", len, POA_CIGAR_STR[op], node_id);
            n[op] += len;
        } else {
            err_fatal(__func__, "Unknown cigar operation: %s\n", op);
        }
    } printf("\n");
    for (i = 0; i < 6; ++i)
        printf("%d%c ", n[i], POA_CIGAR_STR[i]);
    printf("\n");
}

static inline void poa_backtrack(uint64_t *backtrack_z, int best_i, int best_j, int z_col_n, poa_graph_t *graph, int *n_cigar, poa_cigar_t **graph_cigar) {
    int i, j;
    if (n_cigar && graph_cigar) {
        int n_c = 0, m_c = 0, id, which;
        int op_shift[4] = {HOP_OFF_SET, EOP_OFF_SET, FOP_OFF_SET, HOP_OFF_SET}, id_shift[4] = {HID_OFF_SET, EID_OFF_SET, 0, HID_OFF_SET};
        uint64_t d;
        poa_cigar_t *cigar = 0;
        i = best_i, j = best_j, id = poa_graph_index_to_node_id(graph, i), which = 0;
        while (i > 0 && j > 0) {
            d = backtrack_z[(long)(i-1) * z_col_n + j-1];
            which = (d >> op_shift[which]) & 3;
            // printf("(%d,%d) %c\n", i, j, POA_CIGAR_STR[which]);
            if (which == 0) { // match
                cigar = poa_push_cigar(&n_c, &m_c, cigar, POA_CMATCH, 1, id, j-1);
                i = (d >> id_shift[which]) & 0x1fffffff;
                id = poa_graph_index_to_node_id(graph, i);
                j--;
            } else if (which == 3) { // mismatch
                cigar = poa_push_cigar(&n_c, &m_c, cigar, POA_CDIFF, 1, id, j-1);
                i = (d >> id_shift[which]) & 0x1fffffff;
                id = poa_graph_index_to_node_id(graph, i);
                j--;
            } else if (which == 1) { // deletion
                cigar = poa_push_cigar(&n_c, &m_c, cigar, POA_CDEL, 1, id, j-1);
                i = (d >> id_shift[which]) & 0x1fffffff;
                id = poa_graph_index_to_node_id(graph, i);
            } else { // insertion
                cigar = poa_push_cigar(&n_c, &m_c, cigar, POA_CINS, 1, id, j-1);
                j--;
            }
        }
        if (j > 0) cigar = poa_push_cigar(&n_c, &m_c, cigar, POA_CSOFT_CLIP, j, -1, j-1);
        // reverse cigar
        *graph_cigar = poa_reverse_cigar(n_c, cigar);
        *n_cigar = n_c;
#ifdef __DEBUG__
        poa_print_cigar(n_c, *graph_cigar);
#endif
    }
}

#define _set_max_score(best_score, best_i, best_j, score, i, j) { \
    if (score > best_score) {                                     \
        best_score = score; best_i = i; best_j = j;               \
    }                                                             \
}

#ifdef __cplusplus
}
#endif

#endif
