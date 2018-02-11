#include <stdio.h>
#include <stdlib.h>
#include "poa_align.h"
#include "poa.h"
#include "utils.h"

#ifdef __GNUC__
#define LIKELY(x) __builtin_expect((x),1)
#define UNLIKELY(x) __builtin_expect((x),0)
#else
#define LIKELY(x) (x)
#define UNLIKELY(x) (x)
#endif

//poa_graph_cigar_t *poa_init_graph_cigar(void) {
//    poa_graph_cigar_t *poa_cigar = (poa_graph_cigar_t*)_err_malloc(sizeof(poa_graph_cigar_t));
//    poa_cigar->cigar_n = 0; poa_cigar->cigar_m = 1; poa_cigar->cigar = (uint32_t*)_err_malloc(sizeof(uint32_t));
//    poa_cigar->node_n = 0; poa_cigar->node_m = 1; poa_cigar->node_i = (int*)_err_malloc(sizeof(int));
//    return poa_cigar;
//}

//void poa_free_graph_cigar(poa_graph_cigar_t *poa_cigar) {
//    if (poa_cigar->node_m > 0) free(poa_cigar->node_i);
//    if (poa_cigar->cigar_m > 0) free(poa_cigar->cigar);
//    free(poa_cigar);
//}

typedef struct {
	int32_t h, e;
} eh_t;

typedef struct {
    int32_t h, e, f; // h: max, e: vertical, f: horizontal
} dp_matrix_t;


void gen_simple_mat(int m, int8_t *mat, int8_t match, int8_t mismatch) {
	int i, j;
	match = match < 0? -match : match;
	mismatch = mismatch > 0? -mismatch : mismatch;
	for (i = 0; i < m - 1; ++i) {
		for (j = 0; j < m - 1; ++j)
			mat[i * m + j] = i == j? match : mismatch;
		mat[i * m + m - 1] = 0;
	}
	for (j = 0; j < m; ++j)
		mat[(m - 1) * m + j] = 0;
}

/* void poa_push_graph_cigar_id(poa_graph_cigar_t *cigar, int32_t node_id , int32_t query_id) {
    if (cigar->len == cigar->m) {
        cigar->m = cigar->m ? (cigar->m) << 1 : 4;
        cigar->node_id = (int32_t*)_err_realloc(cigar->node_id, (cigar->m) * sizeof(int32_t));
        cigar->query_id = (int32_t*)_err_realloc(cigar->query_id, (cigar->m) * sizeof(int32_t));
    }
    cigar->node_id[cigar->len] = node_id;
    cigar->query_id[cigar->len++] = query_id;
} 

poa_graph_cigar_t *poa_push_graph_cigar(int *n_cigar, int *m_cigar, poa_graph_cigar_t *cigar, int op, int32_t node_id, int32_t query_id) {
    if (*n_cigar == 0 || op != (cigar[(*n_cigar)].op)) {
        if (*n_cigar == *m_cigar) {
            *m_cigar = *m_cigar? (*m_cigar)<<1 : 4;
            cigar = (poa_graph_cigar_t*)_err_realloc(cigar, (*m_cigar) * sizeof(poa_graph_cigar_t));
        }
        cigar[(*n_cigar)].len = cigar[(*n_cigar)].m = 0;
        poa_push_graph_cigar_id(cigar + *n_cigar, node_id, query_id);
        cigar[(*n_cigar)++].op = op;
    } else poa_push_graph_cigar_id(cigar + *n_cigar, node_id, query_id);
    return cigar;
} */

poa_cigar_t *poa_push_cigar(int *n_cigar, int *m_cigar, poa_cigar_t *cigar, int op, int len, int32_t node_id, int32_t query_id) {
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

poa_cigar_t *poa_reverse_cigar(int n_cigar, poa_cigar_t *cigar) {
    int i; poa_cigar_t tmp;
    for (i = 0; i < n_cigar >> 1; ++i) {
        tmp = cigar[i];
        cigar[i] = cigar[n_cigar-1-i];
        cigar[n_cigar-1-i] = tmp;
    }
    return cigar;
}

void poa_print_cigar(int n_cigar, poa_cigar_t *cigar) {
    int i, op, len, node_id, query_id;
    for (i = 0; i < n_cigar; ++i) {
        op = cigar[i] & 0xf; node_id = (int)(cigar[i] >> 34); 
        len = query_id = (int)(cigar[i] >> 4) & 0x3fffffff;
        if (op == POA_CMATCH || op == POA_CDIFF) {
            printf("1%c:%d,%d\n", POA_CIGAR_STR[op], node_id, query_id);
        } else if (op == POA_CDEL) {
            printf("%d%c:%d\n", len, POA_CIGAR_STR[op], node_id);
        } else if (op == POA_CINS || op == POA_CSOFT_CLIP || op == POA_CHARD_CLIP) { 
            printf("%d%c:%d\n", len, POA_CIGAR_STR[op], node_id);
        } else {
            err_fatal(__func__, "Unknown cigar operation: %s\n", op);
        }
    }
}

int poa_ada_global_align_sequence_with_graph(poa_graph_t *graph, uint8_t *query, int qlen, poa_para_t *ppt, int *n_cigar, poa_cigar_t **graph_cigar) {
    return 0;
}
int poa_forefront_global_align_sequence_with_graph(poa_graph_t *graph, uint8_t *query, int qlen, poa_para_t *ppt, int *n_cigar, poa_cigar_t **graph_cigar) {
    return 0;
}


// TODO use gssw; try adaptive banded DP;
// 0. init score table
// 1. traverse graph with topological order; banded width
// 2. backtrack and generate graph_cigar
// h: max score, e: score from vertical, f: score from horizontal
int poa_global_align_sequence_with_graph(poa_graph_t *graph, uint8_t *query, int qlen, poa_para_t *ppt, int *n_cigar, poa_cigar_t **graph_cigar) {
    dp_matrix_t *dp_matrix; // Full: (tlen + 1) * (qlen + 1); Banded: (tlen+1) * (2 * w + 1)
    int8_t *qp, *mat = ppt->mat; // query profile
    int i, j, k, gap_o = ppt->gap_open, gap_e = ppt->gap_ext, gap_oe = ppt->gap_open + ppt->gap_ext, w, matrix_row_n, matrix_col_n, z_col_n;
    int node_i, target_node_n = graph->node_n - 2; // exclude start and end nodes
    uint64_t *z; // backtrack matrix; in each cell: hd << 33 | ed << 2 | fd
                 //                                 h<<62|h_id<<33|e<<31|e_i<<2|f
                 //                                 h_id/e_i: 29 bit XXX cause error when in_edge_n >= pow(2,29)
                 //                                 MATCH:0, DELETION:1, INSERTION:2, MISMATCH:3

    // allocate memory 
    matrix_row_n = target_node_n + 1, matrix_col_n = qlen + 1;
    dp_matrix = (dp_matrix_t*)_err_calloc(matrix_row_n * matrix_col_n, sizeof(dp_matrix_t));
    w = ppt->bw <= 0 ? qlen : ppt->bw; // band width TODO adaptive band width, shift with diagonal
    z_col_n = qlen < 2 * w + 1 ? qlen : 2 * w + 1; // TODO adaptive band width
    z = graph_cigar ? (uint64_t*)_err_malloc(z_col_n * target_node_n * sizeof(uint64_t)) : 0; // FIXME
    qp = (int8_t*)_err_malloc(qlen * ppt->m);

    // generate the query profile
    for (k = i = 0; k < ppt->m; ++k) {
        const int8_t *p = &mat[k * ppt->m];
        for (j = 0; j < qlen; ++j) qp[i++] = p[query[j]];
    }
    // fill the first row
    dp_matrix[0].h = 0; dp_matrix[0].e = INF_32_MIN; // f is useless
    for (j = 1; j < matrix_col_n && j <= w; ++j) 
        dp_matrix[j].h = dp_matrix[j].f = -(gap_o + gap_e * j), dp_matrix[j].e = INF_32_MIN;
    for (; j < matrix_col_n; ++j) dp_matrix[j].h = dp_matrix[j].e = INF_32_MIN; // everything outside the band is -INF
    // DP loop
    for (i = 1; LIKELY(i <= target_node_n); ++i) { // target graph is in the outer loop
        node_i = poa_graph_rank_to_node_id(graph, i); // i: node rank, node_i: node id
        int32_t f = INF_32_MIN, h1, beg, beg_, end;
        int8_t *q = &qp[graph->node[node_i].base * qlen];
        beg = i > w ? i - w : 0; // set band boundary FIXME for multi-in node
        end = i + w + 1 < qlen + 1 ? i + w + 1: qlen + 1; // only loop through [beg,end) of the query sequence 
		h1 = beg == 0 ? -(gap_o + gap_e * i) : INF_32_MIN;
        int pre_n = graph->node[node_i].in_edge_n;
        int *pre_i = graph->node[node_i].in_id;
        int *pre_rank = (int*)_err_malloc(pre_n * sizeof(int));
        for (j = 0; j < pre_n; ++j)
            pre_rank[j] = poa_graph_node_id_to_rank(graph, pre_i[j]);

        dp_matrix_t *pre_matrix;
        dp_matrix_t *cur_matrix = &dp_matrix[i * matrix_col_n];

        if (n_cigar && graph_cigar) { // store backtrack direction
            uint64_t *zi = &z[(long)(i-1) * z_col_n];
            int32_t m, h, e, tmp; uint64_t pre_i, m_i, fd, ed, e_i, hd, mx, m0=0, e1=1, f2=2, x3=3; // h<<62|h_id<<33|e<<31|e_i<<2|f
            // fill the first column
            if (beg == 0) {
                cur_matrix[beg].h = -(gap_o + gap_e * i);
                cur_matrix[beg].e = -(gap_o + gap_e * i);
                cur_matrix[beg].f = INF_32_MIN;
                beg_ = 1;
            } else beg_ = beg;
            for (j = beg_; LIKELY(j < end); ++j) {
                k = 0; // first precursor
                pre_matrix = &dp_matrix[pre_rank[k] * matrix_col_n];
                m_i = e_i = pre_rank[k];
                //cur_matrix[j].f = f = MAX_OF_TWO(cur_matrix[j-1].f, cur_matrix[j-1].h-gap_o) - gap_e;
                tmp = cur_matrix[j-1].h - gap_oe;
                f = cur_matrix[j-1].f - gap_e;
                cur_matrix[j].f = f > tmp ? f : tmp;
                f = cur_matrix[j].f;

                //cur_matrix[j].e = e = MAX_OF_TWO(pre_matrix[j].e, pre_matrix[j].h-gap_o) - gap_e;
                tmp = pre_matrix[j].h - gap_oe;
                e = pre_matrix[j].e - gap_e;
                cur_matrix[j].e = e > tmp ? e : tmp;
                e = cur_matrix[j].e;

                m = pre_matrix[j-1].h + q[j-1];
                mx = q[j-1] == ppt->match ? m0 : x3;
                // hd for current cell
                hd = m >= e ? (mx << 62 | m_i << 33) : (e1 << 62 | e_i << 33);
                h = m >= e ? m : e;
                hd = h >= f ? hd : f2 << 62 ;
                h = h >= f ? h : f;
                cur_matrix[j].h = h;
                
                // ed for next cell
                ed = e - gap_e > m - gap_oe ? (e1 << 31 | e_i << 2) : (mx << 31 | m_i << 2);
                // fd for next cell
                fd = f - gap_e > m - gap_oe ? f2 : 0;
                // if has multiple precursors
                for (k = 1; k < pre_n; ++k) {
                    pre_matrix = &dp_matrix[pre_rank[k] * matrix_col_n];
                    pre_i = pre_rank[k];

                    e = pre_matrix[j].h - gap_oe;
                    if (e > cur_matrix[j].e) {
                        cur_matrix[j].e = e;
                        e_i = pre_i;
                    }
                    e = pre_matrix[j].e - gap_e;
                    if (e > cur_matrix[j].e) {
                        cur_matrix[j].e = e;
                        e_i = pre_i;
                    }
                    e = cur_matrix[j].e;

                    if (pre_matrix[j-1].h + q[j-1] > m) {
                        m = pre_matrix[j-1].h + q[j-1];
                        m_i = pre_i;
                        mx = q[j-1] == ppt->match ? m0 : x3;
                    }
                    // hd for current cell
                    hd = m > cur_matrix[j].h ? (m0 << 62 | m_i << 33) : hd;
                    cur_matrix[j].h = m > cur_matrix[j].h ? m : cur_matrix[j].h;

                    hd = cur_matrix[j].h >= e ? hd : (e1 << 62 | e_i << 33);
                    cur_matrix[j].h = cur_matrix[j].h >= e ? cur_matrix[j].h : e;

                    // ed for next cell
                    ed = e - gap_e > m - gap_oe ? (e1 << 31 | e_i << 2) : (mx << 31 | m_i << 2);
                }
                zi[j - beg_] = hd | ed | fd;
                // printf("direction(%d, %d, %d): hd: %d, h_id: %d, ed: %d, e_id: %d, f_id: %d\n", i, j, beg_, (int)(hd >> 62), (int)((hd >> 33) & 0x1fffffff), (int)(ed >> 31), (int)((ed >> 2) & 0x1fffffff), (int)fd);
            }
            // fill the last column for next row
            if (end <= qlen) {
                cur_matrix[end].e = INF_32_MIN;
                cur_matrix[end].h = cur_matrix[end].f = MAX_OF_TWO(cur_matrix[end-1].e, cur_matrix[end-1].h) - gap_e;
            }
        } else { // only compute score
            int32_t m, h, e;
            // fill the first column
            if (beg == 0) {
                cur_matrix[beg].h = -(gap_o + gap_e * i);
                cur_matrix[beg].e = -(gap_o + gap_e * i);
                cur_matrix[beg].f = INF_32_MIN;
                beg_ = 1;
            } else beg_ = beg;
            for (j = beg_; LIKELY(j < end); ++j) {
                k = 0; // first precursor
                pre_matrix = &dp_matrix[pre_rank[k] * matrix_col_n];
                m = pre_matrix[j-1].h + q[j-1];
                cur_matrix[j].e = e = MAX_OF_TWO(pre_matrix[j].e, pre_matrix[j].h-gap_o) - gap_e;
                cur_matrix[j].f = f = MAX_OF_TWO(cur_matrix[j-1].f, cur_matrix[j-1].h-gap_o) - gap_e;
                h = MAX_OF_THREE(m, e, f);
                cur_matrix[j].h = h;
                for (k = 1; k < pre_n; ++k) {
                    pre_matrix = &dp_matrix[pre_rank[k] * matrix_col_n];
                    m = pre_matrix[j-1].h + q[j-1];
                    e = MAX_OF_TWO(cur_matrix[j].e, MAX_OF_TWO(pre_matrix[j].e, pre_matrix[j].h-gap_o) - gap_e);
                    f = MAX_OF_TWO(cur_matrix[j].f, MAX_OF_TWO(cur_matrix[j-1].f, cur_matrix[j-1].h-gap_o) - gap_e);
                    h = MAX_OF_THREE(m, e, f);

                    cur_matrix[j].e = e;
                    cur_matrix[j].f = f;
                    cur_matrix[j].h = h;
                }
            }
            // fill the last column for next row
            if (end <= qlen) {
                cur_matrix[end].e = INF_32_MIN;
                cur_matrix[end].h = cur_matrix[end].f = MAX_OF_TWO(cur_matrix[end-1].e, cur_matrix[end-1].h) - gap_e;
            }
        }
        /*for (j = 0; j <= qlen; ++j) {
            printf("(%d, %d, %d) ", cur_matrix[j].h, cur_matrix[j].e, cur_matrix[j].f);
        } printf("\n");*/
        free(pre_rank);
    }
    // select best score TODO
    int score = dp_matrix[matrix_row_n * matrix_col_n - 1].h, best_i = matrix_row_n - 1, best_j = matrix_col_n - 1; // TODO based on w
    printf("score: %d\n", score);
    // backtrack from best score
    if (n_cigar && graph_cigar) {
        int n_c = 0, m_c = 0, id, which;
        int op_shift[4] = {62, 31, 0, 62}, id_shift[4] = {33, 2, 0, 33};
        int64_t d;
        poa_cigar_t *cigar = 0;
        i = best_i, j = best_j, id = poa_graph_rank_to_node_id(graph, i), which = 0;
        while (i > 0 && j > 0) {
            d = z[(long)(i-1) * z_col_n + (j-1 - (i > w ? i - w : 0))];
            which = (d >> op_shift[which]) & 3;
            // printf("(%d,%d) %c\n", i, j, POA_CIGAR_STR[which]);
            if (which == 0) { // match
                cigar = poa_push_cigar(&n_c, &m_c, cigar, POA_CMATCH, 1, id/*rank_to_id*/, j-1);
                i = (d >> id_shift[which]) & 0x1fffffff;
                id = poa_graph_rank_to_node_id(graph, i);
                j--;
            } else if (which == 3) { // mismatch
                cigar = poa_push_cigar(&n_c, &m_c, cigar, POA_CDIFF, 1, id, j-1);
                i = (d >> id_shift[which]) & 0x1fffffff;
                id = poa_graph_rank_to_node_id(graph, i);
                j--;
            } else if (which == 1) { // deletion
                cigar = poa_push_cigar(&n_c, &m_c, cigar, POA_CDEL, 1, id, j-1);
                i = (d >> id_shift[which]) & 0x1fffffff;
                id = poa_graph_rank_to_node_id(graph, i);
            } else { // insertion
                cigar = poa_push_cigar(&n_c, &m_c, cigar, POA_CINS, 1, id, j-1);
                j--;
            }
        }
        if (j > 0) cigar = poa_push_cigar(&n_c, &m_c, cigar, POA_CSOFT_CLIP, j, -1, j-1);
        // reverse cigar
        *graph_cigar = poa_reverse_cigar(n_c, cigar);
        *n_cigar = n_c;
        poa_print_cigar(n_c, *graph_cigar);
    }
    free(dp_matrix); free(qp); free(z);
    return score;
}

int poa_local_align_sequence_with_graph(poa_graph_t *graph, uint8_t *query, int qlen, poa_para_t *ppt, int *n_cigar, poa_cigar_t **graph_cigar) {
    return 0;
}

int poa_extend_align_sequence_with_graph(poa_graph_t *graph, uint8_t *query, int qlen, poa_para_t *ppt, int *n_cigar, poa_cigar_t **graph_cigar) {
    return 0;
}
