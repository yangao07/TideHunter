#include <stdio.h>
#include <stdlib.h>
#include "poa_align.h"
#include "poa.h"
#include "utils.h"
#include "simd_instruction.h"

#define TODO_M   0x1
#define TODO_X   0x2
#define TODO_MX  0x3
#define TODO_E   0x4
#define TODO_ME  0x5
#define TODO_F   0x8
#define TODO_MF  0x9
#define TODO_MEF 0xd

#define HAS_M    0x10
#define HAS_X    0x20
#define HAS_MX   0x30
#define HAS_E    0x40
#define HAS_ME   0x50
#define HAS_F    0x80


// xl, el, fl: current extension length of X, E and F
// h: max, e: vertical, f: horizontal
typedef struct {
    int64_t h:32, xl:32;
    int64_t e:32, el:32;
    int64_t f:32, fl:32; 
    uint8_t todo_map; // XXX need to stored seperately
} dp_matrix_t;

void gen_simple_mat(int m, int8_t *mat, int8_t match, int8_t mismatch) {
    int i, j;
    match = match < 0 ? -match : match;
    mismatch = mismatch > 0? -mismatch : mismatch;
    for (i = 0; i < m - 1; ++i) {
        for (j = 0; j < m - 1; ++j)
            mat[i * m + j] = i == j ? match : mismatch;
        mat[i * m + m - 1] = 0;
    }
    for (j = 0; j < m; ++j)
        mat[(m - 1) * m + j] = 0;
}

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
    int best_score = INF_32_MIN, best_i=0, best_j=0;

    // allocate memory 
    matrix_row_n = target_node_n + 1, matrix_col_n = qlen + 1;
    dp_matrix = (dp_matrix_t*)_err_calloc(matrix_row_n * matrix_col_n, sizeof(dp_matrix_t));
    w = ppt->bw <= 0 ? qlen : ppt->bw; // band width TODO adaptive band width, shift with diagonal
    z_col_n = qlen < 2 * w + 1 ? qlen : 2 * w + 1; // TODO adaptive band width
    z = graph_cigar ? (uint64_t*)_err_malloc(z_col_n * target_node_n * sizeof(uint64_t)) : 0;
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
    for (i = 1; i <= target_node_n; ++i) { // target graph is in the outer loop
        node_i = poa_graph_index_to_node_id(graph, i); // i: node index, node_i: node id
        int32_t f = INF_32_MIN, beg, beg_, end;
        int8_t *q = &qp[graph->node[node_i].base * qlen];
        beg = i > w ? i - w : 0; // set band boundary FIXME for multi-in node
        end = i + w + 1 < qlen + 1 ? i + w + 1: qlen + 1; // only loop through [beg,end) of the query sequence 
        int pre_n = graph->node[node_i].in_edge_n;
        int *pre_i = graph->node[node_i].in_id;
        int *pre_index = (int*)_err_malloc(pre_n * sizeof(int));
        for (j = 0; j < pre_n; ++j)
            pre_index[j] = poa_graph_node_id_to_index(graph, pre_i[j]);

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
            for (j = beg_; j < end; ++j) {
                k = 0; // first precursor
                pre_matrix = &dp_matrix[pre_index[k] * matrix_col_n];
                m_i = e_i = pre_index[k];
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

                // ed and fd for next cell
                ed = e - gap_e > m - gap_oe ? (e1 << 31 | e_i << 2) : (mx << 31 | m_i << 2);
                fd = f - gap_e > m - gap_oe ? f2 : 0;
                // if has multiple precursors
                for (k = 1; k < pre_n; ++k) {
                    pre_matrix = &dp_matrix[pre_index[k] * matrix_col_n];
                    pre_i = pre_index[k];

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
                    hd = m > cur_matrix[j].h ? (mx << 62 | m_i << 33) : hd;
                    cur_matrix[j].h = m > cur_matrix[j].h ? m : cur_matrix[j].h;

                    hd = cur_matrix[j].h >= e ? hd : (e1 << 62 | e_i << 33);
                    cur_matrix[j].h = cur_matrix[j].h >= e ? cur_matrix[j].h : e;

                    // ed and fd for next cell
                    ed = e - gap_e > m - gap_oe ? (e1 << 31 | e_i << 2) : (mx << 31 | m_i << 2);
                    fd = f - gap_e > m - gap_oe ? f2 : 0;
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
            for (j = beg_; j < end; ++j) {
                k = 0; // first precursor
                pre_matrix = &dp_matrix[pre_index[k] * matrix_col_n];
                m = pre_matrix[j-1].h + q[j-1];
                cur_matrix[j].e = e = MAX_OF_TWO(pre_matrix[j].e, pre_matrix[j].h-gap_o) - gap_e;
                cur_matrix[j].f = f = MAX_OF_TWO(cur_matrix[j-1].f, cur_matrix[j-1].h-gap_o) - gap_e;
                h = MAX_OF_THREE(m, e, f);
                cur_matrix[j].h = h;
                for (k = 1; k < pre_n; ++k) {
                    pre_matrix = &dp_matrix[pre_index[k] * matrix_col_n];
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
        free(pre_index);
    }
    int in_id, in_index;
    for (i = 0; i < graph->node[POA_SINK_NODE_ID].in_edge_n; ++i) {
        in_id = graph->node[POA_SINK_NODE_ID].in_id[i];
        in_index = poa_graph_node_id_to_index(graph, in_id);
        _set_max_score(best_score, best_i, best_j, dp_matrix[in_index * matrix_col_n + qlen].h, in_index, qlen);
    }

    printf("best_score: (%d, %d) -> %d\n", best_i, best_j, best_score);
    // backtrack from best score
    if (n_cigar && graph_cigar) {
        int n_c = 0, m_c = 0, id, which;
        int op_shift[4] = {62, 31, 0, 62}, id_shift[4] = {33, 2, 0, 33};
        uint64_t d;
        poa_cigar_t *cigar = 0;
        i = best_i, j = best_j, id = poa_graph_index_to_node_id(graph, i), which = 0;
        while (i > 0 && j > 0) {
            d = z[(long)(i-1) * z_col_n + (j-1 - (i-1 > w ? i-1 - w : 0))];
            which = (d >> op_shift[which]) & 3;
            // printf("(%d,%d) %c\n", i, j, POA_CIGAR_STR[which]);
            if (which == 0) { // match
                cigar = poa_push_cigar(&n_c, &m_c, cigar, POA_CMATCH, 1, id/*index_to_id*/, j-1);
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
    free(dp_matrix); free(qp); free(z);
    return best_score;
}

/*
int poa_forefront_global_align_sequence_with_graph(poa_graph_t *graph, uint8_t *query, int qlen, poa_para_t *ppt, int *n_cigar, poa_cigar_t **graph_cigar) {
    dp_matrix_t *dp_matrix; // Full: (tlen + 1) * (qlen + 1); Banded: (tlen+1) * (2 * w + 1)
    int8_t *qp, *mat = ppt->mat; // query profile
    int *forefront_start, *forefront_end, forefront_n;
    int **pre_index, *pre_n, pre_i;
    int rank_i, node_id, index_i, q_i;
    int i, j, k, gap_o = ppt->gap_open, gap_e = ppt->gap_ext, gap_oe = ppt->gap_open + ppt->gap_ext, w, matrix_row_n, matrix_col_n, z_col_n;
    int target_node_n = graph->node_n - 2; // exclude start and end nodes
    uint64_t *backtrack_z; // backtrack matrix; in each cell: hd << 33 | ed << 2 | fd
    //                                 h<<62|h_id<<33|e<<31|e_i<<2|f
    //                                 h_id/e_i: 29 bit XXX cause error when in_edge_n >= pow(2,29)
    //                                 MATCH:0, DELETION:1, INSERTION:2, MISMATCH:3
    int best_score = INF_32_MIN, best_i=0, best_j=0;

    // allocate memory 
    matrix_row_n = target_node_n + 1, matrix_col_n = qlen + 1;
    dp_matrix = (dp_matrix_t*)_err_calloc(matrix_row_n * matrix_col_n, sizeof(dp_matrix_t));
    w = ppt->bw <= 0 ? (MIN_OF_TWO(graph->rank_n-2, qlen)) : (MIN_OF_THREE(graph->rank_n-2, qlen, ppt->bw)); // band width TODO adaptive band width, shift with diagonal
    z_col_n = qlen; // TODO use less memory ???
    backtrack_z = graph_cigar ? (uint64_t*)_err_malloc(z_col_n * target_node_n * sizeof(uint64_t)) : 0;

    qp = (int8_t*)_err_malloc(qlen * ppt->m); // TODO score profile for diagnal
    forefront_n = graph->rank_n + qlen - 1;
    forefront_start = (int*)_err_calloc(forefront_n, sizeof(int)); // [i-start, start]
    forefront_end = (int*)_err_calloc(forefront_n, sizeof(int));   // [end, i-end]
    pre_index = (int**)_err_malloc(graph->node_n * sizeof(int*));
    pre_n = (int*)_err_malloc(graph->node_n * sizeof(int*));

    // generate the query profile
    for (k = i = 0; k < ppt->m; ++k) {
        const int8_t *p = &mat[k * ppt->m];
        for (j = 0; j < qlen; ++j) qp[i++] = p[query[j]];
    }
    // index of pre-node
    for (i = 0; i < graph->node_n; ++i) {
        node_id = poa_graph_index_to_node_id(graph, i); // i: node index
        pre_n[i] = graph->node[node_id].in_edge_n;
        pre_index[i] = (int*)_err_malloc(pre_n[i] * sizeof(int));
        for (j = 0; j < pre_n[i]; ++j) {
            pre_index[i][j] = poa_graph_node_id_to_index(graph, graph->node[node_id].in_id[j]);
        }
    }
    int *rank_index = graph->rank_index;
    // initialize edge for forefront (index)
    for (i = 0; i < qlen; ++i) forefront_start[i] = 0; 
    for (i = 0; i < graph->rank_n-1; ++i) forefront_end[i] = rank_index[i]-1;
    for (i = qlen; i < forefront_n; ++i) forefront_start[i] = rank_index[i-qlen]-1; 
    for (i = graph->rank_n-1; i < forefront_n; ++i) forefront_end[i] = rank_index[graph->rank_n-2]-1;
    // fill the first row and column 
    dp_matrix[0].h = 0; dp_matrix[0].e = INF_32_MIN; dp_matrix[0].f = INF_32_MIN;
    for (i = 1; i <= qlen; ++i) {
        dp_matrix[i].e = INF_32_MIN; dp_matrix[i].f = dp_matrix[i].h = -gap_o - gap_e * i;
    }
    for (i = 1; i <= target_node_n; ++i) {
        //for () // XXX
        //index_i = min_rank_to_index(i);
        index_i = i;
        dp_matrix[index_i * matrix_col_n].f = INF_32_MIN; dp_matrix[index_i * matrix_col_n].e = dp_matrix[index_i * matrix_col_n].h = -gap_o - gap_e * i;
    }

    // DP loop in anti-diagnal direction
    int32_t f, h, e, m, tmp, start, end, sum, pre_sum;
    uint64_t m_i, fd, ed, e_i, hd, mx, m0=0, e1=1, f2=2, x3=3; // h<<62|h_id<<33|e<<31|e_i<<2|f
    dp_matrix_t *cur_dp, *m_pre_dp, *e_pre_dp, *f_pre_dp;
    for (sum = 1; sum < forefront_n; ++sum) {
        start = forefront_start[sum] == 0 ? 1 : forefront_start[sum]; 
        end = poa_graph_index_to_rank(graph, forefront_end[sum]) == sum ? rank_index[sum-1]-1 : forefront_end[sum];
        for (index_i = start; index_i <= end; ++index_i) {
            node_id = poa_graph_index_to_node_id(graph, index_i);
            rank_i = poa_graph_index_to_rank(graph, index_i);
            q_i = sum - rank_i;
            cur_dp = dp_matrix+ index_i * matrix_col_n + q_i;

            // F: from (index_i, q_i-1)
            f_pre_dp = dp_matrix+ index_i * matrix_col_n + q_i - 1;

            if (index_i < forefront_start[sum-1] || index_i > forefront_end[sum-1]) cur_dp->f = INF_32_MIN;
            else {                                                             
                // XXX based on fd
                tmp = f_pre_dp->h - gap_oe;                
                f = f_pre_dp->f - gap_e;              
                cur_dp->f = f > tmp ? f : tmp;       
            }
            pre_i = pre_index[index_i][0]; // index 
            pre_sum = q_i + poa_graph_index_to_rank(graph, pre_i);
            m_i = e_i = pre_i;
            // E: from (pre_i, q_i)
            e_pre_dp = dp_matrix + pre_i * matrix_col_n + q_i;
            if (pre_i < forefront_start[pre_sum] || pre_i > forefront_end[pre_sum]) cur_dp->e = INF_32_MIN;
            else {
                // XXX based on ed
                tmp = e_pre_dp->h - gap_oe;
                e = e_pre_dp->e - gap_e;
                cur_dp->e = e > tmp ? e : tmp;
            }
            // M: from (pre_i, q_i-1)
            m_pre_dp = dp_matrix + pre_i * matrix_col_n + q_i - 1;
            if (pre_sum < 1 || pre_i < forefront_start[pre_sum-1] || pre_i > forefront_end[pre_sum-1]) {
                mx = x3;
                m = INF_32_MIN;
            } else {
                mx = graph->node[node_id].base == query[q_i-1] ? m0 : x3;
                m = m_pre_dp->h + (mx == m0 ? ppt->match : -ppt->mismatch); // XXX score profile for diagnal
            }

            f = cur_dp->f; e = cur_dp->e;

            // hd for current cell
            hd = m >= e ? (mx << 62 | m_i << 33) : (e1 << 62 | e_i << 33);
            h = m >= e ? m : e;
            hd = h >=f ? hd : f2 << 62;
            h = h >= f ? h : f;
            cur_dp->h = h;

            // ed and fd for next cell
            ed = e - gap_e > m - gap_oe ? (e1 << 31 | e_i << 2) : (mx << 31 | m_i << 2);
            fd = f - gap_e > m - gap_oe ? f2 : 0;
            for (j = 1; j < pre_n[index_i]; ++j) {
                pre_i = pre_index[index_i][j];
                pre_sum = q_i + poa_graph_index_to_rank(graph, pre_i);
                // check e TODO USE d of (pre_i,q_i)
                if (pre_i >= forefront_start[pre_sum] && pre_i <= forefront_end[pre_sum]) {
                    e_pre_dp = dp_matrix + pre_i * matrix_col_n + q_i;
                    if (e_pre_dp->h - gap_oe > e) {
                        e = cur_dp->e = e_pre_dp->h - gap_oe;
                        e_i = pre_i;
                    }
                    if (e_pre_dp->e - gap_e > e) {
                        e = cur_dp->e = e_pre_dp->e - gap_e;
                        e_i = pre_i;
                    }
                    if (e > cur_dp->h) {
                        // hd for current cell
                        hd = (e1 << 62 | e_i << 33);
                        cur_dp->h = e;
                    }
                }
                // check m
                if (pre_sum >= 1 && pre_i >= forefront_start[pre_sum-1] && pre_i <= forefront_end[pre_sum-1]) {
                    m_pre_dp = dp_matrix + pre_i * matrix_col_n + q_i - 1;
                    tmp = m_pre_dp->h + (mx == m0 ? ppt->match : -ppt->mismatch); 
                    if (tmp > m) {
                        m = tmp;
                        m_i = pre_i;
                        // hd for current cell
                        if (m > cur_dp->h) {
                            hd =  mx << 62 | m_i << 33;
                            cur_dp->h = m;
                            // fd for next cell
                            fd = f - gap_e > m - gap_oe ? f2 : 0;
                        }
                    }
                }
                // ed for next cell
                ed = e - gap_e > m - gap_oe ? (e1 << 31 | e_i << 2) : (mx << 31 | m_i << 2);

            }
            backtrack_z[(index_i-1) * qlen + q_i-1] = hd | ed | fd;
        }
    }
#ifdef __DEBUG__
    for (j = 0; j <= target_node_n; ++j) {
        printf("rank: %d\t", j);
        for (i = 0; i <= qlen; ++i) {
            printf("%d:(%d,%d,%d)\t", i, dp_matrix[j * matrix_col_n + i].h, dp_matrix[j * matrix_col_n + i].e, dp_matrix[j * matrix_col_n + i].f);
        } printf("\n");
    }
#endif
    int in_id, in_index;
    for (i = 0; i < graph->node[POA_SINK_NODE_ID].in_edge_n; ++i) {
        in_id = graph->node[POA_SINK_NODE_ID].in_id[i];
        in_index = poa_graph_node_id_to_index(graph, in_id);
        _set_max_score(best_score, best_i, best_j, dp_matrix[in_index * matrix_col_n + qlen].h, in_index, qlen);
    }

    printf("best_score: (%d, %d) -> %d\n", best_i, best_j, best_score);
    // backtrack from best score
    if (n_cigar && graph_cigar) {
        int n_c = 0, m_c = 0, id, which;
        int op_shift[4] = {62, 31, 0, 62}, id_shift[4] = {33, 2, 0, 33};
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
    free(dp_matrix); free(qp); free(backtrack_z); free(forefront_start); free(forefront_end);
    for (i = 0; i < graph->node_n; ++i) { free(pre_index[i]); } free(pre_index); free(pre_n);
    return best_score;

    return 0;
}*/
// TODO extend
int poa_ada_extend_align_sequence_with_graph(poa_graph_t *graph, uint8_t *query, int qlen, poa_para_t *ppt, int *n_cigar, poa_cigar_t **graph_cigar) {
    dp_matrix_t *dp_matrix, *cur_line, *pre_line, *next_line, *cur_dp; // Full: (tlen + 1) * (qlen + 1); Banded: (tlen+1) * (2 * w + 1)
    int8_t *qp, *mat = ppt->mat, gap_o = ppt->gap_open, gap_e = ppt->gap_ext, gap_oe = ppt->gap_open + ppt->gap_ext; // query profile
    int **pre_index, *pre_n, pre_i, **next_index, *next_n, next_i;
    int node_id, index_i, min_rank, q_i;
    int i, j, k, w, *dp_beg_cen, *dp_end_cen;
    int target_node_n = graph->node_n - 2; // exclude start and end nodes
    int matrix_row_n = graph->node_n, matrix_col_n = qlen + 2, z_col_n = qlen; // TODO use less memory ???
    uint64_t *backtrack_z, *z; // backtrack matrix; in each cell: hd << 33 | ed << 2 | fd
    //                                 h<<62|h_id<<33|e<<31|e_i<<2|f
    //                                 h_id/e_i: 29 bit XXX cause error when in_edge_n >= pow(2,29)
    //                                 MATCH:0, DELETION:1, INSERTION:2, MISMATCH:3
    int best_score = INF_32_MIN, best_i=0, best_j=0;

    // allocate memory 
    dp_matrix = (dp_matrix_t*)_err_calloc(matrix_row_n * matrix_col_n, sizeof(dp_matrix_t)); // make sure no invalid write 
    dp_beg_cen = (int*)_err_malloc((matrix_row_n) * sizeof(int)); for (i = 1; i < matrix_row_n; ++i) dp_beg_cen[i] = qlen;
    dp_end_cen = (int*)_err_calloc(matrix_row_n, sizeof(int));
    // TODO if no band is used
    // w = ppt->bw <= 0 ? (MIN_OF_TWO(graph->rank_n-2, qlen)) : (MIN_OF_THREE(graph->rank_n-2, qlen, ppt->bw));
    w = 300; // XXX small may cause NO alignment result XXX
    backtrack_z = graph_cigar ? (uint64_t*)_err_malloc(z_col_n * target_node_n * sizeof(uint64_t)) : 0;

    qp = (int8_t*)_err_malloc(qlen * ppt->m); // TODO score profile for diagnal

    pre_index = (int**)_err_malloc(graph->node_n * sizeof(int*));
    pre_n = (int*)_err_malloc(graph->node_n * sizeof(int*));
    next_index = (int**)_err_malloc(graph->node_n * sizeof(int*));
    next_n = (int*)_err_malloc(graph->node_n * sizeof(int*));

    // generate the query profile
    for (k = i = 0; k < ppt->m; ++k) {
        const int8_t *p = &mat[k * ppt->m];
        for (j = 0; j < qlen; ++j) qp[i++] = p[query[j]];
    }

    // index of pre-node
    for (i = 0; i < graph->node_n; ++i) {
        node_id = poa_graph_index_to_node_id(graph, i); // i: node index
        pre_n[i] = graph->node[node_id].in_edge_n;
        pre_index[i] = (int*)_err_malloc(pre_n[i] * sizeof(int));
        for (j = 0; j < pre_n[i]; ++j) {
            pre_index[i][j] = poa_graph_node_id_to_index(graph, graph->node[node_id].in_id[j]);
        }
        
        next_n[i] = graph->node[node_id].out_edge_n;
        next_index[i] = (int*)_err_malloc(next_n[i] * sizeof(int));
        for (j = 0; j < next_n[i]; ++j) {
            next_index[i][j] = poa_graph_node_id_to_index(graph, graph->node[node_id].out_id[j]);
        }
    }

    // TODO no need to store F ???
    // cur_cell: H[i,j], E[i+1,j], F[i,j+1]
    // fill the first row
    cur_line = dp_matrix;
    cur_line->h = 0; cur_line->e = -gap_oe; cur_line->f = -gap_oe;
    cur_line->xl = 0; cur_line->fl = cur_line->el = 1; 
    cur_line->todo_map = HAS_MX; 
    for (i = 1; i <= w; ++i) {
        cur_dp = cur_line + i;
        cur_dp->e = INF_32_MIN; cur_dp->h = -gap_o - gap_e * i; cur_dp->f = cur_dp->h - gap_e;
        cur_dp->fl = i+1; cur_dp->el = cur_dp->xl = 0; 
        cur_dp->todo_map = HAS_MX; 
    }
    for (i = 0; i < next_n[POA_SRC_NODE_ID]; ++i) {
        next_i = next_index[POA_SRC_NODE_ID][i];
        next_line = dp_matrix + next_i * matrix_col_n;
        next_line[0].todo_map = TODO_E;
        for (j = 1; j<= w+1; ++j) {
            next_line[j].todo_map = TODO_MF;
        }
        dp_beg_cen[next_i] = 1;
        dp_end_cen[next_i] = 1;
    }

    // DP loop
    int32_t tmp, h, m, e, f, beg, end, _beg, _end, _next_cen;
    uint64_t m_i, fd, ed, e_i, hd, mx, m0=0x0, e1=0x1, f2=0x2, x3=0x3; // h<<62|h_id<<33|e<<31|e_i<<2|f TODO f<<62|e<<60|e_id<<31|h<<29|h_id
                                                                       // m0=0x0 << 62, x3=0x3<<62
    dp_matrix_t *m_pre_dp, *e_pre_dp, *f_pre_dp;

    for (index_i = 1; index_i < matrix_row_n-1; ++index_i) {
        cur_line = &dp_matrix[index_i * matrix_col_n];
        node_id = poa_graph_index_to_node_id(graph, index_i);
        min_rank = poa_graph_index_to_min_rank(graph, index_i);
        int8_t *q = &qp[graph->node[node_id].base * qlen];
        beg = dp_beg_cen[index_i] >= w ? dp_beg_cen[index_i] - w : 0;
        end = dp_end_cen[index_i] + w <= qlen ? dp_end_cen[index_i] + w : qlen;

        if (beg == 0) {
            cur_line[0].h = -(gap_o + gap_e * min_rank);
            cur_line[0].e = -(gap_o + gap_e * (min_rank+1));
            cur_line[0].f = INF_32_MIN;
            _beg = 1;
        } else _beg = beg;
        for (q_i = _beg; q_i <= end; ++q_i) {
            if (!cur_line[q_i].todo_map) continue;
            cur_dp = cur_line + q_i;
            m = e = f = INF_32_MIN;
            m_i = e_i = -1;
            z = &backtrack_z[(index_i-1) * qlen];

            if (cur_dp->todo_map & TODO_M) { // M: from (pre_i, q_i-1)
                for (i = 0; i < pre_n[index_i]; ++i) {
                    pre_i = pre_index[index_i][i];
                    pre_line = dp_matrix + pre_i * matrix_col_n;
                    m_pre_dp = pre_line + q_i-1;
                    if ((m_pre_dp->todo_map & HAS_M) && m_pre_dp->h > m) {
                        m = m_pre_dp->h; m_i = pre_i;
                    }
                }
            }
            if (cur_dp->todo_map & TODO_E) { // E: from (pre_i, q_i)
                for (j = 0; j < pre_n[index_i]; ++j) {
                    pre_i = pre_index[index_i][j];
                    pre_line = dp_matrix + pre_i * matrix_col_n;
                    e_pre_dp = pre_line + q_i;
                    if ((e_pre_dp->todo_map & HAS_E) && e_pre_dp->e > e) {
                        e = e_pre_dp->e; e_i = pre_i;
                    }
                }
            }
            if (cur_dp->todo_map & TODO_F) { // F: from (index_i, q_i-1)
                f_pre_dp = cur_line + q_i - 1;
                f = f_pre_dp->f;
            }

            m += q[q_i-1]; mx = q[q_i-1] == ppt->match ? m0 : x3;
            // since we have score of M, E and F, now we need to set H for cur cell, set E and F for next cell
            // TODO z-drop
             
            // h,hd for current cell
            hd = m >= e ? (mx << 62 | m_i << 33) : (e1 << 62 | e_i << 33);
            h = m >= e ? m : e;
            hd = h >=f ? hd : f2 << 62;
            h = h >= f ? h : f;
            cur_dp->h = h;

            // e,ed for next cell
            tmp = h - gap_oe;
            e -= gap_e;
            ed = e > tmp ? (e1 << 31 | e_i << 2) : (mx << 31 | m_i << 2);
            e = e > tmp ? e : tmp;
            cur_dp->e = e;

            // f, fd for next cell
            f -= gap_e;
            fd = f > tmp ? f2 : 0;
            f = f > tmp ? f : tmp;
            cur_dp->f = f;

            z[q_i-1] = hd | ed | fd;
        }
        // set TODO_M/E/F for next_line based one cur_line[start] and cur_line[end]
        // set start, end
        if (cur_line[beg].h >= cur_line[end].h) { // move downward
            _next_cen = dp_beg_cen[index_i];
            _beg = beg;
            _end = dp_beg_cen[index_i] + w <= qlen ? dp_beg_cen[index_i] + w : qlen;
            for (i = 0; i < next_n[index_i]; ++i) {
                next_i = next_index[index_i][i];
                next_line = dp_matrix + next_i * matrix_col_n;
                dp_beg_cen[next_i] = dp_beg_cen[next_i] < _next_cen ? dp_beg_cen[next_i] : _next_cen;
                dp_end_cen[next_i] = dp_end_cen[next_i] > _next_cen ? dp_end_cen[next_i] : _next_cen;

                for (j = beg; j <= _end-1; ++j)
                    cur_line[j].todo_map |= HAS_ME;
                cur_line[_end].todo_map |= HAS_E;

                next_line[_beg].todo_map |= TODO_E;
                for (j = _beg+1; j <= _end; ++j) {
                    next_line[j].todo_map |= TODO_MEF;
                }
            }
        } else { // move rightward
            _next_cen = dp_end_cen[index_i] >= qlen ? qlen : dp_end_cen[index_i] + 1;
            _beg = _next_cen >= w ? _next_cen - w : 0;
            _end = end < qlen ? end + 1 : qlen;
            for (i = 0; i < next_n[index_i]; ++i) {
                next_i = next_index[index_i][i];
                next_line = dp_matrix + next_i * matrix_col_n;
                dp_end_cen[next_i] = dp_end_cen[next_i] > _next_cen ? dp_end_cen[next_i] : _next_cen;
                dp_beg_cen[next_i] = dp_beg_cen[next_i] < _next_cen ? dp_beg_cen[next_i] : _next_cen;

                cur_line[_beg-1].todo_map |= HAS_M;
                for (j = _beg; j <= _end-1; ++j) {
                    cur_line[j].todo_map |= HAS_ME;
                }
                if (end == qlen) cur_line[end].todo_map |= HAS_E;

                next_line[_beg].todo_map |= TODO_ME;
                for (j = _beg+1; j < _end; ++j) {
                    next_line[j].todo_map |= TODO_MEF;
                }
                next_line[_end].todo_map |= end == qlen ? TODO_MEF : TODO_MF;
            }
        }
    }

#ifdef __DEBUG__
    for (j = 0; j <= target_node_n; ++j) {
        printf("rank: %d\t", j);
        for (i = 0; i <= qlen; ++i) {
            printf("%d:(%d,%d,%d)\t", i, dp_matrix[j*matrix_col_n+i].h, dp_matrix[j*matrix_col_n+i].e, dp_matrix[j*matrix_col_n+i].f);
        } printf("\n");
    }
#endif
    int in_id, in_index;
    // for global alignment, find best backtrack position
    // TODO semi-global, local
    for (i = 0; i < graph->node[POA_SINK_NODE_ID].in_edge_n; ++i) {
        in_id = graph->node[POA_SINK_NODE_ID].in_id[i];
        in_index = poa_graph_node_id_to_index(graph, in_id);
        _set_max_score(best_score, best_i, best_j, dp_matrix[in_index * matrix_col_n + qlen].h, in_index, qlen);
    }

    printf("best_score: (%d, %d) -> %d\n", best_i, best_j, best_score);
    // backtrack from best score
    if (n_cigar && graph_cigar) {
        int n_c = 0, m_c = 0, id, which;
        int op_shift[4] = {62, 31, 0, 62}, id_shift[4] = {33, 2, 0, 33};
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
    free(dp_matrix); free(qp); free(backtrack_z); free(dp_beg_cen); free(dp_end_cen);
    for (i = 0; i < graph->node_n; ++i) { free(pre_index[i]); } free(pre_index); free(pre_n);
    for (i = 0; i < graph->node_n; ++i) { free(next_index[i]); } free(next_index); free(next_n);
    return best_score;

    return 0;
}

int poa_ada_global_align_sequence_with_graph(poa_graph_t *graph, uint8_t *query, int qlen, poa_para_t *ppt, int *n_cigar, poa_cigar_t **graph_cigar) {
    int8_t *qp, *mat = ppt->mat, gap_o = ppt->gap_open, gap_e = ppt->gap_ext, gap_oe = ppt->gap_open + ppt->gap_ext; // query profile
    int **pre_index, *pre_n, pre_i, **next_index, *next_n;
    int node_id, index_i, q_i;
    int i, j, k, w, *dp_beg, *dp_end;
    int target_node_n = graph->node_n - 2; // exclude start and end nodes
    int matrix_row_n = graph->node_n, matrix_col_n = qlen + 1, z_col_n = qlen; // TODO use less memory? doable, casue band range is known
    int *DP_H, *DP_E, *dp_h, *pre_dp_h, *dp_e, *pre_dp_e, *dp_f;

    uint64_t *backtrack_z, *z; // backtrack matrix; in each cell: hd << 33 | ed << 2 | fd
    uint64_t *hd, *mx, *m_pre_i, *e_pre_i;
    //                                 h<<62|h_id<<33|e<<31|e_i<<2|f
    //                                 h_id/e_i: 29 bit XXX cause error when in_edge_n >= pow(2,29)
    //                                 MATCH:0, DELETION:1, INSERTION:2, MISMATCH:3
    int best_score = INF_32_MIN, best_i=0, best_j=0;

    // allocate memory 
    DP_H = (int*)_err_malloc(matrix_row_n * matrix_col_n * sizeof(int));
    DP_E = (int*)_err_malloc(matrix_row_n * matrix_col_n * sizeof(int));
    dp_f = (int*)_err_malloc(matrix_col_n * sizeof(int));

    dp_beg = (int*)_err_malloc((matrix_row_n) * sizeof(int)); dp_end = (int*)_err_calloc(matrix_row_n, sizeof(int));
    // if w <= 0, do whole global
    // w = ppt->bw <= 0 ? qlen;
    w = 50;
    // calculate band range for each row:
    // have: min_rank, max_rank, min_remain, max_remain
    // then: min_len = min_rank + min_remain, max_len = min_rank + max_remain
    // range: (min_of_two(min_rank, min_rank+qlen-max_len), max_of_two(min_rank+qlen-min_len, max_rank))
    // with w: (min-w, max+w)

    backtrack_z = graph_cigar ? (uint64_t*)_err_malloc(z_col_n * target_node_n * sizeof(uint64_t)) : 0;
    qp = (int8_t*)_err_malloc(qlen * ppt->m);

    hd = (uint64_t*)_err_malloc(qlen * sizeof(int64_t));
    mx = (uint64_t*)_err_malloc(qlen * sizeof(int64_t));
    m_pre_i = (uint64_t*)_err_malloc(matrix_col_n * sizeof(int64_t));
    e_pre_i = (uint64_t*)_err_malloc(matrix_col_n * sizeof(int64_t));

    pre_index = (int**)_err_malloc(graph->node_n * sizeof(int*));
    pre_n = (int*)_err_malloc(graph->node_n * sizeof(int*));
    next_index = (int**)_err_malloc(graph->node_n * sizeof(int*));
    next_n = (int*)_err_malloc(graph->node_n * sizeof(int*));

    // generate the query profile
    for (k = i = 0; k < ppt->m; ++k) {
        const int8_t *p = &mat[k * ppt->m];
        for (j = 0; j < qlen; ++j) qp[i++] = p[query[j]];
    }

    // index of pre-node
    for (i = 0; i < graph->node_n; ++i) {
        node_id = poa_graph_index_to_node_id(graph, i); // i: node index
        pre_n[i] = graph->node[node_id].in_edge_n;
        pre_index[i] = (int*)_err_malloc(pre_n[i] * sizeof(int));
        for (j = 0; j < pre_n[i]; ++j) {
            pre_index[i][j] = poa_graph_node_id_to_index(graph, graph->node[node_id].in_id[j]);
        }
        
        next_n[i] = graph->node[node_id].out_edge_n;
        next_index[i] = (int*)_err_malloc(next_n[i] * sizeof(int));
        for (j = 0; j < next_n[i]; ++j) {
            next_index[i][j] = poa_graph_node_id_to_index(graph, graph->node[node_id].out_id[j]);
        }
    }

    // DP cell: H[i,j], E[i+1,j], F[i,j+1]
    // fill the first row
    dp_beg[0] = GET_DP_BEGIN(graph, w, 0); dp_end[0] = GET_DP_END(graph, w, 0);

    dp_h = DP_H; dp_e = DP_E;
    dp_h[0] = 0; dp_e[0] = -gap_oe;
    for (i = 1; i <= dp_end[0]; ++i) {
        dp_e[i] = INF_32_MIN; dp_h[i] = -gap_o - gap_e * i;
    }

    // DP loop
    int tmp, beg, end, _beg, _end, pre_beg, pre_end;
    uint64_t fd, ed; // f|e|e_id|h|h_id
    uint64_t m0=0x0, e1=0x1, f2=0x2, x3=0x3, he, hf, ee, ff;
                                                     
    he = e1 << HOP_OFF_SET, hf = f2 << HOP_OFF_SET;
    ee = e1 << EOP_OFF_SET; ff = f2 << FOP_OFF_SET;

    for (index_i = 1; index_i < matrix_row_n-1; ++index_i) {
        node_id = poa_graph_index_to_node_id(graph, index_i);

        int8_t *q = &qp[graph->node[node_id].base * qlen];
        dp_h = DP_H + index_i * matrix_col_n; dp_e = DP_E + index_i * matrix_col_n;
        z = &backtrack_z[(index_i-1) * qlen];

        dp_beg[index_i] = GET_DP_BEGIN(graph, w, index_i); dp_end[index_i] = GET_DP_END(graph, w, index_i);

        beg = dp_beg[index_i]; end = dp_end[index_i];
        // init h, e
        for (q_i = beg; q_i <= end; ++q_i) {
            dp_h[q_i] = dp_e[q_i] = INF_32_MIN;
        } 
        
        for (i = 0; i < pre_n[index_i]; ++i) {
            pre_i = pre_index[index_i][i];
            pre_dp_h = DP_H + pre_i * matrix_col_n; pre_dp_e = DP_E + pre_i * matrix_col_n;
            pre_beg = dp_beg[pre_i]; pre_end = dp_end[pre_i];
            // set M from (pre_i, q_i-1)
            _beg = MAX_OF_TWO(beg-1, pre_beg), _end = MIN_OF_TWO(end-1, pre_end);
            for (q_i = _beg; q_i <= _end; ++q_i) { // SIMD parallelization
                m_pre_i[q_i+1] = pre_dp_h[q_i] > dp_h[q_i+1] ? pre_i : m_pre_i[q_i+1];
                dp_h[q_i+1] = MAX_OF_TWO(pre_dp_h[q_i], dp_h[q_i+1]);
            }
            _beg = MAX_OF_TWO(beg, pre_beg), _end = MIN_OF_TWO(end, pre_end);
            // set E from (pre_i, q_i)
            for (q_i = _beg; q_i <= _end; ++q_i) { // SIMD parallelization
                e_pre_i[q_i] = pre_dp_e[q_i] > dp_e[q_i] ? pre_i : e_pre_i[q_i];
                dp_e[q_i] = MAX_OF_TWO(pre_dp_e[q_i], dp_e[q_i]);
            }
        }
        // compare M, E, and F
        for (q_i = beg; q_i <= end; ++q_i) { // SIMD parallelization
            // get M score
            if (q_i) {
                dp_h[q_i] += q[q_i-1]; mx[q_i-1] = q[q_i-1] == ppt->match ? m0 : x3;
                // h,hd for current cell
                hd[q_i-1] = dp_h[q_i] >= dp_e[q_i] ? (mx[q_i-1] << HOP_OFF_SET | m_pre_i[q_i] << HID_OFF_SET) : (he | e_pre_i[q_i] << HID_OFF_SET);
            }
            dp_h[q_i] = MAX_OF_TWO(dp_h[q_i], dp_e[q_i]);
        }

        // set F from (index_i, q_i-1)
        dp_f[beg] = INF_32_MIN;
        for (q_i = beg+1; q_i <= end; ++q_i) { // XXX no SIMD parallelization
            dp_f[q_i] = MAX_OF_TWO(dp_h[q_i-1] - gap_oe, dp_f[q_i-1] - gap_e);
            // fd = h - oe > f - e ? 0 : f2
        }
            
        // since we have score of M, E and F, now we need to set H for cur cell, set E and F for next cell
        for (q_i = beg; q_i <= end; ++q_i) { // SIMD parallelization
            // h,hd for current cell
            if (q_i) hd[q_i-1] = dp_h[q_i] >= dp_f[q_i] ? hd[q_i-1] : hf;
            dp_h[q_i] = MAX_OF_TWO(dp_h[q_i], dp_f[q_i]);

            // e,ed for next cell
            tmp = dp_h[q_i] - gap_oe;
            dp_e[q_i] -= gap_e;
            if (q_i) ed = dp_e[q_i] > tmp ? (ee | e_pre_i[q_i] << EID_OFF_SET) : (mx[q_i-1] << EOP_OFF_SET | m_pre_i[q_i] << EID_OFF_SET);
            dp_e[q_i] = MAX_OF_TWO(dp_e[q_i], tmp);

            // fd for next cell
            if (q_i) {
                dp_f[q_i] -= gap_e;
                fd = dp_f[q_i] > tmp ? ff : (mx[q_i-1] << FOP_OFF_SET);

                z[q_i-1] = hd[q_i-1] | ed | fd;
            }
        }
    }


#ifdef __DEBUG__
    for (j = 0; j <= target_node_n; ++j) {
        printf("index: %d\t", j);
        //for (i = 0; i <= qlen; ++i) {
        for (i = dp_beg[j]; i <= dp_end[j]; ++i) {
            printf("%d:(%d,%d)\t", i, DP_H[j*matrix_col_n+i], DP_E[j*matrix_col_n+i]);
        } printf("\n");
    }
#endif
    int in_id, in_index;
    // for global alignment, find best backtrack position
    // TODO semi-global, local
    for (i = 0; i < graph->node[POA_SINK_NODE_ID].in_edge_n; ++i) {
        in_id = graph->node[POA_SINK_NODE_ID].in_id[i];
        in_index = poa_graph_node_id_to_index(graph, in_id);
        _set_max_score(best_score, best_i, best_j, DP_H[in_index * matrix_col_n + qlen], in_index, qlen);
    }

    printf("best_score: (%d, %d) -> %d\n", best_i, best_j, best_score);
    // backtrack from best score
    poa_backtrack(backtrack_z, best_i, best_j, z_col_n, graph, n_cigar, graph_cigar);
    
    free(DP_H); free(DP_E); free(dp_f); free(m_pre_i); free(e_pre_i); free(hd); free(mx);
    free(qp); free(backtrack_z);
    free(dp_beg); free(dp_end);
    for (i = 0; i < graph->node_n; ++i) { 
        free(pre_index[i]); free(next_index[i]); 
    } 
    free(pre_index); free(pre_n); free(next_index); free(next_n);
    return best_score;

    return 0;
}

// TODO linear gap penalty: gap_o == 0
// TODO W ==> total error in the alignment
// TODO w ==> max consecutive error
/*
int poa_ada_forefront_global_align_sequence_with_graph(poa_graph_t *graph, uint8_t *query, int qlen, poa_para_t *ppt, int *n_cigar, poa_cigar_t **graph_cigar) {
    dp_matrix_t *dp_matrix, **forefront_dp_start, *dp, *cell; // Full: (tlen + 1) * (qlen + 1); Banded: (tlen+1) * (2 * w + 1)
    int8_t *qp, *mat = ppt->mat, gap_o = ppt->gap_open, gap_e = ppt->gap_ext, gap_oe = ppt->gap_open + ppt->gap_ext; // query profile
    int *forefront_start, *forefront_end, *dp_matrix_start, forefront_n, **pre_index, *pre_n, pre_i;
    int rank_i, node_id, index_i, q_i, next_index_i;
    int i, j, k, w;
    int target_node_n = graph->node_n - 2; // exclude start and end nodes
    int matrix_row_n = target_node_n + 1, z_col_n = qlen; // TODO use less memory ???
    uint64_t *backtrack_z; // backtrack matrix; in each cell: hd << 33 | ed << 2 | fd
    //                                 h<<62|h_id<<33|e<<31|e_i<<2|f
    //                                 h_id/e_i: 29 bit XXX cause error when in_edge_n >= pow(2,29)
    //                                 MATCH:0, DELETION:1, INSERTION:2, MISMATCH:3
    int best_score = INF_32_MIN, best_i=0, best_j=0;

    // allocate memory 
    forefront_n = target_node_n + qlen + 1 + 2; // make sure no invalid write
    dp_matrix = (dp_matrix_t*)_err_calloc((matrix_row_n + 1) * forefront_n, sizeof(dp_matrix_t)); // make sure no invalid write 
    // TODO if no band is used
    // w = ppt->bw <= 0 ? (MIN_OF_TWO(graph->rank_n-2, qlen)) : (MIN_OF_THREE(graph->rank_n-2, qlen, ppt->bw));
    w = 10; // XXX small may cause NO alignment result XXX
    backtrack_z = graph_cigar ? (uint64_t*)_err_malloc(z_col_n * target_node_n * sizeof(uint64_t)) : 0;

    qp = (int8_t*)_err_malloc(qlen * ppt->m); // TODO score profile for diagnal

    forefront_start = (int*)_err_malloc(forefront_n * sizeof(int)); // [i-start, start]
    forefront_end = (int*)_err_malloc(forefront_n * sizeof(int));   // [end, i-end]
    dp_matrix_start = (int*)_err_malloc(forefront_n * sizeof(int));
    forefront_dp_start = (dp_matrix_t**)_err_malloc(forefront_n * sizeof(dp_matrix_t*));


    pre_index = (int**)_err_malloc(graph->node_n * sizeof(int*));
    pre_n = (int*)_err_malloc(graph->node_n * sizeof(int*));

    // generate the query profile
    for (k = i = 0; k < ppt->m; ++k) {
        const int8_t *p = &mat[k * ppt->m];
        for (j = 0; j < qlen; ++j) qp[i++] = p[query[j]];
    }

    // index of pre-node
    for (i = 0; i < graph->node_n; ++i) {
        node_id = poa_graph_index_to_node_id(graph, i); // i: node index
        pre_n[i] = graph->node[node_id].in_edge_n;
        pre_index[i] = (int*)_err_malloc(pre_n[i] * sizeof(int));
        for (j = 0; j < pre_n[i]; ++j) {
            pre_index[i][j] = poa_graph_node_id_to_index(graph, graph->node[node_id].in_id[j]);
        }
    }

    // initialize edge for forefront (index)
    for (i = 2; i <= qlen; ++i) forefront_start[i] = 1; 
    for (i = 2; i <= target_node_n; ++i) forefront_end[i] = i-1;
    for (i = qlen+1; i < forefront_n; ++i) forefront_start[i] = i - qlen;
    for (i = target_node_n+1; i < forefront_n; ++i) forefront_end[i] = target_node_n;

    dp_matrix_start[0] = 0;
    for (i = 1; i <= target_node_n+2; ++i) dp_matrix_start[i] = i;
    for (i = target_node_n+3; i < forefront_n; ++i) dp_matrix_start[i] = target_node_n + 2;
    for (i = 1; i < forefront_n; ++i) dp_matrix_start[i] += dp_matrix_start[i-1];
    for (i = 0; i < forefront_n; ++i)  {
#ifdef __DEBUG__
        err_printf("forefront start: %d => %d\n", i, dp_matrix_start[i]);
#endif
        forefront_dp_start[i] = dp_matrix + dp_matrix_start[i];
    }

    // cur_cell: H[i,j], E[i+1,j], F[i,j+1]
    // fill the first row
    cell = dp_matrix;
    cell->h = 0; cell->e = -gap_oe; cell->f = -gap_oe;
    cell->xl = 0; cell->fl = cell->el = 1; cell->todo_map = HAS_MX;
    dp_matrix[4].todo_map = TODO_MX;
    for (i = 1; i <= w; ++i) {
        cell = forefront_dp_start[i];
        cell->e = INF_32_MIN; cell->h = -gap_o - gap_e * i; cell->f = cell->h - gap_e;
        cell->fl = i+1; cell->el = cell->xl = 0; 
        cell->todo_map = HAS_MX; 
        (forefront_dp_start[i+2]+1)->todo_map = TODO_MX;
    }
    // fill the first column
    for (i = 1; i <= w; ++i) { // XXX no need to cal min_rank for all nodes ???
        for (j = graph->min_rank_index[i-1]; j < graph->min_rank_index[i]; ++j) {
            index_i = graph->min_rank_to_index[i];
            node_id = poa_graph_index_to_node_id(graph, index_i);
            cell = forefront_dp_start[i] + i;
            cell->f = INF_32_MIN; cell->h = -gap_o - gap_e * i; cell->e = cell->h - gap_e;
            cell->xl = cell->fl = 0; cell->el = i+1;

            cell->todo_map = HAS_MX;
            for (k = 0; k < graph->node[node_id].out_edge_n; ++k) {
                next_index_i = poa_graph_node_id_to_index(graph, graph->node[node_id].out_id[k]);
                (forefront_dp_start[next_index_i + 1] + next_index_i)->todo_map = TODO_MX;
            }
        }
    }

    // DP loop in anti-diagnal direction
    int32_t h, m, s, xl, e, el, f, fl, start, end, sum;
    uint64_t m_i, fd, ed, e_i, hd, mx, m0=0x0, e1=0x1, f2=0x2, x3=0x3; // h<<62|h_id<<33|e<<31|e_i<<2|f
    dp_matrix_t *cur_dp, *m_pre_dp, *e_pre_dp, *f_pre_dp;

    for (sum = 2; sum < forefront_n; ++sum) {
        start = forefront_start[sum], end = forefront_end[sum];
        dp = forefront_dp_start[sum];
        for (index_i = start; index_i <= end; ++index_i) {
            q_i = sum - index_i;
            cur_dp = dp + index_i;
            if (!cur_dp->todo_map) continue;

            node_id = poa_graph_index_to_node_id(graph, index_i);
            rank_i = poa_graph_index_to_rank(graph, index_i);

            m = e = f = INF_32_MIN;
            m_i = e_i = -1;

            mx = graph->node[node_id].base == query[q_i-1] ? m0 : x3; // XXX score profile for diagnal
            s = (mx == m0 ? ppt->match : -ppt->mismatch);
            xl = fl = el = 0;
            
            if (mx == m0) {
                if (cur_dp->todo_map & TODO_M) { // M: from (pre_i, q_i-1)
                    for (j = 0; j < pre_n[index_i]; ++j) {
                        pre_i = pre_index[index_i][j];
                        m_pre_dp = forefront_dp_start[pre_i + q_i - 1] + pre_i;
                        if (m_pre_dp->todo_map & HAS_M) {
                            if (m_pre_dp->h > m) {
                                m = m_pre_dp->h; m_i = pre_i; xl = m_pre_dp->xl;
                            }
                        }
                    }
                }
            } else { // mx = x3
                if (cur_dp->todo_map & TODO_X) {
                    for (j = 0; j < pre_n[index_i]; ++j) {
                        pre_i = pre_index[index_i][j];
                        m_pre_dp = forefront_dp_start[pre_i + q_i - 1] + pre_i;
                        if (m_pre_dp->todo_map & HAS_X) {
                            if (m_pre_dp->h > m) {
                                m = m_pre_dp->h; m_i = pre_i; xl = m_pre_dp->xl;
                            }
                        }

                    }
                }
            }
            
            if (cur_dp->todo_map & TODO_E) { // E: from (pre_i, q_i)
                for (j = 0; j < pre_n[index_i]; ++j) {
                    pre_i = pre_index[index_i][j];
                    e_pre_dp = forefront_dp_start[pre_i + q_i] + pre_i;
                    if (e_pre_dp->todo_map & HAS_E) {
                        if (e_pre_dp->e > e) {
                            e = e_pre_dp->e; e_i = pre_i; el = e_pre_dp->el;
                        }
                    }
                }
            }
            if (cur_dp->todo_map & TODO_F) { // F: from (index_i, q_i-1)
                f_pre_dp = forefront_dp_start[index_i + q_i - 1] + index_i; // (f_pre_dp->todo_map & HAS_F): always HAS_F
                f = f_pre_dp->f;
                fl = f_pre_dp->fl;
            }
            // since we have score of M, E and F, now we need to set H for cur cell, set E and F for next cell

            // h,hd,xl for current cell
            // TODO z-drop
            m += s;
            hd = m >= e ? (mx << 62 | m_i << 33) : (e1 << 62 | e_i << 33);
            h = m >= e ? m : e;
            hd = h >=f ? hd : f2 << 62;
            h = h >= f ? h : f;
            cur_dp->h = h;
            cur_dp->xl = mx == m0 ? 0 : xl + 1;
            if (cur_dp->xl <= w) {
                cur_dp->todo_map |= HAS_MX;
                for (i = 0; i < graph->node[node_id].out_edge_n; ++i) { // next M: (next_index_i, q_i+1)
                    next_index_i = poa_graph_node_id_to_index(graph, graph->node[node_id].out_id[i]);
#ifdef __DEBUG__
                    err_printf("%d, %d\n", next_index_i, q_i);
#endif
                    (forefront_dp_start[next_index_i + q_i + 1] + next_index_i)->todo_map |= TODO_MX;
                }
            } else {
                cur_dp->todo_map |= HAS_M;
                for (i = 0; i < graph->node[node_id].out_edge_n; ++i) { // next M: (next_index_i, q_i+1)
                    next_index_i = poa_graph_node_id_to_index(graph, graph->node[node_id].out_id[i]);
                    (forefront_dp_start[next_index_i + q_i + 1] + next_index_i)->todo_map |= TODO_M;
                }
            }

            // e,ed,el and f,fd,fl for next cell
            if (e - gap_e > m - gap_oe) {
                ed = (e1 << 31 | e_i << 2);
                cur_dp->e = e - gap_e;
                cur_dp->el = el + 1;
                if (cur_dp->el <= w) {
                    cur_dp->todo_map |= HAS_E;
                    for (i = 0; i < graph->node[node_id].out_edge_n; ++i) { // next E: (next_index_i, q_i)
                        next_index_i = poa_graph_node_id_to_index(graph, graph->node[node_id].out_id[i]);
                        (forefront_dp_start[next_index_i + q_i] + next_index_i)->todo_map |= TODO_E;
                    }
                }
            } else {
                ed = (mx << 31 | m_i << 2);
                cur_dp->e = m - gap_oe;
                cur_dp->el = 1;
                cur_dp->todo_map |= HAS_E;
                for (i = 0; i < graph->node[node_id].out_edge_n; ++i) { // next E: (next_index_i, q_i)
                    next_index_i = poa_graph_node_id_to_index(graph, graph->node[node_id].out_id[i]);
                    (forefront_dp_start[next_index_i + q_i] + next_index_i)->todo_map |= TODO_E;
                }
            }
            if (f - gap_e > m - gap_oe) { // cur_dp->todo_map |= HAS_F;
                fd = f2;
                cur_dp->f = f - gap_e;
                cur_dp->fl = fl + 1;
                if (cur_dp->fl <= w) // next F: (index_i, q_i+1)
                    (forefront_dp_start[index_i + q_i + 1] + index_i)->todo_map |= TODO_F;
            } else {
                fd = 0;
                cur_dp->f = m - gap_oe;
                cur_dp->fl = 1;
                (forefront_dp_start[index_i + q_i + 1] + index_i)->todo_map |= TODO_F;
            }
            backtrack_z[(index_i-1) * qlen + q_i-1] = hd | ed | fd;
        }
        // printf("Sum %d ==> Todo %d\n", sum, todo_i);
    }
#ifdef __DEBUG__
    for (j = 0; j <= target_node_n; ++j) {
        printf("rank: %d\t", j);
        for (i = 0; i <= qlen; ++i) {
            printf("%d:(%d,%d,%d)\t", i, (forefront_dp_start[j + i]+j)->h, (forefront_dp_start[j + i]+j)->e, (forefront_dp_start[j + i]+j)->f);
        } printf("\n");
    }
#endif
    int in_id, in_index;
    // for global alignment, find best backtrack position
    // TODO semi-global, local
    for (i = 0; i < graph->node[POA_SINK_NODE_ID].in_edge_n; ++i) {
        in_id = graph->node[POA_SINK_NODE_ID].in_id[i];
        in_index = poa_graph_node_id_to_index(graph, in_id);
        _set_max_score(best_score, best_i, best_j, (forefront_dp_start[in_index + qlen] + in_index)->h, in_index, qlen);
    }

    printf("best_score: (%d, %d) -> %d\n", best_i, best_j, best_score);
    // backtrack from best score
    if (n_cigar && graph_cigar) {
        int n_c = 0, m_c = 0, id, which;
        int op_shift[4] = {62, 31, 0, 62}, id_shift[4] = {33, 2, 0, 33};
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
    free(dp_matrix); free(dp_matrix_start); free(forefront_dp_start); 
    free(forefront_start); free(forefront_end);
    free(qp); free(backtrack_z);
    for (i = 0; i < graph->node_n; ++i) { free(pre_index[i]); } free(pre_index); free(pre_n);
    return best_score;

    return 0;
}*/
