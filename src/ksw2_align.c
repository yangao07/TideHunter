#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "utils.h"
#include "ksw2.h"

#ifdef __cplusplus
extern "C" {
#endif

int match = 1, mis = 2, gap_open = 2, gap_ext = 1;
int8_t mat[25] = {
     1, -2, -2, -2, 0, 
    -2,  1, -2, -2, 0, 
    -2, -2,  1, -2, 0, 
    -2, -2, -2,  1, 0, 
     0,  0,  0,  0, 0};

static void ksw_gen_simple_mat(int m, int8_t *mat, int8_t a, int8_t b)
{
	int i, j;
	a = a < 0? -a : a;
	b = b > 0? -b : b;
	for (i = 0; i < m - 1; ++i) {
		for (j = 0; j < m - 1; ++j)
			mat[i * m + j] = i == j? a : b;
		mat[i * m + m - 1] = 0;
	}
	for (j = 0; j < m; ++j)
		mat[(m - 1) * m + j] = 0;
}

static void print_seq(const uint8_t *query, int qlen, const uint8_t *target, int tlen) {
    int i;
    printf("Query:  "); for (i = 0; i <qlen; ++i) putchar("ACGT"[query[i]]); putchar('\n');
    printf("Target: "); for (i = 0; i <tlen; ++i) putchar("ACGT"[target[i]]); putchar('\n');
}

static void print_aln(const char *tname, const char *qname, ksw_extz_t *ez)
{
	printf("%s\t%s\t%d", tname, qname, ez->score);
	// printf("\t%d\t%d\t%d", ez->max, ez->max_t, ez->max_q);
	if (ez->n_cigar > 0) {
		int i;
		putchar('\t');
		for (i = 0; i < ez->n_cigar; ++i)
			printf("%d%c", ez->cigar[i]>>4, "MID"[ez->cigar[i]&0xf]);
	}
	putchar('\n');
}

static void print_cigar(int n_cigar, uint32_t *cigar) {
    if (n_cigar > 0) {
        int i;
        printf("cigar:\t");
        for (i = 0; i < n_cigar; ++i)
            printf("%d%c", cigar[i]>>4, "MID"[cigar[i]&0xf]);
        putchar('\n');
    }
}

int *ksw2_get_xid(uint32_t *cigar, int n_cigar, const uint8_t *query, const uint8_t *target) {
    int i, j, *xid = (int*)_err_calloc(4, sizeof(int));
    int qi = 0, ti = 0, op, len;
    for (i = 0; i < n_cigar; ++i) {
        op = cigar[i] & 0xf; len = cigar[i] >> 4;
        if (op == 0) { // MATCH/MISMATCH
            for (j = 0; j < len; ++j) {
                if (query[qi+j] == target[ti+j]) xid[0]++;
                else xid[3]++;
            }
            qi += len; ti += len;
        } else if (op == 1) { // INSERTION
            xid[1] += len;
            qi += len;
        } else if (op == 2) { // DELETION
            xid[2] += len;
            ti += len;
        } else {
            err_fatal_core(__func__, "Unexpected cigar op: %d\n", op);
            free(xid);
            return 0;
        }
    }
    return xid;
}

int backtrack_left_end(int n_cigar, uint32_t *cigar, int qlen, int tlen, int q_left_ext) {
    int t_left_ext = 0, i, j, op, len;
    int q_remain_len = q_left_ext;
    for (i = n_cigar-1; i >=0; --i) {
        op = cigar[i] & 0xf; len = cigar[i] >> 4;
        if (op == 0) { // MATCH/MISMATCH
            if (len >= q_remain_len) {
                t_left_ext += q_remain_len;
                return t_left_ext;
            } else {
                t_left_ext += len;
                q_remain_len -= len;
            }
        } else if (op == 1) { // INSERTION
            if (len >= q_remain_len) {
                return t_left_ext;
            } else {
                q_remain_len -= len;
            }
        } else if (op == 2) { // DELETION
            t_left_ext += len;
        }
    }
    if (q_remain_len > 0) {
        err_fatal_simple("Error: unmatched cigar and q_left_ext.\n");
    }
    return t_left_ext;
}

int ksw2_global(const uint8_t *query, int qlen, const uint8_t *target, int tlen) {
    ksw_extz_t ez; memset(&ez, 0, sizeof(ksw_extz_t));
    int w=-1, zdrop=-1, end_bonus=0, flag = 0;
    ksw_extz2_sse(0, qlen, query, tlen, target, 5, mat, gap_open, gap_ext, w, zdrop, end_bonus, flag, &ez);
    #ifdef __DEBUG__
    print_cigar(ez.n_cigar, ez.cigar);
    #endif
    int *xid = ksw2_get_xid(ez.cigar, ez.n_cigar, query, target);
    int iden_n = 0;
    if (xid != 0) {
        iden_n = xid[0];
        free(xid);
    }
    if (ez.cigar) free(ez.cigar);
    return iden_n; 
}

// TODO use extz2_sse
int ksw2_global_with_cigar(const uint8_t *query, int qlen, const uint8_t *target, int tlen, int *n_cigar, uint32_t **cigar) {
    ksw_extz_t ez; memset(&ez, 0, sizeof(ksw_extz_t));
    int w=-1, zdrop=-1, end_bonus=0, flag = 0;
    ksw_extz2_sse(0, qlen, query, tlen, target, 5, mat, gap_open, gap_ext, w, zdrop, end_bonus, flag, &ez);
    #ifdef __DEBUG__
    print_cigar(ez.n_cigar, ez.cigar);
    #endif
    int *xid = ksw2_get_xid(ez.cigar, ez.n_cigar, query, target);
    int iden_n = 0;
    if (xid != 0) {
        iden_n = xid[0];
        free(xid);
    }
    *n_cigar = ez.n_cigar;
    *cigar = ez.cigar;
    // if (ez.cigar) free(ez.cigar);
    return iden_n; 
}

void ksw2_right_ext(const uint8_t *query, int qlen, const uint8_t *target, int tlen, int *max_q, int *max_t) {
    int iden_n = 0; void *km = 0;
    ksw_extz_t ez; memset(&ez, 0, sizeof(ksw_extz_t));
    int w=-1, zdrop=-1, flag = KSW_EZ_EXTZ_ONLY | KSW_EZ_SCORE_ONLY;
    ksw_extz2_sse(km, qlen, query, tlen, target, 5, mat, gap_open, gap_ext, w, zdrop, 0, flag, &ez);
    *max_q = ez.max_q, *max_t = ez.max_t;
}

void ksw2_left_ext(const uint8_t *query, int qlen, const uint8_t *target, int tlen, int *max_q, int *max_t) {
    uint8_t *rquery = (uint8_t*)_err_malloc(qlen * sizeof(uint8_t));
    uint8_t *rtarget = (uint8_t*)_err_malloc(tlen * sizeof(uint8_t));
    int i, iden_n = 0; void *km = 0;
    for (i = 0; i < qlen; ++i) rquery[i] = query[qlen-i-1];
    for (i = 0; i < tlen; ++i) rtarget[i] = target[tlen-i-1];

    ksw_extz_t ez; memset(&ez, 0, sizeof(ksw_extz_t));
    int w=-1, zdrop=-1, flag = KSW_EZ_EXTZ_ONLY | KSW_EZ_SCORE_ONLY;
    ksw_extz2_sse(km, qlen, rquery, tlen, rtarget, 5, mat, gap_open, gap_ext, w, zdrop, 0, flag, &ez);
    *max_q = ez.max_q, *max_t = ez.max_t;
    free(rquery); free(rtarget);
}

// return target_right_ext_end and query_right_ext_end
//                |: target_end
//|-----target-------
//|-----query-----|
int ksw2_right_extend(const uint8_t *query, int qlen, const uint8_t *target, int tlen, int *mqe_tlen) {
    ksw_extz_t ez; memset(&ez, 0, sizeof(ksw_extz_t));
    int w=-1, zdrop=-1, end_bonus=0, flag = KSW_EZ_RIGHT | KSW_EZ_EXTZ_ONLY;
    ksw_extz2_sse(0, qlen, query, tlen, target, 5, mat, gap_open, gap_ext, w, zdrop, end_bonus, flag, &ez);
    #ifdef __DEBUG__
    printf("max_q: %d, max_t: %d, qlen: %d, tlen: %d\n", ez.max_q, ez.max_t, qlen, tlen);
    print_cigar(ez.n_cigar, ez.cigar);
    #endif
    int *xid = ksw2_get_xid(ez.cigar, ez.n_cigar, query, target);
    int iden_n = 0;
    if (xid != 0) {
        iden_n = xid[0];
        free(xid);
    }
    *mqe_tlen = ez.mqe_t + 1;
    if (ez.cigar) free(ez.cigar);
    return iden_n;
}

// return target_left_ext_end and query_left_ext_end
//  |: target_end
//-------target-----|
//  |-----query-----|
int ksw2_left_extend(const uint8_t *query, int qlen, const uint8_t *target, int tlen, int *mqe_tlen) {
    uint8_t *rquery = (uint8_t*)_err_malloc(qlen * sizeof(uint8_t));
    uint8_t *rtarget = (uint8_t*)_err_malloc(tlen * sizeof(uint8_t));
    int i; ksw_extz_t ez; memset(&ez, 0, sizeof(ksw_extz_t));
    for (i = 0; i < qlen; ++i) rquery[i] = query[qlen-i-1];
    for (i = 0; i < tlen; ++i) rtarget[i] = target[tlen-i-1];
    int w=-1, zdrop=-1, end_bonus=0, flag=KSW_EZ_RIGHT | KSW_EZ_EXTZ_ONLY;
    ksw_extz2_sse(0, qlen, rquery, tlen, rtarget, 5, mat, gap_open, gap_ext, w, zdrop, end_bonus, flag, &ez);
    #ifdef __DEBUG__
    printf("max_q: %d, max_t: %d, qlen: %d, tlen: %d\n", ez.max_q, ez.max_t, qlen, tlen);
    print_cigar(ez.n_cigar, ez.cigar);
    #endif
    int *xid = ksw2_get_xid(ez.cigar, ez.n_cigar, rquery, rtarget);
    int iden_n = 0;
    if (xid != 0) {
        iden_n = xid[0];
        free(xid);
    }
    *mqe_tlen = ez.mqe_t + 1;
    if (ez.cigar) free(ez.cigar);
    free(rquery); free(rtarget);
    return iden_n;
}

#ifdef __cplusplus
}
#endif
