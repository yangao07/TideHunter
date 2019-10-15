#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "tandem_chain.h"
#include "edlib_align.h"
#include "abpoa_cons.h"
#include "ksw2_align.h"
#include "utils.h"

// pos: 0-base
int *partition_seqs_core(char *seq, int seq_len, int8_t *hit_array, int est_period, int chain_start, int chain_end, int *par_n) {
    int i, j;
    int *pos_array = (int*)_err_malloc(sizeof(int) * seq_len), *par_pos = (int*)_err_malloc(sizeof(int) * seq_len);
    int hit_n = 0;
    for (i = 0; i < seq_len; ++i) {
        if (hit_array[i]) {
            #ifdef __DEBUG__
            printf("%d, ", i);
            #endif
            pos_array[hit_n++] = i;
        }
    } 
#ifdef __DEBUG__
    printf("\n");
#endif
    // partition seq into period seperated seqs
    int par_i = 0, l, tot_len; 
    char *query_seq, *target_seq; int copy_num, ed, start, end, target_start;
    //  extend two ends
    l = est_period / 2;
    copy_num = (int)((double)(pos_array[0] - chain_start) / est_period + 0.5);
    if (copy_num >= 1) {
        for (j = copy_num; j >= 1; --j) {
            query_seq = seq + pos_array[0];
            target_start = pos_array[0] - est_period * j - (l<<0);
            // target_start = pos_array[0] - period[pos_array[0]] * j - (l<<1);
            if (target_start < 0) continue;
            target_seq = seq + target_start;
            ed = edlib_align_HW(query_seq, l, target_seq, l<<1, &start, &end);
            // ed = edlib_align_HW(query_seq, l, target_seq, l<<2, &start, &end);
            if (ed >= 0) par_pos[par_i++] = target_start + start;
        }
    }
    par_pos[par_i++] = pos_array[0];
    for (i = 0; i < hit_n-1; ++i) {
        // printf("%d: %d, %d, %d\n", pos_array[i], period[pos_array[i]], period[pos_array[i+1]], ave_p);
        copy_num = (int)((double)(pos_array[i+1] - pos_array[i]) / est_period + 0.5);
        if (copy_num > 1) { // multiple copies: semi-global alignment of prefix l-mer using edlib
            tot_len = pos_array[i+1] - pos_array[i];
            query_seq = seq + pos_array[i];
            for (j = 1; j < copy_num; ++j) {
                target_start = pos_array[i] + est_period * j - (l << 0);
                // target_start = pos_array[i] + ave_p * j - (l << 1);
                target_seq = seq + target_start;
                ed = edlib_align_HW(query_seq, l, target_seq, l<<1, &start, &end);
                // ed = edlib_align_HW(query_seq, l, target_seq, l<<2, &start, &end);
                if (ed < 0) { // no alignment result
                    par_pos[par_i++] = -1; // skip this copy
                } else {
                    par_pos[par_i++] = target_start + start;
                }
            }
        }
        par_pos[par_i++] = pos_array[i+1];
    }
    // extend par_pos on the right
    copy_num = (int)((double)(chain_end - pos_array[hit_n-1]) / est_period + 0.5);
    if (copy_num >= 1) {
        for (j = 1; j <= copy_num; ++j) {
            query_seq = seq + pos_array[hit_n-1];
            target_start = pos_array[hit_n-1] + est_period * j - (l<<0);
            // target_start = pos_array[hit_n-1] + period[pos_array[hit_n-1]] * j - (l<<1);
            if (target_start + (l<<1) >= seq_len) continue;
            // if (target_start + (l<<2) >= seq_len) continue;
            target_seq = seq + target_start;
            ed = edlib_align_HW(query_seq, l, target_seq, l<<1, &start, &end);
            // ed = edlib_align_HW(query_seq, l, target_seq, l<<2, &start, &end);
            if (ed >= 0) par_pos[par_i++] = target_start + start;
        }
    }
    free(pos_array);
    *par_n = par_i; 
    return par_pos;
}

// TODO when (s2) >> (e1), i.e., anchor are too sparse
// start with est_start anchor?
int *get_partition_pos_with_global_alignment(uint8_t *bseq, int seq_len, dp_t **dp, chain_t ch, mini_tandem_para *mtp, int *par_n) {
    int est_ch_i = ch.est_ch_i, est_start = ch.est_start, est_period = ch.est_period;
    int first_end = dp[ch.cell[0].i][ch.cell[0].j].end, last_start = dp[ch.cell[ch.len-1].i][ch.cell[ch.len-1].j].start;
    int *par_pos = (int*)_err_malloc(seq_len * sizeof(int)); *par_n = 0;
    int i, j, k=mtp->k, ch_i1, ch_i2, s1, e1, s2, e2, iden_n, s, e;
    int n_cigar; uint32_t *cigar;
    cell_t c; dp_t d;
    // find non-overlapping adjacent tandem repeat hits
    // global alignment of [s1,s2] and [e1, e2]
    // find best partition position with backtrack
    ch_i2 = est_ch_i, s2 = est_start, e2 = est_start + est_period;
    while (ch_i2 > 0 && s2 >= first_end) {
        for (i = ch_i2 - 1; i >= 0; --i) {
            c = ch.cell[i]; d = dp[c.i][c.j];
            s1 = d.start; e1 = d.end;
            if (e1 == s2) {
                par_pos[(*par_n)++] = s1;
                ch_i2 = i, s2 = s1, e2 = e1;
                break;
            } else if (e1 < s2) { // non-overlapping adjacent tandem repeat hits
                // do global alignment
                iden_n = ksw2_global_with_cigar(bseq+e1-k+1, e2-e1+k, bseq+s1-k+1, s2-s1+k, &n_cigar, &cigar);
                #ifdef __DEBUG__
                printf("iden_n: %d (%d,%d), (%d,%d)\n", iden_n, e2-e1+k, s2-s1+k, s1, s2);
                #endif
                if (iden_n >= MIN_OF_TWO(s2-s1+k, e2-e1+k) * (1-mtp->max_div)) { // extend partition
                    s = s2 - backtrack_left_end(n_cigar, cigar, e2-e1+k, s2-s1+k, e2-s2);
                    par_pos[(*par_n)++] = s; 
                    ch_i2 = i+1; e2 = s2, s2 = s;
                } else { // skip this anchor
                    par_pos[(*par_n)++] = -1; // indicate a separation flag
                    par_pos[(*par_n)++] = e1; 
                    par_pos[(*par_n)++] = s1; 
                    ch_i2 = i; s2 = s1; e2 = e1;
                }
                if (cigar) 
                    free(cigar);
                break;
            }
        }
    }
    // reverse par_pos
    for (i = 0; i < (*par_n)>> 1; ++i) {
        j = par_pos[i]; par_pos[i] = par_pos[*par_n-i-1]; par_pos[*par_n-i-1] = j;
    }
    par_pos[(*par_n)++] = est_start;
    par_pos[(*par_n)++] = est_start + est_period;

    ch_i1 = est_ch_i, s1 = est_start, e1 = est_start + est_period;
    while (ch_i1 < ch.len-1 && e1 <= last_start) {
        for (i = ch_i1+1; i < ch.len; ++i) {
            c = ch.cell[i]; d = dp[c.i][c.j];
            s2 = d.start; e2 = d.end;
            if (s2 == e1) {
                par_pos[(*par_n)++] = e2;
                ch_i1 = i, s1 = s2, e1 = e2;
                break;
            } else if (s2 > e1) { // non-overlapping adjacent tandem repeat hits
                // do global alignment
                iden_n = ksw2_global_with_cigar(bseq+s1-k+1, s2-s1+k, bseq+e1-k+1, e2-e1+k, &n_cigar, &cigar);
                #ifdef __DEBUG__
                printf("iden_n: %d (%d,%d), (%d,%d)\n", iden_n, e2-e1+k, s2-s1+k, s1, s2);
                #endif
                if (iden_n >= MIN_OF_TWO(s2-s1+k, e2-e1+k) * (1-mtp->max_div)) {
                    e = e2 - backtrack_left_end(n_cigar, cigar, s2-s1+k, e2-e1+k, s2-e1);
                    par_pos[(*par_n)++] = e; 
                    ch_i1 = i-1; s1 = e1, e1 = e;
                } else {
                    par_pos[(*par_n)++] = -1; 
                    par_pos[(*par_n)++] = s2; 
                    par_pos[(*par_n)++] = e2; 
                    ch_i1 = i; s1 = s2; e1 = e2;
                }
                if (cigar) 
                    free(cigar);
                break;
            }
        }
    }
    if (*par_n == 0) free(par_pos);
    return par_pos;
}

// TODO do extension first, so that cons would start from the begining of the repeat
int *get_partition_pos_with_narrow_global_alignment(uint8_t *bseq, int seq_len, dp_t **dp, chain_t ch, mini_tandem_para *mtp, int *par_n) {
    int est_ch_i = ch.est_ch_i, est_start = ch.est_start, est_period = ch.est_period;
    int first_end = dp[ch.cell[0].i][ch.cell[0].j].end, last_start = dp[ch.cell[ch.len-1].i][ch.cell[ch.len-1].j].start;
    int *par_pos = (int*)_err_malloc(seq_len * sizeof(int)); *par_n = 0;
    int i, j, k=mtp->k, ch_i, s, e, s1, e1, s2, e2, iden_n;
    int n_cigar; uint32_t *cigar;
    cell_t c; dp_t d;
    // find upstream anchor (s1, e1) and downstream (s2, e2)
    // global alignment of [s1,s2] and [e1, e2]
    // find best partition position by backtracking

    // left extension: find S
    // <=<=<=<=<=<=<=<=<=<=<=<=<=<=<=<=<=<=<=
    // --s1---S---s2----e1---s---e2-------e--
    // --|---------| vs |---------|----------
    ch_i = est_ch_i, s = est_start, e = est_start + est_period;
    while (s >= first_end && ch_i > 0) {
        // find (s1, e1) and (s2, e2)
        s2 = s, e2 = e; s1 = -1; e1 = -1;
        for (i = ch_i - 1; i >= 0; --i) {
            c = ch.cell[i]; d = dp[c.i][c.j];
            s1 = d.start; e1 = d.end;
            if (e1 == s) {
                par_pos[(*par_n)++] = s1;
                ch_i = i, s = s1, e = e1;
                break;
            } else if (e1 < s) {
                // if (s - e1 < k) // overlapped
                // do global alignment
                iden_n = ksw2_global_with_cigar(bseq+e1-k+1, e2-e1+k, bseq+s1-k+1, s2-s1+k, &n_cigar, &cigar);
                #ifdef __DEBUG__
                printf("1: (%d,%d), 2: (%d,%d)\n", s1, e1, s2, e2);
                printf("iden_n: %d (%d,%d), (%d,%d)\n", iden_n, e2-e1+k, s2-s1+k, s1, s2);
                #endif
                if (iden_n >= MIN_OF_TWO(s2-s1+k, e2-e1+k) * (1-mtp->max_div)) { // extend partition
                    e = s; s = s2 - backtrack_left_end(n_cigar, cigar, e2-e1+k, s2-s1+k, e2-s);
                    if (e == s) { // no backtrack
                        ch_i = 0; break;
                    }
                    par_pos[(*par_n)++] = s; 
                    ch_i = i+1;
                } else { // skip this anchor
                    par_pos[(*par_n)++] = -1; // insert a separation flag
                    par_pos[(*par_n)++] = e1; 
                    par_pos[(*par_n)++] = s1; 
                    ch_i = i; s = s1; e = e1;
                }
                if (cigar) free(cigar);
                break;
            } else {
                s2 = s1; e2 = e1;
            }
        }
    }
    // reverse par_pos
    for (i = 0; i < (*par_n)>> 1; ++i) {
        j = par_pos[i]; par_pos[i] = par_pos[*par_n-i-1]; par_pos[*par_n-i-1] = j;
    }
    par_pos[(*par_n)++] = est_start;
    par_pos[(*par_n)++] = est_start + est_period;

    // right extension: find E
    // >=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=
    // --s-------s1---e---s2----e1---E---e2--
    // ----------|---------| vs |---------|--
    ch_i = est_ch_i, s = est_start, e = est_start + est_period;
    while (ch_i < ch.len-1 && e <= last_start) {
        // find (s1,e1) and (s2,e2)
        s1 = s, e1 = e, s2 = e2 = -1;
        for (i = ch_i+1; i < ch.len; ++i) {
            c = ch.cell[i]; d = dp[c.i][c.j];
            s2 = d.start; e2 = d.end;
            if (s2 == e) {
                par_pos[(*par_n)++] = e2;
                ch_i = i; s = s2; e = e2;
                break;
            } else if (s2 > e) { // first (s2, e2)
                // if (s2 - e < k) // overlapped
                iden_n = ksw2_global_with_cigar(bseq+s1-k+1, s2-s1+k, bseq+e1-k+1, e2-e1+k, &n_cigar, &cigar);
                #ifdef __DEBUG__
                printf("1: (%d,%d), 2: (%d,%d)\n", s1, e1, s2, e2);
                printf("iden_n: %d (%d,%d), (%d,%d)\n", iden_n, e2-e1+k, s2-s1+k, s1, s2);
                #endif
                if (iden_n >= MIN_OF_TWO(s2-s1+k, e2-e1+k) * (1-mtp->max_div)) {
                    s = e; e = e2 - backtrack_left_end(n_cigar, cigar, s2-s1+k, e2-e1+k, s2-e);
                    if (e == s) { // no backtrack
                        ch_i = ch.len; break;
                    }
                    par_pos[(*par_n)++] = e; 
                    ch_i = i-1;
                } else {
                    par_pos[(*par_n)++] = -1; 
                    par_pos[(*par_n)++] = s2; 
                    par_pos[(*par_n)++] = e2; 
                    ch_i = i; s = s2; e = e2;
                }
                if (cigar) free(cigar);
                break;
            } else {
                s1 = s2; e1 = e2;
            }
        }
    }
    if (*par_n == 0) free(par_pos);
    return par_pos;
}
