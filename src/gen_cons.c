#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "tide_hunter.h"
#include "edlib_align.h"
#include "abpoa_align.h"
#include "ksw2_align.h"
#include "spoa_align.h"
#include "utils.h"

void write_tandem_cons_seq(tandem_seq_t *tseq, char *cons_seq, int cons_len, int start, int end, double copy_num, mini_tandem_para *mtp, int8_t splint_rotated) {
    if (cons_len < mtp->min_p || cons_len > mtp->max_p) return;
    if (mtp->only_longest && tseq->cons_n == 1) {
        if (end-start > tseq->cons_end[0]-tseq->cons_start[0]) {
            tseq->cons_n = 0; tseq->cons_seq->seq.l = 0;
        } else return;
    }
    if (tseq->cons_seq->seq.l + cons_len >= tseq->cons_seq->seq.m) {
        tseq->cons_seq->seq.m = tseq->cons_seq->seq.l + cons_len + 1;
        tseq->cons_seq->seq.s = (char*)_err_realloc(tseq->cons_seq->seq.s, tseq->cons_seq->seq.m * sizeof(char));
    }
    strcpy(tseq->cons_seq->seq.s + tseq->cons_seq->seq.l, cons_seq); tseq->cons_seq->seq.l += cons_len; 

    if (tseq->cons_n == tseq->cons_m) {
        tseq->cons_m <<= 1;
        tseq->cons_start = (int*)_err_realloc(tseq->cons_start, tseq->cons_m * sizeof(int));
        tseq->cons_end = (int*)_err_realloc(tseq->cons_end, tseq->cons_m * sizeof(int));
        tseq->copy_num = (double*)_err_realloc(tseq->copy_num, tseq->cons_m * sizeof(double));
        tseq->splint_rotated = (int8_t*)_err_realloc(tseq->splint_rotated, tseq->cons_m * sizeof(int8_t));
        tseq->cons_len = (int*)_err_realloc(tseq->cons_len, tseq->cons_m * sizeof(int));
        tseq->cons_score = (int*)_err_realloc(tseq->cons_score, tseq->cons_m * sizeof(int));
    }
    tseq->cons_start[tseq->cons_n] = start; tseq->cons_end[tseq->cons_n] = end; tseq->copy_num[tseq->cons_n] = copy_num;
    tseq->splint_rotated[tseq->cons_n] = splint_rotated;
    tseq->cons_len[tseq->cons_n] = cons_len; tseq->cons_score[tseq->cons_n] = 0; // TODO cons_score
    ++tseq->cons_n;
}

void seqs_msa(int seq_len, uint8_t *bseq, int par_n, int *par_pos, tandem_seq_t *tseq, mini_tandem_para *mtp) {
        int start, end, cons_len=0;
        char *cons_seq = (char*)_err_malloc(seq_len * sizeof(char));
        uint8_t *cons_bseq = (uint8_t*)_err_malloc(seq_len * sizeof(uint8_t));
#ifdef __DEBUG__
        {
            int i, j, seq_i = 0; 
            for (i = 0; i < par_n-1; ++i) {
                if (par_pos[i] >= 0 && par_pos[i+1] >= 0) {
                    start = par_pos[i], end = par_pos[i+1];
                    printf(">seqs_%d:%d-%d\n", end-start, start, end);
                    for (j = start+1; j <= end; ++j) printf("%c", seq[j]);
                    printf("\n");
                }
            }
        }
#endif
        int i = 0, j, s; 
        while (i < par_n-mtp->min_copy) {
            if (par_pos[i] < 0) { 
                i++;
                continue;
            }
            for (j = i+1; j < par_n; ++j) {
                if (par_pos[j] < 0) break;
            }
            if (j - i > mtp->min_copy) { // do multiple sequence alignment and consensus calling for par_pos[i:j]
                cons_len = abpoa_msa(bseq, seq_len, par_pos+i, j-i, cons_bseq);
                for (s = 0; s < cons_len; ++s) cons_seq[s] = "ACGTN"[cons_bseq[s]]; cons_seq[cons_len] = '\0';

                int max_q, max_t, cons_start, cons_end; double copy_num = j-i-1;
                ksw2_left_ext(cons_bseq, cons_len, bseq, par_pos[i]+1, &max_q, &max_t); cons_start = par_pos[i] - max_t;
                // printf("max_q: %d, max_t: %d\n", max_q, max_t);
                copy_num += (max_q + 1.0) / cons_len;
                ksw2_right_ext(cons_bseq, cons_len, bseq+par_pos[j-1]+1, seq_len-par_pos[j-1]-1, &max_q, &max_t); cons_end = par_pos[j-1] + max_t + 1;
                // printf("max_q: %d, max_t: %d\n", max_q, max_t);
                copy_num += (max_q + 1.0) / cons_len;
                // rotate the cons based on splint_seq 
                int splint_rotated = 0;
                if (mtp->splint_seq != NULL && cons_len > mtp->splint_len) { // TODO use ksw
                    int ed, min_ed = mtp->splint_len, sp_len = mtp->splint_len, min_start, min_end, idx = -1;
                    // search splint within a full cons, forward and reverse
                    ed = edlib_align_HW(mtp->splint_seq, mtp->splint_len, cons_seq, cons_len, &start, &end);
                    if (ed >= 0 && sp_len * 0.9 <= end-start && sp_len * 1.1 >= end-start && ed < min_ed) {
                        min_ed = ed; idx = 0; // Forward single cons
                        min_start = start; min_end = end;
                    }
                    ed = edlib_align_HW(mtp->splint_rc_seq, mtp->splint_len, cons_seq, cons_len, &start, &end);
                    if (ed >= 0 && sp_len * 0.9 <= end-start && sp_len * 1.1 >= end-start && ed < min_ed) {
                        min_ed = ed; idx = 1; // Reverse-comp single cons
                        min_start = start; min_end = end;
                    }
                    // search splint within a concatenated 2 cons, forward and reverse
                    char *cons2 = (char*)_err_malloc(((cons_len << 1) + 1) * sizeof(char));
                    strcpy(cons2, cons_seq); strcpy(cons2+cons_len, cons_seq); cons2[cons_len<<1] = '\0'; // concatenated 2 copies
                    ed = edlib_align_HW(mtp->splint_seq, mtp->splint_len, cons2, cons_len<<1, &start, &end);
                    if (ed >= 0 && ed < min_ed && sp_len * 0.9 <= end-start && sp_len * 1.1 >= end-start && start < cons_len && end >= cons_len) {
                        min_ed = ed; idx = 2; // Forward 2 copies cons
                        min_start = start; min_end = end;
                    }
                    ed = edlib_align_HW(mtp->splint_rc_seq, mtp->splint_len, cons2, cons_len<<1, &start, &end);
                    if (ed >= 0 && ed < min_ed && sp_len * 0.9 <= end-start && sp_len * 1.1 >= end-start && start < cons_len && end >= cons_len) {
                        min_ed = ed; idx = 3; // Reverse-comp 2 copies cons
                        min_start = start; min_end = end;
                    }
                    // rotate cons based on ed result
                    switch(idx)
                    {
                        case 0: 
                        case 1:
                            if (cons_len-min_end-1 > 0 && cons_len-min_end-1+min_start >= mtp->min_p && cons_len-min_end-1+min_start <= mtp->max_p) {
                                memcpy(cons_seq, cons2 + min_end + 1, cons_len - min_end - 1);
                                memcpy(cons_seq+cons_len-min_end-1, cons2, min_start);
                                cons_seq[cons_len - min_end - 1 + min_start] = '\0';
                                cons_len = cons_len - min_end - 1 + min_start;
                                splint_rotated = 1;
                            } else cons_len = 0;
                            break;
                        case 2:
                        case 3:
                            min_end -= cons_len;
                            if (min_start-min_end-1 > mtp->min_p && min_start-min_end-1 <= mtp->max_p) {
                                memcpy(cons_seq, cons2 + min_end + 1, min_start - min_end - 1);
                                cons_seq[min_start-min_end-1] = '\0';
                                cons_len = min_start-min_end-1;
                                splint_rotated = 1;
                            } else cons_len = 0;
                            break;
                        default:;
                            // err_printf("%s: no splint sequence found in consensus sequence.\n", read_name);
                    }
                    free(cons2);
                } 
                // TODO find full-length sequence based on 5' and 3' adapter sequences
                // int full_length = 0;
                // if (mtp->five_seq != NULL && mtp->three_seq != NULL && cons_len > mtp->five_len + mtp->three_len) {
                // }
                // write cons_seq into tseq
                write_tandem_cons_seq(tseq, cons_seq, cons_len, cons_start, cons_end, copy_num, mtp, splint_rotated);
            }
            i = j + 1;
        }
        free(cons_seq); free(cons_bseq); free(par_pos);
}
