#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "tide_hunter.h"
#include "edlib_align.h"
#include "abpoa_cons.h"
#include "ksw2_align.h"
#include "utils.h"

void write_tandem_cons_seq(tandem_seq_t *tseq, char *cons_seq, int cons_len, int start, int end, double copy_num, mini_tandem_para *mtp, int8_t full_length, int *par_pos, int pos_n) {
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

    int i;
    if (tseq->cons_n == tseq->cons_m) {
        tseq->cons_m <<= 1;
        tseq->cons_start = (int*)_err_realloc(tseq->cons_start, tseq->cons_m * sizeof(int));
        tseq->cons_end = (int*)_err_realloc(tseq->cons_end, tseq->cons_m * sizeof(int));
        tseq->copy_num = (double*)_err_realloc(tseq->copy_num, tseq->cons_m * sizeof(double));
        tseq->full_length = (int8_t*)_err_realloc(tseq->full_length, tseq->cons_m * sizeof(int8_t));
        tseq->cons_len = (int*)_err_realloc(tseq->cons_len, tseq->cons_m * sizeof(int));
        tseq->cons_score = (int*)_err_realloc(tseq->cons_score, tseq->cons_m * sizeof(int));
        tseq->sub_pos = (int**)_err_realloc(tseq->sub_pos, tseq->cons_m * sizeof(int*));
        tseq->pos_n = (int*)_err_realloc(tseq->pos_n, tseq->cons_m * sizeof(int));
        tseq->pos_m = (int*)_err_realloc(tseq->pos_m, tseq->cons_m * sizeof(int));
        for (i = tseq->cons_n; i < tseq->cons_m; ++i) {
            tseq->pos_m[i] = 0; tseq->pos_n[i] = 0;
            tseq->sub_pos[i] = 0;
        }
    }
    tseq->cons_start[tseq->cons_n] = start; tseq->cons_end[tseq->cons_n] = end; tseq->copy_num[tseq->cons_n] = copy_num;
    tseq->full_length[tseq->cons_n] = full_length;
    tseq->cons_len[tseq->cons_n] = cons_len; tseq->cons_score[tseq->cons_n] = 0; // TODO cons_score
    tseq->pos_n[tseq->cons_n] = pos_n;
    if (pos_n > tseq->pos_m[tseq->cons_n]) {
        tseq->pos_m[tseq->cons_n] = pos_n;
        tseq->sub_pos[tseq->cons_n] = (int*)_err_realloc(tseq->sub_pos[tseq->cons_n], pos_n * sizeof(int));
    }
    for (i = 0; i < pos_n; ++i) tseq->sub_pos[tseq->cons_n][i] = par_pos[i];
    ++tseq->cons_n;
}

void seqs_msa(int seq_len, uint8_t *bseq, int par_n, int *par_pos, tandem_seq_t *tseq, mini_tandem_para *mtp) {
        int cons_len=0;
        char *cons_seq = (char*)_err_malloc(seq_len * sizeof(char));
        uint8_t *cons_bseq = (uint8_t*)_err_malloc(seq_len * sizeof(uint8_t));
#ifdef __DEBUG__
        {
            int i, j, seq_i = 0, start, end; 
            for (i = 0; i < par_n-1; ++i) {
                if (par_pos[i] >= 0 && par_pos[i+1] >= 0) {
                    start = par_pos[i], end = par_pos[i+1];
                    printf(">seqs_%d:%d-%d\n", end-start, start, end);
                    for (j = start+1; j <= end; ++j) printf("%c", "ACGTN"[bseq[j]]);
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
                // find full-length sequence based on 5' and 3' adapter sequences
                int full_length = 0;
                if (mtp->five_seq != NULL && mtp->three_seq != NULL && cons_len > mtp->five_len + mtp->three_len) {
                    char *cons2 = (char*)_err_malloc(((cons_len << 1) + 1) * sizeof(char));
                    strcpy(cons2, cons_seq); strcpy(cons2+cons_len, cons_seq); cons2[cons_len<<1] = '\0'; // concatenated 2 copies
                    int tar_start, tar_end, tot_ed, _5_ed, _3_ed, _5_start, _5_end, _3_start, _3_end;
                    tar_start = tar_end = -1;
                    // |---plus-5'adapter---|---target-sequence---|---minus-3'adapter---|
                    _5_start = _5_end = _3_start = _3_end = -1; tot_ed = INT32_MAX;
                    _5_ed = edlib_align_HW(mtp->five_seq, mtp->five_len, cons2, cons_len<<1, &_5_start, &_5_end);
                    if (_5_ed > mtp->five_len * (1-mtp->ada_match_rat)) goto REV;
                    _3_ed = edlib_align_HW(mtp->three_rc_seq, mtp->three_len, cons2, cons_len<<1, &_3_start, &_3_end);
                    if (_3_ed > mtp->three_len * (1-mtp->ada_match_rat)) goto REV; 
                    if (_3_start <= _5_end) {
                        if (_3_end + cons_len < cons_len << 1 && _3_start + cons_len > _5_end) {
                            tar_start = _5_end + 1;
                            tar_end = _3_start + cons_len - 1;
                            full_length = 1;
                            tot_ed = _5_ed + _3_ed;
                        }
                    } else {
                        tar_start = _5_end + 1;
                        tar_end = _3_start - 1;
                        tot_ed = _5_ed + _3_ed;
                        full_length = 1;
                    }
#ifdef __DEBUG__
                    printf("FOR: %d, %d => %d, %d\n", tar_start, tar_end, _5_ed, _3_ed);
#endif
                    if (tot_ed == 0) goto WRITE_CONS;
                    // |---plus-3'adapter---|---target-sequence---|---minus-5'adapter---|
REV:
                    _5_ed = edlib_align_HW(mtp->five_rc_seq, mtp->five_len, cons2, cons_len<<1, &_5_start, &_5_end);
                    if (_5_ed > mtp->five_len * (1-mtp->ada_match_rat)) goto WRITE_CONS;
                    _3_ed = edlib_align_HW(mtp->three_seq, mtp->three_len, cons2, cons_len<<1, &_3_start, &_3_end);
                    if (_3_ed > mtp->three_len * (1-mtp->ada_match_rat)) goto WRITE_CONS; 
                    if (_5_ed + _3_ed < tot_ed) {
                        if (_5_start <= _3_end) {
                            if (_5_end + cons_len < cons_len << 1 && _5_start + cons_len > _3_end) {
                                tar_start = _3_end + 1;
                                tar_end = _5_start + cons_len - 1;
                                full_length = 2;
                            }
                        } else {
                            tar_start = _3_end + 1;
                            tar_end = _5_start - 1;
                            full_length = 2;
                        }
#ifdef __DEBUG__
                        printf("REV: %d, %d => %d\n", tar_start, tar_end, _5_ed + _3_ed);
#endif
                    }
WRITE_CONS:
                    if (tar_start > 0 && tar_end > tar_start) {
                        memcpy(cons_seq, cons2 + tar_start , tar_end - tar_start + 1);
                        cons_seq[tar_end - tar_start + 1] = '\0';
                        cons_len = tar_end - tar_start + 1;
                    }
                }
                // write cons_seq into tseq
                if (!(mtp->only_full_length) || full_length > 0)
                    write_tandem_cons_seq(tseq, cons_seq, cons_len, cons_start, cons_end, copy_num, mtp, full_length, par_pos+i, j-i);
            }
            i = j + 1;
        }
        free(cons_seq); free(cons_bseq); free(par_pos);
}
