#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "tidehunter.h"
#include "utils.h"
#include "seq.h"
#include "abpoa.h"

#define NAT_E 2.718281828459045

abpoa_para_t *mt_abpoa_init_para(mini_tandem_para *mtp) {
    abpoa_para_t *abpt = abpoa_init_para();
    abpt->cons_agrm = 1; // 0: HB, 1: RC
    abpt->match = mtp->match;     // match score
    abpt->mismatch = mtp->mismatch;  // mismatch penalty
    abpt->gap_open1 = mtp->gap_open1; // first gap open penalty
    abpt->gap_ext1 = mtp->gap_ext1;  // first gap extension penalty 
    abpt->gap_open2 = mtp->gap_open2; // second gap open penalty
    abpt->gap_ext2 = mtp->gap_ext2;  // second gap extension penalty 

    abpt->out_cons = 1;
#ifdef __DEBUG__
    abpt->out_msa = 1; // generate Row-Column multiple sequence alignment(RC-MSA), set 0 to disable
#endif
    abpoa_post_set_para(abpt);
    return abpt;
}

int abpoa_gen_cons(abpoa_t *ab, abpoa_para_t *abpt, uint8_t *bseqs, int seq_len, int *pos, int pos_n, uint8_t *cons_bseq, uint8_t *cons_qual, mini_tandem_para *mtp, int *_n_seqs) {
    int i, n_seqs, cons_len = 0;

    /* clean graph if it is re-used */
    abpoa_reset_graph(ab, abpt, seq_len);
    ab->abs->n_seq = 0;

    int *seq_lens = (int*)malloc(sizeof(int) * (pos_n-1));
    uint8_t **_bseqs = (uint8_t**)malloc(sizeof(uint8_t*) * (pos_n-1));
    /* main graph alignment */
    // |pos|-----|pos|-----pos|
    for (i = n_seqs = 0; i < pos_n-1; ++i) {
        int start = pos[i], end = pos[i+1];
        if (start < 0 || end < 0 || start >= seq_len-1 || end+1 > seq_len) continue;
        // fprintf(stdout, ">%d\n", start);
        seq_lens[n_seqs] = end - start;
        _bseqs[n_seqs] = bseqs + start + 1;
        /*int j; for (j = start; j < end; ++j) fprintf(stdout, "%c", "ACGT"[bseqs[j+1]]); fprintf(stdout, "\n");*/
        ++n_seqs;
    }
#ifdef __DEBUG__
    FILE *outfp = stderr;
#else
    FILE *outfp = NULL;
#endif

    *_n_seqs = n_seqs;
    int min_cov = 0, _min_cov;
    if (mtp->min_frac > 0.0) min_cov = (int)(n_seqs * mtp->min_frac);
    else if (mtp->min_cov > 0) min_cov = mtp->min_cov;
    if (n_seqs <= 2) {
        if (n_seqs <= 1) err_fatal_simple("No enough sequences to perform msa.\n");
        cons_len = seq_lens[0];

        int skip = 0;
        if (min_cov > 0) {
            if (seq_lens[0] != seq_lens[1]) _min_cov = 1;
            else {
                _min_cov = 2;
                for (i = 0; i < cons_len; ++i) {
                    if (_bseqs[0][i] != _bseqs[1][i]) {
                        _min_cov = 1;
                        break;
                    }
                }
            }
            if (_min_cov < min_cov) skip = 1;
        }
        if (skip == 0) {
            for (i = 0; i < cons_len; ++i) {
                cons_bseq[i] = _bseqs[0][i];
                if (cons_qual != NULL) cons_qual[i] = 33;
            }
        } else cons_len = 0;
    } else {
        uint8_t **_cons_bseq; int **_cons_cov, *_cons_l, _cons_n = 0;
        if (min_cov > 0 || cons_qual != NULL) abpoa_msa(ab, abpt, n_seqs, NULL, seq_lens, _bseqs, outfp, &_cons_bseq, &_cons_cov, &_cons_l, &_cons_n, NULL, NULL);
        else abpoa_msa(ab, abpt, n_seqs, NULL, seq_lens, _bseqs, outfp, &_cons_bseq, NULL, &_cons_l, &_cons_n, NULL, NULL);
        if (_cons_n == 1) {
            cons_len = _cons_l[0];
            int skip = 0;
            if (min_cov > 0 || cons_qual != NULL) {
                if (min_cov > 0) {
                    for (i = 0; i < cons_len; ++i) {
                        if (_cons_cov[0][i] < min_cov) {
                            skip = 1; break;
                        }
                    }
                }
                if (cons_qual != NULL) {
                    int phred; double x, p;
                    for (i = 0; i < cons_len; ++i) {
                        // min: 0+33=33, max: 60+33=93
                        x = 13.8 * (1.25 * _cons_cov[0][i] / n_seqs - 0.25);
                        p = 1 - 1.0 / (1.0 + pow(NAT_E, -1 * x));
                        phred = 33 + (int)(-10 * log10(p) + 0.499);
                        cons_qual[i] = phred;
                    }
                }
                free(_cons_cov[0]); free(_cons_cov);
            }
            if (skip == 0) {
                for (i = 0; i < cons_len; ++i) cons_bseq[i] = _cons_bseq[0][i];
            } else cons_len = 0;
            free(_cons_l); free(_cons_bseq[0]); free(_cons_bseq);
        }
    }

    free(seq_lens); free(_bseqs);
    return cons_len;
}
