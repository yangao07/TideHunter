#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "tidehunter.h"
#include "utils.h"
#include "seq.h"

abpoa_para_t *mt_abpoa_init_para(mini_tandem_para *mtp) {
    abpoa_para_t *abpt = abpoa_init_para();
    abpt->m = 5;
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

int abpoa_gen_cons(abpoa_t *ab, abpoa_para_t *abpt, uint8_t *bseqs, int seq_len, int *pos, int pos_n, uint8_t *cons_bseq) {
    int i, seq_n, cons_len = 0;

    /* clean graph if it is re-used */
    abpoa_reset_graph(ab, abpt, seq_len);

    int *seq_lens = (int*)malloc(sizeof(int) * (pos_n-1));
    uint8_t **_bseqs = (uint8_t**)malloc(sizeof(uint8_t*) * (pos_n-1));
    /* main graph alignment */
    // |pos|-----|pos|-----pos|
    for (i = seq_n = 0; i < pos_n-1; ++i) {
        int start = pos[i], end = pos[i+1];
        if (start < 0 || end < 0 || start >= seq_len || end+1 >= seq_len) continue;
        // fprintf(stdout, ">%d\n", start);
        seq_lens[seq_n] = end - start;
        _bseqs[seq_n] = bseqs + start + 1;
        /*int j;
        for (j = start; j < end; ++j)
            fprintf(stdout, "%c", "ACGT"[bseqs[j+1]]);
        fprintf(stdout, "\n");*/
        ++seq_n;
    }
#ifdef __DEBUG__
    FILE *outfp = stderr;
#else
    FILE *outfp = NULL;
#endif
    if (seq_n <= 2) {
        if (seq_n == 0) err_fatal_simple("No enough sequences to perform msa.\n");
        cons_len = seq_lens[0];
        for (i = 0; i < cons_len; ++i) cons_bseq[i] = _bseqs[0][i];
    } else {
        uint8_t **_cons_bseq; int *_cons_l, _cons_n = 0;
        abpoa_msa(ab, abpt, seq_n, NULL, seq_lens, _bseqs, outfp, &_cons_bseq, NULL, &_cons_l, &_cons_n, NULL, NULL);
        if (_cons_n == 1) {
            for (i = 0; i < _cons_l[0]; ++i) cons_bseq[i] = _cons_bseq[0][i];
            cons_len = _cons_l[0];

            free(_cons_l); free(_cons_bseq[0]); free(_cons_bseq);
        }
    }
    // abpoa_msa(ab, abpt, seq_n, NULL, seq_lens, _bseqs, stderr, NULL, NULL, NULL, NULL, NULL, NULL);

    free(seq_lens); free(_bseqs);
    // abpoa_free(ab, abpt); abpoa_free_para(abpt); 
    return cons_len;
}
