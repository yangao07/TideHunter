#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "abpoa.h"

int abpoa_msa(uint8_t *bseqs, int seq_len, int *pos, int pos_n, uint8_t *cons_bseq) {
    int i, seq_n, cons_len = 0;

    abpoa_t *ab = abpoa_init();
    abpoa_para_t *abpt = abpoa_init_para();

    // score parameters
    abpt->m = 5;
    abpt->match = 5;     // match score
    abpt->mismatch = 4;  // mismatch penalty
    abpt->gap_mode = 0;  // linear gap penalty
    abpt->gap_open1 = 0; // gap open penalty
    abpt->gap_ext1 = 8;  // gap extension penalty 
                         // gap_penalty = gap_open + gap_len * gap_ext
    abpt->mat = (int*)malloc(abpt->m * abpt->m * sizeof(int));
    gen_simple_mat(abpt->m, abpt->mat, abpt->match, abpt->mismatch);
    // output options
    #ifdef __DEBUG__
    abpt->out_msa = 1; // generate Row-Column multiple sequence alignment(RC-MSA), set 0 to disable
    //printf("SIMD: %d\n", abpt->simd_flag);
    abpt->use_read_ids = 1;
    set_65536_table();
    #else
    abpt->out_msa = 0;
    #endif
    abpt->out_cons = 1; // generate consensus sequence, set 0 to disable
    // abpt->out_pog = 1; // generate parital order graph using DOT, set 0 to disable

    /* clean graph if it is re-used */
    // abpoa_reset_graph(ab, seq_len, abpt);

    /* main graph alignment */
    // |pos|-----|pos|-----pos|
    for (i = seq_n = 0; i < pos_n-1; ++i) {
        int start = pos[i], end = pos[i+1];
        if (start < 0 || end < 0 || start >= seq_len || end+1 >= seq_len) continue;

        abpoa_cigar_t *abpoa_cigar = 0; int n_cigar = 0;
        abpoa_align_sequence_with_graph(ab, bseqs+start+1, end-start, abpt, &n_cigar, &abpoa_cigar);

        abpoa_add_graph_alignment(ab->abg, abpt, bseqs+start+1, end-start, n_cigar, abpoa_cigar, i, (pos_n-2)/64+1);
        if (n_cigar) free(abpoa_cigar);
        ++seq_n;
    }
    /* generate consensus sequence from graph */
    if (abpt->out_cons && ab->abg->node_n > 2) {
        cons_len = abpoa_generate_consensus(ab->abg, abpt->cons_agrm, 1, 0.0, seq_n, NULL, NULL, NULL, NULL);
        for (i = 0; i < cons_len; ++i) cons_bseq[i] = ab->abg->cons_seq[i];
    }

    /* generate multiple sequence alignment */
    if (abpt->out_msa &&  ab->abg->node_n > 2)
        abpoa_generate_multiple_sequence_alingment(ab->abg, seq_n, stderr);

    abpoa_free(ab, abpt); abpoa_free_para(abpt); 
    return cons_len;
}
