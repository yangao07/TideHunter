#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "poa.h"
#include "poa_graph.h"
#include "utils.h"
#include "poa_align.h"
#include "align.h"

poa_para_t *poa_init_para(void) {
    poa_para_t *ppt = (poa_para_t*)_err_malloc(sizeof(poa_para_t));
    ppt->align_mode = POA_GLOBAL_FLAG;

    // number of residue types
    ppt->m = 5; // nucleotides
    ppt->mat = NULL;

    // score matrix
    ppt->match = POA_MATCH;
    ppt->mismatch = POA_MISMATCH;
    ppt->gap_open = POA_GAP_OPEN;
    ppt->gap_ext = POA_GAP_EXT;

    ppt->bw = 100; // TODO 
    ppt->zdrop = 100;
    ppt->end_bonus = 5;

    return ppt;
}

void poa_free_para(poa_para_t *ppt) {
    if (ppt->mat != NULL) free(ppt->mat);
    free(ppt);
}


// int poa_main(const char *seq_fn, poa_para_t *ppt) { TODO
int poa_main(int seq_n, char (*seq)[100], poa_para_t *ppt){
    int i, j;
    poa_graph_t *graph = poa_init_graph(1);

    uint8_t *bseq = (uint8_t*)_err_malloc(100 * sizeof(uint8_t)); int bseq_m = 100;
    // progressively partial order alignment
    for (i = 0; i < seq_n; ++i) {
        int seq_l = strlen(seq[i]);
        if (seq_l > bseq_m) {
            bseq_m = seq_l;
            bseq = (uint8_t*)_err_realloc(bseq, bseq_m * sizeof(uint8_t));
        }
        for (j = 0; j < seq_l; ++j) bseq[j] = nst_nt4_table[(int)(seq[i][j])];
        poa_cigar_t *poa_cigar=0; int n_cigar=0;
        poa_align_sequence_with_graph(graph, bseq, seq_l, ppt, &n_cigar, &poa_cigar);
        poa_add_graph_alignment(graph, bseq, seq_l, n_cigar, poa_cigar);
        if (n_cigar) free(poa_cigar);
    }
    for (i = 0; i < seq_n; ++i) {
        printf("%s\n", seq[i]);
    }

    // generate consensus from graph
    int cons_seq_l = 0;
    uint8_t *cons_seq = poa_generate_consensus(graph, &cons_seq_l);

    poa_free_graph(graph, 1); free(bseq); if (cons_seq_l > 0) free(cons_seq);
    return 0;
}

int main(int argc, char **argv) {
    int seq_n = 6;
    char seq[100][100] = {
        //"AACAT",
        //"AACT"
        //"ACGTAA",
        //"TTACGGAAGG"
        //"ACGGAA",
        //"ACGTAA"
        //"ACGGTTAA",
        //"ACTTGGAA"
        //"AGGGTTAA",
        //"AGTACCAA",
        //"AGGGCCAA"
        "CATAAAAGAACGTAGGTCGCCCGTCCGTAACCTGTCGGATCACCGGAAAGGACCCGTAAAGTGATAATGAT",
        "ATAAAGGCAGTCGCTCTGTAAGCTGTCGATTCACCGGAAAGATGGCGTTACCACGTAAAGTGATAATGATTAT",
        "ATCAAAGAACGTGTAGCCTGTCCGTAATCTAGCGCATTTCACACGAGACCCGCGTAATGGG",
        "CGTAAATAGGTAATGATTATCATTACATATCACAACTAGGGCCGTATTAATCATGATATCATCA",
        "GTCGCTAGAGGCATCGTGAGTCGCTTCCGTACCGCAAGGATGACGAGTCACTTAAAGTGATAAT",
        "CCGTAACCTTCATCGGATCACCGGAAAGGACCCGTAAATAGACCTGATTATCATCTACAT"
    };
    poa_para_t *ppt = poa_init_para();
    // parse argv TODO

    ppt->mat = (int8_t*)_err_malloc(ppt->m * ppt->m * sizeof(int8_t));
    gen_simple_mat(ppt->m, ppt->mat, ppt->match, ppt->mismatch);
    poa_main(seq_n, seq, ppt);
    poa_free_para(ppt);

    return 0;
}
