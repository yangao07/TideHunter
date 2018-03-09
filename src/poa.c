#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "poa_graph.h"
#include "poa.h"
#include "utils.h"
#include "poa_align.h"
#include "align.h"
#include "poa_graph_visual.h"

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

    ppt->bw = 300; // TODO band width
    ppt->zdrop = 100;
    ppt->end_bonus = 5;

    return ppt;
}

void poa_free_para(poa_para_t *ppt) {
    if (ppt->mat != NULL) free(ppt->mat);
    free(ppt);
}


void cons_to_seq_score(int cons_l, uint8_t *cons_seq, int seq_n, char (*seq)[100], poa_para_t *ppt) {
    int i, j;
    printf("cons_to_seq:\n");
    uint8_t *bseq = (uint8_t*)_err_malloc(100 * sizeof(uint8_t)); int bseq_m = 100;
    // progressively partial order alignment
    for (i = 0; i < seq_n; ++i) {
        poa_graph_t *graph = poa_init_graph(1);
        int seq_l = strlen(seq[i]);
        if (seq_l > bseq_m) {
            bseq_m = seq_l;
            bseq = (uint8_t*)_err_realloc(bseq, bseq_m * sizeof(uint8_t));
        }
        poa_cigar_t *poa_cigar=0; int n_cigar=0;
        poa_align_sequence_with_graph(graph, cons_seq, cons_l, ppt, &n_cigar, &poa_cigar); 
        poa_add_graph_alignment(graph, cons_seq, cons_l, n_cigar, poa_cigar, NULL, NULL); if (n_cigar) free(poa_cigar);

        for (j = 0; j < seq_l; ++j) bseq[j] = nst_nt4_table[(int)(seq[i][j])];
        poa_align_sequence_with_graph(graph, bseq, seq_l, ppt, &n_cigar, &poa_cigar); if (n_cigar) free(poa_cigar);
        poa_free_graph(graph, 1);
    }

    int rest_n = 2, rest_l, rest_i;
    char rest[100][100] = { 
        "GTCGTAAAGAACGTAGGTCGCCCGTCCGTAATCTGTCGGATCACCGGAAAGATGACGACCCGTAAAGTGATAATGATCAT",
        "CATAAAGAACGTAGGTCGCCGTGAGTCCGTAATCCGTACGGATTCACCGGAATGGCGTAGTTACCCGATAAAGTGATAATACAT",
    };
    uint8_t rest_seq[100];
    for (rest_i = 0; rest_i < rest_n; ++rest_i) {
        printf("REST# %d: cons_to_seq:\n", rest_i);
        rest_l = strlen(rest[rest_i]);
        for (i = 0; i < rest_l; ++i) rest_seq[i] = nst_nt4_table[(int)(rest[rest_i][i])];
        // progressively partial order alignment
        for (i = 0; i < seq_n; ++i) {
            poa_graph_t *graph = poa_init_graph(1);
            int seq_l = strlen(seq[i]);
            if (seq_l > bseq_m) {
                bseq_m = seq_l;
                bseq = (uint8_t*)_err_realloc(bseq, bseq_m * sizeof(uint8_t));
            }
            poa_cigar_t *poa_cigar=0; int n_cigar=0;
            poa_align_sequence_with_graph(graph, rest_seq, rest_l, ppt, &n_cigar, &poa_cigar); 
            poa_add_graph_alignment(graph, rest_seq, rest_l, n_cigar, poa_cigar, NULL, NULL); if (n_cigar) free(poa_cigar);

            for (j = 0; j < seq_l; ++j) bseq[j] = nst_nt4_table[(int)(seq[i][j])];
            poa_align_sequence_with_graph(graph, bseq, seq_l, ppt, &n_cigar, &poa_cigar); if (n_cigar) free(poa_cigar);
            poa_free_graph(graph, 1);
        }
    }
    free(bseq);
}

// int poa_main(const char *seq_fn, poa_para_t *ppt) { TODO
int poa_main(int seq_n, char (*seq)[100], poa_para_t *ppt){
    int i, j;
    poa_graph_t *graph = poa_init_graph(1);

    int **seq_node_ids = (int**)_err_malloc(seq_n * sizeof(int*));
    int *seq_node_ids_l = (int*)_err_malloc(seq_n * sizeof(int));
    uint8_t *bseq = (uint8_t*)_err_malloc(100 * sizeof(uint8_t)); int bseq_m = 100;
    // progressively partial order alignment
    for (i = 0; i < seq_n; ++i) {
        int seq_l = strlen(seq[i]);
#ifdef __DEBUG__
        printf("seq(%d): %s\n", seq_l, seq[i]);
#endif
        if (seq_l > bseq_m) {
            bseq_m = seq_l;
            bseq = (uint8_t*)_err_realloc(bseq, bseq_m * sizeof(uint8_t));
        }
        seq_node_ids[i] = (int*)_err_malloc(seq_l * sizeof(int));
        for (j = 0; j < seq_l; ++j) bseq[j] = nst_nt4_table[(int)(seq[i][j])];
        poa_cigar_t *poa_cigar=0; int n_cigar=0;
        poa_align_sequence_with_graph(graph, bseq, seq_l, ppt, &n_cigar, &poa_cigar);
        seq_node_ids_l[i] = 0;
        poa_add_graph_alignment(graph, bseq, seq_l, n_cigar, poa_cigar, seq_node_ids[i], seq_node_ids_l+i);
        //printf("%d %d %d\n", seq_l, graph->rank_n, seq_node_ids_l[i]);
        char poa_dot_fn[100]; sprintf(poa_dot_fn, "./dot_plot/poa_%d.dot", i);
        //poa_graph_visual(graph, poa_dot_fn);
        if (n_cigar) free(poa_cigar);
    }
    for (i = 0; i < seq_n; ++i) {
        printf("%s\n", seq[i]);
    }

    // generate consensus from graph
    poa_generate_consensus(graph);
    /*printf("consensus:\n");
    for (i = 0; i < graph->cons_l; ++i) {
        printf("%c", "ACGTN"[graph->cons_seq[i]]);
    } printf("\n");*/
    // generate multiple sequence alignment
    poa_generate_multiple_sequence_alingment(graph, seq_node_ids, seq_node_ids_l, seq_n, 1, stdout);
    //cons_to_seq_score(graph->cons_l, graph->cons_seq, seq_n, seq, ppt);

    poa_free_graph(graph, 1); free(bseq);
    for (i = 0; i < seq_n; ++i) free(seq_node_ids[i]); free(seq_node_ids); free(seq_node_ids_l);
    return 0;
}

int main(int argc, char **argv) {
    int seq_n = 6;
    char seq[100][100] = {
        //"ACGTAG",
        //"ACGAATAG",
        //"ATCAG",
        //"ATCAG",
        //"ATCAG",
        //"AGTAG",
        //"AGTAG",
        //"AGTCG",
        //"AGTCG",
        //"AGTCG",
        //"AGTCG",

        //"AAAATCGGCCCC",
        //"AAAATCGGCCCC"
        //"AAAATCGGCCCC",
        //"CCCCAGCATTTT",
        //"CCCCAGCATTTT",
        //"GGGGCTACTTTT",
        //"GGGGCTACTTTT"
        //"AACAT",
        //"AACT",
        //"AACT"
        //"AACAT",
        //"ACT",
        //"TTACGGAAGG"
        //"ACGGAA",
        //"ACGTAA"
        //"ACGGTTAA",
        //"ACTTGGAA",
        //"AGGGTTAA"
        //"AGTACCAA",
        //"AGGGCCAA"
        "CCGTAACCTTCATCGGATCACCGGAAAGGACCCGTAAATAGACCTGATTATCATCTACAT",
        "CATAAAAGAACGTAGGTCGCCCGTCCGTAACCTGTCGGATCACCGGAAAGGACCCGTAAAGTGATAATGAT",
        "ATAAAGGCAGTCGCTCTGTAAGCTGTCGATTCACCGGAAAGATGGCGTTACCACGTAAAGTGATAATGATTAT",
        "ATCAAAGAACGTGTAGCCTGTCCGTAATCTAGCGCATTTCACACGAGACCCGCGTAATGGG",
        "CGTAAATAGGTAATGATTATCATTACATATCACAACTAGGGCCGTATTAATCATGATATCATCA",
        "GTCGCTAGAGGCATCGTGAGTCGCTTCCGTACCGCAAGGATGACGAGTCACTTAAAGTGATAAT",
    };
    poa_para_t *ppt = poa_init_para();
    // parse argv TODO

    ppt->mat = (int8_t*)_err_malloc(ppt->m * ppt->m * sizeof(int8_t));
    gen_simple_mat(ppt->m, ppt->mat, ppt->match, ppt->mismatch);
    poa_main(seq_n, seq, ppt);
    poa_free_para(ppt);

    return 0;
}
