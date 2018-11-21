#include <stdlib.h>
#include <stdio.h>
#include "spoa_align.h"

int main(void) {
    char *seqs[10] = {"AACGT", "AAGT", "AACCGT", "ACGT"};
    char *cons_seq = (char*)malloc(10);
    spoa_msa(seqs, 4, cons_seq);
    printf("%s\n", cons_seq);
    free(cons_seq);
    return 0;
}
