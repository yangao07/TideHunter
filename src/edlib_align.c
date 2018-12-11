#include <stdio.h>
#include "edlib_align.h"
#include "edlib.h"

/*
 * para:
 * @query_seq/query_len: sequence and length of query
 * @target_seq/target_len: sequence and length of target
 * @start/end: start and end position on target sequence
 *
 * return:
 * edit distance
 */
#ifdef __cplusplus
extern "C" {
#endif

EdlibEqualityPair additionalEqualities[5] = {
    {'a', 'A'}, 
    {'c', 'C'},
    {'g', 'G'},
    {'t', 'T'},
    {'n', 'N'}
};

int edlib_align_HW(char *query, int qlen, char *target, int tlen, int *start, int *end) {
    int ed = -1;
    EdlibAlignResult result = edlibAlign(query, qlen, target, tlen, edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_LOC, additionalEqualities, 5));
    if (result.status == EDLIB_STATUS_OK) {
        ed = result.editDistance, *start = result.startLocations[0], *end = result.endLocations[0];
        printf("Edlib: %d, %d, %d\n", result.editDistance, result.alignmentLength, result.endLocations[0]);
    }
    edlibFreeAlignResult(result);
    return ed;
}
#ifdef __cplusplus
}
#endif

#ifdef _EDLIB_ALIGN_MAIN
int main(void) {
    int start, end, ed;
    ed = edlib_align_HW("ACGT", 4, "AAAAACGTCGATGCATGCTA", 15, &start, &end);
    return 0;
}
#endif
