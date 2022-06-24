#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "utils.h"
#include "edlib_align.h"
#include "edlib.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * para:
 * @query_seq/query_len: sequence and length of query
 * @target_seq/target_len: sequence and length of target
 * @start/end: start and end position on target sequence
 *
 * return:
 * edit distance
 */
EdlibEqualityPair additionalEqualities[5] = {
    {'a', 'A'}, 
    {'c', 'C'},
    {'g', 'G'},
    {'t', 'T'},
    {'n', 'N'}
};

int* edlibAlignmentToXID(const unsigned char* const alignment, const int alignmentLength) {
    // Maps move code from alignment to char in cigar.
    //                      0    1    2    3
    // moveCodeToChar[] : {'=', 'I', 'D', 'X'};
    int *xid = (int*) calloc(4, sizeof(int));

    int lastMove = -1;  // Char of last move. 0 if there was no previous move.
    int numOfSameMoves = 0, i;
    for (i = 0; i <= alignmentLength; i++) {
        // if new sequence of same moves started
        if (i == alignmentLength || (alignment[i] != lastMove && lastMove != -1)) {
            xid[lastMove] += numOfSameMoves;
            // If not at the end, start new sequence of moves.
            if (i < alignmentLength) {
                // Check if alignment has valid values.
                if (alignment[i] > 3) {
                    free(xid);
                    return 0;
                }
                numOfSameMoves = 0;
            }
        }
        if (i < alignmentLength) {
            lastMove = alignment[i];
            numOfSameMoves++;
        }
    }
    return xid;
}

int edlib_align_NW(char *query, int qlen, char *target, int tlen, int k) {
    int iden_n = 0;
    EdlibAlignResult result = edlibAlign(query, qlen, target, tlen, edlibNewAlignConfig(k, EDLIB_MODE_NW, EDLIB_TASK_PATH, additionalEqualities, 5));
    if (result.status == EDLIB_STATUS_OK) {
        int *xid = edlibAlignmentToXID(result.alignment, result.alignmentLength);
        if (xid != 0) {
             iden_n = xid[0]; 
             free(xid);
        }
    }
    edlibFreeAlignResult(result);
    return iden_n;
}

int edlib_align_HW(char *query, int qlen, char *target, int tlen, int *start, int *end, int k) {
    int ed = -1;
    EdlibAlignResult result = edlibAlign(query, qlen, target, tlen, edlibNewAlignConfig(k, EDLIB_MODE_HW, EDLIB_TASK_LOC, additionalEqualities, 5));
    if (result.status == EDLIB_STATUS_OK) {
        ed = result.editDistance;
        if (ed >= 0) {
            *start = result.startLocations[0];
            *end = result.endLocations[0];
        }
    }
    edlibFreeAlignResult(result);
    return ed;
}

// prefix global-alignment
int edlib_align_SHW(char *query, int qlen, char *target, int tlen, int k) {
    int iden_n = 0;
    EdlibAlignResult result = edlibAlign(query, qlen, target, tlen, edlibNewAlignConfig(k, EDLIB_MODE_SHW, EDLIB_TASK_PATH, additionalEqualities, 5));
    if (result.status == EDLIB_STATUS_OK) {
        int *xid = edlibAlignmentToXID(result.alignment, result.alignmentLength);
        if (xid != 0) {
            iden_n = xid[0]; 
            free(xid);
        }
    }
    edlibFreeAlignResult(result);
    return iden_n;
}

// suffix global-alignment
int edlib_align_UHW(char *query, int qlen, char *target, int tlen, int k) {
    char *rquery = (char*)_err_malloc(qlen * sizeof(char));
    char *rtarget = (char*)_err_malloc(tlen * sizeof(char));
    int i, iden_n=0;
    for (i = 0; i < qlen; ++i) rquery[i] = query[qlen-i-1];
    for (i = 0; i < tlen; ++i) rtarget[i] = target[tlen-i-1];
    iden_n = edlib_align_SHW(rquery, qlen, rtarget, tlen, k);

    free(rquery); free(rtarget);
    return iden_n;
}


#ifdef __cplusplus
}
#endif

#ifdef _EDLIB_ALIGN_MAIN
int main(void) {
    int start, end, ed, iden_n;
    char str1[100] = "TTTTTTT"; //"CCGTCG";
    char str2[100] = "AAAACCGTCGATAACAACCAAACCGGTCG";
    // iden_n = edlib_align_NW(str1, strlen(str1), str2, strlen(str2));
    // printf("iden_n : %d\n", iden_n);
    // iden_n = edlib_align_SHW(str1, strlen(str1), str2, strlen(str2));
    // printf("iden_n : %d\n", iden_n);
    ed = edlib_align_HW(str1, strlen(str1), str2, strlen(str2), &start, &end, -1);
    // printf("ed: %d, start: %d, end: %d\n", ed, start, end);
    return 0;
}
#endif
