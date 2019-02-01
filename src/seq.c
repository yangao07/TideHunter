#include "seq.h"
#include "utils.h"
#include "kseq.h"


#define nt_A 0
#define nt_C 1
#define nt_G 2
#define nt_T 3
#define nt_N 4

KSEQ_INIT(gzFile, gzread)

// ACGTN=>01234
unsigned char nst_nt4_table[256] = {
	0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 5 /*'-'*/, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

// ACGTN=>32104
unsigned char com_nst_nt4_table[256] = {
	3, 2, 1, 0,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4 /*'-'*/, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 3, 4, 2,  4, 4, 4, 1,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  0, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 3, 4, 2,  4, 4, 4, 1,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  0, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

//replace 'N' with 'G':           A  C  G  T  N->G
const uint8_t hash_nt4_table[6] = {0, 1, 2, 3, 2, 2};

char n_char[6] = {'A', 'C', 'G', 'T', 'N' };

uint32_t hash_key(uint8_t *bseq, int seq_len) {
    int i; uint32_t hash_key = 0;
    for (i = 0; i < seq_len; ++i) {
        // if (bseq[i] >= nt_N) err_printf("Error in bseq.\n");
        hash_key = (hash_key << 2) | bseq[i];
    }
    return hash_key;
}

uint32_t hash_shift_key(uint32_t pre_key, uint8_t *bseq, int pre_i, int cur_i, int k) {
    int i; uint32_t hash_key = pre_key;
    for (i = pre_i+1; i <= cur_i; ++i) {
        // if (bseq[i+k] >= nt_N) err_printf("Error in bseq.\n");
        hash_key = (hash_key << 2) | bseq[i];
    }
    return hash_key & ((1 << 2*k) - 1);
}

uint8_t *get_bseq(char *seq, int seq_len) {
    int i;
    uint8_t *bseq = (uint8_t*)_err_malloc(seq_len * sizeof(uint8_t));
    for (i = 0; i < seq_len; ++i) {
        // N(ambiguous base)
        // bseq[i] = hash_nt4_table[nst_nt4_table[(int)seq[i]]];
        bseq[i] = nst_nt4_table[(int)seq[i]];
    }
    return bseq;
}

char *get_rc_seq(char *seq, int seq_len) {
    int i;
    char *rc_seq = (char*)_err_malloc(sizeof(char) * (seq_len+1));
    for (i = 0; i < seq_len; ++i) {
        rc_seq[seq_len-i-1] = "ACGTN"[com_nst_nt4_table[(int)seq[i]]];
    } rc_seq[seq_len] = '\0';
    return rc_seq;
}
