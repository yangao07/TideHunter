#ifndef SEQ_H
#define SEQ_H
#include <stdint.h>

extern unsigned char nst_nt4_table[256];
extern unsigned char com_nst_nt4_table[256];
extern const int8_t hash_nt4_table[6];
extern char n_char[6];

uint32_t hash_key(int8_t *bseq, int seq_len);
uint32_t hash_shift_key(uint32_t pre_key, int8_t *bseq, int pre_i, int cur_i, int k);
int get_bseq(char *seq, int seq_len, int8_t *bseq);

#endif
