#ifndef _MINI_TANDEM_H_
#define _MINI_TANDEM_H_
#include <zlib.h>
#include "abpoa.h"
#include "seq.h"

KSEQ_INIT(gzFile, gzread)

#define THREAD_N 4
#define CHUNK_READ_N 4096 // 1

#define KMER_SIZE 8  // kmer length
#define MAX_KMER_SIZE 16 // uint32_t
#define KMER_WSIZE 1 // window size
#define KMER_SSIZE 1 // step size
//#define KMER_MINM 1  // 1: minimizer

#define REP_RANGE -1 // -1: unlimited
#define MIN_COPY 2
#define MAX_DIV 0.25
#define MIN_PERIOD 2
#define DEF_MIN_PERIOD 30
#define MAX_PERIOD UINT32_MAX
#define DEF_MAX_PERIOD 10000 // UINT16_MAX
#define DEF_MAX_PERIOD_STR "10K" // UINT16_MAX
#define DEF_MIN_LEN 30
#define DEF_MIN_COV 0.0

#define MATCH 2
#define MISMATCH 4
#define GAP_OPEN1 4
#define GAP_EXT1 2
#define GAP_OPEN2 24
#define GAP_EXT2 1

#define FASTA_FMT 1
#define TAB_FMT 2
#define FASTQ_FMT 3
#define TAB_QUAL_FMT 4

#define ADA_MATCH_RAT 0.8

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int64_t k, w; int hpc; // k-mer length, window size, selected minimum m hash values // keep all k-mer when w == m
    int min_copy, min_cov; double max_div, div_exp, min_frac; // max allowed divergence
    int64_t min_p, max_p; // min/max period size
    int64_t max_range; // max range to find tandem repeat, -1 for no limit
    // msa parameters for abPOA
    int match, mismatch, gap_open1, gap_open2, gap_ext1, gap_ext2;
    // char *splint_fn, *splint_seq, *splint_rc_seq; int splint_len;
    float ada_match_rat;
    char *five_fn, *five_seq, *five_rc_seq; int five_len;
    char *three_fn, *three_seq, *three_rc_seq; int three_len;
    int out_fmt, min_len, only_unit, only_longest, only_full_length, single_copy; // only output the cons that spans the longest sequence
    FILE *cons_out, *detail_fp; //char detail_out[1024];
    int n_thread;
} mini_tandem_para;

typedef struct {
    // TODO kstring_t *output_str;
    seq_t *cons_seq;
    int cons_n, cons_m; 
    int *cons_start, *cons_end, *cons_len; double *copy_num, *ave_match;
    // int8_t *splint_rotated; // use cons_len to partition cons_seq when cons_n > 1
    int8_t *full_length;
    int **sub_pos; int *pos_n, *pos_m;
    int *cons_score;
} tandem_seq_t;

typedef struct {
    int tid;
    mini_tandem_para *mtp;
    int n_seqs; kseq_t *read_seq; tandem_seq_t *tseq;
    abpoa_t *ab; abpoa_para_t *abpt;
    // auxiliary memory for computing
    // int8_t **hit_array; int array_m; // for finding max hit
} thread_aux_t;


int mini_tandem(const char *read_fn, mini_tandem_para *mtp);
int tidehunter_core(kseq_t *read_seq, tandem_seq_t *tseq, mini_tandem_para *mtp, abpoa_t *ab, abpoa_para_t *abpt);

mini_tandem_para *mini_tandem_init_para(void);

void mini_tandem_free_para(mini_tandem_para *mtp);

#ifdef __cplusplus
}
#endif

#endif
