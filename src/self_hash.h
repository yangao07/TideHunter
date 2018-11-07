#ifndef SELF_HASH_H
#define SELF_HASH_H

#ifdef __cplusplus
extern "C" {
#endif

typedef int64_t hash_t;

typedef struct {
    int min_p, max_p;
    int start_i, end_i;
    int tot_n;
} hit_bucket_t;

//void radix_sort_hash(hash_t *beg, hash_t *end);

int hash_partition(char *seq, int seq_len, mini_tandem_para *mtp);

#ifdef __cplusplus
}
#endif

#endif
