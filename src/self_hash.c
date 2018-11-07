#include <stdio.h>
#include <stdlib.h>
#include "mini_tandem.h"
#include "self_hash.h"
#include "seq.h"
#include "utils.h"
#include "ksort.h"

#define sort_key_64(x) (x)
KRADIX_SORT_INIT(hash, uint64_t, sort_key_64, 8)

// 0. build Hash index, generate coordinate of all the hits (x,y)
// 1. chaining by partition based on (y-x) values
// 2. keep top N chains, (calcuate density of each chain: tot_N_hits / tot_N_kmer)
// 3. call consensus with each chain
// 4. polish consensus result

int overlap_hash(int8_t *bseq, int seq_len, int k, int w, hash_t *h) {
    int i, hi=0, pos, key;

    pos = 0, key = hash_key(bseq, k);
    h[hi++] = ((hash_t)key << 32) | pos;
    for (i = w; i <= seq_len - k; i += w) {
        pos = i;
        key = hash_shift_key(key, bseq+i-w, 0, w, k); 
        h[hi++] = ((hash_t)key << 32) | pos;
    }
    return hi;
}

int non_overlap_hash(int8_t *bseq, int seq_len, int k, int w, hash_t *h) {
    int i, hi=0, pos, key;

    for (i = 0; i <= seq_len - k; i += w) {
        pos = i;
        key = hash_key(bseq+i, k); 
        h[hi++] = ((hash_t)key << 32) | pos;
    }
    return hi;
}
// return :
// dict { k-mer_int : [pos1, pos2, ... posN] }
int direct_hash(int8_t *bseq, int seq_len, int k, int w, hash_t *h) {
    if (seq_len < w || seq_len < k) return 0;

    if (k > w) return overlap_hash(bseq, seq_len, k, w, h);
    else return non_overlap_hash(bseq, seq_len, k, w, h);
}

//int minimizer_hash(int8_t *bseq, int seq_len, int k, int w) {
//}

int collect_hash_hit(hash_t *h, int hn, hash_t **hit_h) {
    radix_sort_hash(h, h + hn);

    int i, n, n_keys=0;
    int hit_n = 0;

    // for each key, generate C(n,2) hits
    // calculate total hits number
    for (i = 1, n = 1; i < hn; ++i) {
        if (h[i] >> 32 != h[i-1] >> 32) {
            ++n_keys;
            if (n > 1) hit_n += (n * (n-1)) >> 1;
            n = 1;
        } else ++n;
    }
    if (n > 1) hit_n += n * (n-1);

    // generate hits coordinate
    *hit_h = (hash_t*)_err_malloc(hit_n * sizeof(hash_t));
    int start_i = 0, j, k, hi=0;
    for (i = 1, n = 1; i < hn; ++i) {
        if ((h[i] >> 32) != (h[i-1] >> 32)) {
            if (n > 1) {
                for (j = 1; j < n; ++j) {
                    for (k = start_i; k < i-j; ++k) {
                        (*hit_h)[hi++] = ((h[k+j] - h[k]) << 32) | (h[k+j] & 0xffffffff);
                    }
                }
            }
            start_i = i, n = 1;
        } else ++n;
    }
    if (n > 1) {
        for (j = 1; j < n; ++j) {
            for (k = start_i; k < i-j; ++k) {
                (*hit_h)[hi++] = ((h[k+j] - h[k]) << 32) | (h[k+j] & 0xffffffff);
            }
        }
    }
    radix_sort_hash(*hit_h, *hit_h + hit_n);
    return hit_n;
}

int cluster_hash_hit(hash_t *hit_h, int hit_n, double sigma, int bucket_T) {
    hit_bucket_t *hb;
    int i, hb_i=0, hb_m=20;
    double delta_p;

    hb = (hit_bucket_t*)_err_malloc(hb_m * sizeof(hit_bucket_t));
    hb[hb_i].min_p = hb[hb_i].max_p = hit_h[0] >> 32; hb[hb_i].start_i = hb[hb_i].end_i = 0; hb[hb_i].tot_n = 1;
    for (i = 1; i < hit_n; ++i) {
        delta_p = (double)(hit_h[i] >> 32) - (hit_h[i-1] >> 32);
        // printf("delta: %.2f, %.2f\n", delta_p, sigma );
        if (delta_p <= sigma * (hit_h[i] >> 32)) { // update current bucket
            hb[hb_i].max_p = hit_h[i] >> 32;
            hb[hb_i].end_i = i;
            hb[hb_i].tot_n += 1;
        } else { // new delta, add new bucket
            if (hb[hb_i].tot_n < bucket_T) {
                hb[hb_i].min_p = hb[hb_i].max_p = hit_h[i] >> 32;
                hb[hb_i].start_i = hb[hb_i].end_i = i;
                hb[hb_i].tot_n = 1;
            } else {
                if ((hb_i+1) == hb_m) {
                    hb_m <<= 1;
                    hb = (hit_bucket_t*)_err_realloc(hb, hb_m * sizeof(hit_bucket_t));
                }
                ++hb_i;
                hb[hb_i].min_p = hb[hb_i].max_p = hit_h[i] >> 32;
                hb[hb_i].start_i = hb[hb_i].end_i = i;
                hb[hb_i].tot_n = 1;
            }
        }
    }
    for (i = 0; i <= hb_i; ++i) {
        printf("(%3d, %3d) => %4d\n", hb[i].min_p, hb[i].max_p, hb[i].tot_n);
    }
    free(hb);
    return hb_i+1;
}

int hash_partition(char *seq, int seq_len, mini_tandem_para *mtp) {
    int8_t *bseq = (int8_t*)_err_malloc(seq_len * sizeof(int8_t));
    get_bseq(seq, seq_len, bseq);
    hash_t *h = (hash_t*)_err_malloc(seq_len / mtp->w * sizeof(hash_t));
    hash_t *hit_h;
    int hn = direct_hash(bseq, seq_len, mtp->k, mtp->w, h);
    if (hn == 0) return 0;

    int hit_n = collect_hash_hit(h, hn, &hit_h); free(h);

    cluster_hash_hit(hit_h, hit_n, mtp->sigma, mtp->bucket_T);
    free(hit_h); free(bseq);
    return 0;
}
