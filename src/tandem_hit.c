#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "tandem_hit.h"
#include "utils.h"
#include "ksort.h"

#define sort_key_hash(x) (x)
KRADIX_SORT_INIT(hash, hash_t, sort_key_hash, 8)

// invertible integer hash function
static inline uint32_t inv_hash32(uint32_t key)
{
    // return key;
    key = ~key + (key << 15); // key = (key << 15) - key - 1;
    key = key ^ (key >> 12);
    key = key + (key << 2);
    key = key ^ (key >> 4);
    key = (key + (key << 3)) + (key << 11); // key = key * 2057;
    key = key ^ (key >> 16);
    return key;
}

static inline uint64_t inv_hash64(uint64_t key, uint64_t mask)
{
	key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
	key = key ^ key >> 24;
	key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
	key = key ^ key >> 14;
	key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
	key = key ^ key >> 28;
	key = (key + (key << 31)) & mask;
	return key;
}

// hash: kmer-key | rightmost-pos
int direct_hash(uint8_t *bseq, int seq_len, int k, int use_hpc, hash_t *h) {
    int c, l, hi=0, pos; uint32_t key, mask = ((uint64_t)1 << 2*k) - 1;

    for (key = l = pos = 0; pos < seq_len; ++pos) {
        c = bseq[pos];
        if (c >= 4) {
            key = 0; l = 0;
            continue;
        }
        if (use_hpc) while (pos+1 < seq_len && bseq[pos+1] == c) ++pos;
        key = key << 2 | c;
        if (++l >= k) { // get a kmer
            key &= mask; h[hi++] = ((hash_t)key << 32) | pos;
        }
    }
#ifdef __DEBUG__
    for (c = 0; c < hi; ++c) printf("%u: %u\n", (uint32_t)h[c], (uint32_t)(h[c]>>32));
#endif
    return hi;
}

typedef struct { uint32_t x, y; } u64_t;
#define _set_minimizer(minimizer, buf, i, min_buf_pos, equal_min_pos, equal_n) {    \
    if (buf[i].x < minimizer.x) {   \
        minimizer = buf[i]; \
        min_buf_pos = i;    \
        equal_n = 0;    \
    } else if (buf[i].x == minimizer.x) {   \
        equal_min_pos[equal_n++] = buf[i].y;    \
        min_buf_pos = i;    \
    }   \
}

typedef struct { // a simplified version of kdq
	int front, count;
	int a[32];
} tiny_queue_t;

static inline void tq_push(tiny_queue_t *q, int x)
{
	q->a[((q->count++) + q->front) & 0x1f] = x;
}

static inline int tq_shift(tiny_queue_t *q)
{
	int x;
	if (q->count == 0) return -1;
	x = q->a[q->front++];
	q->front &= 0x1f;
	--q->count;
	return x;
}

// para:
// bseq:    seq integer
// seq_len: seq length
// k:       kmer length
// w:       window size (s>=w, non-overlap window)
// use_hpc: use homopolymer-compressed kmer
// h:       kmer-key | rightmost-pos
int minimizer_hash(uint8_t *bseq, int seq_len, int k, int w, int use_hpc, hash_t *h) {
	int i, j, l, c, hi=0, kmer_span=0, buf_pos, min_pos;
    uint32_t key=0, mask = (1ULL << 2*k) - 1;
    u64_t buf[256], min = {UINT32_MAX, UINT32_MAX};
    tiny_queue_t tq; memset(&tq, 0, sizeof(tiny_queue_t));

    for (i = l = buf_pos = min_pos = 0; i < seq_len; ++i) {
        c = bseq[i];
        u64_t info = {UINT32_MAX, UINT32_MAX};
        if (c < 4 ) {
            if (use_hpc) {
                int skip_len = 1;
                if (i + 1 < seq_len && bseq[i+1] == c) {
                    for (skip_len = 2; i + skip_len < seq_len; ++skip_len)
                        if (bseq[i+skip_len] != c) break;
                    i += skip_len - 1;
                }
                tq_push(&tq, skip_len);
                kmer_span += skip_len;
                if (tq.count > k) kmer_span -= tq_shift(&tq);
            } else kmer_span = l+1 < k ? l+1 : k;
            key = (key << 2 | c) & mask;
            ++l;
            if (l >= k && kmer_span < 256) {
                info.x = key; info.y = i;
            }
        } else l = 0, tq.count = tq.front = 0, kmer_span = 0, key = 0;
        buf[buf_pos] = info;

        if (l == w+k-1 && min.x != UINT32_MAX) {
            for (j = buf_pos+1; j < w; ++j)
                if (min.x == buf[j].x && buf[j].y != min.y) h[hi++] = ((hash_t)buf[j].x << 32) | buf[j].y; 
            for (j = 0; j < buf_pos; ++j)
                if (min.x == buf[j].x && buf[j].y != min.y) h[hi++] = ((hash_t)buf[j].x << 32) | buf[j].y; 
        }

        if (info.x <= min.x) {
            if (l >= w+k && min.x != UINT32_MAX) h[hi++] = ((hash_t)min.x << 32) | min.y;
            min = info, min_pos = buf_pos;
        } else if (buf_pos == min_pos) {
            if (l >= w+k-1 && min.x != UINT32_MAX)  h[hi++] = ((hash_t)min.x << 32) | min.y;
            for (j = buf_pos+1, min.x = UINT32_MAX; j < w; ++j)
                if (min.x >= buf[j].x) min = buf[j], min_pos = j;
            for (j = 0; j <= buf_pos; ++j)
                if (min.x >= buf[j].x) min = buf[j], min_pos = j;
            if (l >= w+k-1 && min.x != UINT32_MAX) {
                for (j = buf_pos+1; j < w; ++j)
                    if (min.x == buf[j].x && min.y != buf[j].y) h[hi++] = ((hash_t)buf[j].x << 32) | buf[j].y; 
                for (j = 0; j <= buf_pos; ++j) 
                    if (min.x == buf[j].x && min.y != buf[j].y) h[hi++] = ((hash_t)buf[j].x << 32) | buf[j].y; 
            }
        }
        if (++buf_pos == w) buf_pos = 0;
    }
    if (min.x != UINT32_MAX) h[hi++] = ((hash_t)min.x << 32) | min.y;
#ifdef __DEBUG__
    for (i = 0; i < hi; ++i) printf("%u: %u\n", (uint32_t)h[i], (uint32_t)(h[i]>>32));
#endif
    return hi;

}

// h: kmer-key | rightmost-pos
int build_kmer_hash(uint8_t *bseq, int seq_len, mini_tandem_para *mtp, hash_t *h) {
    if (mtp->w > 1) { // do min hash
        return minimizer_hash(bseq, seq_len, mtp->k, mtp->w, mtp->hpc, h);
    } else { // do direct hash
        return direct_hash(bseq, seq_len, mtp->k, mtp->hpc, h);
    }
    return 0;
}

// TODO use multiple-hit seeds that connects multiple-repeats (n>2)
// discard too close hit, search further hit
int collect_hash_hit(hash_t *h, int hn, uint32_t min_p, uint32_t max_p, hash_t **hash_hit) {
    radix_sort_hash(h, h + hn); // sort h by hash values

    int i, n; int hit_n = 0;

    // calculate total hits number
    for (i = 1, n = 1; i < hn; ++i) {
        if (h[i] >> 32 != h[i-1] >> 32) {
            hit_n += (n-1); // use (n-1): only collect the adjacent hit
            n = 1;
        } else ++n;
    }
    hit_n += (n-1);

    // hash_hit: end:32 | period:32 
    *hash_hit = (hash_t*)_err_malloc(hit_n * sizeof(hash_t));
    int start_i, j, k, hi; uint32_t p;
    for (start_i = 0, hi = 0, i = 1, n = 1; i < hn; ++i) {
        if ((h[i] >> 32) != (h[i-1] >> 32)) {
            if (n > 1) {
                for (j = 1; j < n; ++j) {
                    for (k = j-1; k >= 0; --k) {
                        p = h[start_i+j] - h[start_i+k]; 
                        if (p >= min_p) break; // collect hit whose period >= p
                    }
                    if (p >= min_p && p <= max_p) {
                        (*hash_hit)[hi] = _set_hash_hit(h[start_i+j] & _32mask, p);
                        ++hi;
#ifdef __DEBUG__
                        printf("p: %u, end: %lu, start: %lu\n", p, h[start_i+j] & _32mask, h[start_i+j-1] & _32mask);
#endif
                    } else --hit_n;
                }
            }
            start_i = i, n = 1;
        } else ++n;
    }
    if (n > 1) {
        for (j = 1; j < n; ++j) {
            for (k = j-1; k >= 0; --k) {
                p = h[start_i+j] - h[start_i+k]; 
                    if (p >= min_p) break; // collect hit whose period >= p
            }
            if (p >= min_p && p <= max_p) {
                (*hash_hit)[hi] = _set_hash_hit(h[start_i+j] & _32mask, p);
                ++hi;
#ifdef __DEBUG__
                printf("p: %u, end: %lu, start: %lu\n", p, h[start_i+j] & _32mask, h[start_i+j-1] & _32mask);
#endif
            } else --hit_n;
        }
    }
    radix_sort_hash(*hash_hit, (*hash_hit) + hit_n); // sort hash_hit by end
    return hit_n;
}

int collect_tandem_repeat_hit(uint8_t *bseq, int seq_len, mini_tandem_para *mtp, hash_t **hit_h) {
    // generate hash value for each k-mer
    int hn; hash_t *h = (hash_t*)_err_malloc((seq_len - mtp->w) * sizeof(hash_t));
    if ((hn = build_kmer_hash(bseq, seq_len, mtp, h)) == 0) return 0;
#ifdef __DEBUG__
    printf("hash seed number: %d\n", hn);
#endif
    // collect hash hits
    int hit_n = collect_hash_hit(h, hn, (uint32_t)mtp->min_p, (uint32_t)mtp->max_p, hit_h); free(h);
    return hit_n;
}
