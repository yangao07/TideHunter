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
int direct_hash(uint8_t *bseq, int seq_len, int k, int s, int use_hpc, hash_t *h) {
    int c, l, hi=0, pos; uint32_t key, mask = ((uint64_t)1 << 2*k) - 1;

    for (key = l = pos = 0; pos < seq_len; ++pos) {
        c = bseq[pos];
        if (use_hpc) while (pos+1 < seq_len && bseq[pos+1] == c) ++pos;
        key = key << 2 | c;
        if (++l == k) { // get a kmer
            key &= mask; h[hi++] = ((hash_t)key << 32) | pos;
            if (s > k) {
                l = 0;
                // skip s-k bps
                if (use_hpc) {
                    while (l < s-k) {
                        while (pos + 1 < seq_len && bseq[pos+1] == bseq[pos]) ++pos;
                        ++l; ++pos;
                    }
                    l = 0;
                } else pos += (s-k);
            } else l = k-s;
        }
    }
#ifdef __DEBUG__
    //for (c = 0; c < hi; ++c) printf("%u: %u\n", (uint32_t)h[c], (uint32_t)(h[c]>>32));
#endif
    return hi;
}

int overlap_direct_hash(uint8_t *bseq, int seq_len, int k, int s, int use_hpc, hash_t *h) {
    int c, l, hi=0, pos; uint32_t key, mask = ((uint64_t)1 << 2*k) - 1;

    for (key = pos = l = 0; pos < seq_len; ++pos) {
        c = bseq[pos];
        if (use_hpc) while (pos + 1 < seq_len && bseq[pos+1] == c) ++pos;
        key = key << 2 | c;
        if (++l == k) { // get a kmer
            key &= mask; h[hi++] = ((hash_t)key << 32) | pos;
            l = k-s;
        }
    }
#ifdef __DEBUG__
    //for (c = 0; c < hi; ++c) printf("%u: %u\n", (uint32_t)h[c], (uint32_t)(h[c]>>32));
#endif
    return hi;
}

// para:
// bseq:    seq integer
// seq_len: seq length
// k:       kmer length
// s:       window step
// w:       window size (s>=w, non-overlap window)
// use_hpc: use homopolymer-compressed kmer
// h:       kmer-key | rightmost-pos
// non_overlap: no need to use buffer to store previous kmer hash values
int non_overlap_minimizer_hash(uint8_t *bseq, int seq_len, int k, int s, int w, int use_hpc, hash_t *h) {
    int c, l, j, pos, min_pos, hi=0;
    uint32_t key, minimizer, hash_value, mask = ((uint64_t)1 << 2*k) - 1;
    int *equal_min_pos = (int*)_err_malloc(w * sizeof(int)), equal_n;

    minimizer = UINT32_MAX; equal_n = 0;
    for (key = l = pos = 0; pos < seq_len; ++pos) {
        c = bseq[pos];
        if (use_hpc) while (pos + 1 < seq_len && bseq[pos+1] == c) ++pos;
        key = key << 2 | c;

        if (++l >= k) { // get a kmer
            key &= mask; hash_value = inv_hash32(key);
            if (hash_value < minimizer) {
                min_pos = pos; minimizer = hash_value; equal_n = 0;
            } else if (hash_value == minimizer) {
                equal_min_pos[equal_n++] = pos;
            }
            if (l == k+w-1) { // reach a w
                h[hi++] = ((hash_t)minimizer << 32) | min_pos;
                for (j = 0; j < equal_n; ++j) { 
                    h[hi++] = ((hash_t)minimizer << 32) | equal_min_pos[j];
                }

                if (s >= l) { // total non-overlap
                    // skip s-l bps
                    if (use_hpc) {
                        while (l < s) {
                            while (pos + 1 < seq_len && bseq[pos+1] == bseq[pos] && l < s) ++pos;
                            ++l; ++pos;
                        }
                    } else pos += (s-l);
                    l = 0;
                } else l = l-s; // overlap with next first kmer
                minimizer = UINT32_MAX; equal_n = 0;
            }
        }
    }

#ifdef __DEBUG__
    //for (c = 0; c < hi; ++c) printf("%u: %u\n", (uint32_t)h[c], (uint32_t)(h[c]>>32));
#endif
    free(equal_min_pos);
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
// First w-s: |===w-s===|...
//            |next mini|...
// Next s:    |===w-s===|==s==|...
//            |==s==|next mini|...
// Next s:    |===w-s===|==s==|==s==|...
//            |==s==|==s==|next mini|...
int overlap_minimizer_hash(uint8_t *bseq, int seq_len, int k, int s, int w, int use_hpc, hash_t *h) {
    int i, l, c, hi = 0;
    int pos, buf_pos, min_buf_pos, last_min_buf_pos=-1, equal_n, *equal_min_pos = (int*)_err_malloc(w * sizeof(int));
    uint32_t key, mask = ((uint64_t)1 << 2*k) - 1; 
    u64_t buf[256], minimizer, tmp;

    key = 0; minimizer.x = UINT32_MAX; equal_n = 0; min_buf_pos = -2;
    // first w-s
    for (pos = l = buf_pos = 0; pos < seq_len; ++pos) {
        c = bseq[pos];
        if (use_hpc) while (pos + 1 < seq_len && bseq[pos+1] == c) ++pos;
        key = key << 2 | c;

        if (++l >= k) { // get a kmer
            key &= mask; buf[buf_pos].x = inv_hash32(key); buf[buf_pos].y = (uint32_t)pos;
            _set_minimizer(minimizer, buf, buf_pos, min_buf_pos, equal_min_pos, equal_n);
            ++buf_pos;
            if (l == w-s+k-1) { // get w-s kmers
                ++pos;
                break;
            }
        }
    }
    // with a minimizer(,equal_min_pos, equal_n) in first w-s, calculate general minimizer with additional s
    // [w-s ~ w], [w ~ w+s], [w+s ~ w+2s] ...
    for (l = 0; pos < seq_len; ++pos) {
        c = bseq[pos];
        if (use_hpc) while (pos + 1 < seq_len && bseq[pos+1] == c) ++pos;
        key = (key << 2 | c) & mask; tmp.x = inv_hash32(key); tmp.y = (uint32_t)pos; buf[buf_pos] = tmp;

        if (tmp.x < minimizer.x) { // new minimizer
            if (min_buf_pos == last_min_buf_pos) {// store last old minimizer
                h[hi++] = ((hash_t)minimizer.x << 32) | minimizer.y;
            } minimizer = tmp; min_buf_pos = buf_pos; equal_n = 0;
        } else if (tmp.x == minimizer.x) {
            equal_min_pos[equal_n++] = pos; min_buf_pos = buf_pos;
        }
        if (++l == s) { // reach a new w
            // collect minimizer of w-s for next w
            if ((min_buf_pos <= buf_pos && min_buf_pos > buf_pos - w + s) || min_buf_pos > buf_pos + s) { // current minimizer is in next w
                if (equal_n > 0) {
                    h[hi++] = ((hash_t)minimizer.x << 32) | minimizer.y;
                    for (i = 0; i < equal_n-1; ++i) {
                        h[hi++] = ((hash_t)minimizer.x << 32) | equal_min_pos[i];
                    }
                    minimizer.y = equal_min_pos[equal_n-1]; equal_n = 0;
                }
                last_min_buf_pos = min_buf_pos;
            } else { // out of next w
                h[hi++] = ((hash_t)minimizer.x << 32) | minimizer.y;
                for (i = 0; i < equal_n; ++i) {
                    h[hi++] = ((hash_t)minimizer.x << 32) | equal_min_pos[i];
                }

                int new_start = buf_pos + s + 1;
                if (new_start < w) {
                    minimizer = buf[new_start]; min_buf_pos = new_start; equal_n = 0;
                    for (i = new_start+1; i < w; ++i) {
                        _set_minimizer(minimizer, buf, i, min_buf_pos, equal_min_pos, equal_n);
                    }
                    for (i = 0; i <= buf_pos; ++i) {
                        _set_minimizer(minimizer, buf, i, min_buf_pos, equal_min_pos, equal_n);
                    }
                } else {
                    new_start -= w;
                    minimizer = buf[new_start]; min_buf_pos = new_start; equal_n = 0;
                    for (i = new_start+1; i <= buf_pos; ++i) {
                        _set_minimizer(minimizer, buf, i, min_buf_pos, equal_min_pos, equal_n);
                    }
                }
                last_min_buf_pos = -1;
            }
            l = 0;
        }
        if (++buf_pos == w) buf_pos = 0;
    }

    h[hi++] = ((hash_t)minimizer.x << 32) | minimizer.y;
    for (i = 0; i < equal_n; ++i) { 
        h[hi++] = ((hash_t)minimizer.x << 32) | equal_min_pos[i];
    }

#ifdef __DEBUG__
    //for (c = 0; c < hi; ++c) printf("%u: %u\n", (uint32_t)h[c], (uint32_t)(h[c]>>32));
#endif
    free(equal_min_pos);
    return hi;
}

// TODO non overlap, m minmimal kmers
// m:       collect m minimal kmers (m<w, m==1 for minimizer_hash)
/*int min_hash(uint8_t *bseq, int seq_len, int k, int s, int w, int m, int use_hpc, hash_t *h) {
    return 0;
}*/

// TODO skip N(ambiguous base)
// h: kmer-key | rightmost-pos
int build_kmer_hash(uint8_t *bseq, int seq_len, mini_tandem_para *mtp, hash_t *h) {
    if (mtp->w > 1) { // do min hash
        if (mtp->m == 1) {
            if (mtp->s >= mtp->w) return non_overlap_minimizer_hash(bseq, seq_len, mtp->k, mtp->s, mtp->w, mtp->hpc, h);
            else return overlap_minimizer_hash(bseq, seq_len, mtp->k, mtp->s, mtp->w, mtp->hpc, h);
        } // else return min_hash(bseq, seq_len, mtp->k, mtp->s, mtp->w, mtp->m, mtp->hpc, h);
    } else { // do direct hash
        return direct_hash(bseq, seq_len, mtp->k, mtp->s, mtp->hpc, h);
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
