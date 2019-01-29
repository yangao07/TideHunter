#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mini_tandem.h"
#include "self_chain.h"
#include "edlib_align.h"
#include "abpoa_align.h"
#include "ksw2_align.h"
#include "spoa_align.h"
#include "seq.h"
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

static const char LogTable256[256] = {
#define LT(n) n, n, n, n, n, n, n, n, n, n, n, n, n, n, n, n
    -1, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3,
    LT(4), LT(5), LT(5), LT(6), LT(6), LT(6), LT(6),
    LT(7), LT(7), LT(7), LT(7), LT(7), LT(7), LT(7), LT(7)
};

static inline int ilog2_32(uint32_t v)
{
    uint32_t t, tt;
    if ((tt = v>>16)) return (t = tt>>8) ? 24 + LogTable256[t] : 16 + LogTable256[tt];
    return (t = v>>8) ? 8 + LogTable256[t] : LogTable256[v];
}

// hash: kmer-key | rightmost-pos
int direct_hash(uint8_t *bseq, int seq_len, int k, int s, int use_hpc, hash_t *h) {
    int c, l, hi=0, pos, key; uint32_t mask = (1 << 2*k) - 1;

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
    int c, l, hi=0, pos; uint32_t key, mask = (1 << 2*k) - 1;

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
    uint32_t key, minimizer, hash_value, mask = (1 << 2*k) - 1;
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
    uint32_t key, mask = (1 << 2*k) - 1; 
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

// TODO remove mem_l, it is unnecessary
// TODO use multiple-hit seeds that connects multiple-repeats (n>2)
// discard too close hit, search further hit
// 0. h: hash_value | position
// 1. hash_hit1: period | end
// 2. hit_h2: end | period | K ; K : MEM hit's length
// 3. seed_ids[hit_end] = seed_id
// collect cumulative kmer and hits number for each pos (for estimate of n and e)
int collect_mem_hash_hit(hash_t *h, int hn, int k, int min_p, int max_p, hash_t **hash_hit, int *seed_ids, int *seed_n) {
    radix_sort_hash(h, h + hn); // sort h by hash values

    int i, n; int hit_n = 0; *seed_n = 0;

    // calculate total hits number
    for (i = 1, n = 1; i < hn; ++i) {
        if (h[i] >> 32 != h[i-1] >> 32) {
            hit_n += (n-1); // use (n-1): only collect the adjacent hit
            n = 1;
        } else ++n;
    }
    hit_n += (n-1);

    // generate hash_hit1: period | end
    hash_t *hash_hit1 = (hash_t*)_err_malloc(hit_n * sizeof(hash_t));
    int start_i, j, _k, hi; int32_t p;
    for (start_i = 0, hi = 0, i = 1, n = 1; i < hn; ++i) {
        if ((h[i] >> 32) != (h[i-1] >> 32)) {
            if (n > 1) {
                for (j = 1; j < n; ++j) {
                    for (_k = j-1; _k >= 0; --_k) {
                        p = h[start_i+j] - h[start_i+_k];
                        if (p >= min_p) break; // collect hit whose period >= p
                    }
                    if (p >= min_p && p <= max_p) {
                        hash_hit1[hi] = ((hash_t)p << 32) | (h[start_i+j] & _32mask);
                        ++hi;
                        seed_ids[h[start_i+j] & _32mask] = *seed_n;
                    } else --hit_n;
                    // printf("end: %lld, seed_id: %d\n", h[start_i+j] & _32mask, *seed_n);
                }
                ++(*seed_n);
            }
            start_i = i, n = 1;
        } else ++n;
    }
    if (n > 1) {
        for (j = 1; j < n; ++j) {
            for (_k = j-1; _k >= 0; --_k) {
                p = h[start_i+j] - h[start_i+_k];
                if (p >= min_p) break; // collect hit whose period >= p
            }
            if (p >= min_p && p <= max_p) {
                hash_hit1[hi] = ((hash_t)p << 32) | (h[start_i+j] & _32mask);
                ++hi;
                // printf("end: %lld, seed_id: %d\n", h[start_i+j] & _32mask, *seed_n);
                seed_ids[h[start_i+j] & _32mask] = *seed_n;
            } else --hit_n;
        }
        ++(*seed_n);
    }
    radix_sort_hash(hash_hit1, hash_hit1 + hit_n); // sort hash_hit1 by period


    // generate mem_hash_hit: end:32 | period:16 | K:16 ; K: MEM hit's length
    int mem_hit_n, mem_l;
    for (i = 1, mem_hit_n = 0, n = 0; i < hit_n; ++i) {
        if ((_get_hash_hit1_period(hash_hit1, i) != _get_hash_hit1_period(hash_hit1, i-1)) || (mem_l = _get_hash_hit1_end(hash_hit1, i) - _get_hash_hit1_end(hash_hit1, i-1)) > k) {         

            ++mem_hit_n;
            n = 0;
        } else { // (_get_hash_hit1_end(hash_hit1, i) == _get_hash_hit1_end(hash_hit1, i-1)+1)
            if (n + mem_l > _16mask) {
                n = mem_l - 1;
                ++mem_hit_n;
            } else n += mem_l;
        }
    }
    ++mem_hit_n;

    *hash_hit = (hash_t*)_err_malloc(mem_hit_n * sizeof(hash_t));
    for (hi = 0, i = 1, n = 0; i < hit_n; ++i) {
        if ((_get_hash_hit1_period(hash_hit1, i) != _get_hash_hit1_period(hash_hit1, i-1)) || (mem_l = _get_hash_hit1_end(hash_hit1, i) - _get_hash_hit1_end(hash_hit1, i-1)) > k) {         
            (*hash_hit)[hi++] = _set_mem_hash_hit(hash_hit1[i-1], _get_hash_hit1_period(hash_hit1, i-1), n);
            n = 0;
        } else { // (_get_hash_hit1_end(hash_hit1, i) == _get_hash_hit1_end(hash_hit1, i-1)+1)
            if ((n + mem_l) > _16mask) {
                (*hash_hit)[hi++] = _set_mem_hash_hit(hash_hit1[i-1], _get_hash_hit1_period(hash_hit1, i-1), n);
                n = mem_l-1;
            } else n += mem_l;
        }
    }
    (*hash_hit)[hi++] = _set_mem_hash_hit(hash_hit1[i-1], _get_hash_hit1_period(hash_hit1, i-1), n);
    radix_sort_hash((*hash_hit), (*hash_hit) + mem_hit_n); // sort mem_hash_hit by end
    free(hash_hit1);
    return mem_hit_n;
}

int collect_hash_hit(hash_t *h, int hn, int min_p, int max_p, hash_t **hash_hit, int *seed_ids, int *seed_n) {
    radix_sort_hash(h, h + hn); // sort h by hash values

    int i, n; int hit_n = 0; *seed_n = 0;

    // calculate total hits number
    for (i = 1, n = 1; i < hn; ++i) {
        if (h[i] >> 32 != h[i-1] >> 32) {
            hit_n += (n-1); // use (n-1): only collect the adjacent hit
            n = 1;
        } else ++n;
    }
    hit_n += (n-1);

    // generate hash_hit1: period | end
    // make K always be 0
    // generate hash_hit: end:32 | period:16 | K:16 ; K: MEM hit's length
    *hash_hit = (hash_t*)_err_malloc(hit_n * sizeof(hash_t));
    int start_i, j, k, hi; int32_t p;
    for (start_i = 0, hi = 0, i = 1, n = 1; i < hn; ++i) {
        if ((h[i] >> 32) != (h[i-1] >> 32)) {
            if (n > 1) {
                for (j = 1; j < n; ++j) {
                    for (k = j-1; k >= 0; --k) {
                        p = h[start_i+j] - h[start_i+k]; 
                        if (p >= min_p) break; // collect hit whose period >= p
                    }
                    if (p >= min_p && p <= max_p) {
                        (*hash_hit)[hi] = _set_mem_hash_hit(h[start_i+j] & _32mask, p, 0);
                        ++hi;
                        seed_ids[h[start_i+j] & _32mask] = *seed_n;
#ifdef __DEBUG__
                        printf("p: %d, end: %lld, start: %lld\n", p, h[start_i+j] & _32mask, h[start_i+j-1] & _32mask);
#endif
                    } else --hit_n;
                }
                ++(*seed_n);
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
                (*hash_hit)[hi] = _set_mem_hash_hit(h[start_i+j] & _32mask, p, 0);
                ++hi;
                // printf("end: %lld, seed_id: %d\n", h[start_i+j] & _32mask, *seed_n);
                seed_ids[h[start_i+j] & _32mask] = *seed_n;
#ifdef __DEBUG__
                printf("p: %d, end: %lld, start: %lld\n", p, h[start_i+j] & _32mask, h[start_i+j-1] & _32mask);
#endif
            } else --hit_n;
        }
        ++(*seed_n);
    }
    radix_sort_hash(*hash_hit, (*hash_hit) + hit_n); // sort hash_hit by end
    return hit_n;
}

int dp_score_cmp(const void *a, const void *b) {
    return (((dp_score_t*)b)->score - ((dp_score_t*)a)->score);
}

int triple_i_cmp(const void *a, const void *b) {
    return (*(triple_t*)a).x - (*(triple_t*)b).x;
}

// TODO use heap to only keep top N scores
int sort_dp_score(dp_t **dp, int *array_size, int tot_n, dp_score_t *score_rank) { 
    int i, j, k;
    k = 0;
    for (i = tot_n-1; i >= 0; --i) {
        for (j = 0; j < array_size[i]; ++j) {
            if (dp[i][j].score > 0)
                score_rank[k++] = (dp_score_t){i, j, dp[i][j].score};
        }
    }
    qsort(score_rank, k, sizeof(dp_score_t), dp_score_cmp);
    return k;
}

void reverse_chain(chain_t *ch) {
    int i, tmp_i, tmp_j;
    for (i = 0; i < ch->len >> 1; ++i) {
        tmp_i = ch->cell[i].i, tmp_j = ch->cell[i].j;
        ch->cell[i].i = ch->cell[ch->len-1-i].i, ch->cell[i].j = ch->cell[ch->len-1-i].j; 
        ch->cell[ch->len-1-i].i = tmp_i, ch->cell[ch->len-1-i].j = tmp_j; 
    }
}

int is_overlap_chain1(dp_t **dp, chain_t ch1, chain_t ch2) {
    int s1, e1, s2, e2; cell_t c1, c2;
    c1 = ch1.cell[0]; s1 = dp[c1.i][c1.j].start;
    c1 = ch1.cell[ch1.len-1]; e1 = dp[c1.i][c1.j].start;
    c2 = ch2.cell[0]; s2 = dp[c2.i][c2.j].start;
    c2 = ch2.cell[ch2.len-1]; e2 = dp[c2.i][c2.j].start;

    int min = MIN_OF_TWO(e1-s1, e2-s2);
    int ovlp_len = MIN_OF_TWO(e1, e2) - MAX_OF_TWO(s1, s2);
    if (ovlp_len / (min+0.0) >= 0.5) return 1;
    else return 0;
}

// compare score, keep higher score
int is_overlap_chain(dp_t **dp, chain_t *chain, int ch_n, int ch_i) {
    if (ch_n <= 0 || chain[ch_i].len <= 0) return 0;
    int i; cell_t c = chain[ch_i].cell[chain[ch_i].len-1]; int start = dp[c.i][c.j].start;
    for (i = ch_n-1; i >= 0; --i) {
        if (chain[i].len <= 0) continue;
        if (dp[chain[i].cell[chain[i].len-1].i][chain[i].cell[chain[i].len-1].j].end <= start) break;
        if (is_overlap_chain1(dp, chain[i], chain[ch_i])) {
            if (chain[i].score > chain[ch_i].score) return 1;
            else {
                chain[i].len = 0;
                return 0;
            }
       }
   }
   return 0;
}

// backtrack from (x, y)
int backtrack_dp(dp_t **dp, int x, int y, chain_t *chain, int ch_n, int seq_len) {
    if (dp[x][y].is_tracked) return 0;
    int score = dp[x][y].score;
    int cur_i = x, cur_j = y, pre_i, pre_j, chain_len = 0;
    chain_t *ch = chain+ch_n;
    while (1) {
        // chain_add_hit(cur_i, cur_j);
        dp[cur_i][cur_j].is_tracked = 1;
        ch->cell[chain_len++] = (cell_t){cur_i, cur_j};
        pre_i = dp[cur_i][cur_j].from_i;
        pre_j = dp[cur_i][cur_j].from_j;
        if (pre_i == -1) {
            break;
        } else if (dp[pre_i][pre_j].is_tracked) {
            score -= dp[pre_i][pre_j].score;
            break;
        }
        cur_i = pre_i; cur_j = pre_j;
    } 
    ch->len = chain_len; ch->score = score;
    // printf("start, end: %d, %d, (%d, %d)\n", dp[x][y].start, dp[x][y].end, dp[x][y].end-dp[x][y].start, ch->score);
    reverse_chain(ch);
    if (chain_len > 1 && is_overlap_chain(dp, chain, ch_n, ch_n) == 0) {
        return 1;
    } else return 0;
}

void init_dp(hash_t *hit_h, dp_t **dp, int *hash_index, int *size, int total_n, int k) {
    // dp[total_n][0] = (dp_t){-1, -1, 0, 0, 0, 0, 0};
    int i, j, hash_i;
    int start, end, period, mem_l;
    for (i = 0; i < total_n; ++i) {
        for (j = 0; j < size[i]; ++j) {
            hash_i = hash_index[i] + j;
            end = _get_mem_hash_hit_end(hit_h, hash_i);
            period = _get_mem_hash_hit_period(hit_h, hash_i);
            mem_l = _get_mem_hash_hit_meml(hit_h, hash_i);
            start = end - period;
            dp[i][j] = (dp_t){-1, -1, start, end, mem_l, 2*k+mem_l+0, 0};
        }
    }
}

static inline int get_period_penalty1(int p1, int p2) {
    // int fold, max_p, min_p;
    // max_p = MAX_OF_TWO(p1, p2), min_p = MIN_OF_TWO(p1, p2);
    // fold = max_p / min_p;
    // return MIN_OF_TWO(abs(max_p-fold*min_p) + fold - 1, abs(max_p-(fold+1)*min_p) + fold);
    int delta_period = abs(p1-p2);
    // return delta_period / 2;
    return delta_period * delta_period / 2;
    // int log_d = delta_period ? ilog2_32(delta_period) : 0;
    // return (int)(.01 * delta_period * k) + (log_d >> 1);
}

static inline int get_dis_penalty(int dis1, int dis2) {
    return ilog2_32(dis1 + dis2) / 2;
}


// tend to pick a chain that only consists of periods around the true size
// cur_score = pre_score + match_bases - gap_cost; gap_cost = func(delta_period) + func(dis)
// return value:
#define NO_CON   0 // no connection
#define REG_CON  1 // regular connection
#define SAME_CON 2 // same distance connection
#define OVL_CON  3 // overlapped connection
static inline int get_con_score(dp_t *cur_dp, dp_t *pre_dp, int k, int *con_score) {
    int cur_start = cur_dp->start, pre_start = pre_dp->start, cur_end = cur_dp->end, pre_end = pre_dp->end;
    int cur_period = cur_end-cur_start, pre_period = pre_end-pre_start;
    if (cur_start <= pre_start || cur_period >= pre_period * 1.8 || pre_period >= cur_period * 1.8) {
        return NO_CON;  // crossing hits
    }
    int matched_bases = MIN_OF_TWO(abs(cur_end - pre_end), k + cur_dp->mem_l) + MIN_OF_TWO(abs(cur_start-pre_start), k + cur_dp->mem_l);
    int delta_period = abs(cur_period-pre_period), dis1 = abs(cur_end - pre_end), dis2 = abs(cur_start - pre_start);
    int gap_cost = get_period_penalty1(cur_period, pre_period) + get_dis_penalty(dis1, dis2);
    *con_score = matched_bases - gap_cost;

    if (delta_period == 0) {
        if (matched_bases < 2 * k) return OVL_CON; // overlapped 
        else return SAME_CON;
    } else return REG_CON;
}

// if max potential overlap between new and existing chains >= 1/2 of existing, discard new chain
// chains are sorted by end
int is_in_chain(dp_t **dp, chain_t *ch, int *chain_idx, int ch_n, int cell_i, int cell_j) {
    int i, _i, cell_start = dp[cell_i][0].start, cell_end = dp[cell_i][cell_j].end;
    cell_t c1, c2;
    for (_i = 0; _i < ch_n; ++_i) {
        i = chain_idx[_i];
        if (ch[i].len <= 0) continue;
        c1 = ch[i].cell[0], c2 = ch[i].cell[ch[i].len-1];
        int chain_start = dp[c1.i][c1.j].start, chain_end = dp[c2.i][c2.j].end;

        if (chain_end < cell_start) return 0;
        else if (chain_start > cell_end) continue;
        else if (cell_end - chain_start >= (chain_end - chain_start) / 2) return 1;
        // else if (chain_start <= cell_end && cell_start <= chain_end) return 1;
    }
    return 0;
}

// sort by chain.end
void sort_chain(dp_t **dp, chain_t *chain, int *chain_idx, int ch_n) {
    if (ch_n < 2) return;
    int i, _i, j, _j, ch_end1, ch_end2; cell_t c1, c2;
    for (_i = 0; _i < ch_n-1; ++_i) {
        i = chain_idx[_i];
        if (chain[i].len <= 0) continue;
        c1 = chain[i].cell[chain[i].len-1];
        ch_end1 = dp[c1.i][c1.j].end;
        for (_j = _i+1; _j < ch_n; ++_j) {
            j = chain_idx[_j];
            if (chain[j].len <= 0) continue;
            c2 = chain[j].cell[chain[j].len-1];
            ch_end2 = dp[c2.i][c2.j].end;
            if (ch_end1 < ch_end2) { // switch ch1 and ch2
                chain_idx[_i] = j; chain_idx[_j] = i;
                ch_end1 = ch_end2;
            }
        }
    }
}

int copy_chain(chain_t *src_ch, int seq_len, int start_i, int end_i, chain_t **dest_ch, int *ch_n, int *ch_m) {
    if (start_i < 0 || end_i >= src_ch->len || end_i - start_i < 2) return 0;
    int i;
    if (*ch_n == *ch_m) {
        *ch_m <<= 1;
        *dest_ch = (chain_t*)_err_realloc(*dest_ch, *ch_m * sizeof(chain_t));
        for (i = *ch_n; i < *ch_m; ++i) { 
            (*dest_ch)[i].cell = (cell_t*)_err_malloc(seq_len * sizeof(cell_t));
            (*dest_ch)[i].len = 0;
        }
    }
    for (i = start_i; i <= end_i; ++i) (*dest_ch)[*ch_n].cell[i-start_i] = src_ch->cell[i];
    (*dest_ch)[(*ch_n)++].len = end_i - start_i + 1;
    return 0;
}

int split_chain(char *read_name, char *seq, uint8_t *bseq, int seq_len, dp_t **dp, chain_t *ch, chain_t **post_ch, int *post_ch_n, int *post_ch_m, mini_tandem_para *mtp) {
    copy_chain(ch, seq_len, 0, ch->len-1, post_ch, post_ch_n, post_ch_m); // do not split chain; split chain based on alignment
    return 0;
}

// TODO max allowed fold
int get_period_penalty(int P, triple_t *period, int p_n) {
    int i, pp = 0;
    for (i = 0; i < p_n; ++i) {
        pp += get_period_penalty1(P, period[i].x);
    }
    return pp;
}

// calculate period for each hit
void get_est_period(dp_t **dp, chain_t *ch) {
    int i;
    triple_t *t = (triple_t*)_err_malloc(ch->len * sizeof(triple_t));
    cell_t c; dp_t d;
    // collect period
    ch->max_period = INT32_MIN; ch->min_period = INT32_MAX;
    for (i = 0; i < ch->len; ++i) {
        c = ch->cell[i]; d = dp[c.i][c.j];
        t[i].x = d.end - d.start; t[i].y = d.start; t[i].z = i;
    }
    qsort(t, ch->len, sizeof(triple_t), triple_i_cmp);
    int last_p = -1, pp, min_pp = INT32_MAX;
    ch->est_period = 0;
    for (i = 0; i < ch->len; ++i) {
        if (t[i].x == last_p) continue;
        if (t[i].x > ch->max_period) ch->max_period = t[i].x;
        if (t[i].x < ch->min_period) ch->min_period = t[i].x;
        last_p = t[i].x;
        pp = get_period_penalty(t[i].x, t, ch->len);
#ifdef __DEBUG__
        printf("p: %d, penalty: %d\n", t[i].x, pp);
#endif
        if (pp < min_pp) {
            ch->est_period = t[i].x; ch->est_start = t[i].y; ch->est_ch_i = t[i].z;
            min_pp = pp;
        } 
    }
    free(t);
#ifdef __DEBUG__
    printf("Est_P: %d, %d, %d\n", ch->est_ch_i, ch->est_period, ch->est_start);
#endif
}

// hash_hit: hash table of mem hits
// TODO allocate DP matrix uniformly
int dp_chain(char *read_name, char *seq, uint8_t *bseq, int seq_len, hash_t *hash_hit, int hash_hit_n, mini_tandem_para *mtp, dp_t ***_dp, int *tot_N,chain_t **post_chain, int *post_ch_m) {
    int i, j, k, idx, kmer_k = mtp->k;
    // calculate DP matrix size, allocate DP matrix
    int tot_n = 1, *array_size, *hash_index;
    for (i = 1; i < hash_hit_n; ++i) {
        if (_get_mem_hash_hit_end(hash_hit, i)  != _get_mem_hash_hit_end(hash_hit, i-1))
            tot_n += 1;
    }
    *tot_N = tot_n;
    array_size = (int*)_err_malloc(sizeof(int) * tot_n);
    hash_index = (int*)_err_malloc(sizeof(int) * tot_n);
    *_dp = (dp_t**)_err_calloc((tot_n+1), sizeof(dp_t*));
    dp_t **dp = *_dp;
    dp[tot_n] = (dp_t*)_err_calloc(1, sizeof(dp_t));

    for (i = 1, j = 0, k = 1, idx = 0; i < hash_hit_n; ++i) {
        if (_get_mem_hash_hit_end(hash_hit, i)  != _get_mem_hash_hit_end(hash_hit, i-1)) {
            dp[j] = (dp_t*)_err_calloc(k, sizeof(dp_t));
            hash_index[j] = idx-k+1; array_size[j++] = k;
            k = 1;
        } else ++k;
        ++idx;
    }
    dp[j] = (dp_t*)_err_calloc(k, sizeof(dp_t));
    hash_index[j] = idx-k+1; array_size[j] = k;

    // initialize DP matrix
    // set (tot_n,0) as all cells' precurser
    init_dp(hash_hit, dp, hash_index, array_size, tot_n, kmer_k);

    // main DP process
    int cur_i, cur_j, pre_i, pre_j, max_pre_i, max_pre_j, con_score, score, max_score; 
    int con_res, iter_n, max_h; // max_h: number of meaningless iterations to stop DP
    dp_t *cur_dp, *pre_dp;
    for (cur_i = 1; cur_i < tot_n; ++cur_i) {
        for (cur_j = 0; cur_j < array_size[cur_i]; ++cur_j) {
            cur_dp = dp[cur_i]+cur_j;
            max_score = cur_dp->score;
            max_h = cur_dp->end - cur_dp->start; // TODO use z-drop like strategy
            iter_n = 0;
            for (pre_i = cur_i-1; pre_i >= 0; --pre_i) {
                if (dp[pre_i][0].end < cur_dp->start) goto UPDATE; // TODO max_h is unnecessary when this goto is used
                // if (dp[pre_i][0].end < cur_dp->start - (cur_dp->end-cur_dp->start)) goto UPDATE;
                int gt = 0;
                for (pre_j = 0; pre_j < array_size[pre_i]; ++pre_j) {
                    pre_dp = dp[pre_i]+pre_j;
                    con_res = get_con_score(cur_dp, pre_dp, kmer_k, &con_score);
                    if (con_res == NO_CON) continue;
                    score = dp[pre_i][pre_j].score + con_score;
                    if (score > max_score) {
                        max_score = score; max_pre_i = pre_i, max_pre_j = pre_j;
                        if (con_res == SAME_CON || con_res == OVL_CON) goto UPDATE;
                        gt = 1;
                    } else if (con_res == OVL_CON) {
                        goto UPDATE; 
                    }
                }
                if (gt) iter_n = 0;
                else if (++iter_n >= max_h) goto UPDATE; // only try h iterations
            }
UPDATE:
            if (max_score > cur_dp->score) {
                cur_dp->score = max_score; cur_dp->from_i = max_pre_i; cur_dp->from_j = max_pre_j;
            }
        }
    }

    // TODO backtrack, obtain top N chains, use max-heap
    dp_score_t *score_rank = (dp_score_t*)_err_malloc(hash_hit_n * sizeof(dp_score_t));
    int score_n = sort_dp_score(dp, array_size, tot_n, score_rank);
    int top_N = 1000, ch_n = 0, ch_m = top_N;
    chain_t *chain = (chain_t*)_err_malloc(top_N * sizeof(chain_t)); int *chain_idx = (int*)_err_malloc(sizeof(int) * top_N); // chain_idx[rank]: index in chain
    int post_ch_n = 0; *post_ch_m = top_N; *post_chain = (chain_t*)_err_malloc(top_N * sizeof(chain_t));
    for (i = 0; i < top_N; ++i) {
        chain[i].cell = (cell_t*)_err_malloc(tot_n * sizeof(cell_t)); chain[i].len = 0;
        (*post_chain)[i].cell = (cell_t*)_err_malloc(tot_n * sizeof(cell_t)); (*post_chain)[i].len = 0;
        chain_idx[i] = i;
    }
    for (i = ch_n = 0; i < score_n && ch_n < top_N; ++i) {
        if (is_in_chain(dp, chain, chain_idx, ch_n, score_rank[i].i, score_rank[i].j)) continue;
        if (backtrack_dp(dp, score_rank[i].i, score_rank[i].j, chain, ch_n, seq_len)) ++ch_n;
        sort_chain(dp, chain, chain_idx, ch_n);
    }
#ifdef __DEBUG__
    printf("ch_n: %d\n", ch_n);
    int _i;
    chain_t *ch = chain;
    for (_i = 0; _i < ch_n; ++_i) {
        i = chain_idx[_i];
        if (ch[i].len > 0) {
            int start_i = ch[i].cell[0].i, start_j = ch[i].cell[0].j, end_i = ch[i].cell[ch[i].len-1].i, end_j = ch[i].cell[ch[i].len-1].j;
            int from_i = dp[start_i][start_j].from_i, from_j = dp[start_i][start_j].from_j;
            j = 0;
            printf("\tchain: %d(%d): start: %d, end: %d, p: %d, mem_l: %d, score: %d\n", j+1, ch[i].cell[j].i, dp[ch[i].cell[j].i][ch[i].cell[j].j].start, dp[ch[i].cell[j].i][ch[i].cell[j].j].end, dp[ch[i].cell[j].i][ch[i].cell[j].j].end-dp[ch[i].cell[j].i][ch[i].cell[j].j].start, dp[ch[i].cell[j].i][ch[i].cell[j].j].mem_l, dp[ch[i].cell[j].i][ch[i].cell[j].j].score);
            for (j = 1; j < ch[i].len; ++j) {
                printf("\tchain: %d(%d): start: %d, end: %d, p: %d, mem_l: %d, score: %d, delta: %d\n", j+1, ch[i].cell[j].i, dp[ch[i].cell[j].i][ch[i].cell[j].j].start, dp[ch[i].cell[j].i][ch[i].cell[j].j].end, dp[ch[i].cell[j].i][ch[i].cell[j].j].end-dp[ch[i].cell[j].i][ch[i].cell[j].j].start, dp[ch[i].cell[j].i][ch[i].cell[j].j].mem_l, dp[ch[i].cell[j].i][ch[i].cell[j].j].score, dp[ch[i].cell[j].i][ch[i].cell[j].j].start- dp[ch[i].cell[j-1].i][ch[i].cell[j-1].j].start);
            }
        }
    }
#endif
    // post-process of chains
    post_ch_n = 0;
    for (i = 0; i < ch_n; ++i) {
        split_chain(read_name, seq, bseq, seq_len, dp, chain+i, post_chain, &post_ch_n, post_ch_m, mtp); // split chain in sparse region
    }
    for (i = 0; i < post_ch_n; ++i)
        get_est_period(dp, (*post_chain)+i);
    for (i = 0; i < ch_m; ++i) free(chain[i].cell); free(chain); free(chain_idx);
    free(array_size); free(hash_index); free(score_rank);
    return post_ch_n;
}

void write_tandem_cons_seq(tandem_seq_t *tseq, char *cons_seq, int cons_len, int start, int end, double copy_num, mini_tandem_para *mtp, int8_t splint_rotated) {
    if (cons_len < mtp->min_p || cons_len > mtp->max_p) return;
    if (mtp->only_longest && tseq->cons_n == 1) {
        if (end-start > tseq->cons_end[0]-tseq->cons_start[0]) {
            tseq->cons_n = 0; tseq->cons_seq->seq.l = 0;
        } else return;
    }
    if (tseq->cons_seq->seq.l + cons_len >= tseq->cons_seq->seq.m) {
        tseq->cons_seq->seq.m = tseq->cons_seq->seq.l + cons_len + 1;
        tseq->cons_seq->seq.s = (char*)_err_realloc(tseq->cons_seq->seq.s, tseq->cons_seq->seq.m * sizeof(char));
    }
    strcpy(tseq->cons_seq->seq.s + tseq->cons_seq->seq.l, cons_seq); tseq->cons_seq->seq.l += cons_len; 

    if (tseq->cons_n == tseq->cons_m) {
        tseq->cons_m <<= 1;
        tseq->cons_start = (int*)_err_realloc(tseq->cons_start, tseq->cons_m * sizeof(int));
        tseq->cons_end = (int*)_err_realloc(tseq->cons_end, tseq->cons_m * sizeof(int));
        tseq->copy_num = (double*)_err_realloc(tseq->copy_num, tseq->cons_m * sizeof(double));
        tseq->splint_rotated = (int8_t*)_err_realloc(tseq->splint_rotated, tseq->cons_m * sizeof(int8_t));
        tseq->cons_len = (int*)_err_realloc(tseq->cons_len, tseq->cons_m * sizeof(int));
        tseq->cons_score = (int*)_err_realloc(tseq->cons_score, tseq->cons_m * sizeof(int));
    }
    tseq->cons_start[tseq->cons_n] = start; tseq->cons_end[tseq->cons_n] = end; tseq->copy_num[tseq->cons_n] = copy_num;
    tseq->splint_rotated[tseq->cons_n] = splint_rotated;
    tseq->cons_len[tseq->cons_n] = cons_len; tseq->cons_score[tseq->cons_n] = 0; // TODO cons_score
    ++tseq->cons_n;
}

// pos: 0-base
int *partition_seqs_core(char *seq, int seq_len, int8_t *hit_array, int est_period, int chain_start, int chain_end, int *par_n) {
    int i, j;
    int *pos_array = (int*)_err_malloc(sizeof(int) * seq_len), *par_pos = (int*)_err_malloc(sizeof(int) * seq_len);
    int hit_n = 0;
    for (i = 0; i < seq_len; ++i) {
        if (hit_array[i]) {
            #ifdef __DEBUG__
            printf("%d, ", i);
            #endif
            pos_array[hit_n++] = i;
        }
    } 
#ifdef __DEBUG__
    printf("\n");
#endif
    // partition seq into period seperated seqs
    int par_i = 0, l, tot_len; 
    char *query_seq, *target_seq; int copy_num, ed, start, end, target_start;
    //  extend two ends
    l = est_period / 2;
    copy_num = (int)((double)(pos_array[0] - chain_start) / est_period + 0.5);
    if (copy_num >= 1) {
        for (j = copy_num; j >= 1; --j) {
            query_seq = seq + pos_array[0];
            target_start = pos_array[0] - est_period * j - (l<<0);
            // target_start = pos_array[0] - period[pos_array[0]] * j - (l<<1);
            if (target_start < 0) continue;
            target_seq = seq + target_start;
            ed = edlib_align_HW(query_seq, l, target_seq, l<<1, &start, &end);
            // ed = edlib_align_HW(query_seq, l, target_seq, l<<2, &start, &end);
            if (ed >= 0) par_pos[par_i++] = target_start + start;
        }
    }
    par_pos[par_i++] = pos_array[0];
    for (i = 0; i < hit_n-1; ++i) {
        // printf("%d: %d, %d, %d\n", pos_array[i], period[pos_array[i]], period[pos_array[i+1]], ave_p);
        copy_num = (int)((double)(pos_array[i+1] - pos_array[i]) / est_period + 0.5);
        if (copy_num > 1) { // multiple copies: semi-global alignment of prefix l-mer using edlib
            tot_len = pos_array[i+1] - pos_array[i];
            query_seq = seq + pos_array[i];
            for (j = 1; j < copy_num; ++j) {
                target_start = pos_array[i] + est_period * j - (l << 0);
                // target_start = pos_array[i] + ave_p * j - (l << 1);
                target_seq = seq + target_start;
                ed = edlib_align_HW(query_seq, l, target_seq, l<<1, &start, &end);
                // ed = edlib_align_HW(query_seq, l, target_seq, l<<2, &start, &end);
                if (ed < 0) { // no alignment result
                    par_pos[par_i++] = -1; // skip this copy
                } else {
                    par_pos[par_i++] = target_start + start;
                }
            }
        }
        par_pos[par_i++] = pos_array[i+1];
    }
    // extend par_pos on the right
    copy_num = (int)((double)(chain_end - pos_array[hit_n-1]) / est_period + 0.5);
    if (copy_num >= 1) {
        for (j = 1; j <= copy_num; ++j) {
            query_seq = seq + pos_array[hit_n-1];
            target_start = pos_array[hit_n-1] + est_period * j - (l<<0);
            // target_start = pos_array[hit_n-1] + period[pos_array[hit_n-1]] * j - (l<<1);
            if (target_start + (l<<1) >= seq_len) continue;
            // if (target_start + (l<<2) >= seq_len) continue;
            target_seq = seq + target_start;
            ed = edlib_align_HW(query_seq, l, target_seq, l<<1, &start, &end);
            // ed = edlib_align_HW(query_seq, l, target_seq, l<<2, &start, &end);
            if (ed >= 0) par_pos[par_i++] = target_start + start;
        }
    }
    free(pos_array);
    *par_n = par_i; 
    return par_pos;
}

// multi-hits in one period:
// very close distance of the same kmers causes a negetive score: two identical kmers in one period
static inline int get_max_hit(int8_t *hit_array, int *seed_ids, int seed_n, dp_t **dp, chain_t ch) {
    int i, j, start, end, seed_id, mem_l; dp_t dp_cell;
    int *seed_hits = (int*)_err_calloc(seed_n, sizeof(int));
    int *seed_last_pos = (int*)_err_malloc(seed_n * sizeof(int)); memset(seed_last_pos, -1, seed_n * sizeof(int));
    int max_period = ch.max_period, min_period = ch.min_period;

    for (i = ch.len-1; i >= 0; --i) {
        dp_cell = dp[ch.cell[i].i][ch.cell[i].j];
        start = dp_cell.start, end = dp_cell.end, mem_l = dp_cell.mem_l;
        for (j = 0; j <= mem_l; ++j) {
            seed_id = seed_ids[end-j];

            if (seed_last_pos[seed_id] == -1) { // first
                seed_hits[seed_id] += 2;
                hit_array[start-j] = 1;
                hit_array[end-j] = 1;
            } else {
                if (seed_last_pos[seed_id] == end-j) { // 
                    seed_hits[seed_id] += 1;
                    hit_array[start-j] = 1;
                } else if (seed_last_pos[seed_id] - (end-j) < min_period) {
                    seed_hits[seed_id] = INT32_MIN; // discard impure k-mers
                    hit_array[start-j] = 1;
                    hit_array[seed_last_pos[seed_id]] = 0;
                } else { 
                    seed_hits[seed_id] += 2;
                    hit_array[start-j] = 1;
                    hit_array[end-j] = 1;
                }
            }
            seed_last_pos[seed_id] = start-j;
            // printf("seed_id: %d => start: %d, end: %d\n", seed_id, start-j, end-j);
            // seed_hits[seed_id] += (1-hit_array[seed_id][start-j] + 1-hit_array[seed_id][end-j]);
        }
    }
    
    int max_hit_n = 2, max_id = -1; // max_hit_n >= 3
    for (i = 0; i < seed_n; ++i) {
#ifdef __DEBUG__
        //printf("%d: %d\n", i, seed_hits[i]);
#endif
        if (seed_hits[i] > max_hit_n) {
            max_hit_n = seed_hits[i];
            max_id = i;
       }
    }
#ifdef __DEBUG__
    printf("max_i: %d, max_hit_n: %d\n", max_id, max_hit_n);
#endif
    if (max_id != -1) { // set hit_array as max_id's hit_array
        for (i = ch.len-1; i >= 0; --i) {
            dp_cell = dp[ch.cell[i].i][ch.cell[i].j];
            start = dp_cell.start, end = dp_cell.end, mem_l = dp_cell.mem_l;
            for (j = 0; j <= mem_l; ++j) {
                seed_id = seed_ids[end-j];
                if (seed_id == max_id) {
                    hit_array[end-j] = 1;
                    hit_array[start-j] = 1;
                } else {
                    hit_array[end-j] = 0;
                    hit_array[start-j] = 0;
                }
            }
        }
    }
    free(seed_hits); free(seed_last_pos);
    return max_id;
}

int *partition_seqs(char *seq, int seq_len, dp_t **dp, int *seed_ids, int seed_n, chain_t ch, int *par_n) {
    int *par_pos = NULL, max_id; *par_n = 0;
    int8_t *_hit_array = (int8_t*)_err_calloc(seq_len, sizeof(int8_t));
    int chain_start = dp[ch.cell[0].i][ch.cell[0].j].start, chain_end = dp[ch.cell[ch.len-1].i][ch.cell[ch.len-1].j].end;
    // fill hit array
    if ((max_id = get_max_hit(_hit_array, seed_ids, seed_n, dp, ch)) >= 0)
        par_pos = partition_seqs_core(seq, seq_len, _hit_array, ch.est_period, chain_start, chain_end, par_n);
    free(_hit_array);
    return par_pos;
}

// TODO when (s2) >> (e1), i.e., anchor are too sparse
// TODO Extend by non-full copy, more than one copy
// start with est_start anchor?
int *get_partition_pos_with_global_alignment(uint8_t *bseq, int seq_len, dp_t **dp, chain_t ch, mini_tandem_para *mtp, int *par_n) {
    int est_ch_i = ch.est_ch_i, est_start = ch.est_start, est_period = ch.est_period;
    int first_start = dp[ch.cell[0].i][ch.cell[0].j].start, first_end = dp[ch.cell[0].i][ch.cell[0].j].end;
    int last_start = dp[ch.cell[ch.len-1].i][ch.cell[ch.len-1].j].start, last_end = dp[ch.cell[ch.len-1].i][ch.cell[ch.len-1].j].end;
    int *par_pos = (int*)_err_malloc(seq_len * sizeof(int)); *par_n = 0;
    int i, j, k=mtp->k, ch_i1, ch_i2, s1, e1, s2, e2, iden_n, s, e;
    int n_cigar; uint32_t *cigar;
    cell_t c; dp_t d;
    // find non-overlapping adjacent tandem repeat hits
    // global alignment of [s1,s2] and [e1, e2]
    // find best partition position with backtrack
    ch_i2 = est_ch_i, s2 = est_start, e2 = est_start + est_period;
    while (ch_i2 > 0 && s2 >= first_end) {
        for (i = ch_i2 - 1; i >= 0; --i) {
            c = ch.cell[i]; d = dp[c.i][c.j];
            s1 = d.start; e1 = d.end;
            if (e1 == s2) {
                par_pos[(*par_n)++] = s1;
                ch_i2 = i, s2 = s1, e2 = e1;
                break;
            } else if (e1 < s2) { // non-overlapping adjacent tandem repeat hits
                // do global alignment
                iden_n = ksw2_global_with_cigar(bseq+e1-k+1, e2-e1+k, bseq+s1-k+1, s2-s1+k, &n_cigar, &cigar);
                #ifdef __DEBUG__
                printf("iden_n: %d (%d,%d), (%d,%d)\n", iden_n, e2-e1+k, s2-s1+k, s1, s2);
                #endif
                if (iden_n >= MIN_OF_TWO(s2-s1+k, e2-e1+k) * (1-mtp->max_div)) { // extend partition
                    s = s2 - backtrack_left_end(n_cigar, cigar, e2-e1+k, s2-s1+k, e2-s2);
                    par_pos[(*par_n)++] = s; 
                    ch_i2 = i+1; e2 = s2, s2 = s;
                } else { // skip this anchor
                    par_pos[(*par_n)++] = -1; // indicate a separation flag
                    par_pos[(*par_n)++] = e1; 
                    par_pos[(*par_n)++] = s1; 
                    ch_i2 = i; s2 = s1; e2 = e1;
                }
                if (cigar) 
                    free(cigar);
                break;
            }
        }

    }
    // reverse par_pos
    for (i = 0; i < (*par_n)>> 1; ++i) {
        j = par_pos[i]; par_pos[i] = par_pos[*par_n-i-1]; par_pos[*par_n-i-1] = j;
    }
    par_pos[(*par_n)++] = est_start;
    par_pos[(*par_n)++] = est_start + est_period;

    ch_i1 = est_ch_i, s1 = est_start, e1 = est_start + est_period;
    while (ch_i1 < ch.len-1 && e1 <= last_start) {
        for (i = ch_i1+1; i < ch.len; ++i) {
            c = ch.cell[i]; d = dp[c.i][c.j];
            s2 = d.start; e2 = d.end;
            if (s2 == e1) {
                par_pos[(*par_n)++] = e2;
                ch_i1 = i, s1 = s2, e1 = e2;
                break;
            } else if (s2 > e1) { // non-overlapping adjacent tandem repeat hits
                // do global alignment
                iden_n = ksw2_global_with_cigar(bseq+s1-k+1, s2-s1+k, bseq+e1-k+1, e2-e1+k, &n_cigar, &cigar);
                #ifdef __DEBUG__
                printf("iden_n: %d (%d,%d), (%d,%d)\n", iden_n, e2-e1+k, s2-s1+k, s1, s2);
                #endif
                if (iden_n >= MIN_OF_TWO(s2-s1+k, e2-e1+k) * (1-mtp->max_div)) {
                    e = e2 - backtrack_left_end(n_cigar, cigar, s2-s1+k, e2-e1+k, s2-e1);
                    par_pos[(*par_n)++] = e; 
                    ch_i1 = i-1; s1 = e1, e1 = e;
                } else {
                    par_pos[(*par_n)++] = -1; 
                    par_pos[(*par_n)++] = s2; 
                    par_pos[(*par_n)++] = e2; 
                    ch_i1 = i; s1 = s2; e1 = e2;
                }
                if (cigar) 
                    free(cigar);
                break;
            }
        }
    }
    if (*par_n == 0) free(par_pos);
    return par_pos;
}

// TODO do extension first, so that cons would start from the begining of the repeat
int *get_partition_pos_with_narrow_global_alignment(uint8_t *bseq, int seq_len, dp_t **dp, chain_t ch, mini_tandem_para *mtp, int *par_n) {
    int est_ch_i = ch.est_ch_i, est_start = ch.est_start, est_period = ch.est_period;
    int first_start = dp[ch.cell[0].i][ch.cell[0].j].start, first_end = dp[ch.cell[0].i][ch.cell[0].j].end;
    int last_start = dp[ch.cell[ch.len-1].i][ch.cell[ch.len-1].j].start, last_end = dp[ch.cell[ch.len-1].i][ch.cell[ch.len-1].j].end;
    int *par_pos = (int*)_err_malloc(seq_len * sizeof(int)); *par_n = 0;
    int i, j, k=mtp->k, ch_i, s, e, s1, e1, s2, e2, iden_n;
    int n_cigar; uint32_t *cigar;
    cell_t c; dp_t d;
    // find upstream anchor (s1, e1) and downstream (s2, e2)
    // global alignment of [s1,s2] and [e1, e2]
    // find best partition position by backtracking
    // TODO overlapped hits => partition directly without alignment

    // left extension: find S
    // <=<=<=<=<=<=<=<=<=<=<=<=<=<=<=<=<=<=<=
    // --s1---S---s2----e1---s---e2-------e--
    // --|---------| vs |---------|----------
    ch_i = est_ch_i, s = est_start, e = est_start + est_period;
    while (s >= first_end && ch_i > 0) {
        // find (s1, e1) and (s2, e2)
        s2 = s, e2 = e; s1 = -1; e1 = -1;
        for (i = ch_i - 1; i >= 0; --i) {
            c = ch.cell[i]; d = dp[c.i][c.j];
            s1 = d.start; e1 = d.end;
            if (e1 == s) {
                par_pos[(*par_n)++] = s1;
                ch_i = i, s = s1, e = e1;
                break;
            } else if (e1 < s) {
                // if (s - e1 < k) // overlapped
                // do global alignment
                iden_n = ksw2_global_with_cigar(bseq+e1-k+1, e2-e1+k, bseq+s1-k+1, s2-s1+k, &n_cigar, &cigar);
                #ifdef __DEBUG__
                printf("1: (%d,%d), 2: (%d,%d)\n", s1, e1, s2, e2);
                printf("iden_n: %d (%d,%d), (%d,%d)\n", iden_n, e2-e1+k, s2-s1+k, s1, s2);
                #endif
                if (iden_n >= MIN_OF_TWO(s2-s1+k, e2-e1+k) * (1-mtp->max_div)) { // extend partition
                    e = s; s = s2 - backtrack_left_end(n_cigar, cigar, e2-e1+k, s2-s1+k, e2-s);
                    if (e == s) { // no backtrack
                        ch_i = 0; break;
                    }
                    par_pos[(*par_n)++] = s; 
                    ch_i = i+1;
                } else { // skip this anchor
                    par_pos[(*par_n)++] = -1; // insert a separation flag
                    par_pos[(*par_n)++] = e1; 
                    par_pos[(*par_n)++] = s1; 
                    ch_i = i; s = s1; e = e1;
                }
                if (cigar) free(cigar);
                break;
            } else {
                s2 = s1; e2 = e1;
            }
        }

    }
    // reverse par_pos
    for (i = 0; i < (*par_n)>> 1; ++i) {
        j = par_pos[i]; par_pos[i] = par_pos[*par_n-i-1]; par_pos[*par_n-i-1] = j;
    }
    par_pos[(*par_n)++] = est_start;
    par_pos[(*par_n)++] = est_start + est_period;

    // right extension: find E
    // >=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=
    // --s-------s1---e---s2----e1---E---e2--
    // ----------|---------| vs |---------|--
    ch_i = est_ch_i, s = est_start, e = est_start + est_period;
    while (ch_i < ch.len-1 && e <= last_start) {
        // find (s1,e1) and (s2,e2)
        s1 = s, e1 = e, s2 = e2 = -1;
        for (i = ch_i+1; i < ch.len; ++i) {
            c = ch.cell[i]; d = dp[c.i][c.j];
            s2 = d.start; e2 = d.end;
            if (s2 == e) {
                par_pos[(*par_n)++] = e2;
                ch_i = i; s = s2; e = e2;
                break;
            } else if (s2 > e) { // first (s2, e2)
                // if (s2 - e < k) // overlapped
                iden_n = ksw2_global_with_cigar(bseq+s1-k+1, s2-s1+k, bseq+e1-k+1, e2-e1+k, &n_cigar, &cigar);
                #ifdef __DEBUG__
                printf("1: (%d,%d), 2: (%d,%d)\n", s1, e1, s2, e2);
                printf("iden_n: %d (%d,%d), (%d,%d)\n", iden_n, e2-e1+k, s2-s1+k, s1, s2);
                #endif
                if (iden_n >= MIN_OF_TWO(s2-s1+k, e2-e1+k) * (1-mtp->max_div)) {
                    s = e; e = e2 - backtrack_left_end(n_cigar, cigar, s2-s1+k, e2-e1+k, s2-e);
                    if (e == s) { // no backtrack
                        ch_i = ch.len; break;
                    }
                    par_pos[(*par_n)++] = e; 
                    ch_i = i-1;
                } else {
                    par_pos[(*par_n)++] = -1; 
                    par_pos[(*par_n)++] = s2; 
                    par_pos[(*par_n)++] = e2; 
                    ch_i = i; s = s2; e = e2;
                }
                if (cigar) free(cigar);
                break;
            } else {
                s1 = s2; e1 = e2;
            }
        }
    }
    if (*par_n == 0) free(par_pos);
    return par_pos;
}

void seqs_msa(char *read_name, dp_t **dp, int ch_n, chain_t *chain, int seed_n, int *seed_ids, int seq_len, char *seq, uint8_t *bseq, tandem_seq_t *tseq, mini_tandem_para *mtp) {
    int ch_i; chain_t ch;
    for (ch_i = 0; ch_i < ch_n; ++ch_i) {
        ch = chain[ch_i];
        int par_n, *par_pos;
        par_pos = get_partition_pos_with_narrow_global_alignment(bseq, seq_len, dp, ch, mtp, &par_n);
        // par_pos = get_partition_pos_with_global_alignment(bseq, seq_len, dp, ch, mtp, &par_n);
        // par_pos = get_partition_pos(bseq, seq_len, dp, ch, mtp, &par_n);
        // par_pos = partition_seqs(seq, seq_len, dp, seed_ids, seed_n, ch, &par_n);
        if (par_n < mtp->min_copy+1) {
            free(par_pos);
            continue;
        }
        int start, end, cons_len=0;
        char *cons_seq = (char*)_err_malloc(seq_len * sizeof(char));
        uint8_t *cons_bseq = (uint8_t*)_err_malloc(seq_len * sizeof(uint8_t));
#ifdef __DEBUG__
        {
            int i, j, seq_i = 0; 
            for (i = 0; i < par_n-1; ++i) {
                if (par_pos[i] >= 0 && par_pos[i+1] >= 0) {
                    start = par_pos[i], end = par_pos[i+1];
                    printf(">seqs_%d:%d-%d\n", end-start, start, end);
                    for (j = start+1; j <= end; ++j) printf("%c", seq[j]);
                    printf("\n");
                }
            }
        }
#endif
        int i = 0, j, s, start_par_i, end_par_i; 
        while (i < par_n-mtp->min_copy) {
            if (par_pos[i] < 0) { 
                i++;
                continue;
            }
            for (j = i+1; j < par_n; ++j) {
                if (par_pos[j] < 0) break;
            }
            if (j - i > mtp->min_copy) { // do multiple sequence alignment and consensus calling for par_pos[i:j]
                cons_len = abpoa_msa(bseq, seq_len, par_pos+i, j-i, cons_bseq);
                for (s = 0; s < cons_len; ++s) cons_seq[s] = "ACGTN"[cons_bseq[s]]; cons_seq[cons_len] = '\0';

                int max_q, max_t, cons_start, cons_end; double copy_num = j-i-1;
                ksw2_left_ext(cons_bseq, cons_len, bseq, par_pos[i]+1, &max_q, &max_t); cons_start = par_pos[i] - max_t;
                // printf("max_q: %d, max_t: %d\n", max_q, max_t);
                copy_num += (max_q + 1.0) / cons_len;
                ksw2_right_ext(cons_bseq, cons_len, bseq+par_pos[j-1]+1, seq_len-par_pos[j-1]-1, &max_q, &max_t); cons_end = par_pos[j-1] + max_t + 1;
                // printf("max_q: %d, max_t: %d\n", max_q, max_t);
                copy_num += (max_q + 1.0) / cons_len;
                // rotate the cons based on splint_seq 
                int splint_rotated = 0;
                if (mtp->splint_seq != NULL && cons_len > mtp->splint_len) { // TODO use ksw
                    int ed, min_ed = mtp->splint_len, sp_len = mtp->splint_len, min_start, min_end, idx = -1;
                    // search splint within a full cons, forward and reverse
                    ed = edlib_align_HW(mtp->splint_seq, mtp->splint_len, cons_seq, cons_len, &start, &end);
                    if (ed >= 0 && sp_len * 0.9 <= end-start && sp_len * 1.1 >= end-start && ed < min_ed) {
                        min_ed = ed; idx = 0; // Forward single cons
                        min_start = start; min_end = end;
                    }
                    ed = edlib_align_HW(mtp->splint_rc_seq, mtp->splint_len, cons_seq, cons_len, &start, &end);
                    if (ed >= 0 && sp_len * 0.9 <= end-start && sp_len * 1.1 >= end-start && ed < min_ed) {
                        min_ed = ed; idx = 1; // Reverse-comp single cons
                        min_start = start; min_end = end;
                    }
                    // search splint within a concatenated 2 cons, forward and reverse
                    char *cons2 = (char*)_err_malloc(((cons_len << 1) + 1) * sizeof(char));
                    strcpy(cons2, cons_seq); strcpy(cons2+cons_len, cons_seq); cons2[cons_len<<1] = '\0'; // concatenated 2 copies
                    ed = edlib_align_HW(mtp->splint_seq, mtp->splint_len, cons2, cons_len<<1, &start, &end);
                    if (ed >= 0 && ed < min_ed && sp_len * 0.9 <= end-start && sp_len * 1.1 >= end-start && start < cons_len && end >= cons_len) {
                        min_ed = ed; idx = 2; // Forward 2 copies cons
                        min_start = start; min_end = end;
                    }
                    ed = edlib_align_HW(mtp->splint_rc_seq, mtp->splint_len, cons2, cons_len<<1, &start, &end);
                    if (ed >= 0 && ed < min_ed && sp_len * 0.9 <= end-start && sp_len * 1.1 >= end-start && start < cons_len && end >= cons_len) {
                        min_ed = ed; idx = 3; // Reverse-comp 2 copies cons
                        min_start = start; min_end = end;
                    }
                    // rotate cons based on ed result
                    switch(idx)
                    {
                        case 0: 
                        case 1:
                            if (cons_len-min_end-1 > 0 && cons_len-min_end-1+min_start >= mtp->min_p && cons_len-min_end-1+min_start <= mtp->max_p) {
                                memcpy(cons_seq, cons2 + min_end + 1, cons_len - min_end - 1);
                                memcpy(cons_seq+cons_len-min_end-1, cons2, min_start);
                                cons_seq[cons_len - min_end - 1 + min_start] = '\0';
                                cons_len = cons_len - min_end - 1 + min_start;
                                splint_rotated = 1;
                            } else cons_len = 0;
                            break;
                        case 2:
                        case 3:
                            min_end -= cons_len;
                            if (min_start-min_end-1 > mtp->min_p && min_start-min_end-1 <= mtp->max_p) {
                                memcpy(cons_seq, cons2 + min_end + 1, min_start - min_end - 1);
                                cons_seq[min_start-min_end-1] = '\0';
                                cons_len = min_start-min_end-1;
                                splint_rotated = 1;
                            } else cons_len = 0;
                            break;
                        default:;
                            // err_printf("%s: no splint sequence found in consensus sequence.\n", read_name);
                    }
                    free(cons2);
                } 
                // TODO find full-length sequence based on 5' and 3' adapter sequences
                int full_length = 0;
                if (mtp->five_seq != NULL && mtp->three_seq != NULL && cons_len > mtp->five_len + mtp->three_len) {
                }
                // write cons_seq into tseq
                write_tandem_cons_seq(tseq, cons_seq, cons_len, cons_start, cons_end, copy_num, mtp, splint_rotated);
            }
            i = j + 1;
        }
        free(cons_seq); free(cons_bseq); free(par_pos);
    }
}

// TODO polish everything!!!
// 0. build Hash index, generate coordinate of all the hits (x,y)
// 1. self-chaining based on (y-x) values by dynamic programming
// 2. keep top N chains, (calcuate density of each chain: tot_N_hits / tot_N_kmer)
//    2.1. 2 or more chain may co-exist because of template-switching 
//    2.2. post-analysis of multi-chains: 
//           switch-orientation: reverse complimentary
//           deletion: 
//           insertion:
// 3. call consensus with each chain
// 4. polish consensus result
int mini_tandem_core(kseq_t *read_seq, tandem_seq_t *tseq, mini_tandem_para *mtp) {
    char *seq = read_seq->seq.s; int seq_len = read_seq->seq.l; 

    if (seq_len < mtp->k) return 0;
    uint8_t *bseq = get_bseq(seq, seq_len);

    // generate hash value for each k-mer
    int hn; hash_t *h = (hash_t*)_err_malloc((seq_len - mtp->w) * sizeof(hash_t));
    if ((hn = build_kmer_hash(bseq, seq_len, mtp, h)) == 0) return 0;
#ifdef __DEBUG__
    printf("hash seed number: %d\n", hn);
#endif
    // collect hash hits
    hash_t *hit_h; int *seed_ids = (int*)_err_calloc(seq_len, sizeof(int)), seed_n;
    // int hit_n = collect_mem_hash_hit(h, hn, mtp->k, mtp->min_p, mtp->max_p, &hit_h, seed_ids, &seed_n); free(h);
    int hit_n = collect_hash_hit(h, hn, mtp->min_p, mtp->max_p, &hit_h, seed_ids, &seed_n); free(h);
    // dp chain
    dp_t **dp; int tot_n; chain_t *chain; int ch_m;
    int ch_n = dp_chain(read_seq->name.s, seq, bseq, seq_len, hit_h, hit_n, mtp, &dp, &tot_n, &chain, &ch_m);
    //TODO seperate partition and msa apart
    // partition into seqs and do msa
    seqs_msa(read_seq->name.s, dp, ch_n, chain, seed_n, seed_ids, seq_len, seq, bseq, tseq, mtp);

    int i;
    for (i = 0; i < ch_m; ++i) free(chain[i].cell); free(chain);
    free(seed_ids); free(hit_h); free(bseq);
    for (i = 0; i <= tot_n; ++i) free(dp[i]); free(dp); 
    return 0;
}
