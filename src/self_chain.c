#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mini_tandem.h"
#include "self_chain.h"
#include "edlib_align.h"
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
int direct_hash(uint8_t *bseq, int seq_len, int k, int s, int use_hpc, hash_t *h, int *cu_kmer_m) {
    int c, l, hi=0, pos, key;
    uint32_t mask = (1 << 2*k) - 1;

    for (key = l = pos = 0; pos < seq_len; ++pos) {
        c = bseq[pos];
        if (use_hpc) while (pos+1 < seq_len && bseq[pos+1] == c) ++pos;
        key = key << 2 | c;

        if (++l == k) { // get a kmer
            key &= mask; h[hi++] = ((hash_t)key << 32) | pos;
            cu_kmer_m[pos] = hi;
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
    int c, l, hi=0, pos; 
    uint32_t key, mask = (1 << 2*k) - 1;

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
int non_overlap_minimizer_hash(uint8_t *bseq, int seq_len, int k, int s, int w, int use_hpc, hash_t *h, int *cu_kmer_m) {
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
                cu_kmer_m[min_pos] = hi;
                for (j = 0; j < equal_n; ++j) { 
                    h[hi++] = ((hash_t)minimizer << 32) | equal_min_pos[j];
                    cu_kmer_m[equal_min_pos[j]] = hi;
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
int overlap_minimizer_hash(uint8_t *bseq, int seq_len, int k, int s, int w, int use_hpc, hash_t *h, int *cu_kmer_m) {
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
                cu_kmer_m[minimizer.y] = hi;
            } minimizer = tmp; min_buf_pos = buf_pos; equal_n = 0;
        } else if (tmp.x == minimizer.x) {
            equal_min_pos[equal_n++] = pos; min_buf_pos = buf_pos;
        }
        if (++l == s) { // reach a new w
            // collect minimizer of w-s for next w
            if ((min_buf_pos <= buf_pos && min_buf_pos > buf_pos - w + s) || min_buf_pos > buf_pos + s) { // current minimizer is in next w
                if (equal_n > 0) {
                    h[hi++] = ((hash_t)minimizer.x << 32) | minimizer.y;
                    cu_kmer_m[minimizer.y] = hi;
                    for (i = 0; i < equal_n-1; ++i) {
                        h[hi++] = ((hash_t)minimizer.x << 32) | equal_min_pos[i];
                        cu_kmer_m[minimizer.y] = hi;
                    }
                    minimizer.y = equal_min_pos[equal_n-1]; equal_n = 0;
                }
                last_min_buf_pos = min_buf_pos;
            } else { // out of next w
                h[hi++] = ((hash_t)minimizer.x << 32) | minimizer.y;
                cu_kmer_m[minimizer.y] = hi;
                for (i = 0; i < equal_n; ++i) {
                    h[hi++] = ((hash_t)minimizer.x << 32) | equal_min_pos[i];
                    cu_kmer_m[equal_min_pos[i]] = hi;
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
    cu_kmer_m[minimizer.y] = hi;
    for (i = 0; i < equal_n; ++i) { 
        h[hi++] = ((hash_t)minimizer.x << 32) | equal_min_pos[i];
        cu_kmer_m[equal_min_pos[i]] = hi;
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
int build_kmer_hash(uint8_t *bseq, int seq_len, mini_tandem_para *mtp, hash_t *h, int *cu_kmer_m) {
    if (mtp->w > 1) { // do min hash
        if (mtp->m == 1) {
            if (mtp->s >= mtp->w) return non_overlap_minimizer_hash(bseq, seq_len, mtp->k, mtp->s, mtp->w, mtp->hpc, h, cu_kmer_m);
            else return overlap_minimizer_hash(bseq, seq_len, mtp->k, mtp->s, mtp->w, mtp->hpc, h, cu_kmer_m);
        } // else return min_hash(bseq, seq_len, mtp->k, mtp->s, mtp->w, mtp->m, mtp->hpc, h);
    } else { // do direct hash
        return direct_hash(bseq, seq_len, mtp->k, mtp->s, mtp->hpc, h, cu_kmer_m);
    }
    return 0;
}

// TODO use multiple-hit seeds that connects multiple-repeats (>2)
// 0. h: hash_value | position
// 1. hash_hit1: period | end
// 2. hit_h2: end | period | K ; K : MEM hit's length
// 3. seed_ids[hit_end] = seed_id
// collect cumulative kmer and hits number for each pos (for estimate of n and e)
int collect_hash_hit(hash_t *h, int hn, hash_t **hash_hit, int *seed_ids, int *seed_n) {
    radix_sort_hash(h, h + hn); // sort h by hash values

    int i, n; int hit_n = 0; *seed_n = 0;

    // calculate total hits number
    for (i = 1, n = 1; i < hn; ++i) {
        if (h[i] >> 32 != h[i-1] >> 32) {
#ifdef __ALL_HIT__
            hit_n += (n * (n-1)) >> 1; // for each key, generate C(n,2) hits
#else
            hit_n += (n-1); // use (n-1): only collect the adjacent hit
#endif
            n = 1;
        } else ++n;
    }
#ifdef __ALL_HIT__
    hit_n += (n * (n-1)) >> 1;
#else
    hit_n += (n-1);
#endif

    // generate hash_hit1: period | end
    hash_t *hash_hit1 = (hash_t*)_err_malloc(hit_n * sizeof(hash_t));
    int start_i, j, hi;
    for (start_i = 0, hi = 0, i = 1, n = 1; i < hn; ++i) {
        if ((h[i] >> 32) != (h[i-1] >> 32)) {
            if (n > 1) {
                for (j = 1; j < n; ++j) {
#ifdef __ALL_HIT__
                    int k;
                    for (k = start_i; k < i-j; ++k) {
                        hash_hit1[hi] = _set_hash_hit1(h[k], h[k+j]);
                        ++hi;
                    }
#else
                    hash_hit1[hi] = _set_hash_hit1(h[start_i+j-1], h[start_i+j]);
                    ++hi;
                    // printf("end: %lld, seed_id: %d\n", h[start_i+j] & _32mask, *seed_n);
#endif
                    seed_ids[h[start_i+j] & _32mask] = *seed_n;
                }
                ++(*seed_n);
            }
            start_i = i, n = 1;
        } else ++n;
    }
    if (n > 1) {
        for (j = 1; j < n; ++j) {
#ifdef __ALL_HIT__
            int k;
            for (k = start_i; k < i-j; ++k) {
                hash_hit1[hi] = _set_hash_hit1(h[k], h[k+j]);
                ++hi;
            }
#else
            hash_hit1[hi] = _set_hash_hit1(h[start_i+j-1], h[start_i+j]);
            ++hi;
            // printf("end: %lld, seed_id: %d\n", h[start_i+j] & _32mask, *seed_n);
#endif
            seed_ids[h[start_i+j] & _32mask] = *seed_n;
        }
        ++(*seed_n);
    }
    radix_sort_hash(hash_hit1, hash_hit1 + hit_n); // sort hash_hit1 by period

    // generate hash_hit2: end:32 | period:16 | K:16 ; K: MEM hit's length
    int mem_hit_n;
    for (i = 1, mem_hit_n = 0, n = 0; i < hit_n; ++i) {
        if ((_get_hash_hit1_period(hash_hit1, i) != _get_hash_hit1_period(hash_hit1, i-1)) || (_get_hash_hit1_end(hash_hit1, i) != _get_hash_hit1_end(hash_hit1, i-1)+1)) {         

            ++mem_hit_n;
            n = 0;
        } else { // (_get_hash_hit1_end(hash_hit1, i) == _get_hash_hit1_end(hash_hit1, i-1)+1)
            if (n == _16mask) {
                n = 0;
                ++mem_hit_n;
            } else ++n;
        }
    }
    ++mem_hit_n;

    *hash_hit = (hash_t*)_err_malloc(mem_hit_n * sizeof(hash_t));
    for (hi = 0, i = 1, n = 0; i < hit_n; ++i) {
        // TODO mem dis > 1
        if ((_get_hash_hit1_period(hash_hit1, i) != _get_hash_hit1_period(hash_hit1, i-1)) || (_get_hash_hit1_end(hash_hit1, i) != _get_hash_hit1_end(hash_hit1, i-1)+1)) {         
            (*hash_hit)[hi++] = _set_hash_hit2(hash_hit1[i-1], _get_hash_hit1_period(hash_hit1, i-1), n);
            n = 0;
        } else { // (_get_hash_hit1_end(hash_hit1, i) == _get_hash_hit1_end(hash_hit1, i-1)+1)
            if (n == _16mask) {
                (*hash_hit)[hi++] = _set_hash_hit2(hash_hit1[i-1], _get_hash_hit1_period(hash_hit1, i-1), n);
                n = 0;
            } else ++n;
        }
    }
    (*hash_hit)[hi++] = _set_hash_hit2(hash_hit1[i-1], _get_hash_hit1_period(hash_hit1, i-1), n);
    radix_sort_hash((*hash_hit), (*hash_hit) + mem_hit_n); // sort hash_hit2 by end

    free(hash_hit1);
    return mem_hit_n;
}

int dp_score_cmp(const void *a, const void *b) {
    return (((dp_score_t*)b)->score - ((dp_score_t*)a)->score);
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

static inline void set_cur_hit_n(dp_t **dp, chain_t *ch, int *cu_hit_n, int seq_len) {
    int i, j, cu_n;

    for (j = 0; j < dp[ch->cell[0].i][ch->cell[0].j].end; ++j) cu_hit_n[j] = 0;
    for (i = cu_n = 1; i < ch->len; ++i) {
        for (j = dp[ch->cell[i-1].i][ch->cell[i-1].j].end; j < dp[ch->cell[i].i][ch->cell[i].j].end; ++j) cu_hit_n[j] = cu_n;
        ++cu_n;
    }
    for (; j < seq_len; ++j) cu_hit_n[j] = cu_n;
}

static inline void set_period(dp_t **dp, chain_t *ch, int *period) {
    int i, j, start, end, _period, mem_l;
    cell_t c; dp_t d;
    for (i = ch->len-1; i >= 0; --i) {
        c = ch->cell[i]; d = dp[c.i][c.j];
        start = d.start, end = d.end, mem_l = d.mem_l;
        _period = end - start;
        for (j = 0; j <= mem_l; ++j) {
            period[end-j] = _period;
            period[start-j] = _period;
            // printf("%d: %d, %d\n", end-j, start, _period);
        }
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

int is_overlap_chain(dp_t **dp, chain_t *chain, int ch_n, int ch_i) {
    if (ch_n <= 0) return 0;
    int i; cell_t c = chain[ch_i].cell[chain[ch_i].len-1]; int start = dp[c.i][c.j].start;
    for (i = ch_n-1; i >= 0; ++i) {
        if (dp[chain[i].cell[chain[i].len-1].i][chain[i].cell[chain[i].len-1].j].end <= start) break;
        if (is_overlap_chain1(dp, chain[i], chain[ch_i])) return 1;
    }
    return 0;
}


// backtrack from (x, y)
int backtrack_dp(dp_t **dp, int x, int y, chain_t *chain, int **cu_hit_n, int ch_n, int *period, int seq_len) {
    int cur_i = x, cur_j = y, pre_i, pre_j, chain_len = 0;
    chain_t *ch = chain+ch_n; int *cu = cu_hit_n[ch_n];
    while (1) {
        // chain_add_hit(cur_i, cur_j);
        dp[cur_i][cur_j].is_tracked = 1;
        ch->cell[chain_len++] = (cell_t){cur_i, cur_j};

        pre_i = dp[cur_i][cur_j].from_i;
        pre_j = dp[cur_i][cur_j].from_j;
        if (pre_i == -1 || dp[pre_i][pre_j].is_tracked) break;
        cur_i = pre_i; cur_j = pre_j;
    }
    ch->len = chain_len;
    if (chain_len > 1 && is_overlap_chain(dp, chain, ch_n, ch_n) == 0) {
        reverse_chain(ch);
        // set cu_hit_n
        set_cur_hit_n(dp, ch, cu, seq_len);
        set_period(dp, ch, period);
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
            end = _get_hash_hit2_end(hit_h, hash_i);
            period = _get_hash_hit2_period(hit_h, hash_i);
            mem_l = _get_hash_hit2_meml(hit_h, hash_i);
            start = end - period;
            dp[i][j] = (dp_t){-1, -1, start, end, mem_l, k+mem_l+0, 0};
            // printf("start: %d, end: %d\n", start, end);
        }
    }
}

// cur_score = pre_score + match_bases - gap_cost
// gap_cost = func(delta_period) from minimap2
static inline int get_con_score(dp_t *cur_dp, dp_t *pre_dp, int k, int mem_l, int *con_score) {
    int cur_start = cur_dp->start, pre_start = pre_dp->start;
    if (cur_start <= pre_start) return 0;  // cross-linked hits

    int cur_end = cur_dp->end, pre_end = pre_dp->end;
    int matched_bases = MIN_OF_TWO(cur_end - pre_end, mem_l);
    int delta_period = abs((cur_start - pre_start) - (cur_end - pre_end));
    int log_d = delta_period ? ilog2_32(delta_period) : 0;
    int gap_cost = (int)(.01 * delta_period * k) + (log_d >> 1);

    *con_score = matched_bases - gap_cost;
    return 1;
}

// if max potential overlap between new and existing chains >= 1/2 of existing, discard new chain
// chains are sorted by end
int is_in_chain(dp_t **dp, chain_t *ch, int *chain_idx, int ch_n, int cell_i, int cell_j) {
    int i, _i, cell_start = dp[cell_i][0].start, cell_end = dp[cell_i][cell_j].end;
    cell_t c1, c2;
    for (_i = 0; _i < ch_n; ++_i) {
        i = chain_idx[_i];
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
    // int *ch_idx = (int*)_err_malloc(sizeof(int) * ch_n);
    // for (i = 0; i < ch_n; ++i) ch_idx[i] = i;
    for (_i = 0; _i < ch_n-1; ++_i) {
        i = chain_idx[_i];
        c1 = chain[i].cell[chain[i].len-1];
        ch_end1 = dp[c1.i][c1.j].end;
        for (_j = _i+1; _j < ch_n; ++_j) {
            j = chain_idx[_j];
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
#ifdef __DEBUG__
    printf("Post#%d: %d -> %d\n", *ch_n+1, start_i, end_i);
#endif
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
// TODO only use n that are selected in chain
// calculate expected n for max_div, compare observed n with estimated n
// global alignment of two sequences, determine based on alignment score
int split_chain(char *seq, uint8_t *bseq, int seq_len, dp_t **dp, chain_t *ch, chain_t **post_ch, int *post_ch_n, int *post_ch_m, mini_tandem_para *mtp, int *cu_kmer_m, int *cu_hit_n) {
    int i, start, end, n, m; //, k = mtp->k, w=mtp->w;
    double est_n; dp_t dp_cell, pre_cell;
    int min_i = -1,  min_n = -1, pre_i, iden_n, min_len, from_i, from_j; dp_t min_cell;
    int post_start_i = 0;
    for (i = 0; i < ch->len; ++i) {
        dp_cell = dp[ch->cell[i].i][ch->cell[i].j];
        start = dp_cell.start; end = dp_cell.end;

        n = cu_hit_n[end]-cu_hit_n[start]+1, m = cu_kmer_m[end]-cu_kmer_m[start]+1;
        // without w*w
        est_n = m / mtp->div_exp; // printf("Poisson-split (%d,%d): m: %d, n: %d, e: %.3f\n", start, end, m, n, 1/(1.0*k) * log((m+0.0)/(n)));
        // est_n  = m * pow((1.0 - mtp->max_div), k); printf("Binomial-split(%d,%d): m: %d, n: %d, e: %.3f\n", start, end, m, n, 1.0 - pow((n+0.0)/(m), 1/(1.0*k)));

        if (dp_cell.from_i == -1) continue;

        if (min_n != -1 && start > min_cell.end) { // step out of last sparse valley
            // find precursor whose end <= min_cell.start
            from_i = min_cell.from_i; from_j = min_cell.from_j; iden_n = -1, min_len = min_cell.end-min_cell.start, pre_i = min_i;
            while (from_i != -1 && pre_i > 0) {
                pre_cell = dp[from_i][from_j]; --pre_i;
                if (pre_cell.end == min_cell.start) { // do global alignment
#ifdef __DEBUG__
                    iden_n = ksw2_global(bseq+pre_cell.start+1, min_cell.start-pre_cell.start, bseq+min_cell.start+1, min_cell.end-min_cell.start);
                    printf("G-ksw: %d / %d,%d\n", iden_n, pre_cell.end-pre_cell.start, min_cell.end-min_cell.start);
#endif
                    iden_n = edlib_align_NW(seq+pre_cell.start+1, min_cell.start-pre_cell.start, seq+min_cell.start+1, min_cell.end-min_cell.start);
#ifdef __DEBUG__
                    printf("G-ed:  %d / %d,%d (%d)\n", iden_n, pre_cell.end-pre_cell.start, min_cell.end-min_cell.start, pre_cell.end);
#endif
                    min_len = MIN_OF_TWO(min_cell.end-min_cell.start, min_cell.start-pre_cell.start);
                    break;
                } else if (pre_cell.end < min_cell.start) { // do extension alignment
                    // use edlib prefix (SHW)
#ifdef __DEBUG__
                    iden_n = ksw2_left_ext(bseq+pre_cell.start+1, min_cell.start-pre_cell.start, bseq+min_cell.start+1, min_cell.end-min_cell.start);
                    printf("E-ksw: %d / %d,%d\n", iden_n, min_cell.start-pre_cell.start, min_cell.end-min_cell.start);
#endif
                    iden_n = edlib_align_NW(seq+pre_cell.start+1, min_cell.start-pre_cell.start, seq+min_cell.start+1, min_cell.end-min_cell.start);
#ifdef __DEBUG__
                    printf("G-ed:  %d / %d,%d (%d)\n", iden_n, pre_cell.end-pre_cell.start, min_cell.end-min_cell.start, pre_cell.end);
#endif
                    min_len = MIN_OF_TWO(min_cell.end-min_cell.start, min_cell.start-pre_cell.start);
                    break;
                }
                from_i = pre_cell.from_i; from_j = pre_cell.from_j;
            }
            //min_len = MIN_OF_TWO(min_cell.end-min_cell.start, pre_cell.end-pre_cell.start);
            if (iden_n != -1 && iden_n / (min_len+0.0) < (1-mtp->max_div)) {
                // split chain
#ifdef __DEBUG__
                printf("Remove (%d,%d) (%d,%d) %d\n", min_cell.start, min_cell.end, pre_cell.start, pre_cell.end, pre_i);
#endif
                // from post_start_i to pre_cell
                copy_chain(ch, seq_len, post_start_i, pre_i, post_ch, post_ch_n, post_ch_m);
                post_start_i = min_i;
            }
            min_n = -1;
        }
        if ((double)n < est_n && (min_n == -1 || (start <= min_cell.end && n < min_n))) {
            min_i = i; min_n = n; min_cell = dp_cell;
        }
        // printf("min_cell: %d,%d\n", min_n, min_cell.end);
    }
    if (min_n != -1) {
        // find precursor whose end <= min_cell.start
        from_i = min_cell.from_i; from_j = min_cell.from_j; iden_n = -1, min_len = min_cell.end-min_cell.start, pre_i = min_i;
        while (from_i != -1 && pre_i > 0) {
            pre_cell = dp[from_i][from_j]; --pre_i;
            if (pre_cell.end == min_cell.start) { // do global alignment
#ifdef __DEBUG__
                iden_n = ksw2_global(bseq+pre_cell.start+1, min_cell.start-pre_cell.start, bseq+min_cell.start+1, min_cell.end-min_cell.start);
                printf("G-ksw: %d / %d,%d (%d)\n", iden_n, pre_cell.end-pre_cell.start, min_cell.end-min_cell.start, pre_cell.end);
#endif
                iden_n = edlib_align_NW(seq+pre_cell.start+1, min_cell.start-pre_cell.start, seq+min_cell.start+1, min_cell.end-min_cell.start);
#ifdef __DEBUG__
                printf("G-ed:  %d / %d,%d (%d)\n", iden_n, pre_cell.end-pre_cell.start, min_cell.end-min_cell.start, pre_cell.end);
#endif
                min_len = MIN_OF_TWO(min_cell.end-min_cell.start, pre_cell.end-pre_cell.start);
                break;
            } else if (pre_cell.end < min_cell.start) { // do extension alignment
                // use edlib prefix (SHW)
#ifdef __DEBUG__
                iden_n = ksw2_left_ext(bseq+pre_cell.start+1, min_cell.start-pre_cell.start, bseq+min_cell.start+1, min_cell.end-min_cell.start);
                printf("E-ksw: %d / %d,%d (%d)\n", iden_n, pre_cell.end-pre_cell.start, min_cell.end-min_cell.start, pre_cell.end);
#endif
                iden_n = edlib_align_NW(seq+pre_cell.start+1, min_cell.start-pre_cell.start, seq+min_cell.start+1, min_cell.end-min_cell.start);
#ifdef __DEBUG__
                printf("G-ed:  %d / %d,%d (%d)\n", iden_n, pre_cell.end-pre_cell.start, min_cell.end-min_cell.start, pre_cell.end);
#endif
                min_len = MIN_OF_TWO(min_cell.end-min_cell.start, pre_cell.end-pre_cell.start);
                break;
            }
            from_i = pre_cell.from_i; from_j = pre_cell.from_j;
        }
        //min_len = MIN_OF_TWO(min_cell.end-min_cell.start, pre_cell.end-pre_cell.start);
        if (iden_n != -1 && iden_n / (min_len+0.0) < (1-mtp->max_div)) {
            // split chain
#ifdef __DEBUG__
            printf("Remove (%d,%d) (%d,%d) %d\n", min_cell.start, min_cell.end, pre_cell.start, pre_cell.end, pre_i);
#endif
            // from post_start_i to pre_cell
            copy_chain(ch, seq_len, post_start_i, pre_i, post_ch, post_ch_n, post_ch_m);
            post_start_i = min_i;
        }
    }
    // from post_start_i to ch.len-1
    copy_chain(ch, seq_len, post_start_i, ch->len-1, post_ch, post_ch_n, post_ch_m);
    return 0;
}

// hash_hit: hash table of mem hits
// TODO allocate DP matrix uniformly
int dp_chain(char *seq, uint8_t *bseq, int seq_len, hash_t *hash_hit, int mem_hit_n, mini_tandem_para *mtp, dp_t ***_dp, int *tot_N, int *cu_kmer_m, chain_t **post_chain, int *post_ch_m, int *period) {
    int i, j, k, idx, kmer_k = mtp->k;
    // calculate DP matrix size, allocate DP matrix
    int tot_n = 1, *array_size, *hash_index;
    for (i = 1; i < mem_hit_n; ++i) {
        if (_get_hash_hit2_end(hash_hit, i)  != _get_hash_hit2_end(hash_hit, i-1))
            tot_n += 1;
    }
    *tot_N = tot_n;
    array_size = (int*)_err_malloc(sizeof(int) * tot_n);
    hash_index = (int*)_err_malloc(sizeof(int) * tot_n);
    *_dp = (dp_t**)_err_calloc((tot_n+1), sizeof(dp_t*));
    dp_t **dp = *_dp;
    dp[tot_n] = (dp_t*)_err_calloc(1, sizeof(dp_t));

    for (i = 1, j = 0, k = 1, idx = 0; i < mem_hit_n; ++i) {
        if (_get_hash_hit2_end(hash_hit, i)  != _get_hash_hit2_end(hash_hit, i-1)) {
            dp[j] = (dp_t*)_err_calloc(k, sizeof(dp_t));
            hash_index[j] = idx-k+1; array_size[j++] = k;
            k = 1;
        } else ++k;
        ++idx;
    }
    dp[j] = (dp_t*)_err_calloc(k, sizeof(dp_t));
    hash_index[j] = idx-k; array_size[j] = k;

    // initialize DP matrix
    // set (tot_n,0) as all cells' precurser
    init_dp(hash_hit, dp, hash_index, array_size, tot_n, kmer_k);

    // main DP process
    int cur_i, cur_j, pre_i, pre_j, max_pre_i, max_pre_j, con_score, mem_l, score, max_score; 
    int iter_n, max_h = 50; // TODO number of iterations
    dp_t *cur_dp, *pre_dp;
    for (cur_i = 1; cur_i < tot_n; ++cur_i) {
        for (cur_j = 0; cur_j < array_size[cur_i]; ++cur_j) {
            cur_dp = dp[cur_i]+cur_j;
            max_score = cur_dp->score; mem_l = (int)cur_dp->score;
            iter_n = 0;
            for (pre_i = cur_i-1; pre_i >= 0; --pre_i) {
                if (dp[pre_i][0].end < cur_dp->start) goto UPDATE;
                int gt = 0;
                for (pre_j = 0; pre_j < array_size[pre_i]; ++pre_j) {
                    pre_dp = dp[pre_i]+pre_j;
                    if (get_con_score(cur_dp, pre_dp, kmer_k, mem_l, &con_score) == 0) continue;
                    score = dp[pre_i][pre_j].score + con_score;
                    if (score > max_score) {
                        max_score = score; max_pre_i = pre_i, max_pre_j = pre_j;
                        gt = 1;
                        // iter_n = 0;
                        // if (fabs(con_period - (dp[cur_i][cur_j].end - dp[cur_i][cur_j].start)) == 0) goto UPDATE;
                    } // else if (max_score > 0) ++iter_n;
                    // only try h iterations
                    // if (iter_n >= max_h) goto UPDATE;
                }
                if (gt) iter_n = 0;
                else if (++iter_n >= max_h) goto UPDATE;
            }
UPDATE:
            if (max_score > cur_dp->score) {
                cur_dp->score = max_score; cur_dp->from_i = max_pre_i; cur_dp->from_j = max_pre_j;
            }
        }
    }

    // TODO backtrack, obtain top N chains, use max-heap
    dp_score_t *score_rank = (dp_score_t*)_err_malloc(mem_hit_n * sizeof(dp_score_t));
    int score_n = sort_dp_score(dp, array_size, tot_n, score_rank);
    int top_N = 100, ch_n = 0, ch_m = top_N;
    int **cu_hit_n = (int**)_err_malloc(top_N * sizeof(int*));
    chain_t *chain = (chain_t*)_err_malloc(top_N * sizeof(chain_t)); int *chain_idx = (int*)_err_malloc(sizeof(int) * top_N); // chain_idx[rank]: index in chain
    int post_ch_n = 0; *post_ch_m = top_N; *post_chain = (chain_t*)_err_malloc(top_N * sizeof(chain_t));
    for (i = 0; i < top_N; ++i) {
        chain[i].cell = (cell_t*)_err_malloc(tot_n * sizeof(cell_t)); chain[i].len = 0;
        (*post_chain)[i].cell = (cell_t*)_err_malloc(tot_n * sizeof(cell_t)); (*post_chain)[i].len = 0;
        chain_idx[i] = i;
        cu_hit_n[i] = (int*)_err_calloc(seq_len, sizeof(int));
    }
    for (i = ch_n = 0; i < score_n && ch_n < top_N; ++i) {
        if (is_in_chain(dp, chain, chain_idx, ch_n, score_rank[i].i, score_rank[i].j)) continue;
        if (backtrack_dp(dp, score_rank[i].i, score_rank[i].j, chain, cu_hit_n, ch_n, period, seq_len)) ++ch_n;
        sort_chain(dp, chain, chain_idx, ch_n);
    }
#ifdef __DEBUG__
    printf("ch_n: %d\n", ch_n);
    chain_t *ch = chain;
    for (i = 0; i < ch_n; ++i) {
        if (ch[i].len > 0) {
            int start_i = ch[i].cell[0].i, start_j = ch[i].cell[0].j, end_i = ch[i].cell[ch[i].len-1].i, end_j = ch[i].cell[ch[i].len-1].j;
            int from_i = dp[start_i][start_j].from_i, from_j = dp[start_i][start_j].from_j;
            // printf("%d: score: %d(pre: %d), score_density: %lf, hit_density: %lf, len: %d (%d,%d,%d) -> (%d,%d,%d)\n", i+1, dp[end_i][end_j].score, dp[from_i][from_j].score, (dp[end_i][end_j].score-dp[from_i][from_j].score+0.0)/(dp[end_i][end_j].end-dp[start_i][start_j].start), (ch[i].len+0.0)/(dp[end_i][end_j].end-dp[start_i][start_j].start), ch[i].len, dp[start_i][start_j].start,  dp[start_i][start_j].end-dp[start_i][start_j].start, dp[start_i][start_j].end-dp[start_i][start_j].start,  dp[end_i][end_j].end, dp[end_i][end_j].end-dp[end_i][end_j].start, dp[end_i][end_j].end-dp[end_i][end_j].start);
            j = 0;
            printf("\tchain: %d(%d): start: %d, end: %d, p: %d, mem_l: %d, score: %d\n", j+1, ch[i].cell[j].i, dp[ch[i].cell[j].i][ch[i].cell[j].j].start, dp[ch[i].cell[j].i][ch[i].cell[j].j].end, dp[ch[i].cell[j].i][ch[i].cell[j].j].end-dp[ch[i].cell[j].i][ch[i].cell[j].j].start, dp[ch[i].cell[j].i][ch[i].cell[j].j].mem_l, dp[ch[i].cell[j].i][ch[i].cell[j].j].score);
            for (j = 1; j < ch[i].len; ++j) {
                printf("\tchain: %d(%d): start: %d, end: %d, p: %d, mem_l: %d, score: %d, delta: %d\n", j+1, ch[i].cell[j].i, dp[ch[i].cell[j].i][ch[i].cell[j].j].start, dp[ch[i].cell[j].i][ch[i].cell[j].j].end, dp[ch[i].cell[j].i][ch[i].cell[j].j].end-dp[ch[i].cell[j].i][ch[i].cell[j].j].start, dp[ch[i].cell[j].i][ch[i].cell[j].j].mem_l, dp[ch[i].cell[j].i][ch[i].cell[j].j].score, dp[ch[i].cell[j].i][ch[i].cell[j].j].start- dp[ch[i].cell[j-1].i][ch[i].cell[j-1].j].start);
            }
        }
    }
    /* for (i = 0; i < seq_len; ++i) {
        printf("kmer_hit(%d): %d, %d\n", i, cu_kmer_m[i], cu_hit_n[0][i]);
    }*/
#endif
    // post-process of chains
    post_ch_n = 0;
    for (i = 0; i < ch_n; ++i) {
        split_chain(seq, bseq, seq_len, dp, chain+i, post_chain, &post_ch_n, post_ch_m, mtp, cu_kmer_m, cu_hit_n[i]); // split chain in sparse region
    }
    for (i = 0; i < ch_m; ++i) free(chain[i].cell); free(chain); free(chain_idx);
    for (i = 0; i < top_N; ++i) free(cu_hit_n[i]); free(cu_hit_n);
    free(array_size); free(hash_index); free(score_rank);
    return post_ch_n;
}

int write_tandem_cons_seq(tandem_seq_t *tseq, char *cons_seq, int start, int end) {
    if (tseq->cons_seq->seq.l + strlen(cons_seq) >= tseq->cons_seq->seq.m) {
        tseq->cons_seq->seq.m = tseq->cons_seq->seq.l + strlen(cons_seq) + 1;
        tseq->cons_seq->seq.s = (char*)_err_realloc(tseq->cons_seq->seq.s, tseq->cons_seq->seq.m * sizeof(char));
    }
    strcpy(tseq->cons_seq->seq.s + tseq->cons_seq->seq.l, cons_seq); tseq->cons_seq->seq.l += strlen(cons_seq); 

    if (tseq->cons_n == tseq->cons_m) {
        tseq->cons_m <<= 1;
        tseq->cons_start = (int*)_err_realloc(tseq->cons_start, tseq->cons_m * sizeof(int));
        tseq->cons_end = (int*)_err_realloc(tseq->cons_end, tseq->cons_m * sizeof(int));
        tseq->cons_len = (int*)_err_realloc(tseq->cons_len, tseq->cons_m * sizeof(int));
        tseq->cons_score = (int*)_err_realloc(tseq->cons_score, tseq->cons_m * sizeof(int));
    }
    tseq->cons_start[tseq->cons_n] = start; tseq->cons_end[tseq->cons_n] = end;
    tseq->cons_len[tseq->cons_n] = strlen(cons_seq); tseq->cons_score[tseq->cons_n] = 0; // TODO cons_score
    ++tseq->cons_n;
    return 0;
}
// pos: 0-base
// TODO HPC kmer length or discard MEM  ??
int *partition_seqs_core(char *seq, int seq_len, int *period, int8_t *hit_array, int chain_start, int chain_end, int *par_n) {
    int i, j;
    int *pos_array = (int*)_err_malloc(sizeof(int) * seq_len), *par_pos = (int*)_err_malloc(sizeof(int) * seq_len);
    int hit_n = 0;
    for (i = 0; i < seq_len; ++i) {
        if (hit_array[i]) {
            pos_array[hit_n++] = i;
        }
    }
    // partition seq into period seperated seqs
    int par_i = 0, l = 100, tot_len; // TODO length of l-mer XXX
    char *query_seq, *target_seq; int copy_num, ed, start, end, target_start;
    // extend par_pos on the left
    // TODO extend non-full copy ??? like TRF does
    copy_num = (int)((double)(pos_array[0] - chain_start) / period[pos_array[0]] + 0.5);
    if (copy_num >= 1) {
        for (j = copy_num; j >= 1; --j) {
            query_seq = seq + pos_array[0];
            target_start = pos_array[0] - period[pos_array[0]] * j - (l<<1);
            if (target_start < 0) continue;
            target_seq = seq + target_start;
            ed = edlib_align_HW(query_seq, l, target_seq, l<<2, &start, &end);
            if (ed >= 0) par_pos[par_i++] = target_start + start;
        }
    }
    par_pos[par_i++] = pos_array[0];
    for (i = 0; i < hit_n-1; ++i) {
        int ave_p = (period[pos_array[i+1]] + period[pos_array[i]]) / 2;
        // printf("%d: %d, %d, %d\n", pos_array[i], period[pos_array[i]], period[pos_array[i+1]], ave_p);
        copy_num = (int)((double)(pos_array[i+1] - pos_array[i]) / ave_p + 0.5);
        if (copy_num > 1) { // multiple copies: semi-global alignment of prefix l-mer using edlib
            tot_len = pos_array[i+1] - pos_array[i];
            query_seq = seq + pos_array[i];
            for (j = 1; j < copy_num; ++j) {
                target_start = pos_array[i] + ave_p * j - (l << 1);
                // target_seq = seq + pos_array[i] + tot_len / copy_num * j - 2 * l;
                target_seq = seq + target_start;
                ed = edlib_align_HW(query_seq, l, target_seq, l<<2, &start, &end);
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
    copy_num = (int)((double)(chain_end - pos_array[hit_n-1]) / period[pos_array[hit_n-1]] + 0.5);
    if (copy_num >= 1) {
        for (j = 1; j <= copy_num; ++j) {
            query_seq = seq + pos_array[hit_n-1];
            target_start = pos_array[hit_n-1] + period[pos_array[hit_n-1]] * j - (l<<1);
            if (target_start + (l<<2) >= seq_len) continue;
            target_seq = seq + target_start;
            ed = edlib_align_HW(query_seq, l, target_seq, l<<2, &start, &end);
            if (ed >= 0) par_pos[par_i++] = target_start + start;
        }
    }
    free(pos_array);
    *par_n = par_i; return par_pos;
}

// multi-hits in one period:
// very close distance of the same kmers causes a negetive score: two identical kmers in one period
static inline int get_max_hit(int8_t **hit_array, int *seed_ids, int seed_n, dp_t **dp, chain_t ch) {
    int i, j, start, end, seed_id, mem_l; dp_t dp_cell;
    int *seed_hits = (int*)_err_calloc(seed_n, sizeof(int));
    int *seed_last_pos = (int*)_err_malloc(seed_n * sizeof(int)); memset(seed_last_pos, -1, seed_n * sizeof(int));

    /*for (i = 0; i < ch.len; ++i) {
        dp_cell = dp[ch.cell[i].i][ch.cell[i].j];
        start = dp_cell.start, end = dp_cell.end, mem_l = dp_cell.mem_l;
        for (j = 0; j <= mem_l; ++j) {
            seed_id = seed_ids[end-j];
            // printf("seed_id: %d => start: %d, end: %d\n", seed_id, start-j, end-j);
            seed_hits[seed_id] += (1-hit_array[seed_id][start-j] + 1-hit_array[seed_id][end-j]);

            hit_array[seed_id][start-j] = 1;
            hit_array[seed_id][end-j] = 1;
        }
    }*/

    for (i = ch.len-1; i >= 0; --i) {
        dp_cell = dp[ch.cell[i].i][ch.cell[i].j];
        start = dp_cell.start, end = dp_cell.end, mem_l = dp_cell.mem_l;
        for (j = 0; j <= mem_l; ++j) {
            seed_id = seed_ids[end-j];

            if (seed_last_pos[seed_id] == -1) { // first
                seed_hits[seed_id] += 2;
                hit_array[seed_id][start-j] = 1;
                hit_array[seed_id][end-j] = 1;
                //printf("first\n");
            } else {
                if (seed_last_pos[seed_id] == end-j) { // 
                    seed_hits[seed_id] += 1;
                    hit_array[seed_id][start-j] = 1;
                    // printf("same\n");
                } else if (seed_last_pos[seed_id] - (end-j) < (end - start) / 2) { // too close
                    // seed_hits[seed_id] stay the same;
                    hit_array[seed_id][start-j] = 1;
                    hit_array[seed_id][seed_last_pos[seed_id]] = 0;
                    // printf("close\n");
                } else { 
                    seed_hits[seed_id] += 2;
                    hit_array[seed_id][start-j] = 1;
                    hit_array[seed_id][end-j] = 1;
                    // printf("regular\n");
                }
            }
            seed_last_pos[seed_id] = start-j;
            // printf("seed_id: %d => start: %d, end: %d\n", seed_id, start-j, end-j);
            // seed_hits[seed_id] += (1-hit_array[seed_id][start-j] + 1-hit_array[seed_id][end-j]);
        }
    }
    
    int max_hit_n = 3, max_id = -1; // max_hit_n >= 3
    for (i = 0; i < seed_n; ++i) {
        // printf("%d: %d\n", i, seed_hits[i]);
        if (seed_hits[i] > max_hit_n) {
            max_hit_n = seed_hits[i];
            max_id = i;
       }
    }
#ifdef __DEBUG__
    printf("max_i: %d, max_hit_n: %d\n", max_id, max_hit_n);
#endif
    free(seed_hits); free(seed_last_pos);
    return max_id;
}

int *partition_seqs(char *seq, int seq_len, dp_t **dp, int *period, int *seed_ids, int seed_n, chain_t ch, int chain_start, int chain_end, int *par_n) {
    int i, array_size = seq_len; // (int)(period) * 2;
    int8_t **hit_array = (int8_t**)_err_malloc(sizeof(int8_t*) * array_size);
    for (i = 0; i < array_size; ++i) hit_array[i] = (int8_t*)_err_calloc(seq_len, sizeof(int8_t));
    // fill hit array
    int max_id;
    int *par_pos = NULL; *par_n = 0;
    if ((max_id = get_max_hit(hit_array, seed_ids, seed_n, dp, ch)) >= 0)
        par_pos = partition_seqs_core(seq, seq_len, period, hit_array[max_id], chain_start, chain_end, par_n);
    for (i = 0; i < array_size; ++i) free(hit_array[i]); free(hit_array);
    return par_pos;
}

void seqs_msa(dp_t **dp, int ch_n, chain_t *chain, int *period, int seed_n, int *seed_ids, int seq_len, char *seq, uint8_t *bseq, tandem_seq_t *tseq, mini_tandem_para *mtp) {
    int i, j, ch_i; chain_t ch; int chain_start, chain_end;
    for (ch_i = 0; ch_i < ch_n; ++ch_i) {
        ch = chain[ch_i];
        chain_start = dp[ch.cell[0].i][ch.cell[0].j].start, chain_end = dp[ch.cell[ch.len-1].i][ch.cell[ch.len-1].j].end;
        int par_n, *par_pos;
        par_pos = partition_seqs(seq, seq_len, dp, period, seed_ids, seed_n, ch, chain_start, chain_end, &par_n);
        if (par_n == 0) continue;
        char **seqs = (char**)_err_malloc((par_n - 1) * sizeof(char*)), *cons_seq = (char*)_err_calloc(seq_len, sizeof(char));
        int seq_i = 0, start, end;
        for (i = 0; i < par_n-1; ++i) {
            if (par_pos[i] > 0 && par_pos[i+1] > 0) {
                start = par_pos[i], end = par_pos[i+1];
                seqs[seq_i] = (char*)_err_calloc((end - start + 1), sizeof(char));
                //strncpy(seqs[seq_i++], seq + start, end - start);
                for (j = start; j < end; ++j) {
                    seqs[seq_i][j-start] = "ACGTN"[bseq[j]]; // make sure it's upper case
                } seqs[seq_i++][end-start] = '\0';
#ifdef __DEBUG__
                printf("seqs(%d:%d,%d): %s\n", end-start, start, end, seqs[seq_i-1]);
#endif
            }
        }
        // msa of seqs
        spoa_msa(seqs, seq_i, cons_seq);
        int cons_len = strlen(cons_seq);
#ifdef __DEBUG__
        printf("cons(%d): %s\n", cons_len, cons_seq);
#endif

        // rotate the cons based on splint_seq
        if (mtp->splint_seq != NULL) {
            int ed, min_ed = mtp->splint_len, min_start, min_end, idx = -1;
            // search splint within a full cons, forward and reverse
            ed = edlib_align_HW(mtp->splint_seq, mtp->splint_len, cons_seq, cons_len, &start, &end);
            if (ed >= 0 && ed < min_ed) {
                min_ed = ed; idx = 0; // Forward single cons
                min_start = start; min_end = end;
            }
            ed = edlib_align_HW(mtp->splint_rc_seq, mtp->splint_len, cons_seq, cons_len, &start, &end);
            if (ed >= 0 && ed < min_ed) {
                min_ed = ed; idx = 1; // Reverse-comp single cons
                min_start = start; min_end = end;
            }
            // search splint within a concatenated 2 cons, forward and reverse
            char *cons2 = (char*)_err_malloc(((cons_len << 1) + 1) * sizeof(char));
            strcpy(cons2, cons_seq); strcpy(cons2+cons_len, cons_seq); cons2[cons_len<<1] = '\0'; // concatenated 2 copies
            ed = edlib_align_HW(mtp->splint_seq, mtp->splint_len, cons2, cons_len<<1, &start, &end);
            if (ed >= 0 && ed < min_ed && start < cons_len && end >= cons_len) {
                min_ed = ed; idx = 2; // Forward 2 copies cons
                min_start = start; min_end = end;
            }
            ed = edlib_align_HW(mtp->splint_rc_seq, mtp->splint_len, cons2, cons_len<<1, &start, &end);
            if (ed >= 0 && ed < min_ed && start < cons_len && end >= cons_len) {
                min_ed = ed; idx = 3; // Reverse-comp 2 copies cons
                min_start = start; min_end = end;
            }
            // rotate cons based on ed result
            switch(idx)
            {
                case 0: 
                case 1:
                    memcpy(cons_seq, cons2 + min_end + 1, cons_len - min_end - 1);
                    memcpy(cons_seq+cons_len-min_end-1, cons2, min_start);
                    cons_seq[cons_len - min_end - 1 + min_start] = '\0';
                    break;
                case 2:
                case 3:
                    min_end -= cons_len;
                    memcpy(cons_seq, cons2 + min_end, min_start - min_end - 1);
                    cons_seq[min_start-min_end] = '\0';
                    break;
                default:
                    err_printf("No splint sequence found in consensus sequence.\n");
            }
            free(cons2);
        }

        // write cons_seq into tseq
        write_tandem_cons_seq(tseq, cons_seq, par_pos[0], par_pos[par_n-1]);
        for (i = 0; i < seq_i; ++i) free(seqs[i]); free(seqs); free(cons_seq); free(par_pos);
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
int hash_partition(char *seq, int seq_len, tandem_seq_t *tseq, mini_tandem_para *mtp) {
    if (seq_len < mtp->k) return 0;
    int i;
    uint8_t *bseq = (uint8_t*)_err_malloc(seq_len * sizeof(uint8_t));
    get_bseq(seq, seq_len, bseq);

    // generate hash value for each k-mer
    int hn; hash_t *h = (hash_t*)_err_malloc((seq_len - mtp->w) / mtp->s * mtp->m * sizeof(hash_t));
    int *cu_kmer_m = (int*)_err_calloc(seq_len, sizeof(int));
    if ((hn = build_kmer_hash(bseq, seq_len, mtp, h, cu_kmer_m)) == 0) return 0;
    int last_m = 0;
    for (i = 0; i < seq_len; ++i) {
        if (cu_kmer_m[i] != 0) last_m = cu_kmer_m[i];
        else cu_kmer_m[i] = last_m;
    }
#ifdef __DEBUG__
    printf("hash seed number: %d\n", hn);
#endif

    // collect hash hits
    hash_t *hit_h; int *seed_ids = (int*)_err_calloc(seq_len, sizeof(int)), seed_n;
    int mem_hit_n = collect_hash_hit(h, hn, &hit_h, seed_ids, &seed_n); free(h);

    // dp chain
    dp_t **dp; int tot_n; chain_t *chain; int ch_m;
    int *period = (int*)_err_calloc(seq_len, sizeof(int));
    int ch_n = dp_chain(seq, bseq, seq_len, hit_h, mem_hit_n, mtp, &dp, &tot_n, cu_kmer_m, &chain, &ch_m, period);

    // partition into seqs and do msa
    seqs_msa(dp, ch_n, chain, period, seed_n, seed_ids, seq_len, seq, bseq, tseq, mtp);

    free(period);
    for (i = 0; i < ch_m; ++i) free(chain[i].cell); free(chain);
    free(cu_kmer_m); free(seed_ids); free(hit_h); free(bseq);
    for (i = 0; i <= tot_n; ++i) free(dp[i]); free(dp); 
    return 0;
}
