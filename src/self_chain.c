#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mini_tandem.h"
#include "self_chain.h"
#include "edlib_align.h"
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
int min_hash(uint8_t *bseq, int seq_len, int k, int s, int w, int m, int use_hpc, hash_t *h) {
    return 0;
}

// TODO skip N(ambiguous base)
// h: kmer-key | rightmost-pos
int build_kmer_hash(uint8_t *bseq, int seq_len, mini_tandem_para *mtp, hash_t *h, int *cu_kmer_m) {
    if (mtp->w > 1) { // do min hash
        if (mtp->m == 1) {
            if (mtp->s >= mtp->w) return non_overlap_minimizer_hash(bseq, seq_len, mtp->k, mtp->s, mtp->w, mtp->hpc, h, cu_kmer_m);
            else return overlap_minimizer_hash(bseq, seq_len, mtp->k, mtp->s, mtp->w, mtp->hpc, h, cu_kmer_m);
        } else return min_hash(bseq, seq_len, mtp->k, mtp->s, mtp->w, mtp->m, mtp->hpc, h);
    } else { // do direct hash
        return direct_hash(bseq, seq_len, mtp->k, mtp->s, mtp->hpc, h, cu_kmer_m);
    }
}

// TODO use seed-hit that connects multiple-repeats (>2)
// 0. h: hash_value | position
// 1. hash_hit1: period | end
// 2. hit_h2: end | period | K ; K : MEM hit's length
// 3. seed_ids[hit_end] = seed_id
// TODO collect cumulative kmer and hits number for each pos
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
                    for (k = start_i; k < i-j; ++k) {
                        hash_hit1[hi] = _set_hash_hit1(h[k], h[k+j]);
                        ++hi;
                    }
#else
                    hash_hit1[hi] = _set_hash_hit1(h[start_i+j-1], h[start_i+j]);
                    ++hi;
                    seed_ids[h[start_i+j] & _32mask] = *seed_n;
                    // printf("end: %lld, seed_id: %d\n", h[start_i+j] & _32mask, *seed_n);
#endif
                }
                ++(*seed_n);
            }
            start_i = i, n = 1;
        } else ++n;
    }
    if (n > 1) {
        for (j = 1; j < n; ++j) {
#ifdef __ALL_HIT__
            for (k = start_i; k < i-j; ++k) {
                hash_hit1[hi] = _set_hash_hit1(h[k], h[k+j]);
                ++hi;
            }
#else
            hash_hit1[hi] = _set_hash_hit1(h[start_i+j-1], h[start_i+j]);
            ++hi;
            seed_ids[h[start_i+j] & _32mask] = *seed_n;
            // printf("end: %lld, seed_id: %d\n", h[start_i+j] & _32mask, *seed_n);
#endif
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

// backtrack from (x, y)
int backtrack_dp(dp_t **dp, int tot_n, int x, int y, chain_t *ch, int *cu_hit_n, int seq_len) {
    int cur_i = x, cur_j = y, pre_i, pre_j, chain_len = 0;
    while (1) {
        // chain_add_hit(cur_i, cur_j);
        dp[cur_i][cur_j].is_tracked = 1;
        ch->cell[chain_len++] = (cell_t){cur_i, cur_j};

        pre_i = dp[cur_i][cur_j].from_i;
        pre_j = dp[cur_i][cur_j].from_j;
        if (dp[pre_i][pre_j].is_tracked || pre_i == tot_n) break;
        cur_i = pre_i; cur_j = pre_j;
    }

    if (chain_len > 1) {
        ch->len = chain_len;
        reverse_chain(ch);
        // set cu_hit_n
        set_cur_hit_n(dp, ch, cu_hit_n, seq_len);
        return 1;
    } else return 0;
}

void init_dp(hash_t *hit_h, dp_t **dp, int *hash_index, int *size, int total_n, int k) {
    dp[total_n][0] = (dp_t){-1, -1, 0, 0, 0, 0, 0};
    int i, j, hash_i;
    int start, end, period, mem_l;
    for (i = 0; i < total_n; ++i) {
        for (j = 0; j < size[i]; ++j) {
            hash_i = hash_index[i] + j;
            end = _get_hash_hit2_end(hit_h, hash_i);
            period = _get_hash_hit2_period(hit_h, hash_i);
            mem_l = _get_hash_hit2_meml(hit_h, hash_i);
            start = end - period;
            dp[i][j] = (dp_t){total_n, 0, start, end, mem_l, k+mem_l+0, 0};
            // printf("start: %d, end: %d\n", start, end);
        }
    }
}

// TODO diff_threshold: based on distance and error ratio
// TODO conn_threshold: based on period and error ratio
//
// cur_score = pre_score + match_bases - gap_cost
// gap_cost = func(delta_period) from minimap2
static inline int get_con_score(dp_t *cur_dp, dp_t *pre_dp, int k, int mem_l, int *con_score) {
    int cur_start = cur_dp->start, pre_start = pre_dp->start;
    if (cur_start <= pre_start) return 0;  // cross-linked hits

    int cur_end = cur_dp->end, pre_end = pre_dp->end;
    // double cur_period = cur_dp->period, pre_period = pre_dp->period;
    // TODO disconnect distance
    // G is based on (period, error-rate, k, w)
    // int G = 100;
    // if (cur_start - pre_start >= G) return 0;

    int matched_bases = MIN_OF_TWO(cur_end - pre_end, mem_l);
    int delta_period = abs((cur_start - pre_start) - (cur_end - pre_end));
    int log_d = delta_period ? ilog2_32(delta_period) : 0;
    int gap_cost = (int)(.01 * delta_period * k) + (log_d >> 1);
    // *period = cur_period;

    *con_score = matched_bases - gap_cost;
    return 1;

    /* double f = cur_period / pre_period;
    double s = (double)(MIN_OF_TWO(cur_end - pre_end, k));
    int embed = (cur_start <= pre_start);

    if (MIN_SIM_PERIOD_RAT <= f && f <= MAX_SIM_PERIOD_RAT) { // similar period
        if (embed) return 0;
        *period = cur_period;
    } else if ((int)(f+0.5) == 0) {
        *period = cur_period;
        s = 0;
    } else if (fabs((int)(f+0.5) - f) <= MAX_SIM_PERIOD_RAT - 1.0) { // x-fold period
        return 0; // disconnect multi-period TODO how to incorporate single- and multi-period hits
        *period = pre_period;
        s /= f;
    } else { 
        *period = cur_period;
        s = 0;
    }
    return s; */
}


int is_in_chain(dp_t **dp, chain_t *ch, int ch_n, int cell_i, int cell_j) {
    int i, j, cell_start = dp[cell_i][0].start, cell_end = dp[cell_i][cell_j].end;
    for (i = 0; i < ch_n; ++i) {
        int chain_start = dp[ch[i].cell[0].i][ch[i].cell[0].j].start, chain_end = dp[ch[i].cell[ch[i].len-1].i][ch[i].cell[ch[i].len-1].j].end;
        for (j = ch[i].len-1; j > 0; --j) {
            if (chain_end < cell_start) break;
            else if (chain_start <= cell_end && cell_start <= chain_end) return 1;
        }
    }
    return 0;
}
// a hit should be considered as random-hit and removed:
// if the chain becomes two disconnected chains when the hit is removed
// TODO calculate expected n for 20%-divergence, compare observed n with estimated n
int split_chain(dp_t **dp, chain_t *ch, mini_tandem_para *mtp, int *cu_kmer_m, int *cu_hit_n) {
    int i, start, end, est_n, n, m, k = mtp->k, w=mtp->w;
    dp_t dp_cell;
    for (i = 0; i < ch->len; ++i) {
        dp_cell = dp[ch->cell[i].i][ch->cell[i].j];
        start = dp_cell.start; end = dp_cell.end;
        n = cu_hit_n[end]-cu_hit_n[start]+1, m = cu_kmer_m[end]-cu_kmer_m[start]+1;
        // without w*w
        printf("Poisson-split (%d,%d): m: %d, n: %d, e: %.3f\n", start, end, m, n, 1/(1.0*k) * log((m+0.0)/(n)));
        printf("Binomial-split(%d,%d): m: %d, n: %d, e: %.3f\n", start, end, m, n, 1.0 - pow((n+0.0)/(m), 1/(1.0*k)));
        est_n = (int)(m / mtp->div_exp);

        if (n < est_n) { // too sparse
            // split into two chains
            printf("remove: (%d,%d): m: %d, n: %d, e: %.3f\n", start, end, m, n, 1.0 - pow((n+0.0)/(m), 1/(1.0*k)));
        }
    }

    /*
    if (ch->len < 3) return 0;
    int i;
    int pre_i = ch->cell[0].i, pre_j = ch->cell[0].j; int pre_end = dp[pre_i][pre_j].end;
    int cur_i = ch->cell[1].i, cur_j = ch->cell[1].j; int cur_end = dp[cur_i][cur_j].end;
    int next_i, next_j, next_end, next_start;
    // double cur_p;
    for (i = 1; i < ch->len-1; ++i) {
        next_i = ch->cell[i+1].i, next_j = ch->cell[i+1].j;
        next_end = dp[next_i][next_j].end;
        next_start = dp[next_i][next_j].start;

        // cur_p = dp[cur_i][cur_j].period;
        if (next_start > pre_end) {
        // if ((double)(next_end - pre_end) > cur_p) {
            printf("remove: %d, %d\n", cur_end, cur_end - dp[cur_i][cur_j].start);
            ch->cell[i].i = -1;
        }

        pre_end = cur_end; cur_end = next_end;
        cur_i = next_i; cur_j = next_j;
    }*/

    return 0;
}

// hash_hit: hash table of mem hits
// TODO allocate DP matrix uniformly
int dp_chain(hash_t *hash_hit, int mem_hit_n, int seq_len, mini_tandem_para *mtp, dp_t ***_dp, int *tot_N, int *cu_kmer_m, chain_t **chain, int *ch_m) {
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
    int iter_n, max_h = 100; // TODO number of iterations
    dp_t *cur_dp, *pre_dp;
    for (cur_i = 1; cur_i < tot_n; ++cur_i) {
        for (cur_j = 0; cur_j < array_size[cur_i]; ++cur_j) {
            cur_dp = dp[cur_i]+cur_j;
            max_score = cur_dp->score; mem_l = (int)cur_dp->score;
            iter_n = 0;
            for (pre_i = cur_i-1; pre_i >= 0; --pre_i) {
                // TODO too sparse, set threshold based on error rate profile
                if (dp[pre_i][0].end < cur_dp->start) goto UPDATE;
                for (pre_j = 0; pre_j < array_size[pre_i]; ++pre_j) {
                    pre_dp = dp[pre_i]+pre_j;
                    if (get_con_score(cur_dp, pre_dp, kmer_k, mem_l, &con_score) == 0) continue;
                    score = dp[pre_i][pre_j].score + con_score;
                    if (score > max_score) {
                        max_score = score; max_pre_i = pre_i, max_pre_j = pre_j;
                        iter_n = 0;
                        // if (fabs(con_period - (dp[cur_i][cur_j].end - dp[cur_i][cur_j].start)) == 0) goto UPDATE;
                    } else if (max_score > 0) ++iter_n;
                    // only try h iterations
                    if (iter_n >= max_h) goto UPDATE;
                }
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
    int top_N = 1000, ch_n; *ch_m = top_N;
    int **cu_hit_n = (int**)_err_malloc(top_N * sizeof(int*));
    *chain = (chain_t*)_err_malloc(top_N * sizeof(chain_t));
    for (i = 0; i < top_N; ++i) {
        (*chain)[i].cell = (cell_t*)_err_malloc(tot_n * sizeof(cell_t)); (*chain)[i].len = 0;
        cu_hit_n[i] = (int*)_err_calloc(seq_len, sizeof(int));
    }
    for (i = ch_n = 0; i < score_n && ch_n < top_N; ++i) {
        if (is_in_chain(dp, (*chain), ch_n, score_rank[i].i, score_rank[i].j)) continue;
        if (backtrack_dp(dp, tot_n, score_rank[i].i, score_rank[i].j, (*chain)+ch_n, cu_hit_n[ch_n], seq_len)) ++ch_n;
    }
#ifdef __DEBUG__
    printf("ch_n: %d\n", ch_n);
    chain_t *ch = *chain;
    for (i = 0; i < ch_n; ++i) {
        if (ch[i].len > 0) {
            int start_i = ch[i].cell[0].i, start_j = ch[i].cell[0].j, end_i = ch[i].cell[ch[i].len-1].i, end_j = ch[i].cell[ch[i].len-1].j;
            int from_i = dp[start_i][start_j].from_i, from_j = dp[start_i][start_j].from_j;
            printf("%d: score: %d(pre: %d), score_density: %lf, hit_density: %lf, len: %d (%d,%d,%d) -> (%d,%d,%d)\n", i+1, dp[end_i][end_j].score, dp[from_i][from_j].score, (dp[end_i][end_j].score-dp[from_i][from_j].score+0.0)/(dp[end_i][end_j].end-dp[start_i][start_j].start), (ch[i].len+0.0)/(dp[end_i][end_j].end-dp[start_i][start_j].start), ch[i].len, dp[start_i][start_j].start,  dp[start_i][start_j].end-dp[start_i][start_j].start, dp[start_i][start_j].end-dp[start_i][start_j].start,  dp[end_i][end_j].end, dp[end_i][end_j].end-dp[end_i][end_j].start, dp[end_i][end_j].end-dp[end_i][end_j].start);
            j = 0;
            printf("\tchain: %d(%d): start: %d, end: %d, p: %d, mem_l: %d, score: %d\n", j+1, ch[i].cell[j].i, dp[ch[i].cell[j].i][ch[i].cell[j].j].start, dp[ch[i].cell[j].i][ch[i].cell[j].j].end, dp[ch[i].cell[j].i][ch[i].cell[j].j].end-dp[ch[i].cell[j].i][ch[i].cell[j].j].start, dp[ch[i].cell[j].i][ch[i].cell[j].j].mem_l, dp[ch[i].cell[j].i][ch[i].cell[j].j].score);
            for (j = 1; j < ch[i].len; ++j) {
                printf("\tchain: %d(%d): start: %d, end: %d, p: %d, mem_l: %d, score: %d, delta: %d\n", j+1, ch[i].cell[j].i, dp[ch[i].cell[j].i][ch[i].cell[j].j].start, dp[ch[i].cell[j].i][ch[i].cell[j].j].end, dp[ch[i].cell[j].i][ch[i].cell[j].j].end-dp[ch[i].cell[j].i][ch[i].cell[j].j].start, dp[ch[i].cell[j].i][ch[i].cell[j].j].mem_l, dp[ch[i].cell[j].i][ch[i].cell[j].j].score, dp[ch[i].cell[j].i][ch[i].cell[j].j].start- dp[ch[i].cell[j-1].i][ch[i].cell[j-1].j].start);
            }
        }
    }
#endif

    for (i = 0; i < seq_len; ++i) {
        printf("kmer_hit(%d): %d, %d\n", i, cu_kmer_m[i], cu_hit_n[0][i]);
    }
    // post-process of N chains
    // 1. remove stand-alone hit
    for (i = 0; i < ch_n; ++i) {
        split_chain(dp, (*chain)+i, mtp, cu_kmer_m, cu_hit_n[i]);
    }
    // 2. merge into larger chain that contains INS/DEL

    /*
    chain_t ret_ch; int select_i = 0;
    ret_ch.cell = (cell_t*)_err_malloc((*chain)[select_i].len * sizeof(cell_t));
    ret_ch.len = 0;
    for (i = 0; i < (*chain)[select_i].len; ++i) {
        if ((*chain)[select_i].cell[i].i < 0) break; // TODO break into multiple chains
        ret_ch.cell[i].i = (*chain)[select_i].cell[i].i;
        ret_ch.cell[i].j = (*chain)[select_i].cell[i].j;
        ++(ret_ch.len);
    }*/

    for (i = 0; i < top_N; ++i) free(cu_hit_n[i]); free(cu_hit_n);
    free(array_size); free(hash_index); free(score_rank);
    return ch_n;
}

// TODO HPC kmer length
// TODO multi-hits in one period
// TODO pos: 1-base or 0-base???
// TODO left- and right-search of max-hit query
int partition_seqs_core(char *seq, int8_t *hit_array, int period, int *par_pos) {
    int i, j, seq_len = strlen(seq);
    int *pos_array = (int*)_err_malloc(sizeof(int) * seq_len);
    int hit_n = 0;
    for (i = 0; i < seq_len; ++i) {
        if (hit_array[i]) {
            pos_array[hit_n++] = i;
        }
    }
    // partition seq into period seperated seqs
    // if (pos_array[0] > p) { // 0th copy
    // }
    int par_i = 0, l = 20, tot_len; // TODO length of l-mer
    char *query_seq, *target_seq; int ed, start, end;
    par_pos[par_i++] = pos_array[0];
    for (i = 0; i < hit_n-1; ++i) {
        int copy_num = (int)((double)(pos_array[i+1] - pos_array[i]) / (double)period + 0.5);
        //if (copy_num == 0) 
        //    err_fatal(__func__, "Unexpected copy number. (%d, %d, %lf)\n", pos_array[i], pos_array[i+1], period);

        if (copy_num > 1) { // multiple copies: semi-global alignment of prefix l-mer using edlib
            tot_len = pos_array[i+1] - pos_array[i];
            query_seq = seq + pos_array[i];
            for (j = 1; j < copy_num; ++j) {
                target_seq = seq + pos_array[i] + tot_len / copy_num * j - 2 * l;
                ed = edlib_align_HW(query_seq, l, target_seq, 4 * l, &start, &end);
                if (ed < 0) { // no alignment result
                    par_pos[par_i++] = -1; // skip this copy
                } else {
                    par_pos[par_i++] = pos_array[i] + tot_len / copy_num * j - 2 * l + start;
                }
            }
        }
        par_pos[par_i++] = pos_array[i+1];
    }

    free(pos_array);
    return par_i;
}

static inline int get_max_hit(int8_t **hit_array, int *seed_ids, int seed_n, dp_t **dp, chain_t ch) {
    int i, j, start, end, seed_id, mem_l; dp_t dp_cell;
    int *seed_hits = (int*)_err_calloc(seed_n, sizeof(int));
    for (i = 0; i < ch.len; ++i) {
        dp_cell = dp[ch.cell[i].i][ch.cell[i].j];
        start = dp_cell.start, end = dp_cell.end, mem_l = dp_cell.mem_l;
        for (j = 0; j <= mem_l; ++j) {
            seed_id = seed_ids[end-j];
            // printf("seed_id: %d => start: %d, end: %d\n", seed_id, start-j, end-j);
            seed_hits[seed_id] += (1-hit_array[seed_id][start-j]);
            seed_hits[seed_id] += (1-hit_array[seed_id][end-j]);

            hit_array[seed_id][start-j] = 1;
            hit_array[seed_id][end-j] = 1;
        }
    }
    int max_hit_n = 0, max_id;
    for (i = 0; i < seed_n; ++i) {
        if (seed_hits[i] > max_hit_n) {
            max_hit_n = seed_hits[i];
            max_id = i;
        }
    }
    printf("max_i: %d, max_hit_n: %d\n", max_id, max_hit_n);
    free(seed_hits);
    return max_id;
}

int partition_seqs(char *seq, dp_t **dp, int *seed_ids, int seed_n, chain_t ch, int *par_pos) {
    int i, seq_len;

    int period = dp[ch.cell[0].i][ch.cell[0].j].end - dp[ch.cell[0].i][ch.cell[0].j].start;
    seq_len = strlen(seq);
    int array_size = seq_len; // (int)(period) * 2;
    int8_t **hit_array = (int8_t**)_err_malloc(sizeof(int8_t*) * array_size);
    for (i = 0; i < array_size; ++i) hit_array[i] = (int8_t*)_err_calloc(seq_len, sizeof(int8_t));
    // fill hit array
    int max_id = get_max_hit(hit_array, seed_ids, seed_n, dp, ch);
    int par_n = partition_seqs_core(seq, hit_array[max_id], period, par_pos);

    for (i = 0; i < array_size; ++i) free(hit_array[i]); free(hit_array);
    return par_n;
}

void seqs_msa(dp_t **dp, int ch_n, chain_t *chain, int seed_n, int *seed_ids, int seq_len, char *seq, uint8_t *bseq) {
    int i, j, ch_i; chain_t ch;
    for (ch_i = 0; ch_i < ch_n; ++ch_i) {
        ch = chain[ch_i];
        int p = dp[ch.cell[0].i][ch.cell[0].j].end - dp[ch.cell[0].i][ch.cell[0].j].start; 
        int *par_pos = (int*)_err_malloc(seq_len / p * 2 * sizeof(int));
        int par_n = partition_seqs(seq, dp, seed_ids, seed_n, ch, par_pos);

        char **seqs = (char**)_err_malloc((par_n - 1) * sizeof(char*)), *cons_seq = (char*)_err_malloc(sizeof(char) * p * 2);
        int seq_i = 0, start, end;
        for (i = 0; i < par_n-1; ++i) {
            if (par_pos[i] > 0 && par_pos[i+1] > 0) {
                start = par_pos[i], end = par_pos[i+1];
                seqs[seq_i] = (char*)_err_calloc((end - start + 1), sizeof(char));
                //strncpy(seqs[seq_i++], seq + start, end - start);
                for (j = start; j < end; ++j) {
                    seqs[seq_i][j-start] = "ACGTN"[bseq[j]]; // make sure it's upper case
                } seqs[seq_i++][end-start] = '\0';
                printf("seqs(%d:%d,%d): %s\n", end-start, start, end, seqs[seq_i-1]);
            }
        }
        // msa of seqs
        spoa_msa(seqs, seq_i, cons_seq);
        printf("cons: %s\n", cons_seq);

        free(par_pos);
        for (i = 0; i < seq_i; ++i) free(seqs[i]); 
        free(seqs); free(cons_seq);
    }
}

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
int hash_partition(char *seq, int seq_len, mini_tandem_para *mtp) {
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
    printf("hash seed number: %d\n", hn);

    // collect hash hits
    hash_t *hit_h; int *seed_ids = (int*)_err_calloc(seq_len, sizeof(int)), seed_n;
    int mem_hit_n = collect_hash_hit(h, hn, &hit_h, seed_ids, &seed_n); free(h);

    // dp chain
    dp_t **dp; int tot_n; chain_t *chain; int ch_m;
    int ch_n = dp_chain(hit_h, mem_hit_n, seq_len, mtp, &dp, &tot_n, cu_kmer_m, &chain, &ch_m);

    // partition into seqs and do msa
    seqs_msa(dp, ch_n, chain, seed_n, seed_ids, seq_len, seq, bseq);

    for (i = 0; i < ch_m; ++i) free(chain[i].cell); free(chain);
    free(cu_kmer_m); free(seed_ids); free(hit_h); free(bseq);
    for (i = 0; i <= tot_n; ++i) free(dp[i]); free(dp); 
    return 0;
}
