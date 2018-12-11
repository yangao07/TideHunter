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

// TODO merge consecutive hits into one longer hit
// 0. h: hash_value | position
// 1. hash_hit1: period | end
// 2. hit_h2: end | period | K ; K : MEM hit's length
// 3. seed_ids[hit_end] = seed_id
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

int overlap_rep_reg(int s1, int e1, int s2, int e2) {
    if (e1 < s2 || e2 < s1) return 0;
    else return 1;
}

int reg_cmp(const void *a, const void *b) {
    return ((rep_reg_t*)a)->start - ((rep_reg_t*)b)->start;
}

int int_cmp(const void *a, const void *b) {
    return *(int*)a - *(int*)b;
}

int sort_rep_reg(rep_reg_t **rr, int *rr_n) {
    qsort(*rr, *rr_n, sizeof(rep_reg_t), reg_cmp);
    int i, cur_i;
    for (i=1, cur_i=0; i < *rr_n; ++i) {
        if (overlap_rep_reg((*rr)[cur_i].start, (*rr)[cur_i].end, (*rr)[i].start, (*rr)[i].end)) {
            (*rr)[cur_i].start = MIN_OF_TWO((*rr)[cur_i].start, (*rr)[i].start);
            (*rr)[cur_i].end = MAX_OF_TWO((*rr)[cur_i].end, (*rr)[i].end);
            (*rr)[i].start = -1;
            (*rr)[i].end = -1;
        } else {
            cur_i = i;
        }
    }
    return 0;
}

int dp_score_cmp(const void *a, const void *b) {
    return (((dp_score_t*)a)->score < ((dp_score_t*)b)->score ? 1 : -1);
}

// TODO use heap to only keep top N scores
int sort_dp_score(self_dp_t **dp, int *array_size, int tot_n, dp_score_t *score_rank) { 
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
        tmp_i = ch->chain[i].i, tmp_j = ch->chain[i].j;
        ch->chain[i].i = ch->chain[ch->len-1-i].i, ch->chain[i].j = ch->chain[ch->len-1-i].j; 
        ch->chain[ch->len-1-i].i = tmp_i, ch->chain[ch->len-1-i].j = tmp_j; 
    }
}

// backtrack from (x, y)
void backtrack_dp(self_dp_t **dp, int tot_n, int x, int y, chain_t *_ch, int *ch_n) {
    int cur_i = x, cur_j = y;
    int pre_i, pre_j;
    chain_t *ch = _ch + *ch_n; int chain_len = 0;
    while (1) {
        // chain_add_hit(cur_i, cur_j);
        dp[cur_i][cur_j].is_tracked = 1;
        ch->chain[chain_len++] = (cell_t){cur_i, cur_j};

        pre_i = dp[cur_i][cur_j].from_i;
        pre_j = dp[cur_i][cur_j].from_j;
        if (dp[pre_i][pre_j].is_tracked || pre_i == tot_n) break;
        cur_i = pre_i; cur_j = pre_j;
    }

    if (chain_len > 1) {
        ch->len = chain_len;
        reverse_chain(ch);
        ++(*ch_n);
    }
}

void init_dp(hash_t *hit_h, self_dp_t **dp, int *hash_index, int *size, int total_n, int k) {
    dp[total_n][0] = (self_dp_t){-1, -1, 0, 0, 0, 0.0, 0, 0};
    int i, j, hash_i;
    int start, end, period, mem_l;
    for (i = 0; i < total_n; ++i) {
        for (j = 0; j < size[i]; ++j) {
            hash_i = hash_index[i] + j;
            end = _get_hash_hit2_end(hit_h, hash_i);
            period = _get_hash_hit2_period(hit_h, hash_i);
            mem_l = _get_hash_hit2_meml(hit_h, hash_i);
            start = end - period;
            dp[i][j] = (self_dp_t){total_n, 0, start, end, mem_l, (double)(period), k+mem_l+0, 0};
        }
    }
}

// TODO diff_threshold: based on distance and error ratio
// TODO conn_threshold: based on period and error ratio
//
// TODO INS/DEL connection
// cur_score = pre_score + match_bases - gap_cost
// gap_cost = func(delta_period)
// static inline int get_con_score(int cur_start, int cur_end, double cur_period, int pre_start, int pre_end, double pre_period, double *period, int k, int *con_score) {
static inline int get_con_score(self_dp_t *cur_dp, self_dp_t *pre_dp, double *period, int k, int mem_l, int *con_score) {
    int cur_start = cur_dp->start, pre_start = pre_dp->start;
    if (cur_start <= pre_start) return 0;  // cross-linked hits

    int cur_end = cur_dp->end, pre_end = pre_dp->end;
    double cur_period = cur_dp->period, pre_period = pre_dp->period;
    // TODO disconnect distance
    // G is based on (period, error-rate, k, w)
    // int G = 100;
    // if (cur_start - pre_start >= G) return 0;

    int matched_bases = MIN_OF_TWO(cur_end - pre_end, mem_l);
    int delta_period = abs((cur_start - pre_start) - (cur_end - pre_end));
    int log_d = delta_period ? ilog2_32(delta_period) : 0;
    int gap_cost = (int)(.01 * delta_period * k) + (log_d >> 1);
    *period = cur_period;

    *con_score = matched_bases - gap_cost;
    return 1;

    double f = cur_period / pre_period;
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
    return s;
}

// a hit should be considered as random-hit and removed:
// if the chain becomes two disconnected chains when the hit is removed
int remove_alone_hit(self_dp_t **dp, chain_t ch) {
    if (ch.len < 3) return 0;
    int i;
    int pre_i = ch.chain[0].i, pre_j = ch.chain[0].j; int pre_end = dp[pre_i][pre_j].end;
    int cur_i = ch.chain[1].i, cur_j = ch.chain[1].j; int cur_end = dp[cur_i][cur_j].end;
    int next_i, next_j, next_end, next_start;
    // double cur_p;
    for (i = 1; i < ch.len-1; ++i) {
        next_i = ch.chain[i+1].i, next_j = ch.chain[i+1].j;
        next_end = dp[next_i][next_j].end;
        next_start = dp[next_i][next_j].start;

        // cur_p = dp[cur_i][cur_j].period;
        if (next_start > pre_end) {
        // if ((double)(next_end - pre_end) > cur_p) {
            printf("remove: %d, %d\n", cur_end, cur_end - dp[cur_i][cur_j].start);
            ch.chain[i].i = -1;
        }

        pre_end = cur_end; cur_end = next_end;
        cur_i = next_i; cur_j = next_j;
    }

    return 0;
}

// hash_hit: hash table of mem hits
// TODO allocate DP matrix uniformly
chain_t self_dp_chain(hash_t *hash_hit, int mem_hit_n, int kmer_k, self_dp_t ***_dp, int *tot_N) {
    int i, j, k, idx;
    // calculate DP matrix size, allocate DP matrix
    int tot_n = 1, *array_size;
    int *hash_index;
    for (i = 1; i < mem_hit_n; ++i) {
        if (_get_hash_hit2_end(hash_hit, i)  != _get_hash_hit2_end(hash_hit, i-1))
            tot_n += 1;
    }
    *tot_N = tot_n;
    array_size = (int*)_err_malloc(sizeof(int) * tot_n);
    hash_index = (int*)_err_malloc(sizeof(int) * tot_n);
    *_dp = (self_dp_t**)_err_calloc((tot_n+1), sizeof(self_dp_t*));
    self_dp_t **dp = *_dp;
    dp[tot_n] = (self_dp_t*)_err_calloc(1, sizeof(self_dp_t));

    for (i = 1, j = 0, k = 1, idx = 0; i < mem_hit_n; ++i) {
        if (_get_hash_hit2_end(hash_hit, i)  != _get_hash_hit2_end(hash_hit, i-1)) {
            dp[j] = (self_dp_t*)_err_calloc(k, sizeof(self_dp_t));
            hash_index[j] = idx-k+1; array_size[j++] = k;
            k = 1;
        } else ++k;
        ++idx;
    }
    dp[j] = (self_dp_t*)_err_calloc(k, sizeof(self_dp_t));
    hash_index[j] = idx-k; array_size[j] = k;

    // initialize DP matrix
    // set (tot_n,0) as all cells' precurser
    init_dp(hash_hit, dp, hash_index, array_size, tot_n, kmer_k);

    // main DP process
    int cur_i, cur_j, pre_i, pre_j, max_pre_i, max_pre_j, iter_n;
    int con_score, mem_l;
    int score, max_score; double con_period, max_period;
    int max_h = 100; // TODO
    self_dp_t *cur_dp, *pre_dp;
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

                    if (get_con_score(cur_dp, pre_dp, &con_period, kmer_k, mem_l, &con_score) == 0) continue;
                    score = dp[pre_i][pre_j].score + con_score;
                    if (score > max_score) {
                        max_score = score;
                        max_pre_i = pre_i, max_pre_j = pre_j;
                        max_period = con_period;
                        iter_n = 0;
                        // if (fabs(con_period - (dp[cur_i][cur_j].end - dp[cur_i][cur_j].start)) == 0) goto UPDATE;
                    } else if (max_score > 0) ++iter_n;
                    // only try h iterations
                    if (iter_n >= max_h) goto UPDATE;
                }
            }
UPDATE:
            if (max_score > cur_dp->score) {
                cur_dp->from_i = max_pre_i;
                cur_dp->from_j = max_pre_j;
                cur_dp->period = max_period;
                cur_dp->score = max_score;
            }
        }
    }

    // backtrack, obtain top N chains, use max-heap TODO
    dp_score_t *score_rank = (dp_score_t*)_err_malloc(mem_hit_n * sizeof(dp_score_t));
    int score_n = sort_dp_score(dp, array_size, tot_n, score_rank);
    int top_N = 10000, ch_n;
    chain_t *ch = (chain_t*)_err_malloc(top_N * sizeof(chain_t));
    for (i = 0; i < top_N; ++i) {
        ch[i].chain = (cell_t*)_err_malloc(tot_n * sizeof(cell_t));
        ch[i].len = 0;
    }
    i = 0, ch_n = 0;
    while (i < score_n && ch_n < top_N) {
        backtrack_dp(dp, tot_n, score_rank[i].i, score_rank[i].j, ch, &ch_n);
        ++i;
    }
#ifdef __DEBUG__
    for (i = 0; i < ch_n; ++i) {
        if (ch[i].len > 0) {
            int start_i = ch[i].chain[0].i, start_j = ch[i].chain[0].j, end_i = ch[i].chain[ch[i].len-1].i, end_j = ch[i].chain[ch[i].len-1].j;
            int from_i = dp[start_i][start_j].from_i, from_j = dp[start_i][start_j].from_j;
            printf("%d: score: %d(pre: %d), score_density: %lf, hit_density: %lf, len: %d (%d,%lf,%d) -> (%d,%lf,%d)\n", i+1, dp[end_i][end_j].score, dp[from_i][from_j].score, (dp[end_i][end_j].score-dp[from_i][from_j].score+0.0)/(dp[end_i][end_j].end-dp[start_i][start_j].start), (ch[i].len+0.0)/(dp[end_i][end_j].end-dp[start_i][start_j].start), ch[i].len, dp[start_i][start_j].start,  dp[start_i][start_j].period, dp[start_i][start_j].end-dp[start_i][start_j].start,  dp[end_i][end_j].end, dp[end_i][end_j].period, dp[end_i][end_j].end-dp[end_i][end_j].start);
            j = 0;
            printf("\tchain: %d: start: %d, end: %d, score: %d\n", j+1, dp[ch[i].chain[j].i][ch[i].chain[j].j].start, dp[ch[i].chain[j].i][ch[i].chain[j].j].end, dp[ch[i].chain[j].i][ch[i].chain[j].j].score);
            for (j = 1; j < ch[i].len; ++j) {
                printf("\tchain: %d: start: %d, end: %d, p: %d, score: %d, delta: %d\n", j+1, dp[ch[i].chain[j].i][ch[i].chain[j].j].start, dp[ch[i].chain[j].i][ch[i].chain[j].j].end, dp[ch[i].chain[j].i][ch[i].chain[j].j].end-dp[ch[i].chain[j].i][ch[i].chain[j].j].start, dp[ch[i].chain[j].i][ch[i].chain[j].j].score, dp[ch[i].chain[j].i][ch[i].chain[j].j].start- dp[ch[i].chain[j-1].i][ch[i].chain[j-1].j].start);
            }
        }
    }
#endif

    // post-process of N chains
    // 1. remove stand-alone hit
    ch_n = 1;
    for (i = 0; i < ch_n; ++i) {
        remove_alone_hit(dp, ch[i]);
    }
    // 2. merge into larger chain that contains INS/DEL

    chain_t ret_ch; int select_i = 0;
    // ret_ch.len = ch[select_i].len;
    ret_ch.chain = (cell_t*)_err_malloc(ch[select_i].len * sizeof(cell_t));
    ret_ch.len = 0;
    for (i = 0; i < ch[select_i].len; ++i) {
        if (ch[select_i].chain[i].i < 0) break; // TODO break into multiple chains
        ret_ch.chain[i].i = ch[select_i].chain[i].i;
        ret_ch.chain[i].j = ch[select_i].chain[i].j;
        ++(ret_ch.len);
    }

    // for (i = 0; i <= tot_n; ++i) free(dp[i]); free(dp); 
    free(array_size); free(hash_index);
    for (i = 0; i < top_N; ++i) free(ch[i].chain); free(ch); free(score_rank);
    return ret_ch;
}

// TODO multi-hits in one period
// TODO pos: 1-base or 0-base???
int partition_seqs_core(char *seq, int8_t *hit_array, double period, int *par_pos) {
    int i, j, seq_len = strlen(seq);
    int *pos_array = (int*)_err_malloc(sizeof(int) * seq_len);
    int hit_n = 0;
    for (i = 0; i < seq_len; ++i) {
        if (hit_array[i]) pos_array[hit_n++] = i;
    }
    // partition seq into period seperated seqs
    // if (pos_array[0] > p) { // 0th copy
    // }
    int par_i = 0, l = 20, tot_len; // TODO length of l-mer
    char *query_seq, *target_seq; int ed, start, end;
    par_pos[par_i++] = pos_array[0];
    for (i = 0; i < hit_n-1; ++i) {
        int copy_num = (int)((double)(pos_array[i+1] - pos_array[i]) / period + 0.5);
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

static inline int get_max_hit(int8_t **hit_array, int *seed_ids, int seed_n, self_dp_t **dp, chain_t ch) {
    int i, j, start, end, seed_id, mem_l; self_dp_t dp_cell;
    int *seed_hits = (int*)_err_calloc(seed_n, sizeof(int));
    for (i = 0; i < ch.len; ++i) {
        dp_cell = dp[ch.chain[i].i][ch.chain[i].j];
        start = dp_cell.start, end = dp_cell.end, mem_l = dp_cell.mem_l;
        for (j = 0; j <= mem_l; ++j) {
            // printf("seed_id: %d => start: %d, end: %d\n", seed_id, start, end);
            seed_id = seed_ids[end-j];
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

// TODO reduce time
static inline int get_max_id(int seed_n, int8_t **hit_array, int seq_len) {
    int i, j, hit_n, max_hit_n = 0, max_id;
    for (i = 0; i < seed_n; ++i) {
        hit_n = 0;
        for (j = 0; j < seq_len; ++j) hit_n += hit_array[i][j];
        if (hit_n > max_hit_n) {
            max_hit_n = hit_n;
            max_id = i;
        }
    }
    printf("max_i: %d, max_hit_n: %d\n", max_id, max_hit_n); 
    return max_id;
}

int partition_seqs(char *seq, self_dp_t **dp, int *seed_ids, int seed_n, chain_t ch, int *par_pos) {
    int i, seq_len;

    double period = dp[ch.chain[0].i][ch.chain[0].j].period;
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
    int8_t *bseq = (int8_t*)_err_malloc(seq_len * sizeof(int8_t));
    get_bseq(seq, seq_len, bseq);

    // generate hash value for each k-mer
    hash_t *h = (hash_t*)_err_malloc(seq_len / mtp->w * sizeof(hash_t));
    int hn = direct_hash(bseq, seq_len, mtp->k, mtp->w, h);
    if (hn == 0) return 0;

    hash_t *hit_h; int *seed_ids = (int*)_err_calloc(seq_len, sizeof(int)), seed_n;
    // collect hash hits
    int mem_hit_n = collect_hash_hit(h, hn, &hit_h, seed_ids, &seed_n); free(h);

    // dp chain
    self_dp_t **dp; int tot_n;
    chain_t ch = self_dp_chain(hit_h, mem_hit_n, mtp->k, &dp, &tot_n);

    // partition seqs
    int p = (int)(dp[ch.chain[0].i][ch.chain[0].j].period + 0.5);
    int *par_pos = (int*)_err_malloc(seq_len / p * 2 * sizeof(int));
    int par_n = partition_seqs(seq, dp, seed_ids, seed_n, ch, par_pos);
    char **seqs = (char**)_err_malloc((par_n - 1) * sizeof(char*)); 
    char *cons_seq = (char*)_err_malloc(sizeof(char) * p * 2);
    int i, seq_i = 0, start, end;
    for (i = 0; i < par_n-1; ++i) {
        if (par_pos[i] > 0 && par_pos[i+1] > 0) {
            start = par_pos[i], end = par_pos[i+1];
            seqs[seq_i] = (char*)_err_calloc((end - start + 1), sizeof(char));
            strncpy(seqs[seq_i++], seq + start, end - start);
            printf("seqs(%d:%d,%d): %s\n", end-start, start, end, seqs[seq_i-1]);
        }
    }
    // msa of seqs
    spoa_msa(seqs, seq_i, cons_seq);
    printf("cons: %s\n", cons_seq);
    
    free(seed_ids); free(hit_h); free(bseq); free(ch.chain); free(par_pos);
    for (i = 0; i <= tot_n; ++i) free(dp[i]); free(dp); 
    for (i = 0; i < seq_i; ++i) free(seqs[i]); free(seqs); free(cons_seq);
    return 0;
}
