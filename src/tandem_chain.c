#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "tandem_chain.h"
#include "utils.h"

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
int backtrack_dp(dp_t **dp, int x, int y, chain_t *chain, int ch_n) {
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
    int start, end, period;
    for (i = 0; i < total_n; ++i) {
        for (j = 0; j < size[i]; ++j) {
            hash_i = hash_index[i] + j;
            end = _get_hash_hit_end(hit_h, hash_i);
            period = _get_hash_hit_period(hit_h, hash_i);
            start = end - period;
            dp[i][j] = (dp_t){-1, -1, start, end, k+MIN_OF_TWO(k,period), 0};
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
    int matched_bases = MIN_OF_TWO(abs(cur_end - pre_end), k) + MIN_OF_TWO(abs(cur_start-pre_start), k);
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

// TODO max allowed fold
uint64_t get_period_penalty(int P, triple_t *period, int p_n) {
    int i; uint64_t pp = 0;
    for (i = 0; i < p_n; ++i) {
        pp += get_period_penalty1(P, period[i].x);
    }
    return pp;
}

int get_adj_dis_penalty(dp_t **dp, chain_t *ch, int i) {
    int adj_ave_dis;
    cell_t c1, c2, c;
    c = ch->cell[i];
    if (i == 0) {
        c2 = ch->cell[i+1];
        adj_ave_dis = dp[c2.i][c2.j].start - dp[c.i][c.j].start + dp[c2.i][c2.j].end - dp[c.i][c.j].end;
    } else if (i == ch->len-1) {
        c1 = ch->cell[i-1];
        adj_ave_dis = dp[c.i][c.j].start - dp[c1.i][c1.j].start + dp[c.i][c.j].end - dp[c1.i][c1.j].end;
    } else {
        c1 = ch->cell[i-1]; c2 = ch->cell[i+1];
        adj_ave_dis = (dp[c2.i][c.j].start - dp[c1.i][c1.j].start + dp[c2.i][c2.j].end - dp[c1.i][c1.j].end) / 2;
    }
    return ilog2_32(adj_ave_dis) / 2;
}

// calculate period for each hit
void get_medoid_period(dp_t **dp, chain_t *ch) {
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
    int last_p = -1; uint64_t pp, min_pp = UINT64_MAX;
    ch->est_period = 0;
    for (i = 0; i < ch->len; ++i) {
        if (t[i].x == last_p) continue;
        if (t[i].x > ch->max_period) ch->max_period = t[i].x;
        if (t[i].x < ch->min_period) ch->min_period = t[i].x;
        last_p = t[i].x;
        pp = get_period_penalty(t[i].x, t, ch->len) + get_adj_dis_penalty(dp, ch, t[i].z) * ch->len;
#ifdef __DEBUG__
        printf("p: %d, penalty: %lld\n", t[i].x, pp);
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
int tandem_chain(int seq_len, hash_t *hash_hit, int hash_hit_n, mini_tandem_para *mtp, dp_t ***_dp, int *tot_N, chain_t **post_chain, int *post_ch_m) {
    if (hash_hit_n < 2) return 0;
    int i, j, k, idx, kmer_k = mtp->k;
    // calculate DP matrix size, allocate DP matrix
    int tot_n = 1, *array_size, *hash_index;
    for (i = 1; i < hash_hit_n; ++i) {
        if (_get_hash_hit_end(hash_hit, i)  != _get_hash_hit_end(hash_hit, i-1))
            tot_n += 1;
    }
    *tot_N = tot_n;
    array_size = (int*)_err_malloc(sizeof(int) * tot_n);
    hash_index = (int*)_err_malloc(sizeof(int) * tot_n);
    *_dp = (dp_t**)_err_calloc((tot_n+1), sizeof(dp_t*));
    dp_t **dp = *_dp;
    dp[tot_n] = (dp_t*)_err_calloc(1, sizeof(dp_t));

    for (i = 1, j = 0, k = 1, idx = 0; i < hash_hit_n; ++i) {
        if (_get_hash_hit_end(hash_hit, i)  != _get_hash_hit_end(hash_hit, i-1)) {
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
        if (backtrack_dp(dp, score_rank[i].i, score_rank[i].j, chain, ch_n)) ++ch_n;
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
            printf("\tchain: %d(%d): start: %d, end: %d, p: %d, score: %d\n", j+1, ch[i].cell[j].i, dp[ch[i].cell[j].i][ch[i].cell[j].j].start, dp[ch[i].cell[j].i][ch[i].cell[j].j].end, dp[ch[i].cell[j].i][ch[i].cell[j].j].end-dp[ch[i].cell[j].i][ch[i].cell[j].j].start, dp[ch[i].cell[j].i][ch[i].cell[j].j].score);
            for (j = 1; j < ch[i].len; ++j) {
                printf("\tchain: %d(%d): start: %d, end: %d, p: %d, score: %d, delta: %d\n", j+1, ch[i].cell[j].i, dp[ch[i].cell[j].i][ch[i].cell[j].j].start, dp[ch[i].cell[j].i][ch[i].cell[j].j].end, dp[ch[i].cell[j].i][ch[i].cell[j].j].end-dp[ch[i].cell[j].i][ch[i].cell[j].j].start, dp[ch[i].cell[j].i][ch[i].cell[j].j].score, dp[ch[i].cell[j].i][ch[i].cell[j].j].start- dp[ch[i].cell[j-1].i][ch[i].cell[j-1].j].start);
            }
        }
    }
#endif
    // post-process of chains
    post_ch_n = 0;
    for (i = ch_n-1; i >= 0; --i) {
        copy_chain(chain+chain_idx[i], seq_len, 0, chain[chain_idx[i]].len-1, post_chain, &post_ch_n, post_ch_m);
        //copy_chain(chain+i, seq_len, 0, chain[i].len-1, post_chain, &post_ch_n, post_ch_m); // do not split chain; split chain based on alignment
    }
    for (i = 0; i < post_ch_n; ++i)
        get_medoid_period(dp, (*post_chain)+i);
    for (i = 0; i < ch_m; ++i) free(chain[i].cell); free(chain); free(chain_idx);
    free(array_size); free(hash_index); free(score_rank);
    return post_ch_n;
}
