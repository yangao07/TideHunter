#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mini_tandem.h"
#include "self_chain.h"
#include "seq.h"
#include "utils.h"
#include "ksort.h"

#define sort_key_hash(x) (x)
KRADIX_SORT_INIT(hash, hash_t, sort_key_hash, 8)

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
            // if (n > 1) hit_n += (n-1); // TODO use (n-1)
            n = 1;
        } else ++n;
    }
    if (n > 1) hit_n += (n * (n-1)) >> 1;
    // if (n > 1) hit_n += (n-1);

    // generate hits coordinate
    *hit_h = (hash_t*)_err_malloc(hit_n * sizeof(hash_t));
    int start_i = 0, j, k, hi=0;
    for (i = 1, n = 1; i < hn; ++i) {
        if ((h[i] >> 32) != (h[i-1] >> 32)) {
            if (n > 1) {
                for (j = 1; j < n; ++j) {
                    for (k = start_i; k < i-j; ++k) {
                        //(*hit_h)[hi++] = ((h[k+j] - h[k]) << 32) | (h[k+j] & 0xffffffff); // delta_p | end
                        (*hit_h)[hi++] = ((h[k+j] & 0xffffffff) << 32) | (h[k+j] - h[k]); // end | delta_p
                    }
                }
                //for (j = 1; j < n; ++j) 
                //    (*hit_h)[hi++] = ((h[start_i+j] - h[start_i+j-1]) << 32) | (h[start_i+j] & 0xffffffff);
            }
            start_i = i, n = 1;
        } else ++n;
    }
    if (n > 1) {
        for (j = 1; j < n; ++j) {
            for (k = start_i; k < i-j; ++k) {
                //(*hit_h)[hi++] = ((h[k+j] - h[k]) << 32) | (h[k+j] & 0xffffffff); // delta_p | end
                (*hit_h)[hi++] = ((h[k+j] & 0xffffffff) << 32) | (h[k+j] - h[k]); // end | delta_p
            }
        }
    }
    radix_sort_hash(*hit_h, *hit_h + hit_n);
    return hit_n;
}

int cluster_hash_hit(hash_t *hit_h, int hit_n, hit_bucket_t **hb, int *hb_m, double sigma, int bucket_T) {
    int i, hb_i=0;
    double delta_p;

    (*hb)[hb_i].min_p = (*hb)[hb_i].max_p = hit_h[0] >> 32; (*hb)[hb_i].start_i = (*hb)[hb_i].end_i = 0; (*hb)[hb_i].tot_n = 1;
    for (i = 1; i < hit_n; ++i) {
        delta_p = (double)(hit_h[i] >> 32) - (hit_h[i-1] >> 32);
        // printf("delta: %.2f, %.2f\n", delta_p, sigma );
        if (delta_p <= sigma * (hit_h[i] >> 32)) { // update current bucket
            (*hb)[hb_i].max_p = hit_h[i] >> 32;
            (*hb)[hb_i].end_i = i;
            (*hb)[hb_i].tot_n += 1;
        } else { // new delta, add new bucket
            if ((*hb)[hb_i].tot_n < bucket_T) {
                (*hb)[hb_i].min_p = (*hb)[hb_i].max_p = hit_h[i] >> 32;
                (*hb)[hb_i].start_i = (*hb)[hb_i].end_i = i;
                (*hb)[hb_i].tot_n = 1;
            } else {
                if ((hb_i+1) == *hb_m) {
                    *hb_m <<= 1;
                    *hb = (hit_bucket_t*)_err_realloc(*hb, *hb_m * sizeof(hit_bucket_t));
                }
                ++hb_i;
                (*hb)[hb_i].min_p = (*hb)[hb_i].max_p = hit_h[i] >> 32;
                (*hb)[hb_i].start_i = (*hb)[hb_i].end_i = i;
                (*hb)[hb_i].tot_n = 1;
            }
        }
    }
    int j;
    for (i = 0; i <= hb_i; ++i) {
        //printf("(%3d, %3d) => %4d\n", (*hb)[i].min_p, (*hb)[i].max_p, (*hb)[i].tot_n);
        //for (j = 0; j < (*hb)[i].tot_n; ++j) {
        //    printf("[%ld, %ld] , ", (hit_h[(*hb)[i].start_i+j] & 0xffffffff) - (hit_h[(*hb)[i].start_i+j] >> 32), hit_h[(*hb)[i].start_i+j] & 0xffffffff);
        //}
        //printf("\n");
    }
    return hb_i+1;
}

int overlap_rep_reg(int s1, int e1, int s2, int e2) {
    if (e1 < s2 || e2 < s1) return 0;
    else return 1;
}

int update_rep_reg(rep_reg_t **rr, int *rr_m, int *rr_n, int start, int end) {
    int i, hit = 0;
    for (i = 0; i < *rr_n; ++i) {
        if (overlap_rep_reg((*rr)[i].start, (*rr)[i].end, start, end)) {
            (*rr)[i].start = MIN_OF_TWO((*rr)[i].start, start);
            (*rr)[i].end = MAX_OF_TWO((*rr)[i].end, end);
            hit = 1;
        }
    }
    if (hit == 0) {
        if (*rr_n == *rr_m) {
            *rr_m <<= 1;
            *rr = (rep_reg_t*)_err_realloc(*rr, *rr_m * sizeof(rep_reg_t));
        }
        (*rr)[*rr_n].start = start, (*rr)[*rr_n].end = end;
        *rr_n += 1;
    }
    return 0;
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

// TODO calculate filtered hits projected region
int cal_proj_reg(hash_t *hit_h, int hit_n, hit_bucket_t *hb, int hb_n) {
    int i, j;
    int rr_m = 2, rr_n = 0;
    int *reg = (int*)_err_malloc(hit_n * 2 * sizeof(int));
    int reg_i = 0;
    //rep_reg_t *rr = (rep_reg_t*)_err_malloc(rr_m * sizeof(rep_reg_t));
    for (i = 0; i < hb_n; ++i) {
        for (j = 0; j < hb[i].tot_n; ++j) {
            int start = (hit_h[hb[i].start_i+j] & 0xffffffff) - (hit_h[hb[i].start_i+j] >> 32), end = hit_h[hb[i].start_i+j] & 0xffffffff;
            //reg[reg_i++] = start;
            reg[reg_i++] = end;
        }
    }

    //sort_rep_reg(&rr, &rr_n);
    qsort(reg, reg_i, sizeof(int), int_cmp);

    int w = 300;
    for (i = 0; i < reg_i; ++i) {
        int n = 0;
        for (j = i+1; j < reg_i; ++j) {
            if (reg[j] - reg[i] <= w) n++;
            else break;
        }
        printf("%d : %d\n", reg[i], n);
    }

    free(reg);
    return 0;
}

// TODO
int cal_reg_density(hash_t *hit_h, int hit_n) {
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

// backtrack from (x, y)
int backtrack_dp(self_dp_t **dp, int tot_n, int x, int y, int *start_x, int *start_y) {
    int cur_i = x, cur_j = y;
    int pre_i, pre_j;
    int chain_len = 0;
    while (1) {
        // chain_add_hit(cur_i, cur_j);
        dp[cur_i][cur_j].is_tracked = 1;
        ++chain_len;

        pre_i = dp[cur_i][cur_j].from_i;
        pre_j = dp[cur_i][cur_j].from_j;
        if (dp[pre_i][pre_j].is_tracked || pre_i == tot_n) break;
        cur_i = pre_i; cur_j = pre_j;
    }
    *start_x = cur_i, *start_y = cur_j;
    return chain_len;
}

void init_dp(hash_t *hit_h, self_dp_t **dp, int *hash_index, int *size, int total_n) {
    dp[total_n][0] = (self_dp_t){-1, -1, total_n, 0, 0.0, 0.0, 0};
    int i, j, hash_i;
    for (i = 0; i < total_n; ++i) {
        for (j = 0; j < size[i]; ++j) {
            hash_i = hash_index[i] + j;
            dp[i][j] = (self_dp_t){total_n, 0, i, j, (double)(hit_h[hash_i] & _32mask), 0.0, 0};
        }
    }
}

// TODO INS/DEL connection
double get_con_score(hash_t *hit_h, int pre_i, double pre_p, int cur_i, double *period, int k) {
    double cur_p = (double)(hit_h[cur_i] & _32mask);
    double f = cur_p / pre_p;
    double s = (double)(MIN_OF_TWO((hit_h[cur_i] >> 32) - (hit_h[pre_i] >> 32), k));

    if (0.9 <= f && f <= 1.1) { // similar period
        *period = cur_p;
    } else if ((int)(f+0.5) == 0) {
        *period = cur_p;
        s = 0;
    } else if (fabs((int)(f+0.5) - f) <= 0.1) { // x-fold period
        *period = pre_p;
        s /= (int)(f+0.5);
    } else { 
        *period = cur_p;
        s = 0;
    }
    return s;
}

// TODO allocate DP matrix uniformly
int self_dp_chain(hash_t *hit_h, int hit_n, int kmer_k) {
    int i, j, k;

    // calculate DP matrix size, allocate DP matrix
    int tot_n = 1, *array_size;
    int *hash_index;
    for (i = 1; i < hit_n; ++i) {
        if (hit_h[i] >> 32 != hit_h[i-1] >> 32) {
            tot_n += 1;
        }
    }
    array_size = (int*)_err_malloc(sizeof(int) * tot_n);
    hash_index = (int*)_err_malloc(sizeof(int) * tot_n);
    self_dp_t **dp = (self_dp_t**)_err_calloc((tot_n+1), sizeof(self_dp_t*));
    dp[tot_n] = (self_dp_t*)_err_calloc(1, sizeof(self_dp_t));
    j = 0, k = 1;
    int idx = 0;
    for (i = 1; i < hit_n; ++i) {
        if (hit_h[i] >> 32 != hit_h[i-1] >> 32) {
            dp[j] = (self_dp_t*)_err_calloc(k, sizeof(self_dp_t));
            hash_index[j] = idx-k+1;
            array_size[j++] = k;
            k = 1;
        } else ++k;
        ++idx;
    }
    dp[j] = (self_dp_t*)_err_calloc(k, sizeof(self_dp_t));
    hash_index[j] = idx-k;
    array_size[j] = k;

    // initialize DP matrix
    // set (tot_n,0) as all cells' precurser
    init_dp(hit_h, dp, hash_index, array_size, tot_n);

    // main DP process
    int cur_i, cur_j, pre_i, pre_j, pre_hash_i, cur_hash_i, max_pre_i, max_pre_j, iter_n;
    double score, max_score, con_score, con_period, max_period, pre_p;
    int max_h = 10;
    printf("%d\n", tot_n);
    for (cur_i = 1; cur_i < tot_n; ++cur_i) {
        for (cur_j = 0; cur_j < array_size[cur_i]; ++cur_j) {
            cur_hash_i = hash_index[cur_i] + cur_j;
            max_score = 0, iter_n = 0;
            for (pre_i = cur_i-1; pre_i >= 0; --pre_i) {
                if ((hit_h[hash_index[pre_i]] >> 32) < ((hit_h[cur_hash_i] >> 32) - (hit_h[cur_hash_i] & _32mask))) goto UPDATE;
                for (pre_j = 0; pre_j < array_size[pre_i]; ++pre_j) {
                    pre_hash_i = hash_index[pre_i] + pre_j;
                    pre_p = dp[pre_i][pre_j].period;
                    if ((con_score = get_con_score(hit_h, pre_hash_i, pre_p, cur_hash_i, &con_period, kmer_k)) == 0) continue;
                    score = dp[pre_i][pre_j].score + con_score;
                    if (score > max_score) {
                        max_score = score;
                        max_pre_i = pre_i, max_pre_j = pre_j;
                        max_period = con_period;
                        iter_n = 0;
                        if (fabs(con_period - (hit_h[cur_hash_i] & _32mask)) == 0) goto UPDATE;
                    } else if (max_score > 0) ++iter_n;
                    // only try h iterations
                    if (iter_n >= max_h) goto UPDATE;
                }
            }
UPDATE:
            if (max_score > 0) {
                dp[cur_i][cur_j] = (self_dp_t){max_pre_i, max_pre_j, cur_i, cur_j, max_period, max_score, 0};
            }
        }
    }

    // backtrack, obtain top N chains
    dp_score_t *score_rank = (dp_score_t*)_err_malloc(hit_n * sizeof(dp_score_t));
    int score_n = sort_dp_score(dp, array_size, tot_n, score_rank);
    int top_N = 10000, ch_n, start_i, start_j, chain_len;
    chain_t *ch = (chain_t*)_err_malloc(top_N * sizeof(chain_t));
    i = 0, ch_n = 0;
    while (i < score_n && ch_n < top_N) {
        chain_len = backtrack_dp(dp, tot_n, score_rank[i].i, score_rank[i].j, &start_i, &start_j);
        if (chain_len > 0) {
            ch[ch_n++] = (chain_t){start_i, start_j, score_rank[i].i, score_rank[i].j};
        }
        ++i;
    }
    for (i = 0; i < ch_n; ++i) {
        printf("%d: score: %lf (%ld,%lf,%ld) -> (%ld,%lf,%ld)\n", i+1, dp[ch[i].end_i][ch[i].end_j].score, hit_h[hash_index[ch[i].start_i] + ch[i].start_j] >> 32, dp[ch[i].start_i][ch[i].start_j].period, hit_h[hash_index[ch[i].start_i] + ch[i].start_j] & _32mask, hit_h[hash_index[ch[i].end_i] + ch[i].end_j] >> 32, dp[ch[i].end_i][ch[i].end_j].period, hit_h[hash_index[ch[i].end_i] + ch[i].end_j] & _32mask);
    }

    // post-process of N chains
    // 1. remove stand-alone hit
    for (i = 0; i < ch_n; ++i) {
    }
    // 2. merge into larger chain that contains INS/DEL

    for (i = 0; i <= tot_n; ++i) free(dp[i]); free(dp); free(array_size); free(hash_index);
    free(score_rank); free(ch);
    return 0;
}

// TODO post-process of projected region

// 0. build Hash index, generate coordinate of all the hits (x,y)
// 1. self-chaining based on (y-x) values by dynamic programming
// 2. keep top N chains, (calcuate density of each chain: tot_N_hits / tot_N_kmer)
// 3. call consensus with each chain
// 4. polish consensus result
int hash_partition(char *seq, int seq_len, mini_tandem_para *mtp) {
    int8_t *bseq = (int8_t*)_err_malloc(seq_len * sizeof(int8_t));
    get_bseq(seq, seq_len, bseq);
    hash_t *h = (hash_t*)_err_malloc(seq_len / mtp->w * sizeof(hash_t));
    int hn = direct_hash(bseq, seq_len, mtp->k, mtp->w, h);
    if (hn == 0) return 0;

    hash_t *hit_h;
    int hit_n = collect_hash_hit(h, hn, &hit_h); free(h);
    self_dp_chain(hit_h, hit_n, mtp->k);


    //int hb_m = 20; hit_bucket_t *hb;
    //hb = (hit_bucket_t*)_err_malloc(hb_m * sizeof(hit_bucket_t));
    //int hb_n = cluster_hash_hit(hit_h, hit_n, &hb, &hb_m, mtp->sigma, mtp->bucket_T);
    //cal_proj_reg(hit_h, hit_n, hb, hb_n-1);

    free(hit_h); free(bseq);
    return 0;
}
