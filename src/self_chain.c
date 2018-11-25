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

void reverse_chain(chain_t *ch) {
    int i, tmp_i, tmp_j;
    for (i = 0; i < ch->len >> 1; ++i) {
        tmp_i = ch->chain[i].i, tmp_j = ch->chain[i].j;
        ch->chain[i].i = ch->chain[ch->len-i].i, ch->chain[i].j = ch->chain[ch->len-i].j; 
        ch->chain[ch->len-i].i = tmp_i, ch->chain[ch->len-i].j = tmp_j; 
    }
}

// backtrack from (x, y)
int backtrack_dp(self_dp_t **dp, int tot_n, int x, int y, chain_t *ch) {
    int cur_i = x, cur_j = y;
    int pre_i, pre_j;
    int chain_len = 0;
    while (1) {
        // chain_add_hit(cur_i, cur_j);
        dp[cur_i][cur_j].is_tracked = 1;
        ch->chain[++chain_len] = (cell_t){cur_i, cur_j};

        pre_i = dp[cur_i][cur_j].from_i;
        pre_j = dp[cur_i][cur_j].from_j;
        if (dp[pre_i][pre_j].is_tracked || pre_i == tot_n) break;
        cur_i = pre_i; cur_j = pre_j;
    }
    ch->len = chain_len;
    // reverse chain
    reverse_chain(ch);
    return chain_len;
}

void init_dp(hash_t *hit_h, self_dp_t **dp, int *hash_index, int *size, int total_n) {
    dp[total_n][0] = (self_dp_t){-1, -1, 0, 0, 0.0, 0.0, 0};
    int i, j, hash_i;
    int start, end, period;
    for (i = 0; i < total_n; ++i) {
        for (j = 0; j < size[i]; ++j) {
            hash_i = hash_index[i] + j;
            end = hit_h[hash_i] >> 32;
            period = hit_h[hash_i] & _32mask;
            start = end - period;
            dp[i][j] = (self_dp_t){total_n, 0, start, end, (double)(period), 0.0, 0};
        }
    }
}

// TODO INS/DEL connection
double get_con_score(int cur_start, int cur_end, double cur_period, int pre_start, int pre_end, double pre_period, double *period, int k) {
    double f = cur_period / pre_period;
    double s = (double)(MIN_OF_TWO(cur_end - pre_end, k));
    int embed = (cur_start <= pre_start);

    if (0.9 <= f && f <= 1.1) { // similar period
        if (embed) return 0;
        *period = cur_period;
    } else if ((int)(f+0.5) == 0) {
        *period = cur_period;
        s = 0;
    } else if (fabs((int)(f+0.5) - f) <= 0.1) { // x-fold period
        return 0;
        *period = pre_period;
        s /= f;
    } else { 
        *period = cur_period;
        s = 0;
    }
    return s;
}

int remove_alone_hit(self_dp_t **dp, chain_t ch) {
    if (ch.len < 3) return 0;
    int i;
    int pre_i = ch.chain[0].i, pre_j = ch.chain[0].j; int pre_end = dp[pre_i][pre_j].end;
    int cur_i = ch.chain[1].i, cur_j = ch.chain[1].j; int cur_end = dp[cur_i][cur_j].end;
    int next_i, next_j, next_end;
    double cur_p;
    for (i = 1; i < ch.len-1; ++i) {
        next_i = ch.chain[i+1].i, next_j = ch.chain[i+1].j;
        next_end = dp[next_i][next_j].end;

        cur_p = dp[cur_i][cur_j].period;
        if ((double)(next_end - pre_end) > cur_p) {
            printf("remove: %d, %lf %d\n", cur_end, cur_p, cur_end - dp[cur_i][cur_j].start);
        }

        pre_end = cur_end; cur_end = next_end;
        cur_i = next_i; cur_j = next_j;
    }

    return 0;
}

// TODO allocate DP matrix uniformly
chain_t self_dp_chain(hash_t *hit_h, int hit_n, int kmer_k, self_dp_t ***_dp, int *tot_N) {
    int i, j, k;

    // calculate DP matrix size, allocate DP matrix
    int tot_n = 1, *array_size;
    int *hash_index;
    for (i = 1; i < hit_n; ++i) {
        if (hit_h[i] >> 32 != hit_h[i-1] >> 32) {
            tot_n += 1;
        }
    }
    *tot_N = tot_n;
    array_size = (int*)_err_malloc(sizeof(int) * tot_n);
    hash_index = (int*)_err_malloc(sizeof(int) * tot_n);
    *_dp = (self_dp_t**)_err_calloc((tot_n+1), sizeof(self_dp_t*));
    self_dp_t **dp = *_dp;
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
    int cur_i, cur_j, pre_i, pre_j, max_pre_i, max_pre_j, iter_n;
    double score, max_score, con_score, con_period, max_period, cur_p, pre_p;
    int max_h = 100;
    printf("%d\n", tot_n);
    for (cur_i = 1; cur_i < tot_n; ++cur_i) {
        for (cur_j = 0; cur_j < array_size[cur_i]; ++cur_j) {
            cur_p = dp[cur_i][cur_j].period;
            printf("%d, %d\n", dp[cur_i][cur_j].start, dp[cur_i][cur_j].end);
            if (dp[cur_i][cur_j].start == 3293 && dp[cur_i][cur_j].end == 3718)
                printf("debug");
            max_score = 0, iter_n = 0;
            for (pre_i = cur_i-1; pre_i >= 0; --pre_i) {
                if (dp[pre_i][0].end < dp[cur_i][cur_j].start) goto UPDATE;
                for (pre_j = 0; pre_j < array_size[pre_i]; ++pre_j) {
                    pre_p = dp[pre_i][pre_j].period;
                    if ((con_score = get_con_score(dp[cur_i][cur_i].start, dp[cur_i][0].end, cur_p, dp[pre_i][pre_j].start, dp[pre_i][0].end, pre_p, &con_period, kmer_k)) == 0) continue;
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
            if (max_score > 0) {
                dp[cur_i][cur_j].from_i = max_pre_i;
                dp[cur_i][cur_j].from_j = max_pre_j;
                dp[cur_i][cur_j].period = max_period;
                dp[cur_i][cur_j].score = max_score;
            }
        }
    }

    // backtrack, obtain top N chains
    dp_score_t *score_rank = (dp_score_t*)_err_malloc(hit_n * sizeof(dp_score_t));
    int score_n = sort_dp_score(dp, array_size, tot_n, score_rank);
    int top_N = 10000, ch_n, chain_len;
    chain_t *ch = (chain_t*)_err_malloc(top_N * sizeof(chain_t));
    for (i = 0; i < top_N; ++i) {
        ch[i].chain = (cell_t*)_err_malloc(tot_n * sizeof(cell_t));
    }
    i = 0, ch_n = 0;
    while (i < score_n && ch_n < top_N) {
        chain_len = backtrack_dp(dp, tot_n, score_rank[i].i, score_rank[i].j, ch+ch_n);
        if (chain_len > 0) ++ch_n;
        ++i;
    }
    for (i = 0; i < ch_n; ++i) {
        if (ch[i].len > 50) {
            int start_i = ch[i].chain[0].i, start_j = ch[i].chain[0].j, end_i = ch[i].chain[ch[i].len-1].i, end_j = ch[i].chain[ch[i].len-1].j;
            printf("%d: score: %lf, len: %d (%ld,%lf,%ld) -> (%ld,%lf,%ld)\n", i+1, dp[end_i][end_j].score, ch[i].len, hit_h[hash_index[start_i] + start_j] >> 32, dp[start_i][start_j].period, hit_h[hash_index[start_i] + start_j] & _32mask, hit_h[hash_index[end_i] + end_j] >> 32, dp[end_i][end_j].period, hit_h[hash_index[end_i] + end_j] & _32mask);
        }
    }

    // post-process of N chains
    // 1. remove stand-alone hit
    ch_n = 1;
    for (i = 0; i < ch_n; ++i) {
        remove_alone_hit(dp, ch[i]);
    }
    // 2. merge into larger chain that contains INS/DEL

    chain_t ret_ch;
    ret_ch.len = ch[0].len;
    ret_ch.chain = (cell_t*)_err_malloc(ret_ch.len * sizeof(cell_t));
    for (i = 0; i < ret_ch.len; ++i) {
        ret_ch.chain[i].i = ch[0].chain[i].i;
        ret_ch.chain[i].j = ch[0].chain[i].j;
    }

    // for (i = 0; i <= tot_n; ++i) free(dp[i]); free(dp); 
    free(array_size); free(hash_index);
    for (i = 0; i < top_N; ++i) free(ch[i].chain); free(ch); free(score_rank);
    return ret_ch;
}

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
        if (copy_num == 0) 
            err_fatal(__func__, "Unexpected copy number. (%d, %d, %lf)\n", pos_array[i], pos_array[i+1], period);

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

int partition_seqs(char *seq, self_dp_t **dp, chain_t ch, int *par_pos) {
    int i, j, start, end, seq_len;
    self_dp_t dp_cell;
    double period = dp[ch.chain[0].i][ch.chain[0].j].period;
    seq_len = strlen(seq);
    int array_size = seq_len; // (int)(period) * 2;
    int8_t **hit_array = (int8_t**)_err_malloc(sizeof(int8_t*) * array_size);
    for (i = 0; i < array_size; ++i) hit_array[i] = (int8_t*)_err_calloc(seq_len, sizeof(int8_t));
    int hit_i = 0, hit;
    // fill hit array
    for (i = 0; i < ch.len; ++i) {
        dp_cell = dp[ch.chain[i].i][ch.chain[i].j];
        start = dp_cell.start, end = dp_cell.end;
        printf("start: %d, end: %d\n", start, end);
        hit = 0;
        for (j = 0; j < hit_i; ++j) {
            if (hit_array[j][start]) {
                hit_array[j][end] = 1;
                hit = 1;
                break;
            } else if (hit_array[j][end]) {
                hit_array[j][start] = 1;
                hit = 1;
                break;
            }
        }
        if (hit == 0) {
            printf("%d: (%d, %d)\n", hit_i, start, end);
            hit_array[hit_i][start] = 1;
            hit_array[hit_i++][end] = 1;
        } else {
            printf("%d: (%d, %d)\n", j, start, end);
        }
    }
    // select kmers with max hit 
    int hit_n, max_hit_n = 0, max_i;
    for (i = 0; i < hit_i; ++i) {
        hit_n = 0;
        for (j = 0; j < seq_len; ++j) hit_n += hit_array[i][j];
        if (hit_n > max_hit_n) {
            max_hit_n = hit_n;
            max_i = i;
        }
    }
    int par_n = partition_seqs_core(seq, hit_array[max_i], period, par_pos);

    for (i = 0; i < array_size; ++i) free(hit_array[i]); free(hit_array);
    return par_n;
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

    hash_t *hit_h; self_dp_t **dp; int tot_n;
    int hit_n = collect_hash_hit(h, hn, &hit_h); free(h);
    chain_t ch = self_dp_chain(hit_h, hit_n, mtp->k, &dp, &tot_n);
    int p = (int)(dp[ch.chain[0].i][ch.chain[0].j].period + 0.5);
    int *par_pos = (int*)_err_malloc(seq_len / p * 2 * sizeof(int));
    int par_n = partition_seqs(seq, dp, ch, par_pos);
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
    spoa_msa(seqs, seq_i, cons_seq);
    printf("cons: %s\n", cons_seq);
    
    free(hit_h); free(bseq); free(ch.chain); free(par_pos);
    for (i = 0; i <= tot_n; ++i) free(dp[i]); free(dp); 
    for (i = 0; i < seq_i; ++i) free(seqs[i]); free(seqs); free(cons_seq);
    return 0;
}
