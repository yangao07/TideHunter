#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "tide_hunter.h"
#include "tandem_hit.h"
#include "tandem_chain.h"
#include "partition.h"
#include "gen_cons.h"
#include "utils.h"
#include "seq.h"

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
int tide_hunter_core(kseq_t *read_seq, tandem_seq_t *tseq, mini_tandem_para *mtp) {
    if ((int)(read_seq->seq.l) < mtp->k) return 0;
    int seq_len = read_seq->seq.l; char *seq = read_seq->seq.s; 
    uint8_t *bseq = get_bseq(seq, seq_len);

    // collect tandem repeat hits
    hash_t *hit_h;
    int hit_n = collect_tandem_repeat_hit(bseq, seq_len, mtp, &hit_h);
    // chaining by DP
    dp_t **dp; int tot_n=0; chain_t *chain; int ch_m=0;
    int ch_n = tandem_chain(seq_len, hit_h, hit_n, mtp, &dp, &tot_n, &chain, &ch_m); free(hit_h);
    int ch_i, i; chain_t ch;
    for (ch_i = 0; ch_i < ch_n; ++ch_i) {
        // partition seq into segments
        ch = chain[ch_i];
        int par_n, *par_pos;
        par_pos = get_partition_pos_with_narrow_global_alignment(bseq, seq_len, dp, ch, mtp, &par_n);
        if (par_n < mtp->min_copy+1) {
            free(par_pos); continue;
        }
        // msa and generate consensus
        seqs_msa(seq_len, bseq, par_n, par_pos, tseq, mtp);
    }

    free(bseq);
    if (ch_m > 0) {
        for (i = 0; i < ch_m; ++i) free(chain[i].cell); free(chain);
    }
    if (tot_n > 0) {
        for (i = 0; i <= tot_n; ++i) free(dp[i]); free(dp); 
    }
    return 0;
}
