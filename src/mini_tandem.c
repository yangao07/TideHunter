#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <math.h>
#include "mini_tandem.h"
#include "self_chain.h"
#include "utils.h"
#include "kseq.h"
#include "seq.h"

KSEQ_INIT(gzFile, gzread)

int get_seq_from_fx(gzFile fp, char **seq) {
    kstream_t *fs = ks_init(fp);
    kseq_t *read_seq = (kseq_t*)calloc(1, sizeof(kseq_t));
    read_seq->f = fs;
    int len;
    if (kseq_read(read_seq) > 0) {
        (*seq) = strdup(read_seq->seq.s);
        len = read_seq->seq.l;
        free(read_seq); ks_destroy(fs);
        return len;
    } else {
        err_func_format_printf(__func__, "Warning: No splint sequence found.\n");
        return 0;
    }
}

int mini_tandem_read_seq(kseq_t *read_seq, int chunk_read_n)
{
    kseq_t *s = read_seq;
    int n = 0;
    while (kseq_read(s+n) >= 0) {
        n++;
        if (n >= chunk_read_n) break;
    }
    return n;
}

int COUNT=0;
int THREAD_READ_I;
pthread_rwlock_t RWLOCK;

typedef struct {
    int tid;
    mini_tandem_para *mtp;
    int8_t *hit_array; int array_m;
    int n_seqs; 
    kseq_t *read_seq;
    tandem_seq_t *tseq;
} thread_aux_t;

// 1. output cons.fastq
// 2. output cons.info
void mini_tandem_output(int n_seqs, kseq_t *read_seq, tandem_seq_t *tseq, mini_tandem_para *mtp) {
    int i, seq_i, cons_i, cons_seq_start = 0, cons_seq_end = 0;
    tandem_seq_t *_tseq;
    for (seq_i = 0; seq_i < n_seqs; ++seq_i) {
        _tseq = tseq + seq_i;
        for (cons_i = 0; cons_i < _tseq->cons_n; ++cons_i) { // TODO cons sorted by start,end
            fprintf(stdout, ">%s_cons%d %d-%d:%d\n", (read_seq+seq_i)->name.s, cons_i, _tseq->cons_start[cons_i], _tseq->cons_end[cons_i], _tseq->cons_len[cons_i]);
            cons_seq_end += (tseq+seq_i)->cons_len[cons_i];
            for (i = cons_seq_start; i < cons_seq_end; ++i)  fprintf(stdout, "%c", _tseq->cons_seq->seq.s[i]);
            cons_seq_start += _tseq->cons_len[cons_i];
            fprintf(stdout, "\n");
        }
        _tseq->cons_n = 0;
        _tseq->cons_seq->seq.l = 0;
        cons_seq_start = cons_seq_end = 0;
    }
    if (mtp->detail_fp != NULL) {
    }
}

int mini_tandem_core(kseq_t *read_seq, tandem_seq_t *tandem_seq, int8_t **hit_array, int *array_m, mini_tandem_para *mtp) {
    // 0. build Hash index, generate coordinate of all the hits (x,y)
    // 1. chaining by partition based on (y-x) values
    // 2. keep top N chains, (calcuate density of each chain: tot_N_hits / tot_N_kmer)
    // 3. call consensus with each chain
    // 4. polish consensus result
    hash_partition(read_seq->seq.s, read_seq->seq.l, tandem_seq, hit_array, array_m, mtp);
    return 0;
}

void mini_tandem_main(thread_aux_t *aux)
{
    int i = 0;
    while (1) {
        pthread_rwlock_wrlock(&RWLOCK);
        i = THREAD_READ_I++;
        pthread_rwlock_unlock(&RWLOCK);
        if (i >= aux->n_seqs) break;
        mini_tandem_para *mtp = aux->mtp;
        kseq_t *read_seq = aux->read_seq + i; 
        tandem_seq_t *tandem_seq = aux->tseq + i;
        int8_t *hit_array = aux->hit_array; int *array_m = &(aux->array_m);
        // generate cons_seq from seq , cons_seq may have multiple seqs
        mini_tandem_core(read_seq, tandem_seq, &hit_array, array_m, mtp);
    }
}

static void *mini_tandem_thread_main(void *aux)
{
    thread_aux_t *a = (thread_aux_t*)aux;
    mini_tandem_main(a);
    return aux;
}

tandem_seq_t *alloc_tandem_seq(int n) {
    tandem_seq_t *tseq = (tandem_seq_t*)_err_malloc(n * sizeof(tandem_seq_t));
    int i;
    for (i = 0; i < n; ++i) {
        tseq[i].cons_seq = (seq_t*)calloc(1, sizeof(seq_t));
        tseq[i].cons_n = 0; tseq[i].cons_m = 1;
        tseq[i].cons_start = (int*)_err_malloc(sizeof(int));
        tseq[i].cons_end = (int*)_err_malloc(sizeof(int));
        tseq[i].cons_len = (int*)_err_malloc(sizeof(int));
        tseq[i].cons_score = (int*)_err_malloc(sizeof(int));
    }
    return tseq;
}

mini_tandem_para *mini_tandem_init_para(void) {
    mini_tandem_para *mtp = (mini_tandem_para*)_err_malloc(sizeof(mini_tandem_para));
    mtp->n_thread = THREAD_N;

    mtp->splint_fn = NULL; 
    mtp->splint_seq = NULL; mtp->splint_rc_seq = NULL;
    mtp->splint_len = 0;
    mtp->detail_fp = NULL;
    mtp->k = KMER_SIZE;
    mtp->w = KMER_WSIZE;
    mtp->s = KMER_SSIZE;
    mtp->m = KMER_MINM;
    mtp->hpc = 0;
    mtp->max_div = MAX_DIV;
    mtp->div_exp = exp(mtp->k * MAX_DIV);
    mtp->min_p = MIN_PERIOD;
    mtp->max_p = MAX_PERIOD;


    mtp->max_range = REP_RANGE;

    return mtp;
}

void mini_tandem_free_para(mini_tandem_para *mtp) {
    if (mtp->detail_fp != NULL) err_fclose(mtp->detail_fp);
    if (mtp->splint_fn != NULL) free(mtp->splint_fn);
    if (mtp->splint_seq != NULL) free(mtp->splint_seq);
    if (mtp->splint_rc_seq != NULL) free(mtp->splint_rc_seq);
    free(mtp);
}

int mini_tandem(const char *read_fn, mini_tandem_para *mtp)
{
    int i, n_seqs, THREAD_READ_I = 0;
    gzFile readfp = xzopen(read_fn, "r");
    kstream_t *fs = ks_init(readfp);
    kseq_t *read_seq = (kseq_t*)calloc(CHUNK_READ_N, sizeof(kseq_t));
    for (i = 0; i < CHUNK_READ_N; ++i) read_seq[i].f = fs;
    tandem_seq_t *tseq = alloc_tandem_seq(CHUNK_READ_N);

    if (mtp->splint_fn != NULL) {
        gzFile splint_fp = xzopen(mtp->splint_fn, "r");
        mtp->splint_len = get_seq_from_fx(splint_fp, &(mtp->splint_seq));
        mtp->splint_rc_seq = get_rc_seq(mtp->splint_seq, mtp->splint_len);
        err_gzclose(splint_fp);
    }

    // alloc and initialization for auxiliary data
    thread_aux_t *aux;
    if (mtp->n_thread < 1) mtp->n_thread = 1;
    aux = (thread_aux_t*)calloc(mtp->n_thread, sizeof(thread_aux_t));
    for (i = 0; i < mtp->n_thread; ++i) {
        aux[i].tid = i; 
        aux[i].mtp = mtp;
        aux[i].hit_array = (int8_t*)_err_malloc(sizeof(int8_t)*10000*10000);
        aux[i].array_m = 10000;
    }
    pthread_rwlock_init(&RWLOCK, NULL);

    // core loop
    while ((n_seqs = mini_tandem_read_seq(read_seq, CHUNK_READ_N)) != 0) {
        if (mtp->n_thread <= 1) {
            aux->n_seqs = n_seqs;
            aux->read_seq = read_seq;
            aux->tseq = tseq;
            mini_tandem_main(aux);
        } else {
            pthread_t *tid; pthread_attr_t attr; 
            pthread_attr_init(&attr); pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
            tid = (pthread_t*)calloc(mtp->n_thread, sizeof(pthread_t));
            int j;
            for (j = 0; j < mtp->n_thread; ++j) {
                aux[j].n_seqs = n_seqs; 
                aux[j].read_seq = read_seq;
                aux[j].tseq = tseq;
                pthread_create(&tid[j], &attr, &mini_tandem_thread_main, aux+j);
            }
            for (j = 0; j < mtp->n_thread; ++j) pthread_join(tid[j], 0);
            free(tid);
        }
        // output initial consensus sequences
        mini_tandem_output(n_seqs, read_seq, tseq, mtp);
        double sys_t, usr_t; usr_sys_cputime(&usr_t, &sys_t); err_func_printf(__func__, "User: %.3f sec; Sys: %.3f sec.\n", usr_t, sys_t);
    }
    // free
    for (i = 0; i < mtp->n_thread; ++i) free(aux[i].hit_array);
    pthread_rwlock_destroy(&RWLOCK);
    free(aux);
    for (i = 0; i < CHUNK_READ_N; ++i) {
        free((read_seq+i)->name.s); free((read_seq+i)->comment.s); free((read_seq+i)->seq.s); free((read_seq+i)->qual.s);
        seq_t *cons_seq = tseq[i].cons_seq;
        free(cons_seq->name.s); free(cons_seq->comment.s); free(cons_seq->seq.s); free(cons_seq->qual.s); 
        free(tseq[i].cons_seq); free(tseq[i].cons_start); free(tseq[i].cons_end); free(tseq[i].cons_len); free(tseq[i].cons_score);
    } free(read_seq); free(tseq); ks_destroy(fs); err_gzclose(readfp);
    double sys_t, usr_t; usr_sys_cputime(&usr_t, &sys_t); err_func_printf(__func__, "User: %.3f sec; Sys: %.3f sec.\n", usr_t, sys_t);
    return 0;
}
