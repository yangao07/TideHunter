#include <stdio.h>
#include <stdlib.h>
#include "mini_tandem.h"
#include "utils.h"
#include "kseq.h"
#include <pthread.h>

KSEQ_INIT(gzFile, gzread)

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
    mini_tandem_para *ncp;
    int n_seqs; kseq_t *read_seq, *cons_seq;
} thread_aux_t;

int mini_tandem_main(thread_aux_t *aux)
{
    int i = 0;
    while (1) {
        pthread_rwlock_wrlock(&RWLOCK);
        i = THREAD_READ_I++;
        pthread_rwlock_unlock(&RWLOCK);
        if (i >= aux->n_seqs) break;
        kseq_t *seq = aux->read_seq + i; mini_tandem_para *ncp = aux->ncp;
        kseq_t *cons_seq = aux->cons_seq + i;
    }
    return 0;
}

static void *mini_tandem_thread_main(void *aux)
{
    thread_aux_t *a = (thread_aux_t*)aux;
    return mini_tandem_main(a);
}

int mini_tandem_core(const char *read_fn, mini_tandem_para *ncp)
{
    int i, n_seqs, THREAD_READ_I = 0;
    gzFile readfp = xzopen(read_fn, "r");
    kstream_t *fs = ks_init(readfp);
    kseq_t *read_seq = (kseq_t*)calloc(CHUNK_READ_N, sizeof(kseq_t));
    kseq_t *cons_seq = (kseq_t*)calloc(CHUNK_READ_N, sizeof(kseq_t));
    for (i = 0; i < CHUNK_READ_N; ++i) read_seq[i].f = fs;

    // alloc and initialization for auxiliary data
    thread_aux_t *aux;
    if (ncp->n_thread < 1) ncp->n_thread = 1;
    aux = (thread_aux_t*)calloc(ncp->n_thread, sizeof(thread_aux_t));
    for (i = 0; i < ncp->n_thread; ++i) {
        aux[i].tid = i; 
        aux[i].ncp = ncp;
    }
    pthread_rwlock_init(&RWLOCK, NULL);

    // core loop
    while ((n_seqs = mini_tandem_read_seq(read_seq, CHUNK_READ_N)) != 0) {
        if (ncp->n_thread <= 1) {
            aux->n_seqs = n_seqs;
            aux->read_seq = read_seq;
            aux->cons_seq = cons_seq;
            mini_tandem_main(aux);
        } else {
            pthread_t *tid; pthread_attr_t attr; 
            pthread_attr_init(&attr); pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
            tid = (pthread_t*)calloc(ncp->n_thread, sizeof(pthread_t));
            int j;
            for (j = 0; j < ncp->n_thread; ++j) {
                aux[j].n_seqs = n_seqs; 
                aux[j].read_seq = read_seq;
                aux[j].cons_seq = cons_seq;
                pthread_create(&tid[j], &attr, mini_tandem_thread_main, aux+j);
            }
            for (j = 0; j < ncp->n_thread; ++j) pthread_join(tid[j], 0);
            free(tid);
        }
        // output initial consensus sequences
        mini_tandem_output(n_seqs, ncp);

    }
    // free 
    pthread_rwlock_destroy(&RWLOCK);
    free(aux);
    for (i = 0; i < CHUNK_READ_N; ++i) {
        free((read_seq+i)->name.s); free((read_seq+i)->comment.s); free((read_seq+i)->seq.s); free((read_seq+i)->qual.s);
        free((cons_seq+i)->name.s); free((cons_seq+i)->comment.s); free((cons_seq+i)->seq.s); free((cons_seq+i)->qual.s);
    } free(read_seq); ks_destroy(fs); gzclose(readfp); 
    return 0;
}
