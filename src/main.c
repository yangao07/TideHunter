#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <pthread.h>
#include <math.h>
#include "mini_tandem.h"
#include "self_chain.h"
#include "utils.h"
#include "kseq.h"
#include "seq.h"

const char PROG[20] = "miniTandem";
const char VERSION[20] = "1.0.0";
//const char DATE[20] = "2018-12-01";
const char CONTACT[30] = "yangaoucla@gmail.com";

const struct option mini_tandem_opt [] = {
    { "kmer-length", 1, NULL, 'k' },
    { "window-size", 1, NULL, 'w' },
    { "step-size", 1, NULL, 's' },
    { "minimal-m", 1, NULL, 'm' },
    { "HPC-kmer", 0, NULL, 'H' },
    { "max-diverg", 1, NULL, 'e' },
    { "min-period", 1, NULL, 'p' },
    { "max-period", 1, NULL, 'P' },

    { "rep-range", 1, NULL, 'r' },

    { "splint-seq", 1, NULL, 'S' },
    { "detail-out", 1, NULL, 'd' },

    { "thread", 1, NULL, 't' },
    { 0, 0, 0, 0}
};

static int usage(void)
{
    err_printf("\n");
	err_printf("%s: Tandem repeat detection and consensus calling from noisy concatemeric long-read\n\n", PROG);

    time_t t; time(&t);
    err_printf("Version: %s Build date: %s", VERSION, ctime(&t));
    err_printf("Contact: %s\n\n", CONTACT);

    err_printf("Usage:   %s [options] in.fa/fq > cons_out.fastq\n\n", PROG);

	err_printf("Options: \n");
    err_printf("         -t --thread      [INT]    number of threads to use. [%d]\n", THREAD_N);
    err_printf("         -k --kmer-length [INT]    k-mer length (no larger than 16). [%d]\n", KMER_SIZE); // TODO largest kmer len
    err_printf("         -w --window-size [INT]    window size. [%d]\n", KMER_WSIZE);
    err_printf("         -s --step-size   [INT]    step size. [%d]\n", KMER_SSIZE);
    err_printf("         -m --minimal-m   [INT]    number of minimal k-mer to keep in each window. [%d]\n", KMER_MINM);
    err_printf("         -e --max-diverg  [INT]    maximum allowed divergence rate between two consecutive repeats. [%.2f]\n", MAX_DIV);
    err_printf("         -H --HPC-kmer             use homopolymer-compressed k-mer. [False]\n");
    err_printf("         -p --min-period  [INT]    minimum period size of tandem repeat. (>=%d) [%d]\n", 2, MIN_PERIOD);
    err_printf("         -P --max-period  [INT]    maximum period size of tandem repeat. (<=%d) [%d]\n", MAX_PERIOD, MAX_PERIOD);

//  err_printf("         -r --rep-range   [INT]    maximum range to find tandem repeat. [%d]\n", REP_RANGE); 
//  err_printf("                                   (-1: no limit, tandem repeat can span the whole sequence)\n");

    err_printf("         -S --splint-seq  [STR]    splint sequence in FASTA/FASTQ format. [NULL]\n");
    err_printf("         -d --detail-out  [STR]    detailed information of each consensus. [NULL]\n");
    err_printf("                                   (start, end, score, etc.)\n");

	err_printf("\n");
	return 1;
}
//KSEQ_INIT(gzFile, gzread)

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

thread_aux_t *aux_init(mini_tandem_para *mtp) {
    int i;
    thread_aux_t *aux = (thread_aux_t*)calloc(mtp->n_thread, sizeof(thread_aux_t));
    for (i = 0; i < mtp->n_thread; ++i) {
        aux[i].tid = i; 
        aux[i].mtp = mtp;
    }
    return aux;
}

void aux_free(thread_aux_t *aux, int n_thread) {
    //int i;
    //for (i = 0; i < n_thread; ++i) {
    //}
    free(aux);
}

int COUNT=0;
int THREAD_READ_I;
pthread_rwlock_t RWLOCK;

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

static void *mini_tandem_thread_main(void *aux)
{
    thread_aux_t *a = (thread_aux_t*)aux;
    int i = 0;
    while (1) {
        pthread_rwlock_wrlock(&RWLOCK);
        i = THREAD_READ_I++;
        pthread_rwlock_unlock(&RWLOCK);
        if (i >= a->n_seqs) break;
        mini_tandem_para *mtp = a->mtp;
        kseq_t *read_seq = a->read_seq + i; 
        tandem_seq_t *tandem_seq = a->tseq + i;
        // generate cons_seq from seq , cons_seq may have multiple seqs
        mini_tandem_core(read_seq, tandem_seq, mtp);
    }
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

static inline double get_div_exp(int k, double div) {
    return exp(k * div);
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
    mtp->div_exp = get_div_exp(KMER_SIZE, MAX_DIV);
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
    if (mtp->n_thread < 1) mtp->n_thread = 1;
    thread_aux_t *aux = aux_init(mtp);

    // core loop
    pthread_rwlock_init(&RWLOCK, NULL);
    while ((n_seqs = mini_tandem_read_seq(read_seq, CHUNK_READ_N)) != 0) {
        if (mtp->n_thread <= 1) {
            aux->n_seqs = n_seqs;
            aux->read_seq = read_seq;
            aux->tseq = tseq;
            mini_tandem_thread_main(aux);
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
    pthread_rwlock_destroy(&RWLOCK);

    // free variables
    for (i = 0; i < CHUNK_READ_N; ++i) {
        free((read_seq+i)->name.s); free((read_seq+i)->comment.s); free((read_seq+i)->seq.s); free((read_seq+i)->qual.s);
        seq_t *cons_seq = tseq[i].cons_seq;
        free(cons_seq->name.s); free(cons_seq->comment.s); free(cons_seq->seq.s); free(cons_seq->qual.s); 
        free(tseq[i].cons_seq); free(tseq[i].cons_start); free(tseq[i].cons_end); free(tseq[i].cons_len); free(tseq[i].cons_score);
    } free(read_seq); free(tseq); ks_destroy(fs); err_gzclose(readfp);
    aux_free(aux, mtp->n_thread);
    double sys_t, usr_t; usr_sys_cputime(&usr_t, &sys_t); err_func_printf(__func__, "User: %.3f sec; Sys: %.3f sec.\n", usr_t, sys_t);
    return 0;
}

// TODO add score para 
int main(int argc, char *argv[])
{
    mini_tandem_para *mtp = mini_tandem_init_para();
    int c;
    while ((c = getopt_long(argc, argv, "k:w:m:Hs:r:e:p:P:S:d:t:",mini_tandem_opt, NULL)) >= 0) {
        switch(c)
        {
            case 'k': mtp->k = atoi(optarg); break;
            case 'w': mtp->w = atoi(optarg); break;
            case 'm': mtp->m = atoi(optarg); break;
            case 's': mtp->s = atoi(optarg); break;
            case 'H': mtp->hpc = 1; break;
            case 'e': mtp->max_div = atof(optarg); break;
            case 'p': mtp->min_p = atoi(optarg);
                      if (mtp->min_p < MIN_PERIOD) {
                          err_printf("Error: -p --min-period(%d) needs to be >= %d.\n", mtp->min_p, MIN_PERIOD); 
                          goto End;
                      }
                      break;
            case 'P': mtp->max_p = atoi(optarg); 
                      if (mtp->max_p > MAX_PERIOD) {
                          err_printf("Error: -P --max-period(%d) needs to be <= %d.\n", mtp->max_p, MAX_PERIOD); 
                          goto End;
                      }
                      break;
          //case 'r': mtp->max_range = atoi(optarg); break;
            case 'd': mtp->detail_fp = xopen(optarg, "w"); break;
            case 'S': mtp->splint_fn = strdup(optarg); break;
            case 't': mtp->n_thread = atoi(optarg); break;
            default:
                      err_printf("Error: unknown option: -%c %s.\n", c, optarg);
                      goto End;
        }
    }
	if (argc < 2) return usage();

    mtp->div_exp = get_div_exp(mtp->k, mtp->max_div);
    mini_tandem(argv[optind], mtp);
End:
    mini_tandem_free_para(mtp);
    double sys_t, usr_t; usr_sys_cputime(&usr_t, &sys_t); err_func_printf(__func__, "User: %.3f sec; Sys: %.3f sec.\n", usr_t, sys_t);
    return 0;
}
