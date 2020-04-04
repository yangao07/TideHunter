#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <pthread.h>
#include <math.h>
#include <sys/sysinfo.h>
#include "tide_hunter.h"
#include "utils.h"
#include "kseq.h"

const char PROG[20] = "TideHunter";
const char VERSION[20] = "1.2.2";
const char CONTACT[30] = "yangao07@hit.edu.cn";

const struct option mini_tandem_opt [] = {
	{ "kmer-length", 1, NULL, 'k' },
	{ "window-size", 1, NULL, 'w' },
	{ "step-size", 1, NULL, 's' },
	// { "minimal-m", 1, NULL, 'm' },
	{ "HPC-kmer", 0, NULL, 'H' },

	{ "min-copy", 1, NULL, 'c' },
	{ "max-diverg", 1, NULL, 'e' },
	{ "min-period", 1, NULL, 'p' },
	{ "max-period", 1, NULL, 'P' },

	// { "rep-range", 1, NULL, 'r' },

	// { "splint-seq", 1, NULL, 'S' },
	{ "five-prime", 1, NULL, '5' },
	{ "three-prime", 1, NULL, '3' },
	{ "ada-match-rat", 1, NULL, 'a' },

	{ "cons-out", 1, NULL, 'o' },
	{ "longest", 0, NULL, 'l' },
	{ "full-len", 0, NULL, 'l' },
	{ "out-fmt", 1, NULL, 'f' },
	// { "detail-out", 1, NULL, 'd' },

	{ "thread", 1, NULL, 't' },

    { "help", 0, NULL, 'h' },
    { "version", 0, NULL, 'v' },
	{ 0, 0, 0, 0}
};

static inline int64_t th_parse_num(const char *str)
{
        double x;
        char *p;
        x = strtod(str, &p);
        if (*p == 'G' || *p == 'g') x *= 1e9;
        else if (*p == 'M' || *p == 'm') x *= 1e6;
        else if (*p == 'K' || *p == 'k') x *= 1e3;
        return (int64_t)(x + .499);
}

static int usage(void)
{
	err_printf("\n");
	err_printf("%s: Tandem repeats detection and consensus calling from noisy long reads\n\n", PROG);

	time_t t; time(&t);
	err_printf("Version: %s\n", VERSION);
	//err_printf("Build date: %s", ctime(&t));
	err_printf("Contact: %s\n\n", CONTACT);

	err_printf("Usage:   %s [options] in.fa/fq > cons.fa\n\n", PROG);

	err_printf("Options: \n");
	err_printf("    Seeding:\n");
	err_printf("         -k --kmer-length [INT]    k-mer length (no larger than %d). [%d]\n", MAX_KMER_SIZE, KMER_SIZE);
	err_printf("         -w --window-size [INT]    window size. [%d]\n", KMER_WSIZE);
	err_printf("         -s --step-size   [INT]    step size. [%d]\n", KMER_SSIZE);
	//err_printf("         -m --minimal-m   [INT]    number of minimal k-mer to keep in each window. [%d]\n", KMER_MINM);
	err_printf("         -H --HPC-kmer             use homopolymer-compressed k-mer. [False]\n");

	//err_printf("\n");

	err_printf("    Tandem repeat criteria:\n");
	// TODO min_copy < 2 ???
	err_printf("         -c --min-copy    [INT]    minimum copy number of tandem repeat. [%d]\n", MIN_COPY);
	err_printf("         -e --max-diverg  [INT]    maximum allowed divergence rate between two consecutive repeats. [%.2f]\n", MAX_DIV);
	err_printf("         -p --min-period  [INT]    minimum period size of tandem repeat. (>=%u) [%u]\n", MIN_PERIOD, DEF_MIN_PERIOD);
	err_printf("         -P --max-period  [INT]    maximum period size of tandem repeat. (<=%u) [%s]\n", MAX_PERIOD, DEF_MAX_PERIOD_STR);

//  err_printf("         -r --rep-range   [INT]    maximum range to find tandem repeat. [%d]\n", REP_RANGE); 
//  err_printf("                                   (-1: no limit, tandem repeat can span the whole sequence)\n");
	//err_printf("\n");

	err_printf("    Adapter sequence:\n");
	err_printf("         -5 --five-prime  [STR]    5' adapter sequence (sense strand). [NULL]\n");
	err_printf("         -3 --three-prime [STR]    3' adapter sequence (anti-sense strand). [NULL]\n");
	err_printf("         -a --ada-mat-rat [FLT]    minimum match ratio of adapter sequence. [%.2f]\n", ADA_MATCH_RAT);

	//err_printf("\n");

	err_printf("    Output:\n");
	err_printf("         -o --cons-out    [STR]    output file. [stdout]\n");
	err_printf("         -l --longest              only output the consensus of the longest tandem repeat. [False]\n");
	err_printf("         -F --full-len             only output the consensus that is full-length. [False]\n");
	err_printf("         -f --out-fmt     [INT]    output format. [%d]\n", FASTA_FMT);
	err_printf("                                       %d: FASTA\n", FASTA_FMT);
	err_printf("                                       %d: Tabular\n", TAB_FMT);
	//err_printf("         -S --splint-seq  [STR]    splint sequence in FASTA/FASTQ format. [NULL]\n");
	//err_printf("         -d --detail-out  [STR]    detailed information of each consensus. [NULL]\n");
	//err_printf("                                   (start, end, score, etc.)\n");

	//err_printf("\n");

	err_printf("    Computing resource:\n");
	err_printf("         -t --thread      [INT]    number of threads to use. [%d]\n\n", MIN_OF_TWO(THREAD_N, get_nprocs()));

    err_printf("    General options:\n");
    err_printf("         -h --help                 print this help usage information.\n");
    err_printf("         -v --version              show version number.\n");

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
		err_func_format_printf(__func__, "No sequence found.\n");
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

void read_seq_free(kseq_t *read_seq) {
    free(read_seq->name.s); 
    free(read_seq->comment.s); 
    free(read_seq->seq.s); 
    free(read_seq->qual.s);
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

int COUNT=0;
int THREAD_READ_I;
pthread_rwlock_t RWLOCK;

// TODO cons.fastq
// TODO cons.stats
void mini_tandem_output(int n_seqs, kseq_t *read_seq, tandem_seq_t *tseq, mini_tandem_para *mtp) {
	int i, seq_i, cons_i, cons_seq_start = 0, cons_seq_end = 0;
	tandem_seq_t *_tseq;
	for (seq_i = 0; seq_i < n_seqs; ++seq_i) {
		_tseq = tseq + seq_i;
		for (cons_i = 0; cons_i < _tseq->cons_n; ++cons_i) { // TODO cons sorted by start,end
			if (mtp->out_fmt == FASTA_FMT) {
				// >readName_consN_readLen:start:end:consLen:copyNum:fullLength
				fprintf(mtp->cons_out, ">%s_cons%d_%d_%d_%d_%d_%.1f_%d_", (read_seq+seq_i)->name.s, cons_i, (int)((read_seq+seq_i)->seq.l),  _tseq->cons_start[cons_i]+1, _tseq->cons_end[cons_i]+1, _tseq->cons_len[cons_i], _tseq->copy_num[cons_i], _tseq->full_length[cons_i]);
                fprintf(mtp->cons_out, "%d", _tseq->sub_pos[cons_i][0]+2);
                for (i = 1; i < _tseq->pos_n[cons_i]-1; ++i) fprintf(mtp->cons_out, ",%d", _tseq->sub_pos[cons_i][i]+2);
                fprintf(mtp->cons_out, ",%d\n", _tseq->sub_pos[cons_i][i]+1);
			} else if (mtp->out_fmt == TAB_FMT) {
				// readName/consN/readLen/start/end/consLen/copyNum/fullLength
				fprintf(mtp->cons_out, "%s\tcons%d\t%d\t%d\t%d\t%d\t%.1f\t%d\t", (read_seq+seq_i)->name.s, cons_i, (int)((read_seq+seq_i)->seq.l),  _tseq->cons_start[cons_i]+1, _tseq->cons_end[cons_i]+1, _tseq->cons_len[cons_i], _tseq->copy_num[cons_i], _tseq->full_length[cons_i]);
                fprintf(mtp->cons_out, "%d", _tseq->sub_pos[cons_i][0]+2);
                for (i = 1; i < _tseq->pos_n[cons_i]-1; ++i) fprintf(mtp->cons_out, ",%d", _tseq->sub_pos[cons_i][i]+2);
                fprintf(mtp->cons_out, ",%d\t", _tseq->sub_pos[cons_i][i]+1);
			}
			cons_seq_end += (tseq+seq_i)->cons_len[cons_i];
			for (i = cons_seq_start; i < cons_seq_end; ++i)  
				fprintf(mtp->cons_out, "%c", _tseq->cons_seq->seq.s[i]);
			cons_seq_start += _tseq->cons_len[cons_i];
			fprintf(mtp->cons_out, "\n");
		}
		_tseq->cons_n = 0;
		_tseq->cons_seq->seq.l = 0;
		cons_seq_start = cons_seq_end = 0;
	}
    // if (mtp->detail_fp != NULL) { }
}

static void *mini_tandem_thread_main(void *aux)
{
	thread_aux_t *a = (thread_aux_t*)aux;
	long i = 0;
	while (1) {
		pthread_rwlock_wrlock(&RWLOCK);
		i = THREAD_READ_I++;
		pthread_rwlock_unlock(&RWLOCK);
		if (i >= a->n_seqs) break;
		mini_tandem_para *mtp = a->mtp;
		kseq_t *read_seq = a->read_seq + i; 
		tandem_seq_t *tandem_seq = a->tseq + i;
        // fprintf(stderr, "%s\n", read_seq->name.s); // for debug
		// generate cons_seq from seq , cons_seq may have multiple seqs
		tide_hunter_core(read_seq, tandem_seq, mtp);
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
		tseq[i].copy_num = (double*)_err_malloc(sizeof(double));
		tseq[i].cons_len = (int*)_err_malloc(sizeof(int));
		tseq[i].full_length = (int8_t*)_err_malloc(sizeof(int8_t));
		tseq[i].cons_score = (int*)_err_malloc(sizeof(int));
        tseq[i].pos_n = (int*)_err_calloc(1, sizeof(int));
        tseq[i].pos_m = (int*)_err_calloc(1, sizeof(int));
        tseq[i].sub_pos = (int**)_err_calloc(1, sizeof(int*));
	}
	return tseq;
}
void tandem_seq_free(tandem_seq_t *tseq) {
    seq_t *cons_seq = tseq->cons_seq;
    free(cons_seq->name.s); free(cons_seq->comment.s); free(cons_seq->seq.s); free(cons_seq->qual.s); 
    free(tseq->cons_seq); 
    free(tseq->cons_start); free(tseq->cons_end); 
    free(tseq->cons_len); free(tseq->copy_num); 
    free(tseq->full_length); free(tseq->cons_score); 
    free(tseq->pos_n); free(tseq->pos_m);
    int i;
    for (i = 0; i < tseq->cons_m; ++i) 
        if(tseq->sub_pos[i]) free(tseq->sub_pos[i]); 
    free(tseq->sub_pos);
}
static inline double get_div_exp(int k, double div) {
	return exp(2 * k * div);
}

mini_tandem_para *mini_tandem_init_para(void) {
	mini_tandem_para *mtp = (mini_tandem_para*)_err_malloc(sizeof(mini_tandem_para));
	mtp->n_thread = THREAD_N;

	// mtp->splint_fn = NULL; mtp->splint_seq = NULL; mtp->splint_rc_seq = NULL; mtp->splint_len = 0;
	mtp->ada_match_rat = ADA_MATCH_RAT;
	mtp->five_fn = NULL; mtp->five_seq = NULL; mtp->five_rc_seq = NULL; mtp->five_len = 0;
	mtp->three_fn = NULL; mtp->three_seq = NULL; mtp->three_rc_seq = NULL; mtp->three_len = 0;
	mtp->k = KMER_SIZE;
	mtp->w = KMER_WSIZE;
	mtp->s = KMER_SSIZE;
	mtp->m = KMER_MINM;
	mtp->hpc = 0;
	mtp->min_copy = MIN_COPY;
	mtp->max_div = MAX_DIV;
	mtp->div_exp = get_div_exp(KMER_SIZE, MAX_DIV);
	mtp->min_p = DEF_MIN_PERIOD;
	mtp->max_p = DEF_MAX_PERIOD;

	mtp->cons_out = stdout;
	mtp->only_longest = 0;
	mtp->only_full_length = 0;
	mtp->out_fmt = FASTA_FMT;
	mtp->detail_fp = NULL;

	mtp->max_range = REP_RANGE;

	return mtp;
}

void mini_tandem_free_para(mini_tandem_para *mtp) {
	if (mtp->detail_fp != NULL) err_fclose(mtp->detail_fp);
	if (mtp->cons_out != stdout) err_fclose(mtp->cons_out);
	//if (mtp->splint_fn != NULL) free(mtp->splint_fn); if (mtp->splint_seq != NULL) free(mtp->splint_seq); if (mtp->splint_rc_seq != NULL) free(mtp->splint_rc_seq);
	if (mtp->five_fn != NULL) free(mtp->five_fn); if (mtp->five_seq != NULL) free(mtp->five_seq); if (mtp->five_rc_seq != NULL) free(mtp->five_rc_seq);
	if (mtp->three_fn != NULL) free(mtp->three_fn); if (mtp->three_seq != NULL) free(mtp->three_seq); if (mtp->three_rc_seq != NULL) free(mtp->three_rc_seq);
	free(mtp);
}

int mini_tandem(const char *read_fn, mini_tandem_para *mtp)
{
	int i, n_seqs;
	gzFile readfp = xzopen(read_fn, "r");
	kstream_t *fs = ks_init(readfp);
	kseq_t *read_seq = (kseq_t*)calloc(CHUNK_READ_N, sizeof(kseq_t));
	for (i = 0; i < CHUNK_READ_N; ++i) read_seq[i].f = fs;
	tandem_seq_t *tseq = alloc_tandem_seq(CHUNK_READ_N);
	/* if (mtp->splint_fn != NULL) {
		gzFile splint_fp = xzopen(mtp->splint_fn, "r");
		mtp->splint_len = get_seq_from_fx(splint_fp, &(mtp->splint_seq));
		mtp->splint_rc_seq = get_rc_seq(mtp->splint_seq, mtp->splint_len);
		err_gzclose(splint_fp);
	}*/
	if (mtp->five_fn != NULL && mtp->three_fn != NULL) {
		gzFile five_fp = xzopen(mtp->five_fn, "r"); gzFile three_fp = xzopen(mtp->three_fn, "r");
		mtp->five_len = get_seq_from_fx(five_fp, &(mtp->five_seq)); mtp->five_rc_seq = get_rc_seq(mtp->five_seq, mtp->five_len);
		mtp->three_len = get_seq_from_fx(three_fp, &(mtp->three_seq)); mtp->three_rc_seq = get_rc_seq(mtp->three_seq, mtp->three_len);
		err_gzclose(five_fp); err_gzclose(three_fp);
	}

	// alloc and initialization for auxiliary data
	if (mtp->n_thread < 1) mtp->n_thread = 1;
	thread_aux_t *aux = aux_init(mtp);

	// core loop
	pthread_rwlock_init(&RWLOCK, NULL);
	while ((n_seqs = mini_tandem_read_seq(read_seq, CHUNK_READ_N)) != 0) {
		THREAD_READ_I = 0;
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
		// output consensus sequences
		mini_tandem_output(n_seqs, read_seq, tseq, mtp);
	}
	pthread_rwlock_destroy(&RWLOCK);

	// free variables
	for (i = 0; i < CHUNK_READ_N; ++i) {
        read_seq_free(read_seq+i);
        tandem_seq_free(tseq+i);
	} 
    free(read_seq); free(tseq); ks_destroy(fs); err_gzclose(readfp); free(aux);
	return 0;
}

// TODO add score para 
int main(int argc, char *argv[])
{
	mini_tandem_para *mtp = mini_tandem_init_para();
	int c, op_idx;
	double realtime0 = realtime();
	while ((c = getopt_long(argc, argv, "k:w:m:s:Hhvc:e:p:P:5:3:a:o:lFf:t:", mini_tandem_opt, &op_idx)) >= 0) {
		switch(c)
		{
            case 'h': return usage();
            case 'v': printf("%s\n", VERSION); goto End; break;

			case 'k': mtp->k = atoi(optarg); 
			          if (mtp->k > MAX_KMER_SIZE) {
			          	  err_printf("\n[main] Error: k-mer length can not be larger than %d (%ld).\n", MAX_KMER_SIZE, mtp->k);
			          	  goto End;
				      }
				      break;
			case 'w': mtp->w = atoi(optarg); break;
			case 'm': mtp->m = atoi(optarg); break;
			case 's': mtp->s = atoi(optarg); break;
			case 'H': mtp->hpc = 1; break;

			case 'c': mtp->min_copy = atoi(optarg); break;
			case 'e': mtp->max_div = atof(optarg); break;
			case 'p': mtp->min_p = th_parse_num(optarg);
					  if (mtp->min_p < MIN_PERIOD) {
						  err_printf("Error: -p --min-period needs to be >= %u. (%ld)\n", MIN_PERIOD, mtp->min_p); 
						  goto End;
					  }
					  break;
			case 'P': mtp->max_p = th_parse_num(optarg); 
					  if (mtp->max_p > MAX_PERIOD) {
						  err_printf("Error: -P --max-period needs to be <= %u. (%ld)\n", MAX_PERIOD, mtp->max_p); 
						  goto End;
					  }
					  break;
		  //case 'r': mtp->max_range = th_parse_num(optarg); break;

			//case 'S': mtp->splint_fn = strdup(optarg); break;
			case '5': mtp->five_fn = strdup(optarg); break;
			case '3': mtp->three_fn = strdup(optarg); break;
			case 'a': mtp->ada_match_rat = atof(optarg); break;

			case 'o': mtp->cons_out = xopen(optarg, "w"); break;
			case 'l': mtp->only_longest = 1; break;
			case 'F': mtp->only_full_length = 1; break;
			case 'f': mtp->out_fmt = atoi(optarg); 
					  if (mtp->out_fmt != FASTA_FMT && mtp->out_fmt != TAB_FMT) {
						  err_printf("\n[main] Error: unknown format number: %d. (-%c --%s)\n", atoi(optarg), c, mini_tandem_opt[op_idx].name);
						  goto End;
					  }
					  break;

			case 't': mtp->n_thread = atoi(optarg); break;

			default:
					  err_printf("\n[main] Error: unknown option: -%c %s.\n", c, optarg);
					  goto End;
		}
	}
	if (optind + 1 > argc) {
		err_fprintf(stderr, "\n[main] Error: please specify an input file.\n");
		usage(); goto End;
	}
	if (mtp->only_full_length && (mtp->five_fn == NULL || mtp->three_fn == NULL)) {
		err_printf("\n[main] Error: 5' and 3' adapter sequence need to be provided.\n");
		usage(); goto End;
	}
	if (mtp->five_fn == NULL && mtp->three_fn != NULL) 
		err_printf("\n[main] Warning: only 3' adapter sequence is provided. Full-length sequence cannot be determined.\n");
	if (mtp->five_fn != NULL && mtp->three_fn == NULL) 
		err_printf("\n[main] Warning: only 5' adapter sequence is provided. Full-length sequence cannot be determined.\n");
	mtp->div_exp = get_div_exp(mtp->k, mtp->max_div);
	mini_tandem(argv[optind], mtp);
	err_func_printf(__func__, "Real time: %.3f sec; CPU: %.3f sec; Peak RSS: %.3f GB\n", realtime() - realtime0, cputime(), peakrss() / 1024.0 / 1024.0 / 1024.0);
End:
	mini_tandem_free_para(mtp);
	return 0;
}
