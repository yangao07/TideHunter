import sys, os
import utils as ut
from pyfaidx import Fasta
import pysam as ps
import argparse
import subprocess
import random
from datetime import datetime
from collections import defaultdict as dd
from itertools import permutations as permu
from Bio.Seq import Seq
import trf_cons_align as tca

fxtools = 'fxtools'
miniTandem = 'miniTandem'
minimap2 = 'minimap2'
trf = 'trf.sh'
pbsim = 'pbsim.sh'  # pbsim.sh in.fa out_prefix type(CLR/CCS) depth len-min len-max acc-min acc-max
model = 'CLR'  # 'CLR'
acc_min = 0.85
acc_max = 0.85

def_read_cnt = 1000
read_cnts = [500, 1000, 5000, 10000]
def_par_len = 2000  # pattern size
par_lens = [500, 1000, 2000, 5000, 10000]  # 10000 TODO
def_copy_num = 10
copy_nums = [5, 10, 20, 50, 100]
def_var_type = 'NOR'  # normal TR
var_types = ['NOR', 'INS', 'DEL', 'RC', 'NON']  # insertion/deletion/reverse-comp/chimeric/non-tandem-repeat
def_run_thread = 1
run_threads = [1, 2, 4, 8, 16]
m_k = 8
m_w = 5
m_s = 1
m_H = '-H'  # '-p'


def eval_miniTandem_sam(sam_fn):
    eval_out = dd(lambda: 0)  # (Can_Be_Corr_Map, Can_Be_Map_To_Other) : Num_of_Cons
    corr_map, wrong_map = 0, 0
    with ps.AlignmentFile(sam_fn) as sam_fp:
        rr = ps.AlignedSegment
        for r in sam_fp:
            if r.is_supplementary or r.is_secondary: continue
            qname, tname = r.query_name, r.reference_name

            if '{}_cons'.format(tname) in qname:
                corr_map = 1
            else:
                wrong_map = 1
            eval_out[(corr_map, wrong_map)] += 1
            corr_map = wrong_map = 0
    return eval_out


# tr: 2 copies of tr units
def eval_miniTandem(fa, tr):
    # run miniTandem
    out_fa = fa+'.mini.fa'
    ut.exec_cmd(sys.stderr, 'miniTandem', '{} {} -k{} -w{} -s{} {} > {}'.format(miniTandem, fa, m_k, m_w, m_s, m_H, out_fa))
    # align result with tr
    sam = '{}.mini.sam'.format(fa)
    ut.exec_cmd(sys.stderr, 'minimap2', '{} -a {} {} > {}'.format(minimap2, tr, out_fa, sam))
    # parse BAM
    eval_out = eval_miniTandem_sam(sam)
    print eval_out
    # write eval_out


def eval_trf(fa, tr):
    # run trf
    trf_out = '{}.trf.out'.format(fa)
    ut.exec_cmd(sys.stderr, 'trf', '{} {} {}'.format(trf, fa, trf_out))
    # extract seq, align with tr
    cons_fa, cons_info, cons_all_sam, cons_best_bam = fa + '.trf.fa', fa + '.trf.info', fa + '.trf.sam', fa + '.trf.best.bam'
    tca.cons_align(cons_fa, cons_info, cons_all_sam, cons_best_bam, trf_out, fa, tr, cons_min_len=tca.cons_min_len,
                   cons_min_copy=tca.cons_min_copy, cons_min_frac=tca.cons_min_frac, minimap=minimap2,
                   threads=tca.threads)
    # parse BAM


def get_pbsim_read(pbsim, in_fa, depth, len, out_pre, seed):
    ut.exec_cmd(sys.stderr, 'get_pbsim_read',
                '{} {} {} {} {} {} {} {} {} {}'.format(pbsim, in_fa, out_pre, model, depth, len, len, acc_min, acc_max, seed))


def pbsim_read_core(in_seq, cnt, out_fn, out_seq_pre):
    out_dir = os.path.dirname(os.path.abspath(out_fn)) + '/tmp/'
    length = len(in_seq)

    if not os.path.exists(out_dir):
        ut.exec_cmd(sys.stderr, 'mkdir', 'mkdir {} 2> /dev/null'.format(out_dir))
    out_pre = out_dir + 'tmp'
    # 0. write in_seq to file
    in_fa = out_pre + '.in.fa'
    with open(in_fa, 'w') as in_fp:
        in_fp.write('>seq\n{}\n'.format(in_seq))
    pbsim_pre = out_pre + '.pbsim'
    get_pbsim_read(pbsim, in_fa, cnt, length, pbsim_pre)
    ut.exec_cmd(sys.stderr, 'merge_pbsim', 'cat {}*.fastq > {}.fq'.format(pbsim_pre, pbsim_pre))
    ut.exec_cmd(sys.stderr, 'fastq2fasta', '{} qa {}.fq > {}.fa'.format(fxtools, pbsim_pre, pbsim_pre))
    ut.exec_cmd(sys.stderr, 'rm', 'rm {}.fq {}*.fastq {}*.maf {}*.ref'.format(pbsim_pre, pbsim_pre, pbsim_pre, pbsim_pre))
    ut.exec_cmd(sys.stderr, 'pbsim.fa', 'head -n {} {}.fa > {}; sed -i \'s/S1/{}/g\' {}'.format(cnt * 2, pbsim_pre, out_fn, out_seq_pre, out_fn))
    ut.exec_cmd(sys.stderr, 'rm tmp', 'rm -rf {}'.format(out_dir))


def pbsim_cat_read(in_seq, out_dir, out_fp, read_name, seed):
    cnt = 1  # only one read is needed
    out_dir += '/tmp/'
    length = len(in_seq)
    if not os.path.exists(out_dir):
        ut.exec_cmd(sys.stderr, 'mkdir', 'mkdir {} 2> /dev/null'.format(out_dir))
    out_pre = out_dir + 'tmp'
    # 0. write in_seq to file
    in_fa = out_pre + '.in.fa'
    with open(in_fa, 'w') as in_fp:
        in_fp.write('>seq\n{}\n'.format(in_seq))
    pbsim_pre = out_pre + '.pbsim'
    get_pbsim_read(pbsim, in_fa, cnt, length, pbsim_pre, seed)
    ut.exec_cmd(sys.stderr, 'merge_pbsim', 'cat {}*.fastq > {}.fq'.format(pbsim_pre, pbsim_pre))
    ut.exec_cmd(sys.stderr, 'fastq2fasta', '{} qa {}.fq > {}.fa'.format(fxtools, pbsim_pre, pbsim_pre))
    ut.exec_cmd(sys.stderr, 'rm', 'rm {}.fq {}*.fastq {}*.maf {}*.ref'.format(pbsim_pre, pbsim_pre, pbsim_pre, pbsim_pre))
    cmd = 'head -n {} {}.fa | tail -n1'.format(cnt * 2, pbsim_pre)
    seq = subprocess.check_output(cmd, shell=True)
    out_fp.write('>{}\n{}'.format(read_name, seq))
    # ut.exec_cmd(sys.stderr, 'pbsim.fa',
    #             'head -n {} {}.fa > {}; sed -i \'s/S1/{}/g\' {}'.format(cnt * 2, pbsim_pre, out_fn, out_seq_pre,
    #                                                                     out_fn))
    ut.exec_cmd(sys.stderr, 'rm tmp', 'rm -rf {}'.format(out_dir))


def random_seq(in_fa, length=def_par_len):
    chr = in_fa.keys()  # ['chr1', 'chr2']
    n_chr = len(chr)
    len_chr = [len(in_fa[c][0:].seq) for c in chr]

    while True:
        chr_i = random.randint(0, n_chr - 1)
        if length <= len_chr[chr_i]:
            break
    # print 'random_chr: {}'.format(chr_i)
    while True:
        start = random.randint(1, len_chr[chr_i] - length + 1)
        seq = in_fa.get_seq(chr[chr_i], start, start + length - 1).seq
        if 'N' not in seq:
            break
    # print 'random_start: {}'.format(start)
    return seq, chr_i, start, start + length - 1


def make_tr(in_fa, plen=def_par_len, copy_num=def_copy_num, var_type=def_var_type):
    if var_type == 'NOR':
        seq, chr, start, end = random_seq(in_fa, plen)
        return seq * copy_num, seq, chr, start, end
    elif var_type == 'INS':
        ins_seq, chr, start, end = random_seq(in_fa, plen / 2)
        seq, chr, start, end = random_seq(in_fa, plen)
        ins_pos = int(plen * random.uniform(1.0, copy_num - 1))
        tr_seq = seq * copy_num
        return tr_seq[:ins_pos] + ins_seq + tr_seq[ins_pos:], seq, chr, start, end
    elif var_type == 'DEL':
        seq, chr, start, end = random_seq(in_fa, plen)
        del_pos = int(plen * random.uniform(1.0, copy_num - 1))
        tr_seq = seq * copy_num
        return tr_seq[:del_pos] + tr_seq[del_pos + plen / 2:], seq, chr, start, end
    elif var_type == 'RC':
        seq, chr, start, end = random_seq(in_fa, plen)
        tr_seq = seq * copy_num
        rc_seq = str(Seq(seq).reverse_complement()) * copy_num
        return tr_seq + rc_seq, seq, chr, start, end
    elif var_type == 'NON':
        seq, chr, start, end = random_seq(in_fa, plen * copy_num)
        return seq, seq, chr, start, end
    else:
        ut.fatal_format_time("make_var_type", 'Unknown variance type: {}.'.format(var_type))


def write_seq(seq='', fp=None, name=''):
    fp.write('>{}\n{}\n'.format(name, seq))
    return


def trsim_core(in_fa, out_fa='', out_dir='', cnt=def_read_cnt, plen=def_par_len, copy_num=def_copy_num,
               var_type=def_var_type, seed=0):
    fa_out_fn, tr_out_fn = '{}/{}'.format(out_dir, out_fa), '{}/{}.tr'.format(out_dir, out_fa)
    for out_fn in [fa_out_fn, tr_out_fn]:
        if os.path.exists(out_fn):
            ut.exec_cmd(sys.stderr, 'rm', 'rm {}'.format(out_fn))

    with open(fa_out_fn, 'w') as out_fp, open('{}/{}'.format(out_dir, out_fa + '.tr'), 'w') as tr_fp:
        for read_i in range(cnt):
            tr_seq, unit_seq, chr, start, end = make_tr(in_fa, plen, copy_num, var_type)
            tr_name = '{}_{}_{}:{}-{}'.format(var_type, read_i, in_fa.keys()[chr], start, end)
            pbsim_cat_read(tr_seq, out_dir, out_fp, tr_name, seed)
            tr_fp.write('>{}\n{}\n'.format(tr_name, unit_seq * 2))
    return


def main(args):
    in_fa, sim_fa, out_dir, cnt, plen, copy_num, var_type = Fasta(args.in_fa), args.sim_out, args.out_dir, args.read_count, args.period, args.copy_num, args.var_type
    out_dir = os.path.abspath(out_dir)
    seed = args.random_seed
    # 1. simulate tandem-repeat sequence read
    sim_tr = sim_fa + '.tr'
    trsim_core(in_fa, sim_fa, out_dir, cnt, plen, copy_num + 1, var_type, seed)
    # 2. run miniTandem, TRF
    eval_miniTandem(out_dir+'/'+sim_fa, out_dir+'/'+sim_tr)
    eval_trf(out_dir+'/'+sim_fa, out_dir+'/'+sim_tr)

def parser_argv():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description="TRSim: tandem repeat sequence simulator based on PBSIM")
    parser.add_argument("in_fa", metavar='in.fa', type=str, help='Input read file in FASTA/FASTQ format.')
    parser.add_argument("sim_out", metavar='sim.fa', type=str, help='Simulated reads in FASTA format.')

    general_par = parser.add_argument_group('General options')
    general_par.add_argument('-o', '--out-dir', metavar='./output', type=str, default='.', help='Ouptut folder.')

    sim_par = parser.add_argument_group('Simulation options')
    sim_par.add_argument('-c', '--read-count', metavar='C', type=int, default=def_read_cnt, help='Read count.')
    sim_par.add_argument('-p', '--period', metavar='P', type=int, default=def_par_len, help='Pattern length.')
    sim_par.add_argument('-n', '--copy-num', metavar='N', type=int, default=def_copy_num, help='Copy number.')
    sim_par.add_argument('-v', '--var-type', metavar='V', type=str, default=def_var_type, choices=var_types,
                         help='Variance type.')
    sim_par.add_argument('-s', '--random-seed', metavar='S', type=int, default=0,
                         help='Random seed. (Default: date time)')

    # eval_par = parser.add_argument_group('Evaluation options')
    # eval_par.add_argument('-t', '--threads', type=int, default=def_run_thread, choices=run_threads, help='Number of threads to use.')

    return parser.parse_args()


if __name__ == '__main__':
    # pbsim_read_core('A' * 500, 10, './pbsim.fa', 'new')
    args = parser_argv()
    if not os.path.exists(args.out_dir):
        ut.exec_cmd(sys.stderr, 'mkdir', 'mkdir {}'.format(args.out_dir))
    if args.random_seed == 0:
        random.seed(datetime.now())
    else:
        random.seed(args.random_seed)
        print 'seed: {}'.format(args.random_seed)
    main(args)
