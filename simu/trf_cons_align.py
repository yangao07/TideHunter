import argparse
import os
import sys
import utils as ut

cons_header_ele = ['repeat_start', 'repeat_end', 'period_size', 'copy_number', 'consensus_length', 'match_frac',
                   'indel_frac', 'align_score', 'A', 'C', 'G', 'T', 'entropy', 'cons_seq', 'repeat_seq',
                   'left_flank_seq', 'right_flank_seq']
cons_header_idx = {cons_header_ele[i]: i for i in range(len(cons_header_ele))}

minimap = 'minimap2' #dir_path + '/bin/minimap2'
fxtools = 'fxtools' #dir_path + '/bin/fxtools'
samtools = 'samtools' #dir_path + '/bin/samtools'
cons_min_len = 30
cons_min_copy = 2.0
cons_min_frac = 0.0
threads = 8


def get_read_len(read_fn, out_dir, fxtools=fxtools):
    read_len = {}
    read_len_fn = out_dir + '/' + os.path.basename(read_fn) + '.len'
    ut.exec_cmd(sys.stderr, 'fxtools', '{} lp {} > {}'.format(fxtools, read_fn, read_len_fn))
    with open(read_len_fn, 'r') as len_fp:
        for line in len_fp:
            read_len[line.rsplit()[0]] = int(line.rsplit()[1])
    return read_len


def extract_cons(cons_fa, cons_info, trf_out, read_len, cons_min_len=cons_min_len, cons_min_copy=cons_min_copy, cons_min_frac=cons_min_frac):
    read_name, cons_seq, max_copy_number = '', '', 0
    all_cons_len = 0

    with open(trf_out, 'r') as trf_fp, open(cons_fa, 'w') as cons_fp, open(cons_info, 'w') as info_fp:
        for line in trf_fp:
            ele = line.rsplit()
            if line.startswith('@'):
                if read_name and min_cons_len >= cons_min_len and min_cons_len != 10000:
                    # 1 copy of consensus sequence
                    info_fp.write('{}\t{}\t{}\t{}\t{}\n'.format(read_name, read_len[read_name], len(cons_seq), max_copy_number, all_cons_len / read_len[read_name]))
                    # cons_fp.write('>{}\n{}{}\n'.format(read_name, cons_seq, cons_seq))
                    cons_fp.write('>{}\n{}\n'.format(read_name, cons_seq))
                read_name = line[1:].rsplit()[0]
                min_cons_len = 10000
                cons_seq = ''
            else:  # repeat
                if int(ele[cons_header_idx['consensus_length']]) < min_cons_len \
                        and int(ele[cons_header_idx['consensus_length']]) >= cons_min_len \
                        and read_len[read_name] * cons_min_frac <= (1.0 + float(
                    ele[cons_header_idx['repeat_end']]) - float(ele[cons_header_idx['repeat_start']])) \
                        and float(ele[cons_header_idx['copy_number']]) >= cons_min_copy:
                    all_cons_len = (1.0 + float(ele[cons_header_idx['repeat_end']]) - float(ele[cons_header_idx['repeat_start']]))
                    max_copy_number = float(ele[cons_header_idx['copy_number']])
                    cons_seq = ele[cons_header_idx['cons_seq']]
                    min_cons_len = int(ele[cons_header_idx['consensus_length']])
        if read_name and min_cons_len >= cons_min_len and min_cons_len != 10000:
            info_fp.write('{}\t{}\t{}\t{}\t{}\n'.format(read_name, read_len[read_name], len(cons_seq), max_copy_number, all_cons_len / read_len[read_name]))
            # cons_fp.write('>{}\n{}{}\n'.format(read_name, cons_seq, cons_seq))
            cons_fp.write('>{}\n{}\n'.format(read_name, cons_seq))


# cons_fa: 2 copies of consensus sequence
def align_cons(cons_all_sam, cons_best_bam, cons_fa, ref_fa, minimap=minimap, threads=threads):
    ut.exec_cmd(sys.stderr, 'Mapping', '{} -ax splice -ub --MD {} {} -t {} > {}'.format(minimap, ref_fa, cons_fa, threads, cons_all_sam))
    ut.exec_cmd(sys.stderr, 'FilterBAM', '{} view -bF 4 {} | {} view -bF 2048 | {} view -bF 256 | {} sort -@ {} > {}'.format(samtools, cons_all_sam, samtools, samtools, samtools, threads, cons_best_bam))
    

def cons_align(cons_fa, cons_info, cons_all_sam, cons_best_bam, trf_out, long_read_fn, ref_fa, cons_min_len=cons_min_len, cons_min_copy=cons_min_copy, cons_min_frac=cons_min_frac, minimap=minimap, threads=threads):
    read_len = get_read_len(long_read_fn, os.path.dirname(cons_fa))

    extract_cons(cons_fa, cons_info, trf_out, read_len, cons_min_len, cons_min_copy, cons_min_frac)
    align_cons(cons_all_sam, cons_best_bam, cons_fa, ref_fa, minimap, threads)


def cons_align_core(args):
    cons_align(args.cons_fa, args.cons_info, args.cons_all_sam, args.cons_best_bam, args.trf_out, args.long_read, args.ref_fa, args.min_len, args.min_copy, args.min_frac, args.minimap, args.threads)


def parser_argv():
    # parse command line arguments
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="Extract consensus sequence and align it to genome")
    parser.add_argument("cons_fa", metavar='long.fa', type=str, help='Consensus sequence file.')
    parser.add_argument("cons_info", metavar='long.fa', type=str, help='Consensus information file.')
    parser.add_argument("cons_all_sam", metavar='long.fa', type=str, help='Consensus SAM alignment file.')
    parser.add_argument("cons_best_bam", metavar='long.fa', type=str, help='Sorted best consensus BAM alignment file.')


    parser.add_argument("trf_out", metavar='trf.out', type=str, help='Consensus sequence output file from TRF.')
    parser.add_argument("in_long_read", metavar='long.fa', type=str, help='Original long read.')
    parser.add_argument("ref_fa", type=str, metavar='ref.fa', help='Reference genome sequence file.')

    parser.add_argument('-x', '--fxtools', type=str, help='Path to fxtools.', default=fxtools)
    parser.add_argument('-l', '--min-len', type=int, help='Minimum consensus length to keep.', default=cons_min_len)
    parser.add_argument('-c', '--min-copy', type=float, help='Minimum copy number of consensus to keep.', default=cons_min_copy)
    parser.add_argument('-f', '--min-frac', type=float, help='Minimum fraction of original long read to keep.', default=cons_min_frac)

    parser.add_argument('-m', '--minimap', type=str, help='Path to minimap2.', default=minimap)
    parser.add_argument('-t', '--threads', type=int, help='Number of threads.', default=threads)

    return parser.parse_args()


if __name__ == '__main__':
    args = parser_argv()

    cons_align_core(args)
