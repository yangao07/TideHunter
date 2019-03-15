import sys
import os
import re
import argparse
from pyfaidx import Fasta

def fa_core(in_fa, out_fn, out_type, only_full):
	with Fasta(in_fa) as in_fp, open(out_fn, 'w') as out_fp:
		for r in in_fp:
			name = r.name
			ele = name.rsplit('_')
			if only_full and ele[-1] == '0': continue
			if out_type == 'tab':
				out_fp.write('{}\t{}\t{}\n'.format('_'.join(ele[:-7]), '\t'.join(ele[-7:]), str(r)))
			elif out_type == 'fa':
				out_fp.write('>{}\n{}\n'.format(r.long_name, str(r)))


def tab_core(in_tab, out_fn, out_type, only_full):
	with open(in_tab) as in_fp, open(out_fn, 'w') as out_fp:
		for line in in_fp:
			# print line
			ele = line.rsplit()
			# print len(ele), ele
			if only_full and ele[7] == '0': continue
			if out_type == 'fa':
				out_fp.write('>{}\n{}\n'.format('_'.join(ele[:-1]), ele[-1]))
			elif out_type == 'tab':
				out_fp.write(line)

def th_transform(args):
	if args.in_type == 'tab':
			tab_core(args.in_file, args.out_file, args.out_type, args.only_full)
	elif args.in_type == 'fa':
			fa_core(args.in_file, args.out_file, args.out_type, args.only_full)


# parse command line arguments
def parser_argv():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description="TideHunter output format transformation. (Tabular/FASTA)")
    parser.add_argument('in_file', metavar='in.tab/fa', type=str, help="Tabular/FASTA output file of TideHunber.")
    parser.add_argument('out_file', metavar='out.tab/fa', type=str, help="FASTA/Tabular output file of TideHunber.")
    parser.add_argument('-it', '--in-type', metavar='in_type', default='tab', type=str, choices=['tab','fa'], help="File type of input file.")
    parser.add_argument('-ot', '--out-type', metavar='out_type', default='fa', type=str, choices=['tab','fa'], help="File type of output file.")
    parser.add_argument('-F', '--only-full', default=False, action='store_true', help='Only output the consensus that is full-length.')
    return parser.parse_args()

def main():
    args = parser_argv()
    th_transform(args)

if __name__ == '__main__':
    main()