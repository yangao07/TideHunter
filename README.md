# TideHunter: efficient and sensitive tandem repeat detection from noisy long-reads using seed-and-chain

## Getting started
	git clone https://github.com/yangao07/TideHunter.git --recursive
	cd TideHunter; make
	./bin/TideHunter ./test_data/in.fa > cons.fa

## Table of Contents

- [Introduction](#introduction)
- [Installation](#install)
  - [Operating system](#os)
  - [Cloning and building TideHunter](#build)
- [Getting started with toy example in `test_data`](#start)
- [Commands and options](#cmd)
- [Input and output](#input_output)
  - [Detailed tandem repeat information](#information)
- [Contact](#contact)

## <a name="introduction"></a>Introduction
TideHunter is an efficient and sensitive tandem repeat detection and
consensus calling tool which is designed for tandemly repeated
long-read sequence ([INC-seq](https://doi.org/10.1186/s13742-016-0140-7),
[R2C2](https://doi.org/10.1073/pnas.1806447115)). 

It works with Pacific Biosciences (PacBio) and 
Oxford Nanopore Technologies (ONT) sequencing data at error rates 
up to 20% and is able to detect repeat patterns of any size.

## <a name="install"></a>Installation
### <a name="os"></a>Operating system
TideHunter currently can only be built and run on Linux/Unix systems.

### <a name="build"></a>Cloning and building TideHunter
```
git clone https://github.com/yangao07/TideHunter.git --recursive
cd TideHunter; make
```

## <a name="start"></a>Getting started with toy example in `test_data`
```
./bin/TideHunter ./test_data/test.fa > cons.fa
```

## <a name="cmd"></a>Commands and options
```
Usage:   TideHunter [options] in.fa/fq > cons_out.fa

Options: 
         -t --thread      [INT]    number of threads to use. [1]
         -k --kmer-length [INT]    k-mer length (no larger than 16). [8]
         -s --step-size   [INT]    step size. [1]
         -w --window-size [INT]    window size. [1]
         -H --HPC-kmer             use homopolymer-compressed k-mer. [False]
         -c --min-copy    [INT]    minimum copy number of tandem-repeats. [2]
         -e --max-diverg  [INT]    maximum allowed divergence rate between two consecutive repeats. [0.25]
         -p --min-period  [INT]    minimum period size of tandem repeat. (>=2) [30]
         -P --max-period  [INT]    maximum period size of tandem repeat. (<=65535) [65535]
         -l --longest              only output the consensus of the longest tandem repeat. [False]
         -O --cons-out    [STR]    output consensus sequence in FASTA format. [stdout]

```

## <a name="input_output"></a>Input and output
TideHunter works with FASTA, FASTQ, gzip'd FASTA(.fa.gz) and gzip'd FASTQ(.fq.gz) formats.

The output consensus sequence file is in FASTA format.

### <a name="information"></a>Detailed tandem repeat information 
For each consensus output, TideHunter appends the tandem repeat information 
to the consensus name in the following format:

```
>readName_consN_readLen:start:end:consLen:copyNum
Consensus sequence
```
`readName`: the original read name from input file \
`N`: the ID number of the consensus sequences from the same read, starts from 0\
`readLen`: length of the original long-read\
`start`: start coordinate of the tandem repeat, 1-base\
`end`: end coordinate of the tandem repeat, 1-base\
`consLen`: length of the consensus sequence\
`copyNum`:  copy number of the tandem repeat\
`Consensus sequences`: consensus sequence generated with partial order alignment


## <a name="contact"></a>Contact
Yan Gao yangao07@hit.edu.cn

Yadong Wang ydwang@hit.edu.cn

[github issues](https://github.com/yangao07/TideHunter/issues)
