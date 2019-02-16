# TideHunter: efficient and sensitive tandem repeat detection from noisy long-reads using seed-and-chain

## Getting started
	git clone https://github.com/yangao07/TideHunter.git --recursive
	cd TideHunter; make
	./bin/TideHunter ./test_data/test_50x4.fa > cons.fa

## Table of Contents

- [Introduction](#introduction)
- [Installation](#install)
  - [Operating system](#os)
  - [Cloning and building TideHunter](#build)
- [Getting started with toy example in `test_data`](#start)
- [Usage](#usage)
  - [Generate consensus in FASTA format](#fasta_cons)
  - [Generate consensus in tabular formatr](#tab_cons)
  - [Generate a full-length consensus](#full_cons)
- [Commands and options](#cmd)
- [Input](#input)
  - [Adapter sequence](#adapter)
- [Output](#output)
  - [Tabular format](#tabular)
  - [FASTA format](#fasta)
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
./bin/TideHunter ./test_data/test_1000x10.fa > cons.fa
```

## <a name="usage"></a>Usage
### <a name="fasta_cons"></a>Generate consensus in FASTA format
```
./bin/TideHunter ./test_data/test_1000x10.fa > cons.fa
```
### <a name="tab_cons"></a>Generate consensus in tabular format
```
./bin/TideHunter -f 2 ./test_data/test_1000x10.fa > cons.out
```
### <a name="full_cons"></a>Generate a full-length consensus
```
./bin/TideHunter -5 ./test_data/5prime.fa -3 ./test_data/3prime.fa ./test_data/full_length.fa > cons_full.fa
```

## <a name="cmd"></a>Commands and options
```
Usage:   TideHunter [options] in.fa/fq > cons_out.fa

Options:
    Seeding:
         -k --kmer-length [INT]    k-mer length (no larger than 16). [8]
         -w --window-size [INT]    window size. [1]
         -s --step-size   [INT]    step size. [1]
         -H --HPC-kmer             use homopolymer-compressed k-mer. [False]
    Tandem repeat criteria:
         -c --min-copy    [INT]    minimum copy number of tandem-repeats. [2]
         -e --max-diverg  [INT]    maximum allowed divergence rate between two consecutive repeats. [0.25]
         -p --min-period  [INT]    minimum period size of tandem repeat. (>=2) [30]
         -P --max-period  [INT]    maximum period size of tandem repeat. (<=65535) [65535]
    Adapter sequence:
         -5 --five-prime  [STR]    5' adapter sequence (sense strand). [NULL]
         -3 --three-prime [STR]    3' adapter sequence (anti-sense strand). [NULL]
         -a --ada-mat-rat [FLT]    minimum match ratio of adapter sequence. [0.80]
    Output:
         -o --cons-out    [STR]    output consensus sequence in FASTA format. [stdout]
         -l --longest              only output the consensus of the longest tandem repeat. [False]
         -F --full-len             only output the consensus that is full-length. [False]
         -f --out-fmt     [INT]    output format. [1]
                                       1: FASTA
                                       2: Tabular
    Computing resource:
         -t --thread      [INT]    number of threads to use. [1]

```

## <a name="input_output"></a>Input
TideHunter works with FASTA, FASTQ, gzip'd FASTA(.fa.gz) and gzip'd FASTQ(.fq.gz) formats.

### <a name="adapter"></a>Adapter sequence
Additional adapter sequence files can be provided to TideHunter with `-5` and `-3` options.

TideHunter uses adapter information to search for the full-length sequence from the generated consensus.

Once two adapters are found, TideHunter trims and reorients the consensus sequence.

## <a name="output"></a>Output
TideHunter can output consensus sequence in FASTA format by default, 
it can also provide output in tabular format.

### <a name="tabular"></a>Tabular format
For tabular format, 9 columns will be generated for each consensus sequence:

| id  | Column name | Explanation | 
|:---:|   :---      | ---        |
|  1  | readName    | the original read name |
|  2  | consN       | `N` is the ID number of the consensus sequences from the same read, starts from 0 |
|  3  | readLen     | length of the original long-read |
|  4  | start       | start coordinate of the tandem repeat, 1-base |
|  5  | end         | end coordinate of the tandem repeat, 1-base |
|  6  | consLen     | length of the consensus sequence |
|  7  | copyNum     | copy number of the tandem repeat |
|  8  | fullLen     | 0: not a full-length sequence, 1: sense strand full-length, 2: anti-sense strand full-length |
|  9  | consensus   | consensus sequence |

### <a name="fasta"></a>FASTA format
For FASTA output format, the read name contains detailed information of the detected tandem repeat, 
i.e., the above columns 1 ~ 8.
The sequence is the consensus sequence.

The read name of each consensus sequence has the following format:
```
>readName_consN_readLen_start_end_consLen_copyNum_fullLen
```

## <a name="contact"></a>Contact
Yan Gao yangao07@hit.edu.cn

Yadong Wang ydwang@hit.edu.cn

[github issues](https://github.com/yangao07/TideHunter/issues)
