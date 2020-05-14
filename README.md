# TideHunter: efficient and sensitive tandem repeat detection from noisy long reads using seed-and-chain
[![Latest Release](https://img.shields.io/github/release/yangao07/TideHunter.svg?label=Release)](https://github.com/yangao07/TideHunter/releases/latest)
[![Github All Releases](https://img.shields.io/github/downloads/yangao07/TideHunter/total.svg?label=Download)](https://github.com/yangao07/TideHunter/releases)
[![BioConda Install](https://img.shields.io/conda/dn/bioconda/tidehunter.svg?style=flag&label=BioConda%20install)](https://anaconda.org/bioconda/tidehunter)
[![Published in Bioinformatics](https://img.shields.io/badge/Published%20in-Bioinformatics-blue.svg)](https://doi.org/10.1093/bioinformatics/btz376)
[![GitHub Issues](https://img.shields.io/github/issues/yangao07/TideHunter.svg?label=Issues)](https://github.com/yangao07/TideHunter/issues)
[![Build Status](https://img.shields.io/travis/yangao07/TideHunter/master.svg?label=Master)](https://travis-ci.org/yangao07/TideHunter)
[![License](https://img.shields.io/badge/License-MIT-black.svg)](https://github.com/yangao07/TideHunter/blob/master/LICENSE)
<!--
[![GitHub Downloads](https://img.shields.io/github/downloads/yangao07/TideHunter/total.svg?style=social&logo=github&label=Download)](https://github.com/yangao07/TideHunter/releases)
-->

## Getting started
Download the [latest release](https://github.com/yangao07/TideHunter/releases):
```
wget https://github.com/yangao07/TideHunter/releases/download/v1.3.0/TideHunter-v1.3.0.tar.gz
tar -zxvf TideHunter-v1.3.0.tar.gz && cd TideHunter-v1.3.0
```
Make from source and run with test data:
```
make; ./bin/TideHunter ./test_data/test_50x4.fa > cons.fa
```
Or, install via conda and run with test data:
```
conda install -c bioconda tidehunter
TideHunter ./test_data/test_50x4.fa > cons.fa
```
## Table of Contents

- [Introduction](#introduction)
- [Installation](#install)
  - [Installing TideHunter via conda](#conda)
  - [Building TideHunter from source files](#build)
  - [Pre-built binary executable file for Linux/Unix](#binary)
- [Getting started with toy example in `test_data`](#start)
- [Usage](#usage)
  - [To generate consensus in FASTA format](#fasta_cons)
  - [To generate consensus in tabular format](#tab_cons)
  - [To generate a full-length consensus](#full_cons)
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
 [R2C2](https://doi.org/10.1073/pnas.1806447115), [NanoAmpli-Seq](https://doi.org/10.1093/gigascience/giy140)). 

It works with Pacific Biosciences (PacBio) and 
Oxford Nanopore Technologies (ONT) sequencing data at error rates 
up to 20% and does not have any limitation of the maximal repeat pattern size.

## <a name="install"></a>Installation

### <a name="conda"></a>Installing TideHunter via conda
On Linux/Unix and Mac OS, TideHunter can be installed via
```
conda install -c bioconda tidehunter
```

### <a name="build"></a>Building TideHunter from source files
You can also build TideHunter from source files.
Make sure you have gcc (>=6.4.0) and zlib installed before compiling.
It is recommended to download the latest release of TideHunter 
from the [release page](https://github.com/yangao07/TideHunter/releases).
```
wget https://github.com/yangao07/TideHunter/releases/download/v1.3.0/TideHunter-v1.3.0.tar.gz
tar -zxvf TideHunter-v1.3.0.tar.gz
cd TideHunter-v1.3.0; make
```
Or, you can use `git clone` command to download the source code. 
Don't forget to include the `--recursive` to download the codes of [abPOA](https://github.com/yangao07/abPOA).
This gives you the latest version of TideHunter, which might be still under development.
```
git clone --recursive https://github.com/yangao07/TideHunter.git
cd TideHunter; make
```

### <a name="binary"></a>Pre-built binary executable file for Linux/Unix 
If you meet any compiling issue, please try the pre-built binary file:
```
wget https://github.com/yangao07/TideHunter/releases/download/v1.3.0/TideHunter-v1.3.0_x64-linux.tar.gz
tar -zxvf TideHunter-v1.3.0_x64-linux.tar.gz
```

## <a name="start"></a>Getting started with toy example in `test_data`
```
TideHunter ./test_data/test_1000x10.fa > cons.fa
```

## <a name="usage"></a>Usage
### <a name="fasta_cons"></a>To generate consensus in FASTA format
```
TideHunter ./test_data/test_1000x10.fa > cons.fa
```
### <a name="tab_cons"></a>To generate consensus in tabular format
```
TideHunter -f 2 ./test_data/test_1000x10.fa > cons.out
```
### <a name="full_cons"></a>To generate a full-length consensus sequence
```
TideHunter -5 ./test_data/5prime.fa -3 ./test_data/3prime.fa ./test_data/full_length.fa > cons_full.fa
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
         -P --max-period  [INT]    maximum period size of tandem repeat. (<=4294967295) [10K]
    Adapter sequence:
         -5 --five-prime  [STR]    5' adapter sequence (sense strand). [NULL]
         -3 --three-prime [STR]    3' adapter sequence (anti-sense strand). [NULL]
         -a --ada-mat-rat [FLT]    minimum match ratio of adapter sequence. [0.80]
    Output:
         -o --cons-out    [STR]    output file. [stdout]
         -l --longest              only output the consensus sequence of tandem repeat that covers the longest read sequence. [False]
         -F --full-len             only output the consensus sequence that is full-length. [False]
         -f --out-fmt     [INT]    output format. [1]
                                       1: FASTA
                                       2: Tabular
    Computing resource:
         -t --thread      [INT]    number of threads to use. [1]

    General options:
         -h --help                 print this help usage information.

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

| No. | Column name | Explanation | 
|:---:|   :---      | ---        |
|  1  | readName    | the original read name |
|  2  | consN       | `N` is the ID number of the consensus sequences from the same read, starts from 0 |
|  3  | readLen     | length of the original long read |
|  4  | start       | start coordinate of the tandem repeat, 1-based |
|  5  | end         | end coordinate of the tandem repeat, 1-based |
|  6  | consLen     | length of the consensus sequence |
|  7  | copyNum     | copy number of the tandem repeat |
|  8  | fullLen     | 0: not a full-length sequence, 1: sense strand full-length, 2: anti-sense strand full-length |
|  9  | subPos      | start coordinates of all the tandem repeat unit sequence, followed by the end coordinate of the last tandem repeat unit sequence, separated by `,`, all coordinates are 1-based, see examples below|
| 10  | consensus   | consensus sequence |

For example, here are the output for a non-full-length consensus sequence generated from [test_data/test_50x4.fa](test_data/test_50x4.fa) and the adiagram that illustrates all the coordiantes in the output:
```
test_50x4 cons0 300 51  250 50  4.0 0 59,109,159,208  CGATCGATCGGCATGCATGCATGCTAGTCGATGCATCGGGATCAGCTAGT
```
<!-- ![example](example_50x4.png) -->
<img src="example_50x4.png" width="800">

In this example, TideHunter identifies three consecutive complete tandem repeat units, [59,108], [109,158], [159,208], from the raw read which is 300 bp long.
A consensus sequence with 50 bp is generated from the three complete repeat units. TideHunter further extends the tandem repeat boundary to [51, 250] by aligning the consensus sequence back to the raw read on both sides of the three complete repeat units.

Another example of the output for a full-length consensus sequence generated from [test_data/full_length.fa](test_data/full_length.fa):
```
8f2f7766-4b8e-4c0d-9e2b-caf0e5527b19  cons0 5231  31  5215  203 8.8 1 207,798,1386,1976,2563,3155,3746,4333,4930  ACTAATAAGATCAACAGAATCAGAGTAGATAGTTCCTTGATCGGAACCAAAGGACCCCGTGCCTCAATCTCTATCCTGATGTCATGGGAGTCCTAGCAAAGCTATAGACTCAAGCAAGGCTTGGGGTCCTTTATGGAACCCAAGGATGACTCAGCAATAAAATATTTTGGTTTTGGTTTATAAAAAAAAAAAAAAAAAAAAAA
```
In this example, the `consLen` (i.e., 203) is the length of the full-length consensus sequence excluding the 5' and 3' primer sequences and the `subPos` (i.e., 207,798,1386,1976,2563,3155,3746,4333,4930) is the coordinate information of the identified complete tandem repeat units.

### <a name="fasta"></a>FASTA format
For FASTA output format, the read name contains detailed information of the detected tandem repeat, 
i.e., the above columns 1 ~ 9.
The sequence is the consensus sequence.

The read name of each consensus sequence has the following format:
```
>readName_consN_readLen_start_end_consLen_copyNum_fullLen_subPos
```

## <a name="contact"></a>Contact
Yan Gao yangao07@hit.edu.cn

Yadong Wang ydwang@hit.edu.cn

Yi Xing XINGYI@email.chop.edu

[github issues](https://github.com/yangao07/TideHunter/issues)
